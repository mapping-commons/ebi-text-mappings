#!/usr/bin/env python3
"""Download zooma-format source data and convert to SSSOM/TSV files.

Uses bioregistry for URI-to-CURIE normalisation with preferred prefix casing.
Uses polars for efficient in-memory processing of large datasets.
Follows the SSSOM 1.0 spec for literal mappings (subject_type: rdfs literal).

Usage:
    python scripts/build_sssom.py                  # download all sources
    python scripts/build_sssom.py --local-data data # use local files
"""

import argparse
from datetime import datetime
import hashlib
import html
import io
import sys
import tempfile
from pathlib import Path

import bioregistry
import polars as pl
import requests
import yaml

# Built-in SSSOM prefixes — may be omitted from curie_map per spec
BUILTIN_PREFIXES = frozenset(
    {"owl", "rdf", "rdfs", "semapv", "skos", "sssom", "xsd", "linkml"}
)

MAPPING_SET_ID_BASE = (
    "https://raw.githubusercontent.com/mapping-commons/"
    "ebi-text-mappings/main/mappings"
)

SCRIPT_DIR = Path(__file__).resolve().parent
ROOT_DIR = SCRIPT_DIR.parent
MAPPINGS_DIR = ROOT_DIR / "mappings"
SOURCES_FILE = ROOT_DIR / "sources.yml"


# ---------------------------------------------------------------------------
# helpers
# ---------------------------------------------------------------------------

def load_sources():
    with open(SOURCES_FILE) as fh:
        return yaml.safe_load(fh)


def download_file(url, dest_dir):
    """Download *url* into *dest_dir*, return path."""
    filename = url.rsplit("/", 1)[-1].split("?")[0]
    dest = Path(dest_dir) / filename
    print(f"  Downloading {url}")
    resp = requests.get(url, stream=True, timeout=600)
    resp.raise_for_status()
    with open(dest, "wb") as fh:
        for chunk in resp.iter_content(chunk_size=65536):
            fh.write(chunk)
    return dest


# ---------------------------------------------------------------------------
# CURIE handling via bioregistry
# ---------------------------------------------------------------------------

_curie_cache: dict[str, str | None] = {}
# Maps CURIE prefix → URI expansion actually seen in source data.
# Used to build an accurate curie_map (bioregistry's get_uri_prefix can return
# a different URL than the one used in the source data).
_prefix_uri_expansions: dict[str, str] = {}


def uri_to_curie(uri):
    """Convert a URI to a CURIE with preferred-prefix casing.

    Returns *None* when bioregistry cannot resolve the URI.
    Also records the URI prefix actually used so that the curie_map
    accurately reflects the source data.
    """
    if not uri:
        return None
    uri = uri.strip()
    if not uri:
        return None
    if uri in _curie_cache:
        return _curie_cache[uri]

    curie = bioregistry.curie_from_iri(uri)
    if curie is None:
        _curie_cache[uri] = None
        return None
    prefix, identifier = curie.split(":", 1)
    preferred = bioregistry.get_preferred_prefix(prefix)
    if preferred:
        prefix = preferred
    result = f"{prefix}:{identifier}"
    _curie_cache[uri] = result

    # Record the URI expansion derived from the actual source URI
    if prefix not in _prefix_uri_expansions:
        uri_expansion = uri[: len(uri) - len(identifier)]
        _prefix_uri_expansions[prefix] = uri_expansion

    return result


def get_uri_expansion(prefix):
    """Return the URI expansion for *prefix* from bioregistry (fallback)."""
    result = bioregistry.get_uri_prefix(prefix)
    if result:
        return result
    return bioregistry.get_uri_prefix(prefix.lower())


# ---------------------------------------------------------------------------
# parsing source files
# ---------------------------------------------------------------------------

def detect_delimiter(filepath):
    """Sniff first line: if it contains tabs, assume TSV; otherwise CSV."""
    with open(filepath, encoding="utf-8", errors="replace") as fh:
        first_line = fh.readline()
    return "\t" if "\t" in first_line else ","


def parse_source(filepath, source_config):
    """Parse a zooma-format (or GWAS-format) source file.

    Returns a polars DataFrame with columns: subject_label, semantic_tag,
    and optionally property_type, object_label, zooma:study, zooma:bioentity,
    author_label, mapping_date.
    """
    delimiter = detect_delimiter(filepath)
    column_map = source_config.get("column_map", {})
    object_label_col = source_config.get("object_label_column")

    read_opts = dict(
        separator=delimiter,
        infer_schema_length=0,
        truncate_ragged_lines=True,
        quote_char=None,
        encoding="utf8-lossy",
    )

    # Check for \r line endings; normalize if found (handles CR-only and mixed)
    with open(filepath, "rb") as fh:
        probe = fh.read(65536)
    if b"\r" in probe:
        with open(filepath, "rb") as fh:
            raw = fh.read()
        raw = raw.replace(b"\r\n", b"\n").replace(b"\r", b"\n")
        df = pl.read_csv(io.BytesIO(raw), **read_opts)
    else:
        df = pl.read_csv(filepath, **read_opts)

    df = df.fill_null("")

    if column_map:
        # Custom column mapping (e.g. GWAS)
        exprs = []
        for src_col, dest_col in column_map.items():
            if dest_col == "PROPERTY_VALUE":
                exprs.append(pl.col(src_col).str.strip_chars().alias("subject_label"))
            elif dest_col == "SEMANTIC_TAG":
                exprs.append(pl.col(src_col).str.strip_chars().alias("semantic_tag"))
            elif dest_col == "PROPERTY_TYPE":
                exprs.append(pl.col(src_col).str.strip_chars().alias("property_type"))
        if object_label_col and object_label_col in df.columns:
            exprs.append(
                pl.col(object_label_col).str.strip_chars().alias("object_label")
            )
        df = df.select(exprs)
    else:
        # Standard ZOOMA format
        rename = {"PROPERTY_VALUE": "subject_label", "SEMANTIC_TAG": "semantic_tag"}
        cols = ["PROPERTY_VALUE", "SEMANTIC_TAG"]
        optional = {
            "PROPERTY_TYPE": "property_type",
            "STUDY": "zooma:study",
            "BIOENTITY": "zooma:bioentity",
            "ANNOTATOR": "author_label",
            "ANNOTATION_DATE": "mapping_date",
        }
        for src, dest in optional.items():
            if src in df.columns:
                rename[src] = dest
                cols.append(src)
        df = df.select([c for c in cols if c in df.columns]).rename(rename)
        df = df.with_columns([pl.col(c).str.strip_chars() for c in df.columns])

    # Filter rows missing required fields
    df = df.filter(
        (pl.col("subject_label") != "") & (pl.col("semantic_tag") != "")
    )

    # HTML-unescape subject_label and object_label where needed
    for col in ("subject_label", "object_label"):
        if col in df.columns and df[col].str.contains("&").any():
            unique = df[col].unique()
            mapping = pl.DataFrame({
                col: unique,
                f"_{col}": unique.map_elements(
                    html.unescape, return_dtype=pl.String
                ),
            })
            df = (
                df.join(mapping, on=col, how="left")
                .drop(col)
                .rename({f"_{col}": col})
            )

    # Convert mapping_date to ISO format (YYYY-MM-DD) if date_format given
    date_format = source_config.get("date_format")
    if date_format and "mapping_date" in df.columns:
        unique_dates = df["mapping_date"].unique()
        def _to_iso(val):
            if not val:
                return ""
            try:
                return datetime.strptime(val, date_format).strftime("%Y-%m-%d")
            except ValueError:
                return val  # keep raw if parse fails
        mapping = pl.DataFrame({
            "mapping_date": unique_dates,
            "_mapping_date": unique_dates.map_elements(
                _to_iso, return_dtype=pl.String
            ),
        })
        df = (
            df.join(mapping, on="mapping_date", how="left")
            .drop("mapping_date")
            .rename({"_mapping_date": "mapping_date"})
        )

    # Explode pipe-separated URIs into separate rows
    df = (
        df.with_columns(pl.col("semantic_tag").str.split("|"))
        .explode("semantic_tag")
        .with_columns(pl.col("semantic_tag").str.strip_chars())
        .filter(pl.col("semantic_tag") != "")
    )

    return df


# ---------------------------------------------------------------------------
# conversion to SSSOM
# ---------------------------------------------------------------------------

def convert_to_sssom(df):
    """Convert URIs to CURIEs.  All rows are preserved (no dedup)."""
    # Convert unique URIs only — massive speedup for large files
    unique_uris = df["semantic_tag"].unique().to_list()
    uri_curie_pairs: list[tuple[str, str]] = []
    unresolvable: set[str] = set()
    for uri in unique_uris:
        curie = uri_to_curie(uri)
        if curie is not None:
            uri_curie_pairs.append((uri, curie))
        else:
            unresolvable.add(uri)

    if unresolvable:
        row_count = df.filter(
            pl.col("semantic_tag").is_in(list(unresolvable))
        ).height
        examples = sorted(unresolvable)[:5]
        print(
            f"  Skipped {row_count} rows with unresolvable URIs, "
            f"e.g. {examples}",
            file=sys.stderr,
        )

    if not uri_curie_pairs:
        return pl.DataFrame()

    uri_map = pl.DataFrame({
        "semantic_tag": [p[0] for p in uri_curie_pairs],
        "object_id": [p[1] for p in uri_curie_pairs],
    })

    # Join replaces semantic_tag with object_id; drops unresolvable rows
    df = df.join(uri_map, on="semantic_tag", how="inner").drop("semantic_tag")

    # Add fixed columns
    df = df.with_columns(
        pl.lit("skos:closeMatch").alias("predicate_id"),
        pl.lit("semapv:ManualMappingCuration").alias("mapping_justification"),
    )

    # Rename property_type → subject_category
    if "property_type" in df.columns:
        df = df.rename({"property_type": "subject_category"})

    # Drop zooma extension columns and deduplicate
    for col in ("zooma:study", "zooma:bioentity"):
        if col in df.columns:
            df = df.drop(col)

    # Identify grouping vs aggregation columns
    agg_cols = [c for c in ("author_label", "mapping_date") if c in df.columns]
    group_cols = [c for c in df.columns if c not in agg_cols]

    if agg_cols:
        agg_exprs = []
        if "author_label" in df.columns:
            # Combine unique non-empty authors with semicolon delimiter
            agg_exprs.append(
                pl.col("author_label")
                .filter(pl.col("author_label") != "")
                .unique()
                .sort()
                .str.join(";")
                .alias("author_label")
            )
        if "mapping_date" in df.columns:
            # Keep the earliest non-empty date
            agg_exprs.append(
                pl.col("mapping_date")
                .filter(pl.col("mapping_date") != "")
                .min()
                .alias("mapping_date")
            )
        df = df.group_by(group_cols).agg(agg_exprs)
        # Fill nulls from aggregation (groups with all-empty values)
        df = df.fill_null("")
    else:
        df = df.unique()

    # Sort by all columns for stable diffs
    df = df.sort(df.columns)

    return df


def build_curie_map(df):
    """Build curie_map from prefixes actually used, omitting built-ins.

    Prefers the URI expansion observed in source data (recorded during
    URI-to-CURIE conversion) over bioregistry's canonical expansion.
    Also adds the zooma prefix when extension columns are present.
    """
    prefixes: set[str] = set(
        df["object_id"]
        .str.split_exact(":", 1)
        .struct.field("field_0")
        .unique()
        .to_list()
    )

    curie_map: dict[str, str] = {}
    for prefix in sorted(prefixes):
        if prefix.lower() in BUILTIN_PREFIXES:
            continue
        uri = _prefix_uri_expansions.get(prefix)
        if not uri:
            uri = get_uri_expansion(prefix)
        if uri:
            curie_map[prefix] = uri
        else:
            print(
                f"  Warning: no URI expansion for prefix '{prefix}'",
                file=sys.stderr,
            )

    return curie_map


# ---------------------------------------------------------------------------
# SSSOM/TSV writer
# ---------------------------------------------------------------------------

def _yaml_scalar(value):
    """Format a value for use as a plain YAML scalar, quoting if needed."""
    s = str(value)
    if not s:
        return '""'
    needs_quoting = False
    if s[0] in "{[\"'#&*!|>%@`,":
        needs_quoting = True
    elif ": " in s or " #" in s or s.endswith(":"):
        needs_quoting = True
    elif s.lower() in {"true", "false", "null", "yes", "no", "on", "off"}:
        needs_quoting = True
    if needs_quoting:
        return '"' + s.replace("\\", "\\\\").replace('"', '\\"') + '"'
    return s


def write_sssom_file(filepath, df, source_config, curie_map):
    """Write a single SSSOM/TSV file with embedded YAML metadata block."""
    # Column order: standard SSSOM slots, then extension columns
    standard_order = [
        "subject_label", "subject_category", "predicate_id",
        "object_id", "object_label", "mapping_justification",
        "author_label", "mapping_date",
    ]
    columns = [c for c in standard_order if c in df.columns]

    df_out = df.select(columns)
    name = filepath.name.replace(".sssom.tsv", "")

    with open(filepath, "w", newline="", encoding="utf-8") as fh:
        # --- metadata block ---
        if curie_map:
            fh.write("#curie_map:\n")
            for prefix in sorted(curie_map):
                fh.write(f"#  {prefix}: {_yaml_scalar(curie_map[prefix])}\n")

        fh.write(f"#mapping_set_id: {MAPPING_SET_ID_BASE}/{filepath.relative_to(MAPPINGS_DIR)}\n")

        desc = source_config.get("description", "")
        if desc:
            fh.write(f"#mapping_set_description: {_yaml_scalar(desc)}\n")

        fh.write("#subject_type: rdfs literal\n")

        subject_source = source_config.get("subject_source", "")
        if subject_source:
            fh.write(f"#subject_source: {_yaml_scalar(subject_source)}\n")

        # --- mappings block ---
        fh.write(
            df_out.write_csv(separator="\t", include_header=True, null_value="")
        )

    print(f"  Wrote {df_out.height} mappings → {filepath.name}")


# ---------------------------------------------------------------------------
# splitting (for clinvar-xrefs)
# ---------------------------------------------------------------------------

def split_by_hash(df, split_count=32):
    """Partition rows into *split_count* buckets by hash of (label, id).

    Uses SHA-1 for cross-version determinism. Each bucket is independently sorted.
    """
    def _sha1_bucket(label, oid):
        key = f"{label}\t{oid}".encode()
        return int.from_bytes(hashlib.sha1(key).digest()[:8], "little") % split_count

    df = df.with_columns(
        pl.struct("subject_label", "object_id")
        .map_elements(
            lambda row: _sha1_bucket(row["subject_label"], row["object_id"]),
            return_dtype=pl.Int32,
        )
        .alias("_bucket")
    )
    sort_cols = ["subject_label", "object_id"]
    for ext in ("author_label", "mapping_date"):
        if ext in df.columns:
            sort_cols.append(ext)

    buckets = []
    for i in range(split_count):
        bucket = df.filter(pl.col("_bucket") == i).drop("_bucket")
        if bucket.height > 0:
            bucket = bucket.sort(sort_cols)
            buckets.append((i, bucket))
    return buckets


# ---------------------------------------------------------------------------
# per-source processing
# ---------------------------------------------------------------------------

def process_source(name, source_config, data_path):
    """Parse → convert → write SSSOM file(s) for one source."""
    _curie_cache.clear()
    _prefix_uri_expansions.clear()
    print(f"Processing: {name}")

    df = parse_source(data_path, source_config)
    print(f"  Parsed {df.height} rows")

    df = convert_to_sssom(df)
    print(f"  {df.height} mappings")

    if df.height == 0:
        print(f"  Warning: no mappings for {name}", file=sys.stderr)
        return

    should_split = source_config.get("split", False)
    split_count = source_config.get("split_count", 16)

    if should_split:
        split_dir = MAPPINGS_DIR / name
        split_dir.mkdir(parents=True, exist_ok=True)
        buckets = split_by_hash(df, split_count)
        for i, bucket in buckets:
            cm = build_curie_map(bucket)
            fp = split_dir / f"{name}-{i + 1:02d}.sssom.tsv"
            write_sssom_file(fp, bucket, source_config, cm)
    else:
        cm = build_curie_map(df)
        fp = MAPPINGS_DIR / f"{name}.sssom.tsv"
        write_sssom_file(fp, df, source_config, cm)


# ---------------------------------------------------------------------------
# main
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description="Build SSSOM/TSV files from zooma source data"
    )
    parser.add_argument(
        "--local-data",
        type=Path,
        help="Read source files from this directory instead of downloading",
    )
    args = parser.parse_args()

    config = load_sources()
    sources = config["sources"]

    MAPPINGS_DIR.mkdir(parents=True, exist_ok=True)

    if args.local_data:
        for name, src in sources.items():
            filename = src["import_url"].rsplit("/", 1)[-1].split("?")[0]
            data_path = args.local_data / filename
            if not data_path.exists():
                print(f"Skipping {name}: {data_path} not found")
                continue
            process_source(name, src, data_path)
    else:
        with tempfile.TemporaryDirectory() as tmpdir:
            for name, src in sources.items():
                try:
                    data_path = download_file(src["import_url"], tmpdir)
                    process_source(name, src, data_path)
                except Exception as exc:
                    print(f"  ERROR processing {name}: {exc}", file=sys.stderr)

    write_registry_manifest(sources)


def write_registry_manifest(sources):
    """Generate mappings.yml from the actual output files."""
    refs = []
    for name, src in sources.items():
        if src.get("split"):
            split_dir = MAPPINGS_DIR / name
            if split_dir.is_dir():
                for fp in sorted(split_dir.glob("*.sssom.tsv")):
                    local = f"mappings/{name}/{fp.name}"
                    refs.append({
                        "mapping_set_id": f"{MAPPING_SET_ID_BASE}/{name}/{fp.name}",
                        "local_name": local,
                        "mapping_set_group": name,
                    })
        else:
            fp = MAPPINGS_DIR / f"{name}.sssom.tsv"
            if fp.exists():
                refs.append({
                    "mapping_set_id": f"{MAPPING_SET_ID_BASE}/{fp.name}",
                    "local_name": f"mappings/{fp.name}",
                    "mapping_set_group": name,
                })

    manifest = {
        "mapping_registry_id": "https://w3id.org/sssom/commons/ebi-text-mappings",
        "mapping_registry_title": "EBI Text Mappings",
        "mapping_registry_description": "Curated text-to-ontology mappings from EBI datasources",
        "homepage": "https://github.com/mapping-commons/ebi-text-mappings",
        "issue_tracker": "https://github.com/mapping-commons/ebi-text-mappings/issues",
        "mapping_set_references": refs,
    }

    manifest_path = ROOT_DIR / "mappings.yml"
    with open(manifest_path, "w", encoding="utf-8") as fh:
        yaml.dump(manifest, fh, default_flow_style=False, sort_keys=False)
    print(f"Wrote registry manifest → {manifest_path.name} ({len(refs)} mapping sets)")


if __name__ == "__main__":
    main()
