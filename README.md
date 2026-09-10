# EBI Text Mappings

Curated text-to-ontology mappings from EBI datasources and ClinVar, published as
[SSSOM](https://mapping-commons.github.io/sssom/) TSV files.

These are **literal mappings** (`subject_type: rdfs literal`) — each row maps a
free-text string (`subject_label`) to an ontology term (`object_id`) via
`skos:closeMatch`.  CURIEs are normalised using
[bioregistry](https://bioregistry.io/).

## Sources

Mappings are drawn from the following EBI data sources (configured in
[sources.yml](sources.yml)):

| SSSOM file | Source |
|---|---|
| `sysmicro.sssom.tsv` | [SysMicro / CMPO](https://www.ebi.ac.uk/fg/sym) (static copy in `static/`, the CMPO site is offline) |
| `atlas.sssom.tsv` | [Expression Atlas](https://www.ebi.ac.uk/gxa) |
| `ebisc.sssom.tsv` | [EBiSC](https://cells.ebisc.org/) |
| `cttv.sssom.tsv` | [Open Targets](https://platform.opentargets.org) disease string curation ([opentargets/curation](https://github.com/opentargets/curation/tree/master/mappings/disease)) |
| `opentargets-cancer.sssom.tsv` | [Open Targets](https://platform.opentargets.org) cancer type to EFO mappings |
| `bmb-wp7.sssom.tsv` | [BioMedBridges WP7](https://www.biomedbridges.eu) disease/phenotype mappings (2014, from the Zooma FTP) |
| `uniprot.sssom.tsv` | [UniProt](https://www.ebi.ac.uk/uniprot) |
| `eva-clinvar.sssom.tsv` | [EVA / ClinVar](https://www.ebi.ac.uk/eva) |
| `cbi.sssom.tsv` | [BioSamples (plants)](https://www.ebi.ac.uk/biosamples) |
| `clinvar-xrefs.sssom.tsv` | [ClinVar xrefs](https://www.ncbi.nlm.nih.gov/clinvar) |
| `metabolights.sssom.tsv` | [MetaboLights](https://www.ebi.ac.uk/metabolights) |
| `ukbiobank.sssom.tsv` | [UK Biobank / EFO](https://github.com/EBISPOT/EFO-UKB-mappings) |
| `ebi-biosamples.sssom.tsv` | [BioSamples](https://www.ebi.ac.uk/biosamples/) |
| `faang.sssom.tsv` | [FAANG](https://www.ebi.ac.uk/vg/faang) |
| `hca.sssom.tsv` | [Human Cell Atlas](https://www.ebi.ac.uk/about/collaborations/human-cell-atlas) |
| `gwas.sssom.tsv` | [GWAS Catalog](http://www.ebi.ac.uk/gwas) |
| `pgs.sssom.tsv` | [PGS Catalog](https://www.pgscatalog.org) reported trait to EFO mappings |

Most sources are in Zooma's 7-column format.  GWAS, PGS and the Open Targets
cancer file use their own columns, declared via `column_map` in
[sources.yml](sources.yml); `SEMANTIC_TAG` values may be full URIs or compact
ids such as `MONDO_0004989`, and pipe-separated multi-term values are exploded
into one mapping per term.

## Building locally

```bash
python -m venv .venv && source .venv/bin/activate
pip install -r requirements.txt

# download all sources and build SSSOM files
python scripts/build_sssom.py

# or use local data files
python scripts/build_sssom.py --local-data data
```

Output is written to `mappings/`.

## CI

A GitHub Actions workflow runs daily at 03:00 UTC, downloads fresh source data,
rebuilds SSSOM files, and commits any changes.  See
[.github/workflows/update-mappings.yml](.github/workflows/update-mappings.yml).

If a source fails to download or parse, its previously committed SSSOM file is
kept as-is, the other sources are still updated and committed, and the run is
marked as failed so the problem is visible.
