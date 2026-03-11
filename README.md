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
| `sysmicro.sssom.tsv` | [SysMicro / CMPO](https://www.ebi.ac.uk/fg/sym) |
| `atlas.sssom.tsv` | [Expression Atlas](https://www.ebi.ac.uk/gxa) |
| `ebisc.sssom.tsv` | [EBiSC](https://cells.ebisc.org/) |
| `cttv.sssom.tsv` | [Open Targets](https://www.targetvalidation.org) |
| `uniprot.sssom.tsv` | [UniProt](https://www.ebi.ac.uk/uniprot) |
| `eva-clinvar.sssom.tsv` | [EVA / ClinVar](https://www.ebi.ac.uk/eva) |
| `cbi.sssom.tsv` | [BioSamples (plants)](https://www.ebi.ac.uk/biosamples) |
| `clinvar-xrefs-01…16.sssom.tsv` | [ClinVar xrefs](https://www.ncbi.nlm.nih.gov/clinvar) |
| `metabolights.sssom.tsv` | [MetaboLights](https://www.ebi.ac.uk/metabolights) |
| `ukbiobank.sssom.tsv` | [UK Biobank / EFO](https://github.com/EBISPOT/EFO-UKB-mappings) |
| `ebi-biosamples.sssom.tsv` | [BioSamples](https://www.ebi.ac.uk/biosamples/) |
| `faang.sssom.tsv` | [FAANG](https://www.ebi.ac.uk/vg/faang) |
| `hca.sssom.tsv` | [Human Cell Atlas](https://www.ebi.ac.uk/about/collaborations/human-cell-atlas) |
| `gwas.sssom.tsv` | [GWAS Catalog](http://www.ebi.ac.uk/gwas) |

ClinVar cross-references are split into 16 files by hash to keep file sizes
manageable.

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
