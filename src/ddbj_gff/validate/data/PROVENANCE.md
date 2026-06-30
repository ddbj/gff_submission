# Bundled controlled-vocabulary data — provenance

- `feature-mapping.tsv`: snapshot of `src/main/resources/feature-mapping.tsv` from
  https://github.com/enasequence/gff3tools (branch `main`). Columns:
  SOID, SO term, Definition, Feature (INSDC), Qualifier 1, Qualifier 2 (tab-separated, has header).
  License: Apache-2.0 (see NOTICE). Refresh: re-run the curl in the implementation plan Task 1.
- `dbxref.tsv`: curated subset of the INSDC db_xref vocabulary
  (https://www.insdc.org/submitting-standards/dbxref-qualifier-vocabulary/). One DBTAG per line.
