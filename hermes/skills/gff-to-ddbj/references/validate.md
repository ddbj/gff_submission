# ddbj_gff.validate (reference)

`$PY -m ddbj_gff.validate --gff GFF [--severity CODE=LEVEL]`

Detect-only INSDC profile validator (does not modify the file). Emits tab-separated
`SEVERITY<TAB>code<TAB>message` diagnostic lines, where SEVERITY is `ERROR` / `WARNING` /
`INFO`, followed by a final `summary: ERROR=N` line (and `WARNING=M` when present). The
command exits non-zero (1) iff any ERROR. Fix inputs until there are no `ERROR` lines / exit
0 (`summary: ERROR=0`); warnings may legitimately remain. `--severity CODE=off|error|warning|info`
overrides a rule's level. Unknown/extra GFF attributes are left untouched (not flagged). This
detect-only check should be run after `normalize` and again after `repair`.
