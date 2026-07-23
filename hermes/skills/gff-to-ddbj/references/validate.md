# ddbj_gff.validate (reference)

`$PY -m ddbj_gff.validate --gff GFF [--severity CODE=LEVEL]`

detect-only INSDC profile validator (does not modify the file). Emits diagnostics with a
severity (`ERR:` / `WAR:` / `INFO:`) and a code. Fix inputs until there are no `ERR:`
lines; warnings may legitimately remain. `--severity CODE=off|error|warning|info` overrides
a rule's level. Unknown/extra GFF attributes are left untouched (not flagged). Run after
`normalize` and again after `repair`.
