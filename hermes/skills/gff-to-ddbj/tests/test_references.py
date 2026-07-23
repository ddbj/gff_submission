from pathlib import Path

REF = Path(__file__).resolve().parents[1] / "references"

# each reference file must exist and mention these key tokens (exact flags/ops)
REQUIRED = {
    "normalize.md": ["ddbj_gff.normalize", "--gff", "--out", "--report", "--transl-table"],
    "validate.md":  ["ddbj_gff.validate", "detect-only", "ERR:"],
    "repair.md":    ["ddbj_gff.repair", "--detect", "--apply", "internal-stop-to-misc",
                     "utr-absent-to-partial-mrna", "missing-start-stop-to-partial-cds"],
    "gff2mss.md":   ["gff2mss", "--mss-config", "--common", "--sequence-roles", ".ann"],
    "validator.md": ["ddbj-validator", "-f", "-j 1", "fixed/", "background"],
}


def test_reference_pages_exist_and_document_flags():
    for fname, tokens in REQUIRED.items():
        text = (REF / fname).read_text()
        for tok in tokens:
            assert tok in text, f"{fname} missing {tok!r}"
