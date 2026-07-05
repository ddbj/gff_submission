"""flatfile2gff: DDBJ flatfile -> canonical INSDC GFF3."""
from .molecule import MoleculeInfo, detect_molecule

__all__ = ["MoleculeInfo", "detect_molecule"]
