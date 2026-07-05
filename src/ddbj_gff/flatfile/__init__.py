"""flatfile2gff: DDBJ flatfile -> canonical INSDC GFF3."""
from .molecule import MoleculeInfo, detect_molecule
from .convert import flatfile_to_gff

__all__ = ["MoleculeInfo", "detect_molecule", "flatfile_to_gff"]
