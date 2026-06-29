"""ddbj_gff: INSDC/SO GFF3 parser and object model (phase 1)."""

from .model import Directive, Feature, GffDocument, Span
from .parser import parse

__all__ = ["parse", "GffDocument", "Feature", "Span", "Directive"]
