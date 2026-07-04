from __future__ import annotations

import copy

from ..validate.vocab import load_vocab
from .config import NormalizeConfig
from .passes import (NormalizeContext, pass_directives, pass_coerce_transcript_to_mrna,
                     pass_wrap_cds_in_mrna, pass_circular_origin, pass_so_terms,
                     pass_transl_except, pass_anticodon)
from .report import NormalizationReport

ALL_PASSES = [pass_directives, pass_coerce_transcript_to_mrna, pass_wrap_cds_in_mrna,
              pass_circular_origin, pass_so_terms, pass_transl_except, pass_anticodon]

# actions that represent a clean applied change; everything else needs human attention
_APPLIED = {"add-directive", "rename-type", "add-qualifier", "add-child-feature"}


def normalize(doc, *, seq_lengths=None, config=None) -> tuple:
    config = config or NormalizeConfig()
    work = copy.deepcopy(doc)
    ctx = NormalizeContext(vocab=load_vocab(), seq_lengths=seq_lengths, config=config)
    applied: list = []
    unresolved: list = []
    for run_pass in ALL_PASSES:
        for change in run_pass(work, ctx):
            (applied if change.action in _APPLIED else unresolved).append(change)
    return work, NormalizationReport(applied=applied, unresolved=unresolved)
