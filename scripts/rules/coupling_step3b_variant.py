"""Snakemake rule `coupling_step3b_arv`: Step 3b event data for an arousal-exclusion sensitivity variant.

Runs scripts/rules/coupling_step3b.py unchanged, except that the snore-burst references are read from the
variant's condition labels instead of the primary ones (burst_references defaults to 'primary'). A wrapper keeps the
primary rule's code untouched, so the primary results are not invalidated.
"""
import functools
import runpy

import snoredreamt.snoring as snoring

snoring.burst_references = functools.partial(snoring.burst_references,
                                             labels_variant=snakemake.wildcards.arv)  # noqa: F821
runpy.run_path('scripts/rules/coupling_step3b.py', init_globals={'snakemake': snakemake})  # noqa: F821
