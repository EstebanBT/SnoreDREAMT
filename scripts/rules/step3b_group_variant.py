"""Snakemake rule `step3b_group_arv`: Step 3b/3c group analysis for an arousal-exclusion sensitivity variant.

Runs scripts/rules/step3b_group.py unchanged except for one line: the per-family primary event sets are not
re-selected from the clean data but fixed to those of the primary analysis (input `fixed_primaries`), so that the
sensitivity analysis compares like with like. A wrapper keeps the primary rule's code untouched.
"""
import pandas as pd

PATH = 'scripts/rules/step3b_group.py'
LINE = "PRIMARIES = {fam: tuple(g.iloc[0][['event_set', 'variant']]) for fam, g in cover.groupby('family')}"
cs = pd.read_csv(snakemake.input.fixed_primaries)  # noqa: F821
fixed = {r.family: (r.event_set, r.variant) for r in cs[cs.primary].itertuples()}
src = open(PATH).read()
assert src.count(LINE) == 1, 'step3b_group.py changed: update the wrapper'
src = src.replace(LINE, "PRIMARIES = FIXED_PRIMARIES\nprint('selection by clean z would have been:', "
                        "{fam: tuple(g.iloc[0][['event_set', 'variant']]) for fam, g in cover.groupby('family')})")
exec(compile(src, PATH, 'exec'), {'snakemake': snakemake, 'FIXED_PRIMARIES': fixed,  # noqa: F821
                                  '__name__': '__main__', '__file__': PATH})
