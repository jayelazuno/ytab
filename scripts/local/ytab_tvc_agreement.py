#!/usr/bin/env python3
"""Agreement between treated-vs-control methods (ytab.fitness.agreement): Spearman of effects and signed significance,
and Jaccard / overlap coefficient of the calls per direction, per background x treated condition. Give two or more of
--zscore (a treated_vs_parent result folder), --resampling and --gumbel (ytab treated-vs-control run folders), and
optionally --published (TSV: ORF, call, [effect, score, background, condition]).

  python scripts/local/ytab.py tvc-agreement --zscore output/projects/P/treated_vs_parent/RUN \\
      --resampling output/projects/P/treated_vs_control/RUN --gumbel output/projects/P/treated_vs_control/RUN \\
      --output output/projects/P/treated_vs_control/agreement/NAME
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2]/"src"))
from ytab.fitness import agreement  # noqa: E402


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    for m in ("zscore", "resampling", "gumbel"):p.add_argument(f"--{m}", type=Path)
    p.add_argument("--published", type=Path); p.add_argument("--output", required=True, type=Path)
    a = p.parse_args()
    inputs = {m: getattr(a, m) for m in ("zscore", "resampling", "gumbel") if getattr(a, m)}
    print(json.dumps(agreement.run(inputs, a.output, a.published), indent=1))


if __name__ == "__main__":main()
