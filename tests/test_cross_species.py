import sys
from pathlib import Path

import numpy as np

sys.path.append(str(Path(__file__).resolve().parents[1]/"scripts/local"))
import ytab_cross_species as cs  # noqa: E402


def test_orthologs_with_scer_id_header(tmp_path):
    path = tmp_path/"o.tsv"
    path.write_text("cglab_id\tname\tscer_id\nG1\tA\tYAL001C\nG2\tB\tYBR265W\nG3\tC\tYBR265W\nG4\tD\t\n")
    assert cs.read_orthologs(path) == {"G1":"YAL001C"}   # YBR265W has two genes: not one-to-one


def test_orthologs_without_header(tmp_path):
    path = tmp_path/"o.txt"
    path.write_text("# comment\nC1_001W\tx\ty\tYAL001C\tz\nSPAC1.01\tYBR265W(N)+YDR302W(C)\nSPAC1.02\tYCL030C\nSPAC1.03\tYCL030C\nSPAC2.01\tQ0010\n")
    assert cs.read_orthologs(path) == {"C1_001W":"YAL001C", "SPAC2.01":"Q0010"}


def test_kappa():
    a = np.array([1, 1, 0, 0]); assert cs.kappa(a, a) == 1
    assert abs(cs.kappa(a, np.array([1, 0, 1, 0]))) < 1e-12
