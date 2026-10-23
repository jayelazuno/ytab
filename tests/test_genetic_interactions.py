import numpy as np
import pytest
from scipy import sparse

from ytab.normalization.methods import normalize

from ytab.essentiality import genetic_interactions as gi


def design(rng, genes=400, sites=20, effect_a=None, effect_b=None, pools=False, phi=.2):
    """Strains A and B, each one pool with 2 control + 2 treated cultures (or 2 pools with 1 + 1)."""
    effect_a = np.ones(genes) if effect_a is None else effect_a; effect_b = np.ones(genes) if effect_b is None else effect_b
    blocks, strains, conds, libs, idx, offset = [], [], [], [], [[] for _ in range(genes)], 0
    for strain, eff in (("A", effect_a), ("B", effect_b)):
        for pool in ((1, 2) if pools else (1,)):
            clone = np.exp(np.log(30)+1.5*rng.standard_normal(genes*sites)); g = np.repeat(np.arange(genes), sites)
            draw = lambda mu:rng.negative_binomial(1/phi, (1/phi)/((1/phi)+mu))
            cols = ([draw(clone), draw(clone*eff[g])] if pools else [draw(clone), draw(clone*1.2), draw(clone*eff[g]), draw(clone*eff[g]*.9)])
            block = np.zeros((genes*sites, 0)); blocks.append(np.column_stack(cols))
            n = len(cols); strains += [strain]*n; libs += [f"{strain}{pool}"]*n
            conds += ["control", "treated"] if pools else ["control", "control", "treated", "treated"]
            for k in range(genes):idx[k].append(offset+np.arange(k*sites, (k+1)*sites))
            offset += genes*sites
    raw = np.zeros((offset, len(strains)))
    col = 0; row = 0
    for b in blocks:raw[row:row+len(b), col:col+b.shape[1]] = b; row += len(b); col += b.shape[1]
    genes_meta = [dict(id=f"g{k}", name=f"g{k}") for k in range(genes)]
    raw = raw.astype(np.int64); values = normalize(sparse.csc_matrix(raw), "betageom", libs)[0].toarray()
    return raw, values, strains, conds, libs, genes_meta, [np.concatenate(x) for x in idx]


def test_classes_follow_transit():
    assert gi.classify(-1, -2, 0) == "Aggravating" and gi.classify(1, -.5, -1.5) == "Alleviating" and gi.classify(1, 2, 1) == "Suppressive"


def test_raw_depth_offsets_do_not_bias_calibration():
    rng = np.random.default_rng(2)
    raw, _, s, c, l, genes, idx = design(rng)                  # unnormalized: replicate cultures differ in depth
    rows, _ = gi.analyze(raw.astype(float), s, c, l, genes, idx)
    p = np.array([r["p_value"] for r in rows if r["status"] == "tested"])
    assert .02 <= np.mean(p <= .05) <= .09


def test_pairing_rules():
    with pytest.raises(ValueError, match="equal numbers"):gi.pairs_by_library(["A"]*3, ["control", "control", "treated"], ["L"]*3)
    with pytest.raises(ValueError, match="mixes strains"):gi.pairs_by_library(["A", "B"], ["control", "treated"], ["L", "L"])
    assert gi.pairs_by_library(["A"]*4, ["control", "treated", "control", "treated"], ["L"]*4) == {"A":[("L", 0, 1), ("L", 2, 3)]}


@pytest.mark.parametrize("pools", [False, True])
def test_null_calibrated_and_shared_response_not_called(pools):
    rng = np.random.default_rng(1 if pools else 2)
    shared = np.ones(400); shared[:40] = .25                     # same treatment response in both strains: no GI
    raw, v, s, c, l, genes, idx = design(rng, effect_a=shared, effect_b=shared, pools=pools)
    rows, d = gi.analyze(v, s, c, l, genes, idx, raw=raw)
    p = np.array([r["p_value"] for r in rows if r["status"] == "tested"])
    assert .02 <= np.mean(p <= .05) <= .09 and d["calls"]["No interaction"] >= 395


def test_detects_query_specific_response():
    rng = np.random.default_rng(3)
    eff_b = np.ones(400); eff_b[:40] = .1
    raw, v, s, c, l, genes, idx = design(rng, effect_b=eff_b)
    rows, d = gi.analyze(v, s, c, l, genes, idx)
    called = {r["ORF"] for r in rows if r["call"] in gi.CLASSES}
    assert len(called & {f"g{k}" for k in range(40)}) >= 30 and len(called - {f"g{k}" for k in range(40)}) <= 3
    assert all(r["call"] == "Aggravating" for r in rows if r["ORF"] in called and int(r["ORF"][1:]) < 40)


def test_cli_end_to_end(tmp_path):
    import json, subprocess, sys
    from pathlib import Path
    from test_multicondition import normalization_output
    lost = lambda i:np.where((i >= 60) & (i < 90), 0., 1.)            # gene g02 lost under treatment, query strain only
    one = lambda i:np.ones(len(i))
    names = [("wc1", "W1", one), ("wt1", "W1", one), ("wc2", "W2", one), ("wt2", "W2", one),
             ("mc1", "M1", one), ("mt1", "M1", lost), ("mc2", "M2", one), ("mt2", "M2", lost)]
    meta = [dict(observation_id=f"P::{x}", display_name=x, project="P", library_id=f"P::{l}", pool=l, role="any") for x, l, _ in names]
    gff, reference = normalization_output(tmp_path, [f for *_, f in names], meta)
    strain = lambda x:"WT" if x[0] == "w" else "mut"; cond = lambda x:"control" if x[1] == "c" else "treated"
    args = [sys.executable, str(Path(__file__).resolve().parents[1]/"scripts/local/ytab_genetic_interactions.py"), "--normalized", str(tmp_path/"norm"),
            "--gff", str(gff), "--reference", str(reference), "--output", str(tmp_path/"out"), "--reference-strain", "WT", "--query-strain", "mut"]
    args += [x for n, *_ in names for x in ("--sample", f"P::{n}={strain(n)}:{cond(n)}")]
    done = subprocess.run(args, capture_output=True, text=True)
    assert done.returncode == 0, done.stderr
    out = json.loads(done.stdout); assert out["reference"] == "WT" and out["query"] == "mut"
    import csv
    rows = {r["ORF"]:r for r in csv.DictReader((tmp_path/"out"/"genes.tsv").open(), delimiter="\t")}
    assert float(rows["g02"]["delta_log2FC"]) < -1 and rows["g02"]["interaction"] == "Aggravating"
    bad = subprocess.run(args[:-2]+["--sample", "P::mt2=mut:later"], capture_output=True, text=True)
    assert bad.returncode != 0 and "control|treated" in bad.stderr


def test_holm_matches_definition():
    p = np.array([.01, .04, .03, .005])
    assert np.allclose(gi.holm(p), [.03, .06, .06, .02])


def test_treat_reduces_to_t_test_and_is_conservative_for_rope():
    from scipy import stats
    d, se, df = np.array([.1, -.8, 1.5]), np.array([.2, .3, .4]), np.array([10., 20., 30.])
    assert np.allclose(gi.treat_p(d, se, df, 0.), 2*stats.t.sf(np.abs(d)/se, df))
    assert np.all(gi.treat_p(d, se, df, .5) >= gi.treat_p(d, se, df, 0.))
    assert gi.treat_p(np.array([.5]), np.array([1e-6]), np.array([50.]), .5)[0] > .49   # at the ROPE boundary p ~ 0.5


def test_interval_options_and_holm_calls():
    rng = np.random.default_rng(3)
    eff_b = np.ones(400); eff_b[:40] = .1
    raw, v, s, c, l, genes, idx = design(rng, effect_b=eff_b)
    base, _ = gi.analyze(v, s, c, l, genes, idx)
    rope, _ = gi.analyze(v, s, c, l, genes, idx, settings=gi.Settings(rope=.5))
    fwer, _ = gi.analyze(v, s, c, l, genes, idx, settings=gi.Settings(correction="holm"))
    tested = [r for r in base if r["status"] == "tested"]
    assert all(r["delta_lower"] <= r["delta_log2FC"] <= r["delta_upper"] for r in tested)
    assert all(r["p_value"] >= b["p_value"] for r, b in zip(rope, base) if r["status"] == "tested")
    calls = lambda rows:{r["ORF"] for r in rows if r["call"] in gi.CLASSES}
    assert calls(fwer) <= calls(base) and calls(rope) <= calls(base) and len(calls(fwer)) > 0
    assert all(abs(r["delta_log2FC"]) > .5 for r in rope if r["call"] in gi.CLASSES)
    with pytest.raises(ValueError, match="cannot be combined"):gi.Settings(rope=.5, empirical_null=True).validate()
