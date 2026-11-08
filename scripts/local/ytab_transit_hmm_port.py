#!/usr/bin/env python3
"""Baseline: TRANSIT's HMM (pytransit/analysis/hmm.py, commit 50e23d6) ported as written, run on Hermes data.

The TA-site universe is replaced by the union of insertion sites observed in any of the given libraries
(zeros = sites seen only in other libraries). Each library is TTR-normalized (YTAB ttr), replicates are
averaged (TRANSIT default "Mean"), +1 for the shifted geometric. Emission means, transitions (from pins),
initial probabilities, Viterbi decoding and the gene rule (ES only if every site is ES, else the most
frequent state) follow the TRANSIT source. No parameter is fitted.
Usage: ytab_transit_hmm_port.py OUTPUT [--libraries SRR...]
"""
import argparse
import json
import math
from pathlib import Path
import sys

import numpy as np
import scipy.stats

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT/"src")); sys.path.append(str(ROOT/"scripts/local"))
from ytab.essentiality import control_conditioned_gumbel as ccg
from ytab.normalization.ttr import trim_mean, write_table
import ytab_hermes_hmm as runner

LABELS = {0:"ES", 1:"GD", 2:"NE", 3:"GA"}


def calculate_pins(reads):
    non_ess, temp = [], []
    for rd in reads:
        if rd >= 1:
            if len(temp) < 10:non_ess.extend(temp)
            non_ess.append(rd); temp = []
        else:temp.append(rd)
    return sum(1 for rd in non_ess if rd >= 1)/float(len(non_ess))


def transit_parameters(O):
    reads = O-1; nz = np.sort(reads[reads != 0]); mean_r = np.average(nz[:int(.95*len(nz))])
    mu = np.array([1/.99, .01*mean_r+2, mean_r, mean_r*5.]); L = 1./mu
    pins = calculate_pins(reads); pnon = 1.-pins
    for r in range(100):
        if pnon**r < .01:break
    B1 = scipy.stats.geom(L[2]).pmf(1)
    a = math.log1p(-B1**r); b = r*math.log(B1)+math.log(1./3)
    A = np.full((4, 4), b); np.fill_diagonal(A, a)
    PI = np.zeros(4); PI[0] = .7; PI[1:] = .3/3
    return dict(mu=mu, L=L, pins=pins, r=r, logA=A, PI=PI)


def viterbi(logA, L, PI, O):
    logB = np.column_stack([scipy.stats.geom(l).logpmf(O) for l in L])
    delta = np.log(PI)+logB[0]; back = np.zeros((len(O), 4), dtype=np.int8)
    for t in range(1, len(O)):
        nus = delta[None, :]+logA; back[t] = nus.argmax(1); delta = nus.max(1)+logB[t]
    path = np.zeros(len(O), dtype=int); path[-1] = int(np.argmax(delta))
    for t in range(len(O)-1, 0, -1):path[t-1] = back[t, path[t]]
    return path


def main():
    p = argparse.ArgumentParser(description=__doc__); p.add_argument("output", type=Path); p.add_argument("--libraries", nargs="+", default=runner.GALE)
    a = p.parse_args()
    libs = []
    for s in a.libraries:
        path = ROOT/"output/projects/Gale_2020/create_hit_file"/s/f"{s}_hits.txt"
        chrom = np.loadtxt(path, delimiter="\t", skiprows=1, usecols=(0,), dtype=str, ndmin=1)
        pos, count = np.loadtxt(path, delimiter="\t", skiprows=1, usecols=(11, 12), dtype=np.int64, ndmin=2).T
        d = {}
        for c, q, n in zip(chrom, pos, count):d[(c, int(q))] = d.get((c, int(q)), 0)+int(n)
        libs.append(d)
    sites = sorted(set().union(*libs)); raw = np.array([[lib.get(x, 0) for lib in libs] for x in sites], dtype=float)
    norm = np.zeros_like(raw)
    for j in range(raw.shape[1]):
        pos = raw[:, j][raw[:, j] > 0]; theta = len(pos)/len(raw); norm[:, j] = raw[:, j]*100./(theta*trim_mean(pos, .05))
    O = np.round(norm.mean(1))+1   # TRANSIT combine_replicates("Mean") rounds
    params = transit_parameters(O)
    states = np.zeros(len(sites), dtype=int); chroms = np.array([c for c, _ in sites])
    for c in dict.fromkeys(chroms):   # TRANSIT walks one genome; here each chromosome separately
        m = chroms == c; states[m] = viterbi(params["logA"], params["L"], params["PI"], O[m])
    genes = ccg.annotations(runner.GFF)
    idx = ccg.gene_sites(genes, sites, ccg.Settings(nterm=0., cterm=0., ignore_stop=True))
    rows = []
    for g, i in zip(genes, idx):
        st = states[i]; n = len(i); counts = [int((st == k).sum()) for k in range(4)]
        call = "N/A" if n == 0 else "ES" if counts[0] == n else LABELS[max((counts[k], k) for k in range(4))[1]]
        rows.append(dict(ORF=g["id"], Name=g["name"], sites=n, ES=counts[0], GD=counts[1], NE=counts[2], GA=counts[3],
                         saturation=float((raw[i].sum(1) > 0).mean()) if n else None, call=call))
    a.output.mkdir(parents=True, exist_ok=True); write_table(a.output/"genes.tsv", list(rows[0]), rows)
    summary = dict(samples=a.libraries, sites=len(sites), parameters={k:(v.tolist() if hasattr(v, "tolist") else v) for k, v in params.items()},
                   site_states={LABELS[k]:float((states == k).mean()) for k in range(4)},
                   gene_calls={s:sum(r["call"] == s for r in rows) for s in ["ES", "GD", "NE", "GA", "N/A"]})
    (a.output/"summary.json").write_text(json.dumps(summary, indent=2)+"\n"); print(json.dumps(summary, indent=1))


if __name__ == "__main__":main()
