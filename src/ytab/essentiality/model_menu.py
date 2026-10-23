"""Side-by-side essentiality model menu: one row per gene with every available model's score and call.

Models (see output/essentiality/hermes_hmm/hermes_hmm_methods.md):
- hmm          Hermes window HMM v2 (P_ES_run, call). Unsupervised: needs no labels. Always present.
- hybrid       random forest on SummaryTable gene features + HMM features, trained on labelled genes.
               Present only when labels and gene features (a SummaryTable feature table, or the legacy
               classifier's prediction table, which carries the same columns) are supplied. Labelled training genes get out-of-fold
               scores (folds grouped by chromosome), so no gene is scored by a model that saw its label;
               all other genes are scored by the model trained on every training gene.
- rf_existing  the legacy classifier's score ("RF - G4", trained on S. cerevisiae), for comparison only.
               Present when its prediction table is supplied.
Agreement of the HMM and hybrid ES calls gives a consensus class; nothing here is tuned to labels beyond
the hybrid's own training, whose settings are those of the pre-registered hybrid benchmark.
"""
import numpy as np
import pandas as pd
from sklearn.ensemble import RandomForestClassifier

REVISION = "model_menu_v1"
RF_FEATURES = ["Length", "Hits", "Reads", "Longest free interval", "Freedom index", "Neighborhood index", "Upstream hits 100"]
HMM_FEATURES = ["P_ES_run", "P_ES", "P_GD", "P_NE", "P_GA", "mappable_fraction"]
SEED = 90210
FOLDS = 5
THRESHOLD = .5
AGREEMENT = {(True, True):"consensus essential", (True, False):"HMM only", (False, True):"hybrid only", (False, False):"not essential"}


# SummaryTable feature-table columns -> the gene-feature names the hybrid was benchmarked with
SUMMARY_COLUMNS = {"standard_name":"ORF", "coding_length":"Length", "hits":"Hits", "reads":"Reads", "max_free_region":"Longest free interval",
                   "freedom_index":"Freedom index", "neighborhood_index":"Neighborhood index", "upstream_hits_100":"Upstream hits 100"}


ROUNDING = {"Freedom index":2, "Neighborhood index":3}


def summary_features(table: pd.DataFrame) -> pd.DataFrame:
    """Gene features from a SummaryTable feature table (no classifier needed)."""
    missing = set(SUMMARY_COLUMNS) - set(table.columns)
    if missing:raise ValueError("Feature table is missing columns: " + ", ".join(sorted(missing)))
    out = table[list(SUMMARY_COLUMNS)].rename(columns=SUMMARY_COLUMNS)
    # The validated hybrid was trained on the classifier table's printed precision (%.2f and %.3f); keep it.
    for column, digits in ROUNDING.items():
        out[column] = [float(f"%.{digits}f" % x) if x == x else x for x in out[column].astype(float)]
    return out


def learner(seed=SEED):
    return RandomForestClassifier(n_estimators=500, min_samples_leaf=5, n_jobs=4, random_state=seed)


def chromosome_folds(chromosomes, folds=FOLDS):
    """Fold of each gene: sorted chromosome index mod folds (chromosome-grouped cross-validation)."""
    names = sorted(set(chromosomes))
    if len(names) < 2:raise ValueError("Grouped cross-validation needs at least two chromosomes")
    index = {c:i % min(folds, len(names)) for i, c in enumerate(names)}
    return np.array([index[c] for c in chromosomes])


def hybrid_scores(train, score, seed=SEED):
    """train: rows with features and boolean 'essential'; score: rows to score.
    Returns (out-of-fold scores for train, scores for `score` from the model trained on all of train)."""
    X, y = train[RF_FEATURES+HMM_FEATURES].to_numpy(float), train.essential.to_numpy(bool)
    if y.all() or not y.any():raise ValueError("Training labels must include both classes")
    folds = chromosome_folds(train.chromosome.tolist()); oof = np.zeros(len(train))
    for k in np.unique(folds):
        tr, te = folds != k, folds == k
        oof[te] = learner(seed).fit(X[tr], y[tr]).predict_proba(X[te])[:, 1]
    full = learner(seed).fit(X, y)
    return oof, full.predict_proba(score[RF_FEATURES+HMM_FEATURES].to_numpy(float))[:, 1] if len(score) else np.zeros(0)


def build(hmm, rf=None, labels=None, train_tiers=("silver",), seed=SEED, features=None):
    """hmm: HMM genes table (ORF, Name, chromosome, call, status, HMM_FEATURES).
    features: SummaryTable feature table (gene features for the hybrid) (optional).
    rf: legacy classifier table with 'Standard name', 'RF - G4' and RF_FEATURES (optional; its gene features are
        used only when `features` is not given).
    labels: curated labels (gene_id or cglab_id, label, tier) (optional; needed with gene features for the hybrid).
    Returns (menu DataFrame, summary dict)."""
    d = hmm[["ORF", "Name", "chromosome", "status", "call"]+HMM_FEATURES].rename(columns={"call":"hmm_call", "status":"hmm_status"})
    d["hmm_score"] = d.P_ES_run; d["hmm_essential"] = d.hmm_call == "ES"
    models = ["hmm"]
    if features is not None:
        d = d.merge(summary_features(features), on="ORF", how="left")
    if rf is not None:
        r = rf.rename(columns={"Standard name":"ORF", "RF - G4":"rf_existing_score"})
        r = r[["ORF", "rf_existing_score"]+(RF_FEATURES if features is None else [])]
        d = d.merge(r, on="ORF", how="left"); models.append("rf_existing")
        d["rf_existing_essential"] = np.where(d.rf_existing_score.notna(), d.rf_existing_score >= THRESHOLD, np.nan)
    summary = dict(revision=REVISION, genes=len(d), models=models, threshold=THRESHOLD)
    if (features is not None or rf is not None) and labels is not None:
        lab = labels.rename(columns={"cglab_id":"ORF", "gene_id":"ORF"})[["ORF", "label", "tier"]]
        lab = lab[lab.label.isin(["essential", "nonessential"])]
        usable = (d.hmm_status == "called") & d[RF_FEATURES+HMM_FEATURES].notna().all(axis=1)
        # training genes in label-file order (the order the pre-registered benchmark used)
        train = lab[lab.tier.isin(train_tiers)].merge(d[usable], on="ORF").assign(essential=lambda t:t.label == "essential")
        score = d[usable & ~d.ORF.isin(train.ORF)]
        oof, full = hybrid_scores(train, score, seed)
        d = d.merge(pd.concat([pd.DataFrame(dict(ORF=train.ORF, hybrid_score=oof, hybrid_source="out-of-fold (training gene)")),
                               pd.DataFrame(dict(ORF=score.ORF, hybrid_score=full, hybrid_source="trained model"))]), on="ORF", how="left")
        d["hybrid_essential"] = np.where(d.hybrid_score.notna(), d.hybrid_score >= THRESHOLD, np.nan)
        d["agreement"] = [AGREEMENT[(bool(h), bool(y))] if y == y else None for h, y in zip(d.hmm_essential, d.hybrid_essential)]
        d.loc[d.hmm_status != "called", "agreement"] = None
        d = d.merge(lab.rename(columns={"label":"reference_label", "tier":"reference_tier"}), on="ORF", how="left")
        models.append("hybrid")
        summary |= dict(hybrid_training=dict(tiers=list(train_tiers), genes=len(train), essential=int(train.essential.sum()),
                                             folds=int(len(np.unique(chromosome_folds(train.chromosome.tolist())))), seed=seed),
                        agreement=d.agreement.value_counts().to_dict())
    summary |= dict(hmm_calls=d.hmm_call.value_counts().to_dict())
    lead = ["ORF", "Name", "chromosome", "hmm_score", "hmm_call"]
    lead += [c for c in ("hybrid_score", "hybrid_essential", "agreement", "rf_existing_score", "rf_existing_essential") if c in d]
    rest = [c for c in d.columns if c not in lead]
    return d[lead+rest], summary
