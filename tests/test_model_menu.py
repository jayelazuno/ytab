import numpy as np
import pandas as pd

from ytab.essentiality import model_menu as mm


def tables(n=240, seed=3):
    rng = np.random.default_rng(seed)
    ess = rng.random(n) < .3
    orfs = [f"G{i:04d}" for i in range(n)]
    hmm = pd.DataFrame(dict(ORF=orfs, Name=orfs, chromosome=[f"chr{i % 6}" for i in range(n)], status="called",
                            call=np.where(ess, "ES", "NE"), P_ES_run=np.where(ess, .9, .1)+rng.normal(0, .05, n),
                            P_ES=np.where(ess, .8, .1), P_GD=.05, P_NE=np.where(ess, .1, .8), P_GA=.05, mappable_fraction=1.))
    rf = pd.DataFrame({"Standard name":orfs, "RF - G4":np.where(ess, .8, .2)} | {f:rng.normal(np.where(ess, 1, 0), 1) for f in mm.RF_FEATURES})
    labels = pd.DataFrame(dict(cglab_id=orfs, label=np.where(ess, "essential", "nonessential"), tier=np.where(np.arange(n) < 200, "silver", "gold")))
    return hmm, rf, labels


def test_hmm_only_menu_needs_no_labels():
    hmm, _, _ = tables()
    menu, summary = mm.build(hmm)
    assert summary["models"] == ["hmm"] and "hybrid_score" not in menu and len(menu) == len(hmm)


def test_hybrid_scores_training_genes_out_of_fold():
    hmm, rf, labels = tables()
    menu, summary = mm.build(hmm, rf, labels)
    train = menu.reference_tier == "silver"
    assert (menu.hybrid_source[train] == "out-of-fold (training gene)").all()
    assert (menu.hybrid_source[~train] == "trained model").all()
    assert summary["hybrid_training"]["genes"] == 200
    # an out-of-fold score does not change when the gene's own label is flipped
    flipped = labels.copy(); flipped.loc[0, "label"] = "essential" if labels.label[0] == "nonessential" else "nonessential"
    again, _ = mm.build(hmm, rf, flipped)
    assert again.hybrid_score[0] == menu.hybrid_score[0]


def test_agreement_classes():
    hmm, rf, labels = tables()
    menu, _ = mm.build(hmm, rf, labels)
    expected = [mm.AGREEMENT[(h, bool(y))] for h, y in zip(menu.hmm_call == "ES", menu.hybrid_essential)]
    assert menu.agreement.tolist() == expected


def test_chromosome_folds_group_genes():
    folds = mm.chromosome_folds(["b", "a", "c", "a", "b"], folds=2)
    assert folds.tolist() == [1, 0, 0, 0, 1]


def test_summary_features_match_classifier_table_precision():
    import pandas as pd
    from ytab.essentiality import model_menu as mm
    table = pd.DataFrame(dict(standard_name=["g1"], coding_length=[900], hits=[3], reads=[40], max_free_region=[120],
                              freedom_index=[0.762593], neighborhood_index=[0.1162031], upstream_hits_100=[2]))
    f = mm.summary_features(table)
    assert list(f.columns) == ["ORF"]+mm.RF_FEATURES
    assert f.loc[0, "Freedom index"] == 0.76 and f.loc[0, "Neighborhood index"] == 0.116 and f.loc[0, "Length"] == 900
