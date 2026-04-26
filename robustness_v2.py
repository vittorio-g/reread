"""Robustness analyses for v2 paper.

Runs the holdout experiments (A1 cross-rep, A2 cross-size, A3 cross-pattern),
the realistic-calibration experiment (B1) and the bootstrap CI on the
ablation (C1). Reads sim_article_v2_scores.csv (192k respondents).

Output:
- robustness_holdouts.csv   (A1, A2, A3)
- robustness_calibration.csv (B1)
- robustness_bootstrap.csv   (C1)
"""

import numpy as np
import pandas as pd
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import KFold, GroupKFold
from sklearn.metrics import (matthews_corrcoef, f1_score, roc_auc_score)
import warnings
warnings.filterwarnings("ignore")

# ----------------------------------------------------------------------
# Load + label functions
# ----------------------------------------------------------------------

print("Loading...")
df = pd.read_csv("sim_article_v2_scores.csv")
print(f"  {len(df):,} respondents in {df['size'].nunique()} sizes × {df['rate'].nunique()} rates × {df['rep'].nunique()} reps")

FEATS_FULL = ["z_rr", "z_rr_iter", "irv", "longstring", "d2", "person_tot"]
FEATS_TRIAD = ["z_rr_iter", "irv", "d2"]
FEATS_NORR = ["irv", "longstring", "d2", "person_tot"]

TAU_A = 0.60     # Scenario A: clean (corr=0) vs corrupted (corr > 0.60)
TAU_B = 0.40     # Scenario B: corr > 0.40 vs corr <= 0.40 (everyone in)

def labels_scenario(df_in, scenario):
    if scenario == "A":
        # exclude middle, keep clean (0) and corruption > 0.60
        mask = (df_in["corruption"] == 0) | (df_in["corruption"] > TAU_A)
        d = df_in[mask].copy()
        d["y"] = (d["corruption"] > TAU_A).astype(int)
        return d
    else:  # B
        d = df_in.copy()
        d["y"] = (d["corruption"] > TAU_B).astype(int)
        return d

def calibrate_threshold(probs, y_train, mode="oracle_clean", target_fpr=0.05, rate_hat=None):
    """Pick a threshold from training-fold predictions."""
    probs = np.asarray(probs)
    if mode == "oracle_clean":
        clean = probs[y_train == 0]
        return float(np.quantile(clean, 1 - target_fpr)) if len(clean) > 5 else 0.5
    if mode == "blind95":
        return float(np.quantile(probs, 0.95))
    if mode == "fixed05":
        return 0.5
    if mode == "rate_aware":
        # use empirical positive rate as the proportion to flag
        if rate_hat is None:
            rate_hat = float(y_train.mean())
        return float(np.quantile(probs, 1 - rate_hat))
    raise ValueError(mode)

def fit_predict(d_train, d_test, feats, threshold_mode="oracle_clean"):
    """Fit RF on train, predict on test, return metrics at calibrated threshold."""
    X_tr = d_train[feats].values
    y_tr = d_train["y"].values
    X_te = d_test[feats].values
    y_te = d_test["y"].values
    if y_tr.sum() == 0 or y_tr.sum() == len(y_tr):
        return None
    rf = RandomForestClassifier(n_estimators=200, max_features="sqrt",
                                n_jobs=-1, random_state=42, class_weight=None)
    rf.fit(X_tr, y_tr)
    p_tr = rf.predict_proba(X_tr)[:, 1]
    p_te = rf.predict_proba(X_te)[:, 1]
    rate_hat = float(y_tr.mean())
    thr = calibrate_threshold(p_tr, y_tr, threshold_mode, rate_hat=rate_hat)
    pred = (p_te >= thr).astype(int)
    auc = roc_auc_score(y_te, p_te) if len(np.unique(y_te)) > 1 else np.nan
    mcc = matthews_corrcoef(y_te, pred)
    f1 = f1_score(y_te, pred, zero_division=0)
    tp = int(((pred == 1) & (y_te == 1)).sum()); fp = int(((pred == 1) & (y_te == 0)).sum())
    tn = int(((pred == 0) & (y_te == 0)).sum()); fn = int(((pred == 0) & (y_te == 1)).sum())
    sens = tp / max(tp + fn, 1); spec = tn / max(tn + fp, 1)
    return dict(mcc=mcc, f1=f1, auc=auc, sens=sens, spec=spec, thr=thr,
                n_tr=len(y_tr), n_te=len(y_te), rate_te=float(y_te.mean()))

# ----------------------------------------------------------------------
# Experiment A1: cross-rep holdout (within (size, rate))
# ----------------------------------------------------------------------
print("\n[A1] Cross-rep holdout (8 reps train, 4 reps test)...")
rows_A1 = []
sizes = sorted(df["size"].unique()); rates = sorted(df["rate"].unique())
rng = np.random.RandomState(20260427)
for size in sizes:
    for rate in rates:
        for scenario in ["A", "B"]:
            d = labels_scenario(df[(df["size"] == size) & (df["rate"] == rate)], scenario)
            reps = sorted(d["rep"].unique())
            # 3 random splits of 8 train / 4 test reps
            for split_idx in range(3):
                perm = rng.permutation(reps)
                train_reps = perm[:8]; test_reps = perm[8:]
                d_tr = d[d["rep"].isin(train_reps)]; d_te = d[d["rep"].isin(test_reps)]
                for feat_name, feats in [("full", FEATS_FULL), ("triad", FEATS_TRIAD), ("norr", FEATS_NORR)]:
                    r = fit_predict(d_tr, d_te, feats)
                    if r is None: continue
                    rows_A1.append({"experiment": "A1_cross_rep", "size": size, "rate": rate,
                                    "scenario": scenario, "split": split_idx, "features": feat_name,
                                    **r})
    print(f"  size={size} done")
df_A1 = pd.DataFrame(rows_A1)
print(f"  A1 rows: {len(df_A1)}")

# ----------------------------------------------------------------------
# Experiment A2: cross-size holdout (pool across rates+reps)
# ----------------------------------------------------------------------
print("\n[A2] Cross-size holdout (small <-> large)...")
SMALL = [30, 50, 80, 100]; LARGE = [150, 200, 250, 300]
rows_A2 = []
for direction, train_sizes, test_sizes in [
    ("small_to_large", SMALL, LARGE),
    ("large_to_small", LARGE, SMALL)
]:
    for scenario in ["A", "B"]:
        d_tr = labels_scenario(df[df["size"].isin(train_sizes)], scenario)
        # report per-size MCC on test side
        for test_size in test_sizes:
            d_te = labels_scenario(df[df["size"] == test_size], scenario)
            for feat_name, feats in [("full", FEATS_FULL), ("triad", FEATS_TRIAD), ("norr", FEATS_NORR)]:
                r = fit_predict(d_tr, d_te, feats)
                if r is None: continue
                rows_A2.append({"experiment": "A2_cross_size", "direction": direction,
                                "size_test": test_size, "scenario": scenario,
                                "features": feat_name, **r})
    print(f"  {direction} done")
df_A2 = pd.DataFrame(rows_A2)
print(f"  A2 rows: {len(df_A2)}")

# ----------------------------------------------------------------------
# Experiment A3: cross-pattern holdout (leave-2-patterns-out)
# ----------------------------------------------------------------------
print("\n[A3] Cross-pattern holdout (leave-2-out)...")
careless_patterns = ["random", "longstring", "pure_straight", "acquiescent", "mixed", "fatigue"]
# 5 leave-2-out splits chosen to balance "consistent" vs "inconsistent" in train and test
splits = [
    {"name": "rand+long_held",   "test": ["random", "longstring"]},
    {"name": "mix+fatigue_held", "test": ["mixed", "fatigue"]},
    {"name": "straight+acq_held","test": ["pure_straight", "acquiescent"]},
    {"name": "rand+mix_held",    "test": ["random", "mixed"]},
    {"name": "long+fatigue_held","test": ["longstring", "fatigue"]},
]
rows_A3 = []
for sp in splits:
    test_pat = sp["test"]; train_pat = [p for p in careless_patterns if p not in test_pat]
    for scenario in ["A", "B"]:
        # train: clean + train_patterns
        keep_tr = (df["pattern"] == "clean") | (df["pattern"].isin(train_pat))
        d_tr = labels_scenario(df[keep_tr], scenario)
        # test: clean + test_patterns
        keep_te = (df["pattern"] == "clean") | (df["pattern"].isin(test_pat))
        d_te = labels_scenario(df[keep_te], scenario)
        for feat_name, feats in [("full", FEATS_FULL), ("triad", FEATS_TRIAD), ("norr", FEATS_NORR)]:
            r = fit_predict(d_tr, d_te, feats)
            if r is None: continue
            rows_A3.append({"experiment": "A3_cross_pattern", "split_name": sp["name"],
                            "test_patterns": "+".join(test_pat), "scenario": scenario,
                            "features": feat_name, **r})
    print(f"  split={sp['name']} done")
df_A3 = pd.DataFrame(rows_A3)
print(f"  A3 rows: {len(df_A3)}")

pd.concat([df_A1, df_A2, df_A3], ignore_index=True).to_csv("robustness_holdouts.csv", index=False)
print("Wrote robustness_holdouts.csv")

# ----------------------------------------------------------------------
# Experiment B1: realistic threshold calibration
# ----------------------------------------------------------------------
print("\n[B1] Calibration realism (4 strategies)...")
rows_B1 = []
calib_modes = ["oracle_clean", "blind95", "fixed05", "rate_aware"]
for size in sizes:
    for rate in rates:
        for scenario in ["A", "B"]:
            d = labels_scenario(df[(df["size"] == size) & (df["rate"] == rate)], scenario)
            # 5-fold CV
            kf = KFold(n_splits=5, shuffle=True, random_state=42)
            for fold_idx, (tr, te) in enumerate(kf.split(d)):
                d_tr, d_te = d.iloc[tr], d.iloc[te]
                for mode in calib_modes:
                    r = fit_predict(d_tr, d_te, FEATS_FULL, threshold_mode=mode)
                    if r is None: continue
                    rows_B1.append({"experiment": "B1_calibration", "size": size,
                                    "rate": rate, "scenario": scenario, "fold": fold_idx,
                                    "calib_mode": mode, **r})
    print(f"  size={size} done")
df_B1 = pd.DataFrame(rows_B1)
df_B1.to_csv("robustness_calibration.csv", index=False)
print(f"  B1 rows: {len(df_B1)} -> robustness_calibration.csv")

# ----------------------------------------------------------------------
# Experiment C1: bootstrap CI on Delta MCC (Full vs No-RR) per size
# ----------------------------------------------------------------------
print("\n[C1] Bootstrap CI for Delta MCC by size...")
N_BOOT = 500
rows_C1 = []
for size in sizes:
    for scenario in ["A", "B"]:
        d = labels_scenario(df[df["size"] == size], scenario).reset_index(drop=True)
        if len(d) < 100: continue
        # Train two RFs once on the full pooled data (5-fold OOF predictions)
        kf = KFold(n_splits=5, shuffle=True, random_state=42)
        oof_full = np.zeros(len(d)); oof_norr = np.zeros(len(d))
        for tr, te in kf.split(d):
            d_tr, d_te = d.iloc[tr], d.iloc[te]
            rf_f = RandomForestClassifier(n_estimators=200, max_features="sqrt",
                                          n_jobs=-1, random_state=42).fit(d_tr[FEATS_FULL], d_tr["y"])
            rf_n = RandomForestClassifier(n_estimators=200, max_features="sqrt",
                                          n_jobs=-1, random_state=42).fit(d_tr[FEATS_NORR], d_tr["y"])
            oof_full[te] = rf_f.predict_proba(d_te[FEATS_FULL])[:, 1]
            oof_norr[te] = rf_n.predict_proba(d_te[FEATS_NORR])[:, 1]
        y = d["y"].values
        # threshold = 95th percentile of clean OOF preds for each
        thr_f = float(np.quantile(oof_full[y == 0], 0.95))
        thr_n = float(np.quantile(oof_norr[y == 0], 0.95))
        pred_f = (oof_full >= thr_f).astype(int); pred_n = (oof_norr >= thr_n).astype(int)
        mcc_f_pt = matthews_corrcoef(y, pred_f); mcc_n_pt = matthews_corrcoef(y, pred_n)
        # bootstrap on respondent indices
        rng = np.random.RandomState(20260427 + size + (0 if scenario == "A" else 1))
        boots = []
        N = len(y)
        for b in range(N_BOOT):
            idx = rng.randint(0, N, size=N)
            yb, pf, pn = y[idx], pred_f[idx], pred_n[idx]
            if len(np.unique(yb)) < 2: continue
            mf = matthews_corrcoef(yb, pf); mn = matthews_corrcoef(yb, pn)
            boots.append((mf, mn, mf - mn))
        boots = np.array(boots)
        mcc_full = boots[:, 0]; mcc_norr = boots[:, 1]; delta = boots[:, 2]
        rows_C1.append({"size": size, "scenario": scenario,
                        "mcc_full_pt": mcc_f_pt, "mcc_norr_pt": mcc_n_pt,
                        "delta_pt": mcc_f_pt - mcc_n_pt,
                        "mcc_full_lo": float(np.quantile(mcc_full, 0.025)),
                        "mcc_full_hi": float(np.quantile(mcc_full, 0.975)),
                        "mcc_norr_lo": float(np.quantile(mcc_norr, 0.025)),
                        "mcc_norr_hi": float(np.quantile(mcc_norr, 0.975)),
                        "delta_lo": float(np.quantile(delta, 0.025)),
                        "delta_hi": float(np.quantile(delta, 0.975)),
                        "delta_mean": float(delta.mean()),
                        "delta_p_gt0": float((delta > 0).mean()),
                        "n": N})
    print(f"  size={size} done")
df_C1 = pd.DataFrame(rows_C1)
df_C1.to_csv("robustness_bootstrap.csv", index=False)
print(f"  C1 rows: {len(df_C1)} -> robustness_bootstrap.csv")

print("\nAll done.")
