"""
aggregate_article_sim_v3.py
===========================

Re-trains the Random Forest ensemble on `sim_article_v2_scores.csv` using a
single ReReReRe feature (z_rr) plus the four auxiliary detectors, dropping
the z_rr_iter (iterative) variant entirely.  The paper from this point on
talks about ONE ReReReRe — the standard `ReReReRe()` function — and
nothing else.

Outputs (replacing the v2 pickle):
  - article_assets_v3/article_data.pkl
  - sim_article_v3_metrics.csv         (per size,rate,rep,scenario)
  - sim_article_v3_ablation.csv        (per size,rate,rep,scenario)
  - sim_article_v3_perpat.csv          (per size,rate,rep,scenario,pattern,corruption)
  - robustness_calibration_v3.csv      (per size,rate,scenario,fold,calib_mode)
  - bootstrap_delta_mcc_v3.csv         (per size, scenario)

Plots regenerated (used by build_paper_v3.py):
  - article_assets_v3/fig01_mcc_by_size.png
  - article_assets_v3/fig04_per_pattern.png
  - article_assets_v3/fig08_delta_mcc_by_size.png
  - revision_assets_v3/figR2_calibration_strategies.png
  - revision_assets_v3/figR3_bootstrap_delta_mcc.png

Scenario B only (τ_GT = 0.40) — the deployment-realistic operational
definition; Scenario A is excluded by design.
"""

import os, sys, pickle, time
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl

from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import StratifiedKFold
from sklearn.metrics import (matthews_corrcoef, f1_score, roc_auc_score,
                             average_precision_score, confusion_matrix,
                             cohen_kappa_score)

mpl.rcParams.update({
    "font.size": 10,
    "axes.titlesize": 12,
    "axes.labelsize": 10,
    "legend.fontsize": 9,
    "figure.dpi": 110,
    "savefig.dpi": 150,
    "savefig.bbox": "tight",
})

ROOT       = r"C:\Users\vitto\Desktop\ReReReRe"
SCORES_CSV = os.path.join(ROOT, "sim_article_v2_scores.csv")
OUT_ASSETS = os.path.join(ROOT, "article_assets_v3")
OUT_REV    = os.path.join(ROOT, "revision_assets_v3")
os.makedirs(OUT_ASSETS, exist_ok=True)
os.makedirs(OUT_REV, exist_ok=True)

TAU_GT       = 0.40                                # Scenario B
RR_FEATURES  = ["z_rr", "irv", "longstring", "d2", "person_tot"]   # 5
NORR_FEATURES = ["irv", "longstring", "d2", "person_tot"]          # 4
RF_KWARGS    = dict(n_estimators=200, max_features="sqrt",
                    min_samples_leaf=2, n_jobs=-1, random_state=0)
N_BOOT       = 500


# --------------------------------------------------------------------------
# Helpers
# --------------------------------------------------------------------------
def metrics_at_threshold(y_true, p, thr):
    """Return all 13 metrics at a fixed decision threshold."""
    yhat = (p >= thr).astype(int)
    tn, fp, fn, tp = confusion_matrix(y_true, yhat, labels=[0, 1]).ravel()
    sens = tp / (tp + fn) if (tp + fn) else 0.0
    spec = tn / (tn + fp) if (tn + fp) else 0.0
    ppv  = tp / (tp + fp) if (tp + fp) else 0.0
    npv  = tn / (tn + fn) if (tn + fn) else 0.0
    f1   = f1_score(y_true, yhat, zero_division=0)
    f2   = (5 * tp) / (5 * tp + 4 * fn + fp) if (5 * tp + 4 * fn + fp) else 0.0
    bal  = 0.5 * (sens + spec)
    yj   = sens + spec - 1.0
    gm   = (sens * spec) ** 0.5
    mcc  = matthews_corrcoef(y_true, yhat) if len(np.unique(yhat)) > 1 else 0.0
    kap  = cohen_kappa_score(y_true, yhat)
    if len(np.unique(y_true)) > 1:
        auc   = roc_auc_score(y_true, p)
        auprc = average_precision_score(y_true, p)
    else:
        auc = auprc = float("nan")
    return dict(sens=sens, spec=spec, ppv=ppv, npv=npv, bal_acc=bal,
                youden_j=yj, g_mean=gm, f1=f1, f2=f2, mcc=mcc,
                kappa=kap, auc=auc, auprc=auprc, tau_op=thr)


def calib_thresholds(p, y_true, p_est):
    """Compute thresholds for the four calibration strategies."""
    return {
        # 95th percentile of clean-only predictions  → label-dependent
        "oracle_clean": float(np.quantile(p[y_true == 0], 0.95)) if (y_true == 0).any() else 0.5,
        # 95th percentile of all predictions
        "blind95":      float(np.quantile(p, 0.95)),
        # (1 - p_est)-quantile of all predictions
        "rate_aware":   float(np.quantile(p, 1.0 - p_est)) if 0 < p_est < 1 else 0.5,
        # constant 0.5
        "fixed05":      0.5,
    }


def fit_oof(X, y, seed=0):
    """5-fold stratified OOF predicted probabilities."""
    p_oof = np.zeros(len(y))
    fold_id = np.zeros(len(y), dtype=int)
    skf = StratifiedKFold(n_splits=5, shuffle=True, random_state=seed)
    for k, (tr, te) in enumerate(skf.split(X, y)):
        rf = RandomForestClassifier(**RF_KWARGS)
        rf.fit(X[tr], y[tr])
        p_oof[te] = rf.predict_proba(X[te])[:, 1]
        fold_id[te] = k
    return p_oof, fold_id


# --------------------------------------------------------------------------
# Load data
# --------------------------------------------------------------------------
try:
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
except Exception:
    pass
print(f"[{time.strftime('%H:%M:%S')}] Loading scores ...")
df = pd.read_csv(SCORES_CSV)
# The v2 R simulator stored z_rr with the sign flipped (high = careless) to
# match a "classifier-friendly" direction. Restore the canonical ReReReRe()
# convention where LOW z_rr = careless, so the saved data and the paper's
# narrative agree. (The RF is sign-invariant so this does not change any
# prediction; it only affects how z_rr is described.)
df["z_rr"] = -df["z_rr"]
df["careless"] = (df["corruption"] > TAU_GT).astype(int)
print(f"   shape={df.shape}, careless rate (overall) = {df['careless'].mean():.3f}")
print(f"   z_rr clean median={df.loc[df['careless']==0,'z_rr'].median():.2f}, "
      f"careless median={df.loc[df['careless']==1,'z_rr'].median():.2f} "
      f"(low = careless)")

SIZES = sorted(df["size"].unique())
RATES = sorted(df["rate"].unique())
REPS  = sorted(df["rep"].unique())


# --------------------------------------------------------------------------
# Within-(size,rate,rep) 5-fold CV — Full and No-RR
# --------------------------------------------------------------------------
print(f"[{time.strftime('%H:%M:%S')}] Per-cell 5-fold CV (Full=5 feat & No-RR=4 feat) ...")

metrics_rows  = []
ablation_rows = []
perpat_rows   = []

t0 = time.time()
n_cells = len(SIZES) * len(RATES) * len(REPS)
done = 0

# Persist OOF predictions for downstream pooling
oof_full_all = np.zeros(len(df))
oof_norr_all = np.zeros(len(df))

for size in SIZES:
    for rate in RATES:
        for rep in REPS:
            mask = (df["size"] == size) & (df["rate"] == rate) & (df["rep"] == rep)
            sub  = df.loc[mask]
            y    = sub["careless"].values
            X_f  = sub[RR_FEATURES].values
            X_n  = sub[NORR_FEATURES].values

            p_full, _ = fit_oof(X_f, y, seed=int(rep) * 13 + 1)
            p_norr, _ = fit_oof(X_n, y, seed=int(rep) * 13 + 1)

            # Persist for pooling (Bootstrap, calibration tests)
            oof_full_all[mask.values] = p_full
            oof_norr_all[mask.values] = p_norr

            # Operating point: rate-aware calibration. We use the empirical
            # careless rate (true label rate within the cell) as `rate_hat` —
            # this is the v2.1 convention. In deployment the user supplies it
            # from domain knowledge / pilot data (typical careless rates 5-30%).
            p_est = max(0.005, min(0.995, float(y.mean())))

            # ---- Per-cell aggregated metrics under FPR=5% (clean = y==0) ----
            thr_oracle = float(np.quantile(p_full[y == 0], 0.95)) if (y == 0).any() else 0.5
            mF = metrics_at_threshold(y, p_full, thr_oracle)
            mF.update(size=size, rate=rate, rep=rep, scenario="B",
                      tau=TAU_GT, n_pos=int(y.sum()), n_neg=int(len(y) - y.sum()))
            metrics_rows.append(mF)

            # ---- Ablation: Full vs No-RR at the SAME oracle FPR=5% threshold ----
            thr_oracle_n = float(np.quantile(p_norr[y == 0], 0.95)) if (y == 0).any() else 0.5
            mcc_full = matthews_corrcoef(y, (p_full >= thr_oracle).astype(int))
            mcc_norr = matthews_corrcoef(y, (p_norr >= thr_oracle_n).astype(int))
            f1_full  = f1_score(y, (p_full >= thr_oracle).astype(int), zero_division=0)
            f1_norr  = f1_score(y, (p_norr >= thr_oracle_n).astype(int), zero_division=0)
            ablation_rows.append(dict(
                size=size, rate=rate, rep=rep, scenario="B", tau=TAU_GT,
                mcc_full=mcc_full, f1_full=f1_full,
                mcc_norr=mcc_norr, f1_norr=f1_norr,
                delta_mcc=mcc_full - mcc_norr,
                delta_f1=f1_full - f1_norr,
            ))

            # ---- Per-pattern detection at rate-aware threshold (Full vs No-RR) ----
            thr_ra_full = float(np.quantile(p_full, 1.0 - p_est)) if 0 < p_est < 1 else 0.5
            thr_ra_norr = float(np.quantile(p_norr, 1.0 - p_est)) if 0 < p_est < 1 else 0.5
            for (pat, corr), gsub in sub.groupby(["pattern", "corruption"]):
                idx = gsub.index
                pos = (gsub["careless"] == 1).values
                if pat == "clean" or len(idx) == 0:
                    continue
                mlocal = mask.values
                # Use the OOF predictions slice
                ix = np.where(mlocal)[0]
                # Index back into the cell-local p_full/p_norr arrays:
                # build local positional indexer
                local_pos = np.searchsorted(np.where(mlocal)[0], idx)
                pf = p_full[local_pos]
                pn = p_norr[local_pos]
                rate_full = float((pf >= thr_ra_full).mean())
                rate_norr = float((pn >= thr_ra_norr).mean())
                perpat_rows.append(dict(
                    size=size, rate=rate, rep=rep, scenario="B",
                    pattern=pat, corruption=float(corr),
                    rate_full=rate_full, rate_norr=rate_norr,
                    delta=rate_full - rate_norr,
                ))

            done += 1
            if done % 32 == 0:
                el = time.time() - t0
                print(f"   {done}/{n_cells} cells  elapsed={el:.0f}s")

el = time.time() - t0
print(f"   DONE per-cell CV in {el:.0f}s")

metrics  = pd.DataFrame(metrics_rows)
ablation = pd.DataFrame(ablation_rows)
perpat   = pd.DataFrame(perpat_rows)

metrics.to_csv(os.path.join(ROOT, "sim_article_v3_metrics.csv"),  index=False)
ablation.to_csv(os.path.join(ROOT, "sim_article_v3_ablation.csv"), index=False)
perpat.to_csv(os.path.join(ROOT,   "sim_article_v3_perpat.csv"),   index=False)

# Persist OOF prediction arrays so a downstream failure in the calibration /
# bootstrap / plotting phases does not require redoing the 35-minute CV.
np.save(os.path.join(ROOT, "v3_oof_full.npy"), oof_full_all)
np.save(os.path.join(ROOT, "v3_oof_norr.npy"), oof_norr_all)
print(f"   Saved OOF predictions to v3_oof_full.npy / v3_oof_norr.npy")


# --------------------------------------------------------------------------
# Robustness — calibration comparison, pooled within (size,rate)
# --------------------------------------------------------------------------
print(f"[{time.strftime('%H:%M:%S')}] Calibration comparison at the (size,rate) level ...")

calib_rows = []

for size in SIZES:
    for rate in RATES:
        mask = (df["size"] == size) & (df["rate"] == rate)
        sub  = df.loc[mask].reset_index(drop=True)
        y    = sub["careless"].values
        Xf   = sub[RR_FEATURES].values

        # one 5-fold CV pooled across reps for the calibration sweep
        skf = StratifiedKFold(n_splits=5, shuffle=True, random_state=42)
        for k, (tr, te) in enumerate(skf.split(Xf, y)):
            rf = RandomForestClassifier(**RF_KWARGS)
            rf.fit(Xf[tr], y[tr])
            p_te = rf.predict_proba(Xf[te])[:, 1]
            y_te = y[te]
            # rate_aware: use the empirical training-fold careless rate,
            # matching v2.1 (`rate_hat = y_tr.mean()`). In deployment this is
            # supplied by the user from a pilot or domain prior.
            p_est  = max(0.005, min(0.995, float(y[tr].mean())))
            thrs = calib_thresholds(p_te, y_te, p_est)
            for mode, thr in thrs.items():
                m = metrics_at_threshold(y_te, p_te, thr)
                calib_rows.append(dict(
                    experiment="v3_calibration", size=size, rate=rate,
                    scenario="B", fold=k, calib_mode=mode,
                    mcc=m["mcc"], f1=m["f1"], auc=m["auc"],
                    sens=m["sens"], spec=m["spec"], thr=thr,
                    n_tr=len(tr), n_te=len(te),
                    rate_te=float(y_te.mean()),
                ))

calib = pd.DataFrame(calib_rows)
calib.to_csv(os.path.join(ROOT, "robustness_calibration_v3.csv"), index=False)


# --------------------------------------------------------------------------
# Bootstrap CIs on Δ MCC by size (Scenario B), using OOF predictions
# --------------------------------------------------------------------------
print(f"[{time.strftime('%H:%M:%S')}] Bootstrap 95% CIs on Delta MCC by size ...")

rng = np.random.default_rng(2026)
boot_rows = []

for size in SIZES:
    mask = (df["size"] == size).values
    y    = df.loc[mask, "careless"].values
    pF   = oof_full_all[mask]
    pN   = oof_norr_all[mask]

    # Use the per-cell oracle-FPR=5% threshold consistently
    # Compute ONE pooled threshold for full / no-RR
    thr_F = float(np.quantile(pF[y == 0], 0.95))
    thr_N = float(np.quantile(pN[y == 0], 0.95))

    mcc_F_pt = matthews_corrcoef(y, (pF >= thr_F).astype(int))
    mcc_N_pt = matthews_corrcoef(y, (pN >= thr_N).astype(int))
    delta_pt = mcc_F_pt - mcc_N_pt

    n = len(y)
    deltas = np.empty(N_BOOT)
    for b in range(N_BOOT):
        idx  = rng.integers(0, n, n)
        yb, pFb, pNb = y[idx], pF[idx], pN[idx]
        if len(np.unique(yb)) < 2:
            deltas[b] = np.nan
            continue
        thrFb = float(np.quantile(pFb[yb == 0], 0.95)) if (yb == 0).any() else thr_F
        thrNb = float(np.quantile(pNb[yb == 0], 0.95)) if (yb == 0).any() else thr_N
        mF = matthews_corrcoef(yb, (pFb >= thrFb).astype(int))
        mN = matthews_corrcoef(yb, (pNb >= thrNb).astype(int))
        deltas[b] = mF - mN
    deltas = deltas[~np.isnan(deltas)]
    lo, hi = np.quantile(deltas, [0.025, 0.975])
    boot_rows.append(dict(
        scenario="B", size=size, n=n,
        delta_pt=float(delta_pt),
        delta_lo=float(lo), delta_hi=float(hi),
        delta_p_gt0=float((deltas > 0).mean()),
        mcc_full=float(mcc_F_pt), mcc_norr=float(mcc_N_pt),
    ))

boot = pd.DataFrame(boot_rows)
boot.to_csv(os.path.join(ROOT, "bootstrap_delta_mcc_v3.csv"), index=False)
print(boot.to_string(index=False))


# --------------------------------------------------------------------------
# Aggregations
# --------------------------------------------------------------------------
def agg(d, by, cols):
    g  = d.groupby(by)[cols]
    out = g.mean().reset_index()
    sd = g.std().reset_index()
    se = g.sem().reset_index()
    for c in cols:
        out[c + "_sd"] = sd[c]
        out[c + "_se"] = se[c]
    return out


metric_cols = ["sens", "spec", "ppv", "npv", "bal_acc", "youden_j",
               "g_mean", "f1", "f2", "mcc", "kappa", "auc", "auprc"]
abl_cols    = ["mcc_full", "f1_full", "mcc_norr", "f1_norr",
               "delta_mcc", "delta_f1"]

metrics_by_size      = agg(metrics, ["size", "scenario"], metric_cols)
metrics_by_size_rate = agg(metrics, ["size", "rate", "scenario"], metric_cols)
metrics_by_rate      = agg(metrics, ["rate", "scenario"], metric_cols)
metrics_overall      = agg(metrics, ["scenario"], metric_cols)

abl_by_size      = agg(ablation, ["size", "scenario"], abl_cols)
abl_by_size_rate = agg(ablation, ["size", "rate", "scenario"], abl_cols)
abl_overall      = agg(ablation, ["scenario"], abl_cols)

pp_agg = (perpat.groupby(["scenario", "pattern", "corruption"])
                 [["rate_full", "rate_norr", "delta"]]
                 .mean().reset_index())


# Pooled rate_aware (Scenario B)
ra_pool = (calib[calib.calib_mode == "rate_aware"]
           .groupby("scenario")
           .agg(mcc=("mcc", "mean"), f1=("f1", "mean"),
                sens=("sens", "mean"), spec=("spec", "mean"),
                auc=("auc", "mean")).reset_index())
print("\nRate-aware pooled (Scenario B):")
print(ra_pool.to_string(index=False))


# --------------------------------------------------------------------------
# Plots
# --------------------------------------------------------------------------
print(f"[{time.strftime('%H:%M:%S')}] Plots ...")

# Per-size rate-aware metrics (computed from calib)
ra_size = (calib[calib.calib_mode == "rate_aware"]
           .groupby(["size", "scenario"])
           .agg(mcc=("mcc", "mean"), mcc_sd=("mcc", "std"),
                f1=("f1", "mean"),
                sens=("sens", "mean"), spec=("spec", "mean"),
                auc=("auc", "mean"), n=("mcc", "count"))
           .reset_index())

# Fig 01: MCC by size (rate-aware, Scenario B)
fig, ax = plt.subplots(figsize=(7, 4.4))
sub = ra_size[ra_size.scenario == "B"].sort_values("size")
ax.errorbar(sub["size"], sub["mcc"], yerr=sub["mcc_sd"] / np.sqrt(sub["n"]),
            marker="o", linewidth=2, markersize=7, capsize=4, color="#0a3a6b",
            label="Rate-aware calibration (deployable)")
ax.set_xlabel("Total questionnaire items")
ax.set_ylabel("MCC at rate-aware threshold")
ax.set_title("Detection quality scales smoothly with questionnaire length")
ax.set_xticks(SIZES); ax.grid(alpha=0.3); ax.legend()
ax.set_ylim(0.3, 0.9)
fig.savefig(os.path.join(OUT_ASSETS, "fig01_mcc_by_size.png"))
plt.close(fig)

# Fig 04: per-pattern detection (rate-aware)
PATS = ["pure_straight", "acquiescent", "longstring",
        "fatigue", "mixed", "random"]
fig, axes = plt.subplots(2, 3, figsize=(13.0, 7.0), sharex=True, sharey=True)
for ax, pat in zip(axes.flat, PATS):
    s = pp_agg[(pp_agg.scenario == "B") & (pp_agg.pattern == pat)].sort_values("corruption")
    ax.plot(s.corruption, s.rate_full, marker="o", lw=2, color="#2ca02c",
            label="Full ensemble (RR + 4 aux.)")
    ax.plot(s.corruption, s.rate_norr, marker="s", lw=2, color="#d62728",
            label="Without ReReReRe (4 aux. only)")
    ax.set_title(pat); ax.set_ylim(0, 1.05); ax.grid(alpha=0.3)
axes[0, 0].legend(fontsize=8)
for ax in axes[1]:
    ax.set_xlabel("Corruption fraction")
for ax in axes[:, 0]:
    ax.set_ylabel("Detection rate")
fig.suptitle("Per-pattern detection (Scenario B): RR rescues random/mixed",
             y=1.02, fontsize=13)
fig.tight_layout()
fig.savefig(os.path.join(OUT_ASSETS, "fig04_per_pattern.png"))
plt.close(fig)

# Fig 08: Δ MCC from RR by size (Scenario B)
fig, ax = plt.subplots(figsize=(7.0, 4.3))
df_b = abl_by_size[abl_by_size.scenario == "B"].sort_values("size")
ax.errorbar(df_b["size"], df_b["delta_mcc"], yerr=df_b["delta_mcc_se"],
            marker="o", linewidth=2, capsize=4, color="#0a3a6b")
ax.axhline(0, color="grey", lw=0.7, ls="--")
ax.set_xlabel("Total items")
ax.set_ylabel("Δ MCC from including ReReReRe in the ensemble")
ax.set_title("ReReReRe contribution scales with questionnaire length")
ax.set_xticks(SIZES); ax.grid(alpha=0.3)
fig.savefig(os.path.join(OUT_ASSETS, "fig08_delta_mcc_by_size.png"))
plt.close(fig)

# figR2: calibration strategies
fig, ax = plt.subplots(figsize=(8, 5))
colors = {"oracle_clean": "#1f77b4", "blind95": "#ff7f0e",
          "rate_aware": "#2ca02c",  "fixed05": "#d62728"}
labels = {"oracle_clean": "oracle FPR=5% (label-aware)",
          "blind95":      "blind95",
          "rate_aware":   "rate-aware (deployable)",
          "fixed05":      "fixed 0.5"}
for mode in ["oracle_clean", "blind95", "rate_aware", "fixed05"]:
    s = (calib[(calib.scenario == "B") & (calib.calib_mode == mode)]
         .groupby("size")["mcc"].mean().reset_index().sort_values("size"))
    ax.plot(s["size"], s["mcc"], marker="o", lw=2,
            color=colors[mode], label=labels[mode])
ax.set_xlabel("Questionnaire length (items)")
ax.set_ylabel("MCC")
ax.set_title("Calibration strategies (Scenario B, τ_GT=0.40)")
ax.set_xticks(SIZES); ax.grid(alpha=0.3); ax.legend()
fig.savefig(os.path.join(OUT_REV, "figR2_calibration_strategies.png"))
plt.close(fig)

# figR3: bootstrap Δ MCC
fig, ax = plt.subplots(figsize=(8, 5))
b = boot.sort_values("size")
ax.errorbar(b["size"], b["delta_pt"],
            yerr=[b["delta_pt"] - b["delta_lo"], b["delta_hi"] - b["delta_pt"]],
            marker="o", lw=2, capsize=4, color="#0a3a6b")
ax.axhline(0, color="grey", lw=0.7, ls="--")
ax.set_xlabel("Questionnaire length (items)")
ax.set_ylabel("Δ MCC = MCC(Full 5-feat) − MCC(No-RR 4-feat)")
ax.set_title("Bootstrap 95% CI on the size-conditional ReReReRe contribution")
ax.set_xticks(SIZES); ax.grid(alpha=0.3)
fig.savefig(os.path.join(OUT_REV, "figR3_bootstrap_delta_mcc.png"))
plt.close(fig)


# --------------------------------------------------------------------------
# Pickle for build_paper_v3.py
# --------------------------------------------------------------------------
out = {
    "metrics_overall":      metrics_overall,
    "metrics_by_size":      metrics_by_size,
    "metrics_by_size_rate": metrics_by_size_rate,
    "metrics_by_rate":      metrics_by_rate,
    "abl_overall":          abl_overall,
    "abl_by_size":          abl_by_size,
    "abl_by_size_rate":     abl_by_size_rate,
    "pp_agg":               pp_agg,
    "n_runs":               int(len(metrics)),
    "n_respondents":        int(len(df)),
    "sizes":                SIZES,
    "rates":                RATES,
    "n_reps":               int(max(REPS)),
    "ra_pool":              ra_pool,
    "ra_size":              ra_size,
    "calib":                calib,
    "boot":                 boot,
}
with open(os.path.join(OUT_ASSETS, "article_data.pkl"), "wb") as f:
    pickle.dump(out, f)

print("\nFiles written:")
print(" ", os.path.join(ROOT, "sim_article_v3_metrics.csv"))
print(" ", os.path.join(ROOT, "sim_article_v3_ablation.csv"))
print(" ", os.path.join(ROOT, "sim_article_v3_perpat.csv"))
print(" ", os.path.join(ROOT, "robustness_calibration_v3.csv"))
print(" ", os.path.join(ROOT, "bootstrap_delta_mcc_v3.csv"))
print(" ", os.path.join(OUT_ASSETS, "article_data.pkl"))
print(" ", os.path.join(OUT_ASSETS, "fig01_mcc_by_size.png"))
print(" ", os.path.join(OUT_ASSETS, "fig04_per_pattern.png"))
print(" ", os.path.join(OUT_ASSETS, "fig08_delta_mcc_by_size.png"))
print(" ", os.path.join(OUT_REV,    "figR2_calibration_strategies.png"))
print(" ", os.path.join(OUT_REV,    "figR3_bootstrap_delta_mcc.png"))
