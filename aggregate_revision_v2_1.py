"""Aggregate revision results (A1, A2, A3, B1, C1, D1) and produce
revision figures + a pickle for the v2.1 PDF builder.

Inputs (must exist before running):
  - sim_article_v2_metrics.csv  (v2 baseline)
  - sim_article_v2_ablation.csv (v2 baseline)
  - robustness_holdouts.csv     (A1 + A2 + A3)
  - robustness_calibration.csv  (B1)
  - robustness_bootstrap.csv    (C1)
  - external_rf_ensemble.csv    (D1, optional)
  - external_rf_bootstrap.csv   (D1, optional)

Outputs:
  - revision_assets/revision_data.pkl
  - revision_assets/figR1_holdout_summary.png
  - revision_assets/figR2_calibration_strategies.png
  - revision_assets/figR3_bootstrap_delta_mcc.png
  - revision_assets/figR4_external_rf.png
  - revision_assets/figR5_holdout_size_curve.png
"""

import os, pickle
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

OUT_DIR = "revision_assets"
os.makedirs(OUT_DIR, exist_ok=True)

# ---- baseline (v2) ----
v2_metrics = pd.read_csv("sim_article_v2_metrics.csv")
v2_abl = pd.read_csv("sim_article_v2_ablation.csv")

# ---- revision experiments ----
hold = pd.read_csv("robustness_holdouts.csv")
calib = pd.read_csv("robustness_calibration.csv")
boot = pd.read_csv("robustness_bootstrap.csv")
have_ext = os.path.exists("external_rf_ensemble.csv")
ext_rows = pd.read_csv("external_rf_ensemble.csv") if have_ext else None
ext_ci = pd.read_csv("external_rf_bootstrap.csv") if (have_ext and os.path.exists("external_rf_bootstrap.csv")) else None

# ============================================================
# Figure R1: holdout summary - cross-rep / cross-size / cross-pattern
# Three subpanels: each shows MCC for Full vs No-RR vs Triad
# ============================================================
fig, axes = plt.subplots(1, 3, figsize=(13, 4.2))

# panel 1: A1 cross-rep (per scenario, mean across all sizes/rates)
A1 = hold[hold["experiment"] == "A1_cross_rep"].copy()
piv1 = A1.groupby(["scenario", "features"])["mcc"].mean().unstack()
ax = axes[0]
xs = np.arange(2); w = 0.27
for i, feat in enumerate(["full", "triad", "norr"]):
    if feat in piv1.columns:
        ax.bar(xs + (i - 1) * w, piv1[feat].reindex(["A", "B"]).values,
               w, label=feat.replace("norr", "no-RR"))
ax.set_xticks(xs); ax.set_xticklabels(["Scen A", "Scen B"])
ax.set_ylabel("MCC (cross-rep, test reps)")
ax.set_title("A1 — Cross-rep holdout (8 train / 4 test)")
ax.legend(loc="lower right", fontsize=8)
ax.set_ylim(0, 1.0); ax.grid(axis="y", alpha=0.3)

# panel 2: A2 cross-size, MCC by test size, both directions
A2 = hold[hold["experiment"] == "A2_cross_size"].copy()
ax = axes[1]
for feat, ls in [("full", "-"), ("norr", "--")]:
    sub = A2[(A2["features"] == feat) & (A2["scenario"] == "A")]
    g = sub.groupby("size_test")["mcc"].mean()
    ax.plot(g.index, g.values, marker="o", linestyle=ls, label=f"{feat}")
ax.set_xlabel("Test-size (items)")
ax.set_ylabel("MCC (cross-size holdout, Scen A)")
ax.set_title("A2 — Train on small + large pooled, test by size")
ax.legend(fontsize=8); ax.grid(alpha=0.3)
ax.set_ylim(0, 1.0)

# panel 3: A3 cross-pattern, MCC per split
A3 = hold[hold["experiment"] == "A3_cross_pattern"].copy()
ax = axes[2]
splits = sorted(A3["split_name"].unique())
xs = np.arange(len(splits)); w = 0.27
for i, feat in enumerate(["full", "triad", "norr"]):
    vals = []
    for sp in splits:
        sub = A3[(A3["split_name"] == sp) & (A3["features"] == feat) & (A3["scenario"] == "A")]
        vals.append(float(sub["mcc"].mean()) if len(sub) else np.nan)
    ax.bar(xs + (i - 1) * w, vals, w, label=feat.replace("norr", "no-RR"))
ax.set_xticks(xs); ax.set_xticklabels([s.replace("_held", "") for s in splits], rotation=20, ha="right", fontsize=8)
ax.set_ylabel("MCC (cross-pattern, Scen A)")
ax.set_title("A3 — Leave-2-patterns-out")
ax.set_ylim(0, 1.0); ax.legend(fontsize=8); ax.grid(axis="y", alpha=0.3)

plt.tight_layout()
plt.savefig(f"{OUT_DIR}/figR1_holdout_summary.png", dpi=140)
plt.close()
print("Wrote figR1_holdout_summary.png")

# ============================================================
# Figure R2: calibration strategies — MCC by mode and size
# ============================================================
calib_avg = calib.groupby(["scenario", "calib_mode", "size"])["mcc"].mean().reset_index()
modes = ["oracle_clean", "blind95", "rate_aware", "fixed05"]
fig, axes = plt.subplots(1, 2, figsize=(11, 4.2), sharey=True)
for ax, sc, label in zip(axes, ["A", "B"], ["Scen A (τ=0.60)", "Scen B (τ=0.40)"]):
    for mode in modes:
        sub = calib_avg[(calib_avg["scenario"] == sc) & (calib_avg["calib_mode"] == mode)]
        ax.plot(sub["size"], sub["mcc"], marker="o", label=mode)
    ax.set_xlabel("Items")
    ax.set_title(label)
    ax.grid(alpha=0.3); ax.set_ylim(0, 1.0)
axes[0].set_ylabel("MCC")
axes[0].legend(fontsize=8, loc="lower right")
plt.tight_layout()
plt.savefig(f"{OUT_DIR}/figR2_calibration_strategies.png", dpi=140)
plt.close()
print("Wrote figR2_calibration_strategies.png")

# ============================================================
# Figure R3: bootstrap CI on Delta-MCC per size, both scenarios
# ============================================================
fig, ax = plt.subplots(figsize=(8, 4.2))
boot_a = boot[boot["scenario"] == "A"].sort_values("size")
boot_b = boot[boot["scenario"] == "B"].sort_values("size")
xa = boot_a["size"].values
ya = boot_a["delta_pt"].values
yerr_a = np.vstack([ya - boot_a["delta_lo"].values,
                    boot_a["delta_hi"].values - ya])
xb = boot_b["size"].values
yb = boot_b["delta_pt"].values
yerr_b = np.vstack([yb - boot_b["delta_lo"].values,
                    boot_b["delta_hi"].values - yb])
ax.errorbar(xa, ya, yerr=yerr_a, fmt="o-", capsize=4, label="Scen A (τ=0.60)")
ax.errorbar(xb, yb, yerr=yerr_b, fmt="s--", capsize=4, label="Scen B (τ=0.40)")
ax.axhline(0, color="black", linewidth=0.7, linestyle=":")
ax.set_xlabel("Items")
ax.set_ylabel("Δ MCC (Full ensemble − No-RR)")
ax.set_title("R3 — Bootstrap 95% CI on RR contribution by size (500 resamples)")
ax.legend(); ax.grid(alpha=0.3)
plt.tight_layout()
plt.savefig(f"{OUT_DIR}/figR3_bootstrap_delta_mcc.png", dpi=140)
plt.close()
print("Wrote figR3_bootstrap_delta_mcc.png")

# ============================================================
# Figure R4: external validation MCC (real data, RF on 6 features)
# ============================================================
if ext_rows is not None and len(ext_rows):
    methods_order = ["Full_6", "Triad_3", "NoRR_4", "z_rr_iter_only", "Mah_chi2_001"]
    pal = {"Full_6": "#1f77b4", "Triad_3": "#2ca02c", "NoRR_4": "#ff7f0e",
           "z_rr_iter_only": "#9467bd", "Mah_chi2_001": "#d62728"}
    datasets = list(ext_rows["dataset"].unique())
    fig, ax = plt.subplots(figsize=(9, 4.5))
    xs = np.arange(len(datasets))
    w = 0.16
    for i, m in enumerate(methods_order):
        vals = []
        for ds in datasets:
            sub = ext_rows[(ext_rows["dataset"] == ds) & (ext_rows["method"] == m)]
            vals.append(float(sub["mcc"].iloc[0]) if len(sub) else np.nan)
        ax.bar(xs + (i - 2) * w, vals, w, label=m, color=pal.get(m))
    ax.set_xticks(xs); ax.set_xticklabels(datasets, rotation=15, ha="right", fontsize=9)
    ax.set_ylabel("MCC (5-fold CV on real GT)")
    ax.set_title("R4 — External validation of RF ensemble on real datasets")
    ax.legend(fontsize=8, ncol=3, loc="lower center")
    ax.grid(axis="y", alpha=0.3)
    plt.tight_layout()
    plt.savefig(f"{OUT_DIR}/figR4_external_rf.png", dpi=140)
    plt.close()
    print("Wrote figR4_external_rf.png")
else:
    # placeholder
    fig, ax = plt.subplots(figsize=(8, 3))
    ax.text(0.5, 0.5, "External validation results not available", ha="center", va="center")
    ax.axis("off")
    plt.savefig(f"{OUT_DIR}/figR4_external_rf.png", dpi=140)
    plt.close()
    print("Wrote placeholder figR4_external_rf.png")

# ============================================================
# Figure R5: cross-size holdout, MCC by test size for both directions
# ============================================================
fig, ax = plt.subplots(figsize=(8, 4.2))
A2 = hold[hold["experiment"] == "A2_cross_size"].copy()
for direction, ls, marker in [("small_to_large", "-", "o"),
                               ("large_to_small", "--", "s")]:
    sub = A2[(A2["direction"] == direction) & (A2["features"] == "full") & (A2["scenario"] == "A")]
    g = sub.groupby("size_test")["mcc"].mean().sort_index()
    ax.plot(g.index, g.values, marker=marker, linestyle=ls,
            label=f"Full, {direction.replace('_', '→')}")
    sub_n = A2[(A2["direction"] == direction) & (A2["features"] == "norr") & (A2["scenario"] == "A")]
    gn = sub_n.groupby("size_test")["mcc"].mean().sort_index()
    ax.plot(gn.index, gn.values, marker=marker, linestyle=ls, alpha=0.4,
            label=f"No-RR, {direction.replace('_', '→')}")
ax.set_xlabel("Test-size (items)")
ax.set_ylabel("MCC (Scen A)")
ax.set_title("R5 — Cross-size holdout: small↔large generalisation")
ax.legend(fontsize=8); ax.grid(alpha=0.3)
ax.set_ylim(0, 1.0)
plt.tight_layout()
plt.savefig(f"{OUT_DIR}/figR5_holdout_size_curve.png", dpi=140)
plt.close()
print("Wrote figR5_holdout_size_curve.png")

# ============================================================
# Build pickle for the PDF
# ============================================================
def safe_mean(df, **kw):
    sub = df
    for k, v in kw.items():
        sub = sub[sub[k] == v]
    return float(sub["mcc"].mean()) if len(sub) else float("nan")

revision = {}

# A1 summary: mean MCC of full vs no-RR by scenario
A1 = hold[hold["experiment"] == "A1_cross_rep"]
revision["A1"] = {
    "scenA_full":  safe_mean(A1, scenario="A", features="full"),
    "scenA_triad": safe_mean(A1, scenario="A", features="triad"),
    "scenA_norr":  safe_mean(A1, scenario="A", features="norr"),
    "scenB_full":  safe_mean(A1, scenario="B", features="full"),
    "scenB_triad": safe_mean(A1, scenario="B", features="triad"),
    "scenB_norr":  safe_mean(A1, scenario="B", features="norr"),
}
revision["A1"]["scenA_delta"] = revision["A1"]["scenA_full"] - revision["A1"]["scenA_norr"]
revision["A1"]["scenB_delta"] = revision["A1"]["scenB_full"] - revision["A1"]["scenB_norr"]

# A2 by size (full only, scenario A)
A2 = hold[hold["experiment"] == "A2_cross_size"]
A2_summ = (A2[A2["scenario"] == "A"]
           .groupby(["direction", "features", "size_test"])["mcc"]
           .mean().reset_index())
revision["A2"] = A2_summ.to_dict("records")

# A3 per split (full vs no-RR, scenario A)
A3 = hold[hold["experiment"] == "A3_cross_pattern"]
A3_summ = (A3[A3["scenario"] == "A"]
           .groupby(["split_name", "features"])["mcc"]
           .mean().unstack().reset_index())
revision["A3"] = A3_summ.to_dict("records")

# B1: pooled MCC by mode and scenario
calib_summ = (calib.groupby(["scenario", "calib_mode"])
              [["mcc", "f1", "sens", "spec"]]
              .mean().reset_index())
revision["B1"] = calib_summ.to_dict("records")

# C1
revision["C1"] = boot.to_dict("records")

# D1 (if available)
if ext_rows is not None:
    revision["D1_rows"] = ext_rows.to_dict("records")
    if ext_ci is not None:
        revision["D1_ci"] = ext_ci.to_dict("records")

# Baseline v2 reference for the PDF
mA = v2_metrics[v2_metrics["scenario"] == "A"]
mB = v2_metrics[v2_metrics["scenario"] == "B"]
revision["v2_baseline"] = {
    "scenA_pooled_mcc": float(mA["mcc"].mean()),
    "scenB_pooled_mcc": float(mB["mcc"].mean()),
    "scenA_300_mcc": float(mA[mA["size"] == 300]["mcc"].mean()),
    "scenA_30_mcc":  float(mA[mA["size"] == 30]["mcc"].mean()),
}

with open(f"{OUT_DIR}/revision_data.pkl", "wb") as f:
    pickle.dump(revision, f)
print(f"\nSaved {OUT_DIR}/revision_data.pkl")
print("\nA1 summary:", revision["A1"])
print("Bootstrap CIs (Scenario A):")
for r in revision["C1"]:
    if r["scenario"] == "A":
        print(f"  size={r['size']}: Δ={r['delta_pt']:.3f} CI=[{r['delta_lo']:.3f},{r['delta_hi']:.3f}] P(Δ>0)={r['delta_p_gt0']:.3f}")
