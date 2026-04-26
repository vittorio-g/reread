"""
aggregate_article_sim.py
Reads the four CSVs produced by Sim_Article_Big.R and produces:
  1. Aggregated tables (mean ± SD across reps)
  2. Plots used in the article PDF
  3. A pickle with everything pre-computed for the PDF builder

Run AFTER Sim_Article_Big.R has finished (or after a checkpoint save).
"""

import os
import pickle
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl

mpl.rcParams.update({
    "font.size": 10,
    "axes.titlesize": 12,
    "axes.labelsize": 10,
    "legend.fontsize": 9,
    "figure.dpi": 110,
    "savefig.dpi": 150,
    "savefig.bbox": "tight",
})

ROOT = r"C:\Users\vitto\Desktop\ReReReRe"
OUT_DIR = os.path.join(ROOT, "article_assets")
os.makedirs(OUT_DIR, exist_ok=True)

# ----------------------------------------------------------
# Load
# ----------------------------------------------------------
scores   = pd.read_csv(os.path.join(ROOT, "sim_article_scores.csv"))
metrics  = pd.read_csv(os.path.join(ROOT, "sim_article_metrics.csv"))
ablation = pd.read_csv(os.path.join(ROOT, "sim_article_ablation.csv"))
perpat   = pd.read_csv(os.path.join(ROOT, "sim_article_perpat.csv"))

print(f"scores   {scores.shape}")
print(f"metrics  {metrics.shape}")
print(f"ablation {ablation.shape}")
print(f"perpat   {perpat.shape}")

# ----------------------------------------------------------
# Aggregations
# ----------------------------------------------------------

def agg(df, groupcols, valcols):
    g = df.groupby(groupcols)[valcols]
    out = g.mean().reset_index()
    sd  = g.std().reset_index()
    for c in valcols:
        out[c + "_sd"] = sd[c]
    return out

metric_cols = ["sens", "spec", "ppv", "npv", "bal_acc", "youden_j",
               "g_mean", "f1", "f2", "mcc", "kappa", "auc", "auprc"]

# Headline: per (size, scenario), averaged over rep & rate
metrics_by_size = agg(metrics, ["size", "scenario"], metric_cols)
# Per (size, rate, scenario)
metrics_by_size_rate = agg(metrics, ["size", "rate", "scenario"], metric_cols)
# Overall (per scenario)
metrics_overall = agg(metrics, ["scenario"], metric_cols)

# Ablation per scenario × size
abl_cols = ["mcc_full", "f1_full", "mcc_norr", "f1_norr",
            "mcc_triad", "f1_triad", "delta_mcc", "delta_f1"]
abl_by_size = agg(ablation, ["size", "scenario"], abl_cols)
abl_overall = agg(ablation, ["scenario"], abl_cols)

# Per-pattern (mean detection rates)
pp_agg = (perpat.groupby(["scenario", "pattern", "corruption"])
                 [["rate_full", "rate_norr", "delta"]]
                 .mean().reset_index())

# ----------------------------------------------------------
# Plots
# ----------------------------------------------------------
SIZE_ORDER = sorted(metrics["size"].unique())

# 1. MCC by questionnaire size, both scenarios
fig, ax = plt.subplots(figsize=(7.0, 4.2))
for sc, color, marker in [("A", "#1f77b4", "o"), ("B", "#d62728", "s")]:
    df = metrics_by_size[metrics_by_size.scenario == sc].sort_values("size")
    ax.errorbar(df["size"], df["mcc"], yerr=df["mcc_sd"],
                marker=marker, color=color, lw=2, capsize=3,
                label=f"Scenario {sc} (τ={'0.60' if sc=='A' else '0.40'})")
ax.set_xlabel("Total questionnaire items")
ax.set_ylabel("MCC at FPR=5% (5-fold CV)")
ax.set_title("MCC scales smoothly with questionnaire length")
ax.set_xticks(SIZE_ORDER)
ax.grid(alpha=0.3)
ax.legend()
fig.savefig(os.path.join(OUT_DIR, "fig01_mcc_by_size.png"))
plt.close(fig)

# 2. 13-metric panel for Scenario A across sizes
fig, axes = plt.subplots(3, 4, figsize=(13.0, 8.0), sharex=True)
metric_panel = ["sens", "spec", "ppv", "npv", "f1", "f2",
                "mcc", "kappa", "bal_acc", "youden_j", "auc", "auprc"]
for ax, mc in zip(axes.flat, metric_panel):
    for sc, color in [("A", "#1f77b4"), ("B", "#d62728")]:
        df = metrics_by_size[metrics_by_size.scenario == sc].sort_values("size")
        ax.errorbar(df["size"], df[mc], yerr=df[mc + "_sd"],
                    marker="o", color=color, lw=1.6, capsize=2,
                    label=f"Scenario {sc}")
    ax.set_title(mc.upper().replace("_", " "))
    ax.grid(alpha=0.3)
    ax.set_ylim(0, 1.02)
axes[0, 0].legend(fontsize=8)
for ax in axes[2]:
    ax.set_xlabel("Items")
for ax in axes[:, 0]:
    ax.set_ylabel("Score")
fig.suptitle("All 12 classification metrics across questionnaire sizes",
             y=1.0, fontsize=13)
fig.tight_layout()
fig.savefig(os.path.join(OUT_DIR, "fig02_metric_panel.png"))
plt.close(fig)

# 3. Ablation: full vs no-RR vs triad
fig, axes = plt.subplots(1, 2, figsize=(11.0, 4.2), sharey=True)
for ax, sc in zip(axes, ["A", "B"]):
    df = abl_by_size[abl_by_size.scenario == sc].sort_values("size")
    x = np.arange(len(df))
    w = 0.27
    ax.bar(x - w, df["mcc_full"], w, color="#2ca02c", label="Full ensemble (6 features)")
    ax.bar(x,     df["mcc_triad"], w, color="#1f77b4", label="Triad (iter + IRV + D²)")
    ax.bar(x + w, df["mcc_norr"], w, color="#d62728", label="No RR (4 aux features)")
    ax.set_xticks(x); ax.set_xticklabels(df["size"].astype(int))
    ax.set_xlabel("Items")
    ax.set_title(f"Scenario {sc}  (τ={'0.60' if sc=='A' else '0.40'})")
    ax.grid(alpha=0.3, axis="y")
axes[0].set_ylabel("MCC")
axes[0].legend(loc="lower right", fontsize=8)
fig.suptitle("Ablation: ReReReRe is essential, not redundant", y=1.02, fontsize=13)
fig.tight_layout()
fig.savefig(os.path.join(OUT_DIR, "fig03_ablation.png"))
plt.close(fig)

# 4. Per-pattern detection
PATS = ["pure_straight", "acquiescent", "longstring",
        "fatigue", "mixed", "random"]
fig, axes = plt.subplots(2, 3, figsize=(13.0, 7.0), sharex=True, sharey=True)
for ax, pat in zip(axes.flat, PATS):
    sub = pp_agg[(pp_agg.pattern == pat) & (pp_agg.scenario == "A")]\
              .sort_values("corruption")
    ax.plot(sub.corruption, sub.rate_full, marker="o", lw=2,
            color="#2ca02c", label="Full ensemble")
    ax.plot(sub.corruption, sub.rate_norr, marker="s", lw=2,
            color="#d62728", label="Without ReReReRe")
    ax.set_title(pat)
    ax.set_ylim(0, 1.05)
    ax.grid(alpha=0.3)
axes[0, 0].legend(fontsize=8)
for ax in axes[1]:
    ax.set_xlabel("Corruption fraction")
for ax in axes[:, 0]:
    ax.set_ylabel("Detection rate")
fig.suptitle("Per-pattern detection (Scenario A): RR rescues random/mixed",
             y=1.02, fontsize=13)
fig.tight_layout()
fig.savefig(os.path.join(OUT_DIR, "fig04_per_pattern.png"))
plt.close(fig)

# 5. Score distributions: clean vs careless on z_rr_iter (size=200, rate=0.40)
fig, ax = plt.subplots(figsize=(7.5, 4.2))
sub = scores[(scores["size"] == max(SIZE_ORDER)) & (scores["rate"] == 0.40)]
g_clean = sub[sub.corruption == 0]["z_rr_iter"]
g_low   = sub[(sub.corruption > 0) & (sub.corruption <= 0.40)]["z_rr_iter"]
g_mid   = sub[(sub.corruption > 0.40) & (sub.corruption <= 0.60)]["z_rr_iter"]
g_hi    = sub[sub.corruption > 0.60]["z_rr_iter"]
bins = np.linspace(min(sub["z_rr_iter"].min(), -5),
                   max(sub["z_rr_iter"].max(),  5), 50)
ax.hist(g_clean, bins=bins, alpha=0.55, label="Clean", color="#2ca02c", density=True)
ax.hist(g_low,   bins=bins, alpha=0.55, label="≤40% corrupted", color="#1f77b4", density=True)
ax.hist(g_mid,   bins=bins, alpha=0.55, label="40-60% corrupted", color="#ff7f0e", density=True)
ax.hist(g_hi,    bins=bins, alpha=0.55, label=">60% corrupted", color="#d62728", density=True)
ax.set_xlabel("ReReReRe iterative z-score (sign-flipped: high = more careless)")
ax.set_ylabel("Density")
ax.set_title("Score separation by corruption level (200 items, rate=40%)")
ax.legend()
ax.grid(alpha=0.3)
fig.savefig(os.path.join(OUT_DIR, "fig05_score_distribution.png"))
plt.close(fig)

# 6. Heatmap: MCC vs (size, rate)
piv = (metrics[metrics.scenario == "A"]
        .groupby(["size", "rate"])["mcc"].mean().unstack("rate"))
fig, ax = plt.subplots(figsize=(5.5, 4.5))
im = ax.imshow(piv.values, cmap="viridis", aspect="auto",
               vmin=0, vmax=1, origin="lower")
ax.set_xticks(range(piv.shape[1]))
ax.set_xticklabels([f"{r:.0%}" for r in piv.columns])
ax.set_yticks(range(piv.shape[0]))
ax.set_yticklabels(piv.index.astype(int))
ax.set_xlabel("Sample careless rate")
ax.set_ylabel("Items")
ax.set_title("MCC heatmap (Scenario A, FPR=5%)")
for i in range(piv.shape[0]):
    for j in range(piv.shape[1]):
        ax.text(j, i, f"{piv.values[i,j]:.2f}",
                ha="center", va="center",
                color="white" if piv.values[i, j] < 0.6 else "black")
fig.colorbar(im, ax=ax, label="MCC")
fig.savefig(os.path.join(OUT_DIR, "fig06_heatmap.png"))
plt.close(fig)

# ----------------------------------------------------------
# Pickle for the PDF builder
# ----------------------------------------------------------
with open(os.path.join(OUT_DIR, "article_data.pkl"), "wb") as f:
    pickle.dump({
        "metrics_overall":      metrics_overall,
        "metrics_by_size":      metrics_by_size,
        "metrics_by_size_rate": metrics_by_size_rate,
        "abl_overall":          abl_overall,
        "abl_by_size":          abl_by_size,
        "pp_agg":               pp_agg,
        "n_runs":               len(metrics) // 2,  # /2 because A and B
        "n_respondents":        scores.groupby(["size", "rate", "rep"]).size().sum(),
        "sizes":                SIZE_ORDER,
        "rates":                sorted(metrics["rate"].unique().tolist()),
    }, f)

print("\nWrote outputs to", OUT_DIR)
print("Files:", sorted(os.listdir(OUT_DIR)))
