"""index_map_scree.py — how many components does the careless-index space have?
Scree plot of the pooled index correlation matrix with two decision rules:
  - Kaiser (eigenvalue >= 1), and
  - parallel analysis: the 95th percentile of eigenvalues from random data of the SAME
    size (p = 9 indices, N = pooled respondents), which is the stricter, sampling-aware rule.
Writes figures/fig_scree.png (+ index_map_scree.csv).
"""
import numpy as np, pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

P = pd.read_csv("index_map_pooled_spearman.csv", index_col=0)
idx = [c for c in P.columns if c != "response_time"]      # 9 response-derived indices
R = P.loc[idx, idx].values
m = len(idx)
N = 33000                                                  # pooled respondents (approx.)

w = np.linalg.eigvalsh(R)[::-1]

# ---- parallel analysis: eigenvalues of random uncorrelated data of the same shape ----
rng = np.random.default_rng(20260728)
REPS = 200
sim = np.empty((REPS, m))
for r in range(REPS):
    X = rng.standard_normal((N, m))
    Rr = np.corrcoef(X, rowvar=False)
    sim[r] = np.linalg.eigvalsh(Rr)[::-1]
pa_mean = sim.mean(0)
pa_p95 = np.percentile(sim, 95, axis=0)

k_kaiser = int((w >= 1).sum())
k_pa = int((w > pa_p95).sum())
out = pd.DataFrame({"component": np.arange(1, m + 1), "eigenvalue": w,
                    "prop_var": w / m, "cum_var": np.cumsum(w) / m,
                    "pa_mean": pa_mean, "pa_p95": pa_p95})
out.to_csv("index_map_scree.csv", index=False)
print(out.round(3).to_string(index=False))
print("\nKaiser (>=1): %d components | parallel analysis (> 95th pct of random): %d components" % (k_kaiser, k_pa))

# ---------------------------- figure ----------------------------
plt.rcParams.update({"font.size": 9, "font.family": "serif"})
fig, ax = plt.subplots(figsize=(5.4, 3.5), dpi=200)
x = np.arange(1, m + 1)
ax.plot(x, w, "o-", color="#294060", lw=1.8, ms=5.5, zorder=3, label="observed eigenvalue")
ax.plot(x, pa_p95, "s--", color="#8f3535", lw=1.2, ms=4, zorder=2,
        label="parallel analysis (95th pct, random data)")
ax.axhline(1.0, color="#8a8578", lw=1, ls=":", zorder=1, label="Kaiser criterion (= 1)")
# shade the retained components
ax.axvspan(0.5, k_pa + 0.5, color="#e9eef4", zorder=0)
ax.annotate("%d components retained" % k_pa, xy=(k_pa / 2 + 0.5, w.max() * 0.92),
            ha="center", fontsize=8.5, color="#294060")
ax.set_xlabel("component"); ax.set_ylabel("eigenvalue")
ax.set_xticks(x); ax.set_xlim(0.5, m + 0.5); ax.set_ylim(0, max(w.max(), 1.2) * 1.12)
ax.set_title("Scree plot of the careless-index space (9 indices)", fontsize=10)
ax.legend(frameon=False, fontsize=8, loc="upper right")
for s in ("top", "right"): ax.spines[s].set_visible(False)
fig.tight_layout()
fig.savefig("figures/fig_scree.png", bbox_inches="tight")
print("\nwrote figures/fig_scree.png and index_map_scree.csv")
