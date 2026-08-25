"""GT-free conditional structure: slice the sample by quantile levels of EACH index and
look at how the OTHER indices correlate within each slice. Conditioning on an index partials
out the shared 'carelessness' axis, so within-slice correlations reveal coupling BEYOND
co-elevation. Trend from low slice (careful end) to high slice (careless end) answers:
does the index structure tighten toward the careless extreme? Pooled over all 15 datasets."""
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
plt.rcParams.update({"font.family": "serif", "font.size": 10})

IDX = ["rc","synonyms","antonyms","even_odd","rpr","person_total","d2","longstring","irv"]
LAB = {"rc":"rc","synonyms":"synonyms","antonyms":"antonyms","even_odd":"even-odd","rpr":"RPR",
       "person_total":"person-total","d2":"D²","longstring":"LongString","irv":"IRV"}
long = pd.read_csv("index_map_long.csv"); long[IDX] = long[IDX].apply(pd.to_numeric, errors="coerce")
datasets = list(dict.fromkeys(long.dataset)); m = len(IDX)
Q = 5
fisher = lambda r: np.arctanh(np.clip(r, -0.999, 0.999)); ifisher = np.tanh

def rank_bins(x, q):
    r = pd.Series(x).rank(method="first").values
    n = np.isfinite(x).sum()
    b = np.floor((r - 1) / len(x) * q).astype(int); b[~np.isfinite(x)] = -1
    return np.clip(b, -1, q - 1)

def offmean(C):
    mm = C.shape[0]; return np.nanmean(np.abs(C[~np.eye(mm, dtype=bool)]))

# [stratifier, bin] -> pooled mean|off-diag among OTHERS|
grid_z = np.zeros((m, Q)); grid_n = np.zeros((m, Q))
# pooled 8x8 among-others matrices at low-third / high-third, for two focus stratifiers
focus = {"person_total": {}, "rc": {}}
for fs in focus:
    focus[fs] = {"lo_z": None, "lo_n": None, "hi_z": None, "hi_n": None}

for ds in datasets:
    sub = long[long.dataset == ds]
    for si, S in enumerate(IDX):
        x = sub[S].values
        if np.isfinite(x).sum() < 100: continue
        b = rank_bins(x, Q)
        others = [c for c in IDX if c != S]
        for q in range(Q):
            g = sub[others][b == q]
            if len(g) < 25: continue
            C = g.corr(method="spearman").values
            om = offmean(C)
            if np.isfinite(om): grid_z[si, q] += fisher(om); grid_n[si, q] += 1
        # focus: low third vs high third full among-others matrix
        if S in focus:
            b3 = rank_bins(x, 3)
            for tag, qq in [("lo", 0), ("hi", 2)]:
                g = sub[others][b3 == qq]
                if len(g) < 25: continue
                C = g.corr(method="spearman").values
                Z = fisher(C)
                key_z, key_n = focus[S][tag + "_z"], focus[S][tag + "_n"]
                if key_z is None:
                    focus[S][tag + "_z"] = np.where(np.isfinite(Z), Z, 0.0)
                    focus[S][tag + "_n"] = np.isfinite(Z).astype(float)
                else:
                    focus[S][tag + "_z"] += np.where(np.isfinite(Z), Z, 0.0)
                    focus[S][tag + "_n"] += np.isfinite(Z).astype(float)

grid = ifisher(np.divide(grid_z, grid_n, out=np.zeros_like(grid_z), where=grid_n > 0))
print("=== mean |off-diag corr among the OTHER 8 indices| by quantile of the stratifier ===")
print("stratifier      " + "  ".join(f"Q{q+1}" for q in range(Q)) + "   (Q1=careful end, Q5=careless end)")
for si, S in enumerate(IDX):
    print(f"{LAB[S]:<14} " + "  ".join(f"{grid[si,q]:.2f}" for q in range(Q)) + f"   trend {grid[si,Q-1]-grid[si,0]:+.2f}")
print(f"\nAVERAGE over stratifiers: " + "  ".join(f"{grid[:,q].mean():.2f}" for q in range(Q)) +
      f"   trend {grid[:,Q-1].mean()-grid[:,0].mean():+.2f}")

# ---- figure ----
fig = plt.figure(figsize=(15, 5))
axA = fig.add_subplot(1, 3, 1)
fam = {"rc":"#8f3535","synonyms":"#b45a5a","antonyms":"#b45a5a","even_odd":"#b45a5a","rpr":"#b45a5a",
       "person_total":"#33475f","d2":"#33475f","longstring":"#2f6b3a","irv":"#2f6b3a"}
for si, S in enumerate(IDX):
    axA.plot(range(1, Q+1), grid[si], marker="o", ms=4, color=fam[S], lw=1.4 if S in ("rc","person_total") else 0.9,
             alpha=1 if S in ("rc","person_total") else 0.6, label=LAB[S])
axA.plot(range(1, Q+1), grid.mean(0), marker="s", color="black", lw=2.2, label="mean", zorder=5)
axA.set_xlabel("quantile of stratifying index (1=careful → 5=careless)"); axA.set_ylabel("mean |corr| among the OTHER indices")
axA.set_title("(a) Does the index structure tighten\ntoward the careless end?", fontsize=10)
axA.legend(fontsize=6.5, ncol=2, loc="upper left"); axA.grid(alpha=.25); axA.set_xticks(range(1, Q+1))

# two heatmaps: among-others correlation at low vs high third of person-total
others = [c for c in IDX if c != "person_total"]
Plo = ifisher(np.divide(focus["person_total"]["lo_z"], focus["person_total"]["lo_n"],
                        out=np.zeros_like(focus["person_total"]["lo_z"]), where=focus["person_total"]["lo_n"] > 0)); np.fill_diagonal(Plo, 1)
Phi = ifisher(np.divide(focus["person_total"]["hi_z"], focus["person_total"]["hi_n"],
                        out=np.zeros_like(focus["person_total"]["hi_z"]), where=focus["person_total"]["hi_n"] > 0)); np.fill_diagonal(Phi, 1)
for ax, (Pm, ttl) in zip([fig.add_subplot(1,3,2), fig.add_subplot(1,3,3)],
                         [(Plo, "low person-total (careful end)"), (Phi, "high person-total (careless end)")]):
    im = ax.imshow(Pm, vmin=-1, vmax=1, cmap="RdBu_r"); mm = len(others)
    ax.set_xticks(range(mm)); ax.set_yticks(range(mm))
    ax.set_xticklabels([LAB[c] for c in others], rotation=45, ha="right", fontsize=7); ax.set_yticklabels([LAB[c] for c in others], fontsize=7)
    for i in range(mm):
        for j in range(mm):
            ax.text(j, i, f"{Pm[i,j]:.2f}", ha="center", va="center", color="white" if abs(Pm[i,j])>.55 else "black", fontsize=6)
    ax.set_title(f"(among the other 8 indices)\n{ttl}\nmean|off-diag|={offmean(Pm):.2f}", fontsize=8.5)
plt.tight_layout(); plt.savefig("index_gradient_fig.png", dpi=150); plt.savefig("C:/Users/vitto/Downloads/rerere_gradient.png", dpi=150)
print("\nwrote index_gradient_fig.png")

# specific pairs across person-total slices
print("\n=== key pairs among-others, low vs high person-total third ===")
oi = {c: i for i, c in enumerate(others)}
for a, b in [("rc","synonyms"),("rc","longstring"),("rc","irv"),("synonyms","antonyms"),("even_odd","rpr"),("d2","irv")]:
    print(f"  {LAB[a]}~{LAB[b]:<12} low {Plo[oi[a],oi[b]]:+.2f}   high {Phi[oi[a],oi[b]]:+.2f}")
