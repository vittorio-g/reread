"""Two clean, paper-ready figures for the rc-vs-PT complementarity argument.
Fig 1: rc and PT are distinct constructs (pooled correlation heatmap + factor plot).
Fig 2: distinct validity — normativity(PT) predicts every external criterion (why it wins
       benchmarks), pair-consistency(rc) predicts induced inconsistency (factor x criterion)."""
import numpy as np, pandas as pd, os
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
plt.rcParams.update({"font.family": "serif", "font.size": 11})

RC = os.environ.get("RCLABEL", "rc")   # display label for the rr/rc index (data column stays 'rc')
SUF = os.environ.get("SUF", "")         # output filename suffix (e.g. '_rr' for the thesis)
IDX = ["rc","synonyms","antonyms","even_odd","rpr","person_total","d2","longstring","irv"]
LAB = {"rc":RC,"synonyms":"synonyms","antonyms":"antonyms","even_odd":"even-odd","rpr":"RPR",
       "person_total":"person-total","d2":"Mahalanobis D²","longstring":"LongString","irv":"IRV"}
FAM = {"rc":"#8f3535","synonyms":"#b45a5a","antonyms":"#b45a5a","even_odd":"#b45a5a","rpr":"#b45a5a",
       "person_total":"#2f5c8f","d2":"#5a7fa8","longstring":"#2f6b3a","irv":"#4a8a55"}

P = pd.read_csv("index_map_pooled_spearman.csv", index_col=0).loc[IDX, IDX]
Ld = pd.read_csv("index_map_loadings.csv", index_col=0)
F1, F2 = Ld.columns[0], Ld.columns[1]  # pair-consistency, normativity

# ================= FIGURE 1 =================
fig, ax = plt.subplots(1, 2, figsize=(13, 5.4))
Pm = P.values; m = len(IDX)
im = ax[0].imshow(Pm, vmin=-1, vmax=1, cmap="RdBu_r")
ax[0].set_xticks(range(m)); ax[0].set_yticks(range(m))
ax[0].set_xticklabels([LAB[c] for c in IDX], rotation=42, ha="right", fontsize=9)
ax[0].set_yticklabels([LAB[c] for c in IDX], fontsize=9)
for i in range(m):
    for j in range(m):
        ax[0].text(j, i, f"{Pm[i,j]:.2f}", ha="center", va="center",
                   color="white" if abs(Pm[i,j])>.55 else "black", fontsize=7.5)
# box rc (0) and person_total (5)
for k in (0, 5):
    ax[0].add_patch(plt.Rectangle((k-.5,-.5), 1, m, fill=False, ec="black", lw=1.8))
    ax[0].add_patch(plt.Rectangle((-.5,k-.5), m, 1, fill=False, ec="black", lw=1.8))
ax[0].set_title("(a) Pooled correlation of careless indices\n(Spearman, 15 datasets)  —  %s ↔ PT = %.2f" % (RC, P.loc["rc","person_total"]), fontsize=11)
fig.colorbar(im, ax=ax[0], fraction=0.046, pad=0.04)

for c in IDX:
    x, y = Ld.loc[c, F1], Ld.loc[c, F2]
    ax[1].scatter(x, y, s=170 if c in ("rc","person_total") else 80, color=FAM[c],
                  edgecolor="black" if c in ("rc","person_total") else "none", linewidth=1.8, zorder=3)
    ax[1].annotate(LAB[c], (x, y), fontsize=10 if c in ("rc","person_total") else 8.5,
                   fontweight="bold" if c in ("rc","person_total") else "normal",
                   xytext=(6,5), textcoords="offset points")
ax[1].axhline(0, color="#ccc", lw=.8); ax[1].axvline(0, color="#ccc", lw=.8)
ax[1].set_xlabel("Factor 1: pair-consistency", fontsize=11)
ax[1].set_ylabel("Factor 2: normativity", fontsize=11)
ax[1].set_title("(b) Index space (varimax)\nred = consistency · blue = normativity/outlier · green = response-style", fontsize=10)
ax[1].grid(alpha=.25)
plt.tight_layout()
plt.savefig("fig_rc_pt_distinct%s.png" % SUF, dpi=200); plt.savefig("C:/Users/vitto/Downloads/fig_rc_pt_distinct%s.png" % SUF, dpi=200)

# ================= FIGURE 2 =================
gt = pd.read_csv("index_map_gt.csv")
FAC = [c for c in gt.columns if c.startswith("F")]
order = ["induced","infrequency","attention_check","overclaiming","speeding_RT"]
CRITLAB = {"induced":"induced\ninconsistency","infrequency":"infrequency /\nbogus items",
           "attention_check":"attention\nchecks","overclaiming":"over-\nclaiming","speeding_RT":"response\nspeed"}
fisher = lambda r: np.arctanh(np.clip(r,-.999,.999))
M = np.full((len(FAC), len(order)), np.nan)
for cj, ct in enumerate(order):
    d = gt[gt.criterion == ct]
    for ri, f in enumerate(FAC):
        v = d[f].values; v = v[np.isfinite(v)]
        if len(v): M[ri, cj] = np.tanh(np.nanmean(fisher(v)))
FACLAB = [c.split(":")[1] if ":" in c else c for c in FAC]

fig2, axg = plt.subplots(figsize=(8.2, 4.6))
im = axg.imshow(M, vmin=-0.8, vmax=0.8, cmap="RdBu_r", aspect="auto")
axg.set_xticks(range(len(order))); axg.set_yticks(range(len(FAC)))
axg.set_xticklabels([CRITLAB[c] for c in order], fontsize=9.5)
axg.set_yticklabels(FACLAB, fontsize=10)
for i in range(len(FAC)):
    for j in range(len(order)):
        if np.isfinite(M[i,j]):
            axg.text(j, i, f"{M[i,j]:+.2f}", ha="center", va="center",
                     color="white" if abs(M[i,j])>.45 else "black",
                     fontsize=11, fontweight="bold" if abs(M[i,j])>.3 else "normal")
# separator: external criteria vs negative controls
axg.axvline(2.5, color="#444", lw=1.6, ls="--")
axg.text(1.0, -0.26, "── external careless criteria ──", ha="center", va="top", transform=axg.get_xaxis_transform(),
         fontsize=9, style="italic", color="#444")
axg.text(3.5, -0.26, "── negative controls ──", ha="center", va="top", transform=axg.get_xaxis_transform(),
         fontsize=9, style="italic", color="#444")
axg.set_title("Which factor tracks which criterion (Spearman)\nnormativity (PT) predicts every benchmark; %s-coherence predicts induced inconsistency" % RC,
              fontsize=10.5, pad=12)
fig2.colorbar(im, ax=axg, fraction=0.046, pad=0.04, label="correlation")
plt.tight_layout()
plt.savefig("fig_rc_pt_validity%s.png" % SUF, dpi=200); plt.savefig("C:/Users/vitto/Downloads/fig_rc_pt_validity%s.png" % SUF, dpi=200)
print("wrote fig_rc_pt_distinct.png and fig_rc_pt_validity.png")
print("\nfactor x criterion (pooled):")
print(pd.DataFrame(M, index=FACLAB, columns=order).round(2).to_string())
