"""Nomological map of rc (9 response-derived indices + RT as supplementary variable).
Pooled Spearman (10x10 incl. response_time), varimax factor structure on the 9 indices,
RT projected into factor space via its correlation with per-respondent factor scores.
3D loading map + pairwise 2D projections. Reads index_map_long.csv."""
import numpy as np, pandas as pd, itertools
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D  # noqa
plt.rcParams.update({"font.family": "serif", "font.size": 10})

IDX = ["rc","synonyms","antonyms","even_odd","rpr","person_total","d2","longstring","irv"]  # factor indices
ALL = IDX + ["response_time"]
LAB = {"rc":"rc","synonyms":"synonyms","antonyms":"antonyms","even_odd":"even-odd","rpr":"RPR",
       "person_total":"person-total","d2":"Mahalanobis D²","longstring":"LongString","irv":"IRV","response_time":"response time"}
FAM = {"rc":"#8f3535","synonyms":"#b45a5a","antonyms":"#b45a5a","even_odd":"#b45a5a","rpr":"#b45a5a",
       "person_total":"#33475f","d2":"#33475f","longstring":"#2f6b3a","irv":"#2f6b3a","response_time":"#c8862a"}
df = pd.read_csv("index_map_long.csv"); df[ALL] = df[ALL].apply(pd.to_numeric, errors="coerce")
datasets = list(dict.fromkeys(df.dataset)); mA = len(ALL); m = len(IDX)
fisher = lambda r: np.arctanh(np.clip(r, -0.999, 0.999)); ifisher = np.tanh

# pooled Spearman (10x10, pairwise Fisher-z)
zsum = np.zeros((mA, mA)); zcnt = np.zeros((mA, mA)); per = {}
for ds in datasets:
    C = df[df.dataset == ds][ALL].corr(method="spearman").values; per[ds] = C
    for a in range(mA):
        for b in range(mA):
            if np.isfinite(C[a, b]): zsum[a, b] += fisher(C[a, b]); zcnt[a, b] += 1
pooled = ifisher(np.divide(zsum, zcnt, out=np.zeros_like(zsum), where=zcnt > 0)); np.fill_diagonal(pooled, 1.0)
P = pd.DataFrame(pooled, index=ALL, columns=ALL); P.to_csv("index_map_pooled_spearman.csv")
print("=== Pooled Spearman (10 indices, %d datasets; RT on %d) ===" % (len(datasets), int(zcnt[ALL.index('response_time'),0])))
print(P.round(2).to_string())
print("\nrc profile:")
for k, v in P["rc"].drop("rc").sort_values(key=lambda s: -s.abs()).items(): print(f"  rc – {LAB[k]:<14} {v:+.2f}")

# PCA / varimax on the 9 response-derived indices
P9 = P.loc[IDX, IDX].values
w, V = np.linalg.eigh(P9); o = np.argsort(w)[::-1]; w = w[o]; V = V[:, o]
load = V * np.sqrt(np.clip(w, 0, None))
def varimax(L, it=200, tol=1e-7):
    L = L.copy(); d = 0
    for _ in range(it):
        d0 = d; u, s, vt = np.linalg.svd(L.T @ (L**3 - L @ np.diag((L**2).sum(0)) / L.shape[0]))
        L = L @ (u @ vt); d = s.sum()
        if d0 and abs(d - d0) / d < tol: break
    return L
ev = pd.DataFrame({"eigenvalue": w[:m], "prop_var": w[:m]/m, "cum": np.cumsum(w[:m])/m})
k = int((w[:m] >= 1.0).sum()); print("\nKaiser factors: %d" % k); print(ev.round(3).head(6).to_string())
Lr = varimax(load[:, :k].copy())
for f in range(k):
    if Lr[np.argmax(np.abs(Lr[:, f])), f] < 0: Lr[:, f] *= -1
FAMLAB = {"person_total":"normativity","rc":"pair-consistency","synonyms":"pair-consistency","antonyms":"pair-consistency",
          "even_odd":"split-half","rpr":"split-half","d2":"outlier/variance","irv":"outlier/variance","longstring":"response-style"}
labels = [FAMLAB[IDX[int(np.argmax(np.abs(Lr[:, f])))]] for f in range(k)]
prio = {"pair-consistency":0,"normativity":1,"outlier/variance":2,"split-half":3,"response-style":4}
ordf = sorted(range(k), key=lambda f: (prio.get(labels[f],9), -(Lr[:,f]**2).sum()))
Lr = Lr[:, ordf]; labels = [labels[f] for f in ordf]
Ld = pd.DataFrame(Lr, index=IDX, columns=[f"F{f+1}:{labels[f]}" for f in range(k)]); Ld.to_csv("index_map_loadings.csv")
print("\n=== Varimax loadings (%d factors, %.0f%% var) ===" % (k, 100*ev.cum.iloc[k-1])); print(Ld.round(2).to_string())

# project RT into factor space via per-respondent factor scores
from scipy.stats import spearmanr
zf = np.zeros(k); zc = np.zeros(k)
for ds in datasets:
    sub = df[df.dataset == ds]
    if sub["response_time"].notna().sum() < 20: continue
    Z = sub[IDX].copy(); Z = (Z - Z.mean())/Z.std(ddof=0); Z = Z.fillna(0.0)
    F = Z.values @ Lr; rt = sub["response_time"].values
    for f in range(k):
        ok = np.isfinite(F[:, f]) & np.isfinite(rt)
        if ok.sum() > 20:
            rho = spearmanr(F[ok, f], rt[ok]).correlation
            if np.isfinite(rho): zf[f] += fisher(rho); zc[f] += 1
rt_coord = ifisher(np.divide(zf, zc, out=np.zeros_like(zf), where=zc > 0))
print("\nRT projected onto factors (corr with factor scores):", {labels[f]: round(rt_coord[f],2) for f in range(k)})

# ---------- FIGURE 1: heatmap (10x10) + 3D ----------
prim = [i for i,l in enumerate(labels) if l in ("pair-consistency","normativity","outlier/variance")][:3]
while len(prim) < 3: prim.append([i for i in range(k) if i not in prim][0])
fig = plt.figure(figsize=(14.5, 6.0))
ax0 = fig.add_subplot(1, 2, 1); Pm = P.loc[ALL, ALL].values
im = ax0.imshow(Pm, vmin=-1, vmax=1, cmap="RdBu_r")
ax0.set_xticks(range(mA)); ax0.set_yticks(range(mA))
ax0.set_xticklabels([LAB[c] for c in ALL], rotation=42, ha="right", fontsize=7.5)
ax0.set_yticklabels([LAB[c] for c in ALL], fontsize=7.5)
for a in range(mA):
    for b in range(mA):
        val = Pm[a,b]; txt = f"{val:.2f}" if np.isfinite(val) else "–"
        ax0.text(b, a, txt, ha="center", va="center", color="white" if abs(val)>0.55 else "black", fontsize=6.3)
ax0.set_title("(a) Pooled Spearman, 9 indices + response time\n(%d datasets; RT on 3)" % len(datasets), fontsize=9.5)
fig.colorbar(im, ax=ax0, fraction=0.046, pad=0.04)
ax1 = fig.add_subplot(1, 2, 2, projection="3d")
xs, ys, zs = Lr[:, prim[0]], Lr[:, prim[1]], Lr[:, prim[2]]
for i, c in enumerate(IDX):
    ax1.scatter(xs[i], ys[i], zs[i], s=120 if c=="rc" else 55, color=FAM[c], edgecolor="black" if c=="rc" else "none", linewidth=1.6, depthshade=True)
    ax1.text(xs[i], ys[i], zs[i]+0.05, LAB[c], fontsize=8, fontweight="bold" if c=="rc" else "normal")
ax1.scatter(rt_coord[prim[0]], rt_coord[prim[1]], rt_coord[prim[2]], s=95, marker="D",
            facecolor="none", edgecolor=FAM["response_time"], linewidth=2, depthshade=True)
ax1.text(rt_coord[prim[0]], rt_coord[prim[1]], rt_coord[prim[2]]+0.05, "response time", fontsize=8, color=FAM["response_time"], style="italic")
ax1.set_xlabel(labels[prim[0]], fontsize=8.5); ax1.set_ylabel(labels[prim[1]], fontsize=8.5); ax1.set_zlabel(labels[prim[2]], fontsize=8.5)
ax1.set_title("(b) Index space — 3 primary factors (varimax)\nred=consistency navy=normativity/outlier green=response-style\namber ◇ = response time (supplementary) · 4th factor=split-half", fontsize=8)
ax1.view_init(elev=20, azim=-56)
plt.tight_layout(); plt.savefig("index_map_fig.png", dpi=150); plt.savefig("C:/Users/vitto/Downloads/rerere_index_map.png", dpi=150)

# ---------- FIGURE 2: pairwise 2D projections with RT supplementary ----------
pairs = list(itertools.combinations(range(k), 2)); ncol = 3; nrow = int(np.ceil(len(pairs)/ncol))
fig2, axs = plt.subplots(nrow, ncol, figsize=(5*ncol, 4.4*nrow)); axs = np.array(axs).reshape(-1)
for ax, (i, j) in zip(axs, pairs):
    for t, c in enumerate(IDX):
        x, y = Lr[t, i], Lr[t, j]
        ax.scatter(x, y, s=120 if c=="rc" else 55, color=FAM[c], edgecolor="black" if c=="rc" else "none", linewidth=1.5, zorder=3)
        ax.annotate(LAB[c], (x, y), fontsize=8, xytext=(5,4), textcoords="offset points", fontweight="bold" if c=="rc" else "normal")
    ax.scatter(rt_coord[i], rt_coord[j], s=95, marker="D", facecolor="none", edgecolor=FAM["response_time"], linewidth=2, zorder=4)
    ax.annotate("resp. time", (rt_coord[i], rt_coord[j]), fontsize=8, color=FAM["response_time"], style="italic", xytext=(5,4), textcoords="offset points")
    ax.axhline(0, color="#ccc", lw=.8); ax.axvline(0, color="#ccc", lw=.8)
    ax.set_xlabel(f"F{i+1}: {labels[i]}", fontsize=9); ax.set_ylabel(f"F{j+1}: {labels[j]}", fontsize=9); ax.grid(alpha=.25)
for ax in axs[len(pairs):]: ax.axis("off")
fig2.suptitle("Index space — pairwise 2D projections (varimax, %d factors) + response time (amber ◇, supplementary)" % k, fontsize=10)
plt.tight_layout(); plt.savefig("index_map_fig2d.png", dpi=150); plt.savefig("C:/Users/vitto/Downloads/rerere_index_map_2d.png", dpi=150)
print("\nwrote index_map_fig.png and index_map_fig2d.png")

print("\n=== rc vs neighbours, per dataset (Spearman) — real-key even-odd/RPR where present ===")
print("dataset          rc~syn rc~EO  rc~RPR rc~PT  rc~RT")
for ds in datasets:
    C = pd.DataFrame(per[ds], index=ALL, columns=ALL); g = lambda a,b: C.loc[a,b]
    rt = g('rc','response_time'); rts = f"{rt:+.2f}" if np.isfinite(rt) else "  –"
    print(f"{ds:<16} {g('rc','synonyms'):+.2f}  {g('rc','even_odd'):+.2f}  {g('rc','rpr'):+.2f}  {g('rc','person_total'):+.2f}  {rts}")
