"""Close-up: index correlations within the CARELESS of Study 1 (induced, n=87, clean labels).
Study 1 careless carry a severity gradient (corruption 50/75/100). Raw within-careless
correlation mixes (a) shared severity gradient + (b) residual coupling. We separate them with
a partial correlation controlling severity, and correlate each index with severity."""
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
plt.rcParams.update({"font.family": "serif", "font.size": 10})

IDX = ["rc","synonyms","antonyms","even_odd","rpr","person_total","d2","longstring","irv"]
LAB = {"rc":"rc","synonyms":"synonyms","antonyms":"antonyms","even_odd":"even-odd","rpr":"RPR",
       "person_total":"person-total","d2":"D²","longstring":"LongString","irv":"IRV"}
long = pd.read_csv("index_map_long.csv"); long[IDX] = long[IDX].apply(pd.to_numeric, errors="coerce")
L = pd.read_csv("_study_labels.csv")
sub = long[long.dataset == "Study1_induced"].reset_index(drop=True)
assert len(sub) == len(L), f"{len(sub)} vs {len(L)}"
sub = sub.copy(); sub["careless"] = L["careless"].values; sub["group"] = L["group"].values
sev_map = {"careless_50": 50, "careless_75": 75, "careless_100": 100}
cl = sub[sub.careless == 1].copy(); cf = sub[sub.careless == 0].copy()
cl["sev"] = cl["group"].map(sev_map)
m = len(IDX)
def offmean(C): return np.nanmean(np.abs(C[~np.eye(m, dtype=bool)]))

# rank-transform for Spearman-based partials
def rankz(df, cols): return df[cols].rank()
R_cl = rankz(cl, IDX + []); R_cl["sev"] = cl["sev"].rank().values
Craw = cl[IDX].corr(method="spearman").values
Ccf = cf[IDX].corr(method="spearman").values
# each index ~ severity
sev_corr = {c: cl[c].corr(cl["sev"], method="spearman") for c in IDX}
# partial correlation controlling severity (on ranks)
rr = R_cl[IDX].corr().values  # pearson on ranks == spearman
rs = np.array([R_cl[c].corr(R_cl["sev"]) for c in IDX])  # index~sev on ranks
Cpar = np.zeros((m, m))
for a in range(m):
    for b in range(m):
        num = rr[a, b] - rs[a] * rs[b]
        den = np.sqrt(max(1e-9, (1 - rs[a]**2) * (1 - rs[b]**2)))
        Cpar[a, b] = num / den
np.fill_diagonal(Cpar, 1)

print("=== Study 1 careless (n=%d) ===" % len(cl))
print("mean|off-diag r|:  careful %.3f | careless RAW %.3f | careless PARTIAL(sev) %.3f" %
      (offmean(Ccf), offmean(Craw), offmean(Cpar)))
print("\nindex ~ severity (Spearman, within careless):")
for c in IDX: print(f"  {LAB[c]:<12} {sev_corr[c]:+.2f}")
print("\nrc row: careful | careless RAW | careless PARTIAL(sev)")
for c in IDX:
    if c == "rc": continue
    j = IDX.index(c)
    print(f"  rc~{LAB[c]:<12} {Ccf[0,j]:+.2f}   {Craw[0,j]:+.2f}   {Cpar[0,j]:+.2f}")

# figure: 3 heatmaps + severity bars
fig, ax = plt.subplots(1, 4, figsize=(19, 4.8), gridspec_kw={"width_ratios":[1,1,1,0.6]})
for a, (Pm, ttl) in zip(ax[:3], [(Ccf,"CAREFUL (n=%d)"%len(cf)),(Craw,"CARELESS raw (n=%d)"%len(cl)),(Cpar,"CARELESS | severity controlled")]):
    im = a.imshow(Pm, vmin=-1, vmax=1, cmap="RdBu_r")
    a.set_xticks(range(m)); a.set_yticks(range(m))
    a.set_xticklabels([LAB[c] for c in IDX], rotation=45, ha="right", fontsize=7.5); a.set_yticklabels([LAB[c] for c in IDX], fontsize=7.5)
    for i in range(m):
        for j in range(m):
            a.text(j,i,f"{Pm[i,j]:.2f}",ha="center",va="center",color="white" if abs(Pm[i,j])>.55 else "black",fontsize=6.3)
    a.set_title(f"{ttl}\nmean|off-diag|={offmean(Pm):.2f}", fontsize=9)
fig.colorbar(im, ax=ax[2], fraction=0.046, pad=0.04)
sv = [sev_corr[c] for c in IDX]
ax[3].barh(range(m), sv, color=["#8f3535" if c in("rc","synonyms","antonyms","even_odd","rpr") else "#33475f" if c in("person_total","d2") else "#2f6b3a" for c in IDX])
ax[3].set_yticks(range(m)); ax[3].set_yticklabels([LAB[c] for c in IDX], fontsize=7.5); ax[3].invert_yaxis()
ax[3].axvline(0,color="#888",lw=.8); ax[3].set_xlabel("corr with severity"); ax[3].set_title("index ~ corruption\nlevel (50/75/100)", fontsize=9); ax[3].grid(alpha=.25,axis="x")
plt.tight_layout(); plt.savefig("index_study1_careless_fig.png", dpi=150); plt.savefig("C:/Users/vitto/Downloads/rerere_study1_careless.png", dpi=150)
print("\nwrote index_study1_careless_fig.png")
