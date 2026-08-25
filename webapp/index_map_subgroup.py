"""Within-group index correlations: are the careless indices uncorrelated in CAREFUL
respondents but correlated in CARELESS ones? Split each GT dataset by its careless label,
compute the 9-index Spearman matrix within each subgroup, pool (Fisher-z), and compare.
Also report the standardized mean-shift (careless - careful) per index."""
import numpy as np, pandas as pd, itertools
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
plt.rcParams.update({"font.family": "serif", "font.size": 10})

IDX = ["rc","synonyms","antonyms","even_odd","rpr","person_total","d2","longstring","irv"]
LAB = {"rc":"rc","synonyms":"synonyms","antonyms":"antonyms","even_odd":"even-odd","rpr":"RPR",
       "person_total":"person-total","d2":"D²","longstring":"LongString","irv":"IRV"}
long = pd.read_csv("index_map_long.csv"); long[IDX] = long[IDX].apply(pd.to_numeric, errors="coerce")
D = "../Dataset/gt_benchmark_candidates"
GT = {  # dataset -> (gt_path, column)  [binary careless label, careless=1]
  "Study1_induced": ("_study_labels.csv", "careless"),
  "Kay_S1": ("bogus_bench/kay_s1_gt.csv", "gt"),
  "Kay_S2": ("bogus_bench/kay_s2_gt.csv", "gt"),
  "Kay_S5": ("bogus_bench/kay_s5_gt.csv", "gt"),
  "warning_IPIP300": ("bogus_bench/warning_gt.csv", "gt"),
  "krause": ("bogus_bench/krause_gt.csv", "gt"),
  "douglas": (D+"/douglas2023_dataquality/douglas_labels.csv", "y_attn_ge2"),
}
m = len(IDX)
fisher = lambda r: np.arctanh(np.clip(r, -0.999, 0.999)); ifisher = np.tanh

def offmean(C):  # mean |off-diagonal|
    return np.nanmean(np.abs(C[~np.eye(m, dtype=bool)]))

zc0 = np.zeros((m,m)); nc0 = np.zeros((m,m)); zc1 = np.zeros((m,m)); nc1 = np.zeros((m,m))
shift = np.zeros(m); shn = 0
print("dataset          n_cf n_cl  mean|r|_careful  mean|r|_careless  (pooled)")
for ds,(gp,col) in GT.items():
    sub = long[long.dataset==ds].reset_index(drop=True)
    g = pd.to_numeric(pd.read_csv(gp)[col], errors="coerce").values
    if len(g)!=len(sub): print(f"  ! {ds}: len {len(g)}!={len(sub)}"); continue
    cf = sub[g==0][IDX]; cl = sub[g==1][IDX]
    if len(cf)<30 or len(cl)<30: print(f"  ! {ds}: too few (cf {len(cf)}, cl {len(cl)})"); continue
    C0 = cf.corr(method="spearman").values; C1 = cl.corr(method="spearman").values
    Cp = sub[IDX].corr(method="spearman").values
    for a in range(m):
        for b in range(m):
            if np.isfinite(C0[a,b]): zc0[a,b]+=fisher(C0[a,b]); nc0[a,b]+=1
            if np.isfinite(C1[a,b]): zc1[a,b]+=fisher(C1[a,b]); nc1[a,b]+=1
    # standardized mean-shift (robust): median diff / pooled MAD
    for j,c in enumerate(IDX):
        x0=cf[c].dropna(); x1=cl[c].dropna()
        s=(x0.std()+x1.std())/2
        if s>0: shift[j]+=(x1.median()-x0.median())/s
    shn+=1
    print(f"{ds:<16} {len(cf):>4} {len(cl):>4}    {offmean(C0):.3f}            {offmean(C1):.3f}         ({offmean(Cp):.3f})")

P0 = ifisher(np.divide(zc0,nc0,out=np.zeros_like(zc0),where=nc0>0)); np.fill_diagonal(P0,1)
P1 = ifisher(np.divide(zc1,nc1,out=np.zeros_like(zc1),where=nc1>0)); np.fill_diagonal(P1,1)
shift/=max(shn,1)
D0=pd.DataFrame(P0,index=IDX,columns=IDX); D1=pd.DataFrame(P1,index=IDX,columns=IDX)
D0.to_csv("index_subgroup_careful.csv"); D1.to_csv("index_subgroup_careless.csv")

print("\n=== POOLED mean |off-diagonal r| ===")
print(f"  within CAREFUL : {offmean(P0):.3f}")
print(f"  within CARELESS: {offmean(P1):.3f}")
print("\n=== rc's correlations: careful vs careless ===")
for c in IDX:
    if c=="rc": continue
    print(f"  rc~{LAB[c]:<12} careful {P0[0,IDX.index(c)]:+.2f}   careless {P1[0,IDX.index(c)]:+.2f}   d= {P1[0,IDX.index(c)]-P0[0,IDX.index(c)]:+.2f}")
print("\n=== standardized mean-shift (careless - careful), per index (SD units) ===")
for j,c in enumerate(IDX): print(f"  {LAB[c]:<12} {shift[j]:+.2f}")

# figure: two heatmaps
fig,ax=plt.subplots(1,3,figsize=(16,5))
for a,(Pm,ttl) in zip(ax[:2],[(P0,"within CAREFUL"),(P1,"within CARELESS")]):
    im=a.imshow(Pm,vmin=-1,vmax=1,cmap="RdBu_r")
    a.set_xticks(range(m)); a.set_yticks(range(m))
    a.set_xticklabels([LAB[c] for c in IDX],rotation=45,ha="right",fontsize=8); a.set_yticklabels([LAB[c] for c in IDX],fontsize=8)
    for i in range(m):
        for j in range(m):
            a.text(j,i,f"{Pm[i,j]:.2f}",ha="center",va="center",color="white" if abs(Pm[i,j])>.55 else "black",fontsize=6.5)
    a.set_title(f"{ttl}\nmean|off-diag r| = {offmean(Pm):.2f}",fontsize=10)
fig.colorbar(im,ax=ax[1],fraction=0.046,pad=0.04)
# bar: mean-shift
ax[2].barh(range(m),shift,color=["#8f3535" if c in("rc","synonyms","antonyms","even_odd","rpr") else "#33475f" if c in("person_total","d2") else "#2f6b3a" for c in IDX])
ax[2].set_yticks(range(m)); ax[2].set_yticklabels([LAB[c] for c in IDX],fontsize=8); ax[2].invert_yaxis()
ax[2].axvline(0,color="#888",lw=.8); ax[2].set_xlabel("careless − careful (SD)"); ax[2].set_title("mean elevation in careless",fontsize=10); ax[2].grid(alpha=.25,axis="x")
plt.tight_layout(); plt.savefig("index_subgroup_fig.png",dpi=150); plt.savefig("C:/Users/vitto/Downloads/rerere_subgroup.png",dpi=150)
print("\nwrote index_subgroup_fig.png")
