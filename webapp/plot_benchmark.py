import pandas as pd, numpy as np
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
plt.rcParams.update({"font.family":"serif","font.size":11})
OUT="C:/Users/vitto/Downloads/rerere_overleaf/figures/fig_benchmark.png"
NAVY="#33475f"; RUST="#8f3535"; GREY="#8a8a8a"; GREEN="#2f6b3a"
d=pd.read_csv("benchmark_study.csv")
d=d[d["AUC"].notna()].copy().sort_values("AUC")   # drop antonyms (NA)
d["method"]=d["method"].str.replace(r"^ReReRe","CaReReRe",regex=True)
d["method"]=d["method"].str.replace("rr (index)","rc (index)",regex=False)
y=np.arange(len(d))
def col(m):
    if m.startswith("CaReReRe"): return RUST
    if m.startswith("rc"): return NAVY
    return GREY
colors=[col(m) for m in d["method"]]
fig,ax=plt.subplots(figsize=(7.6,4.8))
ax.axvline(0.5,color="#bbb",lw=1,ls="--",zorder=0)
for i,(_,r) in enumerate(d.iterrows()):
    ax.plot([r.lo95,r.hi95],[i,i],color=colors[i],lw=2.2,zorder=2)
    ax.plot(r.AUC,i,"o",color=colors[i],ms=7,zorder=3)
    # significance vs rr (DeLong)
    p=r.p_vs_rr
    star = "" if pd.isna(p) else ("***" if p<1e-3 else "**" if p<1e-2 else "*" if p<.05 else "n.s.")
    if r.method not in ("rc (index)",):
        ax.text(r.hi95+0.008,i,star,va="center",ha="left",fontsize=9,color="#444")
ax.set_yticks(y); ax.set_yticklabels(d["method"])
ax.set_xlabel("AUC (95% CI, DeLong)"); ax.set_xlim(0.5,1.03)
ax.set_title("Study 1 head-to-head: CaReReRe vs established careless indices",fontsize=11.5)
# legend
from matplotlib.lines import Line2D
ax.legend(handles=[Line2D([0],[0],color=RUST,lw=3,label="CaReReRe (ensemble)"),
                   Line2D([0],[0],color=NAVY,lw=3,label="rc (index)"),
                   Line2D([0],[0],color=GREY,lw=3,label="established indices (careless pkg)")],
          loc="lower right",fontsize=9)
ax.grid(axis="x",alpha=.3)
plt.tight_layout(); plt.savefig(OUT,dpi=150); plt.close()
print("wrote "+OUT)
