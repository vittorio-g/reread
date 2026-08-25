# Forest plot: incremental contribution of rr (ΔoracleMCC of {Long+PT+rr} vs {Long+PT})
# on the current authoritative ext_bench data. Style matches fig_benchmark.py.
import pandas as pd, numpy as np
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
plt.rcParams.update({"font.family":"serif","font.size":11})
RUST="#8f3535"; GREY="#9a9a9a"; ORANGE="#c07a1e"
d=pd.read_csv("ext_bench/rr_incremental_current.csv").set_index("dataset")

# explicit top->bottom order: envelope (J>=60, by J desc) then boundary (J<60, by J desc)
env=["warning300","opsy_16pf","smarvus","mastroianni22"]
bnd=["duckworth","douglas","krause","ivanov21","pennycook20","ogrady19","moss23","alvarez19"]
order=env+bnd
d=d.loc[order]
N=len(order)
ypos=list(range(N))                      # 0..N-1 top->bottom after invert
labels=[f"{ds}   J={int(d.loc[ds,'J'])}" for ds in order]

fig,ax=plt.subplots(figsize=(8.0,5.6))
ax.axvline(0,color="#888",lw=1,ls="--",zorder=1)
for i,ds in enumerate(order):
    r=d.loc[ds]; sig=r.P_dMCC_gt0>=0.975; art=(ds=="moss23")
    c=ORANGE if art else (RUST if sig else GREY)
    ax.plot([r.dMCC_lo,r.dMCC_hi],[i,i],color=c,lw=2.3,zorder=2)
    ax.plot(r.dMCC,i,"o",color=c,ms=8.5 if (sig or art) else 6,
            mfc=c if (sig or art) else "white",mec=c,zorder=3)
ax.set_yticks(ypos); ax.set_yticklabels(labels,fontsize=9)
ax.invert_yaxis()
# divider + group headers (in empty right area)
ax.axhline(len(env)-0.5,color="#ccc",lw=1,zorder=1)
ax.text(0.285,(len(env)-1)/2,"envelope\n(J≥60, multi-item)",fontsize=8.5,style="italic",
        color="#666",ha="right",va="center")
ax.text(0.285,len(env)+ (len(bnd)-1)/2,"boundary\n(short / few-item)",fontsize=8.5,style="italic",
        color="#666",ha="right",va="center")
# moss23 artifact annotation
ai=order.index("moss23")
ax.annotate("spurious: J=8, rr inverted (AUC .29)\n→ not mechanism evidence",
            xy=(d.loc['moss23','dMCC_lo'],ai),xytext=(-0.075,ai-0.55),fontsize=8,color=ORANGE,
            ha="left",arrowprops=dict(arrowstyle="->",color=ORANGE,lw=1))
ax.set_xlim(-0.08,0.30)
ax.set_xlabel(r"$\Delta$ oracle-MCC   (rr added to LongString + Person-Total; 95% bootstrap CI)")
ax.set_title("Incremental contribution of the rr index on external data",fontsize=12.5)
from matplotlib.lines import Line2D
ax.legend(handles=[Line2D([0],[0],marker="o",color=RUST,lw=2.3,label="CI excludes 0  (P>.975)"),
                   Line2D([0],[0],marker="o",color=GREY,mfc="white",lw=2.3,label="CI includes 0"),
                   Line2D([0],[0],marker="o",color=ORANGE,lw=2.3,label="artifact (out of envelope)")],
          loc="lower right",fontsize=8.5,framealpha=.92)
ax.grid(axis="x",alpha=.25)
plt.tight_layout()
OUT="ext_bench/fig_rr_incremental.png"
plt.savefig(OUT,dpi=160); plt.close()
print("wrote "+OUT)
