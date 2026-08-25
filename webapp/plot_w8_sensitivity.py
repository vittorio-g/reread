# Point 8 figure: (A) bootstrap CIs on the ensemble coefficients (false precision),
# (B) AUC robustness to coarsening / perturbing / refitting the weights.
import pandas as pd, numpy as np
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
plt.rcParams.update({"font.family":"serif","font.size":11})
NAVY="#33475f"; RUST="#8f3535"; GREY="#9a9a9a"; GREEN="#2f6b3a"
co=pd.read_csv("ext_bench/w8_coef_ci.csv")
co["param"]=co["param"].replace({"rr":"rc"})
va=pd.read_csv("ext_bench/w8_variant_auc.csv")
va["label"]=va["label"].str.replace("rr + equal aux","rc + equal aux",regex=False)
pe=pd.read_csv("ext_bench/w8_perturb.csv")
lo=pd.read_csv("ext_bench/w8_loo.csv")

fig,(axA,axB)=plt.subplots(1,2,figsize=(11,4.6))

# ---- Panel A: coefficient bootstrap CIs ----
co=co[co.param!="b0"].reset_index(drop=True)   # show the 3 slopes
yy=np.arange(len(co))[::-1]
for i,(_,r) in enumerate(co.iterrows()):
    yv=yy[i]
    axA.plot([r.boot_lo,r.boot_hi],[yv,yv],color=NAVY,lw=2.4,alpha=.55,zorder=2)
    axA.plot(r.boot_md,yv,"o",color=NAVY,ms=6,mfc="white",zorder=3)      # bootstrap median
    axA.plot(r.shipped,yv,"D",color=RUST,ms=8,zorder=4)                   # shipped value
axA.set_yticks(yy); axA.set_yticklabels(co.param)
axA.set_xlabel("logistic coefficient (log-odds per SD)")
axA.set_title("(A) The 4-decimal weights are not identifiable",fontsize=11.5)
axA.axvline(0,color="#ccc",lw=1,ls="--")
axA.grid(axis="x",alpha=.25)
from matplotlib.lines import Line2D
axA.legend(handles=[Line2D([0],[0],marker="D",color=RUST,lw=0,label="shipped weight"),
                    Line2D([0],[0],marker="o",color=NAVY,mfc="white",lw=2.4,label="bootstrap median + 95% CI")],
           loc="lower right",fontsize=8.5)

# ---- Panel B: AUC robustness ----
labels=list(va.label); aucs=list(va.AUC)
ypos=np.arange(len(labels))[::-1]
axB.barh(ypos,[a-0.95 for a in aucs],left=0.95,color=[RUST if "shipped" in l else GREY for l in labels],
         height=0.55,zorder=2)
for y,a in zip(ypos,aucs): axB.text(a+0.0006,y,f"{a:.3f}",va="center",fontsize=8.5)
axB.set_yticks(ypos); axB.set_yticklabels(labels,fontsize=8.5)
axB.set_xlim(0.95,1.0)
# perturbation ±50% band + LOO points annotated as text
p50=pe[pe.pct==50].iloc[0]
axB.set_xlabel("AUC (Study 1, N=157)")
axB.set_title("(B) …and the ensemble is insensitive to them",fontsize=11.5)
axB.axvline(va.AUC[va.label.str.contains("shipped")].iloc[0],color=RUST,lw=1,ls=":",alpha=.7,zorder=1)
axB.grid(axis="x",alpha=.25)
note=(f"±50% random weight jitter: AUC {p50['mean']:.3f} [{p50['lo']:.3f},{p50['hi']:.3f}], min {p50['min']:.3f}\n"
      f"LOO fixed {lo.AUC[0]:.3f} vs refit {lo.AUC[1]:.3f} → freezing wins")
axB.text(0.9515,-0.85,note,fontsize=7.6,color="#333",va="top")
axB.set_ylim(-1.6,len(labels)-0.4)

plt.tight_layout()
OUT="ext_bench/fig_w8_weights.png"; plt.savefig(OUT,dpi=160); plt.close(); print("wrote "+OUT)
