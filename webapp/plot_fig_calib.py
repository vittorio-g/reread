import pandas as pd, numpy as np
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
plt.rcParams.update({"font.family":"serif","font.size":11})
OUTS=["C:/Users/vitto/Downloads/rerere_overleaf/figures/fig_rate_calibration_fine.png",
      "fig_rate_calibration_fine.png"]
NAVY="#33475f"; GREEN="#2f6b3a"; GREY="#888"; AMBER="#8a6d1f"
d=pd.read_csv("fig_calib_pi.csv")
x=100*d["true"].values
m=100*d["pi_mean"].values; lo=100*d["pi_p5"].values; hi=100*d["pi_p95"].values
fm=100*d["flag_mean"].values; flo=100*d["flag_p5"].values; fhi=100*d["flag_p95"].values
fig,ax=plt.subplots(figsize=(7.2,4.9))
# the whole tested range is calibrated under the one-sided estimator
ax.axvspan(x.min(),x.max(),color=NAVY,alpha=.04,zorder=0)
lim=max(x.max(),hi.max())+3
ax.plot([0,lim],[0,lim],ls="--",color=GREY,lw=1.3,label="Identity (perfect calibration)",zorder=1)
# (a) prevalence estimate: how many careless are PRESENT
ax.fill_between(x,lo,hi,color=NAVY,alpha=.16,zorder=2,label="90% coverage band ($\\hat{\\pi}$, 5–95th pct)")
ax.plot(x,m,marker="o",lw=2,color=NAVY,zorder=4,label="Estimated prevalence $\\hat{\\pi}$ (mean)")
# (b) flagged share: how many are actually REMOVED by the shipped Standard cut
ax.plot(x,flo,ls=(0,(3,2)),lw=1.1,color=GREEN,alpha=.75,zorder=3)
ax.plot(x,fhi,ls=(0,(3,2)),lw=1.1,color=GREEN,alpha=.75,zorder=3)
ax.plot(x,fm,marker="s",ms=5,lw=2.2,color=GREEN,zorder=5,
        label="Flagged share, Standard cut (mean; dashed 5–95th pct)")
bias=float(np.max(np.abs(m-x)))
ax.text(lim-0.5,1.5,"$\\hat{\\pi}$ tracks the truth throughout\n(max $|$bias$| = %.1f$ pp over 5--45%%)" % bias,
        color=AMBER,fontsize=8.5,va="bottom",ha="right")
ax.set_xlim(0,lim); ax.set_ylim(0,lim)
ax.set_xlabel("True careless rate (%)")
ax.set_ylabel("Share of the sample (%)")
ax.set_title("What is present vs what is removed:\nprevalence estimate $\\hat{\\pi}$ and flagged share vs true rate",fontsize=11.5)
ax.grid(alpha=.3); ax.legend(loc="upper left",fontsize=8.8,framealpha=.9)
plt.tight_layout()
for p in OUTS: plt.savefig(p,dpi=150)
plt.close()
print("wrote "+" and ".join(OUTS))
