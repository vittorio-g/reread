import pandas as pd, numpy as np, collections
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
plt.rcParams.update({"font.family":"serif","font.size":10})
OUT="C:/Users/vitto/Downloads/rerere_overleaf/figures/fig_robustness.png"
NAVY="#33475f"; RUST="#8f3535"; GREEN="#2f6b3a"; GREY="#8a8a8a"; AMBER="#8a6d1f"
fig,ax=plt.subplots(1,3,figsize=(12,3.9))

# Panel A: AUC vs cor_prop (3 datasets)
h=pd.read_csv("sens_hyper.csv"); cp=h[h.param=="cor_prop"]
cols={"study":NAVY,"smarvus":GREEN,"warning":RUST}
for ds,c in cols.items():
    d=cp[cp.dataset==ds].sort_values("value")
    ax[0].plot(d.value*100, d.auc, marker="o", color=c, label=ds)
ax[0].axvline(3,color=AMBER,ls=":",lw=1)
ax[0].set_xlabel("top-share cor_prop (%)"); ax[0].set_ylabel("Ensemble AUC")
ax[0].set_title("(a) Robust to the pair top-share\n(floor 15 & gate constants: exactly flat)",fontsize=10)
ax[0].legend(fontsize=8,loc="lower right"); ax[0].grid(alpha=.3); ax[0].set_ylim(0.65,1.0)

# Panel B: pair contamination + rr AUC vs prevalence
p=pd.read_csv("sens_pairs.csv")
ag=p.groupby("prevalence").agg(sf=("samefactor_frac_top3","mean"),
    au=("rr_auc",lambda s: np.nanmean(pd.to_numeric(s,errors="coerce")))).reset_index()
ax[1].plot(ag.prevalence*100, ag.sf, marker="s", color=NAVY, label="same-factor fraction of top-3%")
ax[1].plot(ag.prevalence*100, ag.au, marker="o", color=RUST, label="rc AUC")
ax[1].set_xlabel("careless prevalence (%)"); ax[1].set_ylabel("fraction / AUC")
ax[1].set_title("(b) Pair-selection contaminates,\nbut rc AUC holds",fontsize=10)
ax[1].legend(fontsize=8,loc="lower left"); ax[1].grid(alpha=.3); ax[1].set_ylim(0.55,1.02)

# Panel C: one-sided pi-hat calibration (sim overlap) with real-study reference
pi=pd.read_csv("sens_pi.csv")
cg=pi.groupby("true_prev").pi_hat.mean().reset_index()
ax[2].plot([0,32],[0,32],ls="--",color=GREY,label="identity")
ax[2].plot(cg.true_prev*100, cg.pi_hat*100, marker="o", color=NAVY, label="sim (strong overlap)")
# real-study calibration points: read from the (regenerated) Fig-3 data, rates <= 30%
try:
    rc=pd.read_csv("fig_calib_pi.csv"); rc=rc[rc["true"]<=0.31]
    ax[2].plot(100*rc["true"], 100*rc["pi_mean"], marker="^", color=GREEN, label="real study (Fig. 3)")
except Exception: pass
ax[2].axvspan(12,32,color=NAVY,alpha=.05)
ax[2].set_xlabel("true careless rate (%)"); ax[2].set_ylabel(r"estimated $\hat{\pi}$ (%)")
ax[2].set_title("(c) one-sided $\\hat{\\pi}$: tracks the truth even\nunder overlap; conservative at low rates",fontsize=10)
ax[2].legend(fontsize=8,loc="upper left"); ax[2].grid(alpha=.3); ax[2].set_xlim(0,32); ax[2].set_ylim(0,32)

plt.tight_layout(); plt.savefig(OUT,dpi=150); plt.close(); print("wrote "+OUT)
