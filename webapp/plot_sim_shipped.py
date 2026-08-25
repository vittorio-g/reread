import sys, pandas as pd
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
plt.rcParams.update({"font.family":"serif","font.size":11})
FIGDIR=sys.argv[1] if len(sys.argv)>1 else "."
d=pd.read_csv("sim_shipped.csv").sort_values("items")
x=d["items"].values
NAVY="#33475f"; GREEN="#2f6b3a"; RUST="#8f3535"

# ---- Fig 4: detection scales with length (shipped model, at the shipped cut) ----
fig,ax=plt.subplots(figsize=(7.4,4.6))
ax.errorbar(x,d.om_full,yerr=d.om_full_se,marker="o",lw=2,capsize=4,color=NAVY,label="Severity-weighted MCC (automatic cut)")
ax.errorbar(x,d.auc_full,yerr=d.auc_full_se,marker="s",lw=2,capsize=4,color=GREEN,label="AUC")
ax.set_xlabel("Questionnaire length (items)"); ax.set_ylabel("Detection quality")
ax.set_xticks(x); ax.set_ylim(0.45,1.0); ax.grid(alpha=.3); ax.legend(loc="lower right",fontsize=10)
ax.set_title("Shipped ensemble: detection scales with questionnaire length",fontsize=11.5)
plt.tight_layout(); plt.savefig(f"{FIGDIR}/fig_scale_rate_aware.png",dpi=150); plt.close()
print("wrote fig_scale_rate_aware.png")

# ---- Fig 5: rc's contribution vs length (Delta weighted MCC at the shipped cut, bootstrap 95% CI) ----
fig,ax=plt.subplots(figsize=(7.4,4.6))
ax.axhline(0,color="#999",lw=1,ls="--")
ax.fill_between(x,d.d_om_lo,d.d_om_hi,color=RUST,alpha=.15)
ax.plot(x,d.d_om,marker="o",lw=2,color=RUST,label="Δ weighted MCC at the automatic cut (full − no-rc)")
ax.set_xlabel("Questionnaire length (items)"); ax.set_ylabel("rc contribution (Δ weighted MCC, automatic cut)")
ax.set_xticks(x); ax.grid(alpha=.3); ax.legend(loc="upper left",fontsize=10)
ax.set_title("Contribution of the rc index vs length (bootstrap 95% CI)",fontsize=11.5)
plt.tight_layout(); plt.savefig(f"{FIGDIR}/fig_ablation_curves.png",dpi=150); plt.close()
print("wrote fig_ablation_curves.png")
