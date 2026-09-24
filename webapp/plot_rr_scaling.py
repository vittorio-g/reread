import os
import pandas as pd, numpy as np
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
plt.rcParams.update({"font.family":"serif","font.size":11})
OUT="C:/Users/vitto/AppData/Local/Temp/claude/C--Users-vitto-Desktop-ReReReRe/8f426989-5c28-4009-a838-215f9a1e2ca0/scratchpad/tesi/figures/fig_rr_scaling.png"
NAVY="#33475f"; GREY="#8a8a8a"; RUST="#8f3535"
RC=os.environ.get("RCLABEL","rc")   # display label; the CSV column stays auc_rr
d=pd.read_csv("sim_rr_scaling.csv").sort_values("items")
x=d["items"].values
fig,ax=plt.subplots(figsize=(7.2,4.6))
ax.axhline(0.5,color="#ccc",lw=1,ls="--",zorder=0)
ax.fill_between(x,d.auc_ens-d.auc_ens_se,d.auc_ens+d.auc_ens_se,color=GREY,alpha=.12,zorder=1)
ax.plot(x,d.auc_ens,marker="s",lw=1.8,color=GREY,zorder=2,label=f"ensemble ({RC} + partners), for reference")
ax.fill_between(x,d.auc_rr-d.auc_rr_se,d.auc_rr+d.auc_rr_se,color=NAVY,alpha=.16,zorder=2)
ax.plot(x,d.auc_rr,marker="o",lw=2.2,color=NAVY,zorder=3,label=f"${RC}$ index alone")
for xi,yi in zip(x,d.auc_rr): ax.annotate(f"{yi:.2f}",(xi,yi),textcoords="offset points",xytext=(0,-14),ha="center",fontsize=8,color=NAVY)
ax.set_xlabel("Questionnaire length (items)"); ax.set_ylabel("AUC vs.\\ careless label")
ax.set_xticks(x); ax.set_ylim(0.5,1.02); ax.grid(alpha=.3)
ax.legend(loc="lower right",fontsize=9.5)
ax.set_title(f"The ${RC}$ index scales with questionnaire length",fontsize=11.5)
plt.tight_layout(); plt.savefig(OUT,dpi=150); plt.close()
print("wrote "+OUT)
