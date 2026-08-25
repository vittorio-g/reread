import pandas as pd, numpy as np
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
plt.rcParams.update({"font.family":"serif","font.size":11})
OUT="C:/Users/vitto/Downloads/rerere_overleaf/figures/fig_bench_sim.png"
d=pd.read_csv("sim_bench/sim_bench_meanauc.csv").sort_values("items")
x=d["items"].values
series=[("ReReRe","#8f3535","o","-",2.6,"CaReReRe (ensemble)"),
        ("rr","#33475f","o","-",2.6,"rc (index)"),
        ("PersRel","#2f6b3a","s","--",1.8,"Resampled pers. reliability"),
        ("PsychSyn","#6a8caf","^","--",1.5,"Psych. synonyms"),
        ("EvenOdd","#b08a3a","v","--",1.5,"Even-odd"),
        ("IRV","#9a9a9a","D","--",1.5,"IRV"),
        ("LongString","#bdbdbd","P",":",1.4,"LongString"),
        ("Mahalanobis","#c07a7a","x","--",1.6,"Mahalanobis $D^2$")]
fig,ax=plt.subplots(figsize=(7.6,5.0))
for col,c,mk,ls,lw,lab in series:
    if col not in d.columns: continue
    ax.plot(x,d[col],marker=mk,ls=ls,lw=lw,color=c,ms=6,label=lab)
ax.set_xlabel("Questionnaire length (items)"); ax.set_ylabel("AUC vs injected careless label")
ax.set_xticks(x); ax.set_ylim(0.45,1.01); ax.grid(alpha=.3)
ax.set_title("Study 3: competitive AUC vs questionnaire length (20% careless)",fontsize=11.5)
ax.legend(loc="lower right",fontsize=8.6,ncol=2,framealpha=.92)
plt.tight_layout(); plt.savefig(OUT,dpi=150); plt.close()
print("wrote "+OUT)
