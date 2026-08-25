import pandas as pd, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
plt.rcParams.update({"font.family":"serif","font.size":11})
d=pd.read_csv("real_data_test.csv").sort_values("true_rate"); x=d["true_rate"]*100
M={"std":("Standard (z=2.5)","#33475f","o"),"high":("High (z=1.5)","#8f3535","s"),
   "A":("A: pi top-share","#2f6b3a","^"),"B":("B: Bayes","#8a6d1f","D"),
   "C":("C: FDR 10%","#5a5a8f","v"),"D":("D: adaptive z","#c46210","P"),"E":("E: 3-component","#b02f8a","X")}
fig,ax=plt.subplots(figsize=(8.2,5.0))
ax.plot(x,d["oracle"],"--",color="#999",lw=1.6,label="oracle ceiling")
for m,(lab,c,mk) in M.items():
    lw=2.4 if m in("D","E") else 1.6
    ax.plot(x,d[m+"_mcc"],marker=mk,color=c,lw=lw,label=lab)
ax.set_xlabel("Careless rate (%)"); ax.set_ylabel("MCC")
ax.set_ylim(0.2,0.95); ax.grid(alpha=.3); ax.legend(fontsize=8.5,ncol=2,loc="lower left")
ax.set_title("REAL collected study data (n=157, extreme careless): method vs prevalence",fontsize=11.5)
plt.tight_layout(); plt.savefig("v2_fig6_realdata.png",dpi=150); plt.close()
print("wrote v2_fig6_realdata.png")
