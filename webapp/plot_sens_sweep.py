import pandas as pd, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
plt.rcParams.update({"font.family":"serif","font.size":11})
d=pd.read_csv("eval_sens_sweep.csv"); x=d.true_rate*100
fig,ax=plt.subplots(figsize=(8,5.2))
ax.axvspan(10,25,color="#2f6b3a",alpha=.08,label="green range (10–25%)")
ax.plot(x,d.low,   "o-",color="#33475f",lw=2.1,ms=6,label="Low sensitivity")
ax.plot(x,d.medium,"s-",color="#2f6b3a",lw=2.1,ms=6,label="Medium (default)")
ax.plot(x,d["high"],"^-",color="#8f3535",lw=2.1,ms=6,label="High sensitivity")
# region annotations for who wins where
ax.annotate("Low wins\n(clean data)",(8,0.90),ha="center",fontsize=9,color="#33475f",style="italic")
ax.annotate("Medium wins",(17.5,0.905),ha="center",fontsize=9,color="#2f6b3a",style="italic")
ax.annotate("High wins\n(heavy contamination)",(41,0.83),ha="center",fontsize=9,color="#8f3535",style="italic")
ax.set_xlabel("True careless rate (%)"); ax.set_ylabel("MCC (automatic flags)")
ax.set_title("Which sensitivity setting to use depends on how contaminated the data are\n(real careful + careless mixes, 30 resamples/point)",fontsize=11)
ax.set_xticks(range(5,55,5)); ax.set_ylim(0.25,0.95); ax.grid(alpha=.25)
ax.legend(loc="lower left",fontsize=9.5,framealpha=.92)
plt.tight_layout()
for p in ["fig_sensitivity_sweep.png", r"C:/Users/vitto/Downloads/ReReRe_sensitivity_sweep.png"]:
    plt.savefig(p,dpi=160)
print("wrote fig_sensitivity_sweep.png")
