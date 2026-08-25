import pandas as pd, numpy as np, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
plt.rcParams.update({"font.family":"serif","font.size":11})
df=pd.read_csv("eval_mix.csv")
t=df.true_rate.values*100; e=df.est_mean.values*100; s=df.est_sd.values*100
fig,ax=plt.subplots(figsize=(7.2,6.4))
ax.axvspan(10,30,color="#2f6b3a",alpha=0.07,label="realistic range (10–30%)")
ax.plot([0,52],[0,52],"--",color="#8a8577",lw=1.4,label="perfect calibration (y = x)")
ax.errorbar(t,e,yerr=s,marker="o",ms=7,lw=2,capsize=4,color="#33475f",
            label="ReReRe automatic estimate (mean ± SD, 30 samples)")
for xi,yi in zip(t,e): ax.annotate(f"{yi:.0f}%",(xi,yi),textcoords="offset points",xytext=(7,-11),fontsize=8.5,color="#33475f")
ax.set_xlabel("True careless rate"); ax.set_ylabel("Estimated careless rate (automatic)")
ax.set_title("Automatic careless-rate estimation vs ground truth\nReal study data: 38 careful respondents + sampled real careless",fontsize=11.5)
ax.set_xlim(0,52); ax.set_ylim(0,52)
ax.xaxis.set_major_formatter(matplotlib.ticker.PercentFormatter())
ax.yaxis.set_major_formatter(matplotlib.ticker.PercentFormatter())
ax.grid(alpha=0.25); ax.legend(loc="upper left",fontsize=9.5,framealpha=0.9)
ax.text(0.98,0.03,"tracks truth up to ~25%, then plateaus\n(careless-as-minority assumption)",
        transform=ax.transAxes,ha="right",va="bottom",fontsize=8.5,style="italic",color="#6b6457")
plt.tight_layout()
for p in ["fig_rate_calibration.png", r"C:/Users/vitto/Downloads/ReReRe_rate_calibration.png"]:
    plt.savefig(p,dpi=150)
print("wrote fig_rate_calibration.png and Downloads/ReReRe_rate_calibration.png")
