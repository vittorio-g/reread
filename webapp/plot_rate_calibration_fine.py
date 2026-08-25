import pandas as pd, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
plt.rcParams.update({"font.family":"serif","font.size":11})
df=pd.read_csv("eval_mix_fine.csv").sort_values("true_rate")
t=df.true_rate.values*100
reps=df.reps.values if "reps" in df.columns else 300
def ci(col): return 1.96*df[col].values*100/(reps**0.5)   # 95% CI of the mean
sm=df.std_mean.values*100; hm=df.high_mean.values*100
fig,ax=plt.subplots(figsize=(7.6,6.4))
ax.plot([8,32],[8,32],"--",color="#8a8577",lw=1.4,label="perfect calibration (y = x)")
ax.errorbar(t,hm,yerr=ci("high_sd"),marker="s",ms=6,lw=2,capsize=4,color="#8f3535",
            label="High setting (z = 1.5)")
ax.errorbar(t,sm,yerr=ci("std_sd"),marker="o",ms=6,lw=2,capsize=4,color="#2f6b3a",
            label="Standard setting (z = 2.5, default)")
ax.set_xlabel("True careless rate"); ax.set_ylabel("Estimated careless rate (automatic)")
ax.set_title("Automatic careless-rate estimation vs ground truth\n"
             "Real study data, realistic 10–30% range · 300 resamples/point",fontsize=11.5)
ax.set_xlim(9,31); ax.set_ylim(6,32)
ax.xaxis.set_major_formatter(matplotlib.ticker.PercentFormatter())
ax.yaxis.set_major_formatter(matplotlib.ticker.PercentFormatter())
ax.grid(alpha=0.25); ax.legend(loc="upper left",fontsize=9.5,framealpha=0.9)
plt.tight_layout()
for p in ["fig_rate_calibration_fine.png", r"C:/Users/vitto/Downloads/ReReRe_rate_calibration_fine.png"]:
    plt.savefig(p,dpi=150)
print("wrote fig_rate_calibration_fine.png and Downloads copy")
