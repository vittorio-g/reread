import pandas as pd, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
plt.rcParams.update({"font.family":"serif","font.size":11})
d=pd.read_csv("eval_sens_final.csv").sort_values("true_rate"); x=d.true_rate*100
fig,ax=plt.subplots(figsize=(8.6,5.3))
ax.axvspan(10,25,color="#2f6b3a",alpha=.08,label="green range (10–25%)")
ax.plot(x,d.low,   "-o",color="#33475f",lw=1.8,ms=3.8,label="Low sensitivity")
ax.plot(x,d.medium,"-s",color="#2f6b3a",lw=1.8,ms=3.8,label="Medium (default)")
ax.plot(x,d["high"],"-^",color="#8f3535",lw=1.8,ms=4,label="High sensitivity")
ax.annotate("Low wins\n(clean data)",(8.5,0.905),ha="center",fontsize=9,color="#33475f",style="italic")
ax.annotate("Medium wins",(17.5,0.915),ha="center",fontsize=9,color="#2f6b3a",style="italic")
ax.annotate("High wins\n(heavy contamination)",(40,0.86),ha="center",fontsize=9,color="#8f3535",style="italic")
ax.set_xlabel("Careless rate (%)"); ax.set_ylabel("MCC (automatic flags)")
ax.set_title("Which sensitivity setting to use depends on how contaminated the data are\n"
             "(58 fixed attentive respondents; independent resampling per level, 120/point)",fontsize=11)
ax.set_xticks(range(5,51,5)); ax.set_xlim(2.5,52.5); ax.set_ylim(0.25,0.95); ax.grid(alpha=.25)
ax.legend(loc="lower left",fontsize=9.5,framealpha=.92)
plt.tight_layout()
for p in ["fig_sensitivity_final.png", r"C:/Users/vitto/Downloads/ReReRe_sensitivity_final.png"]:
    plt.savefig(p,dpi=160)
print("wrote fig_sensitivity_final.png")
