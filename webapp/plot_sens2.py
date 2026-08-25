import pandas as pd, numpy as np, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
plt.rcParams.update({"font.family":"serif","font.size":11})
OUTS=["C:/Users/vitto/Downloads/rerere_overleaf/figures/fig_sensitivity.png","fig_sensitivity.png"]
GREEN="#2f6b3a"; RUST="#8f3535"; GREY="#8a8577"
d=pd.read_csv("eval_sens2.csv").sort_values("true")
x=100*d["true"].values; s=d["std_mcc"].values; h=d["high_mcc"].values
# crossover: first rate where High beats Standard (linear interp)
cross=None
for i in range(1,len(x)):
    if (s[i-1]>=h[i-1]) and (s[i]<h[i]):
        f=(s[i-1]-h[i-1])/((s[i-1]-h[i-1])-(s[i]-h[i])); cross=x[i-1]+f*(x[i]-x[i-1]); break
fig,ax=plt.subplots(figsize=(7.6,5.0))
if cross:
    ax.axvspan(x[0]-1,cross,color=GREEN,alpha=.07)
    ax.axvspan(cross,x[-1]+1,color=RUST,alpha=.06)
    ax.axvline(cross,color=GREY,lw=1.0,ls=":")
ax.plot(x,s,"-o",color=GREEN,lw=2.2,ms=5,label="shipped cut (z = 2.5)")
ax.plot(x,h,"-s",color=RUST,lw=2.2,ms=5,label="more liberal cut (z = 1.5, not shipped)")
if cross:
    ymin=min(s.min(),h.min())
    ax.annotate("shipped cut better",((x[0]+cross)/2,ymin+0.015),ha="center",fontsize=9,color=GREEN,style="italic")
    ax.annotate("liberal cut better",((cross+x[-1])/2,ymin+0.015),ha="center",fontsize=9,color=RUST,style="italic")
ax.set_xlabel("True careless rate (%)"); ax.set_ylabel("MCC (automatic flags)")
ax.set_title("Why one operating point is enough\n(real careful/careless mixes, 200 resamples per rate; one-sided attentive-anchored cut)",fontsize=10.5)
ax.grid(alpha=.3); ax.legend(loc="upper left",fontsize=9.5,framealpha=.92)
ax.set_xlim(x[0]-1,x[-1]+1)
plt.tight_layout()
for p in OUTS: plt.savefig(p,dpi=160)
print("wrote fig_sensitivity; crossover at %.1f%%" % (cross if cross else float('nan')))
