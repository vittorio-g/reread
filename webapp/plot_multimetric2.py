import pandas as pd, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
plt.rcParams.update({"font.family":"serif","font.size":11})
OUTS=["C:/Users/vitto/Downloads/rerere_overleaf/figures/fig1_multimetric.png","fig1_multimetric.png"]
GREEN="#2f6b3a"; RUST="#8f3535"
d=pd.read_csv("eval_sens2.csv").sort_values("true"); x=100*d["true"].values
fig,ax=plt.subplots(figsize=(7.6,5.0))
ax.axvspan(10,20,color=GREEN,alpha=.07)
ax.plot(x,d["std_prec"], "-o",color=GREEN,lw=2.0,ms=4.5,label="Precision — shipped cut (z = 2.5)")
ax.plot(x,d["std_rec"],  "--o",color=GREEN,lw=1.8,ms=4.5,mfc="none",label="Recall — shipped cut")
ax.plot(x,d["high_prec"],"-s",color=RUST,lw=2.0,ms=4.5,label="Precision — liberal cut (z = 1.5, not shipped)")
ax.plot(x,d["high_rec"], "--s",color=RUST,lw=1.8,ms=4.5,mfc="none",label="Recall — liberal cut")
ax.set_xlabel("True careless rate (%)"); ax.set_ylabel("Precision / Recall")
ax.set_title("Automatic-flag precision and recall\n(one-sided attentive-anchored cut; 200 resamples per rate)",fontsize=10.5)
ax.set_ylim(0,1.02); ax.grid(alpha=.3); ax.legend(fontsize=9,loc="lower left",framealpha=.9)
plt.tight_layout()
for p in OUTS: plt.savefig(p,dpi=160)
print("wrote fig1_multimetric")
