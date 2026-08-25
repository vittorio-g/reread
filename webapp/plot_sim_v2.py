import sys, pandas as pd
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
plt.rcParams.update({"font.family":"serif","font.size":11})
OUT=sys.argv[1] if len(sys.argv)>1 else "."
sc=pd.read_csv("sim_v2_scale.csv").sort_values("items")
pv=pd.read_csv("sim_v2_prev.csv").sort_values("true_rate")
# method -> (label, color, marker)
M={"std":("Standard (z=2.5)","#33475f","o"),"high":("High (z=1.5)","#8f3535","s"),
   "A":("A: pi top-share","#2f6b3a","^"),"B":("B: Bayes boundary","#8a6d1f","D"),
   "C":("C: FDR 10%","#5a5a8f","v"),"D":("D: adaptive z","#c46210","P"),
   "E":("E: 3-component","#b02f8a","X")}

# Fig 1: MCC vs length + oracle ceiling
fig,ax=plt.subplots(figsize=(7.6,4.8)); x=sc["items"]
ax.plot(x,sc["oracle"],"--",color="#999",lw=1.6,label="oracle ceiling")
for m,(lab,c,mk) in M.items():
    ax.errorbar(x,sc[m+"_mcc"],yerr=sc[m+"_mcc_se"],marker=mk,color=c,lw=1.9,capsize=3,label=lab)
ax.set_xlabel("Questionnaire length (items)"); ax.set_ylabel("MCC (new careless definition)")
ax.set_xticks(x); ax.set_ylim(0,1.0); ax.grid(alpha=.3); ax.legend(fontsize=9,loc="lower right")
ax.set_title("Flagging method vs length: MCC",fontsize=11.5)
plt.tight_layout(); plt.savefig(f"{OUT}/v2_fig1_mcc_length.png",dpi=150); plt.close()

# Fig 2: recall & precision vs length (two panels)
fig,axes=plt.subplots(1,2,figsize=(10.2,4.4))
for ax,metric,ttl in [(axes[0],"rec","Recall"),(axes[1],"prec","Precision")]:
    for m,(lab,c,mk) in M.items(): ax.plot(x,sc[m+"_"+metric],marker=mk,color=c,lw=1.9,label=lab)
    ax.set_xlabel("Questionnaire length (items)"); ax.set_ylabel(ttl); ax.set_xticks(x); ax.set_ylim(0,1.02); ax.grid(alpha=.3); ax.set_title(ttl,fontsize=11)
axes[0].legend(fontsize=8.5,loc="lower right")
plt.tight_layout(); plt.savefig(f"{OUT}/v2_fig2_recall_precision.png",dpi=150); plt.close()

# Fig 3: false-positive rate on distractors (acquiescent/fatigue we do NOT want to flag)
fig,ax=plt.subplots(figsize=(7.6,4.6))
for m,(lab,c,mk) in M.items(): ax.plot(x,sc[m+"_fprD"],marker=mk,color=c,lw=1.9,label=lab)
ax.set_xlabel("Questionnaire length (items)"); ax.set_ylabel("Fraction of distractors flagged (lower=better)")
ax.set_xticks(x); ax.set_ylim(0,None); ax.grid(alpha=.3); ax.legend(fontsize=9)
ax.set_title("Wrongly flagging response styles (acquiescent / fatigue)",fontsize=11.5)
plt.tight_layout(); plt.savefig(f"{OUT}/v2_fig3_distractor_fpr.png",dpi=150); plt.close()

# Fig 4: calibration — estimated vs true rate (prevalence sweep, length 120)
fig,ax=plt.subplots(figsize=(7.0,5.4)); tr=pv["true_rate"]*100
ax.plot([0,32],[0,32],"--",color="#8a8577",lw=1.3,label="perfect (y=x)")
for m,(lab,c,mk) in M.items(): ax.plot(tr,pv[m+"_est"]*100,marker=mk,color=c,lw=1.9,label=lab)
ax.set_xlabel("True careless rate (%)"); ax.set_ylabel("Estimated / flagged rate (%)")
ax.set_xlim(3,32); ax.set_ylim(0,34); ax.grid(alpha=.3); ax.legend(fontsize=9,loc="upper left")
ax.set_title("Rate calibration @120 items (+15% distractors)",fontsize=11.5)
plt.tight_layout(); plt.savefig(f"{OUT}/v2_fig4_calibration.png",dpi=150); plt.close()

# Fig 5: MCC vs true rate (prevalence sweep)
fig,ax=plt.subplots(figsize=(7.0,4.6))
ax.plot(tr,pv["oracle"],"--",color="#999",lw=1.6,label="oracle ceiling")
for m,(lab,c,mk) in M.items(): ax.plot(tr,pv[m+"_mcc"],marker=mk,color=c,lw=1.9,label=lab)
ax.set_xlabel("True careless rate (%)"); ax.set_ylabel("MCC"); ax.set_xlim(3,32); ax.set_ylim(0,1.0); ax.grid(alpha=.3); ax.legend(fontsize=9,loc="lower right")
ax.set_title("Flagging method vs prevalence @120 items",fontsize=11.5)
plt.tight_layout(); plt.savefig(f"{OUT}/v2_fig5_mcc_prevalence.png",dpi=150); plt.close()
print("wrote v2_fig1..5")

# ---- weighted-MCC (U-shape: extremes heavy, ambiguous middle light) ----
try:
    # Fig 7: weighted-MCC vs length
    fig,ax=plt.subplots(figsize=(7.8,4.8)); x=sc["items"]
    for m,(lab,c,mk) in M.items():
        lw=2.4 if m=="D" else 1.7
        ax.errorbar(x,sc[m+"_wmcc"],yerr=sc[m+"_wmcc_se"],marker=mk,color=c,lw=lw,capsize=3,label=lab)
    ax.set_xlabel("Questionnaire length (items)"); ax.set_ylabel("Weighted MCC (extremes heavy)")
    ax.set_xticks(x); ax.set_ylim(0.3,0.85); ax.grid(alpha=.3); ax.legend(fontsize=8.5,ncol=2,loc="lower right")
    ax.set_title("Weighted MCC vs length (ambiguous middle down-weighted)",fontsize=11.5)
    plt.tight_layout(); plt.savefig("v2_fig7_wmcc_length.png",dpi=150); plt.close()
    # Fig 8: weighted-MCC vs prevalence
    fig,ax=plt.subplots(figsize=(7.8,4.8)); tr=pv["true_rate"]*100
    for m,(lab,c,mk) in M.items():
        lw=2.4 if m=="D" else 1.7
        ax.plot(tr,pv[m+"_wmcc"],marker=mk,color=c,lw=lw,label=lab)
    ax.set_xlabel("True careless rate (%)"); ax.set_ylabel("Weighted MCC")
    ax.set_xlim(3,32); ax.set_ylim(0.3,0.85); ax.grid(alpha=.3); ax.legend(fontsize=8.5,ncol=2,loc="lower left")
    ax.set_title("Weighted MCC vs prevalence @120 items",fontsize=11.5)
    plt.tight_layout(); plt.savefig("v2_fig8_wmcc_prevalence.png",dpi=150); plt.close()
    print("wrote v2_fig7, v2_fig8 (weighted MCC)")
except Exception as e:
    print("weighted-MCC plots skipped:",e)
