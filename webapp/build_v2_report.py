# -*- coding: utf-8 -*-
"""build_v2_report.py — PDF report on the realistic-careless-definition simulation
and the flexible flagging methods (Standard/High/A/B/C/D). Data-driven tables + figures."""
import pandas as pd, os
from reportlab.lib.pagesizes import A4
from reportlab.lib.styles import getSampleStyleSheet, ParagraphStyle
from reportlab.lib.units import cm
from reportlab.lib import colors
from reportlab.lib.enums import TA_JUSTIFY, TA_CENTER
from reportlab.platypus import (SimpleDocTemplate, Paragraph, Spacer, Image as RLImage,
                                Table, TableStyle, PageBreak)
WEB="."; OUT=r"C:\Users\vitto\Downloads\ReReRe_Flagging_Report_2026-07-13.pdf"
sc=pd.read_csv(f"{WEB}/sim_v2_scale.csv").sort_values("items")
pv=pd.read_csv(f"{WEB}/sim_v2_prev.csv").sort_values("true_rate")
ML={"std":"Standard z2.5","high":"High z1.5","A":"A: pi top-share","B":"B: Bayes","C":"C: FDR 10%","D":"D: adaptive z","E":"E: 3-component"}
METHODS=list(ML.keys())
rd=pd.read_csv(f"{WEB}/real_data_test.csv").sort_values("true_rate") if os.path.exists(f"{WEB}/real_data_test.csv") else None

ss=getSampleStyleSheet()
H1=ParagraphStyle("H1",parent=ss["Heading1"],fontName="Times-Bold",fontSize=16,spaceAfter=8)
H2=ParagraphStyle("H2",parent=ss["Heading2"],fontName="Times-Bold",fontSize=12,spaceBefore=10,spaceAfter=5)
BODY=ParagraphStyle("BODY",parent=ss["BodyText"],fontName="Times-Roman",fontSize=10.3,leading=14.5,alignment=TA_JUSTIFY)
CAP=ParagraphStyle("CAP",parent=BODY,fontSize=8.8,alignment=TA_CENTER,textColor=colors.grey)
story=[]
def P(t,s=BODY): story.append(Paragraph(t,s))
def fig(fn,w=15.5,cap=None):
    p=f"{WEB}/{fn}"
    if os.path.exists(p):
        from PIL import Image as PILImage
        iw,ih=PILImage.open(p).size; story.append(RLImage(p,width=w*cm,height=w*cm*ih/iw))
        if cap: P(cap,CAP)
        story.append(Spacer(1,6))
def tbl(data,fs=8.6):
    t=Table(data,hAlign="CENTER")
    t.setStyle(TableStyle([("FONT",(0,0),(-1,-1),"Times-Roman",fs),("FONT",(0,0),(-1,0),"Times-Bold",fs),
      ("BACKGROUND",(0,0),(-1,0),colors.HexColor("#33475f")),("TEXTCOLOR",(0,0),(-1,0),colors.white),
      ("GRID",(0,0),(-1,-1),0.4,colors.HexColor("#bbbbbb")),("ROWBACKGROUNDS",(0,1),(-1,-1),[colors.white,colors.HexColor("#f2f2f2")]),
      ("ALIGN",(1,0),(-1,-1),"CENTER"),("TOPPADDING",(0,0),(-1,-1),2.5),("BOTTOMPADDING",(0,0),(-1,-1),2.5)]))
    story.append(t); story.append(Spacer(1,6))

def row300():
    r=sc.iloc[-1]
    best=max(METHODS,key=lambda m:r[m+"_mcc"])
    return int(r["items"]),best,r
J300,best300,r300=row300()

P("ReReRe — Realistic careless definition &amp; flexible flagging",H1)
P("Simulation report · 13 July 2026 · shipped ensemble (rr + LongString + Person-Total)",CAP)
P("<b>What changed.</b> We redefined the ground truth to match what we actually want to exclude: "
  "a respondent is <b>careless (to exclude)</b> only if they gave <b>genuinely junk answers on more than "
  "half</b> the questionnaire — i.e. an inconsistent/straightlining pattern (random, mixed, longstring, "
  "pure-straight) at &gt;50% corruption. Two other perturbations are now <b>injected but kept in the data "
  "with label 0</b>: <b>acquiescent</b> responders (yea-saying — a response <i>style</i>, real if biased data) "
  "and <b>fatigue</b> (partial inattention in ≤50% of items, mostly-good data). This makes the negative class "
  "realistic: the detector must separate real junk from response styles it should NOT exclude.")
P("<b>Methods compared.</b> Six ways to turn the ranking (eta) into flags: <b>Standard</b> (cut = m1+2.5·sd1), "
  "<b>High</b> (m1+1.5·sd1), <b>A</b> (flag the top-π share, π = mixture careless proportion), <b>B</b> (Bayes "
  "boundary, posterior&gt;0.5), <b>C</b> (Benjamini-Hochberg FDR at 10%), <b>D</b> (adaptive z that lowers as the "
  "two clusters separate — i.e. as the questionnaire lengthens). All share the same BIC gate (flag nobody when no "
  "distinct careless subpopulation exists). <b>AUC</b> and <b>oracle-MCC</b> are the ranking / ceiling references.")

P("1 · Detection scales with length (all methods)",H2)
fig("v2_fig1_mcc_length.png",15.0,"Fig 1 — MCC vs questionnaire length for each flagging method; dashed = oracle ceiling.")
P("2 · Recall and precision",H2)
fig("v2_fig2_recall_precision.png",16.0,"Fig 2 — Recall (left) and precision (right) vs length.")
P("3 · Do methods wrongly flag response styles?",H2)
fig("v2_fig3_distractor_fpr.png",14.5,"Fig 3 — Fraction of <i>distractors</i> (acquiescent/fatigue, label 0) that get flagged. Lower is better.")
story.append(PageBreak())
P("4 · Rate calibration and prevalence",H2)
fig("v2_fig4_calibration.png",13.0,"Fig 4 — Estimated/flagged rate vs true careless rate (120 items).")
fig("v2_fig5_mcc_prevalence.png",14.0,"Fig 5 — MCC vs true careless rate (120 items).")

# table: methods at 300 items
P("5 · Numbers at "+str(J300)+" items",H2)
head=["Method","MCC","Recall","Precision","Est.rate","Distractor-FPR"]
data=[head]
for m in METHODS:
    data.append([ML[m],f"{r300[m+'_mcc']:.2f}",f"{r300[m+'_rec']:.2f}",f"{r300[m+'_prec']:.2f}",
                 f"{r300[m+'_est']:.2f}",f"{r300[m+'_fprD']:.2f}"])
data.append(["oracle ceiling",f"{r300['oracle']:.2f}","—","—","—","—"])
tbl(data)
story.append(Spacer(1,4))

# ---- data-driven narrative ----
r120=sc[sc["items"]==120].iloc[0]
gapHi=r120["oracle"]-r120["high_mcc"]
oc_lo,oc_hi=sc["oracle"].min(),sc["oracle"].max()
P("6 · Reading the results",H2)
P(f"<b>a) A realistic definition lowers — and is honest about — the ceiling.</b> With extreme "
  f"scrambled careless (the collected study) the ranking was near-perfect. With realistic, "
  f"heterogeneous junk (&gt;50% corruption) and response-style distractors in the negative class, the "
  f"<b>oracle ceiling is {oc_lo:.2f}–{oc_hi:.2f}</b> (AUC ~0.94, plateauing): even the best possible cut "
  f"cannot exceed this, because mild junk and response styles genuinely overlap attentive respondents. "
  f"This is the honest number — the earlier ~0.8–0.9 came from artificially extreme careless.")
P("<b>b) MCC now scales smoothly</b> with length (no erratic dips), rising to about 180 items and then "
  "plateauing as the ceiling saturates. The earlier wild behaviour was largely the old, over-inclusive "
  "definition (which counted acquiescent and mild cases as careless).")
P(f"<b>c) Aggressiveness, not the threshold formula, is what moves MCC.</b> <b>High</b> (z=1.5) wins at "
  f"every length (MCC up to {sc['high_mcc'].max():.2f}) because the realistic careless rate needs a deeper "
  f"cut — but it <b>wrongly flags ~40% of the response-style distractors</b> we decided to keep. The four "
  f"'flexible' methods we tested (A π-share, B Bayes, C FDR, D adaptive-z) all land near the conservative "
  f"<b>Standard</b> and <b>do not beat the existing options</b>. There is no free lunch in the threshold rule.")
P(f"<b>d) The gap to the oracle is real.</b> At 120 items the best method (High) reaches "
  f"{r120['high_mcc']:.2f} while the oracle is {r120['oracle']:.2f} — a gap of {gapHi:.2f}. We hypothesised "
  f"the distractors form a <b>third</b> group that a 2-group mixture can't model, and tested a fix (method E).")

I=ParagraphStyle("I",parent=BODY,leftIndent=12)
P("6b · The three-component model (E): tested and refuted",H2)
P(f"<b>E</b> fits a <b>three-component</b> mixture (clean / response-style / careless) and flags the top "
  f"component. It does deliver the one thing it was designed for — it <b>almost stops flagging response "
  f"styles</b> (distractor-FPR ≈ {sc['E_fprD'].mean():.2f} vs High's {sc['high_fprD'].mean():.2f}) — but its "
  f"<b>MCC is the lowest of all methods</b> (≈{sc['E_mcc'].mean():.2f}) because recall collapses. <b>Why it "
  f"fails:</b> a typical fit gives component means around [−1.6, 6, 52] with proportions [0.65, 0.33, 0.02]. "
  f"The fit finds the clean group, then <b>merges the distractors with most of the careless</b> into the "
  f"middle component, and isolates only a ~2% extreme tail as 'careless'. The careless are <b>not a distinct "
  f"cluster</b> — they are a heavy right tail spanning mild (overlapping the response styles) to extreme. No "
  f"3-group split can separate mild junk from response styles, because on a single score they are the same. "
  f"<b>The 3-component idea is refuted by the data.</b>")

P("7 · The REAL collected data (n = 157)",H2)
fig("v2_fig6_realdata.png",15.0,"Fig 6 — MCC vs careless rate on the real study (70 careful + 87 real "
    "scrambled careless). Real careless are extreme, so they DO form a cluster and everything works better.")
if rd is not None:
    r25=rd.iloc[(rd['true_rate']-0.25).abs().argmin()]; r50=rd.iloc[(rd['true_rate']-0.50).abs().argmin()]
    P(f"On the real data (extreme careless, oracle ~0.90) two things stand out. <b>(i) D (adaptive-z) is the "
      f"most robust method</b> — best from ~20% up: at 25% D={r25['D_mcc']:.2f} vs Standard {r25['std_mcc']:.2f} "
      f"/ High {r25['high_mcc']:.2f}, and at 50% D={r50['D_mcc']:.2f} vs High {r50['high_mcc']:.2f}. Lowering the "
      f"cut as the two clusters separate pays off when careless genuinely cluster. <b>(ii) E (3-component) is "
      f"again the worst</b> ({r25['E_mcc']:.2f} at 25%): with extreme careless it over-splits the careless "
      f"cluster into 'moderate' and 'extreme' and flags only the extreme, so it loses recall even here.")

P("8 · Bottom line &amp; recommendation",H2)
P("<b>The 3-component model does not work</b> — refuted on both simulated and real data. Careless are not a "
  "cluster you can peel off; they are a continuum tail entangled with response styles on a single score.")
P("<b>The bright spot is D (adaptive-z):</b> the most robust method across prevalence on the real data, "
  "beating Standard and High above ~20% careless. It is worth <b>proper tuning and LOO-validation as a "
  "candidate default</b> (the current slope was set by hand).")
P("<b>To separate response styles from real junk you need a new FEATURE, not a new threshold</b> — a second "
  "dimension the ensemble currently lacks: an acquiescence index (high mean-level agreement with low "
  "within-person variance) and a positional/fatigue index. That is the only lever that can lift the ceiling "
  "on realistic, heterogeneous careless.")
P("<b>For deployment today:</b> keep <b>Standard</b> as the safe default; consider promoting <b>D</b> after "
  "validation; document that <b>High</b> trades ~40% response-style false-exclusions for recall. And the "
  "paper's realistic-careless headline should be MCC ~0.55–0.60 with an oracle ceiling ~0.70 — the 0.8–0.9 "
  "holds only for uniformly extreme careless (as in the collected study).")
# save
doc=SimpleDocTemplate(OUT,pagesize=A4,leftMargin=1.6*cm,rightMargin=1.6*cm,topMargin=1.4*cm,bottomMargin=1.4*cm)
doc.build(story)
print("wrote",OUT)
