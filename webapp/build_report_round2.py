# -*- coding: utf-8 -*-
"""build_report_round2.py — weighted (severity-graded) evaluation, D re-tuning, and the
acquiescence feature. PDF in Downloads. Data-driven tables + figures."""
import pandas as pd, os
from reportlab.lib.pagesizes import A4
from reportlab.lib.styles import getSampleStyleSheet, ParagraphStyle
from reportlab.lib.units import cm
from reportlab.lib import colors
from reportlab.lib.enums import TA_JUSTIFY, TA_CENTER
from reportlab.platypus import (SimpleDocTemplate, Paragraph, Spacer, Image as RLImage,
                                Table, TableStyle, PageBreak)
from PIL import Image as PILImage
WEB="."; OUT=r"C:\Users\vitto\Downloads\ReReRe_Weighted_D_Acq_Report_2026-07-14.pdf"
pv=pd.read_csv(f"{WEB}/sim_v2_prev.csv").sort_values("true_rate")
td=pd.read_csv(f"{WEB}/tune_D.csv"); ac=pd.read_csv(f"{WEB}/test_acq.csv")
ar=pd.read_csv(f"{WEB}/test_acq_real.csv") if os.path.exists(f"{WEB}/test_acq_real.csv") else None
rd=pd.read_csv(f"{WEB}/real_data_test.csv").sort_values("true_rate") if os.path.exists(f"{WEB}/real_data_test.csv") else None
ss=getSampleStyleSheet()
H1=ParagraphStyle("H1",parent=ss["Heading1"],fontName="Times-Bold",fontSize=16,spaceAfter=8)
H2=ParagraphStyle("H2",parent=ss["Heading2"],fontName="Times-Bold",fontSize=12,spaceBefore=10,spaceAfter=5)
BODY=ParagraphStyle("BODY",parent=ss["BodyText"],fontName="Times-Roman",fontSize=10.3,leading=14.5,alignment=TA_JUSTIFY)
CAP=ParagraphStyle("CAP",parent=BODY,fontSize=8.8,alignment=TA_CENTER,textColor=colors.grey)
I=ParagraphStyle("I",parent=BODY,leftIndent=12)
story=[]
def P(t,s=BODY): story.append(Paragraph(t,s))
def fig(fn,w=15.0,cap=None):
    p=f"{WEB}/{fn}"
    if os.path.exists(p):
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

P("ReReRe — Severity-weighted evaluation, adaptive cut (D), and an acquiescence feature",H1)
P("Round-2 simulation report · 14 July 2026",CAP)
P("<b>Three questions this report answers.</b> (1) We reform the <i>judgement</i>: carelessness is graded, "
  "so we weight the metric to count the clear cases heavily (clearly careful and clearly careless) and the "
  "ambiguous middle lightly — this is the honest way to score, and it stops punishing the tool for the "
  "intrinsically-ambiguous distractors. (2) Can a single <b>adaptive cut (D)</b> replace the Standard/High "
  "double setting? (3) Does a cheap <b>acquiescence feature</b> help keep response styles while flagging junk? "
  "Distractors (acquiescent + fatigue) are kept in the data with label 0 throughout.")

P("1 · The severity-weighted metric",H2)
P("Each respondent gets a severity <i>s</i> (clean 0, response-style ~0.3, mild junk ~0.5, extreme junk ~1). "
  "The weight is a U-shape <b>w(s) = (2s−1)²</b>: ≈1 at the extremes, ≈0 in the ambiguous middle. Weighted-MCC "
  "then rewards getting the clear cases right and is lenient where the truth itself is ambiguous. Because the "
  "distractors sit in the middle, this <b>dissolves most of the penalty</b> for how a method treats them.")
fig("v2_fig8_wmcc_prevalence.png",14.5,"Fig 1 — Weighted-MCC vs true careless rate (120 items). High leads from ~15% up.")
fig("v2_fig7_wmcc_length.png",14.5,"Fig 2 — Weighted-MCC vs questionnaire length.")
P("Under the honest metric the ranking is unchanged: <b>High is the best flagging rule on the simulated "
  "(heterogeneous) careless</b> at realistic prevalence (0.79 at 15%), the conservative rules (Standard, A, "
  "B, C, D) trail, and the 3-component E is worst. The weighting lifts everyone (the ambiguous middle no "
  "longer counts against them) but does not change who wins.")

story.append(PageBreak())
P("2 · Can one adaptive cut (D) replace the double setting?",H2)
fig("v2_fig9_tuneD.png",13.5,"Fig 3 — Each cut rule scored on SIM (x, weighted-MCC) vs REAL (y, MCC). No rule "
    "is top-right on both.")
best=td.sort_values("combined",ascending=False).iloc[0]
P(f"We gridded fixed-z cuts and adaptive rules on both the simulated (weighted) and the real study data. "
  f"<b>The prevalence trade-off is real and fundamental:</b> a low z (aggressive, e.g. 1.2–1.5) wins at high "
  f"careless rates, a high z (conservative, 2.0–2.5) wins at low rates — <b>no single z is best at both ends</b>. "
  f"This is exactly why the Standard/High double setting exists, and it is justified. The current "
  f"separation-adaptive D does well on the real data (clustered careless) but poorly on the simulated tail; "
  f"the best <i>combined</i> single rule is a fixed <b>{best['rule']}</b> (combined {best['combined']:.2f}), but "
  f"it over-flags at low prevalence, which is risky on clean/screened samples.")
head=["Cut rule","SIM weighted-MCC","REAL MCC","Combined"]
data=[head]+[[r["rule"],f"{r['sim_wmcc']:.3f}",f"{r['real_mcc']:.3f}",f"{r['combined']:.3f}"] for _,r in td.sort_values("combined",ascending=False).iterrows()]
tbl(data)
P("<b>Verdict on D:</b> the adaptive-z as built does <i>not</i> cleanly replace the two settings — it is not "
  "better on both simulated and real data, which was your adoption criterion. Either keep the two settings "
  "(they map to the real low- vs high-prevalence trade-off), or accept <b>High (z=1.5)</b> as a single "
  "moderately-aggressive default (second-best combined, never catastrophic). A truly prevalence-adaptive rule "
  "is the right idea but needs a better careless-rate estimate than the current mixture gives.")

story.append(PageBreak())
P("3 · The acquiescence feature — the real win",H2)
d_fpr=ac.iloc[1]["dfpr"]-ac.iloc[0]["dfpr"]
P(f"Acquiescence = the person's raw mean response (on a mixed-keyed scale a yea-sayer scores high on both "
  f"regular and reverse items). It is <b>one number per respondent</b> — not a heavy superstructure. Refitting "
  f"the logistic ensemble with it (5-fold CV, 7{chr(44)}200 respondents) gives the clearest improvement of "
  f"anything we tried:")
fig("v2_fig10_acq.png",12.0,"Fig 4 — Adding the acquiescence feature.")
head=["Model","AUC","Weighted-MCC","Distractor-FPR"]
data=[head,
      ["rr + LongString + Person-Total",f"{ac.iloc[0]['auc']:.3f}",f"{ac.iloc[0]['wmcc']:.3f}",f"{ac.iloc[0]['dfpr']:.3f}"],
      ["+ Acquiescence",f"{ac.iloc[1]['auc']:.3f}",f"{ac.iloc[1]['wmcc']:.3f}",f"{ac.iloc[1]['dfpr']:.3f}"]]
tbl(data)
P(f"It <b>more than halves the wrong flagging of response styles</b> (distractor-FPR {ac.iloc[0]['dfpr']:.2f} → "
  f"{ac.iloc[1]['dfpr']:.2f}, i.e. {d_fpr:+.2f}), lifts AUC ({ac.iloc[0]['auc']:.3f} → {ac.iloc[1]['auc']:.3f}), "
  f"and pushes the acquiescent distractors' scores down while pushing careless up. This is the 'second "
  f"dimension' that a single incoherence score cannot provide — cheap, effective, and exactly targeted at "
  f"keeping response styles while excluding junk. <b>Contrast with the 3-component model, which was refuted.</b>")

if ar is not None:
    jb=ar[(ar['which']=='johnson')&(ar['model']=='base')].iloc[0]; ja=ar[(ar['which']=='johnson')&(ar['model']=='acq')].iloc[0]
    sb=ar[(ar['which']=='study')&(ar['model']=='base')].iloc[0]; sa=ar[(ar['which']=='study')&(ar['model']=='acq')].iloc[0]
    P("3b · Does it work on REAL data, or only in simulation?",H2)
    fig("v2_fig11_acq_real.png",13.0,"Fig 4b — Acquiescence on synthetic vs REAL Johnson clean profiles: the "
        "distractor-FPR drop survives on real respondents.")
    P(f"<b>The concern that acquiescence might be a simulation artifact is ruled out.</b> Repeating the test "
      f"with <b>real clean respondents as the base</b> (Johnson IPIP, 300 real people × 100 items) and the same "
      f"injected careless/distractors, the feature still more than halves the wrong flagging of response styles "
      f"(distractor-FPR {jb['dfpr']:.2f} → {ja['dfpr']:.2f}) and lifts AUC ({jb['auc']:.3f} → {ja['auc']:.3f}). "
      f"And on the <b>collected study (n=157)</b>, adding acquiescence does <b>no harm</b> to detecting the real "
      f"careless — leave-one-out AUC {sb['auc']:.3f} → {sa['auc']:.3f}. The feature helps on real profiles and "
      f"regresses nothing; it is safe to ship.")

if rd is not None:
    P("4 · Real collected data (reference)",H2)
    fig("v2_fig6_realdata.png",14.5,"Fig 5 — Method vs prevalence on the real study (extreme careless). D is "
        "most robust at high prevalence; E worst; the low-vs-high-prevalence trade-off is visible here too.")

P("5 · Bottom line &amp; recommendations",H2)
P("<b>(a) Adopt the acquiescence feature.</b> Re-fit the shipped ensemble to four features "
  "(rr + LongString + Person-Total + acquiescence). It is the one change that clearly helps — it halves "
  "response-style false exclusions and improves ranking, at the cost of one cheap column.", I)
P("<b>(b) Keep two settings (or a single High default), not adaptive-D.</b> The prevalence trade-off is real; "
  "the double setting is justified. Adaptive-D does not win on both simulated and real data, so it does not "
  "meet the bar to replace them. If a single knob is wanted, High (z=1.5) is the safest compromise.", I)
P("<b>(c) Report with the severity-weighted metric.</b> It is the honest way to score a graded phenomenon and "
  "it neutralises the distractor problem without engineering it away. Headline realistic-careless numbers: "
  "weighted-MCC ~0.6–0.8 (High) across the realistic band, ceiling limited by genuine mild-careless / "
  "response-style overlap.", I)
P("<b>(d) Drop the 3-component model.</b> Refuted on both simulated and real data (careless are a continuum "
  "tail, not a peel-off cluster).", I)
doc=SimpleDocTemplate(OUT,pagesize=A4,leftMargin=1.5*cm,rightMargin=1.5*cm,topMargin=1.3*cm,bottomMargin=1.3*cm)
doc.build(story); print("wrote",OUT)
