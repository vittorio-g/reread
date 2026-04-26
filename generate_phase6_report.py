"""Phase 6 — Threshold sweep report.

Documents the systematic GT-threshold sweep that answers:
'Is corruption > 80% the optimal operational definition of careless?'
"""

from reportlab.lib.pagesizes import A4
from reportlab.lib.styles import ParagraphStyle
from reportlab.lib.units import cm
from reportlab.lib import colors
from reportlab.platypus import (
    SimpleDocTemplate, Paragraph, Spacer, Table, TableStyle, PageBreak, Image
)
from reportlab.lib.enums import TA_CENTER, TA_JUSTIFY, TA_LEFT
import os

OUT = r"C:\Users\vitto\Downloads\ReReReRe_ThresholdSweep_2026-04-26.pdf"
ROOT = r"C:\Users\vitto\Desktop\ReReReRe"

title_style = ParagraphStyle("Title", fontName="Times-Bold",
    fontSize=20, leading=24, alignment=TA_CENTER, spaceAfter=8)
sub_style = ParagraphStyle("Sub", fontName="Times-Italic",
    fontSize=11, leading=14, alignment=TA_CENTER, spaceAfter=18, textColor=colors.grey)
abstract_h = ParagraphStyle("AbsH", fontName="Times-Bold",
    fontSize=10, leading=12, alignment=TA_CENTER, spaceAfter=4)
abstract = ParagraphStyle("Abs", fontName="Times-Roman",
    fontSize=9.5, leading=12, alignment=TA_JUSTIFY,
    spaceAfter=12, leftIndent=0.8*cm, rightIndent=0.8*cm)
h1 = ParagraphStyle("H1", fontName="Times-Bold",
    fontSize=14, leading=17, spaceBefore=18, spaceAfter=8,
    textColor=colors.HexColor("#1a3a6c"))
h2 = ParagraphStyle("H2", fontName="Times-Bold",
    fontSize=12, leading=15, spaceBefore=12, spaceAfter=5,
    textColor=colors.HexColor("#2c5282"))
body = ParagraphStyle("Body", fontName="Times-Roman",
    fontSize=10.5, leading=14, alignment=TA_JUSTIFY,
    spaceAfter=7, firstLineIndent=0.5*cm)
caption = ParagraphStyle("Cap", fontName="Times-Italic",
    fontSize=9, leading=11, alignment=TA_CENTER, spaceAfter=10)
note = ParagraphStyle("Note", fontName="Times-Italic",
    fontSize=8.5, leading=10.5, alignment=TA_LEFT,
    textColor=colors.grey, spaceAfter=4)

cell = ParagraphStyle("Cell", fontName="Times-Roman", fontSize=9.5,
                      leading=11.5, alignment=TA_CENTER)
cell_left = ParagraphStyle("CellL", fontName="Times-Roman", fontSize=9.5,
                            leading=11.5, alignment=TA_LEFT)
cell_bold = ParagraphStyle("CellB", fontName="Times-Bold", fontSize=9.5,
                            leading=11.5, alignment=TA_CENTER)
cell_red = ParagraphStyle("CellR", fontName="Times-Bold", fontSize=9.5,
                            leading=11.5, alignment=TA_CENTER,
                            textColor=colors.HexColor("#a62e2e"))
cell_green = ParagraphStyle("CellG", fontName="Times-Bold", fontSize=9.5,
                            leading=11.5, alignment=TA_CENTER,
                            textColor=colors.HexColor("#2e7d32"))

def wrap(val, left=False, bold=False, color=None):
    if isinstance(val, str):
        if color == "red": s = cell_red
        elif color == "green": s = cell_green
        elif bold: s = cell_bold
        elif left: s = cell_left
        else: s = cell
        return Paragraph(val, s)
    return val

def mktable(data, cols, left_cols=None, highlights=None):
    """highlights: dict {(row_idx, col_idx): 'red'|'green'} for special-coloured cells"""
    if left_cols is None: left_cols = set()
    if highlights is None: highlights = {}
    wrapped = []
    for r, row in enumerate(data):
        out = []
        for i, c in enumerate(row):
            color = highlights.get((r, i))
            out.append(wrap(c, left=(r > 0 and i in left_cols),
                             bold=(r == 0), color=color))
        wrapped.append(out)
    t = Table(wrapped, colWidths=cols, repeatRows=1)
    t.setStyle(TableStyle([
        ("VALIGN", (0, 0), (-1, -1), "MIDDLE"),
        ("LINEABOVE", (0, 0), (-1, 0), 1.0, colors.black),
        ("LINEBELOW", (0, 0), (-1, 0), 0.7, colors.black),
        ("LINEBELOW", (0, -1), (-1, -1), 1.0, colors.black),
        ("PADDING", (0, 0), (-1, -1), 3),
    ]))
    return t

def fitted_image(path, max_w=16*cm, max_h=20*cm):
    if not os.path.exists(path):
        return Paragraph(f"[Missing image: {os.path.basename(path)}]", note)
    img = Image(path)
    iw, ih = img.imageWidth, img.imageHeight
    sx, sy = max_w / iw, max_h / ih
    s = min(sx, sy, 1.0)
    img.drawWidth = iw * s
    img.drawHeight = ih * s
    return img

story = []

# ================ COVER ================
story.append(Spacer(1, 1.5*cm))
story.append(Paragraph(
    "ReReReRe Threshold Sweep:<br/>"
    "What is the optimal operational definition of careless?",
    title_style))
story.append(Paragraph(
    "Systematic search over GT corruption thresholds<br/>"
    "with 13 metrics across two scenarios",
    sub_style))
story.append(Paragraph(
    "Vittorio Guerri &middot; April 26, 2026",
    ParagraphStyle("Author", parent=sub_style, textColor=colors.black,
                   fontSize=10)))

story.append(Spacer(1, 0.8*cm))
story.append(Paragraph("Abstract", abstract_h))
story.append(Paragraph(
    "Earlier phases adopted <i>corruption &gt; 80%</i> as the operational "
    "definition of \"careless\" — chosen as the level above which the index "
    "MCC exceeds 0.7. This threshold was a starting assumption, not an "
    "optimum derived from data. This document reports a systematic sweep "
    "over candidate thresholds tau_GT in {0.10, 0.20, ..., 0.90}, in two "
    "scenarios (clean vs &gt; tau, and &le; tau vs &gt; tau), evaluated on "
    "13 classification metrics. <b>The data shows that &gt; 80% is NOT the "
    "optimal definition under any metric.</b> The argmax depends on the "
    "scenario: tau = 0.60 in Scenario A (MCC = 0.811, F1 = 0.852, AUPRC = "
    "0.942), tau = 0.30-0.40 in Scenario B (MCC = 0.704, F1 = 0.776). "
    "We discuss the trade-off and recommend reporting both scenarios in "
    "the paper.",
    abstract))

# ================ Section 1 ================
story.append(Paragraph("1. The question", h1))
story.append(Paragraph(
    "When labelling respondents as \"careless\" for benchmarking, we must "
    "choose a threshold tau_GT on the underlying corruption level: "
    "respondents with corruption &gt; tau_GT are positive (truly "
    "careless), the rest are negative. Earlier work fixed tau_GT = 0.80, "
    "motivated by the observation that the index's MCC exceeds 0.7 above "
    "this point. But this is a <i>floor</i> derived from the index's "
    "behaviour, not a peak. The right question is:",
    body))
story.append(Paragraph(
    "<i>Across all candidate thresholds, which one maximises classifier "
    "quality? Does this argmax agree across multiple metrics? Does it "
    "change under different label scenarios?</i>",
    body))

# ================ Section 2 ================
story.append(Paragraph("2. Two scenarios", h1))
story.append(Paragraph(
    "We test two ways of defining the negative class:",
    body))
sc_tbl = [
    ["Scenario", "Positives", "Negatives", "When realistic?"],
    ["A — Exclude middle",
     "corruption &gt; tau_GT",
     "clean only (corruption = 0)",
     "When the analyst can drop ambiguous mid-corruption respondents"],
    ["B — All-in",
     "corruption &gt; tau_GT",
     "all others (clean + 0 &lt; corruption &le; tau_GT)",
     "When every respondent must be classified (deployment)"],
]
story.append(mktable(sc_tbl, [3.0*cm, 3.5*cm, 5.5*cm, 4.5*cm], left_cols={0, 1, 2, 3}))

story.append(Paragraph(
    "Scenario A is closer to the early Phase 1-5 setup: it gives the "
    "classifier a clean training signal by excluding the middle. Scenario "
    "B is closer to real deployment: every respondent gets a label, "
    "including those with intermediate corruption (10-80%).",
    body))

# ================ Section 3 ================
story.append(Paragraph("3. Method", h1))
story.append(Paragraph(
    "Random Forest (200 trees, mtry=sqrt(6)) trained on the 6-detector "
    "feature set used throughout this work (z_RR, z_RR_iter_EFA, IRV, "
    "LongString, D&sup2;, Person-Total). 5-fold cross-validated "
    "out-of-fold probabilities. Operating point calibrated to FPR = 5% "
    "(95th percentile of negative-class predictions). At each tau_GT, all "
    "13 metrics are computed at the calibrated operating point. "
    "Threshold-free metrics (AUC, AUPRC) computed once. The same 24,000 "
    "simulated respondents (sim_robust_scores.csv) are used throughout — "
    "no new simulation, just relabelling.",
    body))

# ================ Section 4 — Scenario A ================
story.append(PageBreak())
story.append(Paragraph("4. Results — Scenario A (clean vs &gt; tau)", h1))

# Build full Scenario A table with highlight
data_A = [
    ["tau_GT", "MCC", "F1", "AUPRC", "Kappa", "AUC", "Sens", "TP", "n_pos"],
    ["0.10", "0.671", "0.767", "0.886", "0.656", "0.904", "0.672", "5808", "8637"],
    ["0.20", "0.704", "0.791", "0.900", "0.695", "0.924", "0.713", "5583", "7833"],
    ["0.30", "0.755", "0.826", "0.918", "0.752", "0.944", "0.779", "5232", "6715"],
    ["0.40", "0.787", "0.846", "0.930", "0.787", "0.959", "0.824", "4745", "5758"],
    ["0.50", "0.802", "0.852", "0.938", "0.802", "0.969", "0.852", "4093", "4804"],
    ["0.60", "0.811", "0.852", "0.942", "0.810", "0.976", "0.880", "3381", "3842"],
    ["0.70", "0.799", "0.835", "0.939", "0.798", "0.979", "0.886", "2701", "3050"],
    ["0.80", "0.783", "0.804", "0.939", "0.774", "0.984", "0.920", "1775", "1930"],
    ["0.90", "0.702", "0.703", "0.925", "0.678", "0.985", "0.926", "895",  "967"],
]
# Highlight argmax (row 6 = tau=0.6) and the user's previously-chosen 0.8 (row 8)
hl_A = {
    (6, 1): "green", (6, 2): "green", (6, 3): "green", (6, 4): "green",
    (8, 1): "red",   (8, 2): "red",   (8, 4): "red",
}
story.append(mktable(data_A, [1.4*cm, 1.4*cm, 1.4*cm, 1.6*cm, 1.5*cm, 1.5*cm,
                                 1.5*cm, 1.5*cm, 1.5*cm], highlights=hl_A))
story.append(Paragraph(
    "Table. Scenario A metrics across tau_GT. Green = argmax (tau = 0.60). "
    "Red = previously-assumed 0.80 (suboptimal by ~0.03 MCC).",
    caption))

story.append(Paragraph(
    "<b>Scenario A summary:</b>",
    body))
story.append(Paragraph(
    "• <b>Argmax for MCC, F1, AUPRC, Kappa: tau = 0.60</b> "
    "(MCC = 0.811, F1 = 0.852, AUPRC = 0.942).<br/>"
    "• The user's previous tau = 0.80: MCC = 0.783, F1 = 0.804, "
    "AUPRC = 0.939 — close but not optimal (Δ MCC ≈ 0.03).<br/>"
    "• AUC and Bal. Acc. peak at tau = 0.90 because they reward easy-to-rank "
    "cases (extreme positives only); they ignore PPV. <b>These are not "
    "appropriate criteria here</b>.<br/>"
    "• Sensitivity rises monotonically with tau because higher-corruption "
    "positives are easier; this is mechanical, not informative about the "
    "right operational threshold.",
    body))

# ================ Section 5 — Scenario B ================
story.append(Paragraph("5. Results — Scenario B (≤ tau vs &gt; tau)", h1))

data_B = [
    ["tau_GT", "MCC", "F1", "AUPRC", "Kappa", "AUC", "Sens", "TP", "n_pos"],
    ["0.10", "0.667", "0.760", "0.877", "0.653", "0.900", "0.667", "5760", "8637"],
    ["0.20", "0.678", "0.763", "0.881", "0.668", "0.917", "0.679", "5319", "7833"],
    ["0.30", "0.704", "0.776", "0.878", "0.699", "0.930", "0.716", "4807", "6715"],
    ["0.40", "0.704", "0.769", "0.864", "0.701", "0.938", "0.723", "4161", "5758"],
    ["0.50", "0.690", "0.749", "0.843", "0.689", "0.942", "0.717", "3443", "4804"],
    ["0.60", "0.661", "0.714", "0.805", "0.661", "0.940", "0.700", "2689", "3842"],
    ["0.70", "0.603", "0.653", "0.754", "0.603", "0.937", "0.648", "1977", "3050"],
    ["0.80", "0.548", "0.585", "0.689", "0.545", "0.937", "0.649", "1253", "1930"],
    ["0.90", "0.479", "0.483", "0.664", "0.454", "0.938", "0.689", "666",  "967"],
]
hl_B = {
    (3, 1): "green", (3, 2): "green",
    (4, 1): "green",                      (4, 4): "green",
    (8, 1): "red",   (8, 2): "red",       (8, 4): "red",
}
story.append(mktable(data_B, [1.4*cm, 1.4*cm, 1.4*cm, 1.6*cm, 1.5*cm, 1.5*cm,
                                 1.5*cm, 1.5*cm, 1.5*cm], highlights=hl_B))
story.append(Paragraph(
    "Table. Scenario B metrics across tau_GT. Green = argmax region "
    "(tau = 0.30-0.40). Red = previously-assumed 0.80 (severely "
    "suboptimal: Δ MCC ≈ 0.16).",
    caption))

story.append(Paragraph(
    "<b>Scenario B summary:</b>",
    body))
story.append(Paragraph(
    "• <b>Argmax for MCC, Kappa: tau = 0.40</b> (MCC = 0.704). "
    "F1 argmax: tau = 0.30 (F1 = 0.776). AUPRC argmax: tau = 0.20 "
    "(AUPRC = 0.881). The optimum cluster is tau = 0.30-0.40.<br/>"
    "• <b>The user's previous tau = 0.80: MCC = 0.548, F1 = 0.585, "
    "AUPRC = 0.689</b> — far from optimal. Δ MCC ≈ 0.16.<br/>"
    "• Why? In Scenario B, raising tau_GT moves people with corruption "
    "70-80% INTO the negative class. These are nearly-careless and "
    "indistinguishable from the positive class (corruption &gt; 80%). The "
    "classifier confuses them, and MCC collapses.",
    body))

# ================ Section 6 — Plot ================
story.append(PageBreak())
story.append(Paragraph("6. Visual summary", h1))
story.append(fitted_image(os.path.join(ROOT, "plot_phase6_headline_metrics.png"),
                            max_w=16*cm, max_h=12*cm))
story.append(Paragraph(
    "Figure 1. Headline metrics (MCC, F1, AUPRC, Kappa) across tau_GT in "
    "both scenarios. Scenario A (red) peaks at tau = 0.60 then declines. "
    "Scenario B (teal) peaks at tau = 0.30-0.40 then declines steeply. "
    "tau = 0.80 is in the post-peak decline in both scenarios.",
    caption))

story.append(fitted_image(os.path.join(ROOT, "plot_phase6_all_metrics.png"),
                            max_w=16*cm, max_h=11*cm))
story.append(Paragraph(
    "Figure 2. All 13 metrics across tau_GT in both scenarios. Sensitivity "
    "and Bal. Acc. rise with tau_GT (mechanical: easy positives), MCC, F1, "
    "AUPRC, Kappa peak in the middle.",
    caption))

# ================ Section 7 — Trade-off ================
story.append(PageBreak())
story.append(Paragraph("7. The precision/coverage trade-off", h1))
story.append(Paragraph(
    "The user originally framed the choice as a trade-off between index "
    "precision and the amount of careless caught. The data shows this "
    "trade-off explicitly:",
    body))

trade_tbl = [
    ["Choice of tau_GT", "Index quality", "Coverage (TP at FPR=5%)"],
    ["tau = 0.10 (very inclusive)", "MCC ≈ 0.67",   "5,808 careless caught"],
    ["tau = 0.30",                  "MCC ≈ 0.75",   "5,232 caught (90%)"],
    ["tau = 0.40 (Scenario B opt.)", "MCC ≈ 0.70 (B), 0.79 (A)", "4,161 (B) / 4,745 (A)"],
    ["tau = 0.60 (Scenario A opt.)", "MCC = 0.81 (A), 0.66 (B)", "3,381 (A) / 2,689 (B)"],
    ["tau = 0.80 (previous choice)", "MCC = 0.78 (A), 0.55 (B)", "1,775 (A) / 1,253 (B)"],
    ["tau = 0.90 (very strict)",    "MCC ≈ 0.70 (A), 0.48 (B)", "895 / 666"],
]
story.append(mktable(trade_tbl, [4.5*cm, 4.5*cm, 6.0*cm], left_cols={0, 1, 2}))
story.append(Paragraph(
    "Table. Trade-off table: each tau_GT is a different operationalisation "
    "of \"careless\" with different consequences for both quality and "
    "coverage.",
    caption))

story.append(fitted_image(os.path.join(ROOT, "plot_phase6_tradeoff.png"),
                            max_w=16*cm, max_h=10*cm))
story.append(Paragraph(
    "Figure 3. MCC and sensitivity across tau_GT in both scenarios.",
    caption))

# ================ Section 8 — Recommendation ================
story.append(PageBreak())
story.append(Paragraph("8. Recommendations", h1))

story.append(Paragraph(
    "<b>The user's previous choice tau_GT = 0.80 is suboptimal under "
    "every relevant metric in both scenarios.</b> The data argues for one "
    "of two operational definitions, depending on the use case:",
    body))

rec_tbl = [
    ["Use case", "Best tau_GT", "Best MCC", "Best F1", "Best AUPRC"],
    ["Full-careless detection (with middle excluded, Scenario A)",
     "0.60", "0.811", "0.852", "0.942"],
    ["All-respondent classification (Scenario B)",
     "0.30-0.40", "0.704", "0.776", "0.881"],
    ["Previous (tau = 0.80, Scenario A)", "0.80", "0.783", "0.804", "0.939"],
    ["Previous (tau = 0.80, Scenario B)", "0.80", "0.548", "0.585", "0.689"],
]
story.append(mktable(rec_tbl, [7.5*cm, 1.8*cm, 1.8*cm, 1.8*cm, 1.8*cm],
                      left_cols={0}))

story.append(Paragraph(
    "<b>Paper recommendation:</b> report <b>both</b> scenarios. Scenario "
    "A (tau = 0.60) for the published headline numbers — clean, "
    "high-quality detection of full carelessness with the analyst "
    "controlling the labelling. Scenario B (tau = 0.30-0.40) for the "
    "robustness section — what to expect if the index is deployed in a "
    "setting where every respondent must be classified.",
    body))

story.append(Paragraph(
    "<b>The choice of tau_GT is not just an optimisation choice — it is "
    "a substantive operational definition.</b> tau = 0.60 means \"a "
    "respondent is careless if &gt; 60% of their item answers are "
    "noise.\" tau = 0.30 means \"if &gt; 30% are noise.\" Both are "
    "defensible. tau = 0.80 (\"only if &gt; 80% are noise\") is overly "
    "conservative and leaves real careless respondents unflagged.",
    body))

# ================ Section 9 — Caveats ================
story.append(Paragraph("9. Caveats", h1))
story.append(Paragraph(
    "• The argmax at tau = 0.60 (Scenario A) and 0.30-0.40 (Scenario B) "
    "is specific to the simulation parameters: 6 conditions (8 nF x ipf "
    "pairs), N = 4000 respondents per condition, 6 careless patterns, "
    "corruption levels in steps of 0.10. The argmax could shift with "
    "different mixtures.<br/>"
    "• In Scenario A, n_pos shrinks dramatically with tau_GT (8,637 at "
    "tau = 0.10 down to 967 at tau = 0.90). Confidence intervals widen "
    "accordingly. The MCC = 0.81 at tau = 0.60 is built from 3,842 "
    "positives and 14,400 negatives — comfortably large.<br/>"
    "• AUC ranges from 0.90 (tau = 0.10) to 0.985 (tau = 0.90) and is "
    "essentially flat across the relevant region (~0.97-0.98 for tau = "
    "0.50-0.80). AUC is a poor discriminator between candidate operational "
    "thresholds.<br/>"
    "• 5-fold CV on 24K respondents with random feature subsets means "
    "results are stable across re-seedings.",
    body))

story.append(Paragraph("10. Bottom line", h1))
story.append(Paragraph(
    "The systematic threshold sweep, evaluated on 13 classification "
    "metrics across two scenarios, <b>does not support tau_GT = 0.80</b> "
    "as the operational definition of \"careless.\" The argmax is "
    "tau = 0.60 in Scenario A and tau = 0.30-0.40 in Scenario B. We "
    "recommend abandoning tau = 0.80 and reporting both alternative "
    "operational definitions in the paper, with their distinct "
    "interpretations and use cases.",
    body))

story.append(Spacer(1, 0.5*cm))
story.append(Paragraph(
    "<i>Code:</i> Phase6_ThresholdSweep.R &middot; "
    "<i>Data:</i> phase6_threshold_sweep.csv &middot; "
    "<i>Compiled:</i> 2026-04-26",
    note))

doc = SimpleDocTemplate(OUT, pagesize=A4,
                         leftMargin=2.4*cm, rightMargin=2.4*cm,
                         topMargin=2.2*cm, bottomMargin=2.2*cm,
                         title="ReReReRe Threshold Sweep",
                         author="Vittorio Guerri")
doc.build(story)
print(f"Saved: {OUT}")
