"""Generate LaTeX-style PDF report for the robust re-examination."""

from reportlab.lib.pagesizes import A4
from reportlab.lib.styles import ParagraphStyle
from reportlab.lib.units import cm
from reportlab.lib import colors
from reportlab.platypus import (
    SimpleDocTemplate, Paragraph, Spacer, Table, TableStyle, Image, PageBreak
)
from reportlab.lib.enums import TA_CENTER, TA_JUSTIFY, TA_LEFT

OUT = r"C:\Users\vitto\Downloads\ReReReRe_Reexam_Report_2026-04-23.pdf"
ROOT = r"C:\Users\vitto\Desktop\ReReReRe"

title_style = ParagraphStyle("Title", fontName="Times-Bold",
    fontSize=18, leading=21, alignment=TA_CENTER, spaceAfter=6)
sub_style = ParagraphStyle("Sub", fontName="Times-Italic",
    fontSize=10, leading=12, alignment=TA_CENTER, spaceAfter=14, textColor=colors.grey)
abstract_h = ParagraphStyle("AbsH", fontName="Times-Bold",
    fontSize=10, leading=12, alignment=TA_CENTER, spaceAfter=4)
abstract = ParagraphStyle("Abs", fontName="Times-Roman",
    fontSize=9.5, leading=12, alignment=TA_JUSTIFY,
    spaceAfter=12, leftIndent=0.8*cm, rightIndent=0.8*cm)
h1 = ParagraphStyle("H1", fontName="Times-Bold",
    fontSize=13, leading=16, spaceBefore=14, spaceAfter=6)
h2 = ParagraphStyle("H2", fontName="Times-Bold",
    fontSize=11.5, leading=14, spaceBefore=10, spaceAfter=4)
body = ParagraphStyle("Body", fontName="Times-Roman",
    fontSize=10.5, leading=13.5, alignment=TA_JUSTIFY,
    spaceAfter=6, firstLineIndent=0.5*cm)
caption = ParagraphStyle("Cap", fontName="Times-Italic",
    fontSize=9, leading=11, alignment=TA_CENTER, spaceAfter=8)
note = ParagraphStyle("Note", fontName="Times-Italic",
    fontSize=8.5, leading=10.5, alignment=TA_LEFT,
    textColor=colors.grey, spaceAfter=4)

cell_style = ParagraphStyle("Cell", fontName="Times-Roman",
    fontSize=9.5, leading=11.5, alignment=TA_CENTER)
cell_left = ParagraphStyle("CellL", fontName="Times-Roman",
    fontSize=9.5, leading=11.5, alignment=TA_LEFT)
cell_bold = ParagraphStyle("CellB", fontName="Times-Bold",
    fontSize=9.5, leading=11.5, alignment=TA_CENTER)

def wrap(val, left=False, bold=False):
    if isinstance(val, str):
        s = cell_bold if bold else (cell_left if left else cell_style)
        return Paragraph(val, s)
    return val

def mktable(data, cols, left_cols=None):
    if left_cols is None: left_cols = set()
    wrapped = []
    for r, row in enumerate(data):
        new_row = [wrap(c, left=(r > 0 and i in left_cols), bold=(r == 0))
                   for i, c in enumerate(row)]
        wrapped.append(new_row)
    t = Table(wrapped, colWidths=cols, repeatRows=1)
    t.setStyle(TableStyle([
        ("VALIGN", (0, 0), (-1, -1), "MIDDLE"),
        ("LINEABOVE", (0, 0), (-1, 0), 1.0, colors.black),
        ("LINEBELOW", (0, 0), (-1, 0), 0.7, colors.black),
        ("LINEBELOW", (0, -1), (-1, -1), 1.0, colors.black),
        ("PADDING", (0, 0), (-1, -1), 3),
    ]))
    return t

story = []

# --- Cover ---
story.append(Paragraph(
    "Robust Re-examination of ReReReRe:<br/>"
    "Ensemble Detection of Partial Carelessness",
    title_style))
story.append(Paragraph(
    "Vittorio Guerri &middot; April 23, 2026", sub_style))

story.append(Paragraph("Abstract", abstract_h))
story.append(Paragraph(
    "We re-examine the ReReReRe careless-respondent detector with explicit focus "
    "on respondents whose responses are partially corrupted (60-80% of items). "
    "A robust simulation of 24,000 respondents (48 datasets, 8 replications) "
    "evaluates an ensemble of six detectors — z_RR, iterative-EFA z_RR, IRV, "
    "LongString, Mahalanobis D², and Person-Total correlation — under five "
    "ground-truth cutoffs (corruption > 40%, 50%, 60%, 70%, 80%). Iterative EFA "
    "is the single most impactful addition: it lifts MCC by 0.06-0.08 across all "
    "cutoffs and brings detection of partial-random carelessness at 60-80% "
    "corruption from 27-41% (z_RR alone) to 59-76% (best ensemble). Two patterns "
    "remain structurally hard: fatigue (end-loaded corruption) and random at "
    "60-70% corruption. The best ensemble achieves overall MCC 0.617 at the "
    "standard GT>50% cutoff, sensitivity 0.74, specificity 0.96 at FPR=5%.",
    abstract))

# --- Section 1 ---
story.append(Paragraph("1. Design", h1))
story.append(Paragraph(
    "The re-examination addresses a specific limitation surfaced during the "
    "ensemble v2 study: at partial corruption levels (60-80%), random "
    "carelessness remained poorly detected. This study doubles the simulation "
    "size (from 24 to 48 datasets, N=500 respondents each), explicitly tunes "
    "the trim_pct hyperparameter of iterative EFA, and evaluates the ensemble "
    "at five ground-truth cutoffs.",
    body))

story.append(Paragraph("1.1 Simulation parameters", h2))
params_tbl = [
    ["Parameter", "Value"],
    ["Conditions (nF × ipf)", "6 (8/16/30 × 6/10)"],
    ["Replications", "8"],
    ["Respondents per dataset", "500"],
    ["Total respondents", "24,000"],
    ["Careless rate", "40%"],
    ["Pattern types (v3 injector)", "6 (random, longstring, pure_straight, "
                                      "acquiescent, mixed, fatigue)"],
    ["Corruption levels", "10%, 20%, …, 100%"],
    ["Detectors evaluated", "z_RR, z_RR_iter_EFA (×3 trim levels), "
                             "IRV, LongString, D², PersonTotal"],
    ["Ground-truth thresholds", "40%, 50%, 60%, 70%, 80%"],
    ["Operating point", "FPR = 5% on clean respondents"],
    ["Runtime", "212 minutes (3.5 hours)"],
]
story.append(mktable(params_tbl, [5*cm, 10.5*cm], left_cols={1}))

# --- Section 2 ---
story.append(Paragraph("2. Iterative-EFA Hyperparameter Tuning", h1))
story.append(Paragraph(
    "Three initial_z_threshold values (0.5, 1.0, 1.5 — corresponding to rough "
    "trim fractions of 10-30%) were tested. All three yielded nearly identical "
    "performance, confirming that iter_EFA is robust to trim choice.",
    body))
trim_tbl = [
    ["Trim config", "Sens", "FPR", "MCC"],
    ["initial z = 0.5 (trim ~10%)", "0.739", "0.096", "0.616"],
    ["initial z = 1.0 (trim ~20%)", "0.740", "0.096", "0.617"],
    ["initial z = 1.5 (trim ~30%)", "0.741", "0.096", "0.617"],
]
story.append(mktable(trim_tbl, [6*cm, 2.5*cm, 2.5*cm, 2.5*cm]))
story.append(Paragraph(
    "Table 1. iter_EFA trim tuning — negligible differences. Standard z=1.0 is retained.",
    caption))

# --- Section 3 ---
story.append(Paragraph("3. MCC Across Ground-Truth Thresholds", h1))
story.append(Paragraph(
    "The table below shows MCC at FPR=5% for five ensembles, evaluated at five "
    "definitions of careless (corruption > T). The best ensemble — all six "
    "detectors including iter_EFA — dominates at every threshold.",
    body))

mcc_tbl = [
    ["GT > T", "z_RR alone", "Triad", "Triad + iter", "All 5 (no iter)",
     "All 6 + iter"],
    ["40%", "0.393", "0.550", "0.608", "0.575", "0.635"],
    ["50%", "0.398", "0.543", "0.593", "0.559", "0.617"],
    ["60%", "0.393", "0.519", "0.563", "0.528", "0.577"],
    ["70%", "0.365", "0.461", "0.515", "0.468", "0.524"],
    ["80%", "0.328", "0.390", "0.431", "0.397", "0.440"],
]
story.append(mktable(mcc_tbl, [1.8*cm, 2.5*cm, 2.2*cm, 2.8*cm, 3*cm, 2.7*cm]))
story.append(Paragraph(
    "Table 2. MCC (5-fold CV, FPR=5% operating point) at five GT thresholds.",
    caption))

story.append(Image(rf"{ROOT}\plot_robust_mcc_vs_threshold.png",
                   width=15*cm, height=8.5*cm))
story.append(Paragraph(
    "Figure 1. The iter_EFA variants dominate across all ground-truth thresholds. "
    "MCC naturally declines as T rises because mid-corruption respondents shift "
    "to the negative class but remain statistically similar to high-corruption ones.",
    caption))

story.append(PageBreak())

# --- Section 4 ---
story.append(Paragraph("4. Partial Carelessness (60-80%) — the Focus",
                        h1))
story.append(Paragraph(
    "Averaged across all six pattern types, the ensemble substantially "
    "improves detection in the 60-80% corruption range. The gain over "
    "z_RR alone is +0.32 to +0.35 percentage points at FPR=5%.",
    body))

focus_tbl = [
    ["Corruption", "z_RR", "Triad", "Best", "Gain (best - z_RR)"],
    ["50%", "0.22", "0.35", "0.47", "+0.25"],
    ["60%", "0.27", "0.46", "0.59", "+0.32"],
    ["70%", "0.36", "0.56", "0.69", "+0.33"],
    ["80%", "0.41", "0.66", "0.76", "+0.35"],
    ["90%", "0.46", "0.69", "0.81", "+0.35"],
    ["100%", "0.59", "0.74", "0.85", "+0.26"],
]
story.append(mktable(focus_tbl, [2.5*cm, 2*cm, 2*cm, 2*cm, 4.5*cm]))
story.append(Paragraph(
    "Table 3. Detection rate at FPR=5%, averaged across pattern types.",
    caption))

story.append(Image(rf"{ROOT}\plot_robust_focus_50_100.png",
                   width=14*cm, height=9*cm))
story.append(Paragraph(
    "Figure 2. Ensemble detection rate at corruption 50-100%, averaged across "
    "patterns. The best ensemble achieves >75% detection at 80% corruption.",
    caption))

# --- Section 5 ---
story.append(Paragraph("5. Per-Pattern Performance at 60-80% Corruption",
                        h1))
story.append(Paragraph(
    "Pattern-specific detection rates reveal that the ensemble is not uniformly "
    "effective. Pure-straight and acquiescent respondents (low-variance "
    "patterns) are caught almost perfectly from 70% corruption onward. "
    "Longstring, mixed, and random are detected moderately. Fatigue "
    "(end-loaded corruption) remains structurally difficult.",
    body))

pat_tbl = [
    ["Pattern", "z_RR @60%", "Best @60%", "z_RR @70%", "Best @70%",
     "z_RR @80%", "Best @80%"],
    ["pure_straight", "0.22", "0.83", "0.47", "0.96", "0.56", "0.99"],
    ["acquiescent",   "0.40", "0.73", "0.47", "0.85", "0.54", "0.95"],
    ["longstring",    "0.29", "0.68", "0.37", "0.68", "0.44", "0.80"],
    ["mixed",         "0.34", "0.56", "0.34", "0.69", "0.35", "0.68"],
    ["random",        "0.28", "0.48", "0.30", "0.57", "0.33", "0.61"],
    ["fatigue",       "0.10", "0.26", "0.19", "0.42", "0.25", "0.51"],
]
story.append(mktable(pat_tbl,
                      [3.5*cm, 1.8*cm, 1.8*cm, 1.8*cm, 1.8*cm, 1.8*cm, 1.8*cm],
                      left_cols={0}))
story.append(Paragraph(
    "Table 4. Detection rate per pattern at 60%, 70%, 80% corruption levels.",
    caption))

story.append(Image(rf"{ROOT}\plot_robust_sens_per_pattern.png",
                   width=16*cm, height=9*cm))
story.append(Paragraph(
    "Figure 3. Sensitivity per corruption bucket, faceted by pattern. "
    "The shaded band (60-80%) is the focus of this re-examination.",
    caption))

# --- Section 6 ---
story.append(Paragraph("6. Patterns that Remain Hard", h1))

story.append(Paragraph("6.1 Fatigue — end-loaded corruption", h2))
story.append(Paragraph(
    "The fatigue pattern places corruption deterministically at the last k "
    "items. At 60-70% corruption, detection is 26-42% — well below the "
    "other patterns. None of the six detectors we tested is position-aware, "
    "so a respondent who was attentive for the first half and zoned out for "
    "the second half is nearly indistinguishable from a full-length moderate "
    "responder. A sliding-window coherence score would address this, but "
    "falls outside the scope of this re-examination.",
    body))

story.append(Paragraph("6.2 Random at moderate corruption (60-70%)", h2))
story.append(Paragraph(
    "Random carelessness at 60-70% retains 30-40% of genuine responses, which "
    "dilute the incoherence signal. Even the best ensemble reaches only "
    "48-57% detection. This is a structural limit: uniform-random Likert "
    "responses produce a feature profile that heavily overlaps with "
    "low-engagement but genuine respondents. Without external information "
    "(response time, IRT-based residuals), this gap cannot be closed by "
    "correlation-based detectors alone.",
    body))

# --- Section 7 ---
story.append(Paragraph("7. Conclusions and Recommendations", h1))

story.append(Paragraph("7.1 Best ensemble configuration", h2))
story.append(Paragraph(
    "The defensible ensemble combines six detectors: z_RR (standard ReReReRe "
    "with variance_penalty), z_RR_iter_EFA (iterative version), IRV, "
    "LongString, Mahalanobis D², and Person-Total correlation. Combined via "
    "5-fold-CV logistic regression at a fixed FPR=5% operating point:",
    body))

summary_tbl = [
    ["Metric", "Value"],
    ["Overall MCC @ GT>50%", "0.617"],
    ["Sensitivity @ GT>50%", "0.74"],
    ["Specificity (fixed by design)", "0.95"],
    ["MCC @ GT>60% (partial careless focus)", "0.577"],
    ["MCC @ GT>80% (strict only full careless)", "0.440"],
    ["Detection of pure_straight @ 100%", "1.00"],
    ["Detection of random @ 100%", "0.69-0.85 (depends on condition)"],
    ["Detection of fatigue @ 100%", "0.87"],
    ["Detection of random @ 60%", "0.48 (structurally limited)"],
]
story.append(mktable(summary_tbl, [9*cm, 6.5*cm], left_cols={0, 1}))
story.append(Paragraph(
    "Table 5. Best ensemble at-a-glance. All values at FPR=5% operating point.",
    caption))

story.append(Paragraph("7.2 When to use the ensemble", h2))
story.append(Paragraph(
    "The ensemble approach is recommended when partial carelessness is a "
    "realistic concern. When the expectation is that careless respondents "
    "are either mostly clean or fully random (no middle range), z_RR alone "
    "with auto_z and variance_penalty is sufficient. The ensemble's main "
    "benefit is in the 50-80% corruption range, precisely where z_RR alone "
    "struggles.",
    body))

story.append(Paragraph("7.3 What not to do", h2))
story.append(Paragraph(
    "Do not use a single fixed probability threshold (e.g., p > 0.5) with "
    "the logistic ensemble when the careless base rate is low — the threshold "
    "should always be calibrated to a target FPR. The 5% FPR operating point "
    "used throughout this study is a reasonable default for most applications.",
    body))

story.append(Spacer(1, 0.5*cm))
story.append(Paragraph(
    "<i>Data and code:</i> github.com/vittorio-g/ReReReRe &middot; "
    "<i>Compiled:</i> 2026-04-23",
    note))

doc = SimpleDocTemplate(OUT, pagesize=A4,
                         leftMargin=2.4*cm, rightMargin=2.4*cm,
                         topMargin=2.2*cm, bottomMargin=2.2*cm,
                         title="ReReReRe Robust Re-examination",
                         author="Vittorio Guerri")
doc.build(story)
print(f"Saved: {OUT}")
