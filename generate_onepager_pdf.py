"""Generate 1-page summary PDF of the current state of the ReReReRe research."""

from reportlab.lib.pagesizes import A4
from reportlab.lib.styles import getSampleStyleSheet, ParagraphStyle
from reportlab.lib.units import cm
from reportlab.lib import colors
from reportlab.platypus import (
    SimpleDocTemplate, Paragraph, Spacer, Table, TableStyle, KeepTogether
)
from reportlab.lib.enums import TA_CENTER, TA_JUSTIFY, TA_LEFT

OUT = r"C:\Users\vitto\Downloads\ReReReRe_OnePager_2026-04-21.pdf"

styles = getSampleStyleSheet()

title_style = ParagraphStyle(
    "Title", parent=styles["Title"],
    fontSize=16, leading=18, alignment=TA_CENTER, spaceAfter=4,
    textColor=colors.HexColor("#1a3a6c"),
)
subtitle_style = ParagraphStyle(
    "Subtitle", parent=styles["Normal"],
    fontSize=9, leading=11, alignment=TA_CENTER, spaceAfter=6,
    textColor=colors.grey,
)
h1_style = ParagraphStyle(
    "H1", parent=styles["Heading2"],
    fontSize=11, leading=13, spaceBefore=6, spaceAfter=3,
    textColor=colors.HexColor("#1a3a6c"),
)
body_style = ParagraphStyle(
    "Body", parent=styles["Normal"],
    fontSize=8.5, leading=11, alignment=TA_JUSTIFY, spaceAfter=3,
)
caption_style = ParagraphStyle(
    "Caption", parent=styles["Italic"],
    fontSize=7.5, leading=9, alignment=TA_LEFT,
    textColor=colors.grey, spaceAfter=4,
)

def make_table(data, col_widths, header_color="#1a3a6c", font_size=8):
    t = Table(data, colWidths=col_widths, repeatRows=1)
    t.setStyle(TableStyle([
        ("BACKGROUND", (0, 0), (-1, 0), colors.HexColor(header_color)),
        ("TEXTCOLOR", (0, 0), (-1, 0), colors.white),
        ("FONTNAME", (0, 0), (-1, 0), "Helvetica-Bold"),
        ("FONTSIZE", (0, 0), (-1, -1), font_size),
        ("ALIGN", (0, 0), (-1, -1), "CENTER"),
        ("VALIGN", (0, 0), (-1, -1), "MIDDLE"),
        ("GRID", (0, 0), (-1, -1), 0.4, colors.grey),
        ("ROWBACKGROUNDS", (0, 1), (-1, -1),
         [colors.white, colors.HexColor("#f5f5f5")]),
        ("PADDING", (0, 0), (-1, -1), 2),
    ]))
    return t

# ============================================================
# Story
# ============================================================
story = []

# --- Header ---
story.append(Paragraph("ReReReRe — Research Status One-Pager", title_style))
story.append(Paragraph(
    "Permutation-based careless-respondent detection &middot; "
    "2026-04-21 &middot; github.com/vittorio-g/ReReReRe",
    subtitle_style))

# --- What it is ---
story.append(Paragraph("What the method does", h1_style))
story.append(Paragraph(
    "<b>ReReReRe</b> flags careless respondents in Likert-scale questionnaires by comparing "
    "their observed cross-item coherence against a per-respondent permutation baseline. "
    "For each respondent, a z-score is computed from coupled item pairs (high |r|) vs "
    "random pairs; low z = likely careless. The method is self-calibrating (no tuning on "
    "external data), handles reverse-coded items automatically via "
    "<i>align_signs</i>, and requires only the data matrix as input.",
    body_style))

# --- Architecture ---
story.append(Paragraph("Package architecture (since 2026-04-01)", h1_style))
arch = [
    ["Function", "Engine", "Status"],
    ["ReReReRe()", "Auto: weighted (<=60 items) or coupled (>60 items)", "Validated default"],
    ["ReReReRe_F()", "Per-factor via EFA-D (within-factor pairs)", "Experimental"],
]
story.append(make_table(arch, col_widths=[3.2*cm, 9*cm, 4*cm]))

# --- Multiverse summary ---
story.append(Paragraph("Multiverse evidence (Mar-Apr 2026)", h1_style))
story.append(Paragraph(
    "Across a full multiverse (85K rows) and two dedicated simulations, the method's "
    "key drivers are <b>items_per_factor</b> (31% of variance) and <b>nFactors</b> (26%). "
    "Sample size plateaus at n=300. ReReReRe overtakes Mahalanobis in AUC at nF=8 and in "
    "oracle MCC at nF=10. In the sweet spot (nF>=15, n>=300) fixed defaults give MCC=0.43.",
    body_style))
mv = [
    ["Condition", "ReReReRe MCC", "Practical Mahalanobis MCC", "Winner"],
    ["All conditions (fixed defaults)", "0.271", "0.056", "RR by 4.8x"],
    ["Recommended (nF>=10, n>=100)", "0.336", "0.056-0.060", "RR"],
    ["Sweet spot (nF>=15, n>=300)", "0.435", "~0.065", "RR decisively"],
]
story.append(make_table(mv, col_widths=[5.5*cm, 3.3*cm, 4.5*cm, 2.9*cm]))

# --- Variants tested ---
story.append(Paragraph("All improvement candidates tested (buone_idee.md + Weighted V2)", h1_style))
v = [
    ["Variant", "Oracle MCC", "Delta vs std", "Decision"],
    ["iterative EFA (idea 4)", "0.398", "+0.015", "Worth a flag, not default"],
    ["efa_d (existing)", "0.394", "+0.011", "Kept as ReReReRe_F()"],
    ["std (default)", "0.383", "-", "Unchanged"],
    ["pf_mean (idea 1)", "0.379", "-0.004", "Reject"],
    ["cross-factor baseline (idea 2)", "0.380", "-0.003", "Reject"],
    ["weighted_v2 (proper Pearson)", "0.320", "-0.025", "Reject"],
    ["proplow / pf_comb (idea 1)", "<=0.280", "<=-0.103", "Reject"],
]
story.append(make_table(v, col_widths=[6*cm, 2.5*cm, 2.8*cm, 4.9*cm]))

# --- External validation ---
story.append(Paragraph("External validation on real data", h1_style))
story.append(Paragraph(
    "Six datasets with ground truth + Johnson IPIP-300 inject-and-detect. RR AUC beats "
    "Mahalanobis AUC on 4/6; RR MCC(z=1.5) beats practical Mahalanobis MCC on 4/6. "
    "Johnson IPIP-300 (300 items, inject 5-20% careless): auto-z MCC=0.726, "
    "specificity=1.000. On Pennycook &amp; Rand (30 items), ReReReRe_F() improves the "
    "paper's key correlation r(CRT, Discernment) by +0.047 (paper r=.27 &rarr; .316).",
    body_style))

# --- Current status & next ---
story.append(Paragraph("Current status", h1_style))
status = [
    ["Area", "State"],
    ["Algorithm", "Frozen. Two-function architecture, corProp=0.03, z=1.5, align_signs=TRUE"],
    ["Simulation evidence", "Complete: full multiverse + 2D calibration + all variants + weighted_v2"],
    ["External validation", "6 datasets with GT + Johnson inject-and-detect + Pennycook replication"],
    ["Paper", "Not yet drafted. All figures, tables, and results ready to write up"],
    ["Package", "R functions only; no CRAN release planned in this phase"],
]
story.append(make_table(status, col_widths=[4*cm, 12.2*cm]))

story.append(Spacer(1, 0.15*cm))
story.append(Paragraph(
    "<b>Bottom line:</b> after testing every candidate improvement, the current two-function "
    "default configuration is the best trade-off between simplicity and performance. "
    "Next step is writing up the paper; no further algorithmic exploration is expected.",
    body_style))

# --- Footer ---
story.append(Spacer(1, 0.1*cm))
story.append(Paragraph(
    "Sources: Final Multiverse (57.6K RR calls), Full Multiverse (5,670 RR calls, 85K rows), "
    "2D Calibration (60 cells, 30 reps), All Variants (12,525 rows), Weighted V2 (45 cells, 3 reps). "
    "Commit f6b6dcf.",
    caption_style))

# ============================================================
# Build (single page, tight margins)
# ============================================================
doc = SimpleDocTemplate(OUT, pagesize=A4,
                         leftMargin=1.5*cm, rightMargin=1.5*cm,
                         topMargin=1.2*cm, bottomMargin=1.2*cm,
                         title="ReReReRe Research Status One-Pager",
                         author="Vittorio Guerri")
doc.build(story)
print(f"Saved: {OUT}")
