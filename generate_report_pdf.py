"""Generate PDF report: ReReReRe variants + Weighted V2 comparison."""

from reportlab.lib.pagesizes import A4
from reportlab.lib.styles import getSampleStyleSheet, ParagraphStyle
from reportlab.lib.units import cm
from reportlab.lib import colors
from reportlab.platypus import (
    SimpleDocTemplate, Paragraph, Spacer, PageBreak, Image, Table, TableStyle
)
from reportlab.lib.enums import TA_CENTER, TA_JUSTIFY

OUT = r"C:\Users\vitto\Downloads\ReReReRe_Report_2026-04-21.pdf"
ROOT = r"C:\Users\vitto\Desktop\ReReReRe"

styles = getSampleStyleSheet()

# Custom styles
title_style = ParagraphStyle(
    "Title", parent=styles["Title"],
    fontSize=22, leading=26, alignment=TA_CENTER, spaceAfter=16,
)
subtitle_style = ParagraphStyle(
    "Subtitle", parent=styles["Normal"],
    fontSize=12, leading=14, alignment=TA_CENTER, spaceAfter=20,
    textColor=colors.grey,
)
h1_style = ParagraphStyle(
    "H1", parent=styles["Heading1"],
    fontSize=16, leading=20, spaceBefore=18, spaceAfter=8,
    textColor=colors.HexColor("#1a3a6c"),
)
h2_style = ParagraphStyle(
    "H2", parent=styles["Heading2"],
    fontSize=13, leading=16, spaceBefore=12, spaceAfter=6,
    textColor=colors.HexColor("#2c5282"),
)
body_style = ParagraphStyle(
    "Body", parent=styles["Normal"],
    fontSize=10.5, leading=14, alignment=TA_JUSTIFY, spaceAfter=6,
)
caption_style = ParagraphStyle(
    "Caption", parent=styles["Italic"],
    fontSize=9, leading=11, alignment=TA_CENTER,
    textColor=colors.grey, spaceAfter=12,
)

def table_style(header_color="#1a3a6c"):
    return TableStyle([
        ("BACKGROUND", (0, 0), (-1, 0), colors.HexColor(header_color)),
        ("TEXTCOLOR", (0, 0), (-1, 0), colors.white),
        ("FONTNAME", (0, 0), (-1, 0), "Helvetica-Bold"),
        ("FONTSIZE", (0, 0), (-1, -1), 9),
        ("ALIGN", (0, 0), (-1, -1), "CENTER"),
        ("VALIGN", (0, 0), (-1, -1), "MIDDLE"),
        ("GRID", (0, 0), (-1, -1), 0.5, colors.grey),
        ("ROWBACKGROUNDS", (0, 1), (-1, -1),
         [colors.white, colors.HexColor("#f5f5f5")]),
        ("PADDING", (0, 0), (-1, -1), 4),
    ])

def make_table(data, col_widths=None, highlight_col=None, header_color="#1a3a6c"):
    t = Table(data, colWidths=col_widths, repeatRows=1)
    style = table_style(header_color)
    if highlight_col is not None:
        # Highlight the row with the max value in highlight_col
        col_idx = highlight_col
        values = []
        for r in range(1, len(data)):
            try:
                values.append(float(data[r][col_idx]))
            except (ValueError, TypeError):
                values.append(float("-inf"))
        if values:
            best_row = values.index(max(values)) + 1
            style.add("BACKGROUND", (0, best_row), (-1, best_row),
                      colors.HexColor("#fff3cd"))
            style.add("FONTNAME", (0, best_row), (-1, best_row), "Helvetica-Bold")
    t.setStyle(style)
    return t

# ============================================================
# Build story
# ============================================================
story = []

# --- Cover ---
story.append(Spacer(1, 3 * cm))
story.append(Paragraph("ReReReRe — Variants Comparison Report", title_style))
story.append(Paragraph(
    "All buone_idee.md variants + Weighted V2 head-to-head", subtitle_style))
story.append(Spacer(1, 1 * cm))
story.append(Paragraph(
    "<b>Date:</b> 2026-04-21<br/>"
    "<b>Author:</b> Vittorio Guerri<br/>"
    "<b>Repository:</b> github.com/vittorio-g/ReReReRe",
    ParagraphStyle("Meta", parent=body_style, alignment=TA_CENTER,
                   fontSize=11, leading=18)))

story.append(Spacer(1, 1.5 * cm))
story.append(Paragraph("Executive Summary", h1_style))
story.append(Paragraph(
    "Two major simulations were run to stress-test all candidate improvements to "
    "the ReReReRe method for careless-respondent detection. The <b>All Variants</b> "
    "simulation (13 methods × 75 conditions) tested all three ideas from "
    "<i>buone_idee.md</i>: iterative EFA (idea 4) emerged as the best overall "
    "variant (+4% MCC over std), while cross-factor baseline (idea 2) showed no "
    "benefit and the per-factor aggregation schemes (idea 1) were substantially "
    "worse. A follow-up <b>Weighted V2</b> simulation (4 methods × 45 conditions) "
    "tested a proper weighted-Pearson formulation and confirmed it is inferior to "
    "the existing top three methods in most of the parameter space.",
    body_style))

story.append(PageBreak())

# --- Section 1: Architecture ---
story.append(Paragraph("1. ReReReRe Architecture", h1_style))
story.append(Paragraph(
    "Since 2026-04-01 the package exposes two functions for careless-respondent "
    "detection. Both return a per-respondent z-score (lower = more likely "
    "careless) computed by comparing observed coupled correlations against a "
    "permutation baseline.", body_style))

story.append(Paragraph("1.1 Standard function <i>ReReReRe()</i>", h2_style))
story.append(Paragraph(
    "Auto-switches between two scoring engines: <b>weighted</b> (all item pairs, "
    "each weighted by its sample |r|) for questionnaires ≤60 items, and "
    "<b>coupled</b> (top-k% pairs by |r|, equal-weighted) for longer "
    "questionnaires. This is the validated, recommended method for all primary "
    "analyses (AUC=0.635 on six external datasets).", body_style))

story.append(Paragraph("1.2 Per-factor function <i>ReReReRe_F()</i>", h2_style))
story.append(Paragraph(
    "Wrapper that forces <code>mode='efa_d'</code>: runs an EFA, selects only "
    "within-factor pairs (never across constructs), and weights them by observed "
    "|r|. Experimental — wins on simulated data and short real questionnaires "
    "(Pennycook), loses on external validation sets with noisy ground truth.",
    body_style))

story.append(Paragraph("1.3 Why two functions?", h2_style))
story.append(Paragraph(
    "EFA-based methods consistently win on simulated data but lose on real "
    "datasets with contaminated ground truth. Keeping both exposed covers both "
    "use cases: <i>ReReReRe()</i> for primary analyses, <i>ReReReRe_F()</i> as a "
    "secondary check, particularly on short questionnaires.", body_style))

story.append(PageBreak())

# --- Section 2: All Variants Simulation ---
story.append(Paragraph("2. All Variants Simulation (2026-04-02)", h1_style))
story.append(Paragraph(
    "Comprehensive test of every idea in <i>buone_idee.md</i>. 75 conditions "
    "(5 nF × 3 ipf × 5 reps), N=300, 15% careless, corruption levels 10-100%, "
    "ground truth = respondent with &gt;50% corruption. 13 methods compared. "
    "Runtime: 123 minutes on a local Windows machine.", body_style))

story.append(Paragraph("2.1 Methods tested", h2_style))
methods_tbl = [
    ["Method", "Description", "Source idea"],
    ["std", "Standard ReReReRe (auto weighted/coupled)", "—"],
    ["std_cf", "Standard + cross-factor permutation baseline", "Idea 2"],
    ["efa_d", "EFA within-factor pairs + |r| weights", "Existing"],
    ["efa_d_cf", "EFA-D + cross-factor permutation baseline", "Idea 2"],
    ["iterative", "EFA -> flag -> re-EFA on clean -> re-score", "Idea 4"],
    ["pf_mean", "Mean of per-factor z-scores", "Idea 1"],
    ["proplow (x4)", "% of factors with z < threshold", "Idea 1"],
    ["pf_comb (x3)", "mean(z) - lambda * sqrt(var(z))", "Idea 1"],
]
story.append(make_table(methods_tbl, col_widths=[3 * cm, 9 * cm, 3 * cm]))
story.append(Spacer(1, 0.4 * cm))

story.append(Paragraph("2.2 Overall Oracle MCC ranking", h2_style))
ranking_tbl = [
    ["#", "Method", "Oracle MCC"],
    ["1", "iterative", "0.398"],
    ["2", "efa_d", "0.394"],
    ["3", "efa_d_cf", "0.390"],
    ["4", "std", "0.383"],
    ["5", "std_cf", "0.380"],
    ["6", "pf_mean", "0.379"],
    ["7", "proplow (best)", "0.280"],
    ["8", "pf_comb (best)", "0.277"],
]
story.append(make_table(ranking_tbl,
                        col_widths=[1.5 * cm, 6 * cm, 3 * cm],
                        highlight_col=2))
story.append(Paragraph(
    "<i>Highlighted row: winner. Iterative EFA beats standard by +0.015 MCC "
    "(+4%).</i>", caption_style))

story.append(Paragraph("2.3 Performance by questionnaire length", h2_style))
bins_tbl = [
    ["Total items", "efa_d", "iterative", "std", "pf_mean", "proplow"],
    ["<30", "0.207", "0.196", "0.196", "0.198", "0.181"],
    ["30-60", "0.281", "0.279", "0.261", "0.258", "0.199"],
    ["60-100", "0.379", "0.407", "0.324", "0.355", "0.272"],
    ["100-200", "0.577", "0.580", "0.576", "0.571", "0.383"],
    [">200", "0.724", "0.725", "0.832", "0.705", "0.585"],
]
story.append(make_table(bins_tbl,
                        col_widths=[2.5 * cm, 2 * cm, 2 * cm, 2 * cm, 2 * cm, 2 * cm]))
story.append(Paragraph(
    "<i>efa_d wins on short questionnaires (&lt;60 items), iterative on medium "
    "(60-200), standard on long (&gt;200).</i>", caption_style))

story.append(PageBreak())

story.append(Paragraph("2.4 Iterative EFA — visual comparison", h2_style))
story.append(Paragraph(
    "The three plots below illustrate where iterative EFA provides a meaningful "
    "improvement over the existing methods. The delta plot is the most "
    "diagnostic: positive values mean iterative beats standard.",
    body_style))

for img_name, caption in [
    ("plot_iterative_vs_total_items.png",
     "Figure 1 — Oracle MCC vs total items. All four methods scale upward with "
     "more items. Standard (blue) takes the lead at ≥300 items; iterative (red) "
     "and efa_d (green) cluster together in the 60-200 range."),
    ("plot_iterative_vs_nF.png",
     "Figure 2 — Oracle MCC by number of factors, faceted by items-per-factor. "
     "With ipf=10 and nF=30, standard reaches MCC=0.83 — no EFA method catches "
     "up. At ipf=3 all methods stay flat (too few items)."),
    ("plot_iterative_gain.png",
     "Figure 3 — Delta MCC (iterative − standard). Peak gain of +0.07 at 72-80 "
     "total items. Gain is positive across 40-120 items, then collapses at 300 "
     "items (−0.10)."),
]:
    story.append(Image(rf"{ROOT}\{img_name}", width=15 * cm, height=9 * cm))
    story.append(Paragraph(caption, caption_style))

story.append(PageBreak())

story.append(Paragraph("2.5 Verdict per idea", h2_style))
verdict_tbl = [
    ["Idea", "Outcome", "Δ MCC vs std"],
    ["Cross-factor baseline (2)", "No benefit", "−0.003"],
    ["proplow: % factors z<threshold (1)", "Much worse", "−0.103"],
    ["pf_comb: mean - lambda*sqrt(var) (1)", "Very bad", "-0.106"],
    ["pf_mean: mean of per-factor z (1)", "Equivalent", "−0.004"],
    ["Iterative EFA (4)", "Winner", "+0.015"],
]
story.append(make_table(verdict_tbl,
                        col_widths=[6 * cm, 4 * cm, 3 * cm]))
story.append(Paragraph(
    "<i>Only iterative EFA provides a genuine improvement. Deltas averaged over "
    "75 conditions.</i>", caption_style))

story.append(PageBreak())

# --- Section 3: Weighted V2 ---
story.append(Paragraph("3. Weighted V2 Head-to-Head (2026-04-21)", h1_style))
story.append(Paragraph(
    "Weighted V2 is a proper weighted-Pearson correlation: both the mean and "
    "the variance used to center items are weighted by |r|. The original "
    "<i>weighted</i> mode used a simpler formulation (sum of weighted products). "
    "The V2 variant was pushed by a collaborator; we ran a focused simulation "
    "(45 conditions: 5 nF × 3 ipf × 3 reps, same parameters as §2) to decide "
    "whether to adopt it.", body_style))
story.append(Paragraph(
    "Only the three top performers from §2 (std, efa_d, iterative) were included "
    "as competitors. Runtime: 65 minutes.", body_style))

story.append(Paragraph("3.1 Overall mean Oracle MCC", h2_style))
ov_tbl = [
    ["std", "efa_d", "iterative", "weighted_v2"],
    ["0.345", "0.344", "0.342", "0.320"],
]
story.append(make_table(ov_tbl,
                        col_widths=[3 * cm, 3 * cm, 3 * cm, 3 * cm]))
story.append(Paragraph(
    "<i>Weighted V2 is ~0.025 MCC behind the top three, a 7% relative loss.</i>",
    caption_style))

story.append(Paragraph("3.2 Breakdown by questionnaire length", h2_style))
wv2_bins = [
    ["Total items", "std", "efa_d", "iterative", "weighted_v2", "Δ vs best"],
    ["<30", "0.170", "0.159", "0.150", "0.102", "−0.068"],
    ["30-60", "0.210", "0.212", "0.215", "0.204", "−0.011"],
    ["60-100", "0.299", "0.347", "0.337", "0.370", "+0.023"],
    ["100-200", "0.541", "0.516", "0.514", "0.487", "−0.054"],
    [">200", "0.768", "0.732", "0.752", "0.625", "−0.143"],
]
story.append(make_table(wv2_bins,
                        col_widths=[2.5 * cm, 2 * cm, 2 * cm, 2 * cm, 2.5 * cm, 2.3 * cm]))
story.append(Paragraph(
    "<i>Weighted V2 has a narrow win zone at 60-100 items (+0.023) and "
    "collapses at both extremes: −0.068 below 30 items, −0.143 above 200.</i>",
    caption_style))

story.append(Paragraph("3.3 Winners per cell (nF × ipf)", h2_style))
cells_tbl = [
    ["nF", "ipf", "items", "std", "efa_d", "iterative", "weighted_v2", "winner"],
    ["4", "3", "12", "0.145", "0.155", "0.141", "0.070", "efa_d"],
    ["4", "6", "24", "0.177", "0.169", "0.157", "0.151", "std"],
    ["4", "10", "40", "0.180", "0.204", "0.207", "0.210", "weighted_v2"],
    ["8", "3", "24", "0.188", "0.154", "0.152", "0.085", "std"],
    ["8", "6", "48", "0.244", "0.190", "0.207", "0.193", "std"],
    ["8", "10", "80", "0.374", "0.396", "0.396", "0.329", "efa_d"],
    ["12", "3", "36", "0.177", "0.145", "0.169", "0.156", "std"],
    ["12", "6", "72", "0.265", "0.359", "0.340", "0.383", "weighted_v2"],
    ["12", "10", "120", "0.454", "0.487", "0.495", "0.400", "iterative"],
    ["20", "3", "60", "0.239", "0.308", "0.279", "0.255", "efa_d"],
    ["20", "6", "120", "0.438", "0.441", "0.450", "0.405", "iterative"],
    ["20", "10", "200", "0.633", "0.578", "0.555", "0.621", "std"],
    ["30", "3", "90", "0.260", "0.286", "0.275", "0.398", "weighted_v2"],
    ["30", "6", "180", "0.640", "0.560", "0.557", "0.523", "std"],
    ["30", "10", "300", "0.768", "0.732", "0.752", "0.625", "std"],
]
cells_style = table_style()
for r in range(1, len(cells_tbl)):
    if cells_tbl[r][-1] == "weighted_v2":
        cells_style.add("BACKGROUND", (0, r), (-1, r),
                        colors.HexColor("#fff3cd"))
t = Table(cells_tbl,
          colWidths=[1.2 * cm, 1.2 * cm, 1.5 * cm, 1.8 * cm, 1.8 * cm,
                     2.2 * cm, 2.5 * cm, 2.5 * cm],
          repeatRows=1)
t.setStyle(cells_style)
story.append(t)
story.append(Paragraph(
    "<i>Highlighted cells: weighted_v2 wins 3/15 (40 items, 72 items, 90 "
    "items). Outside its narrow zone it is dominated.</i>", caption_style))

story.append(PageBreak())

# --- Conclusions ---
story.append(Paragraph("4. Conclusions and Recommendations", h1_style))

story.append(Paragraph("4.1 For the package default", h2_style))
story.append(Paragraph(
    "Keep the current two-function architecture unchanged. No candidate variant "
    "is worth promoting to the default:",
    body_style))
story.append(Paragraph(
    "• <b>Iterative EFA</b> is the best-ranked single variant but only by +4%, "
    "and its win zone (60-200 items) is already well covered by existing methods.<br/>"
    "• <b>Weighted V2</b> is strictly inferior to the top three across most of "
    "the parameter space and collapses on very short or very long questionnaires.<br/>"
    "• <b>Cross-factor baseline</b> adds complexity (mandatory EFA) without "
    "gain — reject.<br/>"
    "• <b>Per-factor aggregation schemes</b> (proplow, pf_comb) lose 25-30% "
    "MCC — reject.",
    body_style))

story.append(Paragraph("4.2 For advanced users", h2_style))
story.append(Paragraph(
    "Iterative EFA could be exposed as an opt-in flag on <i>ReReReRe_F()</i> "
    "for the 60-200 items use case. The current standalone function "
    "<i>ReReReRe_F_iterative()</i> is already available but undocumented.",
    body_style))

story.append(Paragraph("4.3 Fixed recommendations (unchanged)", h2_style))
fixed_tbl = [
    ["Parameter", "Value", "Rationale"],
    ["corProp", "0.03", "Full multiverse MCC 0.352 vs 0.327 at 0.05"],
    ["z_threshold", "1.5", "Robust universal default; auto_z equivalent"],
    ["align_signs", "TRUE", "Mandatory for reverse-coded items"],
    ["auto_z", "TRUE", "Convenience; looks up z from total_items LOESS"],
    ["iterations", "100", "Standard permutation count"],
]
story.append(make_table(fixed_tbl, col_widths=[3 * cm, 2 * cm, 10 * cm]))

story.append(Spacer(1, 0.6 * cm))
story.append(Paragraph(
    "<b>Bottom line:</b> after exhaustively testing every improvement candidate "
    "in <i>buone_idee.md</i> plus the Weighted V2 variant, none crosses the "
    "threshold for adoption as a new default. The two-function architecture "
    "with its current defaults (corProp=0.03, z=1.5, align_signs=TRUE) remains "
    "the recommended configuration.",
    body_style))

story.append(Spacer(1, 1 * cm))
story.append(Paragraph(
    "Generated 2026-04-21 from <i>sim_variants_results.csv</i> (12,525 rows) "
    "and <i>sim_weighted_v2_results.csv</i>. Source: "
    "github.com/vittorio-g/ReReReRe, commit f6b6dcf.",
    caption_style))

# ============================================================
# Build
# ============================================================
doc = SimpleDocTemplate(OUT, pagesize=A4,
                         leftMargin=2 * cm, rightMargin=2 * cm,
                         topMargin=2 * cm, bottomMargin=2 * cm,
                         title="ReReReRe Variants Report",
                         author="Vittorio Guerri")
doc.build(story)
print(f"Saved: {OUT}")
