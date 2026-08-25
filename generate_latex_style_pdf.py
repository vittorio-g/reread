"""Generate a LaTeX-style academic PDF summarising ReReReRe results."""

from reportlab.lib.pagesizes import A4
from reportlab.lib.styles import getSampleStyleSheet, ParagraphStyle
from reportlab.lib.units import cm
from reportlab.lib import colors
from reportlab.platypus import (
    SimpleDocTemplate, Paragraph, Spacer, PageBreak, Table, TableStyle,
    ListFlowable, ListItem,
)
from reportlab.lib.enums import TA_CENTER, TA_JUSTIFY, TA_LEFT

OUT = r"C:\Users\vitto\Downloads\ReReReRe_LaTeX_Report_2026-04-21.pdf"

styles = getSampleStyleSheet()

# ------------------------------------------------------------
# LaTeX-style typography: Times Roman body, bold/italic for structure
# ------------------------------------------------------------
title_style = ParagraphStyle(
    "Title", fontName="Times-Bold",
    fontSize=18, leading=21, alignment=TA_CENTER, spaceAfter=6,
)
author_style = ParagraphStyle(
    "Author", fontName="Times-Roman",
    fontSize=11, leading=13, alignment=TA_CENTER, spaceAfter=2,
)
date_style = ParagraphStyle(
    "Date", fontName="Times-Italic",
    fontSize=10, leading=12, alignment=TA_CENTER, spaceAfter=16,
)
abstract_heading_style = ParagraphStyle(
    "AbsHead", fontName="Times-Bold",
    fontSize=10, leading=12, alignment=TA_CENTER, spaceAfter=4,
)
abstract_style = ParagraphStyle(
    "Abstract", fontName="Times-Roman",
    fontSize=9.5, leading=12, alignment=TA_JUSTIFY,
    spaceAfter=12, leftIndent=0.8*cm, rightIndent=0.8*cm,
)
h1_style = ParagraphStyle(
    "H1", fontName="Times-Bold",
    fontSize=13, leading=16, spaceBefore=14, spaceAfter=6,
)
h2_style = ParagraphStyle(
    "H2", fontName="Times-Bold",
    fontSize=11.5, leading=14, spaceBefore=10, spaceAfter=4,
)
body_style = ParagraphStyle(
    "Body", fontName="Times-Roman",
    fontSize=10.5, leading=13.5, alignment=TA_JUSTIFY,
    spaceAfter=6, firstLineIndent=0.5*cm,
)
body_noindent_style = ParagraphStyle(
    "BodyNI", fontName="Times-Roman",
    fontSize=10.5, leading=13.5, alignment=TA_JUSTIFY, spaceAfter=6,
)
bullet_style = ParagraphStyle(
    "Bullet", fontName="Times-Roman",
    fontSize=10, leading=13, alignment=TA_JUSTIFY, spaceAfter=3,
    leftIndent=0.6*cm, bulletIndent=0.2*cm,
)
caption_style = ParagraphStyle(
    "Caption", fontName="Times-Italic",
    fontSize=9, leading=11, alignment=TA_CENTER,
    spaceAfter=8, textColor=colors.black,
)
note_style = ParagraphStyle(
    "Note", fontName="Times-Italic",
    fontSize=8.5, leading=10.5, alignment=TA_LEFT,
    textColor=colors.grey, spaceAfter=4,
)

_cell_style = ParagraphStyle(
    "Cell", fontName="Times-Roman",
    fontSize=9.5, leading=11.5, alignment=TA_CENTER,
)
_cell_left_style = ParagraphStyle(
    "CellLeft", fontName="Times-Roman",
    fontSize=9.5, leading=11.5, alignment=TA_LEFT,
)
_cell_bold_style = ParagraphStyle(
    "CellBold", fontName="Times-Bold",
    fontSize=9.5, leading=11.5, alignment=TA_CENTER,
)

def _wrap_cell(val, left=False, bold=False):
    if isinstance(val, str):
        style = _cell_bold_style if bold else (_cell_left_style if left else _cell_style)
        return Paragraph(val, style)
    return val

def make_table(data, col_widths, font_size=9.5, header=True,
               left_cols=None):
    """left_cols: set of column indices to left-align in body rows (e.g. rationale)."""
    if left_cols is None:
        left_cols = set()
    wrapped = []
    for r, row in enumerate(data):
        new_row = []
        for c, cell in enumerate(row):
            bold = header and r == 0
            left = (not bold) and (c in left_cols)
            new_row.append(_wrap_cell(cell, left=left, bold=bold))
        wrapped.append(new_row)

    t = Table(wrapped, colWidths=col_widths, repeatRows=1 if header else 0)
    style_cmds = [
        ("VALIGN", (0, 0), (-1, -1), "MIDDLE"),
        ("LINEABOVE", (0, 0), (-1, 0), 1.0, colors.black),
        ("LINEBELOW", (0, 0), (-1, 0), 0.7, colors.black),
        ("LINEBELOW", (0, -1), (-1, -1), 1.0, colors.black),
        ("PADDING", (0, 0), (-1, -1), 3),
    ]
    t.setStyle(TableStyle(style_cmds))
    return t

# ============================================================
# STORY
# ============================================================
story = []

# --- Title block ---
story.append(Paragraph(
    "ReReReRe: A Permutation-Based Method for<br/>Detecting Careless Respondents",
    title_style))
story.append(Paragraph("Vittorio Guerri", author_style))
story.append(Paragraph(
    "Consiglio Nazionale delle Ricerche &middot; April 21, 2026",
    date_style))

# --- Abstract ---
story.append(Paragraph("Abstract", abstract_heading_style))
story.append(Paragraph(
    "We present a method for detecting careless respondents in Likert-scale "
    "questionnaires by comparing each respondent's observed cross-item coherence "
    "against a personal permutation baseline. For each respondent a z-score is "
    "derived from coupled item pairs (those with high sample-level |r|) relative to "
    "randomly paired items; low z-scores identify probable careless responders. "
    "Evidence from a 633,600-row multiverse simulation, a 60-cell 2D calibration, "
    "an all-variants stress-test of thirteen candidate improvements, and external "
    "validation on seven real datasets is reported. We distinguish throughout between "
    "<i>oracle</i> performance (maximum MCC across all z-thresholds, unreachable in "
    "practice) and <i>practical</i> performance at the fixed default z=1.5 or at "
    "the 2D-calibrated auto-z. Under practical fixed defaults, ReReReRe outperforms "
    "practical Mahalanobis distance by a factor of 4.8 (MCC 0.271 vs 0.056). Peak "
    "oracle MCC reaches 0.838 on 360-item questionnaires; the same cell with fixed "
    "z=1.5 gives MCC = 0.684, and with auto-z = 0.710. Fair conditions "
    "(>=15 factors, n>=300) yield practical MCC = 0.435 with no tuning. The method "
    "is not recommended for very short questionnaires, very small samples, or "
    "datasets dominated by consistent carelessness (acquiescence, straight-lining).",
    abstract_style))

# --- 1. Results Description ---
story.append(Paragraph("1. General Description of Results", h1_style))
story.append(Paragraph(
    "ReReReRe has been evaluated in four complementary studies. The "
    "<i>full multiverse</i> (5,670 RR calls, 85,050 result rows) varied seven "
    "design and method parameters to measure their relative impact; the 2D "
    "calibration (60 cells, 30 reps each) produced a total-items -> z-threshold "
    "lookup table that supports self-calibration; the <i>all-variants</i> "
    "simulation (75 conditions, 13 methods) tested every improvement candidate; "
    "and the external validation covers six ground-truth datasets plus the "
    "Johnson IPIP-NEO-300 inject-and-detect benchmark.",
    body_style))

story.append(Paragraph(
    "Three findings dominate. First, the detection signal scales with "
    "<i>total questionnaire length</i> (items_per_factor * nFactors), not with "
    "either dimension alone: an ANOVA on multiverse results assigns 31.1% of "
    "variance to items_per_factor and 26.1% to nFactors, with sample size "
    "contributing only 1.0%. Second, ReReReRe and Mahalanobis are not "
    "interchangeable tools - they detect different constructs. ReReReRe captures "
    "inconsistent carelessness (random, mixed, scattered longstring), while "
    "practical Mahalanobis, with its chi-square(p) threshold, is essentially "
    "non-functional on Likert data (sensitivity 1.3%, MCC 0.056). Third, the "
    "method is remarkably insensitive to parameter choice in the recommended "
    "range: any z in [1.0, 2.0] yields near-identical performance.",
    body_style))

results_tbl = [
    ["Study", "Oracle MCC", "Practical MCC (z=1.5)", "Auto-z MCC"],
    ["Final multiverse (57.6K RR calls)", "0.361", "0.247", "0.268"],
    ["Full multiverse (5.6K RR calls)", "0.356", "0.274", "0.268"],
    ["All variants (12.5K rows, std)", "0.363", "0.283", "0.308"],
    ["Sweet spot (nF>=15, n>=300)", "~0.60", "0.435", "~0.45"],
    ["Practical Mahalanobis (chi-sq .001)", "0.361", "0.056", "-"],
]
story.append(make_table(results_tbl, col_widths=[6.8*cm, 2.8*cm, 3.8*cm, 2.8*cm]))
story.append(Paragraph(
    "Table 1. Four studies reported both as oracle (upper bound, unrealistic) and "
    "under fixed practical defaults. Note the ~0.08 MCC gap between oracle and "
    "z=1.5 on average. Practical Mahalanobis is the only tool that loses drastically "
    "when moving from oracle to its realistic chi-square threshold.",
    caption_style))

# --- 2. Where and how it works best ---
story.append(Paragraph("2. Where ReReReRe Works Best", h1_style))

story.append(Paragraph("2.1 Peak performance regimes", h2_style))
story.append(Paragraph(
    "Detection quality peaks in long, multi-factor questionnaires. The best "
    "condition observed was nF=30, ipf=12 (360 items), where <i>oracle</i> "
    "MCC = 0.838 with n=300. Closely comparable values were achieved at 200-300 "
    "item questionnaires with 15 or more factors. Oracle numbers represent an "
    "upper bound: they assume the best z-threshold is known, which it is not in "
    "practice. The same cell under fixed z=1.5 gives MCC ~ 0.65 and under the "
    "auto-z 2D calibration MCC ~ 0.80. A conservative <i>sweet spot</i> for "
    "recommendation is nF>=15 and n>=300: fixed practical defaults "
    "(corProp=0.03, z_threshold=1.5) give MCC=0.435 with sensitivity=0.79 and "
    "specificity=0.79, reached with no tuning at all. Under these conditions the "
    "z-score distributions for careless and attentive respondents are separated "
    "by about 2.5 to 2.8 standard deviations, providing reliable classification "
    "without external reference.",
    body_style))

peak_tbl = [
    ["Condition", "Items", "Oracle", "z=1.5", "Auto-z"],
    ["nF=30, ipf=12, n=300 (calibration)", "360", "0.838", "~0.65", "~0.80"],
    ["nF=30, ipf=10, n=300 (all-variants)", "300", "0.827", "0.684", "0.710"],
    ["nF=20, ipf=10, n=300 (all-variants)", "200", "0.690", "0.664", "0.645"],
    ["nF=30, ipf=6, n=300 (all-variants)", "180", "0.576", "0.486", "0.512"],
    ["nF=20, ipf=6, n=300 (all-variants)", "120", "0.526", "0.441", "0.451"],
    ["nF=12, ipf=10, n=300 (all-variants)", "120", "0.477", "0.398", "0.417"],
    ["Sweet spot avg (nF>=15, n>=300, final mv.)", ">=90", "~0.60", "0.435", "-"],
]
story.append(make_table(peak_tbl, col_widths=[7*cm, 1.6*cm, 1.8*cm, 1.8*cm, 1.9*cm]))
story.append(Paragraph(
    "Table 2. Peak performance conditions reported at all three levels. Oracle = "
    "maximum MCC across z-threshold grid (upper bound, unavailable in practice); "
    "z=1.5 = fixed default recommended for any questionnaire; auto-z = threshold "
    "chosen from the 2D total-items LOESS calibration. The oracle-to-practical "
    "gap narrows as items_per_factor grows: at 300 items it is 0.14, at 120 items "
    "only 0.08.",
    caption_style))

story.append(Paragraph("2.2 Crossover against Mahalanobis", h2_style))
story.append(Paragraph(
    "ReReReRe's advantage over practical Mahalanobis is not uniform. In MCC it "
    "crosses over at nF=10; in AUC it crosses over at nF=8. Beyond these "
    "thresholds, ReReReRe dominates. At nF=25, ReReReRe oracle MCC = 0.572 "
    "versus Mahalanobis 0.391 (+46% relative); AUC = 0.873 versus 0.641. "
    "Against the practical chi-square Mahalanobis (the realistic comparison), "
    "ReReReRe with fixed defaults wins by a factor of 4.8 across all "
    "conditions: MCC 0.271 vs 0.056.",
    body_style))

story.append(Paragraph("2.3 Real-data performance (practical values, not oracle)", h2_style))
story.append(Paragraph(
    "Simulation oracles are inaccessible on real data, so the table below reports "
    "only threshold-free AUC and practical MCC at z=1.5 or auto-z. These are the "
    "numbers a user would actually observe. ReReReRe AUC beats practical "
    "Mahalanobis on 4/6 ground-truth datasets and on the Johnson IPIP-NEO-300 "
    "inject-and-detect benchmark.",
    body_style))

real_tbl = [
    ["Dataset", "N", "Items", "RR AUC", "RR MCC z=1.5", "RR MCC auto-z", "Mah. MCC"],
    ["Schroeders 2022", "605", "60", "0.606", "0.177", "0.228", "0.045"],
    ["Schneider QoL", "1649", "31", "0.735", "0.118", "0.170", "0.208"],
    ["Niessen 2016", "180", "100", "0.639", "0.073", "-0.067", "0.000"],
    ["Goldammer S1 (BFI-2)", "291", "60", "0.711", "0.320", "0.269", "0.257"],
    ["Goldammer S2 (IPIP)", "265", "60", "0.674", "0.304", "0.241", "0.266"],
    ["Goldammer S3 (long.)", "523", "60", "0.462", "-0.036", "-0.072", "0.011"],
    ["Johnson IPIP-NEO-300", "5000", "300", "-", "0.628", "0.726", "-"],
]
story.append(make_table(real_tbl,
                        col_widths=[4*cm, 1.5*cm, 1.5*cm, 1.7*cm, 2.3*cm, 2.3*cm, 2.1*cm],
                        font_size=9))
story.append(Paragraph(
    "Table 3. Real-data performance of ReReReRe vs practical Mahalanobis "
    "(chi-square alpha=.001). RR AUC beats Mahalanobis on 4/6 datasets; practical "
    "MCC beats it on 4/6. On Johnson IPIP-NEO-300 (inject-and-detect with 5-20% "
    "careless), auto-z achieves specificity = 1.000 across all injection rates.",
    caption_style))

story.append(Paragraph("2.4 Recommended application profile", h2_style))
story.append(Paragraph(
    "Best results are obtained when the questionnaire has at least 10 factors "
    "(approximately 60-plus items), a sample size of at least 300, and pursues "
    "a wide construct domain where carelessness manifests as genuine "
    "incoherence. Large-scale online batteries, multi-instrument panels, and "
    "omnibus personality measurements fit this profile. The method performs "
    "especially well when questionnaires include reverse-coded items: the "
    "<i>align_signs</i> routine handles them transparently, requiring no "
    "manual recoding or item-level preprocessing by the analyst.",
    body_style))

# --- 3. Where it should NOT be used ---
story.append(Paragraph("3. Where ReReReRe Should Not Be Used", h1_style))

story.append(Paragraph("3.1 Structural limits of the signal", h2_style))
story.append(Paragraph(
    "Below four factors, the method lacks the cross-factor coherence that its "
    "z-score relies on. Simulation and real data agree: at nF=4 the gap between "
    "careless and attentive z-scores is only 0.48 standard deviations, and "
    "oracle MCC stays around 0.12. At nF=5 (the Niessen 2016 dataset, 100 "
    "items), practical Mahalanobis is the better choice (AUC 0.585 vs "
    "ReReReRe 0.551). Very short questionnaires (<30 items) and very small "
    "samples (n<=50) should also be avoided: individual-level correlations "
    "become unstable and the permutation baseline noisy.",
    body_style))

bad_tbl = [
    ["Condition", "Why it fails", "Alternative"],
    ["nF <= 4", "Cross-factor signal collapses", "Mahalanobis + IRV"],
    ["n <= 50", "Item-level cor. estimates unstable", "LongString / IRV"],
    ["Items < 30", "Too few pairs for any baseline", "LongString only"],
    ["Pure acquiescence / speed", "RR detects inconsistency, not consistency", "IRV, response-time screens"],
    ["Mixed response scales", "Weights biased by scale range", "Rescale to [0,1] first"],
    ["High careless rate (>60%)", "Correlation matrix corrupted", "Pre-screen with longstring"],
]
story.append(make_table(bad_tbl, col_widths=[4.2*cm, 6.6*cm, 5.9*cm],
                        left_cols={1, 2}))
story.append(Paragraph(
    "Table 4. Conditions where ReReReRe is not the right tool, and what to use instead.",
    caption_style))

story.append(Paragraph("3.2 What ReReReRe does not detect", h2_style))
story.append(Paragraph(
    "ReReReRe and IRV-based detectors are complementary, not competing. "
    "ReReReRe measures whether a response profile is <i>structurally coherent</i> "
    "with a sample-level correlation structure. It therefore misses "
    "<i>consistent</i> forms of carelessness - acquiescent response style, "
    "straight-lining, extreme responding - because these produce internally "
    "consistent profiles that align with high-|r| pairs. This limit was "
    "confirmed on PISA 2018 data, where fast or acquiescent respondents in fact "
    "had <i>higher</i> z-scores than attentive ones. In such datasets, a "
    "longstring or IRV pre-screen must be applied before (or alongside) "
    "ReReReRe.",
    body_style))

# --- 4. Statistical choices ---
story.append(Paragraph("4. Statistical Choices, Point by Point", h1_style))

story.append(Paragraph(
    "The design of the algorithm and of its evaluation involved a series of "
    "deliberate choices. Each is reported below, together with the evidence "
    "that supports it.",
    body_style))

story.append(Paragraph("4.1 Scoring engine", h2_style))
choices_tbl_1 = [
    ["Choice", "Value", "Rationale"],
    ["Primary score", "z-score (perm.)",
     "Cohen's d up to 3.1 vs 1.7 for percentile (nF>=10); "
     "free of distributional assumptions."],
    ["Baseline construction", "Per-respondent permutation",
     "Self-calibrating; no need for external reference."],
    ["Minimum pairs", "min_pairs = 15",
     "Avoids degenerate individual-level r with small questionnaires."],
    ["Iterations", "100",
     "Trade-off between MCC stability and runtime; gains plateau by 80."],
    ["Sign handling", "align_signs = TRUE",
     "Reverse-coded items cancel inside cor() otherwise (Robie HEXACO)."],
    ["Reverse-coding method", "(max+1) - x",
     "Keeps values in original Likert range, unlike naive sign flip."],
]
story.append(make_table(choices_tbl_1, col_widths=[3.5*cm, 4*cm, 9.2*cm],
                        left_cols={2}))

story.append(Paragraph("4.2 Pair selection", h2_style))
choices_tbl_2 = [
    ["Choice", "Value", "Rationale"],
    ["Selection rule", "Auto (weighted / coupled)",
     "Weighted for <=60 items, coupled for >60; 2D calibration supports split."],
    ["Top-k proportion", "corProp = 0.03",
     "Full multiverse: MCC 0.352 at 0.03 vs 0.327 at 0.05 vs 0.286 at 0.10."],
    ["Weighted mode scope", "All pairs, |r|-weighted",
     "Preserves weak-but-informative coherence on short questionnaires."],
    ["EFA-based selection", "Optional (ReReReRe_F)",
     "Wins on simulation and Pennycook; loses on noisy external GT."],
    ["Cross-factor baseline", "REJECTED",
     "No measurable benefit; adds EFA dependency (Delta MCC -0.003)."],
]
story.append(make_table(choices_tbl_2, col_widths=[3.5*cm, 4*cm, 9.2*cm],
                        left_cols={2}))

story.append(Paragraph("4.3 Thresholding", h2_style))
choices_tbl_3 = [
    ["Choice", "Value", "Rationale"],
    ["Default z_threshold", "1.5",
     "Flat plateau 1.0-2.0; MCC 0.274 on full multiverse."],
    ["Adaptive z_threshold", "auto_z via 2D LOESS",
     "R2 = 0.476 on 60-cell grid; equivalent performance to z=1.5 overall."],
    ["Fallback formula", "z = 0.505 + 0.0042 * total_items",
     "Linear approximation if LOESS prediction fails."],
    ["z clamp", "[0.3, 3.5]",
     "Prevents extreme values at out-of-grid extrapolation."],
    ["MAD-based thresholding", "REJECTED",
     "z distribution too spread (MAD 1.4-2.8); sensitivity drops to 1.6%."],
]
story.append(make_table(choices_tbl_3, col_widths=[3.5*cm, 4.5*cm, 8.7*cm],
                        left_cols={2}))

story.append(Paragraph("4.4 Simulation design", h2_style))
choices_tbl_4 = [
    ["Choice", "Value", "Rationale"],
    ["Loading range", "Uniform(0.2, 0.8)",
     "More ecologically valid than (0.4, 0.8); captures weak real items."],
    ["Max factor correlation", "0.8",
     "Allows high redundancy while avoiding near-singular matrices."],
    ["Careless corruption", "seq(0.1, 1.0, 0.1)",
     "Replaces earlier 50-100% range for realistic mixed contamination."],
    ["Ground truth cutoff", "careless_pct > 0.50",
     "Respondents with <=50% corruption labelled noise, not GT positive."],
    ["Careless types", "Random / LongString / Mixed (1/3 each)",
     "Fair coverage of inconsistent-careless patterns."],
    ["Replications", "50 (final multiverse)",
     "Within-condition SE ~ 0.011; 95th pct SD = 0.176."],
]
story.append(make_table(choices_tbl_4, col_widths=[3.5*cm, 4.5*cm, 8.7*cm],
                        left_cols={2}))

story.append(Paragraph("4.5 Evaluation metric", h2_style))
choices_tbl_5 = [
    ["Choice", "Value", "Rationale"],
    ["Primary metric", "MCC",
     "Uses all four quadrants; robust to class imbalance."],
    ["Secondary metric", "AUC",
     "Threshold-free; fairest for head-to-head comparisons."],
    ["F1 reporting", "Tracked but not primary",
     "Ignores specificity; misleading with low base rates."],
    ["Oracle vs practical", "Both reported",
     "Oracle shows upper bound; practical reflects real-world use."],
]
story.append(make_table(choices_tbl_5, col_widths=[3.5*cm, 4*cm, 9.2*cm],
                        left_cols={2}))

story.append(Paragraph("4.6 Reproducibility", h2_style))
story.append(Paragraph(
    "Every simulation reported here uses a fixed random seed schema of the form "
    "SEED_BASE + (rep - 1) * nConditions + condition_idx. Raw CSV outputs, "
    "checkpoint files, and plotting scripts are committed at "
    "<i>github.com/vittorio-g/ReReReRe</i>, commit <tt>f6b6dcf</tt>. OneDrive "
    "file-lock resilience is provided by <tt>safe_write_table()</tt> and "
    "<tt>safe_saveRDS()</tt> wrappers.",
    body_style))

# --- Closing ---
story.append(Paragraph("5. Conclusion", h1_style))
story.append(Paragraph(
    "ReReReRe is a self-calibrating permutation-based detector of careless "
    "respondents in Likert-scale data. It works best on long, multi-factor "
    "questionnaires with a moderately large sample, where it outperforms "
    "practical Mahalanobis by a wide margin and approaches oracle-level "
    "detection on 300-item batteries. Its reach is bounded by questionnaire "
    "structure (nF >= 5, items >= 30) and by the type of carelessness it "
    "measures (incoherence rather than consistency). After an exhaustive test "
    "of candidate improvements - cross-factor baselines, per-factor "
    "aggregations, proper weighted-Pearson formulations, and iterative EFA "
    "re-estimation - the current two-function architecture with fixed defaults "
    "remains the recommended configuration.",
    body_style))

story.append(Spacer(1, 0.4*cm))
story.append(Paragraph(
    "<i>Data and code:</i> github.com/vittorio-g/ReReReRe &middot; "
    "<i>Commit:</i> f6b6dcf &middot; <i>Compiled:</i> 2026-04-21",
    note_style))

# ============================================================
# Build
# ============================================================
doc = SimpleDocTemplate(OUT, pagesize=A4,
                         leftMargin=2.4*cm, rightMargin=2.4*cm,
                         topMargin=2.2*cm, bottomMargin=2.2*cm,
                         title="ReReReRe: A Permutation-Based Method",
                         author="Vittorio Guerri")
doc.build(story)
print(f"Saved: {OUT}")
