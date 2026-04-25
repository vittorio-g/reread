"""Generate a step-by-step methodological document explaining the statistical
procedure behind ReReReRe."""

from reportlab.lib.pagesizes import A4
from reportlab.lib.styles import ParagraphStyle
from reportlab.lib.units import cm
from reportlab.lib import colors
from reportlab.platypus import (
    SimpleDocTemplate, Paragraph, Spacer, Table, TableStyle, PageBreak,
    KeepTogether, Image
)
from reportlab.lib.enums import TA_CENTER, TA_JUSTIFY, TA_LEFT

OUT = r"C:\Users\vitto\Downloads\ReReReRe_Methodology_2026-04-25.pdf"

# ---- Styles ----
title_style = ParagraphStyle("Title", fontName="Times-Bold",
    fontSize=20, leading=24, alignment=TA_CENTER, spaceAfter=8)
sub_style = ParagraphStyle("Sub", fontName="Times-Italic",
    fontSize=11, leading=14, alignment=TA_CENTER, spaceAfter=18, textColor=colors.grey)
toc_style = ParagraphStyle("TOC", fontName="Times-Roman",
    fontSize=10.5, leading=14, alignment=TA_LEFT, leftIndent=0.5*cm, spaceAfter=2)
h1 = ParagraphStyle("H1", fontName="Times-Bold",
    fontSize=14, leading=17, spaceBefore=18, spaceAfter=8,
    textColor=colors.HexColor("#1a3a6c"))
h2 = ParagraphStyle("H2", fontName="Times-Bold",
    fontSize=12, leading=15, spaceBefore=12, spaceAfter=5,
    textColor=colors.HexColor("#2c5282"))
h3 = ParagraphStyle("H3", fontName="Times-Bold-Italic",
    fontSize=11, leading=14, spaceBefore=8, spaceAfter=4)
body = ParagraphStyle("Body", fontName="Times-Roman",
    fontSize=10.5, leading=14, alignment=TA_JUSTIFY,
    spaceAfter=7, firstLineIndent=0.5*cm)
body_no_indent = ParagraphStyle("BodyNI", fontName="Times-Roman",
    fontSize=10.5, leading=14, alignment=TA_JUSTIFY, spaceAfter=7)
formula_style = ParagraphStyle("Formula", fontName="Times-Italic",
    fontSize=11, leading=14, alignment=TA_CENTER,
    spaceAfter=10, spaceBefore=8, leftIndent=1*cm, rightIndent=1*cm)
numbered_step = ParagraphStyle("Step", fontName="Times-Roman",
    fontSize=10.5, leading=14, alignment=TA_LEFT,
    spaceAfter=4, leftIndent=0.8*cm, bulletIndent=0.2*cm)
caption = ParagraphStyle("Cap", fontName="Times-Italic",
    fontSize=9, leading=11, alignment=TA_CENTER, spaceAfter=10)
note = ParagraphStyle("Note", fontName="Times-Italic",
    fontSize=8.5, leading=10.5, alignment=TA_LEFT,
    textColor=colors.grey, spaceAfter=4)

cell = ParagraphStyle("Cell", fontName="Times-Roman",
    fontSize=9.5, leading=11.5, alignment=TA_CENTER)
cell_left = ParagraphStyle("CellL", fontName="Times-Roman",
    fontSize=9.5, leading=11.5, alignment=TA_LEFT)
cell_bold = ParagraphStyle("CellB", fontName="Times-Bold",
    fontSize=9.5, leading=11.5, alignment=TA_CENTER)

def wrap(val, left=False, bold=False):
    if isinstance(val, str):
        s = cell_bold if bold else (cell_left if left else cell)
        return Paragraph(val, s)
    return val

def mktable(data, cols, left_cols=None):
    if left_cols is None: left_cols = set()
    wrapped = [[wrap(c, left=(r > 0 and i in left_cols), bold=(r == 0))
                for i, c in enumerate(row)] for r, row in enumerate(data)]
    t = Table(wrapped, colWidths=cols, repeatRows=1)
    t.setStyle(TableStyle([
        ("VALIGN", (0, 0), (-1, -1), "MIDDLE"),
        ("LINEABOVE", (0, 0), (-1, 0), 1.0, colors.black),
        ("LINEBELOW", (0, 0), (-1, 0), 0.7, colors.black),
        ("LINEBELOW", (0, -1), (-1, -1), 1.0, colors.black),
        ("PADDING", (0, 0), (-1, -1), 3),
    ]))
    return t

def code_block(text):
    return Paragraph(
        f"<font face='Courier' size='9'>{text}</font>",
        ParagraphStyle("Code", fontName="Courier", fontSize=9, leading=11,
                       leftIndent=1*cm, rightIndent=1*cm,
                       spaceAfter=8, spaceBefore=4,
                       backColor=colors.HexColor("#f0f0f0"))
    )

# =====================================================================
# Story
# =====================================================================
story = []

# Cover
story.append(Spacer(1, 1.5*cm))
story.append(Paragraph(
    "ReReReRe<br/>"
    "Statistical Procedure, Step by Step",
    title_style))
story.append(Paragraph(
    "A methodological guide to permutation-based<br/>"
    "careless-respondent detection",
    sub_style))
story.append(Paragraph(
    "Vittorio Guerri &middot; April 25, 2026",
    ParagraphStyle("Author", parent=sub_style, textColor=colors.black,
                   fontSize=10)))

story.append(Spacer(1, 1*cm))
story.append(Paragraph("Contents", h1))
toc = [
    "1. Conceptual foundation",
    "2. Notation",
    "3. Step 1 — Compute the sample correlation matrix",
    "4. Step 2 — Select the coupled pairs (top-k% by |r|)",
    "5. Step 3 — Sign alignment for reverse-coded items",
    "6. Step 4 — Compute per-respondent observed coherence",
    "7. Step 5 — Build the per-respondent permutation baseline",
    "8. Step 6 — Convert to a z-score",
    "9. Step 7 — Optional variance penalty",
    "10. Step 8 — Iterative EFA refinement",
    "11. Step 9 — Auto-z calibration of the threshold",
    "12. Step 10 — Multi-detector ensemble",
    "13. Step 11 — Random Forest combiner",
    "14. Step 12 — Operating-point calibration (FPR=5%)",
    "15. Worked example — small numerical illustration",
    "16. Pseudo-code summary",
]
for line in toc:
    story.append(Paragraph(line, toc_style))

story.append(PageBreak())

# =====================================================================
# 1. Conceptual foundation
# =====================================================================
story.append(Paragraph("1. Conceptual foundation", h1))
story.append(Paragraph(
    "Likert-scale questionnaires often contain item pairs that are highly "
    "correlated at the sample level — typically because both items measure "
    "the same latent construct. A respondent who answers the questionnaire "
    "with attention should produce a response profile that is internally "
    "consistent with this sample-level structure: high values on one item "
    "tend to co-occur with high values on its strongly-correlated partner.",
    body))
story.append(Paragraph(
    "<b>The intuition behind ReReReRe is that this internal consistency "
    "becomes a per-respondent diagnostic when compared against a "
    "self-permuted baseline.</b> If a respondent's profile is structurally "
    "coherent on the genuinely high-correlation pairs but shows no "
    "particular coherence on randomly chosen pairs, they are likely "
    "attentive. If their profile is equally (in)coherent on both, the "
    "sample-level structure conveys no information for them — they are "
    "likely careless.",
    body))
story.append(Paragraph(
    "The procedure below converts this intuition into a formal statistic "
    "in twelve steps. The first six produce the core z-score; steps 7-9 "
    "are optional refinements; steps 10-12 wrap the score in a "
    "supervised-learning ensemble for publication-grade accuracy.",
    body))

# =====================================================================
# 2. Notation
# =====================================================================
story.append(Paragraph("2. Notation", h1))
notation_tbl = [
    ["Symbol", "Meaning"],
    ["X", "N × J data matrix (rows = respondents, columns = items)"],
    ["x_{ij}", "value of item j for respondent i"],
    ["x_i", "response vector of respondent i (length J)"],
    ["N",  "number of respondents"],
    ["J",  "number of items"],
    ["R",  "J × J sample correlation matrix"],
    ["k",  "number of coupled pairs selected"],
    ["P",  "set of coupled pairs {(a, b)}"],
    ["P*", "set of random pairs sampled at one permutation iteration"],
    ["B",  "number of permutation iterations (default 100)"],
    ["c_i^obs", "observed coherence score for respondent i"],
    ["c_i^*",   "random coherence score for respondent i (one iteration)"],
    ["μ_i, σ_i", "mean and SD of c_i^* across the B iterations"],
    ["z_i", "final z-score: (c_i^obs − μ_i) / σ_i"],
]
story.append(mktable(notation_tbl, [2.5*cm, 13*cm], left_cols={1}))

# =====================================================================
# 3. Step 1 — Correlation matrix
# =====================================================================
story.append(Paragraph("3. Step 1 — Compute the sample correlation matrix",
                        h1))
story.append(Paragraph(
    "The first step is computing the J × J sample correlation matrix R from "
    "the data matrix X. We use Pearson correlation with pairwise-complete "
    "observations:",
    body))
story.append(Paragraph(
    "R<sub>jk</sub> = corr(X<sub>·j</sub>, X<sub>·k</sub>)  for all j, k ∈ {1, …, J}",
    formula_style))
story.append(Paragraph(
    "We retain both the signed correlation matrix R (used to handle reverse-"
    "coded items in Step 3) and the matrix |R| of absolute values, used to "
    "identify the strongly-correlated pairs in Step 2. Only the lower-"
    "triangular portion of R is used downstream to avoid double-counting.",
    body))
story.append(Paragraph(
    "<b>Why this works:</b> the sample-level correlation captures the "
    "expected pattern of co-occurrence under the assumption that most "
    "respondents are attentive. Even with up to 40-50% careless contamination, "
    "the matrix is dominated by the attentive majority and reliably ranks "
    "true within-construct pairs above random pairs.",
    body))

# =====================================================================
# 4. Step 2 — Coupled pairs
# =====================================================================
story.append(Paragraph("4. Step 2 — Select the coupled pairs (top-k% by |r|)",
                        h1))
story.append(Paragraph(
    "From the absolute-correlation matrix |R|, we identify the top "
    "<i>corProp</i> fraction of pairs with the highest |r|. The default "
    "<i>corProp</i> = 0.03 was chosen via multiverse simulation: smaller "
    "fractions select fewer but more selective pairs, producing a cleaner "
    "signal.",
    body))
story.append(Paragraph(
    "P = {(a, b) : a &lt; b, |R<sub>ab</sub>| ≥ Q<sub>1−corProp</sub>(|R|)}",
    formula_style))
story.append(Paragraph(
    "where Q<sub>q</sub>(·) denotes the q-th quantile. Let k = |P|.",
    body_no_indent))
story.append(Paragraph(
    "<b>Safeguard:</b> if k &lt; min_pairs (default 15), the threshold is "
    "lowered until at least min_pairs are included. This avoids degenerate "
    "individual-level correlations on very short or weakly-structured "
    "questionnaires.",
    body))
story.append(Paragraph(
    "<b>Mode auto-switch:</b> for short questionnaires (≤ 60 items) "
    "ReReReRe uses an alternative \"weighted\" engine that includes ALL "
    "item pairs, weighted by |R|. The selection logic above is the "
    "\"coupled\" engine, used by default for longer questionnaires.",
    body))

# =====================================================================
# 5. Step 3 — Sign alignment
# =====================================================================
story.append(Paragraph("5. Step 3 — Sign alignment for reverse-coded items",
                        h1))
story.append(Paragraph(
    "Some pairs in P may have negative sample correlations, typically "
    "because one item is reverse-coded relative to its partner. If left "
    "unaligned, such pairs cancel against positively-correlated pairs in "
    "the per-respondent coherence computation, driving the score toward 0 "
    "for everyone.",
    body))
story.append(Paragraph(
    "For each pair (a, b) ∈ P with R<sub>ab</sub> &lt; 0, we reverse the "
    "second item's values:",
    body))
story.append(Paragraph(
    "x'<sub>ib</sub> = (max<sub>b</sub> + 1) − x<sub>ib</sub>",
    formula_style))
story.append(Paragraph(
    "where max<sub>b</sub> is the maximum observed value of item b. This "
    "preserves the Likert range (e.g., 1-7 stays 1-7) and inverts the "
    "ordering. After this transformation all pairs in P have positive "
    "expected correlation, which is necessary for the row-correlation in "
    "Step 4 to behave correctly.",
    body))
story.append(Paragraph(
    "The same alignment is applied to the random pairs sampled in Step 5, "
    "to keep the observed and random coherence scores on the same scale.",
    body))

# =====================================================================
# 6. Step 4 — Per-respondent observed coherence
# =====================================================================
story.append(Paragraph("6. Step 4 — Compute per-respondent observed coherence",
                        h1))
story.append(Paragraph(
    "For each respondent i, we compute a single number c<sub>i</sub><sup>obs</sup> "
    "summarizing how internally consistent their responses are with the "
    "structure imposed by P. We use the absolute value of the row-correlation "
    "between two derived vectors A<sub>i</sub> and B<sub>i</sub>:",
    body))
story.append(Paragraph(
    "A<sub>i</sub> = (x<sub>i, a<sub>1</sub></sub>, x<sub>i, a<sub>2</sub></sub>, …, "
    "x<sub>i, a<sub>k</sub></sub>)",
    formula_style))
story.append(Paragraph(
    "B<sub>i</sub> = (x'<sub>i, b<sub>1</sub></sub>, x'<sub>i, b<sub>2</sub></sub>, …, "
    "x'<sub>i, b<sub>k</sub></sub>)",
    formula_style))
story.append(Paragraph(
    "c<sub>i</sub><sup>obs</sup> = |corr(A<sub>i</sub>, B<sub>i</sub>)|",
    formula_style))
story.append(Paragraph(
    "Concretely: A<sub>i</sub> stacks the respondent's values on the "
    "<i>first</i> element of each coupled pair; B<sub>i</sub> stacks the "
    "<i>second</i> element (already sign-aligned). The Pearson correlation "
    "between these two vectors of length k is computed for that single "
    "respondent. High values mean the respondent's profile mirrors the "
    "pattern that the sample-level pairs imply; low values mean it does not.",
    body))
story.append(Paragraph(
    "<b>Edge cases:</b> if A<sub>i</sub> or B<sub>i</sub> has zero variance "
    "(e.g., a pure straight-liner whose values are identical) the "
    "correlation is undefined and is set to 0.",
    body))

# =====================================================================
# 7. Step 5 — Permutation baseline
# =====================================================================
story.append(Paragraph("7. Step 5 — Build the per-respondent permutation baseline",
                        h1))
story.append(Paragraph(
    "The observed coherence c<sub>i</sub><sup>obs</sup> alone is not "
    "enough — it does not say whether this value is high or low <i>for "
    "this respondent</i>, given their idiosyncratic response pattern. "
    "To answer that we build a personal null distribution.",
    body))
story.append(Paragraph(
    "We repeat B times (default B = 100):",
    body))
story.append(Paragraph(
    "1. Sample k random pairs P* uniformly from all (J choose 2) item pairs, "
    "ensuring no self-pairs.<br/>"
    "2. Apply the same sign alignment used for the coupled pairs.<br/>"
    "3. Compute c<sub>i</sub><sup>*</sup> = |corr(A<sub>i</sub><sup>*</sup>, "
    "B<sub>i</sub><sup>*</sup>)| using the random pairs P*.",
    numbered_step))
story.append(Paragraph(
    "After B iterations we have, per respondent, a vector "
    "(c<sub>i,1</sub><sup>*</sup>, …, c<sub>i,B</sub><sup>*</sup>) "
    "approximating the distribution of coherence under the null "
    "hypothesis that the pair structure carries no information for this "
    "respondent.",
    body))
story.append(Paragraph(
    "<b>Why this is the right null:</b> the random pairs use the same "
    "respondent's own responses, so their idiosyncratic mean and variance "
    "are absorbed. The only thing that differs between c<sub>i</sub><sup>obs</sup> "
    "and the c<sub>i,b</sub><sup>*</sup> values is whether the pairs were "
    "informative (top-k%) or random.",
    body))

# =====================================================================
# 8. Step 6 — Z-score
# =====================================================================
story.append(Paragraph("8. Step 6 — Convert to a z-score", h1))
story.append(Paragraph(
    "We compute the per-respondent mean and standard deviation of the "
    "permutation distribution and standardize the observed coherence:",
    body))
story.append(Paragraph(
    "μ<sub>i</sub> = (1/B) Σ<sub>b</sub> c<sub>i,b</sub><sup>*</sup>",
    formula_style))
story.append(Paragraph(
    "σ<sub>i</sub> = SD(c<sub>i,1</sub><sup>*</sup>, …, c<sub>i,B</sub><sup>*</sup>)",
    formula_style))
story.append(Paragraph(
    "z<sub>i</sub> = (c<sub>i</sub><sup>obs</sup> − μ<sub>i</sub>) / σ<sub>i</sub>",
    formula_style))
story.append(Paragraph(
    "When σ<sub>i</sub> = 0 (a degenerate respondent with identical "
    "coherence across all permutations) we set z<sub>i</sub> = 0.",
    body))
story.append(Paragraph(
    "<b>Interpretation:</b> z<sub>i</sub> measures, in standard deviations, "
    "how much more (or less) coherent the respondent looks on the genuine "
    "high-r pairs than on randomly chosen pairs from the same questionnaire.",
    body))
story.append(Paragraph(
    "&nbsp;&nbsp;&nbsp;• Attentive respondent: c<sub>obs</sub> &gt;&gt; μ → z &gt;&gt; 0<br/>"
    "&nbsp;&nbsp;&nbsp;• Careless respondent (random/longstring): c<sub>obs</sub> ≈ μ → z ≈ 0<br/>"
    "&nbsp;&nbsp;&nbsp;• Pure straight-liner: c<sub>obs</sub> = μ = 0 → z is set to 0",
    body_no_indent))
story.append(Paragraph(
    "A simple decision rule flags any respondent with z below a chosen "
    "threshold (default 1.5).",
    body))

# =====================================================================
# 9. Step 7 — Variance penalty
# =====================================================================
story.append(Paragraph("9. Step 7 — Optional variance penalty", h1))
story.append(Paragraph(
    "The basic z-score returns 0 for pure straight-liners (zero-variance "
    "responses), which puts them at the baseline rather than confidently "
    "below it. The optional variance penalty subtracts a function of the "
    "respondent's within-person standard deviation s<sub>i</sub>:",
    body))
story.append(Paragraph(
    "z<sub>i</sub><sup>adj</sup> = z<sub>i</sub> − α · exp(−s<sub>i</sub> / β)",
    formula_style))
story.append(Paragraph(
    "with default α = 3, β = 0.5. A respondent with s<sub>i</sub> = 0 (pure "
    "straight-liner) loses the full α from their z; a respondent with "
    "normal response variance (s<sub>i</sub> ≥ 1.5) loses essentially "
    "nothing. Tuned on simulated data, this penalty adds +0.05 sensitivity "
    "on pure-straight and acquiescent patterns at zero cost on other "
    "patterns.",
    body))

# =====================================================================
# 10. Step 8 — Iterative EFA
# =====================================================================
story.append(Paragraph("10. Step 8 — Iterative EFA refinement", h1))
story.append(Paragraph(
    "Heavy careless contamination (above 30-40%) corrupts the sample "
    "correlation matrix used in Step 1. The high-|r| pairs identified "
    "in Step 2 may include spurious correlations driven by the careless "
    "respondents themselves. The iterative-EFA refinement recomputes the "
    "correlation structure on a cleaner subset.",
    body))
story.append(Paragraph(
    "1. Run the basic z-score procedure (Steps 1-6) on all N respondents → preliminary z<sub>i</sub><sup>(1)</sup>.<br/>"
    "2. Identify the suspected careless: respondents with z<sub>i</sub><sup>(1)</sup> ≤ z<sub>cutoff</sub> (default 1.0 → ~20% trimmed).<br/>"
    "3. On the trimmed (clean) subset, run an exploratory factor analysis with parallel-analysis-determined number of factors.<br/>"
    "4. Use the clean factor loadings to identify within-factor item pairs.<br/>"
    "5. Score ALL N respondents using these clean within-factor pairs → refined z<sub>i</sub><sup>(2)</sup>.",
    numbered_step))
story.append(Paragraph(
    "The refined z-score is more discriminative because the pair selection "
    "is no longer contaminated by careless respondents themselves. In our "
    "simulations this single refinement raises ensemble MCC from 0.55 to "
    "0.62 at GT&gt;50%, with the largest gains on hard-to-detect patterns "
    "(random and fatigue).",
    body))

# =====================================================================
# 11. Step 9 — Auto-z calibration
# =====================================================================
story.append(Paragraph("11. Step 9 — Auto-z calibration of the threshold", h1))
story.append(Paragraph(
    "The optimal flagging threshold for the z-score depends on the "
    "questionnaire dimensions. We pre-computed a 2D LOESS calibration "
    "table from a 60-cell simulation grid (10 nFactors × 6 items-per-factor "
    "× 30 replications). Given a new dataset, the auto-z procedure:",
    body))
story.append(Paragraph(
    "1. Computes total_items = J.<br/>"
    "2. Optionally runs a parallel analysis to estimate nFactors (used "
    "as secondary input to the LOESS).<br/>"
    "3. Predicts the optimal z-threshold via the calibration LOESS, with "
    "fallback to the linear approximation z = 0.505 + 0.0042 · J if the "
    "LOESS fails.<br/>"
    "4. Clamps the result to [0.3, 3.5] to prevent extreme values.",
    numbered_step))
story.append(Paragraph(
    "Auto-z gives roughly the same overall MCC as fixed z = 1.5 but adapts "
    "to the questionnaire scale: on a 48-item form auto-z chooses ~0.7; on "
    "a 300-item battery it chooses ~1.8. This makes the method usable out "
    "of the box without manual threshold tuning.",
    body))

# =====================================================================
# 12. Step 10 — Multi-detector ensemble
# =====================================================================
story.append(Paragraph("12. Step 10 — Multi-detector ensemble", h1))
story.append(Paragraph(
    "The z-score above captures structurally inconsistent responses. To "
    "also detect consistent forms of carelessness (acquiescence, "
    "straight-lining), we combine z<sub>i</sub> with five auxiliary "
    "detectors that each capture an orthogonal signal.",
    body))

aux_tbl = [
    ["Detector", "What it captures", "Sign convention"],
    ["IRV (intra-individual response variability)",
     "Within-person SD of responses",
     "Low = careless"],
    ["LongString",
     "Maximum run of identical consecutive responses",
     "High = careless"],
    ["Mahalanobis D&sup2;",
     "Distance from the sample centroid in item space",
     "High = careless"],
    ["Person-Total correlation",
     "Correlation between respondent's profile and the vector of item means",
     "Low = careless"],
    ["z<sub>RR</sub>",
     "Permutation-based coherence (Steps 1-9)",
     "Low = careless"],
    ["z<sub>RR_iter_EFA</sub>",
     "Iterative-EFA refinement (Step 8)",
     "Low = careless"],
]
story.append(mktable(aux_tbl, [4*cm, 7.5*cm, 4*cm], left_cols={0, 1, 2}))
story.append(Paragraph(
    "Table. The six detectors used in the ensemble. Each captures a "
    "different aspect of careless responding.",
    caption))

story.append(Paragraph(
    "Before combination, we standardize each detector to have mean 0 and "
    "SD 1 across the sample, and we negate those whose natural sign points "
    "to \"low = careless\" so that high values uniformly mean \"more "
    "likely careless\".",
    body))

# =====================================================================
# 13. Step 11 — RF combiner
# =====================================================================
story.append(Paragraph("13. Step 11 — Random Forest combiner", h1))
story.append(Paragraph(
    "We use a Random Forest classifier (200 trees, mtry = √p) trained to "
    "predict the binary careless label from the six standardized detector "
    "scores. The classifier is trained via 5-fold cross-validation on a "
    "labeled simulation set (clean vs corruption&gt;80%), and emits a "
    "predicted probability p<sub>i</sub> ∈ [0, 1] for each respondent.",
    body))
story.append(Paragraph(
    "<b>Why Random Forest:</b> direct comparison against logistic "
    "regression and elastic-net found that RF yields MCC = 0.78 vs "
    "0.71-0.73 for the linear classifiers (a +0.05 MCC improvement). The "
    "advantage comes from RF's ability to model interactions and "
    "non-linearities in the detector space — for example, low IRV and "
    "high D² interact differently than either signal alone.",
    body))
story.append(Paragraph(
    "<b>Feature importance</b> in the trained model shows three primary "
    "contributors (IRV, D², z<sub>RR_iter</sub>) and three secondary "
    "(z<sub>RR</sub>, LongString, Person-Total). The standard z<sub>RR</sub> "
    "becomes redundant once z<sub>RR_iter</sub> is present, but is kept "
    "for diagnostic transparency.",
    body))

# =====================================================================
# 14. Step 12 — Operating point
# =====================================================================
story.append(Paragraph("14. Step 12 — Operating-point calibration (FPR = 5%)",
                        h1))
story.append(Paragraph(
    "The probability output of the Random Forest is converted to a binary "
    "flag by thresholding. Using p<sub>i</sub> &gt; 0.5 is suboptimal "
    "with imbalanced classes — the resulting threshold flags too few or "
    "too many respondents depending on the careless rate.",
    body))
story.append(Paragraph(
    "Instead we calibrate the threshold to a target false-positive rate of "
    "5% on known-clean respondents. Concretely:",
    body))
story.append(Paragraph(
    "1. Compute the predicted probability p<sub>i</sub> for every respondent.<br/>"
    "2. Take the subset of respondents known to be clean (or assumed clean — see note).<br/>"
    "3. Set the flagging threshold τ = 95th percentile of p<sub>i</sub> within that clean subset.<br/>"
    "4. Flag any respondent with p<sub>i</sub> &gt; τ as careless.",
    numbered_step))
story.append(Paragraph(
    "<b>Note on \"clean\" calibration in practice:</b> when the dataset has "
    "no labeled clean subset, a reasonable proxy is the bottom-X% of "
    "respondents on z<sub>RR_iter</sub> (i.e., the most-coherent respondents "
    "by the iterative score). This approximates a clean reference for "
    "threshold setting without external labels.",
    body))
story.append(Paragraph(
    "The 5% target is a conventional choice that balances Type I and Type "
    "II errors for most research applications. Higher targets (10%, 20%) "
    "yield higher sensitivity at the cost of more false positives.",
    body))

# =====================================================================
# 15. Worked example
# =====================================================================
story.append(Paragraph("15. Worked example — small numerical illustration",
                        h1))
story.append(Paragraph(
    "To illustrate the procedure on a tiny scale, consider a J = 4 item "
    "questionnaire administered to N = 4 respondents on a 1-5 Likert. The "
    "data matrix is:",
    body))
example_tbl = [
    ["Resp", "item 1", "item 2", "item 3", "item 4"],
    ["Alice (attentive)",   "5", "5", "1", "1"],
    ["Bob (attentive)",     "4", "4", "2", "2"],
    ["Carol (random)",      "3", "1", "5", "2"],
    ["Dan (straight-liner)", "3", "3", "3", "3"],
]
story.append(mktable(example_tbl,
                      [4*cm, 2.5*cm, 2.5*cm, 2.5*cm, 2.5*cm], left_cols={0}))

story.append(Paragraph(
    "<b>Step 1 (correlation matrix).</b> Computing R yields high positive "
    "correlation between items 1 and 2 (both go high together) and items "
    "3 and 4 (both go low when 1, 2 are high). The inter-block "
    "correlations are strongly negative.",
    body))
story.append(Paragraph(
    "<b>Step 2 (coupled pairs).</b> With corProp = 0.5 (chosen for "
    "illustration; the default 0.03 needs many more items), the top pairs "
    "are {(1, 2), (3, 4)} plus the negatively-correlated cross-block pairs.",
    body))
story.append(Paragraph(
    "<b>Step 3 (sign alignment).</b> The (1, 3) pair has negative R, so we "
    "flip item 3's values: 5 → 1, 1 → 5, etc.",
    body))
story.append(Paragraph(
    "<b>Steps 4-6 (per-respondent z).</b> Alice's row-correlation across "
    "the coupled pairs is high, while her permutation baseline is moderate; "
    "she ends up at z ≈ +2. Carol's observed and permutation correlations "
    "are both moderate-low (her responses don't track the structure); she "
    "lands at z ≈ −0.5. Dan has all-3 responses, so observed and random "
    "are both 0; z = 0 by convention.",
    body))
story.append(Paragraph(
    "<b>Step 7 (variance penalty).</b> Dan has within-person SD = 0, so he "
    "loses the full α = 3 from his z, becoming z = −3. Now confidently "
    "flagged.",
    body))
story.append(Paragraph(
    "<b>Steps 10-12 (ensemble + threshold).</b> A Random Forest trained on "
    "many simulated such datasets combines the z-score with IRV, D², and "
    "the others to produce a probability of carelessness. With FPR=5% "
    "calibration, only Carol and Dan are above the threshold.",
    body))

# =====================================================================
# 16. Pseudo-code summary
# =====================================================================
story.append(Paragraph("16. Pseudo-code summary", h1))
story.append(Paragraph(
    "The full pipeline in compact pseudo-code:",
    body))

pseudo = [
    "function REREREE_DETECT(X, corProp=0.03, B=100, FPR_target=0.05):",
    "    # Step 1",
    "    R = correlation_matrix(X)",
    "    # Step 2",
    "    P = top_corProp_pairs(|R|, corProp)",
    "    # Step 3",
    "    P_aligned = sign_align(P, R)",
    "    # Step 4",
    "    for i in 1..N:",
    "        c_obs[i] = |row_corr(X[i, P_aligned])|",
    "    # Step 5",
    "    for b in 1..B:",
    "        P_random = sample_random_pairs(J, k=|P|, no_self_pairs)",
    "        for i in 1..N:",
    "            c_star[i, b] = |row_corr(X[i, sign_align(P_random, R)])|",
    "    # Step 6",
    "    mu = rowMeans(c_star);  sigma = rowSDs(c_star)",
    "    z_RR = (c_obs - mu) / sigma   # zero where sigma=0",
    "    # Step 7  (optional)",
    "    if variance_penalty:",
    "        s = within_person_SD(X)",
    "        z_RR -= alpha * exp(-s / beta)",
    "    # Step 8  (optional, for high contamination)",
    "    if iterative:",
    "        suspect = z_RR <= z_cutoff",
    "        EFA_clean = exploratory_FA(X[!suspect])",
    "        z_RR_iter = repeat_steps_4_6(X, pairs_from(EFA_clean))",
    "    # Steps 10-11",
    "    features = [z_RR, z_RR_iter, IRV(X), longstring(X), D2(X), person_total(X)]",
    "    probs = random_forest_predict(features)   # trained on simulation",
    "    # Step 12",
    "    tau = quantile(probs[clean_subset], 0.95)",
    "    flagged = probs > tau",
    "    return flagged, probs, z_RR",
]
for line in pseudo:
    story.append(code_block(line))

# Closing
story.append(Spacer(1, 0.3*cm))
story.append(Paragraph(
    "This concludes the step-by-step description of the procedure. The "
    "individual steps can be combined modularly: the basic z-score from "
    "Steps 1-6 alone yields MCC ~ 0.40 in challenging conditions; adding "
    "Steps 7-9 gives ~ 0.55; the full ensemble of Steps 10-12 reaches MCC "
    "= 0.78-0.84 depending on the careless rate. Performance scales "
    "smoothly with questionnaire length, exceeding MCC = 0.95 on 300-item "
    "batteries.",
    body))

story.append(Spacer(1, 0.5*cm))
story.append(Paragraph(
    "<i>Code:</i> github.com/vittorio-g/ReReReRe &middot; "
    "<i>Compiled:</i> 2026-04-25",
    note))

doc = SimpleDocTemplate(OUT, pagesize=A4,
                         leftMargin=2.4*cm, rightMargin=2.4*cm,
                         topMargin=2.2*cm, bottomMargin=2.2*cm,
                         title="ReReReRe Statistical Methodology",
                         author="Vittorio Guerri")
doc.build(story)
print(f"Saved: {OUT}")
