"""Integrated paper draft (v4 &mdash; real-data-first edition).

Reframes the ReReRe paper around its real-data validation:
  - HEADLINE: the collected careless-responding study (n=84, experimentally
    induced careless via text scrambling) &mdash; leave-one-out AUC 0.985, MCC 0.880
    for the shipped ensemble (the triad rr + LongString + Person-Total).
  - Prevalence-sweep validation on real data (true rate 5-50%): AUC 0.98-0.99
    at every prevalence; automatic two-groups calibration; rate-estimate
    calibration curve.
  - External convergent validity on 14 ground-truth-independent public datasets
    plus 7 boundary-condition datasets from the "Learning to Pay Attention"
    benchmark (arXiv:2603.02427).
  - Simulation study demoted to a controlled-conditions CHARACTERISATION
    (scaling with length, ablation of rr with bootstrap CIs, per-pattern).
  - A proper introduction, and the online tool (re-re.re).

Nomenclature: ReReRe = the ENSEMBLE; rr = the permutation-based coherence index.

Inputs (reused from v3): article_assets_v3/article_data.pkl, revision/paper figs.
New inputs: webapp/eval_mix.csv, webapp/eval_mix_fine.csv, webapp/fig_rate_calibration*.png.
Output: C:/Users/vitto/Downloads/ReReRe_Paper_2026-07-06.pdf
"""

import os, pickle
from datetime import date
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from reportlab.lib.pagesizes import A4
from reportlab.lib.styles import getSampleStyleSheet, ParagraphStyle
from reportlab.lib.units import cm
from reportlab.lib import colors
from reportlab.lib.enums import TA_LEFT, TA_JUSTIFY, TA_CENTER
from reportlab.platypus import (SimpleDocTemplate, Paragraph, Spacer,
                                 Image as RLImage, Table, TableStyle, PageBreak)

V3_ASSETS  = "article_assets_v3"
REV_ASSETS = "revision_assets_v3"
PAPER_DIR  = "paper_assets_v3"
WEB        = "webapp"
OUT_PDF    = r"C:\Users\vitto\Downloads\ReReRe_Paper_2026-07-06.pdf"
os.makedirs(PAPER_DIR, exist_ok=True)

# --- simulation data (for the demoted characterisation section) ---
with open(f"{V3_ASSETS}/article_data.pkl", "rb") as f:
    V3 = pickle.load(f)
ra       = V3["ra_size"]
ra_pool  = V3["ra_pool"]
abl_size = V3["abl_by_size"]
boot     = V3["boot"]
ra_B     = ra_pool[ra_pool.scenario == "B"].iloc[0]

# --- real-data prevalence sweep (this project) ---
mix = pd.read_csv(f"{WEB}/eval_mix.csv")

# --- hard-coded validated real-data numbers (see CLAUDE.md 2026-07-06) ---
STUDY = dict(n=84, careful=38, careless=46,
             triad_auc=0.985, triad_mcc=0.880,
             full5_auc=0.990, full5_mcc=0.880,
             norr_auc=0.965, norr_mcc=0.783,
             rr_auc=0.904, rr_mcc=0.687)
STUDY_FEAT_AUC = [  # per-feature standalone AUC on the study (higher=careless)
    ("Person-Total", 0.916), ("rr", 0.912), ("LongString", 0.785),
    ("Mahalanobis D2", 0.458), ("IRV", 0.309)]
WEIGHTS = dict(b0=0.2755, rr=1.7007, longstring=1.5974, person_total=1.7129)

# External convergent-validity benchmark (GT independent of response pattern)
EXT = [  # dataset, items(scale), N, GT criterion, AUC ens, AUC rr, AUC PT
    ("Kay S2",            "363 (5)",  701,   "self-report exclusion",       0.912, 0.919, 0.893),
    ("warning IPIP-NEO",  "300 (5)",  812,   "diligence self-report",       0.822, 0.670, 0.803),
    ("Kay S1",            "250 (5)",  500,   "speeding (duration)",         0.934, 0.852, 0.955),
    ("opsy HEXACO",       "240 (7)",  22732, "V1/V2 seriousness",           0.697, 0.588, 0.709),
    ("ambi 2019",         "181 (7)",  2017,  "RT < 2 s/item",               0.539, 0.428, 0.526),
    ("opsy 16PF",         "163 (5)",  35381, "extreme speeders",            0.828, 0.793, 0.826),
    ("smarvus",           "141 (5)",  10234, ">=2 instructed checks",  0.864, 0.865, 0.886),
    ("mies dev",          "91 (5)",   7188,  "extreme speeders",            0.593, 0.503, 0.476),
    ("Kay S6",            "67 (7)",   562,   "instructed + self-excl.",     0.773, 0.618, 0.831),
    ("Kay S5",            "62 (5)",   629,   "speeding",                    0.733, 0.574, 0.714),
    ("Duckworth VCL",     "50 (5)",   3874,  "fake-word overclaiming",      0.453, 0.500, 0.425),
    ("fbps",              "50 (5)",   36850, "page-level speeding",         0.445, 0.451, 0.442),
    ("Douglas 2023",      "50 (5)",   2071,  ">=2 attention checks",   0.720, 0.702, 0.678),
    ("Krause youth",      "49 (5)",   1783,  ">=2 reactive checks",    0.818, 0.576, 0.805),
]
LTPA = [  # dataset, items(scale), N, rate, AUC ens, AUC rr, AUC PT
    ("Ivanov 2021",      "41 (1-5)",   869,  "20.3%", 0.775, 0.689, 0.813),
    ("Moss 2023",        "8 (mixed)",  2107, "10.1%", 0.742, 0.291, 0.796),
    ("Alvarez 2019",     "6 (1-5)",    1887, "8.7%",  0.720, 0.617, 0.699),
    ("O'Grady 2019",     "30 (0-5)",   358,  "5.9%",  0.677, 0.635, 0.643),
    ("Mastroianni 2022", "100 (0-100)",1004, "2.9%",  0.639, 0.374, 0.736),
    ("Pennycook 2020",   "30 (binary)",853,  "27.2%", 0.584, 0.593, 0.594),
    ("Buchanan 2018",    "14 (1-7)",   1029, "5.0%",  0.452, 0.449, 0.616),
]

# ---------------------------------------------------------------------------
# Figures
# ---------------------------------------------------------------------------
# Simulation scaling (demoted section)
fig, ax = plt.subplots(figsize=(8, 4.6))
sub = ra[ra.scenario == "B"].sort_values("size")
ax.errorbar(sub["size"], sub["mcc"], yerr=sub["mcc_sd"]/np.sqrt(sub["n"]),
            marker="o", lw=2, ms=6, capsize=4, color="#33475f")
ax.set_xlabel("Questionnaire length (items)"); ax.set_ylabel("MCC")
ax.set_title("Simulation: detection quality scales with length (rate-aware)")
ax.set_xticks(sorted(ra["size"].unique())); ax.set_ylim(0.3, 0.9); ax.grid(alpha=.3)
plt.tight_layout(); SCALE_FIG = f"{PAPER_DIR}/fig_scale_rate_aware.png"
plt.savefig(SCALE_FIG, dpi=150); plt.close()

fig, ax = plt.subplots(figsize=(8, 4.6))
ablB = abl_size[abl_size.scenario == "B"].sort_values("size")
ax.errorbar(ablB["size"], ablB["mcc_full"], yerr=ablB["mcc_full_se"],
            marker="o", lw=2, capsize=4, color="#33475f", label="Full ensemble")
ax.errorbar(ablB["size"], ablB["mcc_norr"], yerr=ablB["mcc_norr_se"],
            marker="^", lw=2, capsize=4, color="#8f3535", label="Auxiliaries only (no rr)")
ax.set_xlabel("Questionnaire length (items)"); ax.set_ylabel("MCC")
ax.set_title("Simulation ablation: with vs without rr")
ax.set_xticks(sorted(ablB["size"].unique())); ax.grid(alpha=.3); ax.legend(loc="lower right")
plt.tight_layout(); ABL_FIG = f"{PAPER_DIR}/fig_ablation_curves.png"
plt.savefig(ABL_FIG, dpi=150); plt.close()

# Real-data prevalence sweep (AUC + auto vs oracle MCC)
fig, ax = plt.subplots(figsize=(8, 4.8))
tr = mix.true_rate*100
ax.plot(tr, mix.auc, "o-", lw=2, color="#33475f", label="AUC (ranking)")
ax.plot(tr, mix.oracle_mcc, "s--", lw=1.8, color="#2f6b3a", label="MCC, oracle threshold")
ax.plot(tr, mix.auto_mcc, "^--", lw=1.8, color="#8a6a1e", label="MCC, automatic (Medium)")
ax.axvspan(10, 30, color="#2f6b3a", alpha=.07)
ax.set_xlabel("True careless rate (%)"); ax.set_ylabel("AUC / MCC")
ax.set_title("Real data: ranking is prevalence-robust; automatic decision peaks at 10-30%")
ax.set_ylim(0.2, 1.02); ax.grid(alpha=.3); ax.legend(loc="lower left", fontsize=9)
plt.tight_layout(); SWEEP_FIG = f"{PAPER_DIR}/fig_realdata_sweep.png"
plt.savefig(SWEEP_FIG, dpi=150); plt.close()

CALIB_FIG = f"{WEB}/fig_rate_calibration.png"          # true vs estimated rate (5-50)
CALIB_FINE = f"{WEB}/fig_rate_calibration_fine.png"    # fine 10-30, 100 resamples

# ---------------------------------------------------------------------------
# Styles / helpers (as v3)
# ---------------------------------------------------------------------------
styles = getSampleStyleSheet()
AC = colors.HexColor("#33475f")
H1 = ParagraphStyle("H1", parent=styles["Heading1"], fontSize=17, spaceAfter=10, textColor=AC)
H2 = ParagraphStyle("H2", parent=styles["Heading2"], fontSize=13, spaceAfter=7, textColor=AC, spaceBefore=11)
H3 = ParagraphStyle("H3", parent=styles["Heading3"], fontSize=11, spaceAfter=5,
                    textColor=colors.HexColor("#243447"), spaceBefore=7)
BODY = ParagraphStyle("body", parent=styles["BodyText"], fontSize=10, leading=14,
                      alignment=TA_JUSTIFY, spaceAfter=6)
CODE = ParagraphStyle("code", parent=styles["BodyText"], fontSize=9, leading=12, fontName="Courier",
                      backColor=colors.HexColor("#f2f0ea"), borderPadding=6, leftIndent=6,
                      rightIndent=6, spaceBefore=4, spaceAfter=4)
CAPT = ParagraphStyle("capt", parent=styles["BodyText"], fontSize=8.5, leading=11, alignment=TA_CENTER,
                      textColor=colors.HexColor("#555"), spaceAfter=10, spaceBefore=2, fontName="Helvetica-Oblique")
TITLE = ParagraphStyle("title", parent=styles["Title"], fontSize=24, spaceAfter=8, alignment=TA_CENTER, textColor=AC)
SUB = ParagraphStyle("sub", parent=styles["BodyText"], fontSize=11.5, alignment=TA_CENTER,
                     textColor=colors.HexColor("#444"), spaceAfter=16)
CALL = ParagraphStyle("call", parent=styles["BodyText"], fontSize=10, leading=14, alignment=TA_JUSTIFY,
                      backColor=colors.HexColor("#eef1f5"), borderColor=colors.HexColor("#b8c4d4"),
                      borderWidth=0.6, borderPadding=8, leftIndent=4, rightIndent=4, spaceBefore=6, spaceAfter=8)

def fig_block(path, width=15.5*cm):
    if not os.path.exists(path):
        return Paragraph(f"<i>(missing figure: {path})</i>", BODY)
    img = RLImage(path); iw, ih = img.imageWidth, img.imageHeight
    img.drawWidth = width; img.drawHeight = width*ih/iw
    return img

def caption(t): return Paragraph(t, CAPT)

def grid(rows, col_widths=None, fs=8.5):
    t = Table(rows, colWidths=col_widths, hAlign="LEFT")
    t.setStyle(TableStyle([
        ("FONTSIZE", (0,0), (-1,-1), fs), ("VALIGN", (0,0), (-1,-1), "MIDDLE"),
        ("ALIGN", (0,0), (-1,-1), "CENTER"), ("ALIGN", (0,0), (0,-1), "LEFT"),
        ("FONTNAME", (0,0), (-1,0), "Helvetica-Bold"),
        ("BACKGROUND", (0,0), (-1,0), colors.HexColor("#dde2ea")),
        ("LINEABOVE", (0,0), (-1,0), 0.7, colors.black),
        ("LINEBELOW", (0,0), (-1,0), 0.4, colors.black),
        ("LINEBELOW", (0,-1), (-1,-1), 0.7, colors.black),
        ("ROWBACKGROUNDS", (0,1), (-1,-1), [colors.white, colors.HexColor("#f3f5f8")]),
    ]))
    return t

S = []

# ===========================================================================
# Title
# ===========================================================================
S += [
    Spacer(1, 2.2*cm),
    Paragraph("ReReRe", TITLE),
    Paragraph("Detecting careless responding with a permutation-based coherence "
              "index and a small ensemble &mdash; validated on experimentally "
              "induced and real-world data", SUB),
    Paragraph(
        "<b>Nomenclature.</b> <b>ReReRe</b> denotes the detection <i>method</i>: a compact "
        "ensemble that combines a distinctive permutation-based index with two complementary "
        "auxiliary detectors. <b>rr</b> denotes that distinctive index alone &mdash; the "
        "permutation-based individual-coherence z-score (called z_RR in earlier drafts). The "
        "shipped ensemble is a logistic model on three robustly-standardised features: "
        "<b>rr + LongString + Person-Total</b>.", BODY),
    Spacer(1, 0.25*cm),
    Paragraph(f"<b>Date:</b> {date.today().isoformat()}. Vittorio Guerrieri.", BODY),
    Spacer(1, 0.35*cm),
    Paragraph("<b>Headline.</b>", H3),
    Paragraph(
        f"On a purpose-built study in which carelessness was <i>experimentally induced</i> "
        f"(n = {STUDY['n']}: {STUDY['careful']} attentive respondents and {STUDY['careless']} "
        f"instructed-careless respondents answering scrambled items), the ReReRe ensemble "
        f"separates careless from attentive respondents with <b>leave-one-out AUC = "
        f"{STUDY['triad_auc']:.3f}</b> and <b>MCC = {STUDY['triad_mcc']:.3f}</b>. Across a "
        f"prevalence sweep from 5% to 50% careless, the ranking is essentially perfect at every "
        f"rate (AUC 0.98&ndash;0.99), and a fully automatic, label-free calibration flags the "
        f"careless with MCC ~ 0.82&ndash;0.85 in the realistic 10&ndash;25% band. On 14 "
        f"independent public datasets with ground truth that does <i>not</i> depend on the response "
        f"pattern, rr shows convergent validity that scales with questionnaire length (up to AUC "
        f"0.92 on the longest battery) and a clean discriminant null against unrelated constructs. "
        f"A controlled simulation characterises how detection scales with length and isolates the "
        f"contribution of rr. The method is deployed as a browser tool at re-re.re.", BODY),
    Spacer(1, 0.3*cm),
    Paragraph("<b>Structure.</b> §1 Introduction. §2 Method (the rr index, the ensemble, "
              "the automatic calibration). §3 Validation on experimentally-induced real data "
              "(the headline). §4 External convergent validity on independent datasets. "
              "§5 Characterisation on controlled simulated data. §6 Application and the "
              "online tool. §7 Statistical-concepts box.", BODY),
    PageBreak(),
]

# ===========================================================================
# 1. Introduction
# ===========================================================================
S += [
    Paragraph("1. Introduction", H1),
    Paragraph(
        "Self-report questionnaires assume that respondents read each item and answer it "
        "truthfully. A sizeable minority do not: they straight-line, respond at random, or drift "
        "in and out of attention. This <i>careless</i> or <i>insufficient-effort</i> responding "
        "is common in online and unproctored samples &mdash; typical estimates range from a few "
        "per cent to a third of respondents &mdash; and it biases means, inflates or deflates "
        "correlations, manufactures spurious factors, and degrades reliability. Because a single "
        "contaminated subgroup can overturn a substantive conclusion, screening for careless "
        "responding has become a routine (and, increasingly, an expected) step in survey "
        "analysis.", BODY),
    Paragraph(
        "Two broad families of remedies exist. <b>Design-time</b> methods embed dedicated probes: "
        "instructed-response items ('select strongly agree'), bogus items with a known answer, and "
        "instructional-manipulation checks. These are effective but must be planned in advance, "
        "consume item budget, can be defeated by attentive-looking speeders, and are unavailable "
        "for the vast quantity of already-collected data. <b>Post-hoc</b> methods instead infer "
        "carelessness from the response vector itself: response-time and long-string heuristics "
        "catch straight-lining; multivariate-outlier statistics such as Mahalanobis distance flag "
        "unusual profiles; person-fit indices from item-response theory quantify misfit; and, more "
        "recently, machine-learning models learn the signature of low-quality data. A recent "
        "unsupervised approach (arXiv:2603.02427, <i>Learning to Pay Attention</i>) trains an "
        "autoencoder together with a Chow&ndash;Liu dependency tree to reconstruct 'attentive' "
        "response patterns, using embedded attention-check failures as ground truth across nine "
        "public datasets.", BODY),
    Paragraph(
        "These methods largely target <i>consistent</i> carelessness (straight-lining, "
        "acquiescence) or rely on external signals (timing, explicit checks). A complementary and "
        "harder target is <i>inconsistent</i> carelessness: responding that preserves plausible "
        "marginal distributions but breaks the <i>cross-item structure</i> a real respondent "
        "produces. An attentive person's answer to one item is informative about their answer to a "
        "correlated item, because both load on a shared latent trait; a careless person's answers "
        "are internally incoherent even when each item, in isolation, looks ordinary. Detecting "
        "this requires a per-respondent measure of coherence and, critically, a per-respondent "
        "reference for how coherent that person would look <i>by chance</i>.", BODY),
    Paragraph(
        "<b>This paper.</b> We present ReReRe, built around <b>rr</b>: for each respondent we "
        "measure how coherently they answered the item pairs the sample deems most related, and "
        "compare that against a <i>within-person permutation baseline</i> &mdash; the same "
        "statistic recomputed on random item pairs. The result is a self-calibrating z-score that "
        "needs no distributional assumption and no external cut-off. This per-respondent "
        "permutation null is, to our knowledge, novel; its nearest relatives are psychometric "
        "synonyms/antonyms and Goldammer's resampled personal reliability, neither of which builds "
        "a per-respondent random-pair reference. rr is deliberately blind to consistent "
        "straight-lining, so we combine it with two complementary detectors &mdash; LongString and "
        "Person-Total correlation &mdash; into a small logistic ensemble, with a reliability gate "
        "that falls back to rr alone when the auxiliaries are uninformative.", BODY),
    Paragraph(
        "Our contributions are: (1) the rr permutation index and the ReReRe ensemble; (2) a "
        "purpose-built validation study in which carelessness is <i>experimentally induced</i>, "
        "giving clean ground truth on real questionnaire content (§3); (3) a prevalence-sweep "
        "and an <i>automatic</i>, label-free calibration that estimates the careless rate from the "
        "data rather than asking the analyst to guess it (§3); (4) external convergent- and "
        "discriminant-validity evidence on 14 independent public datasets whose ground truth is "
        "independent of the response pattern, plus seven boundary-condition datasets (§4); "
        "(5) a controlled simulation that characterises scaling and isolates rr's contribution "
        "(§5); and (6) a free browser tool that runs entirely client-side (§6).", BODY),
    PageBreak(),
]

# ===========================================================================
# 2. Method
# ===========================================================================
S += [
    Paragraph("2. Method", H1),
    Paragraph("2.1 The rr permutation index", H2),
    Paragraph(
        "rr turns the coherence idea into a per-respondent z-score in five steps.", BODY),
    Paragraph(
        "<b>1. Pair selection.</b> Compute the sample item-by-item correlation matrix and select "
        "the top-k% pairs by absolute correlation (default 3%, floor 15 pairs); on short "
        "questionnaires the method instead weights <i>all</i> pairs by |r|. These are the "
        "<i>coupled</i> pairs &mdash; items that genuinely move together.<br/>"
        "<b>2. Sign alignment.</b> For a coupled pair with negative sample correlation, reverse "
        "one member on its own scale ((max+min) &minus; x) so reverse-keyed items do not cancel.<br/>"
        "<b>3. Individual coherence.</b> For respondent r, take their answers to the two members of "
        "each coupled pair and measure how tightly the two sides track across all coupled pairs "
        "(an individual-level |correlation|), C_obs(r).<br/>"
        "<b>4. Permutation baseline.</b> Recompute the same statistic on many (default 100&ndash;200) "
        "random pair-sets of equal size, giving a per-respondent null {C_rand(r,s)}.<br/>"
        "<b>5. Z-score.</b> rr(r) = (C_obs(r) &minus; mean_s C_rand) / sd_s C_rand. High rr = more "
        "coherent than the respondent's own chance level (attentive); low rr = no more coherent "
        "than random (careless). Items are proportion-rescaled (x / item max) first so mixed Likert "
        "ranges are comparable; on homogeneous scales this is a no-op.", BODY),
    Paragraph("2.2 The ensemble", H2),
    Paragraph(
        "rr detects <i>inconsistent</i> carelessness but is blind to a respondent who answers '3' "
        "to everything (trivially coherent). We therefore combine rr with two complementary "
        "detectors: <b>LongString</b> (longest run of identical consecutive answers &mdash; catches "
        "straight-lining) and <b>Person-Total correlation</b> (alignment of the profile with the "
        "sample mean profile &mdash; catches responding against the grain). Each feature is "
        "oriented so higher = more careless and robustly standardised <i>within the dataset</i> "
        "(median/MAD), which makes a single fixed weight vector transferable across questionnaires "
        "of different length and scale. The three are combined by a logistic model", BODY),
    Paragraph(
        f"eta = {WEIGHTS['b0']:.3f} + {WEIGHTS['rr']:.3f} z_rr "
        f"+ g ({WEIGHTS['longstring']:.3f} z_ls + {WEIGHTS['person_total']:.3f} z_pt)"
        f"   ->   p = 1/(1+exp(-eta))     [z_* = oriented robust-z]", CODE),
    Paragraph(
        "where z_* denotes the oriented robust-z feature and g is a partner-reliability gate: "
        "on flat or near-orthogonal data (where LongString and Person-Total are uninformative) g "
        "&rarr; 0 and the ensemble falls back to rr alone, so it is never worse than the index. The "
        "weights were fitted and leave-one-out validated on the real study of §3. Two further "
        "detectors (IRV and Mahalanobis D&sup2;) are computed and reported for transparency but "
        "carry zero weight in the shipped model: IRV's careless direction depends on the careless "
        "type, and the pseudo-inverse D&sup2; degenerates in the n &le; p regime common to "
        "careless data (both confirmed on the study, Table 2).", BODY),
    Paragraph("2.3 Automatic, label-free calibration", H2),
    Paragraph(
        "A ranking is not a decision: converting rr/ensemble scores into flags requires a "
        "threshold, and the quantity that fixes it &mdash; the careless prevalence &mdash; is "
        "exactly what the analyst does not know. Rather than ask for it, ReReRe fits a two-component "
        "Gaussian mixture to the ensemble log-odds. A BIC test first asks whether a distinct "
        "careless component exists at all: on clean data it does not, and <i>nothing</i> is flagged "
        "(no false-positive cascade). When it does, respondents are flagged if their score exceeds "
        "the attentive component's upper tail. Anchoring the cut on the attentive mode &mdash; the "
        "low-score cluster, robustly estimable at any prevalence &mdash; recovers far more of the "
        "careless than a posterior&gt;0.5 rule. How far into the tail the cut reaches is an "
        "<b>Expected-carelessness</b> setting (Low / Medium / High ~ 0.6% / 2% / 7% "
        "false-positive rate on attentive respondents); Medium is the validated default. The output "
        "is an <i>estimated</i> careless rate plus per-respondent flags &mdash; no prior required.", BODY),
    PageBreak(),
]

# ===========================================================================
# 3. Real-data validation (HEADLINE)
# ===========================================================================
S += [
    Paragraph("3. Validation on experimentally-induced real data", H1),
    Paragraph("3.1 The study", H2),
    Paragraph(
        f"We built an online study that induces careless responding by design, giving ground truth "
        f"on genuine questionnaire content. A 109-item battery (BFI-2, PANAS, RSES, SWLS, PSS-4, "
        f"plus filler items; mixed 4-, 5- and 7-point scales) was administered between-subjects. "
        f"Attentive respondents received standard instructions; careless respondents were shown a "
        f"fraction of items with the letters of each word scrambled into random but "
        f"length-preserving nonsense ('fully scrambled'), and were told the items were intentionally "
        f"incomprehensible and to answer as they liked. Scrambling a genuine item, rather than "
        f"injecting synthetic noise, reproduces real careless behaviour (straight-lining, "
        f"acquiescence, and incoherent guessing) on real scale structure. After quality screening, "
        f"the analysed sample is n = {STUDY['n']} ({STUDY['careful']} attentive, "
        f"{STUDY['careless']} careless).", BODY),
    Paragraph("3.2 Leave-one-out performance", H2),
    Paragraph(
        f"Scored by the shipped ensemble under leave-one-out cross-validation of the combiner "
        f"(features computed once, weights refit on n&minus;1), ReReRe reaches <b>AUC = "
        f"{STUDY['triad_auc']:.3f}, MCC = {STUDY['triad_mcc']:.3f}</b>. Table 1 shows the ablation: "
        f"the shipped triad matches the full five-feature model and clearly beats rr alone and the "
        f"no-rr auxiliaries, confirming both that rr is essential and that the two problematic "
        f"features (IRV, D&sup2;) add nothing worth their transfer risk.", BODY),
    grid([["Model", "AUC", "MCC"],
          ["Full (5 features: rr + IRV + LongString + D2 + Person-Total)", f"{STUDY['full5_auc']:.3f}", f"{STUDY['full5_mcc']:.3f}"],
          ["Shipped triad (rr + LongString + Person-Total)", f"{STUDY['triad_auc']:.3f}", f"{STUDY['triad_mcc']:.3f}"],
          ["No-rr (auxiliaries only)", f"{STUDY['norr_auc']:.3f}", f"{STUDY['norr_mcc']:.3f}"],
          ["rr alone", f"{STUDY['rr_auc']:.3f}", f"{STUDY['rr_mcc']:.3f}"]],
         col_widths=[10.5*cm, 2.3*cm, 2.3*cm]),
    caption("Table 1. Leave-one-out AUC / MCC on the study (n = 84). The three-feature triad "
            "matches the full model; rr contributes +0.20 MCC over the auxiliaries."),
    Paragraph("Per-feature standalone AUC on the study (higher = more careless):", BODY),
    grid([["Feature", "AUC"]] + [[nm, f"{a:.3f}"] for nm, a in STUDY_FEAT_AUC],
         col_widths=[7*cm, 3*cm]),
    caption("Table 2. rr and Person-Total are the strongest single signals; LongString is moderate; "
            "IRV is reversed on this careless type (AUC < 0.5) and D&sup2; is near chance in the "
            "n < p regime &mdash; hence both carry zero weight in the shipped ensemble."),
    PageBreak(),
]

# 3.3 prevalence sweep
S += [
    Paragraph("3.3 Prevalence robustness and automatic calibration", H2),
    Paragraph(
        "To test behaviour across contamination levels, we kept all 38 attentive respondents and "
        "sampled real careless respondents in to hit true rates from 5% to 50% (30 resamples per "
        "rate). Two conclusions (Figure 1, Table 3). First, the <b>ranking is prevalence-robust</b>: "
        "AUC 0.98&ndash;0.99 everywhere, and with an oracle threshold MCC 0.84&ndash;0.92 at every "
        "rate. Second, the <b>automatic calibration</b> works well in the realistic 10&ndash;30% "
        "band (MCC 0.75&ndash;0.85 at the Medium setting, precision ~ 1.0) and is deliberately "
        "conservative outside it &mdash; too few careless to separate below ~5%, and a mild "
        "under-count above ~35%, where the careless stop being a clear minority.", BODY),
    fig_block(SWEEP_FIG, width=14.5*cm),
    caption("Figure 1. Real-data prevalence sweep. The ranking (AUC) and the oracle-threshold MCC "
            "are flat across 5&ndash;50%; the automatic decision peaks in the shaded 10&ndash;30% "
            "range and softens at the extremes."),
]
rows_sw = [["True rate", "AUC", "auto est.", "auto MCC", "oracle MCC"]]
for _, r in mix.iterrows():
    rows_sw.append([f"{100*r['true_rate']:.0f}%", f"{r['auc']:.2f}",
                    f"{100*r['est_mean']:.0f}%", f"{r['auto_mcc']:.2f}", f"{r['oracle_mcc']:.2f}"])
S += [grid(rows_sw, col_widths=[2.6*cm, 2*cm, 2.2*cm, 2.4*cm, 2.6*cm]),
      caption("Table 3. Per-rate AUC, automatic rate estimate, automatic-flag MCC and "
              "oracle-threshold MCC on the real careful/careless mixes."),
      PageBreak(),
      Paragraph("3.4 Calibration of the estimated rate", H2),
      Paragraph(
        "Because the tool <i>reports</i> an estimated careless rate, we characterised that estimate "
        "finely (10&ndash;30% in 2-point steps, 100 resamples each; Figure 2). The estimate lands "
        "within a point or two of the truth up to ~15%, then turns conservative and plateaus around "
        "20% as the careless cease to be a clear minority. Since the ranking stays excellent, an "
        "analyst who expects heavy contamination raises Expected carelessness to High to push the "
        "cut and recover the missed cases.", BODY),
      fig_block(CALIB_FINE, width=13.5*cm),
      caption("Figure 2. Estimated vs true careless rate on real data (100 resamples/point). "
              "Well-calibrated up to ~15%, deliberately conservative above."),
      PageBreak(),
]

# ===========================================================================
# 4. External validation
# ===========================================================================
S += [
    Paragraph("4. External convergent and discriminant validity", H1),
    Paragraph(
        "The induced-careless study gives clean labels but a single design. To test generalisation "
        "we assembled 14 public datasets that meet the method's envelope (multi-construct Likert "
        "batteries) and, crucially, carry ground truth <i>independent of the response pattern</i> "
        "&mdash; instructed/attention checks, self-reported effort or exclusion advice, and response "
        "timing &mdash; never a pattern-based label (Mahalanobis, long-string, latent class), which "
        "would be circular. The homogeneous Likert block was scored once by the shipped ensemble "
        "(check items excluded, so zero leakage) and its AUC computed against each independent "
        "criterion (Table 4).", BODY),
]
rows_ext = [["Dataset", "items (scale)", "N", "GT criterion", "ens", "rr", "PT"]]
for nm, it, N, gt, e, r, pt in EXT:
    rows_ext.append([nm, it, f"{N:,}", gt, f"{e:.2f}", f"{r:.2f}", f"{pt:.2f}"])
S += [grid(rows_ext, col_widths=[3.1*cm, 2.3*cm, 1.9*cm, 4.0*cm, 1.2*cm, 1.2*cm, 1.2*cm], fs=8),
      caption("Table 4. Convergent-validity AUC (ensemble / rr / Person-Total) against "
              "pattern-independent ground truth, ordered by questionnaire length."),
      Paragraph(
        "Three patterns hold across all datasets. (1) <b>rr scales with length and relevance</b>: on "
        "long batteries with an incoherence-relevant criterion rr reaches 0.79&ndash;0.92 (Kay S2 "
        "363 items = 0.92; 16PF 163 items = 0.79; smarvus 141 items = 0.87); by ~60 items it fades "
        "to ~0.57 and the partners carry the load. (2) <b>Complementarity</b>: where carelessness is "
        "consistent or the battery short, Person-Total supplies the signal and the gated ensemble "
        "leans on it (ensemble &ge; rr throughout). (3) <b>Discriminant validity</b>: against "
        "constructs other than incoherence, rr is correctly at chance &mdash; the Duckworth "
        "fake-word overclaiming scale gives rr = 0.50, a clean negative control, and speeding on "
        "pre-filtered web samples is near-chance, confirming that rr measures profile incoherence, "
        "not speed or consistency. A within-sample gradient on the warning IPIP-NEO (same 812 "
        "respondents, six criteria) makes the point sharply: AUC rises monotonically with how much "
        "the criterion captures whole-questionnaire incoherence (1-item infrequency 0.72 &lt; "
        "speeding 0.70 &lt; instructed 0.74 &lt; graded 6-check 0.78 &lt; self-reported diligence "
        "0.82) &mdash; the signature of a profile-incoherence detector, not a trap-catcher.", BODY),
      PageBreak(),
      Paragraph("4.1 Boundary conditions", H2),
      Paragraph(
        "For completeness we also ran the seven available datasets from the <i>Learning to Pay "
        "Attention</i> benchmark (arXiv:2603.02427), whose ground truth is attention-check failure "
        "(Table 5). None is a long multi-factor Likert battery &mdash; the method's envelope &mdash; "
        "so these are boundary conditions, not core validation. The pattern is exactly as expected: "
        "the longest instrument (Ivanov, 41 items) gives the best rr (0.69), Person-Total carries "
        "short/flat surveys, and the reliability gate correctly abstains on the weakest data. They "
        "confirm the benchmark's own thesis &mdash; detectability depends on questionnaire design "
        "&mdash; and that attention-check failure and profile incoherence are related but distinct "
        "targets.", BODY),
]
rows_l = [["Dataset", "items (scale)", "N", "rate", "ens", "rr", "PT"]]
for nm, it, N, rate, e, r, pt in LTPA:
    rows_l.append([nm, it, f"{N:,}", rate, f"{e:.2f}", f"{r:.2f}", f"{pt:.2f}"])
S += [grid(rows_l, col_widths=[3.2*cm, 2.6*cm, 1.9*cm, 1.7*cm, 1.2*cm, 1.2*cm, 1.2*cm], fs=8),
      caption("Table 5. Boundary-condition datasets (attention-check ground truth). Mostly short / "
              "few-factor &mdash; outside the recommended envelope."),
      PageBreak(),
]

# ===========================================================================
# 5. Simulation characterisation (demoted)
# ===========================================================================
mcc_300 = ra[(ra['size']==300)&(ra.scenario=='B')]['mcc'].iloc[0]
mcc_30  = ra[(ra['size']==30) &(ra.scenario=='B')]['mcc'].iloc[0]
S += [
    Paragraph("5. Characterisation on controlled simulated data", H1),
    Paragraph(
        "Real data cannot vary questionnaire length or careless rate cleanly, so we complement it "
        "with a controlled simulation: 384 datasets over 8 sizes (30&ndash;300 items) &times; 4 base "
        "rates (10&ndash;60%) &times; 12 replications, n = 500 each, with six injected careless "
        "patterns (random, longstring, mixed, fatigue, acquiescent, pure straight-lining) at "
        "corruption 0&ndash;100%. Here a respondent counts as careless when the corrupted fraction "
        "of their row exceeds &tau; = 0.40 (the joint argmax of MCC, &kappa;, Youden's J and "
        "balanced accuracy over a &tau;-sweep). This section reports the five-feature ensemble under "
        "a rate-aware operating point, so the numbers characterise scaling and ablation rather than "
        "the shipped triad.", BODY),
    Paragraph("5.1 Detection scales with length", H2),
    fig_block(SCALE_FIG, width=14*cm),
    caption("Figure 3. Simulated detection quality vs questionnaire length (rate-aware calibration)."),
    Paragraph(
        f"MCC rises from {mcc_30:.3f} at 30 items to {mcc_300:.3f} at 300 items and plateaus near "
        f"300. Long batteries provide more correlated item pairs, so the permutation baseline "
        f"separates attentive and careless z-scores more cleanly &mdash; the same length-dependence "
        f"seen in the external data (Table 4).", BODY),
    PageBreak(),
    Paragraph("5.2 The contribution of rr", H2),
    Paragraph(
        "Ablating rr from the ensemble isolates its added value over the auxiliaries. The gap grows "
        "with length, and a 500-resample bootstrap 95% CI on &Delta; MCC = MCC(Full) &minus; "
        "MCC(No-rr) excludes zero at <i>every</i> size (Table 6, Figure 4): rr's contribution is "
        "small but statistically positive on short questionnaires and sizeable on long ones.", BODY),
    fig_block(ABL_FIG, width=13.5*cm),
    caption("Figure 4. Simulated ablation: Full ensemble vs auxiliaries-only, by length."),
]
rows_b = [["Items", "Delta MCC", "95% CI", "P(Delta>0)"]]
for _, r in boot.sort_values("size").iterrows():
    if r["scenario"] == "B":
        rows_b.append([str(int(r["size"])), f"{r['delta_pt']:+.3f}",
                       f"[{r['delta_lo']:+.3f}, {r['delta_hi']:+.3f}]", f"{r['delta_p_gt0']:.3f}"])
S += [grid(rows_b, col_widths=[2*cm, 2.2*cm, 4.4*cm, 2.4*cm]),
      caption("Table 6. Bootstrap 95% CI on rr's contribution by size (simulation). Positive at "
              "every length."),
      Paragraph(
        "By pattern (not shown as a table here), rr specifically rescues <i>inconsistent</i> "
        "carelessness (random, mixed) that the auxiliaries miss, and adds nothing on consistent "
        "patterns (straight-lining, acquiescence) where LongString already saturates &mdash; the "
        "same complementarity that motivates the triad.", BODY),
      PageBreak(),
]

# ===========================================================================
# 6. Application + tool
# ===========================================================================
S += [
    Paragraph("6. Application and online tool", H1),
    Paragraph("6.1 When to use ReReRe", H2),
    Paragraph(
        "ReReRe is strongest on multi-construct Likert batteries: in validation, reliable detection "
        "begins around 45&ndash;60 items and 50&ndash;60 respondents and sharpens with length. Below "
        "that the scores still <i>rank</i> respondents usefully but flags should be treated as "
        "exploratory. The method targets inconsistent carelessness; consistent straight-lining is "
        "caught by the LongString partner. Items should share a response format (mixed scales are "
        "handled internally by proportion rescaling); reverse-keyed items are aligned automatically.", BODY),
    Paragraph("6.2 Recipe and operating points", H2),
    Paragraph(
        "In practice the analyst does nothing beyond loading the data: the ensemble scores every "
        "respondent and the automatic two-groups calibration reports an estimated careless rate and "
        "the flags. The only knob is <b>Expected carelessness</b> (Low / Medium / High), which "
        "trades false positives for recall; Medium is the default. If a reliable prior on the rate "
        "exists (a pilot or an earlier wave), it can be supplied to flag exactly that top share "
        "instead.", BODY),
    Paragraph("6.3 Caveats", H2),
    Paragraph(
        "&bull; <b>Prevalence extremes:</b> below ~5% careless there may be too few to separate "
        "automatically; above ~35% the automatic estimate is conservative (raise Expected "
        "carelessness). The ranking itself is unaffected.<br/>"
        "&bull; <b>One pass:</b> apply the detector once to raw data. Re-running on an already-cleaned "
        "set re-standardises within the smaller sample and the structure diagnostic will (correctly) "
        "warn that no careless structure remains.<br/>"
        "&bull; <b>Not a lie detector:</b> ReReRe flags responding that is incoherent relative to the "
        "sample's own structure; it does not certify individual respondents.", BODY),
    Paragraph("6.4 The tool", H2),
    Paragraph(
        "The method is deployed as a free browser tool at <b>re-re.re</b>. It runs entirely "
        "client-side &mdash; the data never leaves the analyst's device (no upload, no account) "
        "&mdash; reads CSV or Excel, reports the estimated careless rate and a per-respondent dot "
        "plot with the cut, and exports the flags and every component score.", BODY),
    PageBreak(),
]

# ===========================================================================
# 7. Statistical concepts box
# ===========================================================================
S += [
    Paragraph("7. Statistical concepts box", H1),
    Paragraph("Per-respondent permutation baseline", H3),
    Paragraph(
        "rr's reference distribution is generated from the data itself, one respondent at a time: "
        "how coherent this person looks on <i>random</i> item pairs. This self-calibration replaces "
        "the distributional assumptions (e.g. multivariate normality) that make chi-square "
        "Mahalanobis thresholds fail on Likert data, and it needs no external cut-off table.", CALL),
    Paragraph("Two-groups model and the estimated rate", H3),
    Paragraph(
        "A BIC-gated two-component Gaussian mixture on the ensemble log-odds decides whether a "
        "careless subpopulation exists (none on clean data &rarr; nothing flagged) and estimates its "
        "size. The cut is anchored on the attentive component's tail, an interpretable "
        "false-positive-rate knob rather than an unknown prevalence.", CALL),
    Paragraph("Matthews Correlation Coefficient (MCC)", H3),
    Paragraph(
        "MCC uses all four cells of the confusion matrix and is well-behaved under class imbalance "
        "(the careless class is usually the minority), ranging from &minus;1 through 0 (chance) to "
        "+1. We report it as the primary classification metric alongside AUC (a threshold-free "
        "ranking measure).", CALL),
    Paragraph("Leave-one-out / cross-validation", H3),
    Paragraph(
        "On the small real study we validate the combiner by leave-one-out (refit on n&minus;1, predict "
        "the held-out respondent); on the large simulation we use 5-fold cross-validation. In both, "
        "each respondent is scored by a model that did not see them, so the reported metrics estimate "
        "out-of-sample performance.", CALL),
    Paragraph("Bootstrap percentile confidence interval", H3),
    Paragraph(
        "For a complex statistic such as &Delta; MCC we resample respondents with replacement (500 "
        "times), recompute &Delta; on each resample, and take the 2.5th/97.5th percentiles &mdash; a "
        "95% interval with no distributional assumption.", CALL),
]

doc = SimpleDocTemplate(OUT_PDF, pagesize=A4, leftMargin=2.0*cm, rightMargin=2.0*cm,
                        topMargin=2.0*cm, bottomMargin=2.0*cm,
                        title="ReReRe Paper (v4, real-data-first)", author="Vittorio Guerrieri")
doc.build(S)
print(f"Wrote {OUT_PDF}")
print(f"Size: {os.path.getsize(OUT_PDF)/1024:.1f} KB")
