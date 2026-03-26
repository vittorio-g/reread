---
project_name: "ReReReRe — Careless Respondent Detection via Permutation-Based Individual Correlation"
project_type: research_paper
status: active
priority: 4
urgency: 3
completion_percent: 80
last_updated: "2026-03-27"
description: "Multiverse simulation study evaluating a permutation-based method (ReReReRe) for detecting careless respondents in questionnaire data."
language: en
tags:
  - "careless-responding"
  - "multiverse-analysis"
  - "simulation-study"
  - "r-language"
  - "psychometrics"
collaborators:
  - "Vittorio"
next_steps:
  - "Re-run full multiverse simulation with reverse-coding fix to find new optimal z threshold"
  - "Create publication-quality figures (heatmaps, crossover plot, AUC curves)"
  - "Write paper draft with corrected results and updated default recommendations"
  - "Consider Johnson IPIP-NEO-300 inject-and-detect validation (stress-test reverse coding)"
---

# ReReReRe — Working Notes

## Current Files (root)

| File | Function | Purpose |
|------|----------|---------|
| `Synthetic_Good_Responses_2.R` | `simulated_good_responses()` | Generate clean CFA-based questionnaire data |
| `Careless_machine_2.R` | `inject_careless()` | Inject careless responses (random/longstring/mixed) |
| `A_good_careless_dataset.R` | `careless_corruption()` | SPI-specific wrapper (not used in multiverse) |
| `ReReReRe.R` | `ReReReRe()` + `rowCor_abs()` | Core detection method (vectorized, min_pairs, align_signs) |
| `Multiverse_Runner.R` | — | Pilot simulation runner (single replication) |
| `Explore_Results.R` | — | Loads results, computes summaries, generates plots |
| `Stratified_Analysis.R` | — | Stratified analyses by corruption level, pattern type, metric comparison |
| `Diagnostic_ReReReRe.R` | — | Quick score distribution check (run before multiverse) |
| `Benchmark_Comparison.R` | — | Standalone benchmark comparison (pilot, 72 conditions) |
| `Final_Multiverse_Runner.R` | — | Publication-quality multiverse: R=50, integrated benchmarks, 3 CSVs |
| `analyze_final_results.py` | — | Python analysis of final results (pandas, summary tables) |
| `mahalanobis_practical_analysis.py` | — | Oracle vs practical (chi-square) Mahalanobis threshold analysis |
| `mahalanobis_practical_threshold.R` | — | Re-run Mah with chi-square thresholds on same seeds |
| `External_Validation.R` | — | ReReReRe vs practical Mahalanobis on real datasets with known careless respondents |
| `PISA_Validation.R` | — | PISA 2018 validation (screen-time C/IER ground truth, per-country) |

## External Datasets (`external_datasets/`)

Downloaded 2026-03-25 for validation on real data. Full details in `external_datasets/README.md`.

### Validated datasets (used in External_Validation.R)

| # | Dataset | N | Items | Factors | GT Type | RR AUC | Mah AUC | RR MCC | Mah MCC | Status |
|---|---------|---|-------|---------|---------|--------|---------|--------|---------|--------|
| 06 | Schneider QoL | 1649 | 31 | 4 | Latent class | **0.821** | 0.724 | **0.361** | 0.208 | Best result |
| 01 | Schroeders 2022 | 605 | 60 | ~6 (HEXACO) | Experimental induction | **0.620** | 0.537 | **0.170** | 0.045 | Marginal (nF~6) |
| 08 | Niessen 2016 | 230 | 100 | 5 (Big Five) | Experimental (speed) | 0.551 | **0.585** | 0.066 | **0.159** | Dead zone (nF=5) |

No longstring pre-screening applied in current results. RR wins AUC on 2/3, Mah wins on Niessen.

**AUC direction bug (fixed 2026-03-25):** pROC `direction` parameter was swapped in earlier
runs — `"<"` was used for z_score (should be `">"`) and `">"` for D^2 (should be `"<"`).
Previous AUC values in CLAUDE.md were 1-true_AUC. All values above are corrected.

**Data format notes:**
- Schroeders: `data_mod_resp.csv`, semicolon-delimited. `Careless` (0/1), `HE01_*` items (already
  reverse-recoded, 5-point). 605 rows, 40% careless (experimental). Some items have `_r` suffix.
- Schneider: `carersp.csv`, standard CSV. `c01`=1 means careless (latent class, 7.4%). Items:
  `dep01-dep08`, `pain01-pain08`, `cog01-cog08`, `fat01-fat07` (31 items, 4 domains, 5-point).
- Niessen: `Raw_data.sav` (SPSS, `haven::read_sav()`). Items `E1-E100` (Big Five, interleaved
  by prefix E/C/A/N/O, 5-point). GT: `Conditon` (0/1, speed manipulation). `Use_me` filters
  valid cases.

### Rejected datasets (with rationale)

| # | Dataset | Reason |
|---|---------|--------|
| 05 | Brühlmann 2020 | CrowdFlower "quality" label ≠ response-level carelessness. Both methods at chance. |
| 09 | Robie HEXACO 2022 | GT is multi-method composite including Mahalanobis → circular. But invaluable for discovering sign-alignment bug. |
| 10 | Arias 2020 | Only 36 items across 3 constructs. Too few factors for RR. |
| 07 | gjkev | Only 45 complete cases out of 132 (n < p). |
| 04 | Bloy 2025 | Only 17 gibberish items. |
| 07 | Kuang 2025 | Very short clinical scales (7-14 items each). |

### PISA 2018 (tested, not included — construct mismatch)

| # | Dataset | N | Items | Factors | GT Type | RR AUC | Mah AUC |
|---|---------|---|-------|---------|---------|--------|---------|
| 11 | PISA 2018 ITA | 4078* | 135 | 36 scales | Screen-time C/IER weights | 0.350 | **0.773** |

*After longstring removal. RR signal is **inverted**: careless (fast/acquiescent) respondents
have HIGHER z-scores because they produce more within-scale consistency than attentive
respondents who genuinely vary. RR detects *inconsistent* carelessness (random/mixed), not
*consistent* carelessness (acquiescence/straight-lining). PISA's dominant careless pattern
is the latter. Mahalanobis works because outlier distance captures both types.

Z-scoring items fixes inversion (AUC→0.60) but damages homogeneous-scale datasets
(Schneider 0.83→0.67). Min-max [0,1] rescaling partially helps (AUC→0.43) but is fragile
to data entry errors. Decision: do NOT build rescaling into algorithm. Document "ensure items
use the same response scale" as preprocessing recommendation.

**Files:** `STU/CY07_MSU_STU_QQQ.sav` (1.8GB), `TIM/CY07_MSU_STU_TIM.sav` (581MB).
612K students, 79 countries. ITA: 11,785 rows, 5,048 complete Likert cases.
C/IER weights computed per-scale via `mclust` Gaussian mixture on log response times
(Ulitzsch et al. 2023 method). 48 scales, 135 Likert items found (59 missing from ITA).

### Johnson IPIP-NEO-300 (Strategy 2 — inject and detect)

| # | Dataset | N | Items | Factors | GT |
|---|---------|---|-------|---------|-----|
| 12 | Johnson 2005 | 20,993 | 300 | 30 facets | None (cleaned; 1,455 straightliners already removed) |

**Purpose:** Stress-test align_signs with heavy reverse coding (148/300 items). Inject
careless respondents at various rates, run RR, compare detection against simulation results
for nF=30. Same dataset used by Welz & Alfons (2023) for careless onset detection.
**File:** `ipip20993.sav` (8.6MB SPSS).

## Pipeline

1. `simulated_good_responses(nConstructs, nItems, n)` — lavaan CFA, random loadings (0.2–0.8), random factor cors (0–0.8, PD), → 7-point Likert
2. `inject_careless(data, pct_careless, ...)` — 3 patterns, corruption levels 50%–100%
3. `ReReReRe(data, corProp, cutOff, iterations, min_pairs, align_signs)` — coupled correlation vs random permutations → z-score (primary), percentile, indCors

## Evaluation Design

Old-style validation: all injected careless respondents count as positive class.

- EVAL_THRESHOLD = 0.01 (catches all careless with corruption >= 50%)
- Careless corruption levels: 50%, 60%, 70%, 80%, 90%, 100%
- Clean respondents and non-careless = negative class

**Primary metric: MCC** (Matthews Correlation Coefficient) — uses all four quadrants (TP/TN/FP/FN), handles class imbalance, doesn't ignore specificity like F1 does. F1 also tracked.

## Key Design Choices

- **Custom loop** (not `multiverse` R package) — more control, checkpoint/resume, progress reporting
- **Z-score scoring** — `(coupled_cor - mean_random) / sd_random`; after longstring bug fix, z-score beats percentile for nF>=10 (Cohen's d up to 3.1 vs 1.7)
- **Fixed z_threshold flagging** — `flag if z_score <= z_threshold`; normative MAD-based flagging incompatible (MAD=0 at percentile ceiling, poor results for delta/z)
- **z_threshold applied post-hoc** — ReReReRe returns raw z-scores; different thresholds are free
- **MCC over F1** — F1 ignores specificity; MCC balances all four confusion matrix quadrants
- **Incremental CSV output** — with checkpoint file for crash recovery; safe_write_table/safe_saveRDS wrappers for OneDrive file-lock resilience
- **min_pairs=15** — prevents degenerate individual-level correlations with small questionnaires
- **align_signs=TRUE** (default) — automatically aligns reverse-coded items using the sign of
  the sample-level correlation. Without this, opposite-keyed item pairs cancel each other in
  `rowCor_abs`, producing z≈0 for everyone. Discovered on Robie HEXACO-60 (52% opposite-keyed
  coupled pairs). Emits warning when disabled with negative-sign pairs present.
- **n < p flagged** — worst case is 4 factors x cycling [10,6,3] = 29 items with n=50
- **RR detects inconsistent carelessness only** — random responding, mixed patterns, cross-factor
  incoherence. It does NOT detect consistent carelessness (acquiescence, straight-lining) which
  requires a separate longstring detector as preprocessing. Demonstrated on PISA 2018 where
  fast/acquiescent respondents had HIGHER z-scores than attentive ones.
- **Do not mix response scales** — items with different Likert ranges (e.g., 4-point mixed with
  6-point) can bias the individual-level correlation. Z-scoring fixes this but damages
  homogeneous-scale data. Recommend users ensure uniform response scale as preprocessing.

## Final Multiverse Design (as run, 2026-03-21)

Script: `Final_Multiverse_Runner.R`

| # | Variable | Levels | Count |
|---|----------|--------|-------|
| V1 | nFactors | 4, 6, 8, 10, 12, 15, 20, 25 | 8 |
| V2 | n_respondents | 50, 100, 200, 300, 500, 1000 | 6 |
| V3 | pct_careless | .05, .10, .25, .50 | 4 |
| V4 | corProp | .03, .05, .07, .10, .15, .20 | 6 |
| V5 | z_threshold | 0, 0.5, 1, ..., 5 | 11 (post-hoc) |
| R | replications | 50 | |

**Fixed:** careless types 1/3 each, careless levels 50-100%, iterations=100, min_pairs=15

**Simulation parameters (relaxed, 2026-03-19):**
- Factor loadings: U(0.2, 0.8) — includes weak items common in real questionnaires
- Factor correlations: max 0.8 — allows high redundancy, avoids near-singular matrices

**Data conditions per rep:** 8 x 6 x 4 = 192
**ReReReRe calls per rep:** 192 x 6 (corProp) = 1,152
**Total ReReReRe calls:** 1,152 x 50 = **57,600**
**Total RR result rows:** 57,600 x 11 (z_threshold) = **633,600**
**Benchmark rows:** 192 x 50 x 4 methods = **38,400**

**Output files:**
1. `final_multiverse_results.csv` (131.6 MB) — ReReReRe results per corProp x z_threshold
2. `final_multiverse_results_by_type.csv` (119.3 MB) — detection by careless pattern type
3. `final_benchmark_results.csv` (3.9 MB) — benchmark methods (oracle threshold, once per dataset)
4. `mahalanobis_practical_results.csv` (1.3 MB) — Mahalanobis with oracle + chi-square thresholds

**Design rationale:** Benchmarks integrated directly into the runner (computed once per dataset,
before the corProp loop) so they run on identical simulated data as ReReReRe. Datasets are not
saved, so this was essential for fair comparison.

## Final Multiverse Results (2026-03-26, 50 replications, align_signs=TRUE)

*Current run includes align_signs fix. Previous runs (Mar 21 without align_signs, Mar 12 with
stronger loadings) archived below for comparison.*

### ReReReRe Overall

| Metric | Value | Previous (Mar 21) | Mar 12 |
|--------|-------|----|----|
| Mean MCC | 0.212 | 0.254 | 0.328 |
| Median MCC | 0.189 | 0.236 | 0.307 |
| Max MCC | 0.836 | 1.000 | 1.000 |
| Mean AUC | 0.715 | 0.750 | 0.803 |
| Cells MCC >= 0.3 | 26.5% | 35.1% | 51.7% |
| Cells MCC >= 0.5 | 5.5% | 8.1% | 20.9% |
| Cells MCC >= 0.7 | 0.3% | 0.5% | 4.2% |
| Cells MCC < 0 | 5.5% | 3.3% | 3.1% |

Performance ~0.04 lower than Mar 21 due to align_signs raising the random pair baseline.
The fix is necessary for correctness on real data with reverse-coded items.

### nFactors: strongest driver

| nF | Mean MCC | Median MCC | SD | Previous (Mar 21) |
|----|----------|------------|-----|---|
| 4 | 0.096 | 0.089 | 0.089 | 0.103 |
| 6 | 0.124 | 0.117 | 0.097 | 0.136 |
| 8 | 0.172 | 0.164 | 0.116 | 0.196 |
| 10 | 0.211 | 0.201 | 0.133 | 0.246 |
| 12 | 0.225 | 0.217 | 0.141 | 0.268 |
| 15 | 0.255 | 0.248 | 0.158 | 0.313 |
| 20 | 0.290 | 0.284 | 0.185 | 0.371 |
| 25 | 0.311 | 0.300 | 0.199 | 0.402 |

### n_respondents: second driver

| n | Mean MCC | Previous (Mar 21) |
|---|----------|----------|
| 50 | 0.088 | 0.166 |
| 100 | 0.159 | 0.219 |
| 200 | 0.221 | 0.260 |
| 300 | 0.248 | 0.279 |
| 500 | 0.268 | 0.293 |
| 1000 | 0.280 | 0.297 |

n=1000 adds almost nothing over n=500 (+0.012 MCC). Diminishing returns plateau confirmed.

### Benchmark Comparison (on identical data, oracle best threshold)

| Method | Mean MCC | Mean AUC | Previous MCC |
|--------|----------|----------|---|
| Mahalanobis | 0.361 | 0.642 | 0.443 |
| IRV | 0.290 | 0.471 | 0.246 |
| LongString | 0.132 | 0.517 | 0.089 |
| PersonTotal | 0.025 | 0.318 | 0.026 |

Mahalanobis hurt most by weaker loadings (−0.082). IRV and LongString actually improved (+0.04
each) — weaker items create more variability that IRV picks up, and wider loading ranges make
longstring patterns more distinct against heterogeneous items.

### ReReReRe vs Mahalanobis: crossover at nF=10

Head-to-head using oracle best MCC per cell (optimizing over corProp and z_threshold for RR):

| nF | RR best MCC | Mah MCC | RR wins | RR AUC | Mah AUC |
|----|-------------|---------|---------|--------|---------|
| 4 | 0.217 | **0.328** | 10% | **0.652** | 0.642 |
| 6 | 0.254 | **0.329** | 21% | **0.687** | 0.635 |
| 8 | 0.319 | **0.356** | 37% | **0.741** | 0.642 |
| 10 | **0.370** | 0.367 | 50% | **0.777** | 0.651 |
| 12 | **0.397** | 0.365 | 59% | **0.795** | 0.646 |
| 15 | **0.447** | 0.359 | 76% | **0.823** | 0.633 |
| 20 | **0.528** | 0.395 | 83% | **0.858** | 0.646 |
| 25 | **0.572** | 0.391 | 86% | **0.873** | 0.641 |

**MCC crossover at nF=10** (50% tie). **AUC crossover at nF=8** (0.741 vs 0.642).
Crossover shifted from nF=8 to nF=10 due to align_signs raising the random baseline.
RR is still a better ranker (AUC) earlier than classifier (MCC).

Notable: at nF=4, RR AUC (0.652) already exceeds Mah AUC (0.642) even though MCC favors Mah.

### Detection by careless pattern type (z_threshold 2.0-3.5)

All three patterns detectable with balanced rates:
- Longstring: 82.7% mean detection (was 78.1%)
- Mixed: 80.7% (was 75.9%)
- Random: 78.6% (was 73.6%)

Detection rates actually improved with weaker loadings — likely because the permutation baseline
better calibrates against noisier data, producing cleaner z-score separation.

### Replication stability

- Within-condition SD of MCC: mean=0.077, median=0.061 (was 0.084/0.067)
- SE of condition mean with R=50: ~0.011 (was ~0.012)
- 95th percentile SD: 0.176 (was 0.191)

Slightly more stable than previous run. 50 replications provide tight CIs.

## Mahalanobis: Oracle vs Practical Thresholds

### How Mahalanobis is used in practice

Standard practice (Meade & Craig, 2012; Curran, 2016): flag respondents whose D^2 exceeds the
chi-square critical value at alpha=.001, where df = number of items (p). This is a fixed formula
requiring no data-driven tuning — the Mahalanobis equivalent of ReReReRe's z=2 default.

D^2 ~ chi-square(p) under multivariate normality. But Likert-scale questionnaire data violates
the normality assumption, making the chi-square reference distribution approximate.

### Practical Mahalanobis results (2026-03-21)

Re-ran Mahalanobis on identical datasets (`mahalanobis_practical_threshold.R`, same seeds):

**Overall:**
| Threshold | Mean MCC | Previous (Mar 12) |
|-----------|----------|---|
| Oracle | 0.361 | 0.443 |
| chi-sq(.001) | **0.056** | 0.054 |
| chi-sq(.01) | **0.116** | 0.128 |

**By nFactors (chi-square alpha=.001):**
| nF | Oracle MCC | Practical MCC | Sensitivity | Specificity |
|----|------------|---------------|-------------|-------------|
| 4 | 0.328 | 0.026 | 0.004 | 1.000 |
| 6 | 0.329 | 0.041 | 0.007 | 1.000 |
| 8 | 0.356 | 0.056 | 0.012 | 1.000 |
| 10 | 0.367 | 0.062 | 0.013 | 1.000 |
| 12 | 0.365 | 0.062 | 0.014 | 1.000 |
| 15 | 0.358 | 0.065 | 0.016 | 1.000 |
| 20 | 0.395 | 0.068 | 0.018 | 1.000 |
| 25 | 0.391 | 0.069 | 0.019 | 1.000 |

**Practical Mahalanobis remains essentially non-functional**: sensitivity of 0.4-1.9% means it
catches almost no careless respondents. Results virtually identical to March 12 run — the
chi-square threshold problem is independent of loading strength.

### Three-level comparison (final, 2026-03-26)

| Comparison | RR (fixed) | Mahalanobis | Winner |
|:--|:--|:--|:--|
| Fixed defaults vs practical chi-sq (.001) | 0.271 | 0.056 | **RR by 4.8x** |
| Fixed defaults vs oracle | 0.271 | 0.361 | Mah (but oracle is unrealistic) |
| Oracle vs oracle (MCC) | varies | 0.361 | Crossover at nF=10 |
| AUC vs AUC (threshold-free) | 0.586-0.836 | ~0.642 | RR from nF=8 |

### AUC: the fairest comparison (completely threshold-free)

Since both methods face threshold uncertainty in practice, AUC is the most honest metric:

| nF | RR AUC (fixed) | Mah AUC | Difference |
|----|----------------|---------|------------|
| 4 | 0.586 | **0.642** | -0.056 |
| 6 | 0.632 | **0.635** | -0.004 |
| 8 | **0.699** | 0.642 | **+0.057** |
| 10 | **0.748** | 0.651 | **+0.097** |
| 12 | **0.769** | 0.646 | **+0.123** |
| 15 | **0.797** | 0.633 | **+0.164** |
| 20 | **0.826** | 0.646 | **+0.180** |
| 25 | **0.836** | 0.641 | **+0.195** |

RR AUC with fully fixed defaults (corProp=0.05) beats Mahalanobis AUC from nF=8.
AUC crossover shifted from nF=6 (Mar 21) to nF=8 due to align_signs. Mahalanobis AUC
remains flat at ~0.64. ReReReRe AUC scales from 0.59 to 0.84.

### Paper framing for threshold comparison

The practical threshold results remain the paper's most striking finding. Standard Mahalanobis
practice (chi-square thresholding) is essentially non-functional for Likert data, with
sensitivity of 0.4-1.9%. This is because the chi-square(p) reference distribution assumes
multivariate normality, which Likert data violates badly. The oracle MCC=0.361 that made
Mahalanobis look competitive was entirely dependent on having access to the true labels —
something no researcher ever has.

ReReReRe with fixed defaults (corProp=0.05, z=1.5) outperforms practical Mahalanobis by 4.4x
across all conditions. Even at nF=4 (where RR is weakest), RR MCC=0.108 vs practical Mah=0.026.

## Practical Recommendations

### Fixed defaults for users

The multiverse explored corProp and z_threshold across wide ranges, but **users should not
need to tune these**. The permutation baseline is self-calibrating, so robust defaults work
across conditions. The paper should recommend:

- **corProp = 0.05** — best at nF>=8 where RR is recommended; corProp=0.03 wins for nF>=12
  but 0.05 is more robust across the full range
- **z_threshold = 1.5** — best overall (MCC=0.247) and in recommended zone (0.336) and
  sweet spot (0.435). Previous z=2.0 recommendation was pre-align_signs; the fix raises
  the random baseline, compressing z-scores and shifting the optimal threshold down.
  Per-nF best z: 1.0 for nF=4-12, 1.5 for nF=20, 2.0 only for nF=25.

z=1.5 means "the respondent's coupled correlation is 1.5 SDs above their own random baseline."
Lower than the pre-align_signs z=2.0 because sign alignment raises random pair correlations.

The multiverse validates that these fixed defaults perform well, not that users should navigate
the parameter space. The method is designed for non-technical users who need a single function
call with sensible defaults.

### Fixed-defaults performance (corProp=0.05, z_threshold=1.5)

| Zone | MCC | Sensitivity | Specificity | Previous (Mar 21, z=2.0) |
|------|-----|-------------|-------------|---|
| All conditions | 0.247 | 0.733 | 0.575 | 0.300 |
| Recommended (nF>=10, n>=100) | 0.336 | 0.718 | 0.715 | 0.409 |
| Sweet spot (nF>=15, n>=300) | 0.435 | 0.789 | 0.793 | 0.518 |

Note: the z=2.0 defaults on the new data give: All=0.271, Rec=0.380, Sweet=0.493.
Both z=1.5 and z=2.0 are viable; z=1.5 is slightly better overall but z=2.0 is better
at nF=25. The difference is small enough that either recommendation is defensible.

### When to use ReReReRe

- Questionnaires with >=10 factors (~63+ items): outperforms Mahalanobis on all metrics
- Sweet spot: >=15 factors (~95+ items), n>=300
- Best for: large-scale surveys, multi-instrument batteries, online data collection
- Even for shorter instruments: still vastly outperforms practical Mahalanobis

### When NOT to recommend ReReReRe

- Very small samples (n<=50): insufficient respondent-level precision
- Very short questionnaires (<4 factors): limited signal, but still better than practical Mahalanobis

### Paper framing

Mahalanobis as typically applied (chi-square thresholding) is essentially non-functional for
detecting careless respondents in Likert-scale data (MCC=0.056, sensitivity=1.3%). The oracle
MCC=0.361 that appears in benchmark comparisons depends on post-hoc threshold optimization
using the true labels — something unavailable in practice.

ReReReRe and Mahalanobis operate on fundamentally different principles: Mahalanobis detects
multivariate outliers (distance from centroid), while ReReReRe detects pattern-level corruption
(cross-factor correlation structure). The key advantage of ReReReRe's permutation approach is
self-calibration: the z-score threshold doesn't depend on distributional assumptions.

LongString disadvantage note: our simulated longstrings are non-contiguous (scattered positions),
realistic for scrambled-order questionnaires but disadvantaging the standard LongString detector.
This is fair for our paper's target use case.

## March 12 vs March 21: Impact of Relaxed Loadings

The March 21 run uses more ecologically valid simulation parameters (loadings 0.2-0.8 instead
of 0.4-0.8, factor correlations capped at 0.8 instead of 0.9). Key differences:

1. **All methods drop ~0.07-0.10 MCC** — weaker items = less signal for everyone
2. **Mahalanobis hurt MORE than ReReReRe** — oracle Mah dropped 0.443→0.361 (−0.082); Mah AUC
   dropped ~0.70→~0.64 (more affected by loading quality than RR)
3. **IRV and LongString improved** (+0.04 each) — wider loading range creates more variability
4. **Crossover point unchanged** — MCC crossover at nF~8, AUC crossover at nF=6
5. **Practical Mahalanobis unchanged** — chi-square threshold problem is independent of loadings
6. **Detection rates improved** — by-type rates up ~4-5pp (permutation baseline calibrates better)
7. **Replication stability improved** — lower SDs with weaker loadings
8. **RR's scaling advantage more pronounced** — RR AUC gap over Mah widens further at high nF

The March 21 results are more conservative and more realistic. The paper's narrative is actually
strengthened: under challenging but realistic conditions, ReReReRe's advantage over Mahalanobis
is even more pronounced because Mahalanobis suffers more from weak items.

## External Validation Results (2026-03-25)

### Sign-Alignment Bug Discovery & Fix

**Critical bug found**: when coupled pairs include items with opposite keying directions
(positively-keyed × negatively-keyed within the same factor), the individual-level correlation
computation cancels out. The sample-level `abs(cor())` correctly identifies these as high-|r|
pairs, but `rowCor_abs` computes a single correlation across ALL k pairs. Same-keying pairs
contribute positive values while opposite-keying pairs contribute negative values; these cancel
*inside* `cor()` before `abs()` is applied, giving indCors ≈ 0 for everyone.

**Demonstrated on Robie HEXACO-60**: 52% of coupled pairs had opposite keying. Result: every
respondent got z ≈ 0 (mean indCors = 0.078, mean rand_mean = 0.084). With fix: indCors = 0.299,
z_mean = 3.49, AUC 0.47 → 0.66.

**Fix (align_signs parameter, default TRUE)**: use the sign of the raw (signed) sample
correlation to align all pairs before computing individual-level correlations. For each pair
where the sample-level correlation is negative, flip the B column (multiply by -1). Applied
to both coupled and random pairs. When disabled with negative pairs present, throws a warning.

This is not just a data-preparation issue — it's a fundamental algorithmic requirement for
handling reverse-coded items without manual recoding. The fix is transparent to users.

### Validated Datasets (3 kept, 2 rejected)

| Dataset | Items | Factors | N | GT Type | RR AUC | Mah AUC | RR MCC | Mah MCC |
|---------|-------|---------|---|---------|--------|---------|--------|---------|
| Schneider QoL | 31 | 4 | 1649 | Latent class | **0.821** | 0.724 | **0.361** | 0.208 |
| Schroeders 2022 | 60 | ~6 | 605 | Experimental | **0.620** | 0.537 | **0.170** | 0.045 |
| Niessen 2016 | 100 | 5 | 230 | Experimental (speed) | 0.551 | **0.585** | 0.066 | **0.159** |

**RR wins AUC on 2/3 datasets.** Mah wins both metrics on Niessen (nF=5, dead zone).
Note: These are corrected values after fixing the AUC direction bug (2026-03-25).

**Rejected datasets:**
- Brühlmann 2020: CrowdFlower "quality" label ≠ response-level carelessness. Both methods near chance.
- Robie HEXACO 2022: Ground truth is a multi-method composite that *includes* Mahalanobis → circular.
  (However, the sign-alignment bug was discovered via this dataset, making it invaluable for
  algorithm development.)

### Key observations

1. **Schneider is the strongest result** (AUC=0.821) despite only 4 factors / 31 items.
   The latent-class ground truth from a large internet sample aligns well with what RR detects.

2. **Schroeders** (HEXACO-60, ~6 factors): modest AUC=0.620, consistent with simulated nF=6 zone.
   40% experimentally-instructed careless is unusually high base rate.

3. **Niessen** (IPIP-100, 5 factors): AUC=0.551, Mah wins (0.585). The "Conditon" column
   (speed manipulation) may not be ideal ground truth — instructed fast responders may still
   respond consistently. 100 items but only 5 factors = nF=5 (below RR crossover). RR flagged
   215/230 respondents (93%), suggesting z=2.0 threshold is too liberal for this dataset.

4. **Practical Mahalanobis** (chi-square alpha=.001) remains weak on real data: MCC 0.045-0.208,
   sensitivity 5.5-31%. Confirms simulation finding that chi-square thresholding is too
   conservative for Likert data.

5. **Both AUC bug (direction swap) and sign-alignment bug** were found and fixed during this
   validation. The AUC directions in External_Validation.R were swapped (z_score used "<"
   instead of ">", D^2 used ">" instead of "<").

### Dataset search conclusions

Field-wide data scarcity confirmed: most careless responding papers don't share raw data.
Alfons & Welz (2024) explicitly called for creating an open benchmark repository. Additional
datasets downloaded (Robie HEXACO, Arias 2020) but not suitable for primary validation due to
circularity (Robie) or too few factors (Arias: 3 constructs, 36 items).

## Revision Log

### 2026-03-27 — align_signs: sign flip → proper reverse coding

**Bug identified:** `align_signs` was using sign flip (`B * -1`) instead of proper reverse
coding (`max+1-B`). Sign flip moves values outside the original Likert range, creating
disproportionate deviations from the mean in `rowCor_abs`. This gives negative-correlation
pairs more weight in the individual-level correlation — an artefact of centering, not signal.

**Fix:** replaced `B_coupled * rep(coupled_signs, each=N)` with proper reverse coding
`(max+1) - B` for columns where sample correlation is negative. Uses per-item observed max
(assumes max(observed) = max(scale), valid for Likert data with n >= ~50). Applied to both
coupled and random pairs.

**External validation results (reverse coding, z=2.0 default):**

| Dataset | Items | nF | AUC (rev) | AUC (sign flip) | AUC (no align) | MCC (rev) | MCC (sign flip) |
|---------|-------|----|-----------|-----------------|----------------|-----------|-----------------|
| Schneider QoL | 31 | 4 | 0.696 | 0.821 | 0.822 | 0.068 | 0.361 |
| Schroeders HEXACO | 60 | ~6 | 0.632 | 0.620 | — | 0.213 | 0.170 |
| Niessen IPIP | 100 | 5 | **0.676** | 0.551 | — | **0.364** | 0.066 |

**Niessen dramatically improved** (+0.125 AUC, +0.298 MCC) — IPIP-100 has many reverse-coded
items. Schroeders stable. Schneider dropped because the reverse coding raises rand_mean for
good respondents more than careless (good respondents' random pairs also show structure after
alignment), compressing z-scores. The previous sign-flip result (AUC=0.821) was likely inflated.

**Schneider diagnostic:** 0 negative coupled pairs (no reverse coding in this dataset). The
drop comes entirely from random pair alignment: rand_mean rose from 0.166 to 0.468 for good
respondents vs 0.231 for careless. Cohen's d between groups remains 0.63 — the signal is there
but requires a lower z threshold. Best MCC at z=0.0-0.5 for Schneider (only 4 factors).

**Threshold sensitivity:** optimal z varies by dataset (Schroeders: z=2.0, Schneider: z=0.0,
Niessen: z=1.5). Full multiverse re-run needed to calibrate new default.

### 2026-03-26 — Re-simulation with align_signs, Optimal z Shift, MAD Experiment

**Re-ran full multiverse** with align_signs=TRUE in ReReReRe.R. The fix applies sign alignment
to both coupled AND random pairs, raising the random baseline and compressing z-scores.

**Impact:** Mean MCC dropped 0.254→0.212 (−0.042), AUC 0.750→0.715 (−0.035). MCC crossover
with Mahalanobis shifted from nF=8 to nF=10; AUC crossover from nF=6 to nF=8. Benchmarks
unchanged (they don't use align_signs). The cost is accepted because the fix is necessary
for correctness on real data with reverse-coded items.

**Optimal z_threshold shifted from 2.0 to 1.5.** The mean z for good respondents dropped
(e.g., nF=4: 0.45, nF=25: 5.05 — vs presumably higher before). At nF=4, z=2.0 flags most
good respondents. Per-nF best: z=1.0 for nF≤12, z=1.5 for nF=20, z=2.0 for nF=25. Overall
best is z=1.5 (MCC=0.247 all, 0.336 recommended, 0.435 sweet spot). corProp=0.03 now beats
0.05 for nF>=12, but 0.05 remains more robust across the full range.

**AUC direction bug confirmed and fixed** in External_Validation.R. Both z_score and D^2
directions were swapped (z used "<" instead of ">", D^2 used ">" instead of "<"). All
external validation AUC values corrected.

**MAD-based adaptive thresholding tested and rejected.** The z-score distribution has too
much natural spread (MAD=1.4-2.8) for MAD outlier detection to work. `median(z) - 2*MAD(z)`
lands deeply negative (−2.7 to −0.8), flagging almost nobody (sensitivity 1.6%). Careless
respondents overlap substantially with the lower end of the good respondent distribution,
violating the "tight cluster + distant outliers" assumption. Fixed z=1.5 outperforms MAD
by 9x (MCC=0.213 vs 0.027). Adaptive thresholding remains an open problem for future work.

**Brühlmann dataset dropped.** CrowdFlower "quality" label ≠ response-level carelessness.
Both RR and Mah near chance. Kept Schroeders (RR AUC=0.620), Schneider (0.821), and
Niessen (0.551, dead zone confirmation at nF=5).

### 2026-03-25 — Sign-Alignment Fix, External Validation, PISA Analysis

**Algorithm fix:** Added `align_signs` parameter (default TRUE) to ReReReRe.R. Automatically
detects and corrects reverse-coded item cancellation in the individual-level correlation
computation. When disabled with negative-sign coupled pairs present, emits warning. Fix uses
the sign of the raw sample correlation to align all pair contributions before computing
individual-level cor(). Applied to both coupled and random pairs.

**External validation:** Ran ReReReRe (fixed defaults: corProp=0.05, z=2.0, align_signs=TRUE)
vs practical Mahalanobis (chi-sq alpha=.001) on 3 real datasets. RR wins AUC on all 3.
Fixed two bugs in External_Validation.R: (1) AUC direction swap (z_score and D^2 directions
were backwards), (2) integer overflow in MCC denominator (cast to double).

**Datasets evaluated (corrected AUC values):** Schneider QoL (best: AUC=0.821),
Schroeders HEXACO (AUC=0.620), Niessen IPIP (AUC=0.551 — Mah wins at 0.585, dead zone nF=5).
Rejected: Brühlmann (CrowdFlower "quality" ≠ carelessness, both at chance), Robie (circular GT).
Downloaded but not used: Arias 2020 (too few factors).

**PISA 2018 analysis:** Ran on Italy (ITA) subset: 135 Likert items across 36 scales, N=4078
after longstring removal (>10 consecutive identical). Ground truth: Ulitzsch et al. (2023)
screen-time C/IER attentiveness weights (Gaussian mixture on log response times per scale).
Result: RR signal is **inverted** (AUC=0.35 raw, 0.43 with [0,1] rescaling, 0.60 with
z-scoring). Careless respondents (fast, acquiescent) have HIGHER z-scores because they
produce more within-scale consistency than attentive respondents who genuinely vary answers.
Mahalanobis AUC=0.77 (unaffected by scaling). Conclusion: RR detects *inconsistent*
carelessness (random/mixed responding), not *consistent* carelessness (acquiescence/
straight-lining). PISA's dominant careless pattern is the latter. These are complementary
constructs, not competing measures. PISA not included in final validation — construct
mismatch, not a method failure.

**Scale-mixing investigation:** Tested z-scoring and [0,1] min-max rescaling for mixed
response scales (PISA has 4/5/6-point items). Z-scoring fixes PISA inversion but damages
homogeneous-scale datasets (Schneider AUC 0.83→0.67). Min-max rescaling is less destructive
but fragile to outliers/data entry errors. Decision: do NOT build rescaling into the
algorithm. Document "ensure items use the same response scale" as a preprocessing
recommendation in the paper.

**Johnson IPIP-NEO-300** downloaded (20,993 respondents, 300 items, 30 facets, 148
reverse-coded). No careless labels (cleaned dataset). Suitable for Strategy 2 (inject and
detect) validation, especially to stress-test align_signs with heavy reverse coding.

### 2026-03-21 — Re-simulation Complete (Relaxed Parameters)

Re-ran full multiverse + practical Mahalanobis with relaxed simulation parameters (loadings
0.2-0.8, factor cor max 0.8, n includes 1000). 633,666 RR rows + 38,404 benchmark rows +
9,600 practical Mah rows. All results updated above. Key finding: relative patterns fully
preserved despite ~0.07 MCC drop across board. Crossover unchanged at nF=8. Practical Mah
still non-functional (MCC=0.056). Added safe_write_table/safe_saveRDS wrappers to
Final_Multiverse_Runner.R for OneDrive file-lock resilience.

### 2026-03-19 — Relaxed Simulation Parameters

Changed `Synthetic_Good_Responses_2.R` to use more realistic/challenging data generation:
- Factor loadings: U(0.4, 0.8) → **U(0.2, 0.8)** — includes weak items common in real questionnaires
- Factor correlations: max 0.9 → **max 0.8** — still allows high redundancy but avoids near-singular matrices

Rationale: the 0.4 loading floor was introduced (2026-03-10) before the longstring bug fix
(2026-03-10e). That fix resolved the issue that made low loadings problematic, so the
conservative floor is no longer needed. Lowering to 0.2 and capping factor cor at 0.8 makes
the simulation more ecologically valid (real questionnaires often have weak items and moderate
inter-factor correlations).

### 2026-03-12c — Practical Mahalanobis Threshold Results (CONFIRMED)

Ran `mahalanobis_practical_threshold.R` on identical datasets (same seeds as final multiverse).
Results: chi-square(.001) MCC=0.054 (median 0.000), sensitivity=1.3%, specificity=1.000.
Practical Mahalanobis is essentially non-functional for Likert data. RR fixed defaults (MCC=0.377)
outperform practical Mahalanobis by 7x across all conditions. This is the paper's most striking
finding: the oracle advantage made Mahalanobis look competitive, but no practitioner has that.

### 2026-03-12b — Mahalanobis Practical Threshold Analysis

Created `mahalanobis_practical_analysis.py` to compare oracle vs practical chi-square thresholds.
Key finding: oracle thresholds are 33-39% lower than chi-square(p, .001) in 100% of cells.
Practical Mahalanobis performance would be substantially worse than the oracle MCC=0.443.
AUC comparison (threshold-free) is fairest: RR fixed defaults beat Mah from nF=8.
Mahalanobis AUC is flat at ~0.70 across all nFactors — scaling advantage belongs entirely to RR.

### 2026-03-12 — Final Results Analysis (March 12 run, loadings 0.4-0.8)

Ran `Final_Multiverse_Runner.R` with 8 nFactors, 5 sample sizes (no n=1000), R=50.
Total: 48,000 ReReReRe calls + 32,000 benchmark evaluations. Key findings:
- Mean MCC=0.328, AUC=0.803
- Crossover at nF=8 (MCC), nF=6 (AUC)
- Fixed defaults MCC=0.377 (all), 0.528 (recommended), 0.674 (sweet spot)
- Mahalanobis oracle MCC=0.443, practical MCC=0.054
- All archived as "March 12 run" — superseded by March 21 relaxed-loading run

### 2026-03-10h — Final Multiverse Design

Created `Final_Multiverse_Runner.R` for publication-quality analyses:
- 50 replications per cell (different random datasets) for stable estimates with CIs
- Integrated benchmark methods (Mahalanobis, PersonTotal, LongString, IRV) on identical data
- Three output CSVs: ReReReRe results, by-type detection, benchmark results
- Checkpoint/resume with per-row granularity
- Reproducible seeds: SEED_BASE + (rep-1) x nConditions + condition_idx

### 2026-03-10g — Benchmark Comparison (pilot)

Created `Benchmark_Comparison.R` to compare ReReReRe against standard careless detection methods.
72 conditions (6 nFactors x 4 n x 3 pct_careless), oracle best threshold for each method.
Results showed crossover at nF=10-15 (refined to nF=8 in final run).

### 2026-03-10f — Switch to Z-Score Metric

After fixing the longstring bug (2026-03-10e), re-ran diagnostics and discovered the **metric
comparison reversed**:

| nFactors | d_percentile | d_zscore | d_indCors |
|----------|-------------|----------|-----------|
| 4 | 0.319 | **0.466** | 0.456 |
| 8 | **1.60** | 1.46 | 1.47 |
| 10 | 1.47 | **1.48** | 1.46 |
| 15 | 1.21 | **1.70** | 1.71 |
| 20 (n=300) | 1.44 | **2.18** | 2.21 |
| 20 (n=1000) | 1.69 | **2.40** | 2.37 |
| 25 | 1.73 | **2.96** | 3.11 |

The previous "percentile wins 6/8 conditions" (2026-03-10c-d) was an **artifact of the
contiguous-block longstring bug**. With longstrings implemented correctly (scattered positions),
z-score and indCors dominate for nF>=10.

**