---
project_name: "ReReReRe — Careless Respondent Detection via Permutation-Based Individual Correlation"
project_type: research_paper
status: active
priority: 4
urgency: 3
completion_percent: 80
last_updated: "2026-04-01"
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
  - "Write paper draft with all results (multiverse + external validation)"
  - "Create publication-quality figures (heatmaps, crossover plot, AUC curves)"
  - "Consider additional external datasets if available"
---

# ReReReRe — Working Notes

## Current Files (root)

### Root (core files)

| File | Function | Purpose |
|------|----------|---------|
| `ReReReRe.R` | `ReReReRe()` + `ReReReRe_F()` + helpers | Standard (weighted/coupled auto-switch) + Per-factor variant (EFA-based, experimental) |
| `Synthetic_Good_Responses_2.R` | `simulated_good_responses()` | Generate clean CFA-based questionnaire data |
| `Careless_machine_2.R` | `inject_careless()` | Inject careless responses (random/longstring/mixed) |

### archive/scripts/ (analysis scripts)

| File | Purpose |
|------|---------|
| `Final_Multiverse_Runner.R` | Publication-quality multiverse: R=50, integrated benchmarks, 3 CSVs |
| `Multiverse_Full.R` | Full multiverse: nF×ipf×n×pct×corProp, R=10, auto-z evaluation |
| `Calibration_nF_ipf_Z.R` | 2D nF×ipf→z_threshold calibration (60 cells, 30 reps) |
| `Calibration_nF_to_Z.R` | 1D nF→z_threshold calibration (39 points, 30 reps) |
| `Report_Full_Multiverse.R` | 18-plot comprehensive report of full multiverse |
| `Report_Calibration.R` | 16-plot calibration analysis report |
| `External_Validation_v2.R` | Validation on 3 datasets (corProp=0.03) |
| `External_Validation_Goldammer.R` | Goldammer 2024: 3 experimental studies, BFI-2/IPIP |
| `External_Validation_Johnson.R` | Johnson IPIP-NEO-300: inject-and-detect, 300 items |
| `External_Validation.R` | Original validation (corProp=0.05, superseded by v2) |
| `Benchmark_Comparison.R` | Standalone benchmark comparison (pilot, 72 conditions) |
| `mahalanobis_practical_threshold.R` | Re-run Mah with chi-square thresholds on same seeds |
| `PISA_Validation.R` | PISA 2018 validation (screen-time C/IER ground truth) |
| `Multiverse_Runner.R` | Pilot simulation runner (single replication) |
| `Explore_Results.R` | Loads results, computes summaries, generates plots |
| `Stratified_Analysis.R` | Stratified analyses by corruption level, pattern type |
| `Diagnostic_ReReReRe.R` | Quick score distribution check |
| `A_good_careless_dataset.R` | SPI-specific wrapper (not used in multiverse) |
| `analyze_final_results.py` | Python analysis of final results |
| `mahalanobis_practical_analysis.py` | Oracle vs practical Mahalanobis analysis |
| `Test_Weighted_ReReReRe.R` | Weighted vs standard ReReReRe simulation comparison (18 conditions) |

## External Datasets (`external_datasets/`)

Downloaded 2026-03-25 for validation on real data. Full details in `external_datasets/README.md`.

### Validated datasets (updated 2026-03-28, corProp=0.03)

| # | Dataset | N | Items | Factors | GT Type | RR AUC | Mah AUC | RR MCC z=1.5 | RR MCC auto | Mah MCC | Notes |
|---|---------|---|-------|---------|---------|--------|---------|-------------|-------------|---------|-------|
| 01 | Schroeders 2022 | 605 | 60 | ~6 | Experimental | **0.606** | 0.537 | **0.177** | **0.228** | 0.045 | RR wins both |
| 06 | Schneider QoL | 1649 | 31 | 4 | Latent class | **0.735** | 0.724 | 0.118 | 0.170 | **0.208** | Mah wins MCC |
| 08 | Niessen 2016 | 180 | 100 | 5 | Speed manip | **0.639** | 0.436 | 0.073 | -0.067 | 0.000 | Both weak |
| 13a | Goldammer S1 (BFI-2 bidir) | 291 | 60 | ~8 | Experimental | 0.711 | **0.787** | **0.320** | 0.269 | 0.257 | RR MCC wins |
| 13b | Goldammer S2 (IPIP unidir) | 265 | 60 | ~6 | Experimental | 0.674 | **0.795** | **0.304** | 0.241 | 0.266 | RR MCC wins |
| 13c | Goldammer S3 (longitudinal) | 523 | 60 | ~7 | Experimental | 0.462 | 0.550 | -0.036 | -0.072 | 0.011 | Both fail |
| 12 | Johnson IPIP-300 | 5000 | 300 | 30 | Inject&detect | — | — | 0.63 (mean) | 0.726 | — | Spec=1.000, oracle MCC~0.73 |

**Summary (6 ground-truth datasets + 1 inject-and-detect):** RR AUC > Mah AUC on 4/6 datasets.
RR MCC (z=1.5) > Mah MCC on 4/6. Johnson IPIP-300 (inject-and-detect): oracle MCC=0.73,
auto-z MCC=0.726, specificity=1.000 across all conditions — near-zero false positives with
300 items. Goldammer S3 fails for both methods (67% careless in longitudinal design).
No longstring pre-screening applied.

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

### PISA 2018 (tested, not included — structural limits)

| # | Dataset | N | Items | Factors | GT Type | RR AUC (raw) | RR AUC (0-1) | Mah AUC |
|---|---------|---|-------|---------|---------|--------------|--------------|---------|
| 11 | PISA 2018 ITA | 4078* | 135 | 36 scales | Screen-time C/IER weights | 0.612 | **0.678** | **0.774** |

*After longstring removal. Previous sign-flip implementation gave inverted signal (AUC=0.35);
proper reverse coding fixed this (AUC=0.612). Signal now in correct direction (low z = careless).

**Scale mixing:** PISA has 97 items on 4-point, 27 on 5-point, 5 on 6-point scales. All item
ranges verified correct (min=1, max=expected, computer-based administration). Rescaling items
to [0,1] using theoretical range improves AUC from 0.612 to 0.678 (+0.066), confirming that
mixed scales bias the individual-level correlation.

**Residual gap with Mahalanobis (0.678 vs 0.774):** 135 items across ~36 scales of 2-9 items
each. Many short scales limit the cross-factor signal RR relies on. Additionally, PISA's
dominant careless pattern is acquiescence/speed (consistent carelessness), which RR detects
less well than random/mixed patterns (inconsistent carelessness). These are complementary
constructs. Best RR MCC=0.243 at z=1.0 (raw) vs Mah practical MCC higher.

**Recommendation:** for mixed-scale questionnaires, users should rescale items to a common
range as preprocessing. Do NOT build rescaling into the algorithm — it would damage
homogeneous-scale datasets where it's unnecessary.

**Files:** `STU/CY07_MSU_STU_QQQ.sav` (1.8GB), `TIM/CY07_MSU_STU_TIM.sav` (581MB).
612K students, 79 countries. ITA: 11,785 rows, 5,048 complete Likert cases.
C/IER weights computed per-scale via `mclust` Gaussian mixture on log response times
(Ulitzsch et al. 2023 method). 48 scales, 135 Likert items found (59 missing from ITA).

### Goldammer 2024 (experimental manipulation, BFI-2 / IPIP)

| # | Dataset | N | Items | Factors | GT Type |
|---|---------|---|-------|---------|---------|
| 13a | Study 1 (BFI-2 bidir) | 291 | 60 | ~8 | Experimental (0=honest, 1=33% careless, 2=100%) |
| 13b | Study 2 (IPIP unidir) | 265 | 60 | ~6 | Experimental (same conditions) |
| 13c | Study 3 (longitudinal) | 523 | 60 | ~7 | Experimental (same conditions, t1 only) |

**Source:** ETH Zurich Polybox, downloaded 2026-03-28.
**Results:** RR MCC (z=1.5) beats Mahalanobis on S1 (0.320 vs 0.257) and S2 (0.304 vs 0.266).
S3 fails for both methods — 67% careless rate in longitudinal design, possibly confounded by
practice effects or test-retest variability that masks the careless signal.
**Breakdown S1 (100% vs 33% careless):** 100% careless AUC=0.692, 33% careless AUC=0.738.
Interestingly, 33% partial careless is slightly easier to detect (oracle MCC=0.405 vs 0.398),
possibly because 100% random responding creates more extreme outliers that Mahalanobis catches
better, while partial careless is a more "natural" pattern that RR's correlation-based approach
handles well.

### Johnson IPIP-NEO-300 (Strategy 2 — inject and detect)

| # | Dataset | N | Items | Factors | GT |
|---|---------|---|-------|---------|-----|
| 12 | Johnson 2005 | 20,993 | 300 | 30 facets | None (cleaned; 1,455 straightliners already removed) |

**What it is:** The IPIP-NEO-300 is a 300-item personality questionnaire measuring 30 facets
of the Big Five (6 facets per trait: Extraversion, Agreeableness, Conscientiousness,
Neuroticism, Openness). Collected by John A. Johnson (2005) via a public website where
people could take the personality test for free. Original N=20,993. Johnson already removed
1,455 straightliners, so the dataset is considered "clean" — no careless ground truth exists.

**Why we use it:** (1) **Inject-and-detect strategy** — we take 5,000 clean respondents, inject
artificial careless at various rates, and test if ReReReRe finds them. (2) **Stress-test for
align_signs** — 148/300 items are reverse-coded (~49%), the most extreme case in our validation.
(3) **Largest questionnaire tested** — verifies that simulation predictions for nF=30 hold on
real data. (4) **1-5 Likert scale**, online administration — realistic conditions.
Same dataset used by Welz & Alfons (2023) for careless onset detection.

**File:** `ipip20993.sav` (8.6MB SPSS). Values of 0 = unanswered (recoded to NA).
Complete cases: 7,325 out of 20,993. Sampled to 5,000 for computational feasibility.

**Results (2026-03-30):** Sampled 5,000 complete respondents (0→NA, complete.cases → 7,325,
then random 5,000). Injected careless at 5%, 10%, 20% × 3 reps each. corProp=0.03.

| pct | Mean MCC (oracle, z=3.0) | Mean MCC (z=1.5) | Mean MCC (auto, z=2.37) | Specificity |
|-----|--------------------------|-------------------|-------------------------|-------------|
| 5% | **0.750** | 0.628 | 0.726 | 1.000 |
| 10% | **0.726** | 0.635 | — | 1.000 |
| 20% | **0.702** | 0.607 | — | 1.000 |

Specificity = 1.000 across nearly all conditions — with 300 items, virtually zero false
positives at z≥1.5. Oracle sensitivity ~55-60% at z=3.0. Auto-z chose z=2.37 (appropriate
for 300 items), closing most of the gap to oracle. Results consistent with simulation
predictions for nF=30 (MCC ~0.6-0.8).

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
- **Two-function architecture (2026-04-01)** — `ReReReRe()` is the standard method (weighted
  ≤60 items, coupled >60 items). `ReReReRe_F()` is the per-factor EFA variant (experimental).
  Both EFA-A and EFA-D were tested: they win on simulated data but **lose on real validation
  datasets** (EFA-D AUC=0.561 vs Coupled AUC=0.635, coupled wins 5/6 datasets). However,
  `ReReReRe_F()` improves results on short real questionnaires (Pennycook: r +0.047, R² +0.037).
  Per-factor coherence (per-EFA-factor z-scores) also wins simulation but loses on real data.
  The simulation-to-real gap is caused by parallel analysis producing poor factor assignments
  on real data with cross-loadings and noisy structure.

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

- **corProp = 0.03** — wins overall in full multiverse (MCC=0.352 vs 0.327 for 0.05).
  Fewer but more selective coupled pairs produce cleaner signal. Updated from previous
  recommendation of 0.05 based on full multiverse with items_per_factor varied.
- **z_threshold = 1.5** — best overall (MCC=0.274), essentially tied with z=2.0 (0.274).
  Auto-z calibration (LOESS from parallel analysis) gives identical performance (0.268).
  Any value in the 1.0-2.0 range works well; the method is not threshold-sensitive.

z=1.5 means "the respondent's coupled correlation is 1.5 SDs above their own random baseline."
Lower than the pre-align_signs z=2.0 because sign alignment raises random pair correlations.

**Auto-z option:** available via `auto_z=TRUE` for user convenience. Uses parallel analysis
to estimate nF, then looks up calibrated z from LOESS curve. Performance is identical to
fixed z=1.5 (MCC 0.268 vs 0.274), so it adds no value beyond not requiring a parameter choice.

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

## nF → z_threshold Calibration (2026-03-26, COMPLETED)

### Design rationale

The optimal z_threshold depends strongly on nF. Previous multiverse results showed per-nF best
z ranging from 1.0 (nF≤12) to 2.0 (nF=25), but this was noisy because items/factor varied
via the cycling pattern `c(10, 6, 3)`. To build a clean lookup table for automatic threshold
selection, we run a dedicated calibration with **items/factor fixed at 6**.

### Calibration parameters

| Parameter | Value | Rationale |
|-----------|-------|-----------|
| nF range | 2-40 (step 1) | Full range, fine-grained |
| items/factor | 6 (fixed) | Isolates pure nF effect, realistic middle value |
| z_threshold | 0.1-3.0 (step 0.2) | Post-hoc, covers full useful range |
| n_respondents | 300 | In "recommended" zone |
| pct_careless | 10% | Realistic base rate |
| corProp | 0.05 | Recommended default |
| iterations | 100 | Standard permutation count |
| Replications | 30 | 30 scatter points per nF |
| Careless levels | 50-100% | Standard corruption range |

Total ReReReRe calls: 39 nF × 30 reps = **1,170**. Runtime: ~47 minutes.

### Output files

- `calibration_nF_z_raw.csv` — all reps × nF × z_threshold (full evaluation grid)
- `calibration_nF_z_best.csv` — optimal z per rep per nF (for scatter plot)
- `calibration_nF_z_lookup.csv` — mean optimal z per nF (the actual lookup table)
- `plot_calibration_scatter.png` — scatter cloud with mean overlay, point size = MCC

### Results (2026-03-26)

The optimal z_threshold follows a **U-shaped curve** as a function of nF:

| nF zone | nF range | Mean optimal z | Mean MCC | Interpretation |
|---------|----------|---------------|----------|----------------|
| Dead zone | 2-7 | 0.9-1.2 | 0.09-0.17 | Too few factors, threshold noisy, poor detection |
| Transition | 8-14 | 0.5-0.8 | 0.20-0.34 | z drops as signal emerges but is still weak |
| Sweet spot entry | 15-20 | 0.5-0.7 | 0.35-0.46 | Stable low z, good MCC |
| Scaling zone | 21-28 | 0.7-1.0 | 0.46-0.56 | z rises as z-score distribution separates |
| High nF | 29-40 | 1.3-2.1 | 0.53-0.58 | Strong signal, higher z optimal |

**Key nF breakpoints from lookup table:**
| nF | items | mean_best_z | median_best_z | mean_MCC |
|----|-------|-------------|---------------|----------|
| 4 | 24 | 0.95 | 0.9 | 0.109 |
| 8 | 48 | 0.71 | 0.7 | 0.197 |
| 10 | 60 | 0.79 | 0.8 | 0.255 |
| 12 | 72 | 0.73 | 0.7 | 0.310 |
| 15 | 90 | 0.55 | 0.5 | 0.354 |
| 20 | 120 | 0.63 | 0.7 | 0.455 |
| 25 | 150 | 0.83 | 0.9 | 0.527 |
| 30 | 180 | 1.29 | 1.3 | 0.540 |
| 35 | 210 | 1.83 | 1.8 | 0.542 |
| 40 | 240 | 2.11 | 2.2 | 0.582 |

**U-shape explanation:** At low nF, z-scores are compressed (good and careless respondents
overlap heavily), so a low threshold is needed to catch anyone — but this also catches many
good respondents (low specificity → low MCC). As nF increases, the z-score distributions
separate and the optimal threshold rises because higher z provides better specificity without
losing sensitivity. The minimum of the U (~nF 13-15, z≈0.5) marks where the signal first
becomes reliable enough for meaningful detection.

**Scatter cloud variability:** SD of optimal z ranges from 0.29 (nF=14) to 0.77 (nF=5),
confirming high stochastic noise at low nF. By nF≥25, the scatter tightens (SD≈0.5) and
individual replications cluster more clearly around the mean.

**Scatter cloud plot** (`archive/calibration_nF_z/plot_calibration_scatter.png`): One of the
best figures for the paper. Blue dots = individual replications (jittered, size ∝ MCC), red
dots = mean across reps (size ∝ mean MCC), dashed red line = mean trajectory. The U-shape
is clearly visible: the red line descends from z≈1.0 (nF=2-6) to a minimum of z≈0.5
(nF=13-15), then rises linearly to z≈2.1 (nF=40). Critically, the blue dot SIZE grows
dramatically from left to right — at nF=2-6 all dots are tiny (MCC<0.15) while at nF=35-40
they're large (MCC>0.5). This visually encodes both the threshold AND the quality of
detection improving with more factors.

## Full Multiverse Variable Space

### Dataset variables (simulation)

| # | Variable | Description | Reasonable range |
|---|----------|-------------|------------------|
| D1 | **nFactors** | Number of latent factors | 2-40 |
| D2 | **items_per_factor** | Items per factor | 3, 6, 10, 15, 20 |
| D3 | **n_respondents** | Sample size | 50, 100, 200, 300, 500, 1000 |
| D4 | **pct_careless** | % careless respondents | 5%, 10%, 25%, 50% |
| D5 | **careless_type** | Type (random/longstring/mixed) | 1/3 each (fixed) |
| D6 | **corruption_level** | How much of profile is corrupted | 50%-100% (fixed) |
| D7 | **loading_range** | Factor loading strength | U(0.2, 0.8) fixed |
| D8 | **factor_cor_max** | Max inter-factor correlation | 0.3, 0.5, 0.8 |

### ReReReRe variables

| # | Variable | Description | Reasonable range |
|---|----------|-------------|------------------|
| R1 | **corProp** | Proportion of high-correlation pairs | 0.03, 0.05, 0.07, 0.10, 0.15, 0.20 |
| R2 | **z_threshold** | Z-score flagging threshold | 0.1-3.0 (post-hoc, free) |
| R3 | **auto_z** | Automatic calibration (TRUE/FALSE) | Evaluated post-hoc |
| R4 | **min_pairs** | Minimum coupled pairs | 10, 15, 20, 30 |
| R5 | **iterations** | Number of random permutations | 50, 100, 200 |
| R6 | **align_signs** | Reverse-coded item alignment | TRUE (fixed, mandatory) |

**Notes:**
- z_threshold (R2) is post-hoc: compute once, evaluate at all levels for free
- auto_z (R3) is post-hoc: just look up calibrated z from the LOESS curve
- corProp (R1) is the most expensive: each level requires a separate ReReReRe call
- items_per_factor (D2) was fixed at 6 in calibration — varying it tests LOESS robustness

### Auto-calibration v1: nF-only (2026-03-26)

**Design:** Added `auto_z` parameter to `ReReReRe()`. When `auto_z=TRUE`, runs parallel
analysis to estimate nF, looks up optimal z from LOESS curve (nF=2-40, 39 calibration points).

**Result:** auto-z ≈ fixed z=1.5 (MCC 0.268 vs 0.274). The 1D calibration provides no
meaningful improvement because it ignores items_per_factor, which the full multiverse showed
explains 31% of variance — more than nF itself (26%).

### Auto-calibration v2: total_items (2026-03-27, CURRENT)

**2D calibration** over nF (4-30, 10 levels) × ipf (3-12, 6 levels) = 60 cells, 30 reps each.
Script: `Calibration_nF_ipf_Z.R` | Runtime: 35.5 min | 1,800 RR calls.

**Key finding: total_items is the best single predictor of optimal z.**

Simple linear model: `z = 0.505 + 0.0042 * total_items` (R²=0.476)
vs `z = 0.261 + 0.020*nF + 0.055*ipf` (R²=0.245)

Total_items alone explains twice the variance of nF + ipf separately, confirming that
what matters is the total number of item pairs available for the coupled correlation.

**Performance comparison:**

| Strategy | Mean MCC | vs fixed z=1.5 |
|----------|----------|----------------|
| z fixed = 1.5 | 0.316 | — |
| **Auto-z 2D** | **0.340** | **+0.024** |
| Oracle | 0.380 | +0.064 |

Auto-z 2D closes **38% of the gap** between fixed and oracle. The gain is concentrated
on longer questionnaires where the optimal z diverges most from 1.5:

| Total items | z=1.5 MCC | Auto-z MCC | Gain | Optimal z |
|-------------|-----------|------------|------|-----------|
| 48 | 0.14 | 0.15 | +0.01 | ~0.9 |
| 120 | 0.41 | 0.45 | +0.04 | ~0.6 |
| 200 | 0.66 | 0.66 | ±0.00 | ~1.2 |
| 300 | 0.69 | 0.78 | +0.09 | ~2.4 |
| 360 | 0.65 | 0.82 | +0.17 | ~2.8 |

For questionnaires >200 items, the optimal z rises well above 1.5 and auto-calibration
provides substantial gains. Below 100 items, z≈0.5-1.0 is optimal and the gain is small.

**Optimal z by nF × ipf (from lookup table):**

| nF\ipf | 3 | 6 | 10 | 12 |
|--------|-----|-----|------|------|
| 4 | 1.1 | 0.8 | 1.0 | 0.8 |
| 8 | 1.0 | 1.0 | 0.7 | 0.7 |
| 12 | 1.0 | 0.6 | 0.6 | 0.7 |
| 15 | 0.7 | 0.5 | 0.7 | 1.1 |
| 20 | 0.7 | 0.5 | 1.2 | 1.8 |
| 25 | 0.7 | 0.7 | 1.9 | 2.3 |
| 30 | 0.6 | 1.0 | 2.5 | 2.8 |

The pattern is clear: optimal z forms a **diagonal gradient** from low-z (top-left, few items)
to high-z (bottom-right, many items). This is because with many items, z-score distributions
separate widely between good and careless respondents, allowing a higher threshold for better
specificity without losing sensitivity.

**SD of optimal z** (calibration stability):
- Low total items (<60): SD 0.5-0.7 → very noisy, any z in 0.5-1.5 is roughly equivalent
- Medium (60-150): SD 0.3-0.5 → moderate stability
- High (>200): SD 0.2-0.4 → stable calibration, clear optimal
- nF=30 ipf=12 (360 items): SD=0.17 → extremely stable (z=2.8 ± 0.17)

**MCC at optimal z (best achievable detection):**

| nF\ipf | 3 | 6 | 10 | 12 |
|--------|-------|-------|-------|-------|
| 8 | 0.102 | 0.183 | 0.311 | 0.371 |
| 12 | 0.126 | 0.303 | 0.487 | 0.542 |
| 15 | 0.164 | 0.363 | 0.615 | 0.665 |
| 20 | 0.190 | 0.495 | 0.706 | 0.784 |
| 25 | 0.238 | 0.573 | 0.762 | 0.815 |
| 30 | 0.265 | 0.635 | 0.818 | 0.838 |

At nF=30 ipf=12 (360 items), MCC=0.838 — near-perfect detection.

**Implementation (DONE, 2026-03-28):** Replaced 1D LOESS with 2D total_items-based lookup.
Changes to `ReReReRe.R`:

| Parameter | Before | After | Rationale |
|-----------|--------|-------|-----------|
| corProp default | 0.05 | **0.03** | MCC 0.352 vs 0.327 in full multiverse |
| auto_z engine | 1D (nF → z via LOESS) | **2D (total_items → z)** | R²=0.476 vs 0.245 |
| Calibration table | 39 points (nF=2-40) | **60 points (10 nF × 6 ipf)** | More robust |
| Fallback formula | none | **z = 0.505 + 0.0042 × total_items** | If LOESS fails |
| z clamp range | nF in [2,40] | **z in [0.3, 3.5]** | Prevents extreme thresholds |

When `auto_z=TRUE`, `ReReReRe()` now:
1. Computes total_items = ncol(data)
2. Optionally runs parallel analysis (for nF in output)
3. Looks up z from LOESS fitted on 60-cell 2D calibration grid
4. Falls back to linear formula if LOESS prediction fails

**Output files (in `archive/calibration_nF_ipf_z/`):**
- `calibration_raw.csv` — all reps × cells × z_thresholds
- `calibration_best.csv` — optimal z per rep per cell
- `calibration_lookup.csv` — mean optimal z per nF × ipf
- `calibration_report.txt` — full text report
- `plot_01_heatmap_optimal_z.png` — 2D heatmap of optimal z
- `plot_02_heatmap_mcc.png` — 2D heatmap of MCC at optimal z
- `plot_03_scatter_z_by_total_items.png` — scatter cloud (the key figure)
- `plot_04_scatter_z_by_nF_colored.png` — scatter by nF, colored by ipf
- `plot_05_z_curves_selected.png` — MCC vs z for selected combos
- `plot_06_heatmap_z_variability.png` — SD of optimal z (calibration stability)
- `plot_07_z_vs_total_items_loess.png` — LOESS fit for total_items → z
- `plot_08_fixed_vs_auto_z_mcc.png` — head-to-head comparison

## Full Multiverse Results (2026-03-26, with items_per_factor varied)

Script: `Multiverse_Full.R` | Runtime: 161 minutes | 5,670 RR calls | 85,050 result rows

### Design

| # | Variable | Levels | Count |
|---|----------|--------|-------|
| D1 | nFactors | 4, 8, 12, 16, 20, 25, 30 | 7 |
| D2 | items_per_factor | 3, 6, 10 | 3 |
| D3 | n_respondents | 100, 300, 500 | 3 |
| D4 | pct_careless | .05, .10, .25 | 3 |
| R1 | corProp | 0.03, 0.05, 0.10 | 3 |
| R2 | z_threshold | 0.1-3.0 step 0.2 | 15 (post-hoc) |
| R3 | auto_z | evaluated post-hoc | free |
| Reps | | 10 | |

Fixed: align_signs=TRUE, iterations=100, min_pairs=15, careless types 1/3 each,
corruption 50-100%, loadings U(0.2,0.8).

### Key finding: items_per_factor is a major driver

Oracle best MCC by nFactors × items_per_factor (corProp=0.05, averaged over n and pct):

| nF | ipf=3 | ipf=6 | ipf=10 |
|----|-------|-------|--------|
| 4 | 0.075 | 0.111 | 0.160 |
| 8 | 0.111 | 0.207 | 0.314 |
| 12 | 0.143 | 0.295 | 0.460 |
| 16 | 0.171 | 0.388 | 0.579 |
| 20 | 0.202 | 0.447 | 0.593 |
| 25 | 0.222 | 0.467 | 0.596 |
| 30 | 0.260 | 0.513 | 0.560 |

**Total items matters more than nF alone.** nF=12 with ipf=10 (120 items, MCC=0.46) beats
nF=25 with ipf=3 (75 items, MCC=0.22). The signal comes from the TOTAL number of high-|r|
pairs available, which scales with both nF and ipf.

At nF=30 with ipf=10, MCC slightly drops vs nF=25 — likely because 300 items with n=100-500
approaches the n < p regime where individual-level correlations become noisy.

### corProp: smaller is better

| corProp | Mean MCC (oracle) |
|---------|-------------------|
| 0.03 | **0.352** |
| 0.05 | 0.327 |
| 0.10 | 0.286 |

corProp=0.03 wins overall. Fewer but more selective coupled pairs produce cleaner signal.
This updates the previous recommendation of corProp=0.05.

### Auto-z vs fixed z: the critical comparison

Auto-z (LOESS-calibrated from parallel analysis) vs fixed z thresholds (corProp=0.05):

| nF | auto-z | z=1.0 | z=1.5 | z=2.0 | oracle |
|----|--------|-------|-------|-------|--------|
| 8 | 0.171 | 0.171 | 0.163 | 0.153 | 0.211 |
| 20 | 0.310 | 0.315 | 0.332 | 0.327 | 0.414 |
| 30 | 0.322 | 0.295 | 0.327 | 0.341 | 0.444 |

**Overall means:**
| Strategy | Mean MCC |
|----------|----------|
| auto-z | 0.268 |
| z=1.0 | 0.260 |
| z=1.5 | **0.274** |
| z=2.0 | **0.274** |
| oracle | 0.356 |

**Auto-z ≈ fixed z=1.5 ≈ fixed z=2.0.** The auto-calibration provides no meaningful
improvement over simply using z=1.5 as a universal default. Auto-z correctly adapts (uses
low z at nF=8, matching z=1.0's performance there), but the gain at low nF is offset by
slight losses at high nF.

**Practical implication:** the auto-z feature is "nice to have" for user convenience (no
parameter choice needed), but the paper should emphasize that **z=1.5 is a robust universal
default** and auto-calibration adds complexity without performance gain. The oracle gap
(~0.08 MCC) represents room for improvement, but it requires knowing the true labels —
no threshold selection strategy can close this gap without external information.

### Auto-z by items_per_factor (corProp=0.05)

| nF | ipf=3 | ipf=6 | ipf=10 |
|----|-------|-------|--------|
| 8 | 0.070 | 0.175 | 0.268 |
| 20 | 0.159 | 0.371 | 0.401 |
| 30 | 0.191 | 0.395 | 0.379 |

Auto-z performance tracks items_per_factor closely, confirming that the calibration
(done at ipf=6) generalizes reasonably to ipf=3 and ipf=10.

### Output files (in `archive/multiverse_full/`)

- `multiverse_full_raw.csv` (85,050 rows) — all results per condition × z_threshold
- `multiverse_full_best.csv` (5,670 rows) — oracle best z per cell
- `multiverse_full_auto_z.csv` (2,430 rows) — auto-z performance
- `plot_multiverse_oracle_heatmap.png` — nF × ipf heatmap
- `plot_auto_z_vs_oracle.png` — auto-z vs oracle by nF
- `plot_ipf_effect_auto_z.png` — auto-z by ipf

### Comprehensive Report (18 plots in `archive/multiverse_full/report/`)

Script: `Report_Full_Multiverse.R` — generates full analysis of all variables.

**Overall summary:**
- Mean MCC (oracle best z): 0.322
- Median MCC: 0.272
- Max MCC: 1.000
- Cells MCC >= 0.3: 45.5%
- Cells MCC >= 0.5: 21.5%

**Variable importance (eta-squared from one-way ANOVAs):**

| Variable | eta² | Variance explained |
|----------|------|-------------------|
| **items_per_factor** | **0.311** | **31.1%** |
| **nFactors** | **0.261** | **26.1%** |
| pct_careless | 0.036 | 3.6% |
| corProp | 0.017 | 1.7% |
| n_respondents | 0.010 | 1.0% |

**This is the paper's most important structural finding:** items_per_factor explains MORE
variance than nFactors (31% vs 26%). Previous analyses focused entirely on nF, but ipf is
equally important. The practical recommendation must be framed in terms of **total questionnaire
length** (nF × ipf), not just number of constructs.

**Main effects by variable:**

D1 — nFactors:
| nF | Mean MCC | SD |
|----|----------|-----|
| 4 | 0.117 | 0.063 |
| 8 | 0.207 | 0.108 |
| 12 | 0.294 | 0.158 |
| 16 | 0.368 | 0.191 |
| 20 | 0.405 | 0.199 |
| 25 | 0.423 | 0.217 |
| 30 | 0.438 | 0.220 |

D2 — items_per_factor:
| ipf | Mean MCC | SD |
|-----|----------|-----|
| 3 | 0.169 | 0.090 |
| 6 | 0.341 | 0.174 |
| 10 | 0.455 | 0.221 |

D3 — n_respondents (diminishing returns):
| n | Mean MCC |
|---|----------|
| 100 | 0.290 |
| 300 | 0.335 |
| 500 | 0.340 |

n=500 adds almost nothing over n=300 (+0.005 MCC).

D4 — pct_careless:
| pct | Mean MCC |
|-----|----------|
| 5% | 0.273 |
| 10% | 0.320 |
| 25% | 0.372 |

Higher base rate → higher MCC (more signal to detect, less class imbalance).

R1 — corProp:
| corProp | Mean MCC |
|---------|----------|
| 0.03 | **0.352** |
| 0.05 | 0.327 |
| 0.10 | 0.286 |

R2 — z_threshold (averaged across all conditions):
Best overall: z=1.3 to z=1.5 (MCC≈0.249), very flat plateau from z=0.9 to z=1.9.
The method is remarkably insensitive to z_threshold choice in this range.

**Sensitivity-Specificity trade-off:** at z=1.5, mean sensitivity≈0.73, specificity≈0.58.
Crossover point at approximately z=1.1.

**Z-score separation (good vs careless respondents):**
| nF | Mean z (good) | Mean z (careless) | Gap |
|----|---------------|-------------------|-----|
| 4 | 0.43 | -0.05 | 0.48 |
| 8 | 1.45 | 0.12 | 1.33 |
| 12 | 2.38 | 0.39 | 1.99 |
| 16 | 3.08 | 0.67 | 2.41 |
| 20 | 3.59 | 1.02 | 2.57 |
| 25 | 4.17 | 1.50 | 2.67 |
| 30 | 4.75 | 1.94 | 2.81 |

The gap grows quasi-linearly with nF, explaining the scaling advantage. At nF=4, good and
careless z-scores overlap heavily (gap=0.48). By nF=30, the gap is 2.81 — nearly 3 SD of
separation.

**Total items as unifying predictor:**
| Total items | Mean MCC |
|-------------|----------|
| 12 | 0.075 |
| 48 | 0.188 |
| 72 | 0.291 |
| 120 | 0.443 |
| 160 | 0.558 |
| 200 | 0.580 |
| 300 | 0.550 |

MCC follows a saturating curve: rises steeply to ~150 items, then plateaus. Above 250 items,
MCC can actually decline slightly (n < p regime with n=100-500). The LOESS smooth in
plot 09 shows this clearly.

**Report plots (in `archive/multiverse_full/report/`):**
- `01_mcc_by_nFactors.png` — main effect with errorbars
- `02_mcc_by_items_per_factor.png` — bar chart
- `03_mcc_by_n_respondents.png` — diminishing returns curve
- `04_mcc_by_pct_careless.png` — base rate effect
- `05_mcc_by_corProp.png` — corProp comparison
- `06_mcc_by_z_threshold.png` — z curve with recommended line
- `07_z_curves_by_nFactors.png` — z sensitivity by nF (optimal z shifts)
- `08_heatmap_nF_x_ipf.png` — the key 2D heatmap
- `09_mcc_by_total_items.png` — unifying total_items predictor with LOESS
- `10_nF_x_n_interaction.png` — nF × sample size
- `11_nF_x_pct_interaction.png` — nF × base rate
- `12_nF_x_corProp_interaction.png` — nF × corProp
- `13_auto_z_vs_fixed_vs_oracle.png` — threshold strategy comparison
- `14_sens_spec_tradeoff.png` — sensitivity/specificity crossover
- `15_boxplot_mcc_by_nF.png` — MCC distribution by nF
- `16_heatmap_nF_x_z.png` — full nF × z heatmap with optimal borders
- `17_ipf_within_nF.png` — grouped bars ipf within each nF
- `18_z_separation_by_nF.png` — z-score gap good vs careless

## Revision Log

### 2026-04-01 — Two-function architecture: ReReReRe() + ReReReRe_F()

**Decision:** Split into two public functions instead of a single method with mode switching.

`ReReReRe()` — Standard (default): auto-switches between weighted (≤60 items) and coupled
(>60 items). This is the validated, recommended method for all primary analyses.

`ReReReRe_F()` — Per-factor (experimental): wrapper that calls `ReReReRe(mode="efa_d")`.
Uses EFA to select within-factor pairs, weights by observed |r|. Recommended as secondary
analysis, especially for short questionnaires.

**Motivation:** Per-factor/EFA methods consistently win on simulated data but lose on real
validation datasets with ground truth. However, on the Pennycook & Rand dataset (30 items),
EFA-D improves the paper's key results substantially (r +0.047, R² +0.037) while coupled
makes them worse. The two functions serve complementary purposes.

**Per-factor simulation results (360 conditions, Test_PerFactor_v2.R):**

| Items range | PerFactor mean(z) | Standard coupled | Winner |
|-------------|-------------------|------------------|--------|
| <30 | **0.149** | 0.107 | PF (+39%) |
| 30-60 | **0.256** | 0.197 | PF (+30%) |
| 60-100 | **0.326** | 0.303 | PF (+8%) |
| 100-200 | 0.508 | **0.523** | Std (+3%) |
| >200 | 0.640 | **0.639** | ~Tied |

**Pennycook test (20% flagging cap):**

Study 1 (30 items, N=782): EFA-D r(CRT,Disc) +0.047, R² +0.037. Coupled: −0.009.
Study 2 (24 items, N=2564): Weighted +0.019, EFA-D +0.017. Coupled: −0.002.
All Table 1 correlations improve with EFA-D (e.g., Trump Rep-consistent .197→.270).

**Ground truth validation caveat:** The external datasets have questionable ground truth
quality (63-67% instructed careless in Goldammer, speed-based in Niessen, algorithmic in
Schneider). The Pennycook test uses a different criterion: "do paper results improve?",
which is arguably more ecologically valid.

### 2026-03-31 — EFA-Based Pair Selection: Tested and Rejected (Options A and D)

Tested two EFA-guided pair selection strategies to see if they could replace the coupled/weighted
2-level switch with a single unified method.

**Option A: EFA within-factor pairs + loading product weights (λ_i × λ_j)**

54 conditions (nF=4-20 × ipf=3-10 × 3 reps), N=300, 10% careless.

| Item range | Standard | Weighted | EFA-A | Winner |
|-----------|----------|----------|-------|--------|
| <30 | 0.101 | **0.167** | 0.152 | Weighted |
| 30-60 | 0.183 | **0.247** | 0.240 | Weighted |
| 60-100 | 0.352 | 0.351 | **0.378** | EFA-A |
| 100-200 | **0.523** | 0.400 | 0.489 | Standard |

Overall: Std=0.245, Wt=0.268, EFA-A=0.280. EFA-A wins only at 60-100 items.
**Rejected:** loading products unreliable because parallel analysis underestimates nF (~69%
of true), creating too-large factors that mix constructs.

**Option D: EFA within-factor pairs + observed |r| weights**

Comprehensive simulation (360 conditions): 8 nF (4-30) × 3 ipf (3-10) × 3 pct (5-25%)
× 5 reps, N=300, iterations=50.

| Range | Standard | Weighted | EFA-D | Winner |
|-------|----------|----------|-------|--------|
| <30 items | 0.108 | **0.160** | 0.158 | Wt (by 0.002) |
| 30-60 | 0.198 | 0.244 | **0.261** | EFA-D |
| 60-100 | 0.306 | 0.307 | **0.351** | EFA-D |
| 100-200 | 0.501 | 0.405 | **0.529** | EFA-D |
| >200 | **0.656** | 0.393 | 0.653 | Std (by 0.003) |

Overall on simulated data: EFA-D=0.347 vs Std=0.303 vs Wt=0.286. EFA-D wins 18/24 nF×ipf
cells (75%) on simulated data.

**BUT: external validation on real data showed EFA-D loses badly:**

| Dataset | Items | nF | EFA-D AUC | Coupled AUC | EFA-D MCC₁.₅ | Coupled MCC₁.₅ |
|---------|-------|----|-----------|-------------|-------------|----------------|
| Schroeders | 60 | ~10 | 0.583 | **0.613** | 0.126 | **0.182** |
| Schneider | 31 | ~5 | 0.459 | **0.729** | -0.019 | **0.114** |
| Niessen | 100 | ~5 | 0.587 | **0.637** | 0.037 | **0.073** |
| Goldammer S1 | 60 | ~8 | 0.605 | **0.717** | 0.167 | **0.314** |
| Goldammer S2 | 60 | ~6 | 0.599 | **0.656** | 0.096 | **0.286** |
| Goldammer S3 | 60 | ~7 | **0.532** | 0.456 | 0.096 | -0.032 |

Means: EFA-D AUC=0.561 vs Coupled AUC=0.635. Coupled wins 5/6 datasets.

**Why EFA-D fails on real data:** Parallel analysis on real data with cross-loadings, method
effects, and noisy structure produces poor factor assignments. Simulated data has clean
factor structure that EFA recovers well; real data doesn't. Additionally, datasets with high
careless rates (Goldammer 63-67%) corrupt the correlation matrix, degrading EFA quality.
The coupled method's empirical top-k% by |r| is more robust because the ranking of |r|
survives contamination better than the factor structure does.

**Decision: BOTH REJECTED.** The coupled/weighted 2-level switch remains the default.
The simulation-to-real gap is a cautionary tale: method improvements must be validated on
real data, not just simulations.

Reports in `archive/efa_comparison/` (Option A) and `archive/efa_d_comparison/` (Option D).

### 2026-03-30 — Weighted Mode for Short Questionnaires

**Problem:** Standard ReReReRe selects top-k% item pairs by |r|. With short questionnaires
(≤60 items), few high-|r| pairs exist, leading to weak detection (MCC ~0.03-0.11).

**Solution:** New `rowCor_weighted()` function computes a weighted coherence score using ALL
item pairs, where each pair's contribution is proportional to its sample-level |r|. Strong
correlations count more, weak pairs contribute proportionally less but aren't discarded.

**Algorithm (weighted mode):**
1. Standardize each item across respondents (z-scores)
2. For each pair: compute cross-product z_A × z_B per respondent
3. Weight each pair by its |r_sample|
4. Weighted average across all pairs → "coherence score" per respondent
5. Permutation baseline: same k random pairs, same weights, z-score as usual

**Simulation results (18 conditions: nF=4-20, ipf=3-10, N=300, 10% careless, 3 reps):**

| Item range | Standard (oracle MCC) | Weighted (oracle MCC) | Difference |
|-----------|----------------------|----------------------|------------|
| ≤30 items | 0.110 | **0.151** | **+37%** |
| 30-60 items | 0.193 | **0.255** | **+32%** |
| 60-100 items | **0.315** | 0.300 | −5% |
| 100-200 items | **0.524** | 0.414 | −21% |

**Crossover at ~60 items.** Weighted wins below, coupled wins above.

**Integration into ReReReRe.R:**

| Parameter | Default | Options | Effect |
|-----------|---------|---------|--------|
| `mode` | `"auto"` | `"auto"`, `"coupled"`, `"weighted"` | Auto: weighted if ≤60 items, coupled if >60 |
| `min_r` | 0.0 | 0.0-1.0 | Min |r| to include pair in weighted mode (0 = all) |

New output column: `mode_used` ("coupled" or "weighted").

**Verification test (single seed):**
- 24 items: weighted MCC=0.118 vs coupled MCC=0.031 (auto chose weighted ✓)
- 120 items: coupled MCC=0.442 vs weighted MCC=0.197 (auto chose coupled ✓)

Also moved `Test_Weighted_ReReReRe.R` to `archive/scripts/`.

### 2026-03-30b — EFA-Based Pair Selection: Tested and Rejected (Option A)

Tested EFA-guided pair selection: parallel analysis → EFA (oblimin, minres) → within-factor
pairs only → weights = |λ_i × λ_j| (loading products as shrinkage estimator).

**54 conditions** (nF={4,6,8,10,15,20} × ipf={3,6,10} × 3 reps), N=300, 10% careless.

| Item range | Standard | Weighted | EFA | Winner |
|-----------|----------|----------|-----|--------|
| <30 | 0.101 | **0.167** | 0.152 | Weighted |
| 30-60 | 0.183 | **0.247** | 0.240 | Weighted |
| 60-100 | 0.352 | 0.351 | **0.378** | EFA |
| 100-200 | **0.523** | 0.400 | 0.489 | Standard |

Overall: Standard=0.245, Weighted=0.268, EFA=0.280.

**EFA wins only at 60-100 items** (+0.027 MCC over standard). The loading-product weights
provide a shrinkage advantage (model-implied r is less noisy than observed r). But parallel
analysis systematically underestimates nF (~70% of true), creating too-large factors that
mix items from different constructs. This hurts at <60 items (EFA < Weighted) and provides
no advantage at >100 items (Standard already selects within-factor pairs).

**Decision: rejected.** The gain at 60-100 items doesn't justify adding EFA as a dependency
(model estimation, Heywood case handling, parallel analysis instability). The 2-level
auto-switch (weighted ≤60, coupled >60) captures most of the benefit.

Report with 12 plots saved in `archive/efa_comparison/`. Scripts moved to `archive/scripts/`.

### 2026-03-28b — Comprehensive External Validation (6 datasets + Johnson inject-and-detect)

Ran ReReReRe (corProp=0.03, z=1.5 + auto_z) vs practical Mahalanobis (chi-sq .001) on 6 real
datasets with ground truth, plus inject-and-detect on Johnson IPIP-NEO-300.

**Results summary (corProp=0.03):**

| Dataset | N | Items | nF | RR AUC | Mah AUC | RR MCC z=1.5 | RR auto | Mah MCC |
|---------|---|-------|----|--------|---------|-------------|---------|---------|
| Schroeders 2022 | 605 | 60 | ~10 | **0.606** | 0.537 | **0.177** | **0.228** | 0.045 |
| Schneider QoL | 1649 | 31 | ~5 | **0.735** | 0.724 | 0.118 | 0.170 | **0.208** |
| Niessen 2016 | 180 | 100 | ~6 | **0.639** | 0.436 | 0.073 | -0.067 | 0.000 |
| Goldammer S1 | 291 | 60 | ~8 | 0.711 | **0.787** | **0.320** | 0.269 | 0.257 |
| Goldammer S2 | 265 | 60 | ~6 | 0.674 | **0.795** | **0.304** | 0.241 | 0.266 |
| Goldammer S3 | 523 | 60 | ~7 | 0.462 | 0.550 | -0.036 | -0.072 | 0.011 |

**Key findings:**
1. RR AUC > Mah AUC on 4/6 datasets (all except Goldammer S1-S2 which have 64% careless rate)
2. RR MCC (z=1.5) > Mah MCC on 4/6 datasets — practical Mahalanobis remains weak
3. Auto-z beat z=1.5 on Schroeders (0.228 vs 0.177) — auto chose z=0.72, appropriate for 60 items
4. Goldammer S3 (longitudinal) fails for both methods — 67% careless, t1 data only
5. Niessen: both methods essentially at chance (nF=5, speed manipulation ≠ inconsistent carelessness)
6. Johnson (inject-and-detect, 300 items): mean oracle MCC=0.726-0.750, auto-z MCC=0.726, spec=1.000

**Goldammer S1 breakdown:** 100% careless AUC=0.692, 33% careless AUC=0.738. Partial careless
slightly easier for RR; full random easier for Mahalanobis.

### 2026-03-28 — ReReReRe.R updated with multiverse-informed defaults

Applied all findings from full multiverse and 2D calibration to the algorithm:

1. **corProp default 0.05 → 0.03**: fewer but more selective coupled pairs (MCC +0.025)
2. **auto_z engine replaced**: 1D nF-only LOESS → 2D total_items-based LOESS (R² doubled)
3. **60-point calibration table** embedded (10 nF × 6 ipf, 30 reps each)
4. **Linear fallback formula** added: z = 0.505 + 0.0042 × total_items
5. **z clamped to [0.3, 3.5]** to prevent extreme thresholds

Tested on simulated data: 24 items → z=1.04, 180 items → z=1.00, MCC=0.450 on 30-factor test.
Calibration curve captures the U-shape correctly.

### 2026-03-26d — Full Multiverse with items_per_factor varied

Ran `Multiverse_Full.R`: 7 nF × 3 ipf × 3 n × 3 pct × 3 corProp × R=10 = 5,670 RR calls,
85,050 result rows. Runtime: 161 minutes.

**Three key findings:**

1. **items_per_factor is a major driver** — at nF=12: ipf=3→MCC 0.14, ipf=6→0.30, ipf=10→0.46.
   Total items matters more than nF alone. This means the practical recommendation should be
   framed in terms of total questionnaire length, not just number of constructs.

2. **corProp=0.03 beats 0.05** (MCC 0.352 vs 0.327) — fewer but more selective pairs work better.
   Updates previous recommendation.

3. **Auto-z ≈ fixed z=1.5** (MCC 0.268 vs 0.274) — the LOESS-calibrated threshold provides no
   meaningful improvement over a universal z=1.5 default. Auto-z adapts correctly (low z at low
   nF, high z at high nF) but the aggregate performance is identical. The feature adds
   convenience without performance cost or gain.

### 2026-03-26c — Auto-calibration: nF → z_threshold via parallel analysis

**Added `auto_z` parameter** to `ReReReRe()`. When `auto_z=TRUE`, the function:
1. Runs `psych::fa.parallel()` to estimate nF from the data
2. Looks up the optimal z_threshold via a LOESS-smoothed calibration curve (nF=2-40)
3. Uses the calibrated threshold for flagging

**LOESS curve** (span=0.4, fitted on 39-point calibration means):
- nF=5 → z=1.01, nF=10 → z=0.76, nF=15 → z=0.60 (minimum),
  nF=20 → z=0.69, nF=25 → z=0.83, nF=30 → z=1.25, nF=40 → z=2.10
- Clamped to nF=[2,40] for out-of-range inputs

**Test result** (10-factor simulated data, n=200): parallel analysis detected 7 factors
(typical underestimate with weak loadings), calibrated z=0.95, sensitivity=0.68,
specificity=0.78. Consistent with expected performance at nF~7-10.

**Output changes:** `flagged` column now uses z_score threshold (was percentile-based).
Added `z_threshold_used` and `nFactors_detected` columns to output dataframe.

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

**PISA 2018 re-test with reverse coding:** The sign-flip inversion (AUC=0.35) was partly an
artefact. Reverse coding gives AUC=0.612 (raw) and 0.678 (0-1 rescaled), both in the correct
direction. Mahalanobis AUC=0.774 still wins — PISA has 36 short scales (2-9 items) and the
dominant careless pattern (acquiescence/speed) is consistent, not inconsistent. All 135 item
ranges verified correct (computer-based, no data entry errors). Scale mix: 97x4pt, 27x5pt,
5x6pt. Rescaling to [0,1] helps because it removes the bias from unequal Likert ranges in
the individual-level correlation.

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