# Reviewer point 3 — reframing the contribution (paste-ready prose)

**What the data force us to claim, and nothing more.** All numbers below are from the
*current, reproducible* pipeline (`benchmark_ext.csv` + `rr_incremental_current.csv`,
same frozen 4000-subsample; eta/rr/LongString AUCs reproduce `benchmark_ext.csv` exactly).

---

## 0. The one-paragraph honest position

`rr` and Person–Total are near-equivalent for **ranking** respondents (AUC); the permutation
null does **not** buy ranking power over the simplest profile-typicality index. What it buys
is a **self-calibrating scale** (a per-respondent z against the person's own chance level, so
no external cut-off and comparability across questionnaires) and a **modest, conditional gain
at the classification boundary** (oracle-MCC), visible on long batteries whose ground truth
reflects *profile incoherence* (smarvus ΔoMCC = +.074, 95% CI [+.032,+.100]; Kay S2 +.080
[+.002,+.163]), matching the simulation ablation (Fig 5, +.07–.10 at ≥60 items). Where the
ground truth reflects *speed* or *consistent* careless, or the questionnaire is too short to
build the pair structure, `rr` adds nothing — by design. The ensemble, not `rr` alone, is what
beats established indices: mean external AUC .652 vs PsychSyn .615, rr .573.

---

## 1. Abstract — replace the two overclaims

**DELETE** any phrasing like *"self-calibrating z-score that requires no distributional
assumption or external cut-off"* applied to the whole detector, and *"convergent validity that
scaled with questionnaire length (up to AUC .92)"* as the headline external result.

**REPLACE WITH** (drop-in):

> …yielding a per-respondent *z*-score calibrated against the respondent's own random-pairing
> baseline, so scores are comparable across questionnaires without an external cut-off. Because
> `rr` targets *inconsistent* incoherence, it is combined with two complementary detectors
> (LongString, Person–Total) in a small, reliability-gated logistic ensemble. On experimentally
> induced careless responding the ensemble ranked careless respondents with leave-one-out
> AUC = .985. Against fourteen independent datasets with pattern-independent ground truth, the
> ensemble outperformed established careless indices (Mahalanobis, IRV, LongString, psychometric
> synonyms) — mean AUC .65 vs .62 for the best competitor — while remaining transparent. The
> permutation index contributes over and above the auxiliary detectors specifically on long
> batteries whose ground truth reflects profile incoherence (incremental oracle-MCC up to +.08,
> bootstrap CI excluding zero), and is deliberately inert against speed-based or single-scale
> carelessness — a discriminant property, not a failure.

Note: keep "self-calibrating" only for the *scale/cut-off* property, never for the flag
threshold (which uses a 2-component Gaussian mixture — that IS a distributional assumption).

---

## 2. "The Present Research" — rewrite contribution (c)/(d)

**REPLACE** the external-validity bullet with:

> (c) external validity on 14 independent datasets with pattern-independent ground truth,
> benchmarking the ensemble **head-to-head against four established indices** and isolating the
> permutation index's **incremental** contribution over the auxiliary detectors; (d) a controlled
> simulation localising where `rr` helps. Throughout we distinguish two questions the field
> conflates: **ranking** respondents (where a simple profile index already does well) and
> **flagging** them at a calibrated boundary (where the per-respondent null adds value).

---

## 3. Study 2 — NEW subsection to insert (with the new figure)

**Heading:** *Isolating the permutation index's contribution.*

> Convergent AUC shows the ensemble is competitive, but not whether the permutation index earns
> its place beside the cheaper auxiliaries. We therefore repeated, on real data, the ablation run
> in simulation: within each dataset we fit two logistic detectors under repeated stratified
> cross-validation — auxiliaries only (LongString + Person–Total) versus auxiliaries + `rr` — and
> bootstrapped the difference in oracle-threshold MCC (Figure X). Two facts emerge. First, for
> **ranking**, `rr` and Person–Total are near-identical (e.g. smarvus AUC .834 vs .833), so the
> index does not reorder the sample. Second, at the **classification boundary** `rr` adds a
> modest but reliable gain **exactly where the mechanism predicts** — long, multi-factor batteries
> whose ground truth reflects profile incoherence: smarvus (141 items, instructed-check GT)
> ΔMCC = +.074 [+.032, +.100] and Kay S2 (363 items, self-report exclusion) +.080 [+.002, +.163],
> both bootstrap CIs excluding zero, echoing the simulation (Fig 5). The contribution is null
> where the ground truth is speed-based or the battery too short to construct the pair set — a
> boundary condition, not a failure. (The lone short-scale exception, moss23 at 8 items, is a
> spurious inversion — `rr` there is anti-predictive, AUC .29 — and is flagged, not counted as
> support.)

**Figure X caption:**

> *Incremental contribution of the `rr` index on external data.* Difference in oracle-MCC between
> the full ensemble (LongString + Person–Total + `rr`) and the auxiliaries alone, per dataset,
> with 95% stratified-bootstrap CIs; datasets split into the method's envelope (≥60 multi-item)
> and boundary conditions. `rr` adds reliably only on long, incoherence-relevant batteries
> (smarvus); it is inert elsewhere by design. moss23 (8 items) is a spurious inversion, not
> mechanism evidence.

*(Figure file: `ext_bench/fig_rr_incremental.png`; table: `ext_bench/rr_incremental_current.csv`.)*

---

## 4. Discussion — three edits

**4a. "Relation to existing methods"** — add head-to-head + honest ceiling:

> Benchmarked directly against Mahalanobis distance, IRV, LongString and psychometric synonyms on
> 14 datasets, the ensemble had the highest mean AUC (.65), ahead of psychometric synonyms (.62);
> the permutation index alone was mid-pack (.57). This is the honest ordering: `rr`'s value is not
> as a standalone detector but as a self-calibrating component that sharpens the ensemble's
> decision boundary on the questionnaires it was designed for.

**4b. "A ranking problem before a thresholding problem"** — fold in the AUC/MCC split explicitly:

> A recurring theme is that Person–Total and `rr` rank respondents almost identically, yet differ
> at the flag boundary: the per-respondent null calibrates each score against the respondent's own
> chance level, which is what a *threshold* — not a ranking — needs. This is why `rr`'s external
> contribution surfaces in MCC, not AUC.

**4c. Limitations** — state the envelope plainly:

> `rr`'s incremental value is confined to long, multi-factor batteries whose careless responding is
> inconsistent; on short or single-scale questionnaires, or where carelessness is fast-but-
> consistent, the auxiliaries and design-time checks carry the load and `rr` should be expected to
> add little. The near-ceiling simulation values (AUC → 1.0) reflect a clean corruption measure;
> real external validities (.5–.9) are lower because ground truth is a proxy and real careless is
> partial.

---

## 5. Numbers to standardise across the paper (pick these, kill the stale ones)

| Claim | Value to use | Source |
|---|---|---|
| Study 1 ranking | LOO AUC = .985 | current |
| Ensemble vs best competitor (ext) | mean AUC .652 vs PsychSyn .615 | benchmark_ext.csv |
| rr alone (ext) | mean AUC .573 | benchmark_ext.csv |
| rr incremental (flagship) | smarvus +.074 [.032,.100]; Kay S2 +.080 [.002,.163] | rr_incremental_current.csv |
| rr incremental (simulation) | +.07–.10 at ≥60 items (Fig 5) | sim |
| Discriminant (speed/short) | ΔMCC ≈ 0 (warning, krause, duckworth) | rr_incremental_current.csv |

**⚠️ Reproducibility red flag (fix before submission):** the PDF I reviewed had a Table 3 with
opsy_16PF ens/rr/PT = .828/.793/.826, but the current reproducible pipeline gives .682/.702/—
(uniform gap across all detectors → different scoring generation). Table 3 in the live paper must
be regenerated from ONE pipeline so a reviewer can reproduce it. The current `benchmark_ext.csv`
is that pipeline; make the paper table read from it.
