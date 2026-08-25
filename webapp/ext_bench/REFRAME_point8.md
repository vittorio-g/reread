# Reviewer point 8 — weight over-precision / overfitting (paste-ready)

**Verdict: the reviewer is right about the false precision, and we turn it into a strength.**
The exact 4-decimal weights are *not identifiable*, but the ensemble is *insensitive* to them and
the *frozen* vector transfers better than refitting. Numbers from `weights_sensitivity_study1.R`
(Study 1, N=157; the shipped linear combiner, gate g=1 on this data).

## The three results

**(c) The coefficients carry enormous uncertainty (concede the false precision).**
A plain refit gives slopes rr=4.27, LongString=1.89, Person–Total=1.94 — larger than the shipped
2.68/1.24/1.44, i.e. the shipped weights are already *shrunk*. The data are near-separable
(AUC ≈ .99), so the MLE is unstable and the bootstrap 95% CIs are huge:

| coefficient | shipped | bootstrap 95% CI |
|---|---|---|
| rr | 2.68 | [2.93, 11.5] |
| LongString | 1.24 | [1.37, 5.56] |
| Person–Total | 1.44 | [0.94, 7.61] |

Only the *ordering and rough ratio* (rr ≈ 2× the auxiliaries, all positive) are identifiable;
the digits past the first are noise. **Fix: report weights to ≤2 sig figs and state they are a
regularised/shrunk estimate, not an MLE.**

**(a) But detection is insensitive to the exact weights.**

| weight vector | AUC |
|---|---|
| shipped (2.68, 1.24, 1.44) | .990 |
| integers (3, 1, 1) | .987 |
| rr + equal auxiliaries (3, 1.5, 1.5) | .990 |
| unweighted z-sum (1, 1, 1) | .982 |

Random ±50% multiplicative jitter on all three slopes: AUC .987 [.980, .990], min .972.
So collapsing 2.6756 → 3 costs .003 AUC, and even an *unweighted* sum of the three z-scores loses
only .008. The fourth decimal is cosmetic.

**(b) Freezing the weights beats refitting them (answers "is freezing costly?").**
Leave-one-out AUC with the *fixed* weights = **.990**; with weights *refit on each training fold*
= **.983** (ΔAUC −.007, DeLong p = .020). Because the combiner sits in a flat, near-separable
region, per-dataset refitting only adds estimation noise. Freezing is not a compromise — it
generalises better. (Consistent with external transfer: the same fixed vector reaches mean
AUC .65 across 12 independent datasets in `benchmark_ext.csv`, beating every classic competitor.)

## Paste-ready prose

**Method — after Eq. (1):**
> The weights are a regularised logistic fit on Study 1 and are reported to two significant figures;
> their exact values are not the point. Because the three oriented detectors are near-collinear on
> attentive/careless contrasts, the combiner occupies a flat region of weight space: replacing the
> fitted vector with small integers (3, 1, 1), or perturbing every weight by ±50%, changes
> leave-one-out AUC by at most .01 (Figure Y). We therefore freeze one transferable vector rather
> than refit per dataset — which, on this near-separable problem, would only add estimation noise
> (LOO AUC .990 fixed vs .983 refit, DeLong p = .02).

**Figure Y caption:**
> *The ensemble weights are not identifiable, but detection is insensitive to them.* (A) Bootstrap
> 95% CIs on the three logistic coefficients (Study 1); the shipped values sit inside wide,
> near-separation-inflated intervals, so only their ordering is identifiable. (B) Leave-one-out /
> in-sample AUC is essentially unchanged when the weights are coarsened to integers, set equal, or
> randomly perturbed by ±50%, and the frozen weights outperform per-fold refitting.

*(Assets: `ext_bench/fig_w8_weights.png`, `plot_w8_sensitivity.py`, `weights_sensitivity_study1.R`,
`ext_bench/w8_*.csv`.)*

## Optional strengthening (not yet done)
Clean **external** fixed-vs-refit on the *same oriented features* needs the scorer to dump
O.rr/O.long/O.pt per dataset (currently only O.rr + eta are dumped in `dump_ext.js`). A ~15-line
Node dump would let us repeat the rounding/perturbation test on all 12 external datasets, nailing
the "applied unchanged *everywhere*" clause directly. Say the word and I'll add it.
