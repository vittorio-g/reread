# reread

**Screening careless respondents with a self-calibrating threshold**

`reread` screens questionnaire data for careless or insufficient-effort
responding, and — unlike index-only tools — decides where to cut without labels.
The threshold is anchored on the attentive mode of the score distribution and
reported with a calibrated prevalence estimate and with what the same rule flags
on careless-free data of that size. The score it cuts is **`rc`** (*relative
coherence*), a permutation-based index of individual response coherence that
needs no declared factor structure,
combined with two complementary detectors (longstring and person–total
correlation) in a compact, reliability-gated logistic ensemble. A BIC-gated
BIC-gated, attentive-anchored rule sets the flagging threshold **without labels** and returns a
calibrated prevalence estimate (a one-sided fit of the attentive mode that high
careless rates cannot contaminate).

The package is a faithful R port of the shipped browser tool at
[ca.re-re.re](https://ca.re-re.re) and reproduces its scores.

## Installation

```r
# from a local build until it is on CRAN:
# install.packages("reread_0.1.0.tar.gz", repos = NULL, type = "source")
install.packages("reread")   # once released on CRAN
```

## Quick start

```r
library(reread)

# a data frame / matrix of Likert responses (respondents x items)
fit <- reread(demo_careless$responses)
fit
#> the reread procedure careless-response detection
#>   Respondents: 400   Items: 90
#>   Flagged: ...   Sensitivity: standard
#>   Estimated careless prevalence: ...

# per-respondent scores and flags
head(fit$scores)

# who was flagged
which(fit$flagged)

# score distribution and threshold (needs ggplot2)
plot(fit)
```

There is a single shipped cut, and it is the better choice up to a careless rate
of about 35%. Only past that does a more liberal cut pay off, and then you set it
explicitly:

```r
fl <- auto_flag(ensemble_score(my_data)$eta, z = 1.5)
```

## What it computes

* `reread()` — one-call detection: scores, flags, and a calibrated prevalence.
* `rc_index()` — the permutation coherence index alone.
* `ensemble_score()` — the ensemble log-odds.
* `auto_flag()` — the label-free calibration (one-sided attentive-mode fit, cut at 2.5σ). It is equivariant to the scale of whatever score it is given, so it can be applied to any index.
* `longstring()`, `person_total()`, `irv()`, `mahalanobis_d2()` — auxiliary indices.
* `simulate_clean()` / `inject_careless()` — generate test data.
* `benchmark_indices()` — compare against other indices with DeLong tests.

### Cutting a different score

On external data the ensemble and the person--total correlation are close to
interchangeable, in ranking and in decision quality alike, so the choice between
them is offered rather than hidden:

```r
reread(responses, score = "person_total")   # same cut, single index
```

The ensemble stays the default. It is the more robust of the two — it separates
far better on collected data with clean ground truth, and it removes a steadier
share of the sample across datasets — and the careless probability and the
borderline band are defined only for it (both are `NA` under `"person_total"`).

## When to use it

`reread` targets **inconsistent** careless responding (random and mixed
answering) on **multi-construct** questionnaires; it works best on batteries of
roughly 60 or more items. It does
not detect content-responsive distortions such as faking, and on short or
single-factor scales the auxiliary detectors carry more of the signal.

## Citation

Guerrieri, V., Gallucci, M., & Passarelli, M. (2026). *the reread procedure: A permutation-based
coherence index and a compact ensemble for detecting careless responding in
questionnaire data.*

## License

GPL-3.
