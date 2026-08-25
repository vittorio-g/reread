# reread 0.1.0

* First release.
* `reread()` one-call careless-response detection with automatic, label-free
  flagging and a calibrated prevalence estimate.
* The `rc` index (`rc_index()`), the profile-gated logistic ensemble
  (`ensemble_score()`), and the label-free calibration (`auto_flag()`): a
  one-sided fit of the attentive mode, cut at 2.5 sigma, applied unconditionally.
* `reread(score = "person_total")` applies the identical threshold to the
  person-total correlation alone, for analysts who prefer that index.
* Auxiliary indices: `longstring()`, `person_total()`, `irv()`, `mahalanobis_d2()`.
* Simulation utilities `simulate_clean()` and `inject_careless()`, and
  `benchmark_indices()` for head-to-head comparison with DeLong tests.
* Operational permutation count defaults to 400; per-respondent Monte Carlo standard error
  of the score and a `borderline` flag are reported. Bundled `demo_careless` dataset;
  faithful port of the shipped tool at re-re.re.
