/* gate_tailtest.js — SOLUTION 1: replace the mixture structure-gate with an excess-tail test.
 *
 * Rationale: since the cut is now placed one-sided (m + z*sigma from the left side only), the
 * EM mixture is left with a single job - deciding whether a careless group exists at all. That
 * question is answered more directly, and with one parameter instead of three hand-set thresholds,
 * by asking whether the upper tail is heavier than a clean sample's would be:
 *
 *    H0 (no careless): P(eta > m + z*sigma) = 1 - Phi(z) = .0062 for z = 2.5
 *    observed k = #{eta > m + z*sigma};  gate opens iff  P(X >= k | n, p0) < alpha
 *
 * The test reuses m and sigma from the SAME one-sided fit that places the cut, so the gate and the
 * threshold are finally the same object. Compared here against the shipped mixture gate.
 */
const R = require("./site/rerere.js"); const W = R.WEIGHTS; const fs = require("fs");
const Z = 2.5, P0 = 0.006209665;                       // 1 - Phi(2.5)

/* upper-tail binomial p-value, computed in log space */
function binomUpper(k, n, p) {
  if (k <= 0) return 1;
  if (k > n) return 0;
  const lg = x => { // log-gamma (Lanczos)
    const g = [676.5203681218851,-1259.1392167224028,771.32342877765313,-176.61502916214059,
               12.507343278686905,-0.13857109526572012,9.9843695780195716e-6,1.5056327351493116e-7];
    if (x < 0.5) return Math.log(Math.PI/Math.sin(Math.PI*x)) - lg(1-x);
    x -= 1; let a = 0.99999999999980993, t = x + 7.5;
    for (let i = 0; i < 8; i++) a += g[i]/(x+i+1);
    return 0.5*Math.log(2*Math.PI) + (x+0.5)*Math.log(t) - t + Math.log(a);
  };
  const lchoose = (n,i) => lg(n+1) - lg(i+1) - lg(n-i+1);
  let sum = 0, maxTerm = -Infinity;
  const terms = [];
  for (let i = k; i <= n; i++) {
    const lt = lchoose(n,i) + i*Math.log(p) + (n-i)*Math.log1p(-p);
    terms.push(lt); if (lt > maxTerm) maxTerm = lt;
    if (i > k + 5 && lt < maxTerm - 50) break;          // negligible tail
  }
  for (const lt of terms) sum += Math.exp(lt - maxTerm);
  return Math.min(1, Math.exp(maxTerm) * sum);
}
function tailGate(eta, alpha) {
  const lf = R._leftFit ? R._leftFit(eta) : null;
  let m, sigma;
  if (lf) { m = lf.m; sigma = lf.sigma; }
  else {                                                // fall back: recompute the one-sided fit
    const af = R.autoFlag(eta, "standard"); m = af.attentiveMode; sigma = af.attentiveSigma;
  }
  const tau = m + Z * sigma;
  const k = eta.reduce((s, v) => s + (v > tau ? 1 : 0), 0);
  const p = binomUpper(k, eta.length, P0);
  return { open: p < alpha, k, expected: P0 * eta.length, p, tau, m, sigma,
           flaggedFrac: k / eta.length };
}
module.exports = { tailGate, binomUpper, Z, P0 };
