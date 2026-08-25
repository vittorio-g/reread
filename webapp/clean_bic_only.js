/* clean_bic_only.js — is there a gate variant that actually protects the tail?
 *
 * The A/B showed that the shipped gate (BIC or excess-tail) lowers the AVERAGE
 * false-flag rate on clean data but leaves the worst case untouched: the 95th
 * percentile is identical to the gateless engine in 42 of 48 cells. On a narrow
 * corpus the BIC branch alone had given a perfect 0.00%, so the question is
 * whether that survives the wider grid -- if it does, a genuinely protective
 * option exists and the choice is a real trade rather than a free removal.
 *
 * eta is computed once per clean dataset and the three decision rules are applied
 * to it, so the comparison is exact rather than paired-by-seed.
 *
 *   node clean_bic_only.js   -> clean_bic_only.csv
 */
const fs = require("fs");
const R = require("./site/reread.js");
const Z = 2.5, P_TAIL = 0.0062, TAIL_ALPHA = 1e-6;
function lgamma(x) {
  const c = [76.18009172947146, -86.50532032941677, 24.01409824083091, -1.231739572450155, 0.1208650973866179e-2, -0.5395239384953e-5];
  let y = x, t = x + 5.5; t -= (x + 0.5) * Math.log(t);
  let s = 1.000000000190015; for (let j = 0; j < 6; j++) s += c[j] / ++y;
  return -t + Math.log(2.5066282746310005 * s / x);
}
function binomUpperP(k, n, p) {
  if (k <= 0) return 1; if (k > n) return 0;
  const lch = (n, i) => lgamma(n + 1) - lgamma(i + 1) - lgamma(n - i + 1);
  const T = []; let mx = -Infinity;
  for (let i = k; i <= n; i++) { const lt = lch(n, i) + i * Math.log(p) + (n - i) * Math.log1p(-p); T.push(lt); if (lt > mx) mx = lt; if (i > k + 5 && lt < mx - 50) break; }
  let s = 0; for (const t of T) s += Math.exp(t - mx); return Math.exp(mx) * s;
}
const rng = R._rng(4242);
function gauss() { let u = 0, v = 0; while (u === 0) u = rng(); while (v === 0) v = rng(); return Math.sqrt(-2 * Math.log(u)) * Math.cos(2 * Math.PI * v); }
function qnorm(p) {
  const a = [-3.969683028665376e+01, 2.209460984245205e+02, -2.759285104469687e+02, 1.383577518672690e+02, -3.066479806614716e+01, 2.506628277459239e+00];
  const b = [-5.447609879822406e+01, 1.615858368580409e+02, -1.556989798598866e+02, 6.680131188771972e+01, -1.328068155288572e+01];
  const c = [-7.784894002430293e-03, -3.223964580411365e-01, -2.400758277161838e+00, -2.549732539343734e+00, 4.374664141464968e+00, 2.938163982698783e+00];
  const d = [7.784695709041462e-03, 3.224671290700398e-01, 2.445134137142996e+00, 3.754408661907416e+00];
  const pl = 0.02425; let q, r;
  if (p < pl) { q = Math.sqrt(-2 * Math.log(p)); return (((((c[0] * q + c[1]) * q + c[2]) * q + c[3]) * q + c[4]) * q + c[5]) / ((((d[0] * q + d[1]) * q + d[2]) * q + d[3]) * q + 1); }
  if (p <= 1 - pl) { q = p - 0.5; r = q * q; return (((((a[0] * r + a[1]) * r + a[2]) * r + a[3]) * r + a[4]) * r + a[5]) * q / (((((b[0] * r + b[1]) * r + b[2]) * r + b[3]) * r + b[4]) * r + 1); }
  q = Math.sqrt(-2 * Math.log(1 - p)); return -(((((c[0] * q + c[1]) * q + c[2]) * q + c[3]) * q + c[4]) * q + c[5]) / ((((d[0] * q + d[1]) * q + d[2]) * q + d[3]) * q + 1);
}
const VARIANTS = [
  { loadLo: .4, loadHi: .8, rho: .30, rev: .40, K: 5 }, { loadLo: .2, loadHi: .5, rho: .30, rev: .40, K: 5 },
  { loadLo: .4, loadHi: .8, rho: .60, rev: .40, K: 5 }, { loadLo: .4, loadHi: .8, rho: .30, rev: .00, K: 5 },
  { loadLo: .4, loadHi: .8, rho: .30, rev: .40, K: 7 }
];
function genClean(J, n, v) {
  const ipf = 6, F = Math.max(2, Math.round(J / ipf));
  const load = [], fac = [], tau = [], THR = [];
  for (let c = 1; c < v.K; c++) THR.push(qnorm(c / v.K));
  for (let j = 0; j < J; j++) { let l = v.loadLo + (v.loadHi - v.loadLo) * rng(); if (rng() < v.rev) l = -l; load.push(l); fac.push(j % F); tau.push(0.45 * gauss()); }
  const sq = Math.sqrt(v.rho), sq1 = Math.sqrt(1 - v.rho), M = new Float64Array(n * J);
  for (let i = 0; i < n; i++) {
    const common = gauss(), f = [];
    for (let k = 0; k < F; k++) f.push(sq * common + sq1 * gauss());
    for (let j = 0; j < J; j++) {
      const y = tau[j] + load[j] * f[fac[j]] + Math.sqrt(1 - load[j] * load[j]) * gauss();
      let cat = 1; for (let c = 0; c < v.K - 1; c++) if (y > THR[c]) cat = c + 2;
      M[i * J + j] = cat;
    }
  }
  return M;
}
const JS = [30, 48, 60, 90, 120, 180, 240, 300], NS = [100, 150, 200, 300, 500, 1000], REPS = 30;
const mn = a => a.reduce((s, v) => s + v, 0) / a.length;
const q95 = a => { const s = [...a].sort((x, y) => x - y); return s[Math.floor(0.95 * s.length)]; };
const all = { none: [], both: [], bic: [] };
const out = [["J", "n", "none_mean", "none_p95", "both_mean", "both_p95", "bic_mean", "bic_p95", "bic_open_frac"]];
console.log("clean data: no gate  |  BIC or tail (shipped)  |  BIC only          BIC opens");
console.log("  J    n     mean   p95        mean   p95         mean   p95         on");
for (const J of JS) for (const n of NS) {
  const A = { none: [], both: [], bic: [] }; let bicOpen = 0;
  for (let rep = 0; rep < REPS; rep++) {
    const M = genClean(J, n, VARIANTS[rep % VARIANTS.length]);
    const e = R.ensemble(M, n, J, { seed: 1, iterations: 150, skipD2: true });
    const eta = e.eta, lf = R.leftFit(eta), cut = lf.m + Z * lf.sigma;
    let k = 0; for (const v of eta) if (v > cut) k++;
    const af = R.autoFlag(eta, "standard");
    const tail = binomUpperP(k, n, P_TAIL) < TAIL_ALPHA;
    if (af.bicGate) bicOpen++;
    A.none.push(k / n);
    A.both.push((af.bicGate || tail) ? k / n : 0);
    A.bic.push(af.bicGate ? k / n : 0);
  }
  for (const key of ["none", "both", "bic"]) all[key].push(...A[key]);
  console.log(`${String(J).padStart(3)} ${String(n).padStart(5)}  ${(100 * mn(A.none)).toFixed(2).padStart(5)}% ${(100 * q95(A.none)).toFixed(1).padStart(5)}%    ` +
    `${(100 * mn(A.both)).toFixed(2).padStart(5)}% ${(100 * q95(A.both)).toFixed(1).padStart(5)}%     ` +
    `${(100 * mn(A.bic)).toFixed(2).padStart(5)}% ${(100 * q95(A.bic)).toFixed(1).padStart(5)}%     ${bicOpen}/${REPS}`);
  out.push([J, n, mn(A.none).toFixed(5), q95(A.none).toFixed(5), mn(A.both).toFixed(5), q95(A.both).toFixed(5),
    mn(A.bic).toFixed(5), q95(A.bic).toFixed(5), (bicOpen / REPS).toFixed(3)]);
}
console.log("\noverall (1440 clean datasets):");
for (const key of ["none", "both", "bic"])
  console.log(`  ${key.padEnd(5)} mean ${(100 * mn(all[key])).toFixed(2)}%  p95 ${(100 * q95(all[key])).toFixed(1)}%  max ${(100 * Math.max.apply(null, all[key])).toFixed(1)}%  ` +
    `datasets over 5% flagged: ${all[key].filter(x => x > 0.05).length}/${all[key].length}`);
fs.writeFileSync("clean_bic_only.csv", out.map(r => r.join(",")).join("\n"));
console.log("wrote clean_bic_only.csv");
