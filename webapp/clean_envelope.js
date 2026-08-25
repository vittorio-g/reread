/* clean_envelope.js — what replaces the anti-cascade guarantee.
 *
 * With the structure gate removed the cut is unconditional, so a clean sample no
 * longer flags exactly nobody: it flags whatever its own upper tail puts past
 * m + 2.5*sigma. That is not a defect to hide but a quantity to report, and it
 * depends almost entirely on the two things the analyst already knows -- how many
 * items and how many respondents. This measures it on genuinely clean simulated
 * data across a grid of both, so the tool can tell a user what share of a clean
 * sample of their size the same rule would flag by chance.
 *
 * A flagged share inside that envelope means "indistinguishable from clean"; the
 * old gate asserted the same thing by refusing to act, which threw away the cases
 * where the share was genuinely above it.
 *
 *   node clean_envelope.js   -> clean_envelope.json, clean_envelope.csv
 */
const fs = require("fs");
const R = require("./site/rerere.js");

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
/* clean CFA data; the structural knobs are varied across replications so the
 * envelope is not a property of one generator */
const VARIANTS = [
  { loadLo: .4, loadHi: .8, rho: .30, rev: .40, K: 5 },
  { loadLo: .2, loadHi: .5, rho: .30, rev: .40, K: 5 },
  { loadLo: .4, loadHi: .8, rho: .60, rev: .40, K: 5 },
  { loadLo: .4, loadHi: .8, rho: .30, rev: .00, K: 5 },
  { loadLo: .4, loadHi: .8, rho: .30, rev: .40, K: 7 }
];
function genClean(J, n, v) {
  const ipf = 6, F = Math.max(2, Math.round(J / ipf));
  const load = [], fac = [], tau = [], THR = [];
  for (let c = 1; c < v.K; c++) THR.push(qnorm(c / v.K));
  for (let j = 0; j < J; j++) {
    let l = v.loadLo + (v.loadHi - v.loadLo) * rng(); if (rng() < v.rev) l = -l;
    load.push(l); fac.push(j % F); tau.push(0.45 * gauss());
  }
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

const JS = [30, 48, 60, 90, 120, 180, 240, 300];
const NS = [100, 150, 200, 300, 500, 1000];
const REPS = 30;
const grid = [];
const csv = [["J", "n", "reps", "mean", "p50", "p90", "p95", "max", "share_over_5pct"]];
const t0 = Date.now();
for (const J of JS) {
  for (const n of NS) {
    const shares = [];
    for (let rep = 0; rep < REPS; rep++) {
      const v = VARIANTS[rep % VARIANTS.length];
      const M = genClean(J, n, v);
      const ens = R.ensemble(M, n, J, { seed: 1, iterations: 150, skipD2: true });
      shares.push(ens.nFlagged / n);
    }
    shares.sort((a, b) => a - b);
    const q = p => shares[Math.min(shares.length - 1, Math.floor(p * shares.length))];
    const mean = shares.reduce((s, v) => s + v, 0) / shares.length;
    grid.push({ J, n, mean: +mean.toFixed(5), p95: +q(0.95).toFixed(5) });
    csv.push([J, n, REPS, mean.toFixed(5), q(0.5).toFixed(5), q(0.9).toFixed(5), q(0.95).toFixed(5),
      shares[shares.length - 1].toFixed(5), (shares.filter(x => x > 0.05).length / shares.length).toFixed(3)]);
    console.log(`J=${String(J).padStart(3)} n=${String(n).padStart(4)}  mean ${(100 * mean).toFixed(2)}%  p95 ${(100 * q(0.95)).toFixed(1)}%  max ${(100 * shares[shares.length - 1]).toFixed(1)}%   [${((Date.now() - t0) / 1000).toFixed(0)}s]`);
  }
}
fs.writeFileSync("clean_envelope.csv", csv.map(r => r.join(",")).join("\n"));
fs.writeFileSync("clean_envelope.json", JSON.stringify({ js: JS, ns: NS, reps: REPS, grid }, null, 1));
console.log("\nwrote clean_envelope.json and clean_envelope.csv");
