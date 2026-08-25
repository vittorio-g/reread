/* clean_envelope_ab.js — the cost of removing the gate, measured on identical
 * clean data with both engines side by side.
 *
 * The earlier comparison used one structural variant and mostly large samples,
 * and put the cost at a few tenths of a percentage point. This repeats it on the
 * grid the envelope is calibrated on, which also varies loadings, factor
 * correlation, reverse-keying and response categories -- the conditions under
 * which eta's geometry is least well behaved and the gate would have the most to
 * do. Both engines score the SAME generated matrices.
 *
 *   node clean_envelope_ab.js   -> clean_envelope_ab.csv
 */
const fs = require("fs");
const NEW = require("./site/reread.js");
const OLD = require("./_engine_gated.js");

const rng = NEW._rng(4242);          // same seed as clean_envelope.js
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
const out = [["J", "n", "gateless_mean", "gateless_p95", "gated_mean", "gated_p95", "gated_shut_frac"]];
const mn = a => a.reduce((s, v) => s + v, 0) / a.length;
const q95 = a => { const s = [...a].sort((x, y) => x - y); return s[Math.floor(0.95 * s.length)]; };
console.log("            gateless (now)        with the old gate      gate shut");
console.log("  J    n     mean    p95           mean    p95           on");
const acc = { gl: [], gd: [] };
for (const J of JS) for (const n of NS) {
  const a = [], b = []; let shut = 0;
  for (let rep = 0; rep < REPS; rep++) {
    const M = genClean(J, n, VARIANTS[rep % VARIANTS.length]);
    const eNew = NEW.ensemble(M, n, J, { seed: 1, iterations: 150, skipD2: true });
    const eOld = OLD.ensemble(M, n, J, { seed: 1, iterations: 150, skipD2: true });
    a.push(eNew.nFlagged / n); b.push(eOld.nFlagged / n);
    if (eOld.nFlagged === 0 && eNew.nFlagged > 0) shut++;
  }
  acc.gl.push(...a); acc.gd.push(...b);
  console.log(`${String(J).padStart(3)} ${String(n).padStart(5)}   ${(100 * mn(a)).toFixed(2).padStart(5)}%  ${(100 * q95(a)).toFixed(1).padStart(5)}%        ` +
    `${(100 * mn(b)).toFixed(2).padStart(5)}%  ${(100 * q95(b)).toFixed(1).padStart(5)}%        ${shut}/${REPS}`);
  out.push([J, n, mn(a).toFixed(5), q95(a).toFixed(5), mn(b).toFixed(5), q95(b).toFixed(5), (shut / REPS).toFixed(3)]);
}
console.log(`\noverall on clean data: gateless mean ${(100 * mn(acc.gl)).toFixed(2)}%  p95 ${(100 * q95(acc.gl)).toFixed(1)}%  |  ` +
  `gated mean ${(100 * mn(acc.gd)).toFixed(2)}%  p95 ${(100 * q95(acc.gd)).toFixed(1)}%`);
fs.writeFileSync("clean_envelope_ab.csv", out.map(r => r.join(",")).join("\n"));
console.log("wrote clean_envelope_ab.csv");
