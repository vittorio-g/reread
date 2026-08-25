/* sim_calibrate.js — STAGE A of the simulation-only calibration experiment.
 *
 * Question: does the ensemble have to be fitted on real careless data at all?
 * The shipped weights were estimated on Study 1, which is also where the headline
 * numbers come from. Here everything the analyst would otherwise have to choose —
 * the three combiner weights, the intercept, and the flagging multiplier z — is
 * estimated on SYNTHETIC data only. Stage B (eval_simcal_study1.js) then applies
 * the frozen result to the real collected sample, which the fit has never seen.
 *
 * The corpus deliberately spans more than one generator, so what is learned is
 * not a property of one simulation: questionnaire length, sample size, careless
 * rate, loading strength, factor correlation, reverse-keyed share, number of
 * response categories, and homogeneous vs mixed response scales all vary.
 *
 *   node sim_calibrate.js            -> sim_calib_features.csv, sim_calib.json
 */
const fs = require("fs");
const R = require("./site/rerere.js");

const rng = R._rng(20260805);
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
const thresholds = K => { const t = []; for (let c = 1; c < K; c++) t.push(qnorm(c / K)); return t; };

/* CFA clean data, generalised from sim_shipped.js: per-item intercepts (so item
 * means vary and the partner gate opens as it does on real batteries), a share
 * of reverse-keyed items, and per-item response categories, which lets a single
 * dataset mix 4-, 5- and 7-point scales the way a real multi-instrument battery
 * does (Study 1 mixes 1-4, 0-4, 1-5 and 1-7). */
function genClean(cfg) {
  const { F, ipf, n, loadLo, loadHi, rho, revShare, cats } = cfg;
  const J = F * ipf, load = [], fac = [], tau = [], Kj = [], THR = [];
  for (let j = 0; j < J; j++) {
    let l = loadLo + (loadHi - loadLo) * rng();
    if (rng() < revShare) l = -l;
    load.push(l); fac.push(Math.floor(j / ipf)); tau.push(0.45 * gauss());
    const K = cats[Math.floor(rng() * cats.length)];
    Kj.push(K); THR.push(thresholds(K));
  }
  const sq = Math.sqrt(rho), sq1 = Math.sqrt(1 - rho), rows = [];
  for (let i = 0; i < n; i++) {
    const common = gauss(), f = [];
    for (let k = 0; k < F; k++) f.push(sq * common + sq1 * gauss());
    const row = new Array(J);
    for (let j = 0; j < J; j++) {
      const y = tau[j] + load[j] * f[fac[j]] + Math.sqrt(1 - load[j] * load[j]) * gauss();
      let cat = 1; for (let c = 0; c < Kj[j] - 1; c++) if (y > THR[j][c]) cat = c + 2;
      row[j] = cat;
    }
    rows.push(row);
  }
  return { rows, J, Kj };
}

/* six careless patterns at six corruption levels — the port of the injector the
 * tool's own demo used, but respecting each item's own response range so mixed
 * scales stay legal */
function inject(rows, J, Kj, pct) {
  const N = rows.length;
  const lev = j => { const a = []; const lo = Math.min(...rows.map(r => r[j])); for (let v = 1; v <= Kj[j]; v++) a.push(v); return lo === 0 ? a.map(v => v - 1) : a; };
  const LEV = []; for (let j = 0; j < J; j++) LEV.push(lev(j));
  const ri = m => Math.floor(rng() * m), pick = a => a[ri(a.length)];
  const sampleK = (pool, k) => { const a = pool.slice(); for (let i = 0; i < k; i++) { const j = i + ri(a.length - i); const t = a[i]; a[i] = a[j]; a[j] = t; } return a.slice(0, k); };
  const allIdx = [...Array(J).keys()];
  const hiOf = j => LEV[j][LEV[j].length - 1], loOf = j => LEV[j][0];
  const clampJ = (j, x) => Math.max(loOf(j), Math.min(hiOf(j), x));
  /* a straight-line / acquiescent respondent picks a POSITION on the scale, not a
   * raw number, so the same behaviour maps onto items with different ranges */
  const atPos = (j, pos) => LEV[j][Math.min(LEV[j].length - 1, Math.max(0, Math.round(pos * (LEV[j].length - 1))))];
  const chunks = (row, idx, k, seq) => {
    const nc = k >= 4 ? 2 + ri(Math.min(5, k) - 1) : (k >= 2 ? 2 : 1);
    const asg = []; for (let t = 0; t < k; t++) asg.push(t % nc);
    if (!seq) for (let t = k - 1; t > 0; t--) { const j = ri(t + 1); const tmp = asg[t]; asg[t] = asg[j]; asg[j] = tmp; }
    for (let ch = 0; ch < nc; ch++) { const pos = rng(); for (let t = 0; t < k; t++) if (asg[t] === ch) row[idx[t]] = atPos(idx[t], pos); }
  };
  const patterns = ["random", "longstring", "pure_straight", "acquiescent", "mixed", "fatigue"];
  const levels = [0.5, 0.6, 0.7, 0.8, 0.9, 1.0];
  const M = Math.round(N * pct), order = sampleK([...Array(N).keys()], N), chosen = order.slice(0, M);
  const labels = new Array(N).fill(0), sev = new Array(N).fill(0), out = rows.map(r => r.slice());
  for (let m = 0; m < M; m++) {
    const i = chosen[m], pat = patterns[m % patterns.length], lvl = pick(levels);
    const k = Math.min(J, Math.max(1, Math.round(J * lvl))), row = out[i];
    if (pat === "random") sampleK(allIdx, k).forEach(p => row[p] = pick(LEV[p]));
    else if (pat === "longstring") chunks(row, sampleK(allIdx, k), k, false);
    else if (pat === "pure_straight") { const pos = rng(); sampleK(allIdx, k).forEach(p => row[p] = atPos(p, pos)); }
    else if (pat === "acquiescent") { const pos = 0.7 + 0.3 * rng(); sampleK(allIdx, k).forEach(p => row[p] = clampJ(p, atPos(p, pos) + pick([-1, 0, 0, 0, 1]))); }
    else if (pat === "mixed") { const kl = Math.max(1, Math.round(k / 2)), idx = sampleK(allIdx, k); chunks(row, idx.slice(0, kl), kl, false); idx.slice(kl).forEach(p => row[p] = pick(LEV[p])); }
    else if (pat === "fatigue") { const idx = []; for (let t = J - k; t < J; t++) idx.push(t); chunks(row, idx, k, true); }
    labels[i] = 1; sev[i] = lvl;
  }
  return { out, labels, sev };
}

/* ---------------- corpus design ---------------- */
const SIZES = [[8, 6], [10, 6], [12, 8], [18, 6], [20, 6], [25, 6], [30, 6]];   // 48..180 items
const NS = [150, 300, 500];
const RATES = [0.05, 0.15, 0.30, 0.50];
const REPS = 2;
/* structural variants cycled across datasets so no single set of assumptions
 * dominates the fit */
const VARIANTS = [
  { loadLo: 0.4, loadHi: 0.8, rho: 0.30, revShare: 0.40, cats: [5], tag: "standard" },
  { loadLo: 0.2, loadHi: 0.5, rho: 0.30, revShare: 0.40, cats: [5], tag: "weak-loadings" },
  { loadLo: 0.4, loadHi: 0.8, rho: 0.60, revShare: 0.40, cats: [5], tag: "correlated-factors" },
  { loadLo: 0.4, loadHi: 0.8, rho: 0.30, revShare: 0.00, cats: [5], tag: "no-reverse-keying" },
  { loadLo: 0.4, loadHi: 0.8, rho: 0.30, revShare: 0.40, cats: [7], tag: "7-point" },
  { loadLo: 0.4, loadHi: 0.8, rho: 0.30, revShare: 0.40, cats: [4, 5, 7], tag: "mixed-scales" }
];

/* skipD2: Mahalanobis carries zero ensemble weight and is not used here, and its
 * eigendecomposition is the expensive part at J up to 180 with n as low as 150 */
const OPTS = { seed: 1, iterations: 150, skipD2: true };
const rowsOut = [["dataset", "variant", "F", "ipf", "J", "n", "rate", "profileInfo", "gate", "rr", "ls", "pt", "label", "severity"]];
const pooled = { X: [], y: [], ds: [] };

let d = 0, t0 = Date.now();
for (const [F, ipf] of SIZES)
  for (const n of NS)
    for (const rate of RATES)
      for (let rep = 0; rep < REPS; rep++) {
        const v = VARIANTS[d % VARIANTS.length];
        const cfg = Object.assign({ F, ipf, n }, v);
        const { rows, J, Kj } = genClean(cfg);
        const { out, labels, sev } = inject(rows, J, Kj, rate);
        const mat = new Float64Array(n * J);
        for (let i = 0; i < n; i++) for (let j = 0; j < J; j++) mat[i * J + j] = out[i][j];
        const ef = R.ensembleFeatures(mat, n, J, OPTS), O = ef.oriented;
        const g = Math.max(0, Math.min(1, (ef.profileInfo - 0.01) / 0.05));
        for (let i = 0; i < n; i++) {
          rowsOut.push([d, v.tag, F, ipf, J, n, rate, ef.profileInfo.toFixed(4), g.toFixed(4),
            O.rr[i].toFixed(4), O.longstring[i].toFixed(4), O.person_total[i].toFixed(4), labels[i], sev[i]]);
          /* design matrix in DEPLOYMENT form: the partners enter through the gate,
           * exactly as eta is assembled at scoring time */
          pooled.X.push([1, O.rr[i], g * O.longstring[i], g * O.person_total[i]]);
          pooled.y.push(labels[i]);
          pooled.ds.push(d);
        }
        d++;
        if (d % 20 === 0) console.log(`  ${d} datasets, ${pooled.y.length} respondents, ${((Date.now() - t0) / 1000).toFixed(0)}s`);
      }

console.log(`corpus : ${d} simulated datasets, ${pooled.y.length} respondents, ` +
  `${pooled.y.reduce((s, v) => s + v, 0)} careless (${(100 * pooled.y.reduce((s, v) => s + v, 0) / pooled.y.length).toFixed(1)}%)`);

/* ---------------- fit the combiner (same ridge logistic as the shipped fit) ---------------- */
const P = 4;
function solve(A, b) {
  const m = b.length, Aug = A.map((r, i) => r.concat(b[i]));
  for (let c = 0; c < m; c++) {
    let pv = c; for (let r = c + 1; r < m; r++) if (Math.abs(Aug[r][c]) > Math.abs(Aug[pv][c])) pv = r;
    [Aug[c], Aug[pv]] = [Aug[pv], Aug[c]];
    const dd = Aug[c][c] || 1e-12;
    for (let r = 0; r < m; r++) if (r !== c) { const f = Aug[r][c] / dd; for (let k = c; k <= m; k++) Aug[r][k] -= f * Aug[c][k]; }
  }
  return Aug.map((r, i) => r[m] / (r[i] || 1e-12));
}
function fit(X, y, lambda) {
  let b = new Array(P).fill(0);
  for (let it = 0; it < 60; it++) {
    const g = new Array(P).fill(0), H = Array.from({ length: P }, () => new Array(P).fill(0));
    for (let i = 0; i < X.length; i++) {
      let e = 0; for (let k = 0; k < P; k++) e += b[k] * X[i][k];
      const p = 1 / (1 + Math.exp(-e)), w = Math.max(p * (1 - p), 1e-6);
      for (let k = 0; k < P; k++) { g[k] += (y[i] - p) * X[i][k]; for (let l = 0; l < P; l++) H[k][l] += w * X[i][k] * X[i][l]; }
    }
    for (let k = 1; k < P; k++) { g[k] -= lambda * b[k]; H[k][k] += lambda; }
    const s = solve(H, g); let mx = 0;
    for (let k = 0; k < P; k++) { b[k] += s[k]; mx = Math.max(mx, Math.abs(s[k])); }
    if (mx < 1e-9) break;
  }
  return b;
}
const beta = fit(pooled.X, pooled.y, 1.0);
const Wsim = { b0: +beta[0].toFixed(4), rr: +beta[1].toFixed(4), irv: 0, longstring: +beta[2].toFixed(4), d2: 0, person_total: +beta[3].toFixed(4) };
console.log("weights (simulation-only): " + JSON.stringify(Wsim));

/* ---------------- pick the flagging multiplier z on simulation too ---------------- */
function mcc(flag, lab) {
  let tp = 0, fp = 0, fn = 0, tn = 0;
  for (let i = 0; i < lab.length; i++) { if (lab[i] === 1) { flag[i] ? tp++ : fn++; } else { flag[i] ? fp++ : tn++; } }
  const dd = Math.sqrt((tp + fp) * (tp + fn) * (tn + fp) * (tn + fn));
  return dd ? (tp * tn - fp * fn) / dd : 0;
}
const byDs = new Map();
pooled.ds.forEach((k, i) => { if (!byDs.has(k)) byDs.set(k, []); byDs.get(k).push(i); });
const ZGRID = []; for (let z = 1.0; z <= 3.51; z += 0.25) ZGRID.push(+z.toFixed(2));
const zScores = ZGRID.map(z => {
  const per = [];
  for (const idx of byDs.values()) {
    const eta = idx.map(i => pooled.X[i].reduce((s, v, k) => s + v * beta[k], 0));
    const lab = idx.map(i => pooled.y[i]);
    per.push(mcc(R.autoFlag(eta, z).flagged, lab));
  }
  return { z, mcc: per.reduce((s, v) => s + v, 0) / per.length };
});
zScores.sort((a, b) => b.mcc - a.mcc);
const zBest = zScores[0];
console.log("z sweep (mean per-dataset MCC on simulation):");
zScores.slice().sort((a, b) => a.z - b.z).forEach(r =>
  console.log(`   z=${r.z.toFixed(2)}  MCC=${r.mcc.toFixed(3)}${r.z === zBest.z ? "   <- best" : (r.z === 2.5 ? "   (shipped)" : "")}`));

fs.writeFileSync("sim_calib_features.csv", rowsOut.map(r => r.join(",")).join("\n"));
fs.writeFileSync("sim_calib.json", JSON.stringify({
  weights: Wsim, zBest: zBest.z, zSweep: zScores.slice().sort((a, b) => a.z - b.z),
  corpus: { datasets: d, respondents: pooled.y.length, sizes: SIZES, ns: NS, rates: RATES, reps: REPS, variants: VARIANTS.map(v => v.tag), opts: OPTS }
}, null, 2));
console.log(`\nwrote sim_calib_features.csv (${rowsOut.length - 1} rows) and sim_calib.json  [${((Date.now() - t0) / 1000).toFixed(0)}s]`);
