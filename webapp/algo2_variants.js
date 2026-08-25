/* algo2_variants.js — can Algorithm 2 be repaired so that a simulation-only
 * calibration also transfers its CUT, not just its ranking?
 *
 * The simulation-only experiment (sim_calibrate.js -> eval_simcal_study1.js)
 * found that the ranking transfers essentially free but the automatic cut does
 * not: with simulation-fitted weights the attentive score is right-skewed
 * (skew ~0.9 against ~0.0 for the shipped weights), the one-sided fit reads its
 * sigma off the left side assuming near-symmetry, so sigma is too small, the cut
 * lands too low and pi-hat sticks near 30% whatever the truth.
 *
 * The parallel review session traced a related fragility to a concrete bug in
 * the same estimator: the mode finder takes the FIRST local maximum above 15% of
 * the peak, so a spurious bump in the low tail can win over the dominant mode,
 * after which sigma is estimated from the handful of points below it. Taking the
 * GLOBAL mode instead was measured there as free (identical on six real
 * datasets, better on clean simulated data). That fix is V1 here.
 *
 * V2/V3 additionally attack the skew problem, by estimating the scale from
 * central quantiles of the whole sample rather than from the left half alone.
 *
 * Every variant is evaluated on three things that must hold together:
 *   (a) the anti-cascade guarantee on genuinely clean simulated data,
 *   (b) the shipped configuration must not get worse,
 *   (c) the simulation-only configuration should get better.
 *
 *   node algo2_variants.js        -> algo2_variants.csv
 */
const fs = require("fs");
const R = require("./site/rerere.js");

const PHI1 = 0.8413447;

/* ---------- the estimator, with its two choices made explicit ---------- */
function kde(eta) {
  const n = eta.length, s = Array.from(eta).sort((a, b) => a - b);
  const q1 = s[Math.floor(0.25 * n)], q3 = s[Math.floor(0.75 * n)], iqr = q3 - q1;
  let mean = 0; for (let i = 0; i < n; i++) mean += eta[i]; mean /= n;
  let v = 0; for (let i = 0; i < n; i++) v += (eta[i] - mean) ** 2; v /= n;
  const sd0 = Math.sqrt(v);
  let h = 0.9 * Math.min(sd0 || 1, iqr > 0 ? iqr / 1.34 : (sd0 || 1)) * Math.pow(n, -0.2);
  if (!(h > 0)) h = (sd0 || 1) * 0.3;
  const G = 512, lo = s[0] - h, hi = s[n - 1] + h, grid = [], dens = [];
  for (let g = 0; g < G; g++) {
    const x = lo + (hi - lo) * g / (G - 1);
    let d = 0; for (let i = 0; i < n; i++) d += Math.exp(-0.5 * ((x - eta[i]) / h) ** 2);
    grid.push(x); dens.push(d);
  }
  return { grid, dens, s, sd0 };
}
const quantile = (sorted, p) => sorted[Math.min(sorted.length - 1, Math.max(0, Math.floor(p * sorted.length)))];
/* standard normal quantile (Acklam), for the central-quantile scale estimator */
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
const CQ_LO = 0.10, CQ_HI = 0.60, CQ_SPAN = qnorm(CQ_HI) - qnorm(CQ_LO);

/* The central-quantile scale is only legitimate while the upper quantile it
 * reads is still attentive. A first-pass left-only fit gives a rough prevalence,
 * which is enough to place that upper quantile safely inside the attentive mass
 * — and to fall back to the left-only scale entirely when the careless are too
 * numerous for any right-hand quantile to be trusted. */
function adaptiveHi(eta, m, sigLeft) {
  const c1 = m + sigLeft;
  let nb = 0; for (const v of eta) if (v < c1) nb++;
  const pi0 = Math.max(0, Math.min(1, 1 - (nb / PHI1) / eta.length));
  return { pi0, hi: Math.min(CQ_HI, 0.95 * (1 - pi0)) };
}

function leftFitVariant(eta, modeRule, sigmaRule) {
  const { grid, dens, s, sd0 } = kde(eta), G = grid.length;
  const dmax = Math.max.apply(null, dens);
  let m = s[Math.floor(s.length / 2)];
  if (modeRule === "firstLocal") {
    for (let g = 1; g < G - 1; g++)
      if (dens[g] >= dens[g - 1] && dens[g] >= dens[g + 1] && dens[g] > 0.15 * dmax) { m = grid[g]; break; }
  } else {                                   // "global": the dominant mode
    let best = -1; for (let g = 0; g < G; g++) if (dens[g] > best) { best = dens[g]; m = grid[g]; }
  }
  const below = Array.from(eta).filter(x => x < m).map(x => m - x).sort((a, b) => a - b);
  let sigLeft = below.length ? 1.4826 * below[Math.floor(below.length / 2)] : 0;
  if (!(sigLeft > 0)) sigLeft = sd0 || 1e-9;
  /* central-quantile scale: measures the spread across the attentive BODY rather
   * than the left half only, so a right-skewed attentive cluster no longer
   * shrinks it. Safe while the careless sit above the CQ_HI quantile. */
  const sigCentral = Math.max(1e-9, (quantile(s, CQ_HI) - quantile(s, CQ_LO)) / CQ_SPAN);
  let sigma = sigLeft, pi0 = null;
  if (sigmaRule === "central") sigma = sigCentral;
  else if (sigmaRule === "max") sigma = Math.max(sigLeft, sigCentral);
  else if (sigmaRule === "switch") {          // hard fallback when the careless are many
    const a = adaptiveHi(eta, m, sigLeft); pi0 = a.pi0;
    sigma = pi0 >= 0.30 ? sigLeft : sigCentral;
  } else if (sigmaRule === "adaptive") {      // slide the upper quantile inside the attentive mass
    const a = adaptiveHi(eta, m, sigLeft); pi0 = a.pi0;
    const span = qnorm(a.hi) - qnorm(CQ_LO);
    sigma = span > 0.2
      ? Math.max(1e-9, (quantile(s, a.hi) - quantile(s, CQ_LO)) / span)
      : sigLeft;
  }
  return { m, sigma, sigLeft, sigCentral, pi0, nBelow: below.length };
}
function autoFlagVariant(eta, z, V) {
  const g = R.gaussMix1D(Array.from(eta));
  const pooled = Math.sqrt((g.v1 + g.v2) / 2) || 1;
  const sep = (g.m2 - g.m1) / pooled;
  const two = (g.bic1 - g.bic2 > 2) && g.pi > 0.01 && g.pi < 0.7 && sep > 1.0;
  const lf = leftFitVariant(eta, V.mode, V.sigma);
  const cut = lf.m + z * lf.sigma;
  const flagged = two ? Array.from(eta, v => v > cut) : eta.map(() => false);
  let pi = 0;
  if (two) {
    const c1 = lf.m + lf.sigma;
    let nb = 0; for (const v of eta) if (v < c1) nb++;
    pi = Math.max(0, Math.min(1, 1 - (nb / PHI1) / eta.length));
  }
  return { flagged, pi, two, cut, lf };
}

const VARIANTS = [
  { tag: "V0 shipped (first-local mode, left MAD)", mode: "firstLocal", sigma: "left" },
  { tag: "V1 global mode, left MAD", mode: "global", sigma: "left" },
  { tag: "V2 global mode, central-quantile scale", mode: "global", sigma: "central" },
  { tag: "V3 global mode, max(left, central)", mode: "global", sigma: "max" },
  { tag: "V4 global mode, central below 30% else left", mode: "global", sigma: "switch" },
  { tag: "V5 global mode, adaptive upper quantile", mode: "global", sigma: "adaptive" }
];

/* ---------- weight vectors under test ---------- */
const Wsim = JSON.parse(fs.readFileSync("sim_calib.json", "utf8")).weights;
const Wship = R.WEIGHTS;

/* ---------- Study 1 prevalence sweep ---------- */
const M = R.parseCSV(fs.readFileSync("_study_matrix.csv", "utf8"));
const J = M.header.length;
const y = R.parseCSV(fs.readFileSync("_study_labels.csv", "utf8")).rows.map(r => Number(r[0]));
const careful = [], careless = [];
y.forEach((v, i) => (v === 1 ? careless : careful).push(i));
const rng = R._rng(20260805);
const sampleK = (pool, k) => {
  const a = pool.slice();
  for (let i = 0; i < k; i++) { const j = i + Math.floor(rng() * (a.length - i)); const t = a[i]; a[i] = a[j]; a[j] = t; }
  return a.slice(0, k);
};
const mccOf = (flag, yy) => {
  let tp = 0, fp = 0, fn = 0, tn = 0;
  for (let i = 0; i < yy.length; i++) { if (yy[i] === 1) { flag[i] ? tp++ : fn++; } else { flag[i] ? fp++ : tn++; } }
  const d = Math.sqrt((tp + fp) * (tp + fn) * (tn + fp) * (tn + fn));
  return { mcc: d ? (tp * tn - fp * fn) / d : 0, prec: tp + fp ? tp / (tp + fp) : 0, rec: tp / (tp + fn) };
};

const RATES = [0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.40, 0.50];
const REPS = 60;
const acc = {};   // acc[variant][weights][rate] = {mcc:[], pi:[], prec:[], rec:[]}
for (const V of VARIANTS) { acc[V.tag] = { sim: {}, shipped: {} }; for (const w of ["sim", "shipped"]) for (const r of RATES) acc[V.tag][w][r] = { mcc: [], pi: [], prec: [], rec: [] }; }

let t0 = Date.now();
for (const rate of RATES) {
  const k = Math.round(careful.length * rate / (1 - rate));
  for (let rep = 0; rep < REPS; rep++) {
    const idx = careful.concat(sampleK(careless, k));
    const n = idx.length, mat = new Float64Array(n * J);
    for (let i = 0; i < n; i++) for (let j = 0; j < J; j++) mat[i * J + j] = R.toNum(M.rows[idx[i]][j]);
    const yy = idx.map(i => y[i]);
    const ef = R.ensembleFeatures(mat, n, J, { seed: 1, iterations: 200, skipD2: true }), O = ef.oriented;
    const gg = Math.max(0, Math.min(1, (ef.profileInfo - 0.01) / 0.05));
    const etaOf = W => Array.from({ length: n }, (_, i) =>
      W.b0 + W.rr * O.rr[i] + gg * (W.longstring * O.longstring[i] + W.person_total * O.person_total[i]));
    const E = { sim: etaOf(Wsim), shipped: etaOf(Wship) };
    for (const V of VARIANTS)
      for (const w of ["sim", "shipped"]) {
        const f = autoFlagVariant(E[w], 2.5, V), c = mccOf(f.flagged, yy);
        const A = acc[V.tag][w][rate];
        A.mcc.push(c.mcc); A.pi.push(f.pi); A.prec.push(c.prec); A.rec.push(c.rec);
      }
  }
  console.log(`  rate ${(100 * rate).toFixed(0)}% done [${((Date.now() - t0) / 1000).toFixed(0)}s]`);
}

/* ---------- anti-cascade check on genuinely clean simulated data ---------- */
function gauss(rr) { let u = 0, v = 0; while (u === 0) u = rr(); while (v === 0) v = rr(); return Math.sqrt(-2 * Math.log(u)) * Math.cos(2 * Math.PI * v); }
function genClean(F, ipf, n, rr) {
  const Jc = F * ipf, load = [], fac = [], tau = [], THR = [];
  for (let c = 1; c < 5; c++) THR.push(qnorm(c / 5));
  for (let j = 0; j < Jc; j++) { let l = 0.4 + 0.4 * rr(); if (rr() < 0.4) l = -l; load.push(l); fac.push(Math.floor(j / ipf)); tau.push(0.45 * gauss(rr)); }
  const rho = 0.3, sq = Math.sqrt(rho), sq1 = Math.sqrt(1 - rho), rows = [];
  for (let i = 0; i < n; i++) {
    const common = gauss(rr), f = [];
    for (let kk = 0; kk < F; kk++) f.push(sq * common + sq1 * gauss(rr));
    const row = new Array(Jc);
    for (let j = 0; j < Jc; j++) {
      const yv = tau[j] + load[j] * f[fac[j]] + Math.sqrt(1 - load[j] * load[j]) * gauss(rr);
      let cat = 1; for (let c = 0; c < 4; c++) if (yv > THR[c]) cat = c + 2;
      row[j] = cat;
    }
    rows.push(row);
  }
  return { rows, J: Jc };
}
const CLEAN = [[10, 6, 200], [10, 6, 500], [20, 6, 500], [30, 6, 500]];
const CLEAN_REPS = 30;
const clean = {};
for (const V of VARIANTS) { clean[V.tag] = { sim: [], shipped: [] }; }
const crng = R._rng(777);
for (const [F, ipf, n] of CLEAN) {
  for (let rep = 0; rep < CLEAN_REPS; rep++) {
    const { rows, J: Jc } = genClean(F, ipf, n, crng);
    const mat = new Float64Array(n * Jc);
    for (let i = 0; i < n; i++) for (let j = 0; j < Jc; j++) mat[i * Jc + j] = rows[i][j];
    const ef = R.ensembleFeatures(mat, n, Jc, { seed: 1, iterations: 150, skipD2: true }), O = ef.oriented;
    const gg = Math.max(0, Math.min(1, (ef.profileInfo - 0.01) / 0.05));
    const etaOf = W => Array.from({ length: n }, (_, i) =>
      W.b0 + W.rr * O.rr[i] + gg * (W.longstring * O.longstring[i] + W.person_total * O.person_total[i]));
    const E = { sim: etaOf(Wsim), shipped: etaOf(Wship) };
    for (const V of VARIANTS) for (const w of ["sim", "shipped"]) {
      const f = autoFlagVariant(E[w], 2.5, V);
      clean[V.tag][w].push(f.flagged.filter(Boolean).length / n);
    }
  }
  console.log(`  clean F=${F} n=${n} done [${((Date.now() - t0) / 1000).toFixed(0)}s]`);
}

/* ---------- report ---------- */
const mn = a => a.reduce((s, v) => s + v, 0) / a.length;
const out = [["variant", "weights", "rate", "mcc", "precision", "recall", "pi_hat", "pi_bias"]];
for (const w of ["shipped", "sim"]) {
  console.log(`\n=== ${w === "shipped" ? "SHIPPED weights (fitted on real data)" : "SIMULATION-ONLY weights"} — Study 1 mixtures, z=2.5 ===`);
  console.log("variant                                      " + RATES.map(r => (100 * r).toFixed(0).padStart(6) + "%").join(""));
  for (const V of VARIANTS) {
    let l1 = "  MCC   " + V.tag.padEnd(38), l2 = "  pi-hat" + " ".repeat(38);
    for (const r of RATES) {
      const A = acc[V.tag][w][r];
      l1 += mn(A.mcc).toFixed(3).padStart(7);
      l2 += (100 * mn(A.pi)).toFixed(0).padStart(6) + "%";
      out.push([V.tag, w, r, mn(A.mcc).toFixed(4), mn(A.prec).toFixed(4), mn(A.rec).toFixed(4), mn(A.pi).toFixed(4), (mn(A.pi) - r).toFixed(4)]);
    }
    console.log(l1); console.log(l2);
  }
}
console.log("\n=== anti-cascade: fraction flagged on genuinely CLEAN simulated data (mean / worst of 120) ===");
for (const V of VARIANTS)
  console.log("  " + V.tag.padEnd(42) +
    ["shipped", "sim"].map(w => `${w}: ${(100 * mn(clean[V.tag][w])).toFixed(2)}% / ${(100 * Math.max.apply(null, clean[V.tag][w])).toFixed(1)}%`).join("   "));
for (const V of VARIANTS) for (const w of ["shipped", "sim"])
  out.push([V.tag, w, "CLEAN", "", "", "", mn(clean[V.tag][w]).toFixed(4), Math.max.apply(null, clean[V.tag][w]).toFixed(4)]);

fs.writeFileSync("algo2_variants.csv", out.map(r => r.join(",")).join("\n"));
console.log("\nwrote algo2_variants.csv");
