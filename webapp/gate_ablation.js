/* gate_ablation.js — two engine questions at once, on the same footing.
 *
 * (1) Scale: the shipped sigma is the left half-width at half maximum of the KDE
 *     peak; the previous version used 1.4826 * MAD of the points below the mode.
 *     HWHM was adopted because the MAD is inflated by eta's long coherent left
 *     tail on some datasets, which pushed the cut past everyone. On Study 1 it
 *     goes the other way (sigma 1.51 -> 2.02, 72 flags -> 60). Which is better
 *     overall, and where?
 *
 * (2) Gate: the structure gate exists to guarantee that clean data yields no
 *     flags. It also produces the awkward states -- abstaining on scores that do
 *     carry signal, and (because gate and cut are estimated from different
 *     objects) opening while nobody passes the cut. What if it were removed and
 *     the anti-cascade guarantee were replaced by a warning?
 *
 * Six configurations are scored on three things that have to be weighed together:
 * external decision quality, behaviour across prevalence on the real study, and
 * the price paid on genuinely clean data.
 *
 *   node gate_ablation.js   -> gate_ablation.csv
 */
const fs = require("fs"), path = require("path");
const R = require("./site/rerere.js");
const Z = 2.5, P_TAIL = 0.0062, TAIL_ALPHA = 1e-6;

/* upper-tail binomial p-value in log space (the engine's own test, reimplemented
 * here so the tail branch can be recomputed under a different sigma) */
function lgamma(x) {
  const c = [76.18009172947146, -86.50532032941677, 24.01409824083091,
             -1.231739572450155, 0.1208650973866179e-2, -0.5395239384953e-5];
  let y = x, t = x + 5.5; t -= (x + 0.5) * Math.log(t);
  let s = 1.000000000190015;
  for (let j = 0; j < 6; j++) s += c[j] / ++y;
  return -t + Math.log(2.5066282746310005 * s / x);
}
function binomUpperP(k, n, p) {
  if (k <= 0) return 1;
  if (k > n) return 0;
  const lch = (n, i) => lgamma(n + 1) - lgamma(i + 1) - lgamma(n - i + 1);
  const terms = []; let mx = -Infinity;
  for (let i = k; i <= n; i++) {
    const lt = lch(n, i) + i * Math.log(p) + (n - i) * Math.log1p(-p);
    terms.push(lt); if (lt > mx) mx = lt;
    if (i > k + 5 && lt < mx - 50) break;
  }
  let s = 0; for (const t of terms) s += Math.exp(t - mx);
  return Math.exp(mx) * s;
}

const CONFIGS = [
  { tag: "shipped: HWHM + gate(BIC|tail)", sigma: "hwhm", gate: "both" },
  { tag: "HWHM + gate(BIC only)",          sigma: "hwhm", gate: "bic"  },
  { tag: "HWHM + NO gate",                 sigma: "hwhm", gate: "none" },
  { tag: "previous: MAD + gate(BIC only)", sigma: "mad",  gate: "bic"  },
  { tag: "MAD + gate(BIC|tail)",           sigma: "mad",  gate: "both" },
  { tag: "MAD + NO gate",                  sigma: "mad",  gate: "none" }
];

/* one decision for one score under one configuration */
function decide(score, cfg) {
  const s = Array.from(score), n = s.length;
  const lf = R.leftFit(s);
  const sig = cfg.sigma === "hwhm" ? lf.sigma : lf.sigmaMad;
  const cut = lf.m + Z * sig;
  let kTail = 0; for (const v of s) if (v > cut) kTail++;
  const af = R.autoFlag(s, "standard");        // for its BIC-branch verdict only
  const tailOpen = binomUpperP(kTail, n, P_TAIL) < TAIL_ALPHA;
  const open = cfg.gate === "none" ? true
             : cfg.gate === "bic" ? af.bicGate
             : (af.bicGate || tailOpen);
  const flagged = open ? s.map(v => v > cut) : s.map(() => false);
  return { flagged, open, cut, sigma: sig, m: lf.m };
}
function conf(flag, y) {
  let tp = 0, fp = 0, fn = 0, tn = 0;
  for (let i = 0; i < y.length; i++) { if (y[i] === 1) { flag[i] ? tp++ : fn++; } else { flag[i] ? fp++ : tn++; } }
  const d = Math.sqrt((tp + fp) * (tp + fn) * (tn + fp) * (tn + fn));
  return { mcc: d ? (tp * tn - fp * fn) / d : 0, prec: tp + fp ? tp / (tp + fp) : NaN,
           rec: tp / (tp + fn), share: (tp + fp) / y.length };
}
const mean = a => { const v = a.filter(Number.isFinite); return v.length ? v.reduce((s, x) => s + x, 0) / v.length : NaN; };

/* ================= (a) external datasets ================= */
const DIR = "ext_bench2";
const ORDER = ["kay_s2","warning300","kay_s1","opsy_16pf","smarvus","kay_s6","kay_s5","douglas","krause"];
function readY(nm) {
  const L = fs.readFileSync(path.join(DIR, nm + "_y.csv"), "utf8").trim().split(/\r?\n/);
  const h = L[0].split(",").map(x => x.replace(/"/g, "")), c = {};
  h.forEach(x => c[x] = []);
  for (let i = 1; i < L.length; i++) { const f = L[i].split(","); h.forEach((x, k) => c[x].push(Number(f[k]))); }
  return c;
}
const EXT = {};
for (const cfg of CONFIGS) EXT[cfg.tag] = {};
for (const nm of ORDER) {
  const Y = readY(nm);
  for (const cfg of CONFIGS) {
    const d = decide(Y.eta, cfg);
    EXT[cfg.tag][nm] = Object.assign(conf(d.flagged, Y.y), { open: d.open });
  }
}
console.log("=== (a) nine convergent external datasets, ensemble score ===");
console.log("configuration                       mean MCC   mean %flagged   gate shut   flags nobody");
for (const cfg of CONFIGS) {
  const E = ORDER.map(d => EXT[cfg.tag][d]);
  console.log(`${cfg.tag.padEnd(34)} ${mean(E.map(e => e.mcc)).toFixed(3)}      ` +
    `${(100 * mean(E.map(e => e.share))).toFixed(1)}%          ` +
    `${E.filter(e => !e.open).length}/9        ${E.filter(e => e.share === 0).length}/9`);
}
console.log("\nper dataset (MCC):");
console.log("configuration                     " + ORDER.map(d => d.slice(0, 7).padStart(8)).join(""));
for (const cfg of CONFIGS)
  console.log(cfg.tag.padEnd(34) + ORDER.map(d => EXT[cfg.tag][d].mcc.toFixed(3).padStart(8)).join(""));

/* ================= (b) Study 1 across prevalence ================= */
const M = R.parseCSV(fs.readFileSync("_study_matrix.csv", "utf8"));
const J1 = M.header.length;
const y1 = R.parseCSV(fs.readFileSync("_study_labels.csv", "utf8")).rows.map(r => Number(r[0]));
const careful = [], careless = [];
y1.forEach((v, i) => (v === 1 ? careless : careful).push(i));
const rng = R._rng(20260819);
const sampleK = (pool, k) => {
  const a = pool.slice();
  for (let i = 0; i < k; i++) { const j = i + Math.floor(rng() * (a.length - i)); const t = a[i]; a[i] = a[j]; a[j] = t; }
  return a.slice(0, k);
};
const RATES = [0.05, 0.10, 0.15, 0.20, 0.30, 0.50];
const REPS = 50;
const S1 = {}; for (const cfg of CONFIGS) { S1[cfg.tag] = {}; for (const r of RATES) S1[cfg.tag][r] = { mcc: [], share: [] }; }
for (const rate of RATES) {
  const k = Math.round(careful.length * rate / (1 - rate));
  for (let rep = 0; rep < REPS; rep++) {
    const idx = careful.concat(sampleK(careless, k));
    const n = idx.length, mat = new Float64Array(n * J1);
    for (let i = 0; i < n; i++) for (let j = 0; j < J1; j++) mat[i * J1 + j] = R.toNum(M.rows[idx[i]][j]);
    const yy = idx.map(i => y1[i]);
    const ens = R.ensemble(mat, n, J1, { seed: 1, iterations: 200, skipD2: true });
    for (const cfg of CONFIGS) {
      const c = conf(decide(ens.eta, cfg).flagged, yy);
      S1[cfg.tag][rate].mcc.push(c.mcc); S1[cfg.tag][rate].share.push(c.share);
    }
  }
  process.stdout.write(`  [study1 ${(100 * rate).toFixed(0)}% done]`);
}
console.log("\n\n=== (b) Study 1 mixtures, MCC of the automatic flags (50 resamples per rate) ===");
console.log("configuration                     " + RATES.map(r => ((100 * r).toFixed(0) + "%").padStart(8)).join("") + "     mean");
for (const cfg of CONFIGS) {
  const v = RATES.map(r => mean(S1[cfg.tag][r].mcc));
  console.log(cfg.tag.padEnd(34) + v.map(x => x.toFixed(3).padStart(8)).join("") + "   " + mean(v).toFixed(3));
}

/* ================= (c) genuinely clean data ================= */
function gauss(rr) { let u = 0, v = 0; while (u === 0) u = rr(); while (v === 0) v = rr(); return Math.sqrt(-2 * Math.log(u)) * Math.cos(2 * Math.PI * v); }
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
const CL = {}; for (const cfg of CONFIGS) CL[cfg.tag] = [];
const crng = R._rng(777);
for (const [F, ipf, n] of CLEAN) {
  for (let rep = 0; rep < CLEAN_REPS; rep++) {
    const { rows, J: Jc } = genClean(F, ipf, n, crng);
    const mat = new Float64Array(n * Jc);
    for (let i = 0; i < n; i++) for (let j = 0; j < Jc; j++) mat[i * Jc + j] = rows[i][j];
    const ens = R.ensemble(mat, n, Jc, { seed: 1, iterations: 150, skipD2: true });
    for (const cfg of CONFIGS) {
      const d = decide(ens.eta, cfg);
      CL[cfg.tag].push(d.flagged.filter(Boolean).length / n);
    }
  }
}
console.log("\n=== (c) genuinely CLEAN simulated data (120 datasets): share wrongly flagged ===");
console.log("configuration                       mean     worst    datasets with >5% flagged");
for (const cfg of CONFIGS) {
  const v = CL[cfg.tag];
  console.log(`${cfg.tag.padEnd(34)} ${(100 * mean(v)).toFixed(2)}%   ${(100 * Math.max.apply(null, v)).toFixed(1)}%     ${v.filter(x => x > 0.05).length}/120`);
}

const out = [["configuration","ext_mean_mcc","ext_mean_share","ext_gate_shut","s1_mean_mcc","clean_mean","clean_worst","clean_over5pct"]];
for (const cfg of CONFIGS) {
  const E = ORDER.map(d => EXT[cfg.tag][d]);
  out.push([cfg.tag, mean(E.map(e => e.mcc)).toFixed(4), mean(E.map(e => e.share)).toFixed(4),
    E.filter(e => !e.open).length, mean(RATES.map(r => mean(S1[cfg.tag][r].mcc))).toFixed(4),
    mean(CL[cfg.tag]).toFixed(4), Math.max.apply(null, CL[cfg.tag]).toFixed(4),
    CL[cfg.tag].filter(x => x > 0.05).length]);
}
fs.writeFileSync("gate_ablation.csv", out.map(r => r.join(",")).join("\n"));
console.log("\nwrote gate_ablation.csv");
