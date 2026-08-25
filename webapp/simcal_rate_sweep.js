/* simcal_rate_sweep.js — STAGE C of the simulation-only calibration experiment.
 *
 * Study 1 as collected is 55% careless by design, which is not a deployment
 * condition. This re-runs the comparison across realistic prevalences: keep all
 * 70 careful respondents, draw k real careless into the mixture to hit a target
 * rate, RECOMPUTE the engine on that mixture (pair selection, sign alignment and
 * the robust-z reference all depend on the sample), then score it with each
 * FROZEN weight vector. Nothing is refitted inside the loop.
 *
 * The comparison is deliberately unfair to the simulation-only fit: the shipped
 * weights were estimated on these very respondents' labels, so they have seen
 * every resample. Sim-only has seen none of them.
 *
 *   node simcal_rate_sweep.js        -> simcal_rate_sweep.csv
 */
const fs = require("fs");
const R = require("./site/reread.js");

const CAL = JSON.parse(fs.readFileSync("sim_calib.json", "utf8"));
const Wsim = CAL.weights, Zsim = CAL.zBest, Wship = R.WEIGHTS;

const M = R.parseCSV(fs.readFileSync("_study_matrix.csv", "utf8"));
const J = M.header.length;
const y = R.parseCSV(fs.readFileSync("_study_labels.csv", "utf8")).rows.map(r => Number(r[0]));
const carefulIdx = [], carelessIdx = [];
y.forEach((v, i) => (v === 1 ? carelessIdx : carefulIdx).push(i));

const rng = R._rng(20260805);
const sampleK = (pool, k) => {
  const a = pool.slice();
  for (let i = 0; i < k; i++) { const j = i + Math.floor(rng() * (a.length - i)); const t = a[i]; a[i] = a[j]; a[j] = t; }
  return a.slice(0, k);
};
const auc = (s, yy) => {
  const po = [], ne = [];
  for (let i = 0; i < yy.length; i++) (yy[i] ? po : ne).push(s[i]);
  let c = 0; for (const a of po) for (const b of ne) c += a > b ? 1 : a === b ? 0.5 : 0;
  return c / (po.length * ne.length);
};
const mccOf = (flag, yy) => {
  let tp = 0, fp = 0, fn = 0, tn = 0;
  for (let i = 0; i < yy.length; i++) { if (yy[i] === 1) { flag[i] ? tp++ : fn++; } else { flag[i] ? fp++ : tn++; } }
  const d = Math.sqrt((tp + fp) * (tp + fn) * (tn + fp) * (tn + fn));
  return { mcc: d ? (tp * tn - fp * fn) / d : 0, prec: tp + fp ? tp / (tp + fp) : 0, rec: tp / (tp + fn) };
};

const RATES = [0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.40, 0.50];
const REPS = 100;
const OPTS = { seed: 1, iterations: 200, skipD2: true };
const acc = {};
const t0 = Date.now();

for (const rate of RATES) {
  const k = Math.round(carefulIdx.length * rate / (1 - rate));
  const A = { auc_sim: [], auc_ship: [], mcc_sim: [], mcc_ship: [], prec_sim: [], rec_sim: [], prec_ship: [], rec_ship: [], pi_sim: [], pi_ship: [], auc_rr: [], auc_pt: [], mcc_sim_z25: [], mcc_ship_zsim: [], skew_sim: [], skew_ship: [] };
  for (let rep = 0; rep < REPS; rep++) {
    const idx = carefulIdx.concat(sampleK(carelessIdx, k));
    const n = idx.length, mat = new Float64Array(n * J);
    for (let i = 0; i < n; i++) for (let j = 0; j < J; j++) mat[i * J + j] = R.toNum(M.rows[idx[i]][j]);
    const yy = idx.map(i => y[i]);
    const ef = R.ensembleFeatures(mat, n, J, OPTS), O = ef.oriented;
    const g = Math.max(0, Math.min(1, (ef.profileInfo - 0.01) / 0.05));
    const eta = W => Array.from({ length: n }, (_, i) =>
      W.b0 + W.rr * O.rr[i] + g * (W.longstring * O.longstring[i] + W.person_total * O.person_total[i]));
    const eS = eta(Wsim), eH = eta(Wship);
    const fS = R.autoFlag(eS, Zsim), fH = R.autoFlag(eH, 2.5);
    const cS = mccOf(fS.flagged, yy), cH = mccOf(fH.flagged, yy);
    /* separate the two things simulation had to guess: the weights and the cut.
     * sim weights at the shipped z isolates the weights; shipped weights at the
     * sim z isolates the cut. */
    const fS25 = R.autoFlag(eS, 2.5), fH20 = R.autoFlag(eH, Zsim);
    A.mcc_sim_z25.push(mccOf(fS25.flagged, yy).mcc);
    A.mcc_ship_zsim.push(mccOf(fH20.flagged, yy).mcc);
    /* skewness of the attentive-only score, the quantity the one-sided fit
     * assumes is near-normal when it reads sigma off the left side */
    const att = eS.filter((_, i) => yy[i] === 0), attH = eH.filter((_, i) => yy[i] === 0);
    const skew = a => { const m = a.reduce((s, v) => s + v, 0) / a.length;
      const v2 = a.reduce((s, v) => s + (v - m) ** 2, 0) / a.length;
      const v3 = a.reduce((s, v) => s + (v - m) ** 3, 0) / a.length;
      return v2 > 0 ? v3 / Math.pow(v2, 1.5) : 0; };
    A.skew_sim.push(skew(att)); A.skew_ship.push(skew(attH));
    A.auc_sim.push(auc(eS, yy)); A.auc_ship.push(auc(eH, yy));
    A.auc_rr.push(auc(Array.from(O.rr), yy)); A.auc_pt.push(auc(Array.from(O.person_total), yy));
    A.mcc_sim.push(cS.mcc); A.mcc_ship.push(cH.mcc);
    A.prec_sim.push(cS.prec); A.rec_sim.push(cS.rec);
    A.prec_ship.push(cH.prec); A.rec_ship.push(cH.rec);
    A.pi_sim.push(fS.pi); A.pi_ship.push(fH.pi);
  }
  acc[rate] = A;
  const mn = a => a.reduce((s, v) => s + v, 0) / a.length;
  console.log(`rate ${(100 * rate).toFixed(0).padStart(2)}%  n=${carefulIdx.length + k}  |  ` +
    `AUC sim ${mn(A.auc_sim).toFixed(3)} ship ${mn(A.auc_ship).toFixed(3)}  |  ` +
    `MCC  sim@z${Zsim} ${mn(A.mcc_sim).toFixed(3)}  sim@z2.5 ${mn(A.mcc_sim_z25).toFixed(3)}  ` +
    `ship@z2.5 ${mn(A.mcc_ship).toFixed(3)}  ship@z${Zsim} ${mn(A.mcc_ship_zsim).toFixed(3)}  |  ` +
    `pi-hat sim ${(100 * mn(A.pi_sim)).toFixed(1)}% ship ${(100 * mn(A.pi_ship)).toFixed(1)}%  |  ` +
    `skew(attentive) sim ${mn(A.skew_sim).toFixed(2)} ship ${mn(A.skew_ship).toFixed(2)}`);
}

const mn = a => a.reduce((s, v) => s + v, 0) / a.length;
const sd = a => { const m = mn(a); return Math.sqrt(mn(a.map(v => (v - m) ** 2))); };
const out = [["rate", "n", "k_careless", "reps",
  "auc_sim", "auc_sim_sd", "auc_shipped", "auc_shipped_sd", "d_auc", "d_auc_lo", "d_auc_hi",
  "mcc_sim", "mcc_shipped", "d_mcc", "mcc_sim_z25", "mcc_shipped_zsim",
  "prec_sim", "rec_sim", "prec_shipped", "rec_shipped",
  "pi_sim", "pi_shipped", "skew_att_sim", "skew_att_shipped", "auc_rr", "auc_pt"]];
for (const rate of RATES) {
  const A = acc[rate], k = Math.round(carefulIdx.length * rate / (1 - rate));
  const dA = A.auc_sim.map((v, i) => v - A.auc_ship[i]);
  const s = [...dA].sort((a, b) => a - b);
  out.push([rate, carefulIdx.length + k, k, REPS,
    mn(A.auc_sim).toFixed(4), sd(A.auc_sim).toFixed(4), mn(A.auc_ship).toFixed(4), sd(A.auc_ship).toFixed(4),
    mn(dA).toFixed(4), s[Math.floor(0.025 * s.length)].toFixed(4), s[Math.floor(0.975 * s.length)].toFixed(4),
    mn(A.mcc_sim).toFixed(4), mn(A.mcc_ship).toFixed(4), (mn(A.mcc_sim) - mn(A.mcc_ship)).toFixed(4),
    mn(A.mcc_sim_z25).toFixed(4), mn(A.mcc_ship_zsim).toFixed(4),
    mn(A.prec_sim).toFixed(4), mn(A.rec_sim).toFixed(4), mn(A.prec_ship).toFixed(4), mn(A.rec_ship).toFixed(4),
    mn(A.pi_sim).toFixed(4), mn(A.pi_ship).toFixed(4),
    mn(A.skew_sim).toFixed(4), mn(A.skew_ship).toFixed(4), mn(A.auc_rr).toFixed(4), mn(A.auc_pt).toFixed(4)]);
}
fs.writeFileSync("simcal_rate_sweep.csv", out.map(r => r.join(",")).join("\n"));
console.log("\nwrote simcal_rate_sweep.csv");
