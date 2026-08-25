/* eval_simcal_study1.js — STAGE B of the simulation-only calibration experiment.
 *
 * Takes the combiner calibrated purely on synthetic data (sim_calib.json, from
 * sim_calibrate.js), freezes it, and applies it to the real collected sample
 * (Study 1) that the fit has never seen. Reports ranking and flagging quality
 * against comparators that DID see real data, so the price of never touching it
 * is visible:
 *
 *   sim-only      weights and z from simulation alone            <- the experiment
 *   shipped       weights fitted on all of Study 1               (in-sample, optimistic)
 *   LOO-refit     weights refitted leaving each respondent out   (honest real-data fit)
 *   equal (1,1,1) unweighted z-sum, no fitting at all
 *   integer (3,1,1) the rounded shape the robustness analysis found sufficient
 *   single indices rr / LongString / Person-Total alone
 *
 *   node eval_simcal_study1.js       -> simcal_study1.csv + report on stdout
 */
const fs = require("fs");
const R = require("./site/reread.js");

const CAL = JSON.parse(fs.readFileSync("sim_calib.json", "utf8"));
const Wsim = CAL.weights, Zsim = CAL.zBest;
const Wship = R.WEIGHTS;

/* ---- Study 1 features (same engine, same defaults, computed once) ---- */
const M = R.parseCSV(fs.readFileSync("_study_matrix.csv", "utf8"));
const n = M.rows.length, J = M.header.length;
const mat = new Float64Array(n * J);
for (let i = 0; i < n; i++) for (let j = 0; j < J; j++) mat[i * J + j] = R.toNum(M.rows[i][j]);
const lab = R.parseCSV(fs.readFileSync("_study_labels.csv", "utf8"));
const y = lab.rows.map(r => Number(r[0]));
const grp = lab.rows.map(r => r[1].replace(/"/g, ""));

const ef = R.ensembleFeatures(mat, n, J, { seed: 1, iterations: 400 });
const O = ef.oriented;
const g = Math.max(0, Math.min(1, (ef.profileInfo - 0.01) / 0.05));
const X = Array.from({ length: n }, (_, i) => [1, O.rr[i], g * O.longstring[i], g * O.person_total[i]]);

/* ---- metrics ---- */
const auc = (s, yy) => {
  const po = [], ne = [];
  for (let i = 0; i < yy.length; i++) (yy[i] ? po : ne).push(s[i]);
  let c = 0; for (const a of po) for (const b of ne) c += a > b ? 1 : a === b ? 0.5 : 0;
  return c / (po.length * ne.length);
};
function conf(flag) {
  let tp = 0, fp = 0, fn = 0, tn = 0;
  for (let i = 0; i < y.length; i++) { if (y[i] === 1) { flag[i] ? tp++ : fn++; } else { flag[i] ? fp++ : tn++; } }
  const d = Math.sqrt((tp + fp) * (tp + fn) * (tn + fp) * (tn + fn));
  return { tp, fp, tn, fn, mcc: d ? (tp * tn - fp * fn) / d : 0, prec: tp + fp ? tp / (tp + fp) : 0, rec: tp / (tp + fn) };
}
const oracleMCC = s => {
  const srt = [...s].sort((a, b) => a - b); let best = -1;
  for (const thr of srt) { const c = conf(s.map(v => v >= thr)); if (c.mcc > best) best = c.mcc; }
  return best;
};
const spearman = (a, b) => {
  const rk = v => { const idx = v.map((x, i) => [x, i]).sort((p, q) => p[0] - q[0]); const r = new Array(v.length);
    let i = 0; while (i < idx.length) { let j = i; while (j + 1 < idx.length && idx[j + 1][0] === idx[i][0]) j++;
      const m = (i + j) / 2 + 1; for (let k = i; k <= j; k++) r[idx[k][1]] = m; i = j + 1; } return r; };
  const ra = rk(a), rb = rk(b), m = ra.length;
  const ma = ra.reduce((s, v) => s + v, 0) / m, mb = rb.reduce((s, v) => s + v, 0) / m;
  let num = 0, da = 0, db = 0;
  for (let i = 0; i < m; i++) { num += (ra[i] - ma) * (rb[i] - mb); da += (ra[i] - ma) ** 2; db += (rb[i] - mb) ** 2; }
  return num / Math.sqrt(da * db);
};

/* ---- ridge logistic, identical to the one that produced the shipped weights ---- */
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
function fit(rows, yy) {
  let b = new Array(P).fill(0);
  for (let it = 0; it < 60; it++) {
    const gr = new Array(P).fill(0), H = Array.from({ length: P }, () => new Array(P).fill(0));
    for (let i = 0; i < rows.length; i++) {
      let e = 0; for (let k = 0; k < P; k++) e += b[k] * rows[i][k];
      const p = 1 / (1 + Math.exp(-e)), w = Math.max(p * (1 - p), 1e-6);
      for (let k = 0; k < P; k++) { gr[k] += (yy[i] - p) * rows[i][k]; for (let l = 0; l < P; l++) H[k][l] += w * rows[i][k] * rows[i][l]; }
    }
    for (let k = 1; k < P; k++) { gr[k] -= 1.0 * b[k]; H[k][k] += 1.0; }
    const s = solve(H, gr); let mx = 0;
    for (let k = 0; k < P; k++) { b[k] += s[k]; mx = Math.max(mx, Math.abs(s[k])); }
    if (mx < 1e-9) break;
  }
  return b;
}
const etaOf = W => X.map(r => W.b0 + W.rr * r[1] + W.longstring * r[2] + W.person_total * r[3]);

/* LOO on Study 1: the honest real-data fit */
const etaLOO = new Array(n);
for (let i = 0; i < n; i++) {
  const Xi = X.filter((_, k) => k !== i), yi = y.filter((_, k) => k !== i);
  const b = fit(Xi, yi);
  etaLOO[i] = b[0] + b[1] * X[i][1] + b[2] * X[i][2] + b[3] * X[i][3];
}

const MODELS = [
  ["sim-only (weights + z from simulation)", etaOf(Wsim), Zsim],
  ["shipped (fit on all of Study 1)", etaOf(Wship), 2.5],
  ["LOO-refit on Study 1", etaLOO, 2.5],
  ["equal weights (1, 1, 1)", X.map(r => r[1] + r[2] + r[3]), 2.5],
  ["integer weights (3, 1, 1)", X.map(r => 3 * r[1] + r[2] + r[3]), 2.5],
  ["rr alone", Array.from(O.rr), 2.5],
  ["LongString alone", Array.from(O.longstring), 2.5],
  ["Person-Total alone", Array.from(O.person_total), 2.5]
];

console.log(`Study 1: n=${n}, J=${J}, careless=${y.reduce((s, v) => s + v, 0)} (${(100 * y.reduce((s, v) => s + v, 0) / n).toFixed(1)}%), gate g=${g.toFixed(3)}`);
console.log(`simulation corpus: ${CAL.corpus.datasets} datasets, ${CAL.corpus.respondents} respondents, variants ${CAL.corpus.variants.join("/")}`);
console.log(`weights  sim-only : ${JSON.stringify(Wsim)}   z=${Zsim}`);
console.log(`weights  shipped  : ${JSON.stringify(Wship)}   z=2.5\n`);
console.log("model                                     AUC   oracleMCC |  autoflag: n   MCC   prec   rec   pi-hat");
console.log("-".repeat(103));

const out = [["model", "auc", "oracle_mcc", "z", "n_flagged", "mcc", "precision", "recall", "pi_hat"]];
const detail = {};
for (const [name, eta, z] of MODELS) {
  const a = auc(eta, y), om = oracleMCC(eta);
  const fl = R.autoFlag(eta, z), c = conf(fl.flagged);
  const nf = fl.flagged.filter(Boolean).length;
  console.log(`${name.padEnd(41)} ${a.toFixed(3)}   ${om.toFixed(3)}     |  ` +
    `${String(nf).padStart(3)}  ${c.mcc.toFixed(3)}  ${c.prec.toFixed(3)}  ${c.rec.toFixed(3)}  ${(100 * fl.pi).toFixed(1)}%`);
  out.push([name, a.toFixed(4), om.toFixed(4), z, nf, c.mcc.toFixed(4), c.prec.toFixed(4), c.rec.toFixed(4), fl.pi.toFixed(4)]);
  detail[name] = { eta, flagged: fl.flagged };
}

const eSim = detail["sim-only (weights + z from simulation)"].eta;
const eShip = detail["shipped (fit on all of Study 1)"].eta;
console.log(`\nagreement sim-only vs shipped : Spearman ${spearman(eSim, eShip).toFixed(4)}   ` +
  `same flags on ${detail["sim-only (weights + z from simulation)"].flagged
    .filter((v, i) => v === detail["shipped (fit on all of Study 1)"].flagged[i]).length}/${n} respondents`);

/* paired bootstrap over respondents: is the AUC gap real, or resampling noise?
 * Both scores are recomputed on the same resample, so the pairing is preserved. */
const brng = R._rng(20260805);
function bootDeltaAUC(a, b, B) {
  const d = [];
  for (let r = 0; r < B; r++) {
    const idx = Array.from({ length: n }, () => Math.floor(brng() * n));
    const ya = idx.map(i => y[i]);
    if (ya.every(v => v === ya[0])) { r--; continue; }
    d.push(auc(idx.map(i => a[i]), ya) - auc(idx.map(i => b[i]), ya));
  }
  d.sort((p, q) => p - q);
  return { mean: d.reduce((s, v) => s + v, 0) / d.length, lo: d[Math.floor(0.025 * B)], hi: d[Math.floor(0.975 * B)] };
}
const B = 2000;
for (const [other, name] of [[etaLOO, "LOO-refit"], [eShip, "shipped (in-sample)"]]) {
  const c = bootDeltaAUC(eSim, other, B);
  console.log(`  delta AUC  sim-only - ${name.padEnd(20)} = ${c.mean >= 0 ? "+" : ""}${c.mean.toFixed(4)}` +
    `  95% CI [${c.lo.toFixed(4)}, ${c.hi.toFixed(4)}]${(c.lo <= 0 && c.hi >= 0) ? "   (includes zero)" : ""}`);
}

console.log("\nsim-only flags by condition (the collected ground truth):");
const by = {};
grp.forEach((gr_, i) => {
  by[gr_] = by[gr_] || { n: 0, f: 0 };
  by[gr_].n++; if (detail["sim-only (weights + z from simulation)"].flagged[i]) by[gr_].f++;
});
for (const k of Object.keys(by).sort())
  console.log(`   ${k.padEnd(14)} ${String(by[k].f).padStart(3)}/${String(by[k].n).padStart(3)} (${(100 * by[k].f / by[k].n).toFixed(0)}%)`);

/* also: how far off would the sim-only cut be if the shipped z were used instead */
const flSimAtShipped = R.autoFlag(eSim, 2.5), cShipZ = conf(flSimAtShipped.flagged);
console.log(`\nsim-only weights with the shipped z=2.5 : ${flSimAtShipped.flagged.filter(Boolean).length} flagged, ` +
  `MCC ${cShipZ.mcc.toFixed(3)}, precision ${cShipZ.prec.toFixed(3)}, recall ${cShipZ.rec.toFixed(3)}`);
out.push(["sim-only weights @ shipped z=2.5", auc(eSim, y).toFixed(4), oracleMCC(eSim).toFixed(4), 2.5,
  flSimAtShipped.flagged.filter(Boolean).length, cShipZ.mcc.toFixed(4), cShipZ.prec.toFixed(4), cShipZ.rec.toFixed(4), flSimAtShipped.pi.toFixed(4)]);

fs.writeFileSync("simcal_study1.csv", out.map(r => r.join(",")).join("\n"));
console.log("\nwrote simcal_study1.csv");
