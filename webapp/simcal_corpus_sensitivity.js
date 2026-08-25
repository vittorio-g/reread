/* simcal_corpus_sensitivity.js — STAGE D of the simulation-only calibration
 * experiment: how much of the transfer depends on how the simulation was built?
 *
 * Refits the combiner from sim_calib_features.csv under restrictions of the
 * corpus — one structural variant at a time, one careless rate at a time, one
 * length band at a time — and applies each frozen result to Study 1. If the
 * transfer only works for one particular corpus design, that is a weakness of
 * the claim; if it survives every restriction, the claim is about the method.
 *
 *   node simcal_corpus_sensitivity.js   -> simcal_corpus_sensitivity.csv
 */
const fs = require("fs");
const R = require("./site/reread.js");

/* ---- simulated corpus, respondent level ---- */
const raw = fs.readFileSync("sim_calib_features.csv", "utf8").split(/\r?\n/).filter(Boolean);
const head = raw[0].split(",");
const col = {}; head.forEach((h, i) => col[h] = i);
const corpus = raw.slice(1).map(l => {
  const f = l.split(",");
  return {
    variant: f[col.variant], J: +f[col.J], n: +f[col.n], rate: +f[col.rate], gate: +f[col.gate],
    x: [1, +f[col.rr], +f[col.gate] * +f[col.ls], +f[col.gate] * +f[col.pt]], y: +f[col.label]
  };
});

/* ---- Study 1 ---- */
const M = R.parseCSV(fs.readFileSync("_study_matrix.csv", "utf8"));
const n1 = M.rows.length, J1 = M.header.length;
const mat = new Float64Array(n1 * J1);
for (let i = 0; i < n1; i++) for (let j = 0; j < J1; j++) mat[i * J1 + j] = R.toNum(M.rows[i][j]);
const y1 = R.parseCSV(fs.readFileSync("_study_labels.csv", "utf8")).rows.map(r => Number(r[0]));
const ef = R.ensembleFeatures(mat, n1, J1, { seed: 1, iterations: 400 });
const O = ef.oriented, g1 = Math.max(0, Math.min(1, (ef.profileInfo - 0.01) / 0.05));
const X1 = Array.from({ length: n1 }, (_, i) => [1, O.rr[i], g1 * O.longstring[i], g1 * O.person_total[i]]);

const auc = (s, yy) => {
  const po = [], ne = [];
  for (let i = 0; i < yy.length; i++) (yy[i] ? po : ne).push(s[i]);
  let c = 0; for (const a of po) for (const b of ne) c += a > b ? 1 : a === b ? 0.5 : 0;
  return c / (po.length * ne.length);
};
const conf = flag => {
  let tp = 0, fp = 0, fn = 0, tn = 0;
  for (let i = 0; i < y1.length; i++) { if (y1[i] === 1) { flag[i] ? tp++ : fn++; } else { flag[i] ? fp++ : tn++; } }
  const d = Math.sqrt((tp + fp) * (tp + fn) * (tn + fp) * (tn + fn));
  return { mcc: d ? (tp * tn - fp * fn) / d : 0, prec: tp + fp ? tp / (tp + fp) : 0, rec: tp / (tp + fn), n: tp + fp };
};

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
function fit(rows) {
  let b = new Array(P).fill(0);
  for (let it = 0; it < 60; it++) {
    const gr = new Array(P).fill(0), H = Array.from({ length: P }, () => new Array(P).fill(0));
    for (const rw of rows) {
      let e = 0; for (let k = 0; k < P; k++) e += b[k] * rw.x[k];
      const p = 1 / (1 + Math.exp(-e)), w = Math.max(p * (1 - p), 1e-6);
      for (let k = 0; k < P; k++) { gr[k] += (rw.y - p) * rw.x[k]; for (let l = 0; l < P; l++) H[k][l] += w * rw.x[k] * rw.x[l]; }
    }
    for (let k = 1; k < P; k++) { gr[k] -= 1.0 * b[k]; H[k][k] += 1.0; }
    const s = solve(H, gr); let mx = 0;
    for (let k = 0; k < P; k++) { b[k] += s[k]; mx = Math.max(mx, Math.abs(s[k])); }
    if (mx < 1e-9) break;
  }
  return b;
}

const RESTRICTIONS = [["full corpus", () => true]];
for (const v of [...new Set(corpus.map(r => r.variant))]) RESTRICTIONS.push([`variant = ${v}`, r => r.variant === v]);
for (const rt of [...new Set(corpus.map(r => r.rate))].sort((a, b) => a - b)) RESTRICTIONS.push([`rate = ${(100 * rt).toFixed(0)}% only`, r => r.rate === rt]);
RESTRICTIONS.push(["short only (J <= 60)", r => r.J <= 60]);
RESTRICTIONS.push(["long only (J >= 120)", r => r.J >= 120]);
RESTRICTIONS.push(["small n only (n = 150)", r => r.n === 150]);

console.log("corpus restriction              rows     b0     w_rr   w_ls   w_pt |  Study 1: AUC   MCC   prec   rec   flagged");
console.log("-".repeat(112));
const out = [["restriction", "rows", "b0", "w_rr", "w_ls", "w_pt", "auc", "mcc", "precision", "recall", "n_flagged"]];
const Z = JSON.parse(fs.readFileSync("sim_calib.json", "utf8")).zBest;
for (const [name, keep] of RESTRICTIONS) {
  const rows = corpus.filter(keep);
  if (rows.length < 500 || !rows.some(r => r.y === 1) || !rows.some(r => r.y === 0)) continue;
  const b = fit(rows);
  const eta = X1.map(x => b.reduce((s, v, k) => s + v * x[k], 0));
  const a = auc(eta, y1), c = conf(R.autoFlag(eta, Z).flagged);
  console.log(`${name.padEnd(30)} ${String(rows.length).padStart(5)}  ` +
    `${b[0].toFixed(2).padStart(6)} ${b[1].toFixed(2).padStart(6)} ${b[2].toFixed(2).padStart(6)} ${b[3].toFixed(2).padStart(6)} |  ` +
    `   ${a.toFixed(3)} ${c.mcc.toFixed(3)} ${c.prec.toFixed(3)} ${c.rec.toFixed(3)}   ${c.n}`);
  out.push([name, rows.length, b[0].toFixed(4), b[1].toFixed(4), b[2].toFixed(4), b[3].toFixed(4),
    a.toFixed(4), c.mcc.toFixed(4), c.prec.toFixed(4), c.rec.toFixed(4), c.n]);
}
fs.writeFileSync("simcal_corpus_sensitivity.csv", out.map(r => r.join(",")).join("\n"));
console.log("\nwrote simcal_corpus_sensitivity.csv");
