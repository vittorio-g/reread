/* envelope.js — operating envelope of ReReRe on the collected study data.
 * Subsamples J items × n respondents (stratified careful/careless) from the
 * 84×109 ground-truth dataset, runs the tool engine, and records the AUC of
 * z_rr for separating careful vs careless. Maps where the method holds up.
 *
 * Run:  node envelope.js
 */
const fs = require("fs");
const R = require("./site/rerere.js");

const M = R.parseCSV(fs.readFileSync("_study_matrix.csv", "utf8"));
const L = R.parseCSV(fs.readFileSync("_study_labels.csv", "utf8"));
const items = M.header;                       // 109 item names
const rows = M.rows.map(r => r.map(Number));  // numeric matrix rows
const y = L.rows.map(r => Number(r[0]));      // 1 = careless
const N0 = rows.length, Jall = items.length;

const carefulIdx = y.map((v, i) => v === 0 ? i : -1).filter(i => i >= 0);
const carelessIdx = y.map((v, i) => v === 1 ? i : -1).filter(i => i >= 0);

const rng = R._rng(12345);
function sampleK(arr, k) {                     // k without replacement
  const a = arr.slice();
  for (let i = 0; i < k; i++) {
    const j = i + Math.floor(rng() * (a.length - i));
    const t = a[i]; a[i] = a[j]; a[j] = t;
  }
  return a.slice(0, k);
}
function aucLowPos(scores, labels) {           // low score => positive(careless)
  const arr = scores.map((v, i) => ({ v: -v, y: labels[i] })).sort((a, b) => a.v - b.v);
  let r = 0, np = 0, nn = 0;
  arr.forEach((o, i) => { if (o.y) { r += i + 1; np++; } else nn++; });
  return (np && nn) ? (r - np * (np + 1) / 2) / (np * nn) : NaN;
}

const J_GRID = [10, 20, 30, 45, 60, 80, 109];
const N_GRID = [20, 30, 40, 50, 60, 70, 76];   // total, stratified 50/50
const REPS = 30, ITERS = 100;

function subrun(J, nTot) {
  const half = Math.floor(nTot / 2);
  if (half > carefulIdx.length || half > carelessIdx.length) return null;
  const cols = sampleK([...Array(Jall).keys()], J);
  const pick = [...sampleK(carefulIdx, half), ...sampleK(carelessIdx, nTot - half)];
  const n = pick.length;
  const sub = new Float64Array(n * J);
  const lab = new Array(n);
  for (let i = 0; i < n; i++) {
    lab[i] = y[pick[i]];
    for (let j = 0; j < J; j++) sub[i * J + j] = rows[pick[i]][cols[j]];
  }
  try {
    const res = R.run(sub, n, J, { iterations: ITERS, seed: 1 + Math.floor(rng() * 1e6) });
    return { auc: aucLowPos(res.z, lab), level: res.diagnostic.level };
  } catch (e) { return { auc: NaN, level: "err" }; }
}

console.log("ReReRe operating envelope — AUC (careful vs careless), mean of " + REPS + " subsamples");
console.log("dataset: " + N0 + " respondents (" + carefulIdx.length + " careful / " +
            carelessIdx.length + " careless), " + Jall + " items\n");

const grid = {};
for (const nTot of N_GRID) {
  for (const J of J_GRID) {
    const aucs = [];
    let okStruct = 0;
    for (let r = 0; r < REPS; r++) {
      const res = subrun(J, nTot);
      if (res && Number.isFinite(res.auc)) { aucs.push(res.auc); if (res.level === "ok" || res.level === "marginal") okStruct++; }
    }
    aucs.sort((a, b) => a - b);
    grid[nTot + "_" + J] = {
      mean: aucs.reduce((s, v) => s + v, 0) / aucs.length,
      p10: aucs[Math.floor(aucs.length * 0.1)] || NaN,
      pGood: aucs.filter(a => a >= 0.8).length / aucs.length
    };
  }
}

// AUC table
function fmtTable(field, label) {
  console.log("\n=== " + label + " ===");
  let head = "n \\ J".padEnd(8);
  for (const J of J_GRID) head += ("J=" + J).padStart(8);
  console.log(head);
  for (const nTot of N_GRID) {
    let line = ("n=" + nTot).padEnd(8);
    for (const J of J_GRID) {
      const v = grid[nTot + "_" + J][field];
      line += (Number.isFinite(v) ? v.toFixed(2) : "  -").padStart(8);
    }
    console.log(line);
  }
}
fmtTable("mean", "Mean AUC");
fmtTable("pGood", "Fraction of subsamples with AUC >= 0.80 (reliability)");

// frontier: smallest J reaching AUC 0.8 / 0.7 at each n
console.log("\n=== Breakdown frontier (mean AUC) ===");
for (const thr of [0.8, 0.7]) {
  console.log("  AUC >= " + thr + ":");
  for (const nTot of N_GRID) {
    const js = J_GRID.filter(J => grid[nTot + "_" + J].mean >= thr);
    console.log("    n=" + String(nTot).padEnd(3) + " -> " +
      (js.length ? "needs J >= " + js[0] + " items" : "not reached at any J"));
  }
}
