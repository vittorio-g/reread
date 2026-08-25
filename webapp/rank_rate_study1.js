/* rank_rate_study1.js — Study 1: is the RANKING (Algorithm 1, the rc index) robust as the
 * careless share grows? Cutoff-free: AUC only.
 * Why this can break: rc selects the coupled pairs from the SAMPLE correlation matrix and
 * Person-Total uses the SAMPLE mean profile, so at high careless rates the careless
 * contaminate the very reference they are judged against. The engine is therefore
 * RECOMPUTED on every mixture (never reuse full-sample scores).
 * Pools: 70 careful, 87 careless -> rates up to ~90% by thinning the careful side. */
const fs = require("fs");
const R = require("./site/reread.js");

const M = R.parseCSV(fs.readFileSync("_study_matrix.csv", "utf8"));
const y = R.parseCSV(fs.readFileSync("_study_labels.csv", "utf8")).rows.map(r => Number(r[0]));
const J = M.header.length;
const rows = M.rows.map(r => r.map(v => R.toNum(v)));
const carefulIdx = y.map((v, i) => v === 0 ? i : -1).filter(i => i >= 0);
const carelessIdx = y.map((v, i) => v === 1 ? i : -1).filter(i => i >= 0);

const rng = R._rng(20260728);
function sampleK(arr, k) { const a = arr.slice(); for (let i = 0; i < k; i++) { const j = i + Math.floor(rng() * (a.length - i)); const t = a[i]; a[i] = a[j]; a[j] = t; } return a.slice(0, k); }
function auc(sc, lab) {
  const idx = [...sc.keys()].filter(i => Number.isFinite(sc[i])).sort((a, b) => sc[a] - sc[b]);
  const rk = []; idx.forEach((id, r) => rk[id] = r + 1);
  let n1 = 0, n0 = 0, s = 0;
  for (const i of idx) { if (lab[i] === 1) { n1++; s += rk[i]; } else n0++; }
  return n1 && n0 ? (s - n1 * (n1 + 1) / 2) / (n1 * n0) : NaN;
}
const mean = v => v.reduce((s, x) => s + x, 0) / v.length;
const sd = v => { const m = mean(v); return Math.sqrt(mean(v.map(x => (x - m) ** 2))); };

const REPS = 30;
console.log("Study 1 ranking robustness (AUC only, engine recomputed per mixture, " + REPS + " reps/rate)");
console.log("rate  N    (nCareful/nCareless) | AUC rc         | AUC ensemble   | AUC LongString | AUC PersonTotal");
const out = [["true_rate", "N", "n_careful", "n_careless", "auc_rc", "sd_rc", "auc_ens", "sd_ens", "auc_ls", "auc_pt"]];

for (const target of [0.05, 0.10, 0.20, 0.30, 0.40, 0.50, 0.60, 0.70, 0.80, 0.90]) {
  // keep the larger side whole, thin the other, so N stays as large as possible
  let nCareful, nCareless;
  if (target <= carelessIdx.length / (carefulIdx.length + carelessIdx.length)) {
    nCareful = carefulIdx.length; nCareless = Math.round(nCareful * target / (1 - target));
  } else {
    nCareless = carelessIdx.length; nCareful = Math.round(nCareless * (1 - target) / target);
  }
  if (nCareless < 3 || nCareful < 3) continue;
  const arc = [], aen = [], als = [], apt = [];
  for (let rep = 0; rep < REPS; rep++) {
    const ids = sampleK(carefulIdx, nCareful).concat(sampleK(carelessIdx, nCareless));
    const n = ids.length, m = new Float64Array(n * J);
    for (let i = 0; i < n; i++) for (let j = 0; j < J; j++) m[i * J + j] = rows[ids[i]][j];
    const lab = ids.map(i => y[i]);
    const res = R.ensemble(m, n, J, { iterations: 200, seed: 1000 + rep });
    arc.push(auc(res.rr.map(v => -v), lab));          // low rc = careless -> negate
    aen.push(auc(res.p, lab));
    als.push(auc(res.oriented.longstring, lab));
    apt.push(auc(res.oriented.person_total, lab));
  }
  const N = nCareful + nCareless, tr = nCareless / N;
  console.log(
    (100 * tr).toFixed(0).padStart(3) + "%  " + String(N).padStart(3) +
    "  (" + String(nCareful).padStart(2) + "/" + String(nCareless).padStart(2) + ")        | " +
    mean(arc).toFixed(3) + " (" + sd(arc).toFixed(3) + ") | " +
    mean(aen).toFixed(3) + " (" + sd(aen).toFixed(3) + ") | " +
    mean(als).toFixed(3) + "          | " + mean(apt).toFixed(3));
  out.push([tr.toFixed(3), N, nCareful, nCareless, mean(arc).toFixed(4), sd(arc).toFixed(4),
            mean(aen).toFixed(4), sd(aen).toFixed(4), mean(als).toFixed(4), mean(apt).toFixed(4)]);
}
fs.writeFileSync("rank_rate_study1.csv", out.map(r => r.join(",")).join("\n"));
console.log("\nwrote rank_rate_study1.csv");
