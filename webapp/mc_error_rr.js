/* mc_error_rr.js — Monte-Carlo error of the rr z-score itself (reviewer #4).
 * For the study data, vary the number of permutations and the seed; measure
 *  (a) per-respondent SD of z across seeds at fixed iters (the MC error),
 *  (b) AUC(rr) stability across seeds,
 *  (c) how many respondents hit the degenerate (zero-variance) rule,
 *  (d) the number of coupled pairs k (is the 15-pair floor active?).
 */
const fs = require("fs");
const R = require("./site/rerere.js");

const M = R.parseCSV(fs.readFileSync("_study_matrix.csv", "utf8"));
const y = R.parseCSV(fs.readFileSync("_study_labels.csv", "utf8")).rows.map(r => Number(r[0]));
const n = M.rows.length, J = M.header.length;
const mat = new Float64Array(n * J);
for (let i = 0; i < n; i++) for (let j = 0; j < J; j++) mat[i * J + j] = R.toNum(M.rows[i][j]);

function auc(sc, yy) {
  const idx = [...sc.keys()].sort((a, b) => sc[a] - sc[b]); const rk = [];
  idx.forEach((id, r) => rk[id] = r + 1);
  let n1 = 0, n0 = 0, s = 0;
  for (let i = 0; i < sc.length; i++) { if (yy[i] === 1) { n1++; s += rk[i]; } else n0++; }
  return (s - n1 * (n1 + 1) / 2) / (n1 * n0);
}
const mean = v => v.reduce((s, x) => s + x, 0) / v.length;
const sd = v => { const m = mean(v); return Math.sqrt(mean(v.map(x => (x - m) ** 2))); };
const quantile = (v, q) => { const s = [...v].sort((a, b) => a - b); return s[Math.min(s.length - 1, Math.floor(q * s.length))]; };

const SEEDS = 30;
const base = R.run(mat, n, J, { iterations: 200, seed: 1 });
const k = base.params.k, nPairs = base.params.nPairs;
console.log(`Study: n=${n}, J=${J} items, careless=${y.reduce((s,v)=>s+v,0)}/${n}`);
console.log(`Coupled pairs k=${k} of ${nPairs} (corProp 0.03); 15-pair floor active? ${k <= 15 ? "YES" : "no"}\n`);

console.log("iters | mean per-resp SD(z) | 95th pct SD(z) | AUC(rr) mean | AUC(rr) SD across seeds");
for (const iters of [100, 200, 400, 800]) {
  const Zruns = [];      // Zruns[s][i]
  const aucs = [];
  for (let s = 0; s < SEEDS; s++) {
    const z = R.run(mat, n, J, { iterations: iters, seed: 1000 + s }).z;
    Zruns.push(z);
    aucs.push(auc(z.map(v => -v), y));   // low z = careless -> negate for AUC
  }
  const perRespSD = [];
  for (let i = 0; i < n; i++) perRespSD.push(sd(Zruns.map(z => z[i])));
  console.log(
    String(iters).padEnd(6) + "| " +
    mean(perRespSD).toFixed(4).padEnd(20) + "| " +
    quantile(perRespSD, 0.95).toFixed(4).padEnd(15) + "| " +
    mean(aucs).toFixed(4).padEnd(13) + "| " +
    sd(aucs).toFixed(4)
  );
}

// degenerate (zero-variance) respondents: how many, and what z they receive
const zc = base.z, minZ = Math.min(...zc);
let atFloor = 0; for (let i = 0; i < n; i++) if (Math.abs(zc[i] - minZ) < 1e-9) atFloor++;
console.log(`\nRespondents pinned to the most-careless z by the zero-variance rule (approx): ${atFloor}`);
console.log(`(These are pure straight-liners: C_obs and C_rand both 0 -> z undefined -> assigned min z=${minZ.toFixed(3)}.)`);
