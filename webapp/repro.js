/* repro.js — reproducibility & idempotence of ReReRe (tool engine).
 *  A. Determinism: same data + same seed -> identical z.
 *  B. Monte-Carlo stability: across many seeds, how stable are z and the
 *     flag decision? Reported per iteration count -> recommended default.
 *  C. Idempotence: clean a dataset (remove flagged), reload it, and check
 *     whether the algorithm flags NEW respondents. Iterate to convergence.
 *
 * Run:  node repro.js
 */
const fs = require("fs");
const R = require("./site/rerere.js");

const M = R.parseCSV(fs.readFileSync("_study_matrix.csv", "utf8"));
const L = R.parseCSV(fs.readFileSync("_study_labels.csv", "utf8"));
const J = M.header.length;
const allRows = M.rows.map(r => r.map(Number));
const y = L.rows.map(r => Number(r[0]));          // 1 = careless (ground truth)
const N = allRows.length;
const ZTHR = 1.5;

function mat(rowsIdx) {
  const n = rowsIdx.length, m = new Float64Array(n * J);
  for (let i = 0; i < n; i++) for (let j = 0; j < J; j++) m[i * J + j] = allRows[rowsIdx[i]][j];
  return { m, n };
}
function runIdx(rowsIdx, opts) {
  const { m, n } = mat(rowsIdx);
  return R.run(m, n, J, Object.assign({ iterations: 200, seed: 1 }, opts));
}
const mean = v => v.reduce((s, x) => s + x, 0) / v.length;
const sd = v => { const m = mean(v); return Math.sqrt(mean(v.map(x => (x - m) ** 2))); };

/* ---------- A. determinism ---------- */
console.log("== A. Determinism (same data, same seed) ==");
{
  const a = runIdx([...Array(N).keys()], { seed: 7 }).z;
  const b = runIdx([...Array(N).keys()], { seed: 7 }).z;
  let maxDiff = 0;
  for (let i = 0; i < a.length; i++) maxDiff = Math.max(maxDiff, Math.abs(a[i] - b[i]));
  console.log("  max |z1 - z2| over " + N + " respondents = " + maxDiff.toExponential(2) +
              (maxDiff === 0 ? "   -> BIT-IDENTICAL" : "   -> NOT identical!"));
}

/* ---------- B. Monte-Carlo stability across seeds ---------- */
console.log("\n== B. Monte-Carlo stability across 40 random seeds (full data) ==");
{
  const idx = [...Array(N).keys()];
  for (const iters of [50, 100, 200, 500]) {
    const K = 40;
    const zRuns = [];
    const flagRuns = [];
    for (let s = 0; s < K; s++) {
      const z = runIdx(idx, { iterations: iters, seed: 1000 + s }).z;
      zRuns.push(z);
      flagRuns.push(z.map(v => v <= ZTHR));
    }
    // per-respondent z SD across seeds
    const zsd = [];
    for (let i = 0; i < N; i++) zsd.push(sd(zRuns.map(z => z[i])));
    // flag-flip: respondents whose flag is not unanimous across seeds
    let flips = 0;
    for (let i = 0; i < N; i++) {
      const votes = flagRuns.map(f => f[i]);
      if (votes.some(v => v) && votes.some(v => !v)) flips++;
    }
    console.log("  iterations=" + String(iters).padEnd(4) +
      " | mean per-resp z SD = " + mean(zsd).toFixed(3) +
      " | flag-unstable respondents = " + flips + "/" + N +
      " (" + Math.round(100 * flips / N) + "%)");
  }
  console.log("  (flag-unstable = flagged in some seeds but not others at z<=" + ZTHR + ")");
}

/* ---------- C. idempotence ---------- */
console.log("\n== C. Idempotence: clean, reload, re-flag? ==");
{
  let keep = [...Array(N).keys()];
  console.log("  round 0: n=" + keep.length);
  for (let round = 1; round <= 5; round++) {
    const res = runIdx(keep, { iterations: 200, seed: 42 });
    const flaggedLocal = [];
    keep.forEach((gi, li) => { if (res.z[li] <= ZTHR) flaggedLocal.push(gi); });
    const nCarefulFlagged = flaggedLocal.filter(gi => y[gi] === 0).length;
    console.log("  round " + round + ": on n=" + keep.length +
      " -> flags " + flaggedLocal.length +
      " (" + (flaggedLocal.length - nCarefulFlagged) + " true-careless, " +
      nCarefulFlagged + " careful) | diagnostic=" + res.diagnostic.level);
    if (flaggedLocal.length === 0) { console.log("  -> CONVERGED: a cleaned set produces no new flags."); break; }
    keep = keep.filter(gi => flaggedLocal.indexOf(gi) < 0);
  }
}

/* ---------- C2. false-positive cascade on a GROUND-TRUTH-clean set ---------- */
console.log("\n== C2. Reload a set of ONLY verified-careful respondents ==");
{
  const carefulOnly = [...Array(N).keys()].filter(i => y[i] === 0);
  const res = runIdx(carefulOnly, { iterations: 300, seed: 9 });
  const flagged = res.z.filter(v => v <= ZTHR).length;
  console.log("  " + carefulOnly.length + " verified-careful respondents, no careless present.");
  console.log("  ReReRe flags " + flagged + " of them (" +
    Math.round(100 * flagged / carefulOnly.length) + "%) at z<=" + ZTHR +
    " | diagnostic=" + res.diagnostic.level);
  console.log("  (this is the false-positive floor: how many clean people a cleaned dataset still flags)");
}
