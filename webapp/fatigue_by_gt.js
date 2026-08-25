/* fatigue_by_gt.js — is the late-questionnaire decay a PERSON effect (fatigue) or a CONTENT
 * effect (the second half simply contains more homogeneous scales)?
 *
 * Decisive contrast: split respondents by the archived ground truth and compute the same
 * first-vs-second-half deltas within each group.
 *   - CONTENT effect  -> the delta is the same for attentive and careless respondents.
 *   - FATIGUE effect  -> the delta is larger (worse) for the careless / low-effort group,
 *     because disengagement is what grows over the session.
 * Metrics: rc (permutation z; lower = less coherent) and LongString (raw longest run;
 * higher = more straightlining).
 */
const fs = require("fs");
const path = require("path");
const R = require("./site/rerere.js");

const ITERS = 100, NMAX = 4000;
const EXT = "ext_bench2";
const idx = R.parseCSV(fs.readFileSync(path.join(EXT, "index.csv"), "utf8"));
const NAME = idx.header.indexOf("name"), JJ = idx.header.indexOf("J");

const mean = v => v.length ? v.reduce((s, x) => s + x, 0) / v.length : NaN;
const sd = v => { const m = mean(v); return Math.sqrt(v.reduce((s, x) => s + (x - m) * (x - m), 0) / Math.max(v.length - 1, 1)); };

function scoreCols(rows, cols) {
  const n = rows.length, J = cols.length;
  const m = new Float64Array(n * J);
  for (let i = 0; i < n; i++) for (let j = 0; j < J; j++) m[i * J + j] = rows[i][cols[j]];
  const res = R.run(m, n, J, { iterations: ITERS, seed: 1 });
  return { z: res.z, ls: res.longstring };
}

console.log("FATIGUE: person effect or content effect? First-vs-second-half deltas BY ground truth");
console.log("(rc delta negative = less coherent late; LS delta positive = more straightlining late)\n");
console.log("dataset          J   nAtt nCar | rc: dAtt   dCar   diff   | LS: dAtt   dCar   diff");
const out = [["dataset","J","n_att","n_car","rc_d_att","rc_d_car","rc_diff","ls_d_att","ls_d_car","ls_diff"]];

for (const r of idx.rows) {
  const tag = r[NAME], J = Number(r[JJ]);
  if (J < 60) continue;
  const mf = path.join(EXT, tag + "_M.csv"), yf = path.join(EXT, tag + "_y.csv");
  if (!fs.existsSync(mf) || !fs.existsSync(yf)) { console.log("  (missing) " + tag); continue; }

  const P = R.parseCSV(fs.readFileSync(mf, "utf8"));
  let rows = [P.header].concat(P.rows).map(rr => rr.map(v => R.toNum(v)));   // header is data
  const Y = R.parseCSV(fs.readFileSync(yf, "utf8"));
  const yi = Y.header.indexOf("y");
  let lab = Y.rows.map(rr => Number(rr[yi]));
  if (lab.length !== rows.length) { console.log("  (row mismatch) " + tag + " " + rows.length + " vs " + lab.length); continue; }

  if (rows.length > NMAX) {                       // keep it fast, preserve label balance
    const keep = [...Array(rows.length).keys()].filter((_, i) => i % Math.ceil(rows.length / NMAX) === 0);
    rows = keep.map(i => rows[i]); lab = keep.map(i => lab[i]);
  }
  const half = Math.floor(J / 2);
  const first = [...Array(half).keys()], second = [...Array(J - half).keys()].map(i => i + half);
  const A = scoreCols(rows, first), B = scoreCols(rows, second);

  const dz = [], dls = [];
  for (let i = 0; i < rows.length; i++) { dz.push(B.z[i] - A.z[i]); dls.push(B.ls[i] - A.ls[i]); }
  const att = i => lab[i] === 0, car = i => lab[i] === 1;
  const pick = (arr, f) => arr.filter((_, i) => f(i)).filter(Number.isFinite);
  const dzA = pick(dz, att), dzC = pick(dz, car), dlA = pick(dls, att), dlC = pick(dls, car);
  if (dzC.length < 10) { console.log(tag.padEnd(15) + "  (too few careless: " + dzC.length + ")"); continue; }

  const rcDiff = mean(dzC) - mean(dzA), lsDiff = mean(dlC) - mean(dlA);
  console.log(
    tag.padEnd(15) + String(J).padStart(4) + String(dzA.length).padStart(5) + String(dzC.length).padStart(5) + " | " +
    (mean(dzA) >= 0 ? "+" : "") + mean(dzA).toFixed(2).padStart(6) + " " +
    (mean(dzC) >= 0 ? "+" : "") + mean(dzC).toFixed(2).padStart(6) + " " +
    (rcDiff >= 0 ? "+" : "") + rcDiff.toFixed(2).padStart(6) + "  | " +
    (mean(dlA) >= 0 ? "+" : "") + mean(dlA).toFixed(2).padStart(6) + " " +
    (mean(dlC) >= 0 ? "+" : "") + mean(dlC).toFixed(2).padStart(6) + " " +
    (lsDiff >= 0 ? "+" : "") + lsDiff.toFixed(2).padStart(6));
  out.push([tag, J, dzA.length, dzC.length, mean(dzA).toFixed(4), mean(dzC).toFixed(4), rcDiff.toFixed(4),
            mean(dlA).toFixed(4), mean(dlC).toFixed(4), lsDiff.toFixed(4)]);
}
fs.writeFileSync("fatigue_by_gt.csv", out.map(r => r.join(",")).join("\n"));
console.log("\ndiff = careless minus attentive. Fatigue predicts rc diff NEGATIVE and LS diff POSITIVE.");
console.log("wrote fatigue_by_gt.csv");
