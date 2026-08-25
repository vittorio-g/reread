/* eval_tisp.js — run the shipped CaReReRe ensemble on the TISP dataset and evaluate against
 * the attention-check ground truth.
 *   GT : ATTCHECK_RES != 1  ("please select 'strongly disagree'")  -> careless
 * Matrix = the 57 Likert items administered BEFORE the check (so that failers and passers have
 * equally complete data; scoring later blocks would leak the survey's routing pattern).
 * Reports AUC (ensemble, rc, Person-Total, LongString) and MCC under the deployed automatic
 * rule (Standard and High) plus the oracle best threshold.
 */
const R = require("./site/rerere.js"); const W = R.WEIGHTS; const fs = require("fs");
const DIR = "../Dataset/gt_benchmark_candidates/tisp/";
const P = R.parseCSV(fs.readFileSync(DIR + "_matrix.csv", "utf8"));
const J = P.header.length;
const rows = P.rows.map(r => r.map(v => { const x = Number(v); return Number.isFinite(x) ? x : NaN; }));
const lab = R.parseCSV(fs.readFileSync(DIR + "_labels.csv", "utf8")).rows.map(r => Number(r[0]));
const n = Math.min(rows.length, lab.length);
console.log("n=" + n + "  J=" + J + "  careless=" + lab.slice(0, n).reduce((s, v) => s + v, 0) +
            " (" + (100 * lab.slice(0, n).reduce((s, v) => s + v, 0) / n).toFixed(1) + "%)");

const mat = new Float64Array(n * J);
for (let i = 0; i < n; i++) for (let j = 0; j < J; j++) mat[i * J + j] = rows[i][j];

console.log("scoring ...");
const t0 = Date.now();
const ef = R.ensembleFeatures(mat, n, J, { seed: 1, iterations: 200 });
const O = ef.oriented;
const g = Math.max(0, Math.min(1, (ef.profileInfo - 0.01) / 0.05));
const eta = new Array(n);
for (let i = 0; i < n; i++)
  eta[i] = W.b0 + W.rr * O.rr[i] + g * (W.longstring * O.longstring[i] + W.person_total * O.person_total[i]);
console.log("scored in " + ((Date.now() - t0) / 1000).toFixed(0) + "s  | profileInfo=" +
            ef.profileInfo.toFixed(4) + "  gate g=" + g.toFixed(3));

function auc(score, y) {
  const idx = score.map((s, i) => [s, y[i]]).sort((a, b) => a[0] - b[0]);
  let r = 0, n1 = 0, n0 = 0, i = 0;
  while (i < idx.length) {
    let j = i; while (j + 1 < idx.length && idx[j + 1][0] === idx[i][0]) j++;
    const rank = (i + j) / 2 + 1;
    for (let t = i; t <= j; t++) { if (idx[t][1] === 1) { r += rank; n1++; } else n0++; }
    i = j + 1;
  }
  return n1 && n0 ? (r - n1 * (n1 + 1) / 2) / (n1 * n0) : NaN;
}
function mcc(flag, y) {
  let tp = 0, tn = 0, fp = 0, fn = 0;
  for (let i = 0; i < y.length; i++) { if (flag[i]) { y[i] ? tp++ : fp++; } else { y[i] ? fn++ : tn++; } }
  const d = Math.sqrt((tp + fp) * (tp + fn) * (tn + fp) * (tn + fn));
  return { mcc: d ? (tp * tn - fp * fn) / d : 0, prec: (tp + fp) ? tp / (tp + fp) : 0,
           rec: (tp + fn) ? tp / (tp + fn) : 0, frac: (tp + fp) / y.length };
}
const y = lab.slice(0, n);
console.log("\n=== AUC (threshold-free) ===");
for (const [nm, s] of [["ensemble", eta], ["rc", Array.from(O.rr)],
                       ["Person-Total", Array.from(O.person_total)],
                       ["LongString", Array.from(O.longstring)]])
  console.log("  " + nm.padEnd(14) + auc(s, y).toFixed(3));

console.log("\n=== MCC ===");
for (const lvl of ["standard", "high"]) {
  const af = R.autoFlag(eta, lvl);
  const m = mcc(af.flagged, y);
  console.log("  auto " + lvl.padEnd(9) + "MCC=" + m.mcc.toFixed(3) + "  prec=" + m.prec.toFixed(2) +
              "  rec=" + m.rec.toFixed(2) + "  flagged=" + (100 * m.frac).toFixed(1) + "%" +
              "  gateFired=" + af.twoComp);
}
const srt = [...eta].sort((a, b) => a - b);
let best = { mcc: -1 };
for (let q = 2; q <= 98; q += 1) {
  const thr = srt[Math.floor(q / 100 * (n - 1))];
  const m = mcc(eta.map(v => v > thr), y);
  if (m.mcc > best.mcc) best = { ...m, q };
}
console.log("  oracle       MCC=" + best.mcc.toFixed(3) + "  prec=" + best.prec.toFixed(2) +
            "  rec=" + best.rec.toFixed(2) + "  at top " + (100 - best.q) + "%");
fs.writeFileSync(DIR + "_scores.csv",
  "eta,rc,longstring,person_total,y\n" +
  Array.from({ length: n }, (_, i) => [eta[i], O.rr[i], O.longstring[i], O.person_total[i], y[i]].join(",")).join("\n"));
console.log("\nwrote _scores.csv");
