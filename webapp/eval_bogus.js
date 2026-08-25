/* eval_bogus.js — generic CaReReRe evaluation against bogus-item ground truth.
 * usage: node eval_bogus.js <dir> <labelfile> <label>
 * Scores <dir>/_matrix.csv once with the shipped ensemble and reports AUC (ensemble, rc,
 * Person-Total, LongString) plus MCC under the deployed automatic rule and at the oracle cut.
 */
const R = require("./site/reread.js"); const W = R.WEIGHTS; const fs = require("fs");
const DIR = process.argv[2], LAB = process.argv[3] || "_labels.csv", TAG = process.argv[4] || "";
const MAT = process.argv[5] || "_matrix.csv";
const P = R.parseCSV(fs.readFileSync(DIR + "/" + MAT, "utf8"));
const J = P.header.length;
const rows = P.rows.map(r => r.map(v => { const x = Number(v); return Number.isFinite(x) ? x : NaN; }));
const lab = R.parseCSV(fs.readFileSync(DIR + "/" + LAB, "utf8")).rows.map(r => Number(r[0]));
const n = Math.min(rows.length, lab.length), y = lab.slice(0, n);
const mat = new Float64Array(n * J);
for (let i = 0; i < n; i++) for (let j = 0; j < J; j++) mat[i * J + j] = rows[i][j];
const ef = R.ensembleFeatures(mat, n, J, { seed: 1, iterations: 200 });
const O = ef.oriented, g = Math.max(0, Math.min(1, (ef.profileInfo - 0.01) / 0.05));
const eta = new Array(n);
for (let i = 0; i < n; i++)
  eta[i] = W.b0 + W.rr * O.rr[i] + g * (W.longstring * O.longstring[i] + W.person_total * O.person_total[i]);
function auc(s, y) {
  const idx = s.map((v, i) => [v, y[i]]).sort((a, b) => a[0] - b[0]);
  let r = 0, n1 = 0, n0 = 0, i = 0;
  while (i < idx.length) { let j = i; while (j + 1 < idx.length && idx[j + 1][0] === idx[i][0]) j++;
    const rk = (i + j) / 2 + 1;
    for (let t = i; t <= j; t++) { if (idx[t][1] === 1) { r += rk; n1++; } else n0++; } i = j + 1; }
  return n1 && n0 ? (r - n1 * (n1 + 1) / 2) / (n1 * n0) : NaN;
}
function mcc(flag, y) { let tp=0,tn=0,fp=0,fn=0;
  for (let i = 0; i < y.length; i++) { if (flag[i]) { y[i] ? tp++ : fp++; } else { y[i] ? fn++ : tn++; } }
  const d = Math.sqrt((tp+fp)*(tp+fn)*(tn+fp)*(tn+fn));
  return { mcc: d ? (tp*tn-fp*fn)/d : 0, prec:(tp+fp)?tp/(tp+fp):0, rec:(tp+fn)?tp/(tp+fn):0, frac:(tp+fp)/y.length }; }
const af = R.autoFlag(eta, "standard"), m = mcc(af.flagged, y);
const srt = [...eta].sort((a, b) => a - b); let best = { mcc: -1, q: 0 };
for (let q = 1; q <= 99; q++) { const t = srt[Math.floor(q/100*(n-1))];
  const mm = mcc(eta.map(v => v > t), y); if (mm.mcc > best.mcc) best = { ...mm, q }; }
const prev = 100 * y.reduce((s, v) => s + v, 0) / n;
console.log([TAG.padEnd(22), "n=" + String(n).padStart(6), "J=" + String(J).padStart(4),
  "prev=" + prev.toFixed(1).padStart(5) + "%",
  "| AUC ens=" + auc(eta, y).toFixed(3), "rc=" + auc(Array.from(O.rr), y).toFixed(3),
  "PT=" + auc(Array.from(O.person_total), y).toFixed(3), "LS=" + auc(Array.from(O.longstring), y).toFixed(3),
  "| MCCauto=" + m.mcc.toFixed(3), "(flag " + (100*m.frac).toFixed(1) + "%, gate " + af.twoComp + ")",
  "MCCoracle=" + best.mcc.toFixed(3)].join("  "));
