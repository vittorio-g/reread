/* ls_presentation_order.js -- Study 1 LongString: column order vs the recorded
 * presentation order.
 *
 * Items were presented in a per-participant randomized order (item_order is
 * archived), while LongString -- like every published implementation -- is
 * computed on the canonical column order of the exported matrix. This scores the
 * shipped ensemble both ways to check whether that choice matters.
 *
 *   node ls_presentation_order.js   (expects ls_pres_order.csv: one permutation
 *                                    per study row, canonical indices, csv)
 */
const fs = require("fs");
const R = require("./site/rerere.js");

function readCsv(p) {
  const L = fs.readFileSync(p, "utf8").trim().split(/\r?\n/);
  return L.map(l => l.split(","));
}
const mat = readCsv("_study_matrix.csv"), head = mat[0];
const n = mat.length - 1, J = head.length;
const M = new Float64Array(n * J);
for (let i = 0; i < n; i++) for (let j = 0; j < J; j++) {
  const v = parseFloat(mat[i + 1][j]);
  M[i * J + j] = Number.isFinite(v) ? v : NaN;
}
const y = readCsv("_study_labels.csv").slice(1).map(r => +r[0]);
const perm = readCsv("ls_pres_order.csv").slice(1).map(r => r.map(Number));
if (perm.length !== n) throw new Error("permutation rows != study rows");

const ef = R.ensembleFeatures(M, n, J, { seed: 1, iterations: 400 });
const W = R.WEIGHTS;
const g = Math.max(0, Math.min(1, (ef.profileInfo - 0.01) / 0.05));  // as in ensemble()

// LongString recomputed on each respondent's own presentation order, on the same
// proportion-rescaled, median-imputed values the engine uses.
const colMax = new Float64Array(J), colMed = new Float64Array(J);
for (let j = 0; j < J; j++) {
  const v = [];
  for (let i = 0; i < n; i++) if (Number.isFinite(M[i * J + j])) v.push(M[i * J + j]);
  v.sort((a, b) => a - b);
  colMax[j] = v[v.length - 1]; colMed[j] = v[v.length >> 1];
}
const lsPres = new Float64Array(n);
for (let i = 0; i < n; i++) {
  const seq = perm[i].map(j => {
    const v = M[i * J + j];
    return (Number.isFinite(v) ? v : colMed[j]) / colMax[j];
  });
  let best = 1, run = 1;
  for (let k = 1; k < seq.length; k++) {
    if (seq[k] === seq[k - 1]) { run++; if (run > best) best = run; } else run = 1;
  }
  lsPres[i] = best;
}

// robust z, matching the engine's standardization of the LongString feature
function robustZ(a) {
  const s = [...a].sort((x, y) => x - y), med = s[s.length >> 1];
  const d = a.map(v => Math.abs(v - med)).sort((x, y) => x - y);
  const mad = 1.4826 * d[d.length >> 1] || 1;
  return Array.from(a, v => (v - med) / mad);
}
const zLsPres = robustZ(Array.from(lsPres));
const O = ef.oriented;

function auc(sc, lab) {
  const p = [], q = [];
  sc.forEach((s, i) => (lab[i] ? p : q).push(s));
  let t = 0;
  for (const a of p) for (const b of q) t += a > b ? 1 : a === b ? 0.5 : 0;
  return t / (p.length * q.length);
}
const etaShipped = [], etaPres = [], etaNoRcS = [], etaNoRcP = [];
for (let i = 0; i < n; i++) {
  const rc = W.rr * O.rr[i], pt = W.person_total * O.person_total[i];
  etaShipped.push(W.b0 + rc + g * (W.longstring * O.longstring[i] + pt));
  etaPres.push(W.b0 + rc + g * (W.longstring * zLsPres[i] + pt));
  etaNoRcS.push(W.b0 + g * (W.longstring * O.longstring[i] + pt));
  etaNoRcP.push(W.b0 + g * (W.longstring * zLsPres[i] + pt));
}
const r = (a, b) => {
  const m = x => x.reduce((s, v) => s + v, 0) / x.length, ma = m(a), mb = m(b);
  let sab = 0, sa = 0, sb = 0;
  for (let i = 0; i < a.length; i++) {
    const da = a[i] - ma, db = b[i] - mb; sab += da * db; sa += da * da; sb += db * db;
  }
  return sab / Math.sqrt(sa * sb);
};
console.log("gate g =", g.toFixed(3), " (1 = partners fully in)");
console.log("\nLongString alone      AUC   column order %s | presentation order %s",
  auc(Array.from(O.longstring), y).toFixed(3), auc(zLsPres, y).toFixed(3));
console.log("rc alone              AUC   %s", auc(Array.from(O.rr), y).toFixed(3));
console.log("Person-Total alone    AUC   %s", auc(Array.from(O.person_total), y).toFixed(3));
console.log("\nEnsemble (frozen w.)  AUC   column order %s | presentation order %s",
  auc(etaShipped, y).toFixed(3), auc(etaPres, y).toFixed(3));
console.log("Auxiliaries, no rc    AUC   column order %s | presentation order %s",
  auc(etaNoRcS, y).toFixed(3), auc(etaNoRcP, y).toFixed(3));
console.log("\nrc's AUC increment          column order %s | presentation order %s",
  (auc(etaShipped, y) - auc(etaNoRcS, y)).toFixed(3),
  (auc(etaPres, y) - auc(etaNoRcP, y)).toFixed(3));
console.log("corr(LS column, LS presentation) = %s", r(Array.from(O.longstring), zLsPres).toFixed(3));
