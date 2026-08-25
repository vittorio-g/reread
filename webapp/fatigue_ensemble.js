/* fatigue_ensemble.js — FATIGUE TEST on the full CaReReRe ensemble (not the rc index alone).
 *
 * WHAT CAN BE COMPARED. The ensemble features are robust-z standardised WITHIN the scored data,
 * so mean(eta) is centred in each half by construction and carries no level information. The
 * ensemble-level quantities that DO compare across halves are:
 *   - pi_hat      : the calibration's estimated careless rate (one-sided fit of the attentive mode)
 *   - flag rate   : the share the shipped automatic cut actually flags
 *   - flag asymmetry: respondents flagged in the 2nd half but not the 1st, vs the reverse.
 *                   Fatigue predicts an EXCESS of clean->flagged over flagged->clean.
 *   - AUC vs GT   : does the ensemble detect the archived careless better late than early?
 *
 * CONTROL. Halves differ in content, not only position, so each dataset is also split at random
 * (same sizes) a few times; the random split is the null for every quantity above.
 */
const fs = require("fs");
const path = require("path");
const R = require("./site/reread.js");

const ITERS = 150, NMAX = 3000, RAND = 5;
const rng = R._rng(20260806);
const EXT = "ext_bench2";
const idx = R.parseCSV(fs.readFileSync(path.join(EXT, "index.csv"), "utf8"));
const NAME = idx.header.indexOf("name"), JJ = idx.header.indexOf("J");

const mean = v => v.length ? v.reduce((s, x) => s + x, 0) / v.length : NaN;
function auc(sc, lab) {
  const ok = [...sc.keys()].filter(i => Number.isFinite(sc[i]) && (lab[i] === 0 || lab[i] === 1));
  const ord = ok.slice().sort((a, b) => sc[a] - sc[b]); const rk = [];
  ord.forEach((id, r) => rk[id] = r + 1);
  let n1 = 0, n0 = 0, s = 0;
  for (const i of ok) { if (lab[i] === 1) { n1++; s += rk[i]; } else n0++; }
  return n1 && n0 ? (s - n1 * (n1 + 1) / 2) / (n1 * n0) : NaN;
}
function ens(rows, cols) {
  const n = rows.length, J = cols.length;
  const m = new Float64Array(n * J);
  for (let i = 0; i < n; i++) for (let j = 0; j < J; j++) m[i * J + j] = rows[i][cols[j]];
  return R.ensemble(m, n, J, { iterations: ITERS, seed: 1 });
}

console.log("FATIGUE TEST — CaReReRe ENSEMBLE, first vs second half");
console.log("(mean eta is not comparable across halves: features are robust-z within the scored data)\n");
console.log("dataset          J   n    | pi1st  pi2nd | flag1st flag2nd | ->flag <-flag  asym | AUC1st AUC2nd | rand asym");
const out = [["dataset","J","n","pi_first","pi_second","flag_first","flag_second","to_flag","from_flag","asym","auc_first","auc_second","rand_asym_mean"]];

for (const r of idx.rows) {
  const tag = r[NAME], J = Number(r[JJ]);
  if (J < 60) continue;
  const mf = path.join(EXT, tag + "_M.csv"), yf = path.join(EXT, tag + "_y.csv");
  if (!fs.existsSync(mf)) { console.log("  (missing) " + tag); continue; }
  const P = R.parseCSV(fs.readFileSync(mf, "utf8"));
  let rows = [P.header].concat(P.rows).map(rr => rr.map(v => R.toNum(v)));
  const Y = R.parseCSV(fs.readFileSync(yf, "utf8"));
  const yi = Y.header.indexOf("y");
  let lab = Y.rows.map(rr => Number(rr[yi]));
  if (lab.length !== rows.length) { console.log("  (row mismatch) " + tag); continue; }
  if (rows.length > NMAX) {
    const step = Math.ceil(rows.length / NMAX);
    const keep = [...Array(rows.length).keys()].filter(i => i % step === 0);
    rows = keep.map(i => rows[i]); lab = keep.map(i => lab[i]);
  }
  const half = Math.floor(J / 2);
  const first = [...Array(half).keys()], second = [...Array(J - half).keys()].map(i => i + half);
  const A = ens(rows, first), B = ens(rows, second);

  let toFlag = 0, fromFlag = 0;
  for (let i = 0; i < rows.length; i++) {
    if (!A.flagged[i] && B.flagged[i]) toFlag++;
    else if (A.flagged[i] && !B.flagged[i]) fromFlag++;
  }
  const asym = toFlag - fromFlag;

  // random-split null for the asymmetry
  const rAsym = [];
  for (let k = 0; k < RAND; k++) {
    const perm = [...Array(J).keys()];
    for (let i = J - 1; i > 0; i--) { const j = Math.floor(rng() * (i + 1)); [perm[i], perm[j]] = [perm[j], perm[i]]; }
    const X = ens(rows, perm.slice(0, half)), Z = ens(rows, perm.slice(half));
    let t = 0, f = 0;
    for (let i = 0; i < rows.length; i++) { if (!X.flagged[i] && Z.flagged[i]) t++; else if (X.flagged[i] && !Z.flagged[i]) f++; }
    rAsym.push(t - f);
  }
  const ra = mean(rAsym);

  console.log(
    tag.padEnd(15) + String(J).padStart(4) + String(rows.length).padStart(5) + "  | " +
    (100 * A.estimatedRate).toFixed(1).padStart(5) + "% " + (100 * B.estimatedRate).toFixed(1).padStart(5) + "% | " +
    (100 * A.nFlagged / rows.length).toFixed(1).padStart(6) + "% " + (100 * B.nFlagged / rows.length).toFixed(1).padStart(6) + "% | " +
    String(toFlag).padStart(5) + String(fromFlag).padStart(7) + String(asym).padStart(6) + " | " +
    auc(A.p, lab).toFixed(3) + "  " + auc(B.p, lab).toFixed(3) + " | " + (ra >= 0 ? "+" : "") + ra.toFixed(1));
  out.push([tag, J, rows.length, A.estimatedRate.toFixed(4), B.estimatedRate.toFixed(4),
            (A.nFlagged / rows.length).toFixed(4), (B.nFlagged / rows.length).toFixed(4),
            toFlag, fromFlag, asym, auc(A.p, lab).toFixed(4), auc(B.p, lab).toFixed(4), ra.toFixed(2)]);
}
fs.writeFileSync("fatigue_ensemble.csv", out.map(r => r.join(",")).join("\n"));
console.log("\nasym = (flagged only in 2nd half) - (flagged only in 1st). Fatigue predicts asym > 0, above the random-split null.");
console.log("wrote fatigue_ensemble.csv");
