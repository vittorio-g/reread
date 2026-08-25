/* diag_gate.js — why does the BIC gate refuse to fire on some datasets?
 * The shipped acceptance test is:
 *    (bic1 - bic2 > 2) && pi > 0.01 && pi < 0.7 && separation > 1.0
 * Print each component per dataset so we can see which condition is the binding one. */
const R = require("./site/rerere.js"); const W = R.WEIGHTS; const fs = require("fs");
const D = "../Dataset/gt_benchmark_candidates";
const SETS = [
  ["TISP",        D + "/tisp",                 "_labels.csv"],
  ["TUMI",        D + "/tumi",                 "_labels1.csv"],
  ["warn IPIP",   D + "/warning_ipipneo300",   "_labels1.csv"],
  ["smarvus",     D + "/smarvus",              "_labels1.csv"],
  ["Kay S2",      D + "/kay_idris_idria/s2",   "_labels1.csv"],
  ["Kay S1",      D + "/kay_idris_idria/s1",   "_labels1.csv"],
];
console.log("dataset      n      J    | bic1-bic2   pi      sep    | gate  binding condition");
for (const [name, dir, lf] of SETS) {
  const P = R.parseCSV(fs.readFileSync(dir + "/_matrix.csv", "utf8"));
  const J = P.header.length;
  const rows = P.rows.map(r => r.map(v => { const x = Number(v); return Number.isFinite(x) ? x : NaN; }));
  const n = rows.length;
  const mat = new Float64Array(n * J);
  for (let i = 0; i < n; i++) for (let j = 0; j < J; j++) mat[i * J + j] = rows[i][j];
  const ef = R.ensembleFeatures(mat, n, J, { seed: 1, iterations: 200 });
  const O = ef.oriented, g = Math.max(0, Math.min(1, (ef.profileInfo - 0.01) / 0.05));
  const eta = new Array(n);
  for (let i = 0; i < n; i++)
    eta[i] = W.b0 + W.rr * O.rr[i] + g * (W.longstring * O.longstring[i] + W.person_total * O.person_total[i]);
  const gm = R.gaussMix1D(eta);
  const pooled = Math.sqrt((gm.v1 + gm.v2) / 2) || 1;
  const sep = (gm.m2 - gm.m1) / pooled;
  const dbic = gm.bic1 - gm.bic2;
  const cond = [["dBIC>2", dbic > 2], ["w>.01", gm.pi > 0.01], ["w<.7", gm.pi < 0.7], ["sep>1", sep > 1.0]];
  const ok = cond.every(c => c[1]);
  const fail = cond.filter(c => !c[1]).map(c => c[0]).join(",") || "-";
  const af = R.autoFlag(eta, "standard");     // the SHIPPED rule (one-sided cut + pi by subtraction)
  console.log([name.padEnd(11), String(n).padStart(5), String(J).padStart(5), "|",
    dbic.toFixed(1).padStart(9), gm.pi.toFixed(3).padStart(6), sep.toFixed(2).padStart(6), "|",
    (ok ? "ON " : "OFF"), fail.padEnd(8), "|",
    "piHat=" + (100 * af.pi).toFixed(1).padStart(5) + "%",
    "flagged=" + (100 * af.estRate).toFixed(1).padStart(5) + "%"].join(" "));
}
