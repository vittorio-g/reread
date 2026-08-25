/* diag_pooling.js — is pooling many countries/sites a congenital problem for these datasets?
 * If between-country differences dominate, the sample correlation matrix (which drives pair
 * selection) and the sample mean profile (which drives Person-Total) describe GROUP differences
 * rather than within-person coherence, and eta becomes multimodal for geographic reasons.
 * Test: score TISP pooled vs Germany-only, and decompose eta's variance by country.
 */
const R = require("./site/reread.js"); const W = R.WEIGHTS; const fs = require("fs");
const DIR = "../Dataset/gt_benchmark_candidates/tisp/";
function score(sfx) {
  const P = R.parseCSV(fs.readFileSync(DIR + "_matrix" + sfx + ".csv", "utf8"));
  const J = P.header.length;
  const rows = P.rows.map(r => r.map(v => { const x = Number(v); return Number.isFinite(x) ? x : NaN; }));
  const y = R.parseCSV(fs.readFileSync(DIR + "_labels" + sfx + ".csv", "utf8")).rows.map(r => Number(r[0]));
  const cty = R.parseCSV(fs.readFileSync(DIR + "_country" + sfx + ".csv", "utf8")).rows.map(r => r[0]);
  const n = Math.min(rows.length, y.length);
  const mat = new Float64Array(n * J);
  for (let i = 0; i < n; i++) for (let j = 0; j < J; j++) mat[i * J + j] = rows[i][j];
  const ef = R.ensembleFeatures(mat, n, J, { seed: 1, iterations: 200 });
  const O = ef.oriented, g = Math.max(0, Math.min(1, (ef.profileInfo - 0.01) / 0.05));
  const eta = new Array(n);
  for (let i = 0; i < n; i++)
    eta[i] = W.b0 + W.rr * O.rr[i] + g * (W.longstring * O.longstring[i] + W.person_total * O.person_total[i]);
  return { eta, y: y.slice(0, n), cty: cty.slice(0, n), O, n, J };
}
function auc(s, y) {
  const idx = s.map((v, i) => [v, y[i]]).sort((a, b) => a[0] - b[0]);
  let r = 0, n1 = 0, n0 = 0, i = 0;
  while (i < idx.length) { let j = i; while (j + 1 < idx.length && idx[j + 1][0] === idx[i][0]) j++;
    const rk = (i + j) / 2 + 1;
    for (let t = i; t <= j; t++) { if (idx[t][1] === 1) { r += rk; n1++; } else n0++; } i = j + 1; }
  return n1 && n0 ? (r - n1 * (n1 + 1) / 2) / (n1 * n0) : NaN;
}
/* share of eta variance that is BETWEEN countries (eta-squared) */
function etaSq(x, grp) {
  const by = new Map();
  for (let i = 0; i < x.length; i++) { if (!by.has(grp[i])) by.set(grp[i], []); by.get(grp[i]).push(x[i]); }
  const gm = x.reduce((a, b) => a + b, 0) / x.length;
  let ssb = 0, sst = 0;
  for (const v of x) sst += (v - gm) * (v - gm);
  for (const [, arr] of by) { const m = arr.reduce((a, b) => a + b, 0) / arr.length; ssb += arr.length * (m - gm) * (m - gm); }
  return { r2: ssb / sst, groups: by.size };
}
for (const [tag, sfx] of [["POOLED (68 countries)", ""], ["GERMANY only", "_DEU"]]) {
  const S = score(sfx);
  const af = R.autoFlag(S.eta, "standard");
  const prev = 100 * S.y.reduce((a, b) => a + b, 0) / S.n;
  console.log("\n=== " + tag + " ===");
  console.log("  n=" + S.n + " J=" + S.J + " careless=" + prev.toFixed(1) + "%");
  console.log("  AUC  ensemble=" + auc(S.eta, S.y).toFixed(3) +
              "  rc=" + auc(Array.from(S.O.rr), S.y).toFixed(3) +
              "  PT=" + auc(Array.from(S.O.person_total), S.y).toFixed(3) +
              "  LS=" + auc(Array.from(S.O.longstring), S.y).toFixed(3));
  console.log("  gate=" + (af.twoComp ? "ON" : "OFF") + "  flagged=" + (100 * af.estRate).toFixed(1) +
              "%  piHat=" + (100 * af.pi).toFixed(1) + "%");
  if (sfx === "") {
    const e = etaSq(S.eta, S.cty);
    console.log("  >> between-country share of eta variance: " + (100 * e.r2).toFixed(1) + "%  (" + e.groups + " countries)");
    const rr = etaSq(Array.from(S.O.rr), S.cty), pt = etaSq(Array.from(S.O.person_total), S.cty);
    console.log("  >> same for rc: " + (100 * rr.r2).toFixed(1) + "%   for Person-Total: " + (100 * pt.r2).toFixed(1) + "%");
    /* careless rate by country - does prevalence itself vary by country? */
    const by = new Map();
    for (let i = 0; i < S.n; i++) { if (!by.has(S.cty[i])) by.set(S.cty[i], []); by.get(S.cty[i]).push(S.y[i]); }
    const rates = [...by.entries()].filter(([, a]) => a.length >= 100)
      .map(([c, a]) => [c, a.reduce((x, z) => x + z, 0) / a.length, a.length])
      .sort((a, b) => a[1] - b[1]);
    console.log("  >> careless rate by country (>=100 respondents), lowest and highest:");
    for (const [c, r, m] of rates.slice(0, 3)) console.log("       " + c + " " + (100 * r).toFixed(1) + "%  (n=" + m + ")");
    console.log("       ...");
    for (const [c, r, m] of rates.slice(-3)) console.log("       " + c + " " + (100 * r).toFixed(1) + "%  (n=" + m + ")");
  }
}
