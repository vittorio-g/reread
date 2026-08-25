/* dump_feats_ext.js — dump ALL oriented features (including the two that carry zero weight in the
 * shipped model, IRV and Mahalanobis D^2) for every external dataset, so we can test whether adding
 * them to the ensemble helps. Writes <dir>/_feats.csv with y1,y2 where available. */
const R = require("./site/reread.js"); const fs = require("fs");
const D = "../Dataset/gt_benchmark_candidates";
const SETS = [
  ["TISP",      D + "/tisp",               ["_labels.csv"]],
  ["TUMI",      D + "/tumi",               ["_labels1.csv", "_labels2.csv"]],
  ["warnIPIP",  D + "/warning_ipipneo300", ["_labels1.csv", "_labels2.csv"]],
  ["smarvus",   D + "/smarvus",            ["_labels1.csv", "_labels2.csv"]],
  ["KayS2",     D + "/kay_idris_idria/s2", ["_labels1.csv", "_labels2.csv"]],
  ["KayS1",     D + "/kay_idris_idria/s1", ["_labels1.csv", "_labels2.csv"]],
];
for (const [name, dir, labs] of SETS) {
  const P = R.parseCSV(fs.readFileSync(dir + "/_matrix.csv", "utf8"));
  const J = P.header.length;
  const rows = P.rows.map(r => r.map(v => { const x = Number(v); return Number.isFinite(x) ? x : NaN; }));
  const n = rows.length;
  const mat = new Float64Array(n * J);
  for (let i = 0; i < n; i++) for (let j = 0; j < J; j++) mat[i * J + j] = rows[i][j];
  const ef = R.ensembleFeatures(mat, n, J, { seed: 1, iterations: 200 });
  const O = ef.oriented;
  const ys = labs.map(l => R.parseCSV(fs.readFileSync(dir + "/" + l, "utf8")).rows.map(r => Number(r[0])));
  const head = ["rc", "longstring", "person_total", "irv", "d2"].concat(ys.map((_, i) => "y" + (i + 1)));
  const out = [head.join(",")];
  for (let i = 0; i < n; i++)
    out.push([O.rr[i], O.longstring[i], O.person_total[i], O.irv[i], O.d2[i]]
      .concat(ys.map(y => y[i])).join(","));
  fs.writeFileSync(dir + "/_feats.csv", out.join("\n"));
  console.log(name.padEnd(10) + " n=" + String(n).padStart(5) + " J=" + String(J).padStart(4) +
              " -> _feats.csv (" + ys.length + " label set(s))");
}
