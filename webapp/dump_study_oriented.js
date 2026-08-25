/* dump_study_oriented.js — Study 1 oriented features in EXACTLY the same form as
 * the external dump (ext_bench_oriented): robust-z, higher = more careless.
 * Output: _study_oriented.csv with y,rr,longstring,person_total  */
const fs = require("fs");
const R = require("./site/reread.js");

const M = R.parseCSV(fs.readFileSync("_study_matrix.csv", "utf8"));
const L = R.parseCSV(fs.readFileSync("_study_labels.csv", "utf8"));
const n = M.rows.length, J = M.header.length;
const mat = new Float64Array(n * J);
for (let i = 0; i < n; i++) for (let j = 0; j < J; j++) mat[i * J + j] = R.toNum(M.rows[i][j]);

// label column: first column of _study_labels.csv ("careless")
const y = L.rows.map(r => Number(r[0]));
if (y.length !== n) { console.error("row mismatch", y.length, n); process.exit(1); }

const ef = R.ensembleFeatures(mat, n, J, { iterations: 400, seed: 1 });
const O = ef.oriented;
const lines = ["y,rr,longstring,person_total"];
for (let i = 0; i < n; i++)
  lines.push([y[i], O.rr[i].toFixed(6), O.longstring[i].toFixed(6), O.person_total[i].toFixed(6)].join(","));
fs.writeFileSync("_study_oriented.csv", lines.join("\n"));
console.log("wrote _study_oriented.csv  n=" + n + "  careless=" + y.reduce((s, v) => s + v, 0) +
            "  profileInfo=" + ef.profileInfo.toFixed(3));
