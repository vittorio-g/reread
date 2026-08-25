/* score_waves.js — within-person validation study on COVID-Dynamic 16 waves.
 * For each wave w1..w16 in Dataset/gt_benchmark_candidates/covid_dynamic/waves/:
 * load w{w}_matrix.csv (Likert items, header row, empty cells = missing) and
 * w{w}_meta.csv (pid,flag,nvalid; row-aligned), score with the CaReReRe engine
 * (ensembleFeatures, seed=1, iterations=150), take O.rr (careless-oriented,
 * robust-standardized within sample) and the gated eta, and pool everything
 * into waves/panel.csv (pid,wave,rc,eta,flag,nvalid). */
const R = require("./site/reread.js");
const fs = require("fs");
const path = require("path");

const DIR = path.join(__dirname, "..", "Dataset", "gt_benchmark_candidates", "covid_dynamic", "waves");
const W = R.WEIGHTS;
const out = ["pid,wave,rc,eta,flag,nvalid"];

for (let w = 1; w <= 16; w++) {
  const t0 = Date.now();
  const P = R.parseCSV(fs.readFileSync(path.join(DIR, "w" + w + "_matrix.csv"), "utf8"));
  const J = P.header.length;
  const rows = P.rows.filter(r => r.length > 1 || (r.length === 1 && r[0] !== ""));
  const n = rows.length;
  const mat = new Float64Array(n * J);
  for (let i = 0; i < n; i++)
    for (let j = 0; j < J; j++) {
      const s = rows[i][j];
      // empty cell = missing -> NaN (Number("") would be 0, which is wrong)
      const x = (s == null || s === "") ? NaN : Number(s);
      mat[i * J + j] = Number.isFinite(x) ? x : NaN;
    }

  const M = R.parseCSV(fs.readFileSync(path.join(DIR, "w" + w + "_meta.csv"), "utf8"));
  const meta = M.rows.filter(r => r.length >= 3 && r[0] !== "");
  if (meta.length !== n) throw new Error("w" + w + ": matrix n=" + n + " != meta n=" + meta.length);

  // d2 carries zero ensemble weight; skip its costly eigendecomposition.
  const ef = R.ensembleFeatures(mat, n, J, { seed: 1, iterations: 150, skipD2: true });
  const O = ef.oriented;
  const g = Math.max(0, Math.min(1, (ef.profileInfo - 0.01) / 0.05));

  for (let i = 0; i < n; i++) {
    const eta = W.b0 + W.rr * O.rr[i] +
      g * (W.longstring * O.longstring[i] + W.person_total * O.person_total[i]);
    out.push([meta[i][0], w, O.rr[i].toFixed(6), eta.toFixed(6), meta[i][1], meta[i][2]].join(","));
  }
  console.log("w" + String(w).padStart(2) + "  n=" + n + "  J=" + J +
    "  gate=" + g.toFixed(3) + "  profileInfo=" + ef.profileInfo.toFixed(4) +
    "  [" + ((Date.now() - t0) / 1000).toFixed(1) + "s]");
}

fs.writeFileSync(path.join(DIR, "panel.csv"), out.join("\n") + "\n");
console.log("\nwrote " + path.join(DIR, "panel.csv") + "  (" + (out.length - 1) + " person-waves)");
