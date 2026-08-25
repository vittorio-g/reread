/* study1_gateless.js — Study 1 recomputed with the gate removed, against the
 * gated engine on the SAME resamples.
 *
 * 70 attentive respondents are held fixed and real careless ones are drawn in to
 * hit a target rate; the engine is recomputed on every mixture (pair selection,
 * sign alignment and the robust-z reference all depend on the sample), so nothing
 * is carried over between rates. Both engines see identical data.
 *
 *   node study1_gateless.js   -> study1_gateless.csv
 */
const fs = require("fs");
const NEW = require("./site/reread.js");
const OLD = require("./_engine_gated.js");

const M = NEW.parseCSV(fs.readFileSync("_study_matrix.csv", "utf8"));
const J = M.header.length;
const y = NEW.parseCSV(fs.readFileSync("_study_labels.csv", "utf8")).rows.map(r => Number(r[0]));
const careful = [], careless = [];
y.forEach((v, i) => (v === 1 ? careless : careful).push(i));

const rng = NEW._rng(20260819);
const sampleK = (pool, k) => {
  const a = pool.slice();
  for (let i = 0; i < k; i++) { const j = i + Math.floor(rng() * (a.length - i)); const t = a[i]; a[i] = a[j]; a[j] = t; }
  return a.slice(0, k);
};
function conf(flag, yy) {
  let tp = 0, fp = 0, fn = 0, tn = 0;
  for (let i = 0; i < yy.length; i++) { if (yy[i] === 1) { flag[i] ? tp++ : fn++; } else { flag[i] ? fp++ : tn++; } }
  const d = Math.sqrt((tp + fp) * (tp + fn) * (tn + fp) * (tn + fn));
  return { mcc: d ? (tp * tn - fp * fn) / d : 0, prec: tp + fp ? tp / (tp + fp) : NaN,
           rec: tp / (tp + fn), share: (tp + fp) / yy.length };
}
const mn = a => { const v = a.filter(Number.isFinite); return v.length ? v.reduce((s, x) => s + x, 0) / v.length : NaN; };

const RATES = [0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.40, 0.50];
const REPS = 100;
const out = [["rate", "n", "engine", "mcc", "precision", "recall", "flag_share", "pi_hat", "silent_frac"]];
console.log("Study 1: gateless vs gated, engine recomputed on every mixture, " + REPS + " resamples per rate");
console.log("rate   n    engine     MCC   precision recall  %flagged  pi-hat   flags nobody");
for (const rate of RATES) {
  const k = Math.round(careful.length * rate / (1 - rate));
  const A = { gateless: { m: [], p: [], r: [], s: [], pi: [], zero: 0 }, gated: { m: [], p: [], r: [], s: [], pi: [], zero: 0 } };
  for (let rep = 0; rep < REPS; rep++) {
    const idx = careful.concat(sampleK(careless, k));
    const n = idx.length, mat = new Float64Array(n * J);
    for (let i = 0; i < n; i++) for (let j = 0; j < J; j++) mat[i * J + j] = NEW.toNum(M.rows[idx[i]][j]);
    const yy = idx.map(i => y[i]);
    for (const [tag, ENG] of [["gateless", NEW], ["gated", OLD]]) {
      const e = ENG.ensemble(mat, n, J, { seed: 1, iterations: 200, skipD2: true });
      const c = conf(e.flagged, yy), a = A[tag];
      a.m.push(c.mcc); a.p.push(c.prec); a.r.push(c.rec); a.s.push(c.share);
      a.pi.push(e.estimatedRate != null ? e.params.estimatedRate : NaN);
      if (e.nFlagged === 0) a.zero++;
    }
  }
  for (const tag of ["gateless", "gated"]) {
    const a = A[tag];
    console.log(`${(100 * rate).toFixed(0).padStart(3)}%  ${String(careful.length + k).padStart(3)}  ${tag.padEnd(9)} ` +
      `${mn(a.m).toFixed(3)}   ${mn(a.p).toFixed(3)}    ${mn(a.r).toFixed(3)}   ${(100 * mn(a.s)).toFixed(1)}%     ` +
      `${(100 * mn(a.pi)).toFixed(1)}%      ${a.zero}/${REPS}`);
    out.push([rate, careful.length + k, tag, mn(a.m).toFixed(4), mn(a.p).toFixed(4), mn(a.r).toFixed(4),
      mn(a.s).toFixed(4), mn(a.pi).toFixed(4), (a.zero / REPS).toFixed(3)]);
  }
}
/* the full collected sample, for the record: 55% careless, far outside the range
 * the cut is calibrated for, so it is reported and not used as a headline */
{
  const n = y.length, mat = new Float64Array(n * J);
  for (let i = 0; i < n; i++) for (let j = 0; j < J; j++) mat[i * J + j] = NEW.toNum(M.rows[i][j]);
  console.log("\nfull collected sample (55.4% careless, outside the calibrated range):");
  for (const [tag, ENG] of [["gateless", NEW], ["gated", OLD]]) {
    const e = ENG.ensemble(mat, n, J, { seed: 1, iterations: 400, skipD2: true });
    const c = conf(e.flagged, y);
    console.log(`  ${tag.padEnd(9)} MCC ${c.mcc.toFixed(3)}  precision ${c.prec.toFixed(3)}  recall ${c.rec.toFixed(3)}  ` +
      `flagged ${e.nFlagged} (${(100 * e.nFlagged / n).toFixed(0)}%)`);
    out.push(["full", n, tag, c.mcc.toFixed(4), c.prec.toFixed(4), c.rec.toFixed(4),
      (e.nFlagged / n).toFixed(4), "", ""]);
  }
}
fs.writeFileSync("study1_gateless.csv", out.map(r => r.join(",")).join("\n"));
console.log("\nwrote study1_gateless.csv");
