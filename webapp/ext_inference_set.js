/* ext_inference_set.js -- external head-to-head recomputed on the behavioral-criterion
 * inference set.
 *
 * Two datasets carry a purely SELF-ASSESSED criterion (kay_s2: exclusion advice;
 * opsy_16pf: self-rated accuracy) -- a global judgement of one's own data quality,
 * promptable by straight-lining or speeding, i.e. patterns rc is blind to. This
 * recomputes every summary and every paired test on the seven datasets whose
 * criterion has a behavioral component, alongside the full nine, so the paper can
 * report both. Duckworth (discriminant control) is excluded from every average.
 *
 *   node ext_inference_set.js   -> ext_inference_set.csv
 */
const fs = require("fs");

const rows = fs.readFileSync("ext_deployable.csv", "utf8").trim().split(/\r?\n/).map(l => l.split(","));
const h = rows[0], R = rows.slice(1).map(r => Object.fromEntries(r.map((v, i) => [h[i], v])));

const CONTROL = ["duckworth"];             // discriminant negative control
const SELF_ASSESSED = ["kay_s2", "opsy_16pf"];
const all = [...new Set(R.map(r => r.dataset))];
const nine = all.filter(d => !CONTROL.includes(d));
const seven = nine.filter(d => !SELF_ASSESSED.includes(d));
const INDICES = ["ens", "rc", "PT", "IRV", "PsychSyn", "LongString", "Mahalanobis"];

const val = (d, idx, f) => {
  const r = R.find(x => x.dataset === d && x.index === idx);
  return r ? parseFloat(r[f]) : NaN;
};
const mean = a => a.reduce((s, v) => s + v, 0) / a.length;

/* exact one-sided Wilcoxon signed-rank: P(ensemble > competitor), enumerating all
   2^n sign assignments. Zero differences are dropped (Wilcoxon's convention). */
function exactSignedRank(diffs) {
  const d = diffs.filter(v => Math.abs(v) > 1e-12);
  const n = d.length;
  if (!n) return { p: 1, nEff: 0, wins: 0 };
  const ranks = d.map((v, i) => ({ a: Math.abs(v), s: Math.sign(v), i }))
    .sort((x, y) => x.a - y.a);
  // average ranks for ties in |d|
  const rk = new Array(n);
  for (let i = 0; i < n;) {
    let j = i; while (j + 1 < n && ranks[j + 1].a === ranks[i].a) j++;
    const avg = (i + j) / 2 + 1;
    for (let k = i; k <= j; k++) rk[k] = avg;
    i = j + 1;
  }
  let Wobs = 0;
  for (let i = 0; i < n; i++) if (ranks[i].s > 0) Wobs += rk[i];
  let ge = 0;
  for (let m = 0; m < (1 << n); m++) {
    let w = 0;
    for (let i = 0; i < n; i++) if (m & (1 << i)) w += rk[i];
    if (w >= Wobs) ge++;
  }
  return { p: ge / (1 << n), nEff: n, wins: d.filter(v => v > 0).length };
}

const out = [["set", "index", "n_datasets", "mean_auc", "mean_mcc_auto", "mean_flag_share",
  "n_abstain", "wins_auc", "p_auc", "wins_mcc", "p_mcc"]];

for (const [name, S] of [["nine", nine], ["seven", seven]]) {
  console.log(`\n=== ${name === "seven" ? "INFERENCE SET (behavioral criterion)" : "all convergent"} `
    + `-- ${S.length} datasets: ${S.join(", ")} ===`);
  console.log("index         AUC    MCC   flagged  abstains   vs ens: AUC        MCC");
  for (const idx of INDICES) {
    const auc = S.map(d => val(d, idx, "auc")), mcc = S.map(d => val(d, idx, "mcc"));
    const fl = S.map(d => val(d, idx, "flag_share"));
    const ab = S.filter(d => val(d, idx, "n_flagged") === 0).length;
    let pa = { p: NaN, wins: NaN }, pm = { p: NaN, wins: NaN };
    if (idx !== "ens") {
      pa = exactSignedRank(S.map(d => val(d, "ens", "auc") - val(d, idx, "auc")));
      pm = exactSignedRank(S.map(d => val(d, "ens", "mcc") - val(d, idx, "mcc")));
    }
    console.log("%s %s  %s  %s%%   %s/%d      %s  %s   %s  %s",
      idx.padEnd(12), mean(auc).toFixed(3), mean(mcc).toFixed(3),
      (100 * mean(fl)).toFixed(1).padStart(5), ab, S.length,
      idx === "ens" ? "  -" : `${pa.wins}/${pa.nEff}`.padStart(4),
      idx === "ens" ? "     " : pa.p.toFixed(3),
      idx === "ens" ? "  -" : `${pm.wins}/${pm.nEff}`.padStart(4),
      idx === "ens" ? "     " : pm.p.toFixed(3));
    out.push([name, idx, S.length, mean(auc).toFixed(4), mean(mcc).toFixed(4),
      mean(fl).toFixed(4), ab, pa.wins, pa.p.toFixed(4), pm.wins, pm.p.toFixed(4)]);
  }
}
fs.writeFileSync("ext_inference_set.csv", out.map(r => r.join(",")).join("\n"));
console.log("\nwrote ext_inference_set.csv");
