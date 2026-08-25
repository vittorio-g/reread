/* sweep_injection_curve.js — injection curve for every (sub)dataset.
 * At each target prevalence p we build a sample of REAL respondents: k who failed the attention
 * check and m who passed, with k/(k+m) = p. When the careless pool runs out we keep all of them
 * and REDUCE the attentive side instead, so the proportion stays exact and only n shrinks.
 * >=REPS independent resamples per point; we report the mean estimated prevalence (pi-hat).
 * Output: injection_curve.csv  (dataset, target, true_prev, pi_mean, pi_sd, flag_mean, n)
 */
const R = require("./site/rerere.js"); const W = R.WEIGHTS; const fs = require("fs");
let SEED = 777001; const rng = () => { SEED = (SEED*1103515245+12345)&0x7fffffff; return SEED/0x7fffffff; };
const REPS = 10, ITER = 100, NCAP = 1400, NMIN = 180;
const TARGETS = [0.01,0.03,0.05,0.08,0.12,0.16,0.20,0.25,0.30,0.35,0.40];
const D = "../Dataset/gt_benchmark_candidates";
const SETS = [
  ["TISP Germany",  D+"/tisp",               "_labels_DEU.csv", "_matrix_DEU.csv"],
  ["TISP Poland",   D+"/tisp",               "_labels_POL.csv", "_matrix_POL.csv"],
  ["TISP Australia",D+"/tisp",               "_labels_AUS.csv", "_matrix_AUS.csv"],
  ["smarvus England",D+"/smarvus_England",   "_labels2.csv",    "_matrix.csv"],
  ["smarvus Egypt", D+"/smarvus_Egypt",      "_labels2.csv",    "_matrix.csv"],
  ["smarvus Canada",D+"/smarvus_Canada",     "_labels2.csv",    "_matrix.csv"],
  ["warning control",D+"/warning_control",   "_labels2.csv",    "_matrix.csv"],
  ["Kay S1",        D+"/kay_idris_idria/s1", "_labels2.csv",    "_matrix.csv"],
  /* TUMI enters with GT>=1 (33 careless): with GT>=2 it has only 10, which cannot
   * populate any target rate at a usable n. High rates are skipped automatically. */
  ["TUMI",          D+"/tumi",               "_labels1.csv",    "_matrix.csv"],
];
function shuffle(a){for(let i=a.length-1;i>0;i--){const j=Math.floor(rng()*(i+1));const t=a[i];a[i]=a[j];a[j]=t;}return a;}
const out = [["dataset","target","true_prev","pi_mean","pi_sd","flag_mean","n","piraw_mean","estimable_frac"]];
for (const [nm, dir, lf, mf] of SETS) {
  const P = R.parseCSV(fs.readFileSync(dir + "/" + mf, "utf8")), J = P.header.length;
  const rows = P.rows.map(r => r.map(v => { const x = Number(v); return Number.isFinite(x) ? x : NaN; }));
  const y = R.parseCSV(fs.readFileSync(dir + "/" + lf, "utf8")).rows.map(r => Number(r[0]));
  const n0 = Math.min(rows.length, y.length);
  const att = [], car = [];
  for (let i = 0; i < n0; i++) (y[i] ? car : att).push(rows[i]);
  console.log("\n" + nm + "  J=" + J + "  attentive=" + att.length + "  careless=" + car.length);
  for (const p of TARGETS) {
    /* biggest sample with EXACT proportion p: limited either by the careless or the attentive pool */
    let k = Math.floor(Math.min(car.length, att.length * p / (1 - p)));
    let m = Math.round(k * (1 - p) / p);
    if (m > att.length) { m = att.length; k = Math.round(m * p / (1 - p)); }
    if (k + m > NCAP) { const s = NCAP / (k + m); k = Math.max(1, Math.round(k * s)); m = Math.round(m * s); }
    if (k < 3 || k + m < NMIN) { console.log("  p=" + (100*p).toFixed(0) + "% skipped (pool too small)"); continue; }
    const pis = [], fls = [], raws = [], ests = [];
    for (let r = 0; r < REPS; r++) {
      const A = shuffle(att.slice()).slice(0, m), C = shuffle(car.slice()).slice(0, k);
      const sub = A.concat(C), n = sub.length;
      const mat = new Float64Array(n * J);
      for (let i = 0; i < n; i++) for (let j = 0; j < J; j++) mat[i * J + j] = sub[i][j];
      const ef = R.ensembleFeatures(mat, n, J, { seed: 1 + r, iterations: ITER });
      const O = ef.oriented, g = Math.max(0, Math.min(1, (ef.profileInfo - 0.01) / 0.05));
      const eta = new Array(n);
      for (let i = 0; i < n; i++)
        eta[i] = W.b0 + W.rr*O.rr[i] + g*(W.longstring*O.longstring[i] + W.person_total*O.person_total[i]);
      const af = R.autoFlag(eta, "standard");
      pis.push(af.pi); fls.push(af.estRate);
      raws.push(af.piRaw !== undefined ? af.piRaw : af.pi);
      ests.push(af.piEstimable ? 1 : 0);
    }
    const mean = a => a.reduce((s,v)=>s+v,0)/a.length;
    const sd = a => { const mu = mean(a); return Math.sqrt(mean(a.map(v => (v-mu)*(v-mu)))); };
    const truep = k / (k + m);
    out.push([nm, p.toFixed(2), truep.toFixed(4), mean(pis).toFixed(4), sd(pis).toFixed(4),
              mean(fls).toFixed(4), k + m, mean(raws).toFixed(4), mean(ests).toFixed(2)]);
    console.log("  true=" + (100*truep).toFixed(1).padStart(5) + "%  n=" + String(k+m).padStart(5) +
                "  piHat=" + (100*mean(pis)).toFixed(1).padStart(5) + "% (sd " + (100*sd(pis)).toFixed(1) +
                ")  flagged=" + (100*mean(fls)).toFixed(1) + "%");
  }
}
fs.writeFileSync("injection_curve.csv", out.map(r => r.join(",")).join("\n"));
console.log("\nwrote injection_curve.csv");
