/* sweep_bogus_strictness.js — how does the QUALITY of the prevalence estimate improve as the
 * carelessness criterion gets stricter (>=k failed checks, k rising)?
 *
 * Rationale: the absolute gap between pi-hat and the check-based rate is NOT interpretable,
 * because the "attentive" pool is itself contaminated (the checks miss people). What IS
 * interpretable is the SLOPE of pi-hat against the injected careless fraction: each injected
 * careless respondent should move the estimate by one respondent (slope 1). A stricter criterion
 * selects more definite careless -> prediction: slope rises toward 1 with k.
 *
 * Design per (dataset, k): careless pool = fails>=k, attentive pool = fails<k (the honest
 * complement - it absorbs the milder failers, so the intercept may rise with k; the slope is the
 * quality metric). Injection: fixed-proportion mixing of REAL respondents, REPS resamples per
 * target rate; every resample is written out so the analysis can bootstrap over them.
 */
const R = require("./site/rerere.js"); const W = R.WEIGHTS; const fs = require("fs");
let SEED = 909090; const rng = () => { SEED = (SEED*1103515245+12345)&0x7fffffff; return SEED/0x7fffffff; };
const REPS = 10, ITER = 100, NCAP = 1200, NMIN = 180;
const TARGETS = [0.02,0.05,0.08,0.12,0.16,0.20,0.25,0.30,0.35,0.40];
const D = "../Dataset/gt_benchmark_candidates";
const SETS = [
  ["warning IPIP", D + "/warning_ipipneo300", [1,2,3]],
  ["smarvus",      D + "/smarvus",            [1,2,3,4,5]],
  ["Kay S1",       D + "/kay_idris_idria/s1", [1,2,3]],
  ["Kay S2",       D + "/kay_idris_idria/s2", [1,2]],
];
function shuffle(a){for(let i=a.length-1;i>0;i--){const j=Math.floor(rng()*(i+1));const t=a[i];a[i]=a[j];a[j]=t;}return a;}
const out = [["dataset","k","target","true_prev","rep","n","piHat","piRaw","estimable","flagged","heavy"]];
for (const [nm, dir, KS] of SETS) {
  const P = R.parseCSV(fs.readFileSync(dir + "/_matrix.csv", "utf8")), J = P.header.length;
  const rows = P.rows.map(r => r.map(v => { const x = Number(v); return Number.isFinite(x) ? x : NaN; }));
  const nf = R.parseCSV(fs.readFileSync(dir + "/_failcount.csv", "utf8")).rows.map(r => Number(r[0]));
  const n0 = Math.min(rows.length, nf.length);
  for (const k of KS) {
    const att = [], car = [];
    for (let i = 0; i < n0; i++) (nf[i] >= k ? car : att).push(rows[i]);
    console.log("\n" + nm + "  k>=" + k + "   (J=" + J + ", att=" + att.length + ", car=" + car.length + ")");
    for (const p of TARGETS) {
      let kk = Math.floor(Math.min(car.length, att.length * p / (1 - p)));
      let m = Math.round(kk * (1 - p) / p);
      if (m > att.length) { m = att.length; kk = Math.round(m * p / (1 - p)); }
      if (kk + m > NCAP) { const s = NCAP / (kk + m); kk = Math.max(1, Math.round(kk * s)); m = Math.round(m * s); }
      if (kk < 3 || kk + m < NMIN) continue;
      const truep = kk / (kk + m);
      let piBar = 0, nEst = 0;
      for (let r = 0; r < REPS; r++) {
        const A = shuffle(att.slice()).slice(0, m), C = shuffle(car.slice()).slice(0, kk);
        const sub = A.concat(C), n = sub.length;
        const mat = new Float64Array(n * J);
        for (let i = 0; i < n; i++) for (let j = 0; j < J; j++) mat[i * J + j] = sub[i][j];
        const ef = R.ensembleFeatures(mat, n, J, { seed: 11 + r, iterations: ITER });
        const O = ef.oriented, g = Math.max(0, Math.min(1, (ef.profileInfo - 0.01) / 0.05));
        const eta = new Array(n);
        for (let i = 0; i < n; i++)
          eta[i] = W.b0 + W.rr*O.rr[i] + g*(W.longstring*O.longstring[i] + W.person_total*O.person_total[i]);
        const af = R.autoFlag(eta, "standard");
        out.push([nm, k, p.toFixed(2), truep.toFixed(4), r, n, af.pi.toFixed(4),
                  (af.piRaw !== undefined ? af.piRaw : af.pi).toFixed(4),
                  af.piEstimable ? 1 : 0, af.estRate.toFixed(4), af.heavyContamination ? 1 : 0]);
        piBar += af.pi; nEst += af.piEstimable ? 1 : 0;
      }
      console.log("  true=" + (100*truep).toFixed(1).padStart(5) + "%  n=" + String(kk+m).padStart(5) +
                  "  piHat=" + (100*piBar/REPS).toFixed(1).padStart(5) + "%  estimable " + nEst + "/" + REPS);
    }
  }
}
fs.writeFileSync("strictness_sweep.csv", out.map(r => r.join(",")).join("\n"));
console.log("\nwrote strictness_sweep.csv (" + (out.length - 1) + " resample rows)");
