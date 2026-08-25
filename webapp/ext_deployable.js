/* ext_deployable.js — Study 2, scored the way it would actually be used.
 *
 * The external benchmark reported AUC (threshold-free) and oracle MCC (the best
 * cutoff for each score, chosen with the labels). Neither is what an analyst
 * gets: a real screening decision needs a threshold, and nobody has the labels
 * that would place it. So every index here is put through the SAME label-free
 * calibration as the ensemble -- the structure gate, then m + 2.5*sigma on the
 * one-sided attentive fit -- and scored on what that decision actually yields.
 *
 * Lending the calibration to the competitors is the only way to ask the question
 * at all: none of them ships a cutoff rule. It is a courtesy rather than a
 * neutral test, since the calibration was developed on the ensemble score, and
 * the results are reported with that stated.
 *
 * An index whose gate abstains scores MCC = 0 by construction. That is not a
 * defect of the metric: in deployment it means the index flags nobody, which is
 * exactly the outcome the analyst would face, and it is reported as such.
 *
 *   node ext_deployable.js   -> ext_deployable.csv, ext_deployable_pairs.csv
 */
const fs = require("fs"), path = require("path");
const R = require("./site/rerere.js");
const DIR = "ext_bench2";

const ORDER = ["kay_s2","warning300","kay_s1","opsy_16pf","smarvus","kay_s6","kay_s5","douglas","krause","duckworth"];
const LABEL = { kay_s2:"Kay S2", warning300:"warning IPIP-NEO", kay_s1:"Kay S1", opsy_16pf:"opsy 16PF",
  smarvus:"smarvus", kay_s6:"Kay S6", kay_s5:"Kay S5", douglas:"Douglas 2023", krause:"Krause youth",
  duckworth:"Duckworth VCL" };
const CONV = ORDER.filter(d => d !== "duckworth");     // the negative control never enters a mean

function readCsv(file) {
  const lines = fs.readFileSync(path.join(DIR, file), "utf8").trim().split(/\r?\n/);
  const head = lines[0].split(",").map(h => h.replace(/"/g, ""));
  const cols = {}; head.forEach(h => cols[h] = []);
  for (let i = 1; i < lines.length; i++) {
    const f = lines[i].split(",");
    head.forEach((h, k) => cols[h].push(f[k] === "NA" ? NaN : Number(f[k])));
  }
  return cols;
}
const auc = (s, y) => {
  const po = [], ne = [];
  for (let i = 0; i < y.length; i++) if (Number.isFinite(s[i])) (y[i] ? po : ne).push(s[i]);
  if (!po.length || !ne.length) return NaN;
  let c = 0; for (const a of po) for (const b of ne) c += a > b ? 1 : a === b ? 0.5 : 0;
  return c / (po.length * ne.length);
};
function conf(flag, y) {
  let tp = 0, fp = 0, fn = 0, tn = 0;
  for (let i = 0; i < y.length; i++) { if (y[i] === 1) { flag[i] ? tp++ : fn++; } else { flag[i] ? fp++ : tn++; } }
  const d = Math.sqrt((tp + fp) * (tp + fn) * (tn + fp) * (tn + fn));
  return { mcc: d ? (tp * tn - fp * fn) / d : 0, prec: tp + fp ? tp / (tp + fp) : NaN, rec: tp / (tp + fn), nf: tp + fp };
}
/* exact Wilcoxon signed-rank, one-sided (A > B); n <= 20 pairs */
function wilcoxon(a, b) {
  const d = a.map((v, i) => v - b[i]).filter(v => Number.isFinite(v) && v !== 0);
  const n = d.length; if (!n) return { p: NaN, n: 0, wins: 0 };
  const abs = d.map(Math.abs), idx = abs.map((v, i) => [v, i]).sort((p, q) => p[0] - q[0]);
  const rank = new Array(n); let i = 0;
  while (i < n) { let j = i; while (j + 1 < n && idx[j + 1][0] === idx[i][0]) j++;
    const r = (i + j) / 2 + 1; for (let k = i; k <= j; k++) rank[idx[k][1]] = r; i = j + 1; }
  let W = 0; for (let k = 0; k < n; k++) if (d[k] > 0) W += rank[k];
  let ge = 0; const tot = 1 << n;
  for (let m = 0; m < tot; m++) { let w = 0; for (let k = 0; k < n; k++) if (m & (1 << k)) w += rank[k]; if (w >= W) ge++; }
  return { p: ge / tot, n, wins: d.filter(v => v > 0).length };
}

const NAMES = ["ens", "rc", "PT", "Mahalanobis", "IRV", "LongString", "PsychSyn"];
const rowsOut = [["dataset","n","rate","index","auc","gate_open","n_flagged","flag_share","mcc","precision","recall","pi_hat"]];
const store = {};    // store[index][dataset] = {...}

for (const nm of ORDER) {
  const Y = readCsv(nm + "_y.csv");
  const C = fs.existsSync(path.join(DIR, nm + "_comp.csv")) ? readCsv(nm + "_comp.csv") : null;
  const y = Y.y, n = y.length, rate = y.reduce((s, v) => s + v, 0) / n;
  const S = { ens: Y.eta, rc: Y.rr, PT: Y.pt };
  if (C) for (const k of ["Mahalanobis", "IRV", "LongString", "PsychSyn"]) S[k] = C[k];
  for (const k of NAMES) {
    const s = S[k];
    if (!s || s.every(v => !Number.isFinite(v))) { store[k] = store[k] || {}; store[k][nm] = null; continue; }
    /* a few competitor scores can be non-finite for isolated respondents; the
     * calibration needs a complete vector, so those are held out of the fit and
     * left unflagged, which is what a tool would do with a missing score */
    const ok = s.map(v => Number.isFinite(v));
    const sub = s.filter((v, i) => ok[i]);
    const fl = R.autoFlag(sub, "standard");
    const flagged = new Array(n).fill(false);
    let j = 0; for (let i = 0; i < n; i++) if (ok[i]) flagged[i] = fl.flagged[j++];
    const c = conf(flagged, y), a = auc(s, y);
    store[k] = store[k] || {};
    store[k][nm] = { auc: a, mcc: c.mcc, prec: c.prec, rec: c.rec, nf: c.nf, share: c.nf / n,
                     gate: fl.twoComp, pi: fl.pi };
    rowsOut.push([nm, n, rate.toFixed(3), k, Number.isFinite(a) ? a.toFixed(3) : "NA", fl.twoComp,
      c.nf, (c.nf / n).toFixed(4), c.mcc.toFixed(3),
      Number.isFinite(c.prec) ? c.prec.toFixed(3) : "NA", c.rec.toFixed(3), fl.pi.toFixed(3)]);
  }
}

/* ---------- per-dataset table ---------- */
const f3 = x => (x == null || !Number.isFinite(x)) ? "  --  " : x.toFixed(3).replace(/^0/, " ").replace(/^-0/, "-").padStart(6);
console.log("=== Study 2, deployable: every index cut by the shipped label-free calibration ===\n");
console.log("AUC (threshold-free)");
console.log("dataset            " + NAMES.map(k => k.slice(0, 6).padStart(7)).join(""));
for (const nm of ORDER) console.log(LABEL[nm].padEnd(18) + NAMES.map(k => f3(store[k][nm] && store[k][nm].auc) + " ").join(""));
console.log("\nMCC at the automatic cut");
console.log("dataset            " + NAMES.map(k => k.slice(0, 6).padStart(7)).join(""));
for (const nm of ORDER)
  console.log(LABEL[nm].padEnd(18) + NAMES.map(k => {
    const v = store[k][nm]; if (!v) return "  --   ";
    return f3(v.mcc) + (v.nf === 0 ? "*" : " ");
  }).join(""));
console.log("  (* the calibration flags nobody)");
console.log("\nshare of the sample flagged");
console.log("dataset            " + NAMES.map(k => k.slice(0, 6).padStart(7)).join(""));
for (const nm of ORDER)
  console.log(LABEL[nm].padEnd(18) + NAMES.map(k => {
    const v = store[k][nm]; return v ? ((100 * v.share).toFixed(0) + "%").padStart(6) + " " : "  --   ";
  }).join(""));

/* ---------- means and paired tests over the nine convergent datasets ---------- */
const mean = arr => { const v = arr.filter(Number.isFinite); return v.length ? v.reduce((s, x) => s + x, 0) / v.length : NaN; };
console.log("\n=== mean over the nine convergent datasets (Duckworth VCL excluded: negative control) ===");
console.log("index          mean AUC   mean MCC(auto)   mean flagged   datasets flagging nobody");
const summary = [["index","mean_auc","mean_mcc_auto","mean_flag_share","n_abstain"]];
for (const k of NAMES) {
  const a = mean(CONV.map(d => store[k][d] && store[k][d].auc));
  const m = mean(CONV.map(d => store[k][d] && store[k][d].mcc));
  const sh = mean(CONV.map(d => store[k][d] && store[k][d].share));
  const ab = CONV.filter(d => store[k][d] && store[k][d].nf === 0).length;
  console.log(`${k.padEnd(14)} ${a.toFixed(3)}      ${m.toFixed(3)}          ${(100 * sh).toFixed(1)}%          ${ab}/9`);
  summary.push([k, a.toFixed(4), m.toFixed(4), sh.toFixed(4), ab]);
}
console.log("\n=== paired tests across the nine convergent datasets (one-sided, exact signed-rank) ===");
console.log("comparison                       AUC: wins  p    |  MCC(auto): wins  p");
const pairs = [["ens","Mahalanobis"],["ens","IRV"],["ens","LongString"],["ens","PsychSyn"],["ens","PT"],["ens","rc"],
               ["rc","Mahalanobis"],["rc","IRV"],["rc","LongString"],["rc","PsychSyn"],["rc","PT"]];
const pout = [["a","b","auc_wins","auc_n","auc_p","mcc_wins","mcc_n","mcc_p"]];
for (const [x, z] of pairs) {
  const ax = CONV.map(d => store[x][d] && store[x][d].auc), az = CONV.map(d => store[z][d] && store[z][d].auc);
  const mx = CONV.map(d => store[x][d] && store[x][d].mcc), mz = CONV.map(d => store[z][d] && store[z][d].mcc);
  const wa = wilcoxon(ax, az), wm = wilcoxon(mx, mz);
  console.log(`${(x + " vs " + z).padEnd(32)} ${wa.wins}/${wa.n}   ${wa.p.toFixed(3)}  |  ${wm.wins}/${wm.n}   ${wm.p.toFixed(3)}`);
  pout.push([x, z, wa.wins, wa.n, wa.p.toFixed(4), wm.wins, wm.n, wm.p.toFixed(4)]);
}
fs.writeFileSync("ext_deployable.csv", rowsOut.map(r => r.join(",")).join("\n"));
fs.writeFileSync("ext_deployable_pairs.csv", pout.map(r => r.join(",")).join("\n"));
fs.writeFileSync("ext_deployable_summary.csv", summary.map(r => r.join(",")).join("\n"));
console.log("\nwrote ext_deployable.csv, ext_deployable_pairs.csv, ext_deployable_summary.csv");
