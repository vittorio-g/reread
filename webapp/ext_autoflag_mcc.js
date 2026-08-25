/* ext_autoflag_mcc.js — Table 4 currently reports the ORACLE (best-threshold) MCC,
 * which is the most generous possible treatment of every column: each score is
 * given the cutoff that maximises its own MCC, chosen with the labels in hand.
 * No analyst has that. This adds the DEPLOYABLE counterpart: the MCC each score
 * reaches when it is cut by the shipped label-free calibration (structure gate +
 * one-sided fit of the attentive mode, z = 2.5).
 *
 * Applying that same algorithm to rc and to Person-Total is the fair comparison
 * and also the informative one: it lends the competitors the calibration layer
 * they do not have, so the table separates "which score ranks better" from "which
 * score can be turned into a decision without labels".
 *
 * Scores come from the committed ext_bench2 dump (eta, rc, PT per respondent,
 * oriented high = careless), so nothing is re-derived here.
 *
 *   node ext_autoflag_mcc.js   -> ext_autoflag_mcc.csv
 */
const fs = require("fs"), path = require("path");
const R = require("./site/reread.js");
const DIR = "ext_bench2";

const ORDER = ["kay_s2","warning300","kay_s1","opsy_16pf","smarvus","kay_s6","kay_s5","douglas","krause","duckworth"];
const LABEL = { kay_s2:"Kay S2", warning300:"warning IPIP-NEO", kay_s1:"Kay S1", opsy_16pf:"opsy 16PF",
  smarvus:"smarvus", kay_s6:"Kay S6", kay_s5:"Kay S5", douglas:"Douglas 2023", krause:"Krause youth",
  duckworth:"Duckworth VCL" };

function readY(name) {
  const txt = fs.readFileSync(path.join(DIR, name + "_y.csv"), "utf8").trim().split(/\r?\n/);
  const head = txt[0].split(",");
  const col = {}; head.forEach((h, i) => col[h.replace(/"/g, "")] = i);
  const y = [], rr = [], eta = [], pt = [];
  for (let i = 1; i < txt.length; i++) {
    const f = txt[i].split(",");
    y.push(Number(f[col.y])); rr.push(Number(f[col.rr]));
    eta.push(Number(f[col.eta])); pt.push(Number(f[col.pt]));
  }
  return { y, rr, eta, pt };
}
function conf(flag, y) {
  let tp = 0, fp = 0, fn = 0, tn = 0;
  for (let i = 0; i < y.length; i++) { if (y[i] === 1) { flag[i] ? tp++ : fn++; } else { flag[i] ? fp++ : tn++; } }
  const d = Math.sqrt((tp + fp) * (tp + fn) * (tn + fp) * (tn + fn));
  return { mcc: d ? (tp * tn - fp * fn) / d : 0, prec: tp + fp ? tp / (tp + fp) : 0, rec: tp / (tp + fn), n: tp + fp };
}
function oracleMCC(s, y) {
  const srt = Array.from(new Set(s)).sort((a, b) => a - b);
  const step = Math.max(1, Math.floor(srt.length / 400));
  let best = -1;
  for (let i = 0; i < srt.length; i += step) { const c = conf(s.map(v => v >= srt[i]), y); if (c.mcc > best) best = c.mcc; }
  return best;
}

const out = [["dataset","n","rate","score","auc","oracle_mcc","gate_open","n_flagged","flag_share","auto_mcc","auto_prec","auto_rec","pi_hat"]];
const auc = (s, y) => { const po = [], ne = []; for (let i = 0; i < y.length; i++) (y[i] ? po : ne).push(s[i]);
  let c = 0; for (const a of po) for (const b of ne) c += a > b ? 1 : a === b ? 0.5 : 0; return c / (po.length * ne.length); };

console.log("dataset            n      rate |          AUC          |     MCC (oracle)      |   MCC (automatic cut) |  flagged share");
console.log("                                |  ens    rc     PT     |  ens    rc     PT     |  ens    rc     PT     |  ens    rc     PT");
console.log("-".repeat(126));
for (const nm of ORDER) {
  const D = readY(nm), y = D.y, rate = y.reduce((s, v) => s + v, 0) / y.length;
  const S = { ens: D.eta, rc: D.rr, PT: D.pt };
  const A = {}, O = {}, F = {};
  for (const k of Object.keys(S)) {
    A[k] = auc(S[k], y); O[k] = oracleMCC(S[k], y);
    const fl = R.autoFlag(S[k], "standard"), c = conf(fl.flagged, y);
    F[k] = { fl, c, share: fl.flagged.filter(Boolean).length / y.length };
    out.push([nm, y.length, rate.toFixed(3), k, A[k].toFixed(3), O[k].toFixed(3), fl.twoComp,
      fl.flagged.filter(Boolean).length, F[k].share.toFixed(4), c.mcc.toFixed(3), c.prec.toFixed(3), c.rec.toFixed(3),
      (fl.pi != null ? fl.pi.toFixed(3) : "")]);
  }
  const f3 = x => x.toFixed(3).replace(/^0/, " ").padStart(6);
  const pc = x => (100 * x).toFixed(0).padStart(5) + "%";
  console.log(`${LABEL[nm].padEnd(18)} ${String(y.length).padStart(5)} ${rate.toFixed(3)} |` +
    `${f3(A.ens)}${f3(A.rc)}${f3(A.PT)} |` +
    `${f3(O.ens)}${f3(O.rc)}${f3(O.PT)} |` +
    `${f3(F.ens.c.mcc)}${f3(F.rc.c.mcc)}${f3(F.PT.c.mcc)} |` +
    `${pc(F.ens.share)}${pc(F.rc.share)}${pc(F.PT.share)}` +
    (F.ens.fl.twoComp ? "" : "   [ens gate closed]"));
}
fs.writeFileSync("ext_autoflag_mcc.csv", out.map(r => r.join(",")).join("\n"));
console.log("\nwrote ext_autoflag_mcc.csv");

/* means over the nine convergent-validity datasets (Duckworth VCL is the
 * negative control and is excluded from any average) */
const CONV = ORDER.filter(d => d !== "duckworth");
const rows = out.slice(1);
const pick = (d, k, col) => Number(rows.find(r => r[0] === d && r[3] === k)[col]);
console.log("\nmean over the nine convergent datasets");
console.log("               AUC            MCC (oracle)        MCC (automatic cut)");
for (const k of ["ens", "rc", "PT"]) {
  const a = CONV.map(d => pick(d, k, 4)).reduce((s, v) => s + v, 0) / CONV.length;
  const o = CONV.map(d => pick(d, k, 5)).reduce((s, v) => s + v, 0) / CONV.length;
  const m = CONV.map(d => pick(d, k, 9)).reduce((s, v) => s + v, 0) / CONV.length;
  console.log(`  ${k.padEnd(4)}      ${a.toFixed(3)}            ${o.toFixed(3)}                ${m.toFixed(3)}`);
}
const wins = (k1, k2, col) => CONV.filter(d => pick(d, k1, col) > pick(d, k2, col)).length;
console.log(`\nens beats PT on oracle MCC in ${wins("ens","PT",5)}/9 datasets, on the automatic cut in ${wins("ens","PT",9)}/9`);
console.log(`rc  gate stays closed (nothing flagged) on ${CONV.filter(d => pick(d,"rc",8) === 0).length}/9 datasets`);
