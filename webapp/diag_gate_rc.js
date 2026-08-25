/* diag_gate_rc.js — why does the shipped gate abstain so often when it is handed
 * rc alone, while it opens on the ensemble for every external dataset?
 *
 * The gate has two branches: a BIC mixture test (does a second component exist,
 * with pi in range and separation > 1?) and an excess-tail binomial test (is the
 * mass above m + 2.5*sigma heavier than the 0.62% a careless-free sample would
 * put there?). Abstention means BOTH said no. This prints, per dataset and per
 * score, which clause failed and -- the quantity that actually settles it -- where
 * the true careless respondents sit in the score's own distribution, expressed in
 * sigma units above the attentive mode. If the labelled careless are not in the
 * far tail of a score, no far-tail rule can find them, and abstaining is correct
 * rather than a failure of the gate.
 *
 *   node diag_gate_rc.js   -> diag_gate_rc.csv
 */
const fs = require("fs"), path = require("path");
const R = require("./site/rerere.js");
const DIR = "ext_bench2";
const ORDER = ["kay_s2","warning300","kay_s1","opsy_16pf","smarvus","kay_s6","kay_s5","douglas","krause","duckworth"];
const LABEL = { kay_s2:"Kay S2", warning300:"warning IPIP", kay_s1:"Kay S1", opsy_16pf:"opsy 16PF",
  smarvus:"smarvus", kay_s6:"Kay S6", kay_s5:"Kay S5", douglas:"Douglas", krause:"Krause", duckworth:"Duckworth" };

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
  for (let i = 0; i < y.length; i++) (y[i] ? po : ne).push(s[i]);
  let c = 0; for (const a of po) for (const b of ne) c += a > b ? 1 : a === b ? 0.5 : 0;
  return c / (po.length * ne.length);
};
const conf = (flag, y) => {
  let tp = 0, fp = 0, fn = 0, tn = 0;
  for (let i = 0; i < y.length; i++) { if (y[i] === 1) { flag[i] ? tp++ : fn++; } else { flag[i] ? fp++ : tn++; } }
  const d = Math.sqrt((tp + fp) * (tp + fn) * (tn + fp) * (tn + fn));
  return d ? (tp * tn - fp * fn) / d : 0;
};
/* the threshold that maximises MCC, expressed in sigma units above the attentive
 * mode: how far into the tail the careless actually are for this score */
function bestCutInSigma(s, y, m, sig) {
  const cand = Array.from(new Set(s)).sort((a, b) => a - b);
  const step = Math.max(1, Math.floor(cand.length / 400));
  let best = -1, at = NaN;
  for (let i = 0; i < cand.length; i += step) {
    const v = conf(s.map(x => x >= cand[i]), y);
    if (v > best) { best = v; at = cand[i]; }
  }
  return { mcc: best, z: (at - m) / sig };
}

const out = [["dataset","n","score","auc","dBIC","pi_mix","sep","bic_gate","k_tail","exp_tail","tail_p","tail_gate","gate","n_flagged","mcc_auto","best_cut_z","best_mcc","careless_median_z"]];
console.log("dataset       score |  AUC  | dBIC     pi_mix  sep   BIC? |  kTail(exp)   tailP     tail? | GATE | flagged  MCC  | best cut at  careless median");
console.log("-".repeat(140));
for (const nm of ORDER) {
  const Y = readCsv(nm + "_y.csv"), y = Y.y, n = y.length;
  const S = { ens: Y.eta, rc: Y.rr, PT: Y.pt };
  for (const k of ["ens", "rc", "PT"]) {
    const s = S[k];
    const fl = R.autoFlag(s, "standard");
    const g = R.gaussMix1D(s);
    const pooled = Math.sqrt((g.v1 + g.v2) / 2) || 1;
    const sep = (g.m2 - g.m1) / pooled;
    const lf = R.leftFit(s);
    const bc = bestCutInSigma(s, y, lf.m, lf.sigma);
    /* where the labelled careless actually sit, in sigma units above the mode */
    const cz = s.filter((_, i) => y[i] === 1).map(v => (v - lf.m) / lf.sigma).sort((a, b) => a - b);
    const medZ = cz[Math.floor(cz.length / 2)];
    const expTail = 0.0062 * n;
    const nf = fl.flagged.filter(Boolean).length;
    console.log(`${LABEL[nm].padEnd(13)} ${k.padEnd(4)}  | ${auc(s, y).toFixed(3)} | ` +
      `${(g.bic1 - g.bic2).toFixed(0).padStart(7)}  ${g.pi.toFixed(2)}   ${sep.toFixed(2)}  ${fl.bicGate ? "YES" : "no "} | ` +
      `${String(fl.tailK).padStart(5)}(${expTail.toFixed(1).padStart(5)})  ${fl.tailP.toExponential(1).padStart(8)}  ${fl.tailP < 1e-6 ? "YES" : "no "} | ` +
      `${fl.twoComp ? "OPEN " : "SHUT "}| ${String(nf).padStart(5)}  ${conf(fl.flagged, y).toFixed(2).padStart(5)} | ` +
      `${bc.z.toFixed(2).padStart(5)} sigma   ${medZ.toFixed(2).padStart(6)} sigma`);
    out.push([nm, n, k, auc(s, y).toFixed(4), (g.bic1 - g.bic2).toFixed(1), g.pi.toFixed(4), sep.toFixed(3),
      fl.bicGate, fl.tailK, expTail.toFixed(2), fl.tailP.toExponential(3), fl.tailP < 1e-6, fl.twoComp,
      nf, conf(fl.flagged, y).toFixed(4), bc.z.toFixed(3), bc.mcc.toFixed(4), medZ.toFixed(3)]);
  }
  console.log("");
}
fs.writeFileSync("diag_gate_rc.csv", out.map(r => r.join(",")).join("\n"));
console.log("wrote diag_gate_rc.csv");
