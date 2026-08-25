/* check_demo_study1.js — audit what the demo publishes.
 *
 * Two jobs. (1) Confirm the payload carries results only: no item column ids, no
 * per-respondent array longer than the roster. (2) Score the published results
 * against the ground truth carried in the row ids, so the numbers a visitor sees
 * are checked rather than assumed.
 *
 * Reads site/demo_data.js alone — it needs no access to the responses, which is
 * the point.
 *
 *   node check_demo_study1.js
 */
const fs = require("fs"), path = require("path"), vm = require("vm");
const HERE = __dirname;

const ctx = { self: {}, console };
vm.createContext(ctx);
vm.runInContext(fs.readFileSync(path.join(HERE, "site", "demo_data.js"), "utf8"), ctx);
const D = ctx.self.RERERE_DEMO_RESULT, META = ctx.self.RERERE_DEMO_META;
if (!D || !D.res) throw new Error("demo_data.js does not define RERERE_DEMO_RESULT");
const res = D.res, ids = D.ids;

/* ---- 1. disclosure audit ---- */
const raw = JSON.stringify(D);
const itemIds = fs.readFileSync(path.join(HERE, "_study_items.csv"), "utf8")
  .split(/\r?\n/).slice(1).filter(Boolean).map(l => l.split(",")[0]);
const leaked = itemIds.filter(id => raw.includes(id));
const tooLong = [];
(function walk(v, where) {
  if (Array.isArray(v)) { if (v.length > ids.length) tooLong.push(`${where} (${v.length})`); v.forEach((x, i) => walk(x, `${where}[${i}]`)); }
  else if (v && typeof v === "object") for (const k of Object.keys(v)) walk(v[k], `${where}.${k}`);
})(D, "payload");
console.log(`disclosure : item column ids present = ${leaked.length}` +
            `   arrays longer than the ${ids.length} respondents = ${tooLong.length}` +
            `   payload = ${(raw.length / 1024).toFixed(1)} KB`);
if (leaked.length || tooLong.length) {
  console.log("  LEAK: " + [...leaked, ...tooLong].join(", "));
  process.exitCode = 1;
}

/* ---- 2. results vs the collected ground truth ---- */
const truth = ids.map(id => id[0] === "C");
const auc = (s, y) => {
  const idx = s.map((v, i) => [v, i]).sort((a, b) => a[0] - b[0]);
  const rank = new Array(s.length); let i = 0;
  while (i < idx.length) { let j = i; while (j + 1 < idx.length && idx[j + 1][0] === idx[i][0]) j++;
    const r = (i + j) / 2 + 1; for (let k = i; k <= j; k++) rank[idx[k][1]] = r; i = j + 1; }
  const np = y.filter(Boolean).length, nn = y.length - np;
  let sp = 0; for (let k = 0; k < y.length; k++) if (y[k]) sp += rank[k];
  return (sp - np * (np + 1) / 2) / (np * nn);
};
let tp = 0, fp = 0, tn = 0, fn = 0;
for (let i = 0; i < truth.length; i++) {
  if (res.flagged[i] && truth[i]) tp++; else if (res.flagged[i]) fp++;
  else if (truth[i]) fn++; else tn++;
}
const den = Math.sqrt((tp + fp) * (tp + fn) * (tn + fp) * (tn + fn));
const f3 = x => Number(x).toFixed(3);

console.log(`dataset    : n=${META.n}  items=${META.items}  careless=${META.careless} ` +
            `(${(100 * META.careless / META.n).toFixed(0)}%)  settings ${JSON.stringify(META.opts)}`);
console.log(`ranking    : ensemble AUC ${f3(auc(res.eta, truth))}   rr ${f3(auc(res.rr.map(v => -v), truth))}   ` +
            `LongString ${f3(auc(res.longstring, truth))}   Person-Total ${f3(auc(res.person_total.map(v => -v), truth))}`);
console.log(`diagnostic : ${res.diagnostic.level}  (signal ${res.diagnostic.signal_strength.toFixed(2)}x noise, ` +
            `top-pair mean |r| ${res.diagnostic.top_mean_r.toFixed(2)})   partner gate ${f3(res.partnerGate)}`);
console.log(`flagging   : ${res.nFlagged} flagged (${(100 * res.estimatedRate).toFixed(0)}%)   ` +
            `MCC ${f3(den > 0 ? (tp * tn - fp * fn) / den : 0)}  precision ${f3(tp / (tp + fp))}  recall ${f3(tp / (tp + fn))}` +
            `   (TP ${tp}, FP ${fp}, TN ${tn}, FN ${fn})`);

const byGroup = {};
ids.forEach((id, i) => {
  const g = id[0] === "P" ? "careful" : "careless " + id.slice(1).split("_")[0] + "%";
  byGroup[g] = byGroup[g] || { n: 0, flagged: 0 };
  byGroup[g].n++; if (res.flagged[i]) byGroup[g].flagged++;
});
for (const g of Object.keys(byGroup).sort())
  console.log(`             ${g.padEnd(14)} ${String(byGroup[g].flagged).padStart(3)}/${String(byGroup[g].n).padStart(3)}` +
              ` flagged (${(100 * byGroup[g].flagged / byGroup[g].n).toFixed(0)}%)`);
