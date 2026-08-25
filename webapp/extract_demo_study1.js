/* extract_demo_study1.js — build site/demo_data.js from the Study 1 validation
 * sample (n = 157, 109 items).
 *
 * WHAT IS PUBLISHED, AND WHAT IS NOT. The demo shows the tool's OUTPUT, not its
 * input: this script runs the shipped engine here, offline, and emits only the
 * per-respondent results (careless probability, the ensemble log-odds, the four
 * component indices, flags, the structure diagnostic). The item responses never
 * leave this machine. That is enough for a complete, honest demo — the results
 * view, the ranking plot, the flag table, the CSV downloads and the live re-flag
 * controls all read from the result object alone — while 157 real respondents'
 * answers stay unpublished until the study data is released.
 *
 * The 12 numbers kept per respondent cannot reconstruct their 109 answers.
 *
 * The ground truth travels in the row id (P = careful, C50/C75/C100 = instructed
 * careless on that share of the battery), so a visitor can still check every
 * flag against what was actually collected.
 *
 *   node extract_demo_study1.js
 */
const fs = require("fs");
const path = require("path");
const vm = require("vm");
const HERE = __dirname;

function readCSV(file) {
  const txt = fs.readFileSync(path.join(HERE, file), "utf8").replace(/^﻿/, "");
  const lines = txt.split(/\r?\n/).filter(l => l.length);
  const split = l => l.split(",").map(c => c.replace(/^"|"$/g, ""));
  return { header: split(lines[0]), rows: lines.slice(1).map(split) };
}

const M = readCSV("_study_matrix.csv");     // 157 x 109, item columns only
const L = readCSV("_study_labels.csv");     // careless, group
const I = readCSV("_study_items.csv");      // id, domain, scale, reverse, rmax, rmin

if (M.rows.length !== L.rows.length)
  throw new Error(`matrix/labels length mismatch: ${M.rows.length} vs ${L.rows.length}`);
const items = I.rows.map(r => Object.fromEntries(I.header.map((h, k) => [h, r[k]])));
if (items.length !== M.header.length)
  throw new Error(`items/columns mismatch: ${items.length} vs ${M.header.length}`);

const gi = L.header.indexOf("group");
const idFor = (group, i) => {
  const n = String(i + 1).padStart(3, "0");
  return group === "careful" ? "P" + n : "C" + group.replace("careless_", "") + "_" + n;
};

const ids = [];
const counts = {};
let nMissingCells = 0;
for (let i = 0; i < M.rows.length; i++) {
  const g = L.rows[i][gi];
  counts[g] = (counts[g] || 0) + 1;
  ids.push(idFor(g, i));
  for (const v of M.rows[i]) if (v === "NA") nMissingCells++;
}

/* ---- run the shipped engine on the real responses, here and now ---- */
const ctx = { self: {}, console };
vm.createContext(ctx);
vm.runInContext(fs.readFileSync(path.join(HERE, "site", "reread.js"), "utf8"), ctx);
const R = ctx.self.ReReRe;

/* the app's default Advanced settings, so the published numbers are exactly
 * what a visitor would get if they ran this file through the tool themselves */
const OPTS = { corProp: 0.03, iterations: 200, minPairs: 15, seed: 1, sensitivity: "standard" };

const rows = M.rows.map(r => r.slice());                 // strings; toNum handles "NA"
const B = R.buildMatrix(rows, rows[0].map((_, j) => j));
const res = R.ensemble(B.M, B.n, B.J, OPTS);

/* ---- keep only what the results view reads; drop everything else ---- */
const r4 = a => Array.from(a, v => Number(Number(v).toFixed(4)));
const published = {
  p: r4(res.p),
  eta: r4(res.eta),                    // reflag/sensitivity controls need this
  etaSe: r4(res.etaSe),
  borderline: Array.from(res.borderline, Boolean),
  flagged: Array.from(res.flagged, Boolean),
  nFlagged: res.nFlagged,
  rr: r4(res.rr),
  longstring: Array.from(res.longstring, Number),
  person_total: r4(res.person_total),
  irv: r4(res.irv),
  d2: r4(res.d2),
  nMissing: Array.from(res.nMissing, Number),
  diagnostic: res.diagnostic,
  threshold: res.threshold,
  profileInfo: res.profileInfo,
  partnerGate: res.partnerGate,
  estimatedRate: res.estimatedRate,
  flagMode: res.flagMode,
  twoComponent: res.twoComponent,
  flagStatus: res.flagStatus,
  params: res.params
};
for (const k of ["p", "eta", "flagged", "rr", "person_total"])
  if (!published[k] || published[k].length !== ids.length)
    throw new Error(`published.${k} is missing or the wrong length`);

/* scale inventory, for the caption: this battery is deliberately mixed-scale
 * (1–4, 0–4, 1–5 and 1–7 items side by side). "Random" is the study's decoy
 * block: items with no common factor, kept because respondents answered them. */
const byScale = {};
for (const it of items) {
  const name = it.scale === "Random" ? "unrelated filler" : it.scale;
  const k = name + " (" + it.rmin + "–" + it.rmax + ")";
  byScale[k] = (byScale[k] || 0) + 1;
}
const scaleList = Object.entries(byScale).map(([k, v]) => v + " " + k).join(", ");
const nCareless = M.rows.length - (counts["careful"] || 0);

const meta = {
  n: M.rows.length, items: M.header.length,
  careful: counts["careful"] || 0, careless: nCareless,
  byGroup: counts, scales: byScale, missing: nMissingCells, opts: OPTS
};

const banner =
`/* demo_data.js — PRECOMPUTED results for the CaReReRe validation sample
   (Study 1): ${meta.n} real respondents × ${meta.items} items, ${scaleList}.
   ${meta.careful} answered the battery honestly; ${nCareless} were instructed to answer
   carelessly on a known share of it (${counts["careless_50"]} on 50%, ${counts["careless_75"]} on 75%, ${counts["careless_100"]} on 100%),
   so the ground truth is collected rather than injected or simulated, and it
   travels in the row id (P = careful, C50/C75/C100 = instructed careless).

   THE RESPONSES THEMSELVES ARE NOT PUBLISHED. This file carries only the engine's
   OUTPUT — careless probability, ensemble log-odds, the four component indices,
   flags and the structure diagnostic — computed offline by
   webapp/extract_demo_study1.js at the tool's default settings
   (${JSON.stringify(OPTS)}). Twelve numbers per
   respondent cannot reconstruct 109 answers. The demo therefore shows exactly
   what the tool produces without disclosing the data it was produced from; the
   study data is released separately. Regenerate with extract_demo_study1.js. */
`;

/* Guard, checked BEFORE anything is written: the payload must not carry a single
 * item column id, and must not carry any per-respondent array longer than the
 * roster — either would mean response-level data slipped through. Instrument
 * names in the caption ("60 BFI-2 (1–5)") are deliberate and are not checked. */
const payload = JSON.stringify({ type: "done", res: published, ids, excluded: [] });
for (const it of items)
  if (payload.includes(it.id))
    throw new Error(`REFUSING TO WRITE: item column id "${it.id}" is in the published payload`);
(function assertShapes(v, where) {
  if (Array.isArray(v)) {
    if (v.length > ids.length)
      throw new Error(`REFUSING TO WRITE: array at ${where} has ${v.length} entries (> ${ids.length} respondents)`);
    v.forEach((x, i) => assertShapes(x, `${where}[${i}]`));
  } else if (v && typeof v === "object") {
    for (const k of Object.keys(v)) assertShapes(v[k], `${where}.${k}`);
  }
})(JSON.parse(payload), "payload");

const out =
  banner +
  "self.RERERE_DEMO_RESULT = " + payload + ";\n" +
  "self.RERERE_DEMO_META = " + JSON.stringify(meta) + ";\n";

const dest = path.join(HERE, "site", "demo_data.js");
fs.writeFileSync(dest, out, "utf8");

console.log("wrote " + dest + "  (" + out.length + " bytes)");
console.log("  n=" + meta.n + "  items=" + meta.items + "  careful=" + meta.careful + "  careless=" + meta.careless);
console.log("  groups : " + JSON.stringify(counts));
console.log("  scales : " + scaleList);
console.log("  engine : diagnostic=" + published.diagnostic.level + "  flagged=" + published.nFlagged +
            "  estRate=" + (100 * published.estimatedRate).toFixed(1) + "%");
const tp = ids.filter((id, i) => published.flagged[i] && id[0] === "C").length;
const fp = ids.filter((id, i) => published.flagged[i] && id[0] === "P").length;
console.log("  vs GT  : " + tp + " true positives, " + fp + " false positives out of " + meta.careful + " careful");
console.log("  payload: results only — " + Object.keys(published).length + " fields, no item column ids, " +
            "no array longer than " + ids.length + " (checked before writing)");
