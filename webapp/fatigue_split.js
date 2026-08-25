/* fatigue_split.js — does response quality decay across the questionnaire?
 *
 * DESIGN. For each dataset we split the items in ADMINISTRATION (column) order into a first
 * and a second half, score each half independently, and compare the two WITHIN respondent.
 *   - rc (the permutation z) is the right metric here: unlike the ensemble's robust-z features
 *     it has an absolute meaning (coherence relative to that person's own chance baseline),
 *     so it is comparable across halves. eta/robust-z would be centred within half by
 *     construction and could not show a level difference.
 *   - LongString (raw longest run) is reported too: straightlining is the other fatigue signal.
 *
 * CONFOUND + CONTROL. The two halves differ in CONTENT (different scales, different numbers of
 * factors), not only in position, so a raw first-vs-second difference is not by itself evidence
 * of fatigue. We therefore also compute RANDOM splits of the same size: the random-split delta
 * is the null for "splitting this questionnaire in two changes rc by this much". Fatigue means
 * the positional delta is more negative than the random-split delta.
 *
 * CAVEAT (per dataset): column order is only a proxy for administration order. Where the study
 * randomized item order, the positional split is a null by construction -- such datasets act as
 * negative controls rather than tests.
 */
const fs = require("fs");
const path = require("path");
const R = require("./site/reread.js");

const ITERS = 100, NMAX = 2500, RAND_SPLITS = 10, SEED = 20260806;
const rng = R._rng(SEED);

function readMatrix(file, hasHeader) {
  const txt = fs.readFileSync(file, "utf8");
  const P = R.parseCSV(txt);
  let rows = P.rows;
  if (hasHeader === false) rows = [P.header].concat(P.rows);   // header was actually data
  const J = rows[0].length;
  return { rows: rows.map(r => r.map(v => R.toNum(v))), J };
}
function subsample(rows, nmax) {
  if (rows.length <= nmax) return rows;
  const idx = [...Array(rows.length).keys()];
  for (let i = rows.length - 1; i > 0; i--) { const j = Math.floor(rng() * (i + 1)); [idx[i], idx[j]] = [idx[j], idx[i]]; }
  return idx.slice(0, nmax).map(i => rows[i]);
}
function scoreCols(rows, cols) {
  const n = rows.length, J = cols.length;
  const m = new Float64Array(n * J);
  for (let i = 0; i < n; i++) for (let j = 0; j < J; j++) m[i * J + j] = rows[i][cols[j]];
  const res = R.run(m, n, J, { iterations: ITERS, seed: 1 });
  return { z: res.z, ls: res.longstring };
}
const mean = v => v.reduce((s, x) => s + x, 0) / v.length;
const sd = v => { const m = mean(v); return Math.sqrt(v.reduce((s, x) => s + (x - m) * (x - m), 0) / (v.length - 1)); };
function pairedStats(a, b) {          // b - a, per respondent
  const d = a.map((v, i) => b[i] - v).filter(Number.isFinite);
  const m = mean(d), s = sd(d);
  const t = m / (s / Math.sqrt(d.length));
  return { delta: m, dz: m / s, t, n: d.length };
}

// ---- dataset inventory: only J >= 60 so each half has >= 30 items ----
const EXT = "ext_bench2";
const idx = R.parseCSV(fs.readFileSync(path.join(EXT, "index.csv"), "utf8"));
const NAME = idx.header.indexOf("name"), JJ = idx.header.indexOf("J");
const sets = [];
for (const r of idx.rows) {
  const J = Number(r[JJ]);
  if (J >= 60) sets.push({ tag: r[NAME], file: path.join(EXT, r[NAME] + "_M.csv"), header: false });
}
sets.push({ tag: "study1", file: "_study_matrix.csv", header: true });

console.log("FATIGUE TEST — first vs second half of the questionnaire (rc = permutation z)");
console.log("positional delta = rc(2nd half) - rc(1st half), within respondent; negative = worse coherence late");
console.log("random-split delta = same quantity for " + RAND_SPLITS + " random halvings (content control)\n");
console.log("dataset          n     J   | rc 1st  rc 2nd  posΔ    dz     | randΔ (mean±sd)  | excess | LS 1st LS 2nd");

const out = [["dataset","n","J","rc_first","rc_second","pos_delta","dz","rand_delta_mean","rand_delta_sd","excess","ls_first","ls_second","ls_pos_delta","ls_rand_delta","ls_excess"]];
for (const s of sets) {
  if (!fs.existsSync(s.file)) { console.log("  (missing) " + s.tag); continue; }
  let { rows, J } = readMatrix(s.file, s.header);
  rows = subsample(rows, NMAX);
  const half = Math.floor(J / 2);
  const first = [...Array(half).keys()];
  const second = [...Array(J - half).keys()].map(i => i + half);

  const A = scoreCols(rows, first), B = scoreCols(rows, second);
  const pos = pairedStats(A.z, B.z);
  const lsA = mean(A.ls), lsB = mean(B.ls);

  // random-split control: shuffle columns, take the same sizes
  const rdeltas = [], rlsdeltas = [];
  for (let r = 0; r < RAND_SPLITS; r++) {
    const perm = [...Array(J).keys()];
    for (let i = J - 1; i > 0; i--) { const j = Math.floor(rng() * (i + 1)); [perm[i], perm[j]] = [perm[j], perm[i]]; }
    const ra = perm.slice(0, half), rb = perm.slice(half);
    const X = scoreCols(rows, ra), Y = scoreCols(rows, rb);
    rdeltas.push(pairedStats(X.z, Y.z).delta);
    rlsdeltas.push(mean(Y.ls) - mean(X.ls));
  }
  const rm = mean(rdeltas), rs = sd(rdeltas);
  const rlm = mean(rlsdeltas), rls = sd(rlsdeltas);
  const excess = pos.delta - rm;
  const lsExcess = (lsB - lsA) - rlm;

  console.log(
    s.tag.padEnd(15) + String(rows.length).padStart(5) + String(J).padStart(5) + "   | " +
    mean(A.z).toFixed(2).padStart(6) + "  " + mean(B.z).toFixed(2).padStart(6) + "  " +
    (pos.delta >= 0 ? "+" : "") + pos.delta.toFixed(3) + "  " + (pos.dz >= 0 ? "+" : "") + pos.dz.toFixed(2) +
    "  | " + (rm >= 0 ? "+" : "") + rm.toFixed(3) + "±" + rs.toFixed(3) +
    "    | " + (excess >= 0 ? "+" : "") + excess.toFixed(3) +
    " | " + lsA.toFixed(1).padStart(5) + " " + lsB.toFixed(1).padStart(5) +
    "  LSposΔ " + ((lsB-lsA)>=0?"+":"") + (lsB-lsA).toFixed(2) + "  LSrandΔ " + (rlm>=0?"+":"") + rlm.toFixed(2) + "  LSexcess " + (lsExcess>=0?"+":"") + lsExcess.toFixed(2));
  out.push([s.tag, rows.length, J, mean(A.z).toFixed(4), mean(B.z).toFixed(4), pos.delta.toFixed(4),
            pos.dz.toFixed(4), rm.toFixed(4), rs.toFixed(4), excess.toFixed(4), lsA.toFixed(3), lsB.toFixed(3),
            (lsB-lsA).toFixed(3), rlm.toFixed(3), lsExcess.toFixed(3)]);
}
fs.writeFileSync("fatigue_split.csv", out.map(r => r.join(",")).join("\n"));
console.log("\nwrote fatigue_split.csv");
