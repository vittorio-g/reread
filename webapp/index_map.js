/* index_map.js — Nomological map of rc against the family of careless indices.
 * Per respondent, 9 indices ORIENTED so high = careless:
 *   rc, psych_synonyms, psych_antonyms, even_odd, rpr,
 *   person_total, mahalanobis_d2, longstring, irv
 * Five come oriented from the shipped engine (ensembleFeatures.oriented). Four
 * consistency indices computed here: psych synonyms/antonyms (Meade & Craig 2012),
 * even-odd, and RPR (resampled personal reliability, Goldammer 2024). Even-odd and
 * RPR need an item->scale map, built data-driven from the correlation structure
 * (top-k eigenvectors, item assigned to its dominant factor, keyed by loading sign)
 * so every dataset gets all 9 indices uniformly.
 * Correlations downstream are SPEARMAN, so only orientation sign matters (no z needed). */
const R = require("./site/reread.js");
const fs = require("fs");
const D = "../Dataset/gt_benchmark_candidates";

// In-envelope multi-construct batteries for the label-free nomological map (no GT needed,
// so broader than the AUC benchmark): the 10 external batteries of Table tab:ext + Study 1
// + four further public personality batteries (HEXACO-240, AMBI, MIES, FBPS).
const SETS = [
  ["Study1_induced",   "_study_matrix.csv",                     "in"],
  ["Kay_S2",           D+"/kay_idris_idria/kay_matrix_s2.csv",  "in"],
  ["Kay_S1",           D+"/kay_idris_idria/kay_matrix_s1.csv",  "in"],
  ["warning_IPIP300",  D+"/warning_ipipneo300/_matrix.csv",     "in"],
  ["opsy_hexaco",      D+"/opsy_hexaco_matrix.csv",             "in"],
  ["opsy_16pf",        D+"/opsy_16pf_matrix.csv",               "in"],
  ["ambi_2019",        D+"/ambi_2019_matrix.csv",               "in"],
  ["smarvus",          D+"/smarvus/smarvus_matrix.csv",         "in"],
  ["mies_dev",         D+"/mies_dev_matrix.csv",                "in"],
  ["Kay_S6",           D+"/kay_idris_idria/kay_matrix_s6.csv",  "border"],
  ["Kay_S5",           D+"/kay_idris_idria/kay_matrix_s5.csv",  "border"],
  ["fbps",             D+"/fbps_validation_matrix.csv",         "border"],
  ["duckworth",        D+"/duckworth_grit_vcl/duckworth_matrix.csv", "border"],
  ["krause",           D+"/krause_ier_youth/krause_matrix.csv", "border"],
  ["douglas",          D+"/douglas2023_dataquality/douglas_matrix.csv", "border"],
];

const SUBSAMPLE = 4000;
function readMatrix(path) {
  const P = R.parseCSV(fs.readFileSync(path, "utf8")), J = P.header.length;
  let rows = P.rows.map(r => r.map(v => { const x = Number(v); return Number.isFinite(x) ? x : NaN; }));
  let keep = rows.map((_, i) => i);
  if (rows.length > SUBSAMPLE) {
    const rng = R._rng(42); const idx = rows.map((_, i) => i);
    for (let i = idx.length - 1; i > 0; i--) { const j = Math.floor(rng() * (i + 1)); const t = idx[i]; idx[i] = idx[j]; idx[j] = t; }
    keep = idx.slice(0, SUBSAMPLE); rows = keep.map(i => rows[i]);
  }
  return { rows, J, keep };
}
// build scales + keying sign from a real item->scale key CSV (one scale-id per column, -1 excluded)
function scalesFromKey(path, Xz, J, n) {
  const ids = R.parseCSV(fs.readFileSync(path, "utf8")).rows.map(r => parseInt(r[0], 10));
  const k = Math.max(0, ...ids) + 1;
  const scales = Array.from({ length: k }, () => []);
  for (let j = 0; j < J && j < ids.length; j++) if (ids[j] >= 0) scales[ids[j]].push(j);
  const sign = new Float64Array(J).fill(1);
  for (const items of scales) {
    if (items.length < 2) continue;
    const sm = new Float64Array(n);
    for (let i = 0; i < n; i++) { let s = 0; for (const j of items) s += Xz[i * J + j]; sm[i] = s / items.length; }
    let mm = 0; for (let i = 0; i < n; i++) mm += sm[i]; mm /= n;
    for (const j of items) {
      let sab = 0; for (let i = 0; i < n; i++) sab += Xz[i * J + j] * (sm[i] - mm);
      sign[j] = sab >= 0 ? 1 : -1;
    }
  }
  return { scales, sign };
}
function colMeans(rows, J) {
  const m = new Float64Array(J), c = new Int32Array(J);
  for (const r of rows) for (let j = 0; j < J; j++) if (Number.isFinite(r[j])) { m[j] += r[j]; c[j]++; }
  for (let j = 0; j < J; j++) m[j] = c[j] ? m[j] / c[j] : 0;
  return m;
}
// standardized (mean-imputed) data + correlation matrix
function prep(rows, J, mean) {
  const n = rows.length, X = new Float64Array(n * J);
  for (let i = 0; i < n; i++) for (let j = 0; j < J; j++) { const v = rows[i][j]; X[i * J + j] = Number.isFinite(v) ? v : mean[j]; }
  const mu = new Float64Array(J), sd = new Float64Array(J);
  for (let j = 0; j < J; j++) { let s = 0; for (let i = 0; i < n; i++) s += X[i * J + j]; mu[j] = s / n; }
  for (let j = 0; j < J; j++) { let s = 0; for (let i = 0; i < n; i++) { const d = X[i * J + j] - mu[j]; s += d * d; } sd[j] = Math.sqrt(s / n) || 1e-9; }
  const Xz = new Float64Array(n * J);
  for (let i = 0; i < n; i++) for (let j = 0; j < J; j++) Xz[i * J + j] = (X[i * J + j] - mu[j]) / sd[j];
  const C = new Float64Array(J * J);
  for (let a = 0; a < J; a++) for (let b = a; b < J; b++) {
    let s = 0; for (let i = 0; i < n; i++) s += Xz[i * J + a] * Xz[i * J + b];
    const r = s / n; C[a * J + b] = r; C[b * J + a] = r;
  }
  return { Xz, C };
}
// top-k eigenvectors of C via power iteration + deflation
function topEigvecs(C, J, k) {
  const Cw = Float64Array.from(C), vecs = [];
  for (let f = 0; f < k; f++) {
    const rng = R._rng(1000 + f); let v = new Float64Array(J);
    for (let j = 0; j < J; j++) v[j] = rng() - 0.5;
    let nv = Math.hypot(...v) || 1; for (let j = 0; j < J; j++) v[j] /= nv;
    let lam = 0;
    for (let it = 0; it < 60; it++) {
      const w = new Float64Array(J);
      for (let a = 0; a < J; a++) { let s = 0; const base = a * J; for (let b = 0; b < J; b++) s += Cw[base + b] * v[b]; w[a] = s; }
      let norm = 0; for (let j = 0; j < J; j++) norm += w[j] * w[j]; norm = Math.sqrt(norm) || 1;
      for (let j = 0; j < J; j++) w[j] /= norm; v = w; lam = norm;
    }
    vecs.push(v);
    for (let a = 0; a < J; a++) { const base = a * J; for (let b = 0; b < J; b++) Cw[base + b] -= lam * v[a] * v[b]; }
  }
  return vecs;
}
// assign each item to its dominant factor; keying sign from that loading
function assignFactors(vecs, J) {
  const k = vecs.length, fac = new Int32Array(J), sign = new Float64Array(J);
  for (let j = 0; j < J; j++) {
    let best = 0, bf = 0; for (let f = 0; f < k; f++) { const a = Math.abs(vecs[f][j]); if (a > best) { best = a; bf = f; } }
    fac[j] = bf; sign[j] = vecs[bf][j] >= 0 ? 1 : -1;
  }
  // scale -> item index list
  const scales = Array.from({ length: k }, () => []);
  for (let j = 0; j < J; j++) scales[fac[j]].push(j);
  return { scales, sign };
}
// per-person consistency = cor over scales of (half-A score, half-B score)
function halfScoreConsistency(Xz, J, scales, sign, n, splitter) {
  const out = new Array(n).fill(NaN);
  const usable = [];
  for (const items of scales) if (items.length >= 4) usable.push(items);
  if (usable.length < 3) return out;
  for (let i = 0; i < n; i++) {
    const vA = [], vB = [];
    for (const items of usable) {
      const [A, B] = splitter(items);
      if (A.length < 1 || B.length < 1) continue;
      let a = 0, b = 0; for (const j of A) a += sign[j] * Xz[i * J + j]; for (const j of B) b += sign[j] * Xz[i * J + j];
      vA.push(a / A.length); vB.push(b / B.length);
    }
    if (vA.length < 3) continue;
    const m = arr => arr.reduce((s, v) => s + v, 0) / arr.length;
    const mA = m(vA), mB = m(vB); let sab = 0, saa = 0, sbb = 0;
    for (let t = 0; t < vA.length; t++) { const da = vA[t] - mA, db = vB[t] - mB; sab += da * db; saa += da * da; sbb += db * db; }
    out[i] = (saa > 0 && sbb > 0) ? sab / Math.sqrt(saa * sbb) : 0;
  }
  return out;
}
const evenOddSplit = items => {
  const A = [], B = []; for (let t = 0; t < items.length; t++) (t % 2 ? B : A).push(items[t]); return [A, B];
};
function rprScore(Xz, J, scales, sign, n, B, seed) {
  const acc = new Float64Array(n), cnt = new Int32Array(n);
  for (let b = 0; b < B; b++) {
    const rng = R._rng(seed + b);
    const splitter = items => {
      const sh = items.slice(); for (let i = sh.length - 1; i > 0; i--) { const j = Math.floor(rng() * (i + 1)); const t = sh[i]; sh[i] = sh[j]; sh[j] = t; }
      const h = Math.floor(sh.length / 2); return [sh.slice(0, h), sh.slice(h)];
    };
    const s = halfScoreConsistency(Xz, J, scales, sign, n, splitter);
    for (let i = 0; i < n; i++) if (Number.isFinite(s[i])) { acc[i] += s[i]; cnt[i]++; }
  }
  const out = new Array(n).fill(NaN);
  for (let i = 0; i < n; i++) if (cnt[i] > 0) out[i] = acc[i] / cnt[i];
  return out;
}
// psych synonyms/antonyms pair selection + within-person cor (from v1)
function selectPairs(C, J, thr, minPairs) {
  const off = [];
  for (let a = 0; a < J; a++) for (let b = a + 1; b < J; b++) off.push([a, b, C[a * J + b]]);
  let syn = off.filter(p => p[2] >= thr).map(p => [p[0], p[1]]);
  let anti = off.filter(p => p[2] <= -thr).map(p => [p[0], p[1]]);
  if (syn.length < minPairs) { const s = off.slice().sort((p, q) => q[2] - p[2]); syn = s.slice(0, minPairs).filter(p => p[2] > 0).map(p => [p[0], p[1]]); }
  if (anti.length < 8) { const s = off.slice().sort((p, q) => p[2] - q[2]); anti = s.slice(0, minPairs).filter(p => p[2] < 0).map(p => [p[0], p[1]]); }
  return { syn, anti };
}
function personPairCor(Xz, J, pairs, n) {
  const out = new Array(n).fill(NaN); if (pairs.length < 3) return out;
  for (let i = 0; i < n; i++) {
    const xs = [], ys = []; for (const [a, b] of pairs) { xs.push(Xz[i * J + a]); ys.push(Xz[i * J + b]); }
    const m = arr => arr.reduce((s, v) => s + v, 0) / arr.length; const mx = m(xs), my = m(ys);
    let sxy = 0, sxx = 0, syy = 0; for (let t = 0; t < xs.length; t++) { const dx = xs[t] - mx, dy = ys[t] - my; sxy += dx * dy; sxx += dx * dx; syy += dy * dy; }
    out[i] = (sxx > 0 && syy > 0) ? sxy / Math.sqrt(sxx * syy) : 0;
  }
  return out;
}

const COLS = ["rc", "synonyms", "antonyms", "even_odd", "rpr", "person_total", "d2", "longstring", "irv", "response_time"];
const lines = ["dataset," + COLS.join(",")];
console.log("dataset          n     J    k  #syn #anti #scl  scaleKey  RT");
for (const [name, mp] of SETS) {
  let M; try { M = readMatrix(mp); } catch (e) { console.log(name.padEnd(16) + " MISSING"); continue; }
  const { rows, J, keep } = M, n = rows.length;
  const flat = new Float64Array(n * J);
  for (let i = 0; i < n; i++) for (let j = 0; j < J; j++) flat[i * J + j] = rows[i][j];
  const ef = R.ensembleFeatures(flat, n, J, { iterations: 200, seed: 1 });
  const O = ef.oriented;
  const mean = colMeans(rows, J);
  const { Xz, C } = prep(rows, J, mean);
  // synonyms / antonyms
  const { syn, anti } = selectPairs(C, J, 0.40, 20);
  const oSyn = personPairCor(Xz, J, syn, n).map(v => Number.isFinite(v) ? -v : NaN);
  const oAnti = personPairCor(Xz, J, anti, n).map(v => Number.isFinite(v) ? v : NaN);
  // scales for even-odd + RPR: real key if available, else data-driven eigenvectors
  const keyPath = "keys/" + name + "_scale.csv";
  let scales, sign, scaleSrc, k;
  if (fs.existsSync(keyPath)) { ({ scales, sign } = scalesFromKey(keyPath, Xz, J, n)); scaleSrc = "real"; k = scales.length; }
  else { k = Math.max(3, Math.min(20, Math.round(J / 10))); ({ scales, sign } = assignFactors(topEigvecs(C, J, k), J)); scaleSrc = "data"; }
  const nUsable = scales.filter(s => s.length >= 4).length;
  const eo = halfScoreConsistency(Xz, J, scales, sign, n, evenOddSplit).map(v => Number.isFinite(v) ? -v : NaN);
  const rpr = rprScore(Xz, J, scales, sign, n, 20, 7).map(v => Number.isFinite(v) ? -v : NaN);
  // response time (oriented careless = fast => -rt), aligned via keep indices
  const rtPath = "keys/" + name + "_rt.csv";
  let ort = new Array(n).fill(NaN), rtSrc = "-";
  if (fs.existsSync(rtPath)) {
    const rtFull = R.parseCSV(fs.readFileSync(rtPath, "utf8")).rows.map(r => Number(r[0]));
    ort = keep.map(i => (Number.isFinite(rtFull[i]) && rtFull[i] > 0 ? -rtFull[i] : NaN)); rtSrc = "yes";
  }
  for (let i = 0; i < n; i++) {
    const rec = [name, O.rr[i], oSyn[i], oAnti[i], eo[i], rpr[i], O.person_total[i], O.d2[i], O.longstring[i], O.irv[i], ort[i]];
    lines.push(rec.map(v => (typeof v === "number" ? (Number.isFinite(v) ? v.toFixed(5) : "NA") : v)).join(","));
  }
  console.log(name.padEnd(16) + String(n).padStart(5) + String(J).padStart(5) + String(k).padStart(5) +
    String(syn.length).padStart(5) + String(anti.length).padStart(6) + String(nUsable).padStart(5) +
    scaleSrc.padStart(9) + rtSrc.padStart(5));
}
fs.writeFileSync("index_map_long.csv", lines.join("\n"));
console.log("\nwrote index_map_long.csv (" + (lines.length - 1) + " rows, 10 indices)");
