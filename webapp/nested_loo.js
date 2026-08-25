/* nested_loo.js — quantify residual leakage in the LOO of the CaReReRe ensemble.
 *
 * Reviewer point (#6): the reported LOO refits only the logistic COMBINER on n-1;
 * the FEATURES are computed once on the full sample, so the held-out respondent
 * leaks into (i) the pair-selection correlation matrix, (ii) the sign-alignment,
 * (iii) the Person-Total reference mean profile, and (iv) the robust median/MAD
 * standardization. The permutation null and LongString/IRV are per-respondent and
 * do not leak.
 *
 * This script computes TWO leave-one-out AUCs on the induced-careless study
 * (n=157, 109 items), differing ONLY in where the sample-level quantities come
 * from:
 *   (A) combiner-only LOO  — features from the FULL sample (the reported setup).
 *   (B) full nested   LOO  — for every fold, re-derive pair-selection, signs,
 *                            Person-Total reference, and median/MAD on the n-1
 *                            training rows, then score the held-out respondent.
 * The gap A - B is the optimistic bias the reviewer worries about.
 *
 * The feature engine is reimplemented here VERBATIM from webapp/site/reread.js and
 * validated (VALIDATION block) against the shipped R.run() on the full sample, so
 * the numbers are the shipped method's, not an approximation.
 */
const fs = require("fs");
const path = require("path");
const R = require("./site/reread.js");   // shipped lib: parseCSV, detectColumns, buildMatrix, run, _rng

/* ---- verbatim internals copied from reread.js (kept identical) ---- */
function colStats(M, n, J) {
  const max = new Float64Array(J).fill(-Infinity), min = new Float64Array(J).fill(Infinity), median = new Float64Array(J), buf = [];
  for (let j = 0; j < J; j++) {
    buf.length = 0;
    for (let i = 0; i < n; i++) { const v = M[i * J + j]; if (Number.isFinite(v)) { buf.push(v); if (v > max[j]) max[j] = v; if (v < min[j]) min[j] = v; } }
    buf.sort((a, b) => a - b);
    median[j] = buf.length ? buf[Math.floor(buf.length / 2)] : 0;
    if (!Number.isFinite(max[j])) { max[j] = 1; min[j] = 0; }
    if (max[j] === 0) max[j] = 1;
  }
  return { max, min, median };
}
// prepare rows `rowsIdx` of the full matrix M into rescaled/imputed P, using
// column stats `cs` (max/median) that were computed on the TRAINING rows.
function prepareWith(M, J, rowsIdx, cs) {
  const n = rowsIdx.length, P = new Float64Array(n * J);
  for (let a = 0; a < n; a++) { const i = rowsIdx[a];
    for (let j = 0; j < J; j++) { let v = M[i * J + j]; if (!Number.isFinite(v)) v = cs.median[j]; P[a * J + j] = v / cs.max[j]; } }
  const minP = new Float64Array(J); for (let j = 0; j < J; j++) minP[j] = cs.min[j] / cs.max[j];
  return { P, minP };
}
function correlationPairs(P, n, J) {
  const mean = new Float64Array(J), sd = new Float64Array(J);
  for (let j = 0; j < J; j++) { let s = 0; for (let i = 0; i < n; i++) s += P[i * J + j]; mean[j] = s / n;
    let q = 0; for (let i = 0; i < n; i++) { const d = P[i * J + j] - mean[j]; q += d * d; } sd[j] = Math.sqrt(q / (n - 1)) || 0; }
  const Z = new Float64Array(n * J);
  for (let i = 0; i < n; i++) for (let j = 0; j < J; j++) Z[i * J + j] = sd[j] > 0 ? (P[i * J + j] - mean[j]) / sd[j] : 0;
  const nPairs = (J * (J - 1)) / 2, pa = new Int32Array(nPairs), pb = new Int32Array(nPairs), pr = new Float64Array(nPairs);
  let p = 0;
  for (let a = 0; a < J; a++) for (let b = a + 1; b < J; b++, p++) { let s = 0; for (let i = 0; i < n; i++) s += Z[i * J + a] * Z[i * J + b]; pa[p] = a; pb[p] = b; pr[p] = (sd[a] > 0 && sd[b] > 0) ? s / (n - 1) : 0; }
  return { pa, pb, pr, nPairs };
}
function indCors(P, minP, n, J, pa, pb, pr, idx) {  // swap-invariant ICC (matches reread.js)
  const k = idx.length, out = new Float64Array(n);
  for (let i = 0; i < n; i++) {
    let sAB = 0, sSq = 0, sSum = 0;
    for (let t = 0; t < k; t++) { const p = idx[t]; const a = P[i * J + pa[p]]; const bj = pb[p], bv = P[i * J + bj]; const b = pr[p] < 0 ? (1 + minP[bj]) - bv : bv; sAB += a * b; sSq += a * a + b * b; sSum += a + b; }
    const m = sSum / (2 * k); const num = sAB / k - m * m; const den = sSq / (2 * k) - m * m;
    out[i] = den > 1e-12 ? Math.abs(num / den) : 0;
  }
  return out;
}
function sampleIdx(nPairs, k, rand) {
  const pool = new Int32Array(nPairs); for (let i = 0; i < nPairs; i++) pool[i] = i;
  const out = new Int32Array(k);
  for (let t = 0; t < k; t++) { const r = t + Math.floor(rand() * (nPairs - t)); const tmp = pool[t]; pool[t] = pool[r]; pool[r] = tmp; out[t] = pool[t]; }
  return out;
}
// per-respondent longstring + Person-Total against a GIVEN reference mean profile cm.
function auxWithRef(P, n, J, cm) {
  const personTotal = new Float64Array(n), longstring = new Int32Array(n);
  let mc = 0; for (let j = 0; j < J; j++) mc += cm[j]; mc /= J;
  for (let i = 0; i < n; i++) {
    let m = 0; for (let j = 0; j < J; j++) m += P[i * J + j]; m /= J;
    let sab = 0, sa = 0, sb = 0;
    for (let j = 0; j < J; j++) { const da = P[i * J + j] - m, db = cm[j] - mc; sab += da * db; sa += da * da; sb += db * db; }
    personTotal[i] = (sa > 0 && sb > 0) ? sab / Math.sqrt(sa * sb) : 0;
    let best = 1, run = 1; for (let j = 1; j < J; j++) { if (P[i * J + j] === P[i * J + j - 1]) { run++; if (run > best) best = run; } else run = 1; }
    longstring[i] = best;
  }
  return { personTotal, longstring };
}

/* ---- raw features (rr, longstring, person_total) for ALL rows, with all
 *      sample-level quantities derived from `trainIdx`. Mirrors run()+auxiliaries. */
function rawFeatures(M, n, J, trainIdx, allIdx, iters, seed) {
  const cs = colStats(subMatrix(M, J, trainIdx), trainIdx.length, J);      // train max/min/median
  const { P: Pall, minP } = prepareWith(M, J, allIdx, cs);                  // rescale/impute ALL with train stats
  // train slice of Pall (positions of trainIdx within allIdx == same order since allIdx=0..n-1)
  const Ptr = subRows(Pall, J, trainIdx);
  const { pa, pb, pr, nPairs } = correlationPairs(Ptr, trainIdx.length, J); // train sample correlations
  const k = Math.min(Math.max(15, Math.round(0.03 * nPairs)), nPairs);
  const order = Array.from({ length: nPairs }, (_, i) => i).sort((x, y) => Math.abs(pr[y]) - Math.abs(pr[x]));
  const topIdx = Int32Array.from(order.slice(0, k));
  const coupled = indCors(Pall, minP, n, J, pa, pb, pr, topIdx);
  const rand = R._rng(seed), sum = new Float64Array(n), sumSq = new Float64Array(n);
  for (let b = 0; b < iters; b++) { const rc = indCors(Pall, minP, n, J, pa, pb, pr, sampleIdx(nPairs, k, rand)); for (let i = 0; i < n; i++) { sum[i] += rc[i]; sumSq[i] += rc[i] * rc[i]; } }
  const z = new Float64Array(n); let minValid = Infinity;
  for (let i = 0; i < n; i++) { const mu = sum[i] / iters, va = sumSq[i] / iters - mu * mu, sd = Math.sqrt(Math.max(va, 0)); z[i] = sd > 1e-12 ? (coupled[i] - mu) / sd : NaN; if (Number.isFinite(z[i]) && z[i] < minValid) minValid = z[i]; }
  for (let i = 0; i < n; i++) if (!Number.isFinite(z[i])) z[i] = minValid;
  // Person-Total reference = TRAIN mean profile
  const cm = new Float64Array(J); for (const i of trainIdx) for (let j = 0; j < J; j++) cm[j] += Pall[i * J + j]; for (let j = 0; j < J; j++) cm[j] /= trainIdx.length;
  const aux = auxWithRef(Pall, n, J, cm);
  return { rr: z, longstring: aux.longstring, person_total: aux.personTotal };
}
function subMatrix(M, J, idx) { const m = new Float64Array(idx.length * J); for (let a = 0; a < idx.length; a++) for (let j = 0; j < J; j++) m[a * J + j] = M[idx[a] * J + j]; return m; }
function subRows(P, J, idx) { const m = new Float64Array(idx.length * J); for (let a = 0; a < idx.length; a++) for (let j = 0; j < J; j++) m[a * J + j] = P[idx[a] * J + j]; return m; }

/* ---- robust within-dataset z with med/MAD taken from a TRAINING slice ---- */
function robustParams(x, trainIdx) {
  const v = trainIdx.map(i => x[i]).sort((a, b) => a - b), med = v[Math.floor(v.length / 2)];
  const dev = trainIdx.map(i => Math.abs(x[i] - med)).sort((a, b) => a - b);
  let mad = dev[Math.floor(dev.length / 2)] * 1.4826;
  if (!(mad > 1e-9)) { let m = 0; for (const i of trainIdx) m += x[i]; m /= trainIdx.length; let q = 0; for (const i of trainIdx) { const d = x[i] - m; q += d * d; } const sd = Math.sqrt(q / Math.max(trainIdx.length - 1, 1)) || 1; return { c: m, s: sd }; }
  return { c: med, s: mad };
}
// oriented (higher = more careless) feature rows for ALL, standardized with TRAIN med/MAD
function orient(raw, trainIdx, n) {
  const pRR = robustParams(raw.rr, trainIdx), pLS = robustParams(raw.longstring, trainIdx), pPT = robustParams(raw.person_total, trainIdx);
  const O = [];
  for (let i = 0; i < n; i++) O.push([
    1,
    -((raw.rr[i] - pRR.c) / pRR.s),                 // low rr = careless -> negate
    (raw.longstring[i] - pLS.c) / pLS.s,            // high = careless
    -((raw.person_total[i] - pPT.c) / pPT.s)        // low = careless -> negate
  ]);
  return O;
}

/* ---- logistic fit (IRLS, ridge lambda=1 on non-intercept) — from fit_triad.js ---- */
function solve(A, b) { const m = b.length, Aug = A.map((r, i) => r.concat(b[i]));
  for (let c = 0; c < m; c++) { let pv = c; for (let r = c + 1; r < m; r++) if (Math.abs(Aug[r][c]) > Math.abs(Aug[pv][c])) pv = r; [Aug[c], Aug[pv]] = [Aug[pv], Aug[c]]; const d = Aug[c][c] || 1e-12;
    for (let r = 0; r < m; r++) if (r !== c) { const f = Aug[r][c] / d; for (let kk = c; kk <= m; kk++) Aug[r][kk] -= f * Aug[c][kk]; } }
  return Aug.map((r, i) => r[m] / (r[i] || 1e-12)); }
function fitLogit(rows, yy) { const P = rows[0].length; let b = new Array(P).fill(0);
  for (let it = 0; it < 50; it++) { const g = new Array(P).fill(0), H = Array.from({ length: P }, () => new Array(P).fill(0));
    for (let i = 0; i < rows.length; i++) { let e = 0; for (let kk = 0; kk < P; kk++) e += b[kk] * rows[i][kk]; const p = 1 / (1 + Math.exp(-e)), w = Math.max(p * (1 - p), 1e-6);
      for (let kk = 0; kk < P; kk++) { g[kk] += (yy[i] - p) * rows[i][kk]; for (let l = 0; l < P; l++) H[kk][l] += w * rows[i][kk] * rows[i][l]; } }
    for (let kk = 1; kk < P; kk++) { g[kk] -= 1.0 * b[kk]; H[kk][kk] += 1.0; } const s = solve(H, g); let mx = 0; for (let kk = 0; kk < P; kk++) { b[kk] += s[kk]; mx = Math.max(mx, Math.abs(s[kk])); } if (mx < 1e-8) break; }
  return b; }
const predict = (b, x) => 1 / (1 + Math.exp(-(b[0] * x[0] + b[1] * x[1] + b[2] * x[2] + b[3] * x[3])));

/* ---- metrics ---- */
function auc(sc, y) { const idx = [...sc.keys()].sort((a, b) => sc[a] - sc[b]); const rk = []; idx.forEach((id, r) => rk[id] = r + 1);
  let n1 = 0, n0 = 0, s = 0; for (let i = 0; i < sc.length; i++) { if (y[i] === 1) { n1++; s += rk[i]; } else n0++; } return n1 && n0 ? (s - n1 * (n1 + 1) / 2) / (n1 * n0) : NaN; }
function bestMCC(sc, y) { const ts = [...new Set(sc)].sort((a, b) => a - b); let best = -2, bt = 0;
  for (const t of ts) { let tp = 0, fp = 0, fn = 0, tn = 0; for (let i = 0; i < y.length; i++) { const f = sc[i] >= t; if (y[i] === 1) f ? tp++ : fn++; else f ? fp++ : tn++; } const d = Math.sqrt((tp + fp) * (tp + fn) * (tn + fp) * (tn + fn)), m = d ? (tp * tn - fp * fn) / d : 0; if (m > best) { best = m; bt = t; } } return { mcc: best, thr: bt }; }
const pearson = (a, b) => { const n = a.length; let ma = 0, mb = 0; for (let i = 0; i < n; i++) { ma += a[i]; mb += b[i]; } ma /= n; mb /= n; let sab = 0, sa = 0, sb = 0; for (let i = 0; i < n; i++) { const da = a[i] - ma, db = b[i] - mb; sab += da * db; sa += da * da; sb += db * db; } return sab / Math.sqrt(sa * sb); };

/* ==================== run ==================== */
const ITERS = 200, SEED = 1;
const raw = fs.readFileSync(path.join(__dirname, "_study_matrix.csv"), "utf8");
const parsed = R.parseCSV(raw);
const det = R.detectColumns(parsed.header, parsed.rows);
const { M, n, J } = R.buildMatrix(parsed.rows, det.itemCols);
const y = R.parseCSV(fs.readFileSync(path.join(__dirname, "_study_labels.csv"), "utf8")).rows.map(r => Number(r[0]));
if (y.length !== n) throw new Error("row mismatch " + n + " vs " + y.length);
const ALL = [...Array(n).keys()];
console.log("Study: n=" + n + ", J=" + J + " items, careless=" + y.reduce((s, v) => s + v, 0) + "/" + n + "  (iters=" + ITERS + ")\n");

/* --- VALIDATION: our full-sample rr must equal the shipped R.run() rr --- */
const rawFull = rawFeatures(M, n, J, ALL, ALL, ITERS, SEED);
const shipped = R.run(M, n, J, { corProp: 0.03, iterations: ITERS, minPairs: 15, seed: SEED });
let maxdiff = 0; for (let i = 0; i < n; i++) maxdiff = Math.max(maxdiff, Math.abs(rawFull.rr[i] - shipped.z[i]));
console.log("VALIDATION vs shipped R.run():  Pearson(rr)=" + pearson(rawFull.rr, shipped.z).toFixed(6) +
  "  max|Δrr|=" + maxdiff.toExponential(2) + (maxdiff < 1e-9 ? "  -> bit-identical engine ✓" : "  <- CHECK"));

/* --- (A) combiner-only LOO: features from FULL sample, refit weights on n-1 --- */
const O_full = orient(rawFull, ALL, n);
const pA = new Array(n);
for (let i = 0; i < n; i++) { const tr = ALL.filter(j => j !== i); const b = fitLogit(tr.map(j => O_full[j]), tr.map(j => y[j])); pA[i] = predict(b, O_full[i]); }
const aucA = auc(pA, y), mccA = bestMCC(pA, y);

/* --- (B) full nested LOO: re-derive ALL sample quantities on n-1 each fold --- */
const pB = new Array(n);
for (let i = 0; i < n; i++) {
  const tr = ALL.filter(j => j !== i);
  const rawTr = rawFeatures(M, n, J, tr, ALL, ITERS, SEED);   // sample quantities from train, features for all
  const O = orient(rawTr, tr, n);                              // standardize with train med/MAD
  const b = fitLogit(tr.map(j => O[j]), tr.map(j => y[j]));
  pB[i] = predict(b, O[i]);
  if ((i + 1) % 40 === 0) process.stdout.write("  nested fold " + (i + 1) + "/" + n + "\r");
}
const aucB = auc(pB, y), mccB = bestMCC(pB, y);

console.log("\n\n================ RESULTS (leave-one-out) ================");
console.log("(A) combiner-only LOO  (features from full sample) : AUC=" + aucA.toFixed(3) + "  MCC*=" + mccA.mcc.toFixed(3));
console.log("(B) full nested   LOO  (everything refit on n-1)    : AUC=" + aucB.toFixed(3) + "  MCC*=" + mccB.mcc.toFixed(3));
console.log("---------------------------------------------------------");
console.log("Optimistic bias  ΔAUC = " + (aucA - aucB).toFixed(3) + "   ΔMCC = " + (mccA.mcc - mccB.mcc).toFixed(3));
console.log("(MCC* = best-threshold MCC; A reproduces the reported combiner-only LOO.)");
fs.writeFileSync(path.join(__dirname, "nested_loo_results.json"), JSON.stringify({
  n, J, iters: ITERS, validation_maxdiff: maxdiff,
  A_combiner_only: { auc: aucA, mcc: mccA.mcc }, B_full_nested: { auc: aucB, mcc: mccB.mcc },
  delta: { auc: aucA - aucB, mcc: mccA.mcc - mccB.mcc }
}, null, 2));
