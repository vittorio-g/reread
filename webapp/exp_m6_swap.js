/* exp_m6_swap.js — Reviewer M6: is rc well-defined per respondent?
 * The shipped coherence is |cor(A,B)| across the k coupled pairs, A=left member (raw
 * proportion), B=sign-aligned right member. Swapping which member is "left" vs "right"
 * for a single pair changes the means of A and B and hence cor(A,B) -> rc depends on the
 * (arbitrary) column order. We (1) quantify how much rc and the flag set move under random
 * per-pair orientation flips and full column permutations, and (2) test a swap-INVARIANT
 * coherence: sign-aligned mean cross-product of standardized items, mean_t(zA_t * zB_t),
 * which is symmetric to swapping and to column order. Compare stability AND validity. */
const R = require("./site/rerere.js"); const fs = require("fs");
const auc = (s, y) => { let po = [], ne = []; for (let i = 0; i < y.length; i++)(y[i] ? po : ne).push(s[i]);
  if (!po.length || !ne.length) return NaN; let c = 0; for (const a of po) for (const b of ne) c += a > b ? 1 : a === b ? 0.5 : 0; return c / (po.length * ne.length); };
const spearman = (x, y) => { const n = x.length, rx = rank(x), ry = rank(y); let sx = 0, sy = 0, sxx = 0, syy = 0, sxy = 0;
  for (let i = 0; i < n; i++) { sx += rx[i]; sy += ry[i]; } sx /= n; sy /= n;
  for (let i = 0; i < n; i++) { const dx = rx[i] - sx, dy = ry[i] - sy; sxx += dx * dx; syy += dy * dy; sxy += dx * dy; } return sxy / Math.sqrt(sxx * syy); };
function rank(a) { const idx = a.map((v, i) => [v, i]).sort((p, q) => p[0] - q[0]); const r = new Array(a.length);
  for (let k = 0; k < idx.length;) { let j = k; while (j < idx.length && idx[j][0] === idx[k][0]) j++; const avg = (k + j - 1) / 2 + 1; for (let t = k; t < j; t++) r[idx[t][1]] = avg; k = j; } return r; }
const sd = a => { const m = a.reduce((s, v) => s + v, 0) / a.length; return Math.sqrt(a.reduce((s, v) => s + (v - m) * (v - m), 0) / a.length); };

// ---- faithful pipeline pieces (mirror rerere.js) ----
function prepare(M, n, J) {
  const max = new Float64Array(J), min = new Float64Array(J), med = new Float64Array(J);
  for (let j = 0; j < J; j++) { const col = []; for (let i = 0; i < n; i++) { const v = M[i * J + j]; if (Number.isFinite(v)) col.push(v); }
    col.sort((a, b) => a - b); max[j] = col[col.length - 1]; min[j] = col[0]; med[j] = col[Math.floor(col.length / 2)]; }
  const P = new Float64Array(n * J);
  for (let i = 0; i < n; i++) for (let j = 0; j < J; j++) { let v = M[i * J + j]; if (!Number.isFinite(v)) v = med[j]; P[i * J + j] = v / max[j]; }
  const minP = new Float64Array(J); for (let j = 0; j < J; j++) minP[j] = min[j] / max[j];
  return { P, minP };
}
function corPairs(P, n, J) {
  const mean = new Float64Array(J), sdv = new Float64Array(J);
  for (let j = 0; j < J; j++) { let s = 0; for (let i = 0; i < n; i++) s += P[i * J + j]; mean[j] = s / n;
    let q = 0; for (let i = 0; i < n; i++) { const d = P[i * J + j] - mean[j]; q += d * d; } sdv[j] = Math.sqrt(q / (n - 1)) || 0; }
  const Z = new Float64Array(n * J);
  for (let i = 0; i < n; i++) for (let j = 0; j < J; j++) Z[i * J + j] = sdv[j] > 0 ? (P[i * J + j] - mean[j]) / sdv[j] : 0;
  const nPairs = J * (J - 1) / 2; const pa = new Int32Array(nPairs), pb = new Int32Array(nPairs), pr = new Float64Array(nPairs);
  let p = 0; for (let a = 0; a < J; a++) for (let b = a + 1; b < J; b++, p++) { let s = 0; for (let i = 0; i < n; i++) s += Z[i * J + a] * Z[i * J + b];
    pa[p] = a; pb[p] = b; pr[p] = (sdv[a] > 0 && sdv[b] > 0) ? s / (n - 1) : 0; }
  return { Z, pa, pb, pr, nPairs };
}
// SHIPPED: |cor(A,B)| across pairs, A=P[pa], B=reverse-coded P[pb]. `flip` toggles pa/pb per pair.
function cohCor(P, minP, n, J, pa, pb, pr, idx, flip) {
  const k = idx.length, out = new Float64Array(n), A = new Float64Array(k), B = new Float64Array(k);
  for (let i = 0; i < n; i++) {
    for (let t = 0; t < k; t++) { const p = idx[t]; let ai = pa[p], bi = pb[p]; if (flip && flip[t]) { const tmp = ai; ai = bi; bi = tmp; }
      A[t] = P[i * J + ai]; const bv = P[i * J + bi]; B[t] = pr[p] < 0 ? (1 + minP[bi]) - bv : bv; }
    let ma = 0, mb = 0; for (let t = 0; t < k; t++) { ma += A[t]; mb += B[t]; } ma /= k; mb /= k;
    let sab = 0, sa = 0, sb = 0; for (let t = 0; t < k; t++) { const da = A[t] - ma, db = B[t] - mb; sab += da * db; sa += da * da; sb += db * db; }
    out[i] = (sa > 0 && sb > 0) ? Math.abs(sab / Math.sqrt(sa * sb)) : 0;
  }
  return out;
}
// SWAP-INVARIANT (proper): intraclass-style correlation for exchangeable pairs =
// Pearson on the doubled data (a,b)+(b,a). Uses the SAME per-respondent raw values as cor,
// just symmetrized -> |num/den| with num=mean(a*b)-m^2, den=mean((a^2+b^2)/2)-m^2, m=mean of all members.
function cohSI(P, minP, n, J, pa, pb, pr, idx) {
  const k = idx.length, out = new Float64Array(n);
  for (let i = 0; i < n; i++) {
    let sAB = 0, sSq = 0, sSum = 0;
    for (let t = 0; t < k; t++) { const p = idx[t]; const a = P[i * J + pa[p]]; const bv = P[i * J + pb[p]]; const b = pr[p] < 0 ? (1 + minP[pb[p]]) - bv : bv;
      sAB += a * b; sSq += a * a + b * b; sSum += a + b; }
    const m = sSum / (2 * k); const num = sAB / k - m * m; const den = sSq / (2 * k) - m * m;
    out[i] = den > 1e-12 ? Math.abs(num / den) : 0;
  }
  return out;
}
function topIdx(pr, nPairs, k) { const o = Array.from({ length: nPairs }, (_, i) => i).sort((x, y) => Math.abs(pr[y]) - Math.abs(pr[x])); return Int32Array.from(o.slice(0, k)); }
function nullZ(cohFn, n, k, nPairs, iters, seed) { // returns per-resp mu,sd over random pair sets
  const rng = R._rng(seed); const sum = new Float64Array(n), sq = new Float64Array(n);
  for (let b = 0; b < iters; b++) { const pool = Int32Array.from({ length: nPairs }, (_, i) => i); const samp = new Int32Array(k);
    for (let t = 0; t < k; t++) { const r = t + Math.floor(rng() * (nPairs - t)); const tmp = pool[t]; pool[t] = pool[r]; pool[r] = tmp; samp[t] = pool[t]; }
    const c = cohFn(samp); for (let i = 0; i < n; i++) { sum[i] += c[i]; sq[i] += c[i] * c[i]; } }
  const mu = new Float64Array(n), sdn = new Float64Array(n);
  for (let i = 0; i < n; i++) { mu[i] = sum[i] / iters; sdn[i] = Math.sqrt(Math.max(0, sq[i] / iters - mu[i] * mu[i])); }
  return { mu, sdn };
}
const readM = p => { const P = R.parseCSV(fs.readFileSync(p, "utf8")), J = P.header.length; const rows = P.rows.map(r => r.map(Number)); return { rows, J }; };

const D = "../Dataset/gt_benchmark_candidates";
const SETS = [
  ["Study1", "_study_matrix.csv", "_study_labels.csv", "careless"],
  ["warning", D + "/warning_ipipneo300/_matrix.csv", "bogus_bench/warning_gt.csv", "gt"],
  ["Kay_S2", D + "/kay_idris_idria/kay_matrix_s2.csv", "bogus_bench/kay_s2_gt.csv", "gt"],
];
const ITERS = 800, FLIPS = 100, PERMS = 40, SHARE = 0.15;

for (const [name, mp, gp, ycol] of SETS) {
  const M = readM(mp), n = M.rows.length, J = M.J;
  const flat = new Float64Array(n * J); for (let i = 0; i < n; i++) for (let j = 0; j < J; j++) flat[i * J + j] = M.rows[i][j];
  const L = R.parseCSV(fs.readFileSync(gp, "utf8")); const yi = L.header.indexOf(ycol); const y = L.rows.slice(0, n).map(r => Number(r[yi]));
  const { P, minP } = prepare(flat, n, J); const { Z, pa, pb, pr, nPairs } = corPairs(P, n, J);
  const k = Math.min(Math.max(15, Math.round(0.03 * nPairs)), nPairs); const idx = topIdx(pr, nPairs, k);
  // baseline z for cor and SI (fixed null)
  const nulC = nullZ(s => cohCor(P, minP, n, J, pa, pb, pr, s, null), n, k, nPairs, ITERS, 7);
  const nulS = nullZ(s => cohSI(P, minP, n, J, pa, pb, pr, s), n, k, nPairs, ITERS, 7);
  const zc0raw = cohCor(P, minP, n, J, pa, pb, pr, idx, null), zs0raw = cohSI(P, minP, n, J, pa, pb, pr, idx);
  const zC = i => (nulC.sdn[i] > 1e-12 ? (zc0raw[i] - nulC.mu[i]) / nulC.sdn[i] : 0);
  const zS = i => (nulS.sdn[i] > 1e-12 ? (zs0raw[i] - nulS.mu[i]) / nulS.sdn[i] : 0);
  const zC0 = Array.from({ length: n }, (_, i) => zC(i)), zS0 = Array.from({ length: n }, (_, i) => zS(i));
  // low z = careless => flag the SHARE most-careless = lowest z
  const flagSet = z => { const o = z.map((v, i) => [v, i]).sort((a, b) => a[0] - b[0]).slice(0, Math.round(SHARE * n)).map(a => a[1]); return new Set(o); };
  const jacc = (A, B) => { let inter = 0; for (const x of A) if (B.has(x)) inter++; return inter / (A.size + B.size - inter); };
  const fC0 = flagSet(zC0), fS0 = flagSet(zS0);

  // Test A: random per-pair orientation flips (fixed pair set, fixed null)
  const rng = R._rng(123); let spC = 0, spS = 0, mvC = 0, mvS = 0, jaC = 0, jaS = 0;
  const sdzc = sd(zC0), sdzs = sd(zS0);
  for (let b = 0; b < FLIPS; b++) {
    const flip = Array.from({ length: k }, () => rng() < 0.5);
    const cc = cohCor(P, minP, n, J, pa, pb, pr, idx, flip);
    const zcb = Array.from({ length: n }, (_, i) => nulC.sdn[i] > 1e-12 ? (cc[i] - nulC.mu[i]) / nulC.sdn[i] : 0);
    // SI is invariant to flips by construction; recompute to confirm numerically
    const zsb = zS0.slice();
    spC += spearman(zC0, zcb); spS += spearman(zS0, zsb);
    let dC = 0, dS = 0; for (let i = 0; i < n; i++) { dC += Math.abs(zcb[i] - zC0[i]); dS += Math.abs(zsb[i] - zS0[i]); } mvC += dC / n / sdzc; mvS += dS / n / sdzs;
    jaC += jacc(fC0, flagSet(zcb)); jaS += jacc(fS0, flagSet(zsb));
  }
  // Test B: full column permutation (recompute everything, realistic "banal input permutation")
  const rng2 = R._rng(321); let spCp = 0, mvCp = 0, jaCp = 0, spSp = 0, aucCp = [], aucSp = [];
  for (let b = 0; b < PERMS; b++) {
    const perm = Int32Array.from({ length: J }, (_, i) => i); for (let i = J - 1; i > 0; i--) { const r = Math.floor(rng2() * (i + 1)); const t = perm[i]; perm[i] = perm[r]; perm[r] = t; }
    const Mp = new Float64Array(n * J); for (let i = 0; i < n; i++) for (let j = 0; j < J; j++) Mp[i * J + j] = flat[i * J + perm[j]];
    const pp = prepare(Mp, n, J); const cp = corPairs(pp.P, n, J); const kp = Math.min(Math.max(15, Math.round(0.03 * cp.nPairs)), cp.nPairs); const ip = topIdx(cp.pr, cp.nPairs, kp);
    const nc = nullZ(s => cohCor(pp.P, pp.minP, n, J, cp.pa, cp.pb, cp.pr, s, null), n, kp, cp.nPairs, 300, 7);
    const ns = nullZ(s => cohSI(pp.P, pp.minP, n, J, cp.pa, cp.pb, cp.pr, s), n, kp, cp.nPairs, 300, 7);
    const ccp = cohCor(pp.P, pp.minP, n, J, cp.pa, cp.pb, cp.pr, ip, null), csp = cohSI(pp.P, pp.minP, n, J, cp.pa, cp.pb, cp.pr, ip);
    const zcp = Array.from({ length: n }, (_, i) => nc.sdn[i] > 1e-12 ? (ccp[i] - nc.mu[i]) / nc.sdn[i] : 0);
    const zsp = Array.from({ length: n }, (_, i) => ns.sdn[i] > 1e-12 ? (csp[i] - ns.mu[i]) / ns.sdn[i] : 0);
    spCp += spearman(zC0, zcp); spSp += spearman(zS0, zsp);
    let dC = 0; for (let i = 0; i < n; i++) dC += Math.abs(zcp[i] - zC0[i]); mvCp += dC / n / sdzc;
    jaCp += jacc(fC0, flagSet(zcp));
    aucCp.push(auc(zcp.map(v => -v), y)); aucSp.push(auc(zsp.map(v => -v), y)); // -z: high = careless
  }
  console.log(`\n===== ${name} (n=${n}, J=${J}, k=${k}) =====`);
  console.log(`rank agreement cor~SI (Spearman): ${spearman(zC0, zS0).toFixed(3)}`);
  console.log(`AUC vs label:   cor ${auc(zC0.map(v => -v), y).toFixed(3)}   ICC ${auc(zS0.map(v => -v), y).toFixed(3)}`);
  console.log(`-- Test A: random per-pair orientation flips (${FLIPS} draws) --`);
  console.log(`  Spearman(z, z_flip):   cor ${(spC / FLIPS).toFixed(3)}   ICC ${(spS / FLIPS).toFixed(3)}`);
  console.log(`  mean|Δz| / sd(z):      cor ${(mvC / FLIPS).toFixed(3)}   ICC ${(mvS / FLIPS).toFixed(3)}`);
  console.log(`  flag-set Jaccard(${SHARE * 100}%): cor ${(jaC / FLIPS).toFixed(3)}   ICC ${(jaS / FLIPS).toFixed(3)}`);
  console.log(`-- Test B: full column permutation (${PERMS} perms) --`);
  console.log(`  Spearman(z, z_perm):   cor ${(spCp / PERMS).toFixed(3)}   ICC ${(spSp / PERMS).toFixed(3)}`);
  console.log(`  mean|Δz| / sd(z):      cor ${(mvCp / PERMS).toFixed(3)}`);
  console.log(`  flag-set Jaccard:      cor ${(jaCp / PERMS).toFixed(3)}`);
  console.log(`  AUC across perms:      cor ${(aucCp.reduce((s, v) => s + v, 0) / PERMS).toFixed(3)} ± ${sd(aucCp).toFixed(3)}   ICC ${(aucSp.reduce((s, v) => s + v, 0) / PERMS).toFixed(3)} ± ${sd(aucSp).toFixed(3)}`);
}
