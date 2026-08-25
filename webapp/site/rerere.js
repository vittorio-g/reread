/* rerere.js — ReReRe core library (browser + Node).
 *
 * Canonical algorithm (as in the paper): coupled-pair permutation z-score.
 *  1. Rescale items to proportions (x / item_max) so mixed Likert ranges
 *     are comparable (provably a no-op on homogeneous scales).
 *  2. Sample correlation matrix; take the top corProp% item pairs by |r|
 *     (at least minPairs).
 *  3. Per respondent: swap-invariant coherence across those coupled pairs
 *     (intraclass correlation = Pearson on the symmetrized (a,b)+(b,a) data),
 *     with reverse-coded partners sign-aligned via (max+min) - x. Symmetric
 *     because the coupled pairs have no natural left/right member.
 *  4. Permutation baseline: the same statistic on `iterations` random
 *     pair sets of the same size.
 *  5. z = (coupled - mean_random) / sd_random.  LOW z = careless.
 *
 * Everything is dependency-free and synchronous; run it in a Web Worker.
 */
(function (root, factory) {
  if (typeof module === "object" && module.exports) module.exports = factory();
  else root.ReReRe = factory();
})(typeof self !== "undefined" ? self : this, function () {
  "use strict";

  /* ---------------- seeded RNG (mulberry32) — reproducible runs -------- */
  function rng(seed) {
    let a = seed >>> 0;
    return function () {
      a |= 0; a = (a + 0x6D2B79F5) | 0;
      let t = Math.imul(a ^ (a >>> 15), 1 | a);
      t = (t + Math.imul(t ^ (t >>> 7), 61 | t)) ^ t;
      return ((t ^ (t >>> 14)) >>> 0) / 4294967296;
    };
  }

  /* ---------------- CSV parsing (delimiter-sniffing, quoted fields) ---- */
  function sniffDelimiter(text) {
    const head = text.slice(0, 4000);
    const cand = [",", ";", "\t"];
    let best = ",", bestScore = -1;
    for (const d of cand) {
      const lines = head.split(/\r?\n/).filter(l => l.length).slice(0, 10);
      if (!lines.length) continue;
      const counts = lines.map(l => l.split(d).length);
      const median = counts.sort((a, b) => a - b)[Math.floor(counts.length / 2)];
      const consistent = counts.filter(c => c === median).length / counts.length;
      const score = (median > 1 ? median : 0) * consistent;
      if (score > bestScore) { bestScore = score; best = d; }
    }
    return best;
  }

  function parseCSV(text, delimiter) {
    text = text.replace(/^﻿/, "");                       // strip BOM
    const d = delimiter || sniffDelimiter(text);
    const rows = [];
    let row = [], field = "", inQ = false;
    for (let i = 0; i < text.length; i++) {
      const c = text[i];
      if (inQ) {
        if (c === '"') {
          if (text[i + 1] === '"') { field += '"'; i++; }
          else inQ = false;
        } else field += c;
      } else if (c === '"') inQ = true;
      else if (c === d) { row.push(field); field = ""; }
      else if (c === "\n") { row.push(field); rows.push(row); row = []; field = ""; }
      else if (c === "\r") { /* skip */ }
      else field += c;
    }
    if (field.length || row.length) { row.push(field); rows.push(row); }
    while (rows.length && rows[rows.length - 1].every(f => f === "")) rows.pop();
    if (!rows.length) throw new Error("Empty file");
    return { header: rows[0].map(h => String(h).trim()), rows: rows.slice(1), delimiter: d };
  }

  /* ---------------- data preparation ----------------------------------- */
  /* Numeric-ish parser: accepts "3", "3.0", "3,5" (EU decimal), "" -> NaN */
  function toNum(s) {
    if (s == null) return NaN;
    s = String(s).trim();
    if (s === "" || /^(na|n\/a|nan|null|missing|\.)$/i.test(s)) return NaN;
    const v = Number(s.replace(",", "."));
    return Number.isFinite(v) ? v : NaN;
  }

  /* Decide which columns are items: >=80% numeric among non-empty cells,
   * and at least 2 distinct values. Returns {itemCols, idCol, excluded}. */
  function detectColumns(header, rows) {
    const n = rows.length, J = header.length;
    const itemCols = [], excluded = [];
    let idCol = -1;
    for (let j = 0; j < J; j++) {
      let nonEmpty = 0, numeric = 0;
      const distinct = new Set();
      for (let i = 0; i < n; i++) {
        const raw = rows[i][j];
        if (raw == null || String(raw).trim() === "") continue;
        nonEmpty++;
        const v = toNum(raw);
        if (Number.isFinite(v)) { numeric++; distinct.add(v); }
      }
      const name = header[j] || ("col" + (j + 1));
      const looksId = /^(id|participant|subject|resp|pid|prolific|code|matricola|token)/i.test(name);
      if (nonEmpty === 0) { excluded.push({ col: name, reason: "empty" }); continue; }
      if (numeric / nonEmpty >= 0.8 && distinct.size >= 2 && !looksId) itemCols.push(j);
      else {
        if (idCol < 0 && (looksId || numeric / nonEmpty < 0.8)) idCol = j;
        excluded.push({ col: name, reason: looksId ? "looks like an ID" : "non-numeric or constant" });
      }
    }
    return { itemCols, idCol, excluded };
  }

  /* Build numeric matrix (n x J) with NaN for missing. */
  function buildMatrix(rows, itemCols) {
    const n = rows.length, J = itemCols.length;
    const M = new Float64Array(n * J);
    for (let i = 0; i < n; i++)
      for (let j = 0; j < J; j++)
        M[i * J + j] = toNum(rows[i][itemCols[j]]);
    return { M, n, J };
  }

  /* ---------------- core math ------------------------------------------ */
  function colStats(M, n, J) {
    const max = new Float64Array(J).fill(-Infinity);
    const min = new Float64Array(J).fill(Infinity);
    const median = new Float64Array(J);
    const buf = [];
    for (let j = 0; j < J; j++) {
      buf.length = 0;
      for (let i = 0; i < n; i++) {
        const v = M[i * J + j];
        if (Number.isFinite(v)) {
          buf.push(v);
          if (v > max[j]) max[j] = v;
          if (v < min[j]) min[j] = v;
        }
      }
      buf.sort((a, b) => a - b);
      median[j] = buf.length ? buf[Math.floor(buf.length / 2)] : 0;
      if (!Number.isFinite(max[j])) { max[j] = 1; min[j] = 0; }
      if (max[j] === 0) max[j] = 1;
    }
    return { max, min, median };
  }

  /* Proportion-rescaled, median-imputed working matrix P; missing count. */
  function prepare(M, n, J) {
    const { max, min, median } = colStats(M, n, J);
    const P = new Float64Array(n * J);
    const nMissing = new Int32Array(n);
    for (let i = 0; i < n; i++)
      for (let j = 0; j < J; j++) {
        let v = M[i * J + j];
        if (!Number.isFinite(v)) { nMissing[i]++; v = median[j]; }
        P[i * J + j] = v / max[j];
      }
    const minP = new Float64Array(J);
    for (let j = 0; j < J; j++) minP[j] = min[j] / max[j];
    return { P, minP, nMissing };
  }

  /* Column-standardised copy of P (for the sample correlation matrix). */
  function correlationPairs(P, n, J) {
    const mean = new Float64Array(J), sd = new Float64Array(J);
    for (let j = 0; j < J; j++) {
      let s = 0; for (let i = 0; i < n; i++) s += P[i * J + j];
      mean[j] = s / n;
      let q = 0; for (let i = 0; i < n; i++) { const d = P[i * J + j] - mean[j]; q += d * d; }
      sd[j] = Math.sqrt(q / (n - 1)) || 0;
    }
    const Z = new Float64Array(n * J);
    for (let i = 0; i < n; i++)
      for (let j = 0; j < J; j++)
        Z[i * J + j] = sd[j] > 0 ? (P[i * J + j] - mean[j]) / sd[j] : 0;
    const nPairs = (J * (J - 1)) / 2;
    const pa = new Int32Array(nPairs), pb = new Int32Array(nPairs), pr = new Float64Array(nPairs);
    let p = 0;
    for (let a = 0; a < J; a++)
      for (let b = a + 1; b < J; b++, p++) {
        let s = 0;
        for (let i = 0; i < n; i++) s += Z[i * J + a] * Z[i * J + b];
        pa[p] = a; pb[p] = b; pr[p] = (sd[a] > 0 && sd[b] > 0) ? s / (n - 1) : 0;
      }
    return { pa, pb, pr, nPairs };
  }

  /* Per-respondent coherence across a pair set (sign-aligned partners), as a
   * SWAP-INVARIANT intraclass-style correlation. The coupled pairs are exchangeable
   * (there is no natural "left"/"right" member), so a plain cor(A,B) between the
   * left-vector and right-vector is ill-defined: it depends on the arbitrary column
   * order that fixes which member is A vs B, and swapping one pair changes it. We use
   * the intraclass correlation instead = Pearson on the symmetrized data (a,b)+(b,a):
   * with m the mean of all 2k members, num = mean(a*b) - m^2, den = mean((a^2+b^2)/2) - m^2.
   * This is exactly invariant to per-pair member swaps and to input column order, while
   * ranking respondents essentially identically to the old |cor| (Spearman >= .89). */
  function indCors(P, minP, n, J, pa, pb, pr, idx) {
    const k = idx.length;
    const out = new Float64Array(n);
    for (let i = 0; i < n; i++) {
      let sAB = 0, sSq = 0, sSum = 0;
      for (let t = 0; t < k; t++) {
        const p = idx[t];
        const a = P[i * J + pa[p]];
        const bj = pb[p];
        const bv = P[i * J + bj];
        const b = pr[p] < 0 ? (1 + minP[bj]) - bv : bv;   // reverse-code partner
        sAB += a * b; sSq += a * a + b * b; sSum += a + b;
      }
      const m = sSum / (2 * k);
      const num = sAB / k - m * m;
      const den = sSq / (2 * k) - m * m;
      out[i] = den > 1e-12 ? Math.abs(num / den) : 0;
    }
    return out;
  }

  function sampleIdx(nPairs, k, rand) {
    // partial Fisher-Yates over an index pool
    const pool = new Int32Array(nPairs);
    for (let i = 0; i < nPairs; i++) pool[i] = i;
    const out = new Int32Array(k);
    for (let t = 0; t < k; t++) {
      const r = t + Math.floor(rand() * (nPairs - t));
      const tmp = pool[t]; pool[t] = pool[r]; pool[r] = tmp;
      out[t] = pool[t];
    }
    return out;
  }

  /* ---------------- structure diagnostic ------------------------------- */
  /* Two ingredients:
   *  - separation ratio: mean(|r| of top-k pairs) / mean(|r| of all pairs).
   *    NOTE: under pure noise this ratio self-inflates (~3) because the
   *    top tail of |r| is a sampling artefact — so it cannot stand alone.
   *  - signal strength: mean(|r| of top-k pairs) relative to the H0 noise
   *    expectation. Under independence r ~ N(0, 1/(n-1)); the mean of the
   *    top ~3% of |r| is ≈ 2.4/sqrt(n). strength = top_mean_r / (2.4/sqrt(n)).
   *    Noise-only data gives strength ≈ 1 regardless of J.
   * The level is driven by signal strength (noise-calibrated). */
  function diagnostic(pr, k, n) {
    const abs = Array.from(pr, Math.abs).sort((a, b) => b - a);
    const allMean = abs.reduce((s, v) => s + v, 0) / abs.length;
    const topMean = abs.slice(0, k).reduce((s, v) => s + v, 0) / k;
    const ratio = allMean > 0 ? topMean / allMean : NaN;
    const noiseTop = 2.4 / Math.sqrt(Math.max(n, 2));
    const strength = topMean / noiseTop;
    let level, advice;
    if (!Number.isFinite(strength) || strength < 1.3) {
      level = "no_structure";
      advice = "The strongest item correlations are indistinguishable from sampling noise: this dataset shows no usable multi-construct structure. ReReRe is NOT suitable here — rely on LongString / IRV / response times instead.";
    } else if (strength < 1.8) {
      level = "weak";
      advice = "Weak structure: the coupled pairs barely rise above sampling noise. ReReRe may produce false positives — treat z-scores as a ranking, not a hard flag, and corroborate with other indices.";
    } else if (strength < 2.8) {
      level = "marginal";
      advice = "Marginal structure: ReReRe is usable, but detection quality will be reduced. Consider a stricter threshold and corroborating indices.";
    } else {
      level = "ok";
      advice = "Good multi-construct structure — ReReRe should work well on this dataset.";
    }
    return { separation_ratio: ratio, signal_strength: strength, level, advice,
             top_mean_r: topMean, all_mean_r: allMean, noise_top_r: noiseTop };
  }

  /* ---------------- auxiliary indices ---------------------------------- */
  function auxiliaries(P, n, J) {
    // IRV: within-person SD on the proportion scale
    const irv = new Float64Array(n);
    // person-total: r between the person's profile and the sample mean profile
    const cm = new Float64Array(J);
    for (let j = 0; j < J; j++) {
      let s = 0; for (let i = 0; i < n; i++) s += P[i * J + j];
      cm[j] = s / n;
    }
    const personTotal = new Float64Array(n);
    // longstring: longest run of identical raw answers in column order
    const longstring = new Int32Array(n);
    for (let i = 0; i < n; i++) {
      let m = 0; for (let j = 0; j < J; j++) m += P[i * J + j];
      m /= J;
      let q = 0; for (let j = 0; j < J; j++) { const d = P[i * J + j] - m; q += d * d; }
      irv[i] = Math.sqrt(q / (J - 1));
      let mc = 0; for (let j = 0; j < J; j++) mc += cm[j]; mc /= J;
      let sab = 0, sa = 0, sb = 0;
      for (let j = 0; j < J; j++) {
        const da = P[i * J + j] - m, db = cm[j] - mc;
        sab += da * db; sa += da * da; sb += db * db;
      }
      personTotal[i] = (sa > 0 && sb > 0) ? sab / Math.sqrt(sa * sb) : 0;
      let best = 1, run = 1;
      for (let j = 1; j < J; j++) {
        if (P[i * J + j] === P[i * J + j - 1]) { run++; if (run > best) best = run; }
        else run = 1;
      }
      longstring[i] = best;
    }
    return { irv, personTotal, longstring };
  }

  /* ---------------- Mahalanobis D^2 (pseudo-inverse, N<p safe) ---------- */
  /* Symmetric eigendecomposition via cyclic Jacobi (m up to a few hundred).
   * Only squared eigenvector entries / squared projections are used downstream,
   * so eigenvector sign conventions are irrelevant. */
  function jacobiEigen(Ain, m) {
    const A = Float64Array.from(Ain), V = new Float64Array(m * m);
    for (let i = 0; i < m; i++) V[i * m + i] = 1;
    for (let sweep = 0; sweep < 100; sweep++) {
      let off = 0;
      for (let p = 0; p < m; p++) for (let q = p + 1; q < m; q++) { const a = A[p * m + q]; off += a * a; }
      if (off < 1e-18) break;
      for (let p = 0; p < m; p++) for (let q = p + 1; q < m; q++) {
        const apq = A[p * m + q];
        if (Math.abs(apq) < 1e-300) continue;
        const app = A[p * m + p], aqq = A[q * m + q];
        const phi = 0.5 * Math.atan2(2 * apq, app - aqq);
        const c = Math.cos(phi), s = Math.sin(phi);
        for (let i = 0; i < m; i++) { const aip = A[i * m + p], aiq = A[i * m + q]; A[i * m + p] = c * aip - s * aiq; A[i * m + q] = s * aip + c * aiq; }
        for (let i = 0; i < m; i++) { const api = A[p * m + i], aqi = A[q * m + i]; A[p * m + i] = c * api - s * aqi; A[q * m + i] = s * api + c * aqi; }
        for (let i = 0; i < m; i++) { const vip = V[i * m + p], viq = V[i * m + q]; V[i * m + p] = c * vip - s * viq; V[i * m + q] = s * vip + c * viq; }
      }
    }
    const val = new Float64Array(m);
    for (let i = 0; i < m; i++) val[i] = A[i * m + i];
    return { val, vec: V };
  }

  /* Mahalanobis D^2 in the leading principal-component subspace.
   *
   * Full pseudo-inverse Mahalanobis is degenerate when n<=p: the data fill the
   * (n-1)-dim centred subspace, so every respondent is equidistant and D^2 is
   * constant. Instead we keep the components explaining `keepVar` of the total
   * variance (default 90%) and whiten within that subspace. This yields a
   * genuine outlyingness score in BOTH regimes: with n>p and keepVar->1 it
   * reproduces the exact Mahalanobis distance; with n<=p it drops the noisy
   * minor directions that make the full pseudo-inverse constant.
   *
   * Eigendecompose whichever of the Gram (n×n) or covariance (p×p) matrix is
   * smaller. Per-component contribution to D^2_i:
   *   covariance route: (Xc_i · v_c)^2 / lambda_c
   *   Gram route:       (n-1) · U_ic^2         (identical, via Xc V = U Sigma) */
  function mahalanobisD2(P, n, J, keepVar) {
    keepVar = keepVar != null ? keepVar : 0.90;
    const mean = new Float64Array(J), keepCol = [];
    for (let j = 0; j < J; j++) {
      let s = 0; for (let i = 0; i < n; i++) s += P[i * J + j]; mean[j] = s / n;
      let q = 0; for (let i = 0; i < n; i++) { const d = P[i * J + j] - mean[j]; q += d * d; }
      if (q > 1e-12) keepCol.push(j);
    }
    const p = keepCol.length;
    const Xc = new Float64Array(n * p);
    for (let i = 0; i < n; i++) for (let c = 0; c < p; c++) Xc[i * p + c] = P[i * J + keepCol[c]] - mean[keepCol[c]];
    const d2 = new Float64Array(n), REL = 1e-9;

    // lambda_c = per-component variance; term(i,c) = contribution to D^2_i.
    let m, lambda, termAt;
    if (n <= p) {
      const G = new Float64Array(n * n);
      for (let i = 0; i < n; i++) for (let j = i; j < n; j++) {
        let s = 0; for (let c = 0; c < p; c++) s += Xc[i * p + c] * Xc[j * p + c];
        G[i * n + j] = s; G[j * n + i] = s;
      }
      const { val, vec } = jacobiEigen(G, n); m = n;
      lambda = Float64Array.from(val, g => g / (n - 1));      // variance
      termAt = (i, c) => { const u = vec[i * n + c]; return (n - 1) * u * u; };
    } else {
      const S = new Float64Array(p * p);
      for (let a = 0; a < p; a++) for (let b = a; b < p; b++) {
        let s = 0; for (let i = 0; i < n; i++) s += Xc[i * p + a] * Xc[i * p + b];
        s /= (n - 1); S[a * p + b] = s; S[b * p + a] = s;
      }
      const { val, vec } = jacobiEigen(S, p); m = p;
      lambda = val;
      termAt = (i, c) => {
        const lam = lambda[c]; if (!(lam > 0)) return 0;
        let proj = 0; for (let a = 0; a < p; a++) proj += Xc[i * p + a] * vec[a * p + c];
        return proj * proj / lam;
      };
    }
    // rank components by variance, keep the leading ones up to keepVar
    let lmax = 0, total = 0;
    for (let c = 0; c < m; c++) { if (lambda[c] > lmax) lmax = lambda[c]; if (lambda[c] > 0) total += lambda[c]; }
    const tol = lmax * REL;
    const idx = Array.from({ length: m }, (_, c) => c)
      .filter(c => lambda[c] > tol).sort((a, b) => lambda[b] - lambda[a]);
    const kept = []; let cum = 0;
    for (const c of idx) { kept.push(c); cum += lambda[c]; if (cum >= keepVar * total) break; }
    for (let i = 0; i < n; i++) { let acc = 0; for (const c of kept) acc += termAt(i, c); d2[i] = acc; }
    return d2;
  }

  /* ---------------- ensemble (ReReRe) ---------------------------------- */
  /* Robust within-dataset z (median / MAD, SD fallback) — makes each feature
   * dataset-invariant so a single fixed weight vector transfers across
   * questionnaires of different length and scale. */
  function robustZ(x, n) {
    const s = Float64Array.from(x).sort((a, b) => a - b);
    const med = s[Math.floor(n / 2)];
    const dev = Float64Array.from(x, v => Math.abs(v - med)).sort((a, b) => a - b);
    let mad = dev[Math.floor(n / 2)] * 1.4826;
    if (!(mad > 1e-9)) {
      let m = 0; for (let i = 0; i < n; i++) m += x[i]; m /= n;
      let q = 0; for (let i = 0; i < n; i++) { const d = x[i] - m; q += d * d; }
      const sd = Math.sqrt(q / Math.max(n - 1, 1)) || 1;
      return Float64Array.from(x, v => (v - m) / sd);
    }
    return Float64Array.from(x, v => (v - med) / mad);
  }

  /* Compute the five ensemble features and orient them so HIGHER = more
   * careless, each as a robust within-dataset z. Reuses the rr engine. */
  function ensembleFeatures(M, n, J, opts) {
    const base = run(M, n, J, opts);
    const { P } = prepare(M, n, J);
    // d2 carries zero ensemble weight; skip its (costly, N>p) eigendecomposition
    // when an analysis does not need the diagnostic column.
    const d2 = (opts && opts.skipD2) ? new Float64Array(n) : mahalanobisD2(P, n, J);
    // Profile informativeness = Var(item means) / mean(item variance). When the
    // item-mean profile is flat (near 0), the partner detectors (Person-Total,
    // LongString) are unreliable — Person-Total correlates against a flat mean
    // profile (noise) and attentive respondents can show incidental runs. This
    // gates the partners so the ensemble falls back to rr on such data.
    const cm = new Float64Array(J), cv = new Float64Array(J);
    for (let j = 0; j < J; j++) {
      let s = 0; for (let i = 0; i < n; i++) s += P[i * J + j]; cm[j] = s / n;
      let q = 0; for (let i = 0; i < n; i++) { const d = P[i * J + j] - cm[j]; q += d * d; } cv[j] = q / (n - 1);
    }
    let gm = 0; for (let j = 0; j < J; j++) gm += cm[j]; gm /= J;
    let vm = 0; for (let j = 0; j < J; j++) vm += (cm[j] - gm) * (cm[j] - gm); vm /= J;
    let mv = 0; for (let j = 0; j < J; j++) mv += cv[j]; mv /= J;
    const profileInfo = mv > 0 ? vm / mv : 0;
    const o_rr = robustZ(base.z, n);           // low rr  = careless -> negate
    const o_irv = robustZ(base.irv, n);        // high    = careless
    const o_ls = robustZ(base.longstring, n);  // high    = careless
    const o_d2 = robustZ(d2, n);               // high    = careless
    const o_pt = robustZ(base.personTotal, n); // low     = careless -> negate
    const O = { rr: [], irv: [], longstring: [], d2: [], person_total: [] };
    for (let i = 0; i < n; i++) {
      O.rr.push(-o_rr[i]); O.irv.push(o_irv[i]); O.longstring.push(o_ls[i]);
      O.d2.push(o_d2[i]); O.person_total.push(-o_pt[i]);
    }
    return { base, oriented: O, d2: Array.from(d2), profileInfo };
  }

  /* Fixed logistic weights on the oriented robust-z features, fitted and
   * leave-one-out validated on the real careless-study data (n=84, 109 items).
   * Shipped model = the transferable triad rr + longstring + person_total
   * (LOO AUC 0.985, MCC 0.880 — matching the full 5-feature RF, AUC 0.965).
   * irv and d2 are still computed and shown as diagnostic columns but carry
   * zero weight: irv's careless direction depends on the careless type, and the
   * pseudo-inverse d2 degenerates in the n<=p regime common to careless data.
   *
   * rr (low = careless) is the robust anchor. The two partners (longstring high
   * = careless, person_total low = careless) are helpful when the item-mean
   * profile is informative but mislead on flat/orthogonal data, so their joint
   * contribution is multiplied by a per-dataset reliability gate derived from
   * profileInfo (see ensembleFeatures). On flat data the gate -> 0 and the
   * ensemble degrades gracefully to rr alone (never worse); on real batteries
   * the gate -> 1 and the ensemble reaches LOO AUC 0.984 on the study data. */
  var WEIGHTS = { b0: -0.1057, rr: 2.7645, irv: 0, longstring: 1.2689, d2: 0, person_total: 1.4452 };

  /* ---------------- automatic calibration -------------------------------- */
  /* 1-D two-component Gaussian mixture via EM on the ensemble log-odds. Used
   * ONLY as a structure gate: to decide, data-drivenly, whether a careless
   * subpopulation exists at all (BIC of the 1- vs 2-component fits). It is NOT
   * used to place the cut or to estimate the rate — both come from the
   * one-sided fit of the attentive mode (leftFit below), which the careless
   * cannot contaminate. A full two-component fit collapses at high prevalence:
   * its low component absorbs careless observations and the estimate shrinks. */
  function gaussMix1D(x) {
    const n = x.length;
    let mean = 0; for (let i = 0; i < n; i++) mean += x[i]; mean /= n;
    let varr = 0; for (let i = 0; i < n; i++) varr += (x[i] - mean) * (x[i] - mean); varr /= n;
    const sd = Math.sqrt(varr) || 1, vfloor = Math.max(1e-4, varr * 0.02);
    let m1 = mean - 0.7 * sd, m2 = mean + 0.7 * sd, v1 = varr, v2 = varr, pi = 0.3;
    const norm = (v, m, s2) => Math.exp(-((v - m) * (v - m)) / (2 * s2)) / Math.sqrt(2 * Math.PI * s2);
    const post = new Array(n).fill(0);
    for (let it = 0; it < 300; it++) {
      let sp = 0;
      for (let i = 0; i < n; i++) { const a = pi * norm(x[i], m2, v2), b = (1 - pi) * norm(x[i], m1, v1); post[i] = a / (a + b + 1e-300); sp += post[i]; }
      pi = Math.min(0.999, Math.max(1e-3, sp / n));
      let s2 = 0, s1 = 0; for (let i = 0; i < n; i++) { s2 += post[i] * x[i]; s1 += (1 - post[i]) * x[i]; }
      m2 = s2 / (sp + 1e-9); m1 = s1 / (n - sp + 1e-9);
      let q2 = 0, q1 = 0; for (let i = 0; i < n; i++) { q2 += post[i] * (x[i] - m2) * (x[i] - m2); q1 += (1 - post[i]) * (x[i] - m1) * (x[i] - m1); }
      v2 = Math.max(vfloor, q2 / (sp + 1e-9)); v1 = Math.max(vfloor, q1 / (n - sp + 1e-9));
    }
    if (m2 < m1) { const tm = m1; m1 = m2; m2 = tm; const tv = v1; v1 = v2; v2 = tv; pi = 1 - pi; for (let i = 0; i < n; i++) post[i] = 1 - post[i]; }
    let ll2 = 0; for (let i = 0; i < n; i++) ll2 += Math.log(pi * norm(x[i], m2, v2) + (1 - pi) * norm(x[i], m1, v1) + 1e-300);
    let ll1 = 0; for (let i = 0; i < n; i++) ll1 += Math.log(norm(x[i], mean, varr) + 1e-300);
    return { pi, post, m1, m2, v1, v2, bic1: -2 * ll1 + 2 * Math.log(n), bic2: -2 * ll2 + 5 * Math.log(n) };
  }

  /* One-sided fit of the ATTENTIVE mode: the DOMINANT KDE mode m of the score
   * distribution, and a scale estimated from the LEFT side only (MAD of the
   * points below m). The left side is careless-free by construction at any
   * realistic prevalence, so neither quantity is contaminated when the careless
   * are numerous. The careless, by contrast, form a heavy multi-type TAIL with
   * no mode of their own (verified empirically up to a 45% rate), so the
   * attentive mode is the only anchor the score distribution offers.
   *
   * m is the GLOBAL maximum of the density, not the first local maximum clearing
   * 15% of the peak as in earlier versions: that rule let a spurious bump in the
   * low tail win over the dominant mode, after which sigma was estimated from the
   * handful of points below it and collapsed, taking the cut with it. The global
   * mode is free — identical flags on Study 1 and on six labelled external
   * datasets — and strictly better on clean simulated data with few respondents,
   * where the old rule occasionally flagged a fifth of a clean sample. */
  function leftFit(eta) {
    const n = eta.length, s = Array.from(eta).sort(function (a, b) { return a - b; });
    const q1 = s[Math.floor(0.25 * n)], q3 = s[Math.floor(0.75 * n)], iqr = q3 - q1;
    let mu = 0; for (let i = 0; i < n; i++) mu += eta[i]; mu /= n;
    let va = 0; for (let i = 0; i < n; i++) va += (eta[i] - mu) * (eta[i] - mu); va /= n;
    const sd = Math.sqrt(va);
    let h = 0.9 * Math.min(sd || 1, (iqr / 1.34) || (sd || 1)) * Math.pow(n, -0.2);
    if (!(h > 0)) h = (sd || 1) * 0.3;
    const G = 512, lo = s[0] - h, hi = s[n - 1] + h, dens = new Float64Array(G);
    for (let g = 0; g < G; g++) {
      const xg = lo + (hi - lo) * g / (G - 1); let d = 0;
      for (let i = 0; i < n; i++) { const zz = (xg - eta[i]) / h; d += Math.exp(-0.5 * zz * zz); }
      dens[g] = d;
    }
    let dmax = -1, gm = -1;
    for (let g = 0; g < G; g++) if (dens[g] > dmax) { dmax = dens[g]; gm = g; }
    const m = gm >= 0 ? lo + (hi - lo) * gm / (G - 1) : s[Math.floor(n / 2)];
    /* Scale of the attentive mode from its LEFT half-width at half maximum.
     * The MAD of the whole left side is inflated when eta has a long coherent
     * tail (|z_rc| grows with battery length), which pushed the 2.5-sigma cut
     * beyond everyone and drove the one-sided prevalence negative. The HWHM
     * measures the shoulder of the peak itself: the far tail lies below
     * half-maximum and cannot influence it, and the careless mass sits on the
     * right side, which is never used. The KDE kernel widens the apparent
     * shoulder, so its bandwidth is deconvolved (floored at half the raw
     * width to keep the estimate from collapsing when hw ~ h). */
    let sigHW = 0;
    if (gm > 0) {
      let gl = gm; while (gl > 0 && dens[gl] > 0.5 * dmax) gl--;
      if (dens[gl] <= 0.5 * dmax) {
        const xHalf = lo + (hi - lo) * gl / (G - 1);
        const hw = (m - xHalf) / 1.17741;
        sigHW = Math.sqrt(Math.max(hw * hw - h * h, 0.25 * hw * hw));
      }
    }
    const below = [];
    for (let i = 0; i < n; i++) if (eta[i] < m) below.push(m - eta[i]);
    below.sort(function (a, b) { return a - b; });
    let sigMad = below.length ? 1.4826 * below[Math.floor(below.length / 2)] : 0;
    if (!(sigMad > 0)) sigMad = (sd || 1e-9);
    const sig = (sigHW > 0) ? sigHW : sigMad;   // MAD only as fallback (degenerate KDE)
    return { m: m, sigma: sig, sigmaMad: sigMad, nBelow: below.length };
  }

  /* Decide flags automatically from the log-odds. First a BIC mixture test asks
   * whether a careless subpopulation exists at all (so clean data flags ~nobody,
   * no false-positive cascade) — the mixture serves ONLY as this structure gate.
   * When the gate opens, both the flagging cut and the prevalence estimate come
   * from the one-sided fit of the attentive mode (leftFit): the cut is the
   * z-sigma upper tail of that mode, and the prevalence is estimated by counting
   * the mass below m + 1*sigma (careless-free by construction) and dividing by
   * Phi(1) — so neither quantity is distorted when the careless are numerous,
   * which is where a full two-component fit collapses (its low component absorbs
   * careless and the estimate shrinks).
   *
   * A SINGLE cut is shipped: z=2.5 (~0.6% FPR-on-attentive). Anchoring on the
   * one-sided fit extends its useful range to a ~35% careless rate — past any
   * realistic prevalence — so the second ("High", z=1.5) preset that the earlier
   * two-setting scheme offered no longer earns its place and was removed. Legacy
   * names (high/low/medium) still resolve, to z=2.5, so old calls keep working.
   * Beyond ~35% careless, pass a raw numeric z or a manual prevalence instead. */
  var SENS_Z = { standard: 2.5, high: 2.5, low: 2.5, medium: 2.5 };
  var PHI1 = 0.8413447;
  /* upper tail of the standard normal, and the upper-tail binomial p-value (log space).
   * Used by the excess-tail branch of the structure gate: under "no careless subpopulation"
   * the share of respondents above m + z*sigma is 1 - Phi(z), so a tail much heavier than
   * that is itself evidence that a careless group exists - a test the EM mixture is too
   * weak to pass on small samples even when the group is obvious. */
  function normUpper(z) {
    const t = 1 / (1 + 0.2316419 * Math.abs(z));
    const d = 0.3989422804014327 * Math.exp(-z * z / 2);
    const p = d * t * (0.319381530 + t * (-0.356563782 + t * (1.781477937 +
              t * (-1.821255978 + t * 1.330274429))));
    return z >= 0 ? p : 1 - p;
  }
  function lgamma(x) {
    const c = [676.5203681218851, -1259.1392167224028, 771.32342877765313,
               -176.61502916214059, 12.507343278686905, -0.13857109526572012,
               9.9843695780195716e-6, 1.5056327351493116e-7];
    if (x < 0.5) return Math.log(Math.PI / Math.sin(Math.PI * x)) - lgamma(1 - x);
    x -= 1; let a = 0.99999999999980993; const t = x + 7.5;
    for (let i = 0; i < 8; i++) a += c[i] / (x + i + 1);
    return 0.5 * Math.log(2 * Math.PI) + (x + 0.5) * Math.log(t) - t + Math.log(a);
  }
  function binomUpperP(k, n, p) {
    if (k <= 0) return 1;
    if (k > n) return 0;
    const lch = function (n, i) { return lgamma(n + 1) - lgamma(i + 1) - lgamma(n - i + 1); };
    const terms = []; let mx = -Infinity;
    for (let i = k; i <= n; i++) {
      const lt = lch(n, i) + i * Math.log(p) + (n - i) * Math.log1p(-p);
      terms.push(lt); if (lt > mx) mx = lt;
      if (i > k + 5 && lt < mx - 50) break;
    }
    let s = 0; for (let j = 0; j < terms.length; j++) s += Math.exp(terms[j] - mx);
    return Math.min(1, Math.exp(mx) * s);
  }
  /* The tail branch must be overwhelming, not merely significant. An n-adaptive relaxation
   * (1e-3 below n=400) was tried and REJECTED by validation on 2026-08-06: it did not recover
   * availability at low prevalence (power, not threshold, is binding at small n) while raising
   * the clean-data fire rate at n=350 from 4% to 16%. Keep the flat conservative value. */
  function tailAlpha(n) { return 1e-6; }

  /* ---------------- clean-data false-flag envelope ---------------------- */
  /* With the gate removed the cut is unconditional, so a clean sample flags
   * whatever its own upper tail puts past m + z*sigma rather than exactly nobody.
   * That share is almost entirely a function of the two things the analyst
   * already knows -- battery length and sample size -- so it can be measured once
   * and reported: CLEAN_ENV holds the mean and 95th percentile of the flagged
   * share on genuinely clean simulated data over a grid of both (30 replications
   * per cell, structural parameters varied across replications; see
   * webapp/clean_envelope.js). A flagged share inside this envelope is
   * indistinguishable from what clean data produces by chance. The old structure
   * gate asserted the same conclusion by refusing to act, which also discarded
   * the samples whose share was genuinely above it. */
  var CLEAN_ENV = { js: [30, 48, 60, 90, 120, 180, 240, 300], ns: [100, 150, 200, 300, 500, 1000],
    cells: [
      [0.0247,0.2000], [0.0047,0.0200], [0.0025,0.0100], [0.0073,0.0500], [0.0101,0.0880], [0.0040,0.0250], 
      [0.0293,0.1300], [0.0184,0.0933], [0.0035,0.0300], [0.0128,0.0400], [0.0039,0.0200], [0.0033,0.0130], 
      [0.0370,0.1600], [0.0196,0.1533], [0.0143,0.0700], [0.0164,0.1167], [0.0139,0.0640], [0.0063,0.0240], 
      [0.0620,0.2400], [0.0249,0.1400], [0.0317,0.1200], [0.0206,0.1333], [0.0127,0.0660], [0.0065,0.0480], 
      [0.0380,0.1500], [0.0256,0.1067], [0.0118,0.0500], [0.0148,0.0833], [0.0159,0.0680], [0.0053,0.0230], 
      [0.0543,0.1900], [0.0311,0.1533], [0.0255,0.1200], [0.0074,0.0333], [0.0135,0.0920], [0.0103,0.0650], 
      [0.0303,0.1700], [0.0304,0.1200], [0.0287,0.0750], [0.0138,0.0633], [0.0167,0.0680], [0.0124,0.0400], 
      [0.0393,0.1300], [0.0358,0.1133], [0.0215,0.0950], [0.0189,0.0900], [0.0143,0.0520], [0.0122,0.0370]
    ] };
  function cleanEnvelope(n, J) {
    const gj = CLEAN_ENV.js, gn = CLEAN_ENV.ns;
    /* bilinear interpolation in log-log, clamped at the grid edges: outside the
     * measured range the nearest cell is the honest answer, not an extrapolation */
    const pos = (v, arr) => {
      if (v <= arr[0]) return [0, 0, 0];
      if (v >= arr[arr.length - 1]) return [arr.length - 1, arr.length - 1, 0];
      let k = 0; while (k < arr.length - 2 && arr[k + 1] < v) k++;
      const t = (Math.log(v) - Math.log(arr[k])) / (Math.log(arr[k + 1]) - Math.log(arr[k]));
      return [k, k + 1, t];
    };
    const [j0, j1, tj] = pos(J, gj), [n0, n1, tn] = pos(n, gn);
    const at = (a, b) => CLEAN_ENV.cells[a * gn.length + b];
    const mix = key => {
      const v00 = at(j0, n0)[key], v01 = at(j0, n1)[key],
            v10 = at(j1, n0)[key], v11 = at(j1, n1)[key];
      return (1 - tj) * ((1 - tn) * v00 + tn * v01) + tj * ((1 - tn) * v10 + tn * v11);
    };
    return { mean: mix(0), p95: mix(1),
             outsideGrid: J < gj[0] || J > gj[gj.length - 1] || n < gn[0] || n > gn[gn.length - 1] };
  }

  /* Flagging is now UNCONDITIONAL: every sample gets the cut m + z*sigma, and the
   * structure gate that used to precede it has been removed.
   *
   * The gate answered "does a distinct careless subpopulation exist?" and, when it
   * said no, flagged nobody — which bought a guarantee that clean data yields no
   * flags. Measured across six configurations, that guarantee was already gone:
   * it belonged to the BIC branch alone (0.00% flagged on 180 clean simulated
   * datasets, every size), and the excess-tail branch added later gives it up
   * (it opens on 2-6 of 30 clean datasets at short batteries and small n, flagging
   * up to 19% of a clean sample). Keeping the conjunction therefore protected
   * almost nothing while still refusing to act: on real mixtures at a 5% careless
   * rate the gate shut on 13 of 50 samples that genuinely contained careless
   * respondents, costing recall (.585 against .735) with no gain in precision
   * (.892 against .910). Removing it changes nothing at all on the external
   * benchmark (mean MCC .190 either way) and recovers the low-prevalence loss.
   *
   * What replaces the guarantee is an honest number rather than a promise: the
   * caller is told what share of a CLEAN sample of this size the same rule would
   * flag by chance (see cleanEnvelope), so an unremarkable flagged share can be
   * recognised as noise instead of being silently suppressed. */
  function autoFlag(eta, level) {
    // level may be the preset name "standard" (legacy names resolve to it) OR a raw
    // numeric z (power-user override; flagged as out-of-envelope by the UI).
    const z = (typeof level === "number" && isFinite(level)) ? level
              : (SENS_Z[level] || SENS_Z.standard);
    const n = eta.length;
    const lf = leftFit(eta);
    const cut = lf.m + z * lf.sigma;
    const flagged = new Array(n).fill(false);
    for (let i = 0; i < n; i++) flagged[i] = eta[i] > cut;
    const nFlagged = flagged.filter(Boolean).length;
    const estRate = nFlagged / n;
    // one-sided prevalence: attentive mass read below m + 1*sigma, Phi(1)-corrected
    const c1 = lf.m + lf.sigma;
    let nBelow = 0; for (let i = 0; i < n; i++) if (eta[i] < c1) nBelow++;
    /* A clearly negative raw value means the Gaussian-shoulder mapping is out of
     * its domain (skewed eta): report "not estimable" instead of silently 0. */
    const piRaw = 1 - (nBelow / PHI1) / n;
    const piEstimable = piRaw >= -0.02;
    const pi = Math.max(0, Math.min(1, piRaw));
    /* Retained as REPORTED DIAGNOSTICS, no longer as a veto: how far apart a
     * two-component fit puts the groups, and how much heavier the upper tail is
     * than a careless-free sample's would be. Both are informative about whether
     * the flagged set is a real group; neither now decides anything. */
    const g = gaussMix1D(eta);
    const pooledSd = Math.sqrt((g.v1 + g.v2) / 2) || 1;
    const sep = (g.m2 - g.m1) / pooledSd;
    const bicGate = (g.bic1 - g.bic2 > 2) && g.pi > 0.01 && g.pi < 0.7 && sep > 1.0;
    const tailP = binomUpperP(nFlagged, n, normUpper(z));
    const status = nFlagged === 0 ? "none-past-cut" : "flagged";
    /* Above ~40% estimated contamination the geometry degrades (the careless mass
     * competes with the attentive mode for the KDE peak - observed on real data):
     * the estimate is then a rough indication, not a measurement. */
    const heavyContamination = (pi > 0.40 || estRate > 0.40);
    return { status, flagged, estRate, etaThreshold: cut, pi: pi,
             piRaw: piRaw, piEstimable: piEstimable,
             heavyContamination: heavyContamination, separation: sep,
             bicGate: bicGate, tailP: tailP, tailK: nFlagged,
             /* deprecated: the gate no longer decides. Kept so callers that read it
              * still see the two-component diagnostic rather than crashing. */
             twoComp: bicGate || tailP < tailAlpha(n),
             attentiveMode: lf.m, attentiveSigma: lf.sigma,
             sensitivity: (typeof level === "number") ? level : (level || "standard") };
  }

  /* Full ReReRe ensemble score. Returns a per-respondent careless probability
   * plus every component feature. Flagging is AUTOMATIC by default (two-groups
   * model: estimates the careless rate and flags its members, flags ~none on
   * clean data). A manual expected rate (opts.prevalence, a fraction) overrides
   * it with a top-share cut. */
  function ensemble(M, n, J, opts) {
    opts = opts || {};
    const W = opts.weights || WEIGHTS;
    const ef = ensembleFeatures(M, n, J, opts);
    const O = ef.oriented;
    // partner reliability gate (0 = fall back to rr, 1 = full ensemble)
    const g = Math.max(0, Math.min(1, (ef.profileInfo - 0.01) / 0.05));
    const prob = new Float64Array(n), eta = new Float64Array(n);
    for (let i = 0; i < n; i++) {
      const e = W.b0 + W.rr * O.rr[i] + W.irv * O.irv[i] +
                g * (W.longstring * O.longstring[i] + W.person_total * O.person_total[i]) +
                W.d2 * O.d2[i];
      eta[i] = e; prob[i] = 1 / (1 + Math.exp(-e));
    }
    // per-respondent Monte Carlo SE of the log-odds. rr is the only stochastic
    // term: SE(eta) = |w_rr| * SE(z) / MAD(z).
    const zArr = ef.base.z, zSe = ef.base.zSe || zArr.map(() => 0);
    const sz = Float64Array.from(zArr).sort((a, b) => a - b);
    const medz = sz[Math.floor(n / 2)];
    const madz = (Float64Array.from(zArr, v => Math.abs(v - medz)).sort((a, b) => a - b)[Math.floor(n / 2)] * 1.4826) || 1;
    const etaSe = Array.from({ length: n }, (_, i) => Math.abs(W.rr) * zSe[i] / madz);
    let thr, flagged, flagMode, twoComp = null, estRate, flagStatus = null;
    if (opts.prevalence != null) {                       // manual override: top share
      const srt = Array.from(prob).sort((a, b) => a - b);
      const cut = Math.min(n - 1, Math.max(0, Math.floor((1 - opts.prevalence) * n)));
      thr = srt[cut];
      flagged = Array.from(prob, v => v >= thr);
      flagMode = "manual";
      estRate = flagged.filter(Boolean).length / n;
      flagStatus = "manual";
    } else {                                              // automatic attentive-mode cut
      const af = autoFlag(Array.from(eta), opts.sensitivity);
      flagged = af.flagged;
      thr = 1 / (1 + Math.exp(-af.etaThreshold));         // always defined: the cut is unconditional
      flagMode = "auto";
      twoComp = af.twoComp;                               // reported diagnostic only
      estRate = af.estRate;
      flagStatus = af.status;
    }
    // borderline: score within ~2 SE of the flagging threshold (uncertain flag)
    const etaThr = (thr > 0 && thr < 1) ? Math.log(thr / (1 - thr)) : Infinity;
    const borderline = Array.from(eta, (e, i) => isFinite(etaThr) && Math.abs(e - etaThr) < 2 * etaSe[i]);
    return {
      p: Array.from(prob), eta: Array.from(eta), etaSe, borderline,
      nBorderline: borderline.filter(Boolean).length,
      flagged,
      nFlagged: flagged.filter(Boolean).length,
      rr: ef.base.z, irv: ef.base.irv, longstring: ef.base.longstring,
      d2: ef.d2, person_total: ef.base.personTotal, nMissing: ef.base.nMissing,
      oriented: O, diagnostic: ef.base.diagnostic, threshold: thr,
      profileInfo: ef.profileInfo, partnerGate: g,
      estimatedRate: estRate, flagMode: flagMode, twoComponent: twoComp, flagStatus: flagStatus,
      cleanEnvelope: cleanEnvelope(n, J),
      params: Object.assign({}, ef.base.params, {
        weights: W, prevalence: opts.prevalence != null ? opts.prevalence : null,
        profileInfo: ef.profileInfo, partnerGate: g,
        flagMode: flagMode, estimatedRate: estRate, twoComponent: twoComp,
        sensitivity: flagMode === "auto" ? (opts.sensitivity || "standard") : null
      })
    };
  }

  /* ---------------- main entry ----------------------------------------- */
  /* opts: corProp (0.03), iterations (200), minPairs (15), zThreshold (1.5),
   *       seed (1), onProgress (fn 0..1) */
  function run(M, n, J, opts) {
    opts = opts || {};
    const corProp   = opts.corProp    != null ? opts.corProp    : 0.03;
    const iters     = opts.iterations != null ? opts.iterations : 400;
    const minPairs  = opts.minPairs   != null ? opts.minPairs   : 15;
    const zThr      = opts.zThreshold != null ? opts.zThreshold : 1.5;
    const rand      = rng(opts.seed != null ? opts.seed : 1);
    const progress  = opts.onProgress || function () {};

    if (J < 6)  throw new Error("Need at least 6 item columns (found " + J + ").");
    if (n < 10) throw new Error("Need at least 10 respondents (found " + n + ").");

    const { P, minP, nMissing } = prepare(M, n, J);
    progress(0.05);
    const { pa, pb, pr, nPairs } = correlationPairs(P, n, J);
    progress(0.15);

    const k = Math.min(Math.max(minPairs, Math.round(corProp * nPairs)), nPairs);
    const diag = diagnostic(pr, k, n);

    // coupled pairs = top-k by |r|
    const order = Array.from({ length: nPairs }, (_, i) => i)
      .sort((x, y) => Math.abs(pr[y]) - Math.abs(pr[x]));
    const topIdx = Int32Array.from(order.slice(0, k));
    const coupled = indCors(P, minP, n, J, pa, pb, pr, topIdx);

    // permutation baseline
    const sum = new Float64Array(n), sumSq = new Float64Array(n);
    for (let b = 0; b < iters; b++) {
      const rc = indCors(P, minP, n, J, pa, pb, pr, sampleIdx(nPairs, k, rand));
      for (let i = 0; i < n; i++) { sum[i] += rc[i]; sumSq[i] += rc[i] * rc[i]; }
      if ((b & 7) === 0) progress(0.15 + 0.8 * (b / iters));
    }
    const z = new Float64Array(n);
    let minValid = Infinity;
    for (let i = 0; i < n; i++) {
      const mu = sum[i] / iters;
      const va = sumSq[i] / iters - mu * mu;
      const sd = Math.sqrt(Math.max(va, 0));
      z[i] = sd > 1e-12 ? (coupled[i] - mu) / sd : NaN;
      if (Number.isFinite(z[i]) && z[i] < minValid) minValid = z[i];
    }
    // degenerate (near-constant) respondents: no variance -> no coherence;
    // assign the most-careless observed z (they are straight-liners, caught
    // by LongString anyway) instead of NaN.
    for (let i = 0; i < n; i++) if (!Number.isFinite(z[i])) z[i] = minValid;
    // per-respondent Monte Carlo standard error of the permutation z-score
    const zSe = Array.from(z, v => Math.sqrt((1 + v * v / 2) / iters));

    const aux = auxiliaries(P, n, J);
    const flagged = Array.from(z, v => v <= zThr);
    progress(1);

    return {
      z: Array.from(z),
      zSe,
      coupled: Array.from(coupled),
      flagged,
      nFlagged: flagged.filter(Boolean).length,
      irv: Array.from(aux.irv),
      personTotal: Array.from(aux.personTotal),
      longstring: Array.from(aux.longstring),
      nMissing: Array.from(nMissing),
      diagnostic: diag,
      params: { corProp, iterations: iters, minPairs, zThreshold: zThr, k, nPairs, n, J }
    };
  }

  return { parseCSV, detectColumns, buildMatrix, run, ensemble, ensembleFeatures, leftFit,
           autoFlag, gaussMix1D, mahalanobisD2, robustZ, toNum, _rng: rng, _jacobiEigen: jacobiEigen,
           get WEIGHTS() { return WEIGHTS; }, set WEIGHTS(w) { WEIGHTS = w; } };
});
