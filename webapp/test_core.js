/* test_core.js — Node validation of reread.js
 * 1. Synthetic: factor-structured data + injected random careless -> z must
 *    separate them (AUC > 0.8) and diagnostic must say ok/marginal.
 * 2. No-structure data -> diagnostic must warn (weak/no_structure).
 * 3. Fidelity: run on the real study matrix (exported CSV) and correlate with
 *    the R pipeline's z_rr from features.csv (Spearman > 0.9 expected).
 */
const fs = require("fs");
const path = require("path");
const R = require("./site/reread.js");

function mulberry(seed) { return R._rng(seed); }

/* ---------- 1. synthetic factor data + careless ---------- */
function synth(n, nF, ipf, pctCareless, seed) {
  const rand = mulberry(seed);
  const J = nF * ipf;
  const gauss = () => {
    let u = 0, v = 0;
    while (u === 0) u = rand();
    while (v === 0) v = rand();
    return Math.sqrt(-2 * Math.log(u)) * Math.cos(2 * Math.PI * v);
  };
  const M = new Float64Array(n * J);
  const careless = new Array(n).fill(false);
  for (let i = 0; i < n; i++) {
    const f = Array.from({ length: nF }, gauss);
    const isCare = rand() < pctCareless;
    careless[i] = isCare;
    for (let j = 0; j < J; j++) {
      let x;
      if (isCare) x = 1 + Math.floor(rand() * 7);           // uniform random 1..7
      else {
        const load = 0.65;
        const latent = load * f[Math.floor(j / ipf)] + Math.sqrt(1 - load * load) * gauss();
        x = Math.max(1, Math.min(7, Math.round(4 + 1.5 * latent)));
      }
      M[i * J + j] = x;
    }
  }
  return { M, n, J, careless };
}

function auc(scoreLowIsPositive, labels) {
  const s = scoreLowIsPositive.map((v, i) => [-v, labels[i] ? 1 : 0])
    .sort((a, b) => a[0] - b[0]);
  let rank = 1, sumR = 0, n1 = 0;
  for (const [, y] of s.map((x, i) => [i, x[1]])) { /* placeholder */ }
  // simple rank-based AUC
  const arr = scoreLowIsPositive.map((v, i) => ({ v: -v, y: labels[i] ? 1 : 0 }))
    .sort((a, b) => a.v - b.v);
  let r = 0, npos = 0, nneg = 0;
  arr.forEach((o, i) => { if (o.y) { r += i + 1; npos++; } else nneg++; });
  return (r - npos * (npos + 1) / 2) / (npos * nneg);
}

let pass = 0, fail = 0;
function check(name, cond, detail) {
  if (cond) { pass++; console.log("  PASS " + name + (detail ? "  (" + detail + ")" : "")); }
  else { fail++; console.log("  FAIL " + name + (detail ? "  (" + detail + ")" : "")); }
}

console.log("== 1. synthetic structured data (n=300, 10 factors x 8 items, 20% careless) ==");
{
  const { M, n, J, careless } = synth(300, 10, 8, 0.20, 42);
  const res = R.run(M, n, J, { iterations: 100, seed: 7 });
  const a = auc(res.z, careless);
  check("AUC > 0.80", a > 0.80, "AUC=" + a.toFixed(3));
  check("diagnostic ok/marginal", ["ok", "marginal"].includes(res.diagnostic.level),
        "level=" + res.diagnostic.level + " ratio=" + res.diagnostic.separation_ratio.toFixed(2));
  const zc = res.z.filter((_, i) => careless[i]);
  const zg = res.z.filter((_, i) => !careless[i]);
  const mean = v => v.reduce((s, x) => s + x, 0) / v.length;
  check("mean z careless < mean z careful", mean(zc) < mean(zg),
        mean(zc).toFixed(2) + " vs " + mean(zg).toFixed(2));
  check("flags mostly careless at z<=1.5",
        res.flagged.filter((f, i) => f && careless[i]).length >
        res.flagged.filter((f, i) => f && !careless[i]).length,
        res.nFlagged + " flagged");
}

console.log("== 2. structureless data -> diagnostic must warn ==");
{
  const rand = mulberry(3);
  const n = 200, J = 40;
  const M = new Float64Array(n * J);
  for (let i = 0; i < n * J; i++) M[i] = 1 + Math.floor(rand() * 7);
  const res = R.run(M, n, J, { iterations: 60, seed: 7 });
  check("level weak/no_structure", ["weak", "no_structure"].includes(res.diagnostic.level),
        "level=" + res.diagnostic.level + " ratio=" + res.diagnostic.separation_ratio.toFixed(2));
}

console.log("== 3. CSV parser ==");
{
  const t1 = R.parseCSV('id,a,b\np1,"1,5",2\np2,3,4\n');
  check("comma+quotes", t1.rows.length === 2 && t1.rows[0][1] === "1,5");
  const t2 = R.parseCSV("id;a;b\np1;1;2\np2;3;4");
  check("semicolon sniff", t2.delimiter === ";" && t2.rows.length === 2);
  const t3 = R.parseCSV("a\tb\n1\t2");
  check("tab sniff", t3.delimiter === "\t");
  check("EU decimal", R.toNum("3,5") === 3.5 && R.toNum("4.2") === 4.2);
}

console.log("== 4. fidelity vs R pipeline on real study data ==");
{
  const mPath = path.join(__dirname, "_study_matrix.csv");
  const fPath = path.join(__dirname, "_study_zrr.csv");
  if (fs.existsSync(mPath) && fs.existsSync(fPath)) {
    const parsed = R.parseCSV(fs.readFileSync(mPath, "utf8"));
    const det = R.detectColumns(parsed.header, parsed.rows);
    const { M, n, J } = R.buildMatrix(parsed.rows, det.itemCols);
    const res = R.run(M, n, J, { iterations: 200, seed: 11 });
    const rz = R.parseCSV(fs.readFileSync(fPath, "utf8")).rows.map(r => Number(r[0]));
    // Spearman
    const rankv = v => {
      const idx = v.map((x, i) => [x, i]).sort((a, b) => a[0] - b[0]);
      const rk = new Array(v.length);
      idx.forEach(([, i], r) => rk[i] = r);
      return rk;
    };
    const ra = rankv(res.z), rb = rankv(rz);
    const mean = v => v.reduce((s, x) => s + x, 0) / v.length;
    const ma = mean(ra), mb = mean(rb);
    let sab = 0, sa = 0, sb = 0;
    for (let i = 0; i < ra.length; i++) {
      sab += (ra[i] - ma) * (rb[i] - mb);
      sa += (ra[i] - ma) ** 2; sb += (rb[i] - mb) ** 2;
    }
    const rho = sab / Math.sqrt(sa * sb);
    check("Spearman(z_JS, z_R) > 0.9", rho > 0.9, "rho=" + rho.toFixed(3) + " n=" + n + " J=" + J);
  } else {
    console.log("  SKIP (export _study_matrix.csv/_study_zrr.csv first)");
  }
}

console.log("\n" + pass + " passed, " + fail + " failed");
process.exit(fail ? 1 : 0);
