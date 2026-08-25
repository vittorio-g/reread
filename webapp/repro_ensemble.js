/* repro_ensemble.js — reproducibility & idempotence of the ReReRe ENSEMBLE.
 *  A. Determinism: same data+seed -> identical ensemble probability.
 *  B. Monte-Carlo stability across seeds (p SD; flag-set stability).
 *  C. Idempotence: clean (remove flagged), reload, re-flag -> converge?
 *  C2. False-positive floor on a verified-careful-only set.
 * Flagging uses an ABSOLUTE eta cutoff calibrated once on the full study at its
 * careless-rate operating point (analogous to the rr tool's fixed z=1.5). */
const fs = require("fs");
const R = require("./site/reread.js");
const M = R.parseCSV(fs.readFileSync("_study_matrix.csv","utf8"));
const L = R.parseCSV(fs.readFileSync("_study_labels.csv","utf8"));
const allRows = M.rows.map(r=>r.map(Number)), y = L.rows.map(r=>Number(r[0]));
const N = allRows.length, J = M.header.length;
const RATE = y.reduce((s,v)=>s+v,0)/N;
function mat(idx){const n=idx.length,m=new Float64Array(n*J);for(let i=0;i<n;i++)for(let j=0;j<J;j++)m[i*J+j]=allRows[idx[i]][j];return {m,n};}
function ens(idx,opts){const {m,n}=mat(idx);return R.ensemble(m,n,J,Object.assign({iterations:300,seed:1,skipD2:true},opts));}
const mean=v=>v.reduce((s,x)=>s+x,0)/v.length, sd=v=>{const m=mean(v);return Math.sqrt(mean(v.map(x=>(x-m)**2)));};
const quantile=(v,q)=>{const s=[...v].sort((a,b)=>a-b);return s[Math.min(s.length-1,Math.floor(q*s.length))];};

// absolute eta cutoff from full study at the careless-rate operating point
const full = ens([...Array(N).keys()],{seed:1});
const TAU = quantile(full.eta, 1-RATE);
console.log("Ensemble flag rule: eta >= "+TAU.toFixed(3)+" (full-study operating point at rate "+RATE.toFixed(3)+")\n");

console.log("== A. Determinism (same data, same seed) ==");
{ const a=ens([...Array(N).keys()],{seed:7}).p, b=ens([...Array(N).keys()],{seed:7}).p;
  let mx=0; for(let i=0;i<N;i++) mx=Math.max(mx,Math.abs(a[i]-b[i]));
  console.log("  max |p1-p2| over "+N+" respondents = "+mx.toExponential(2)+(mx===0?"  -> BIT-IDENTICAL":"  -> NOT identical")); }

console.log("\n== B. Monte-Carlo stability across 40 seeds (full data) ==");
{ const idx=[...Array(N).keys()], K=40;
  for(const iters of [50,100,200,300]){
    const pRuns=[],flagRuns=[];
    for(let s=0;s<K;s++){const p=ens(idx,{iterations:iters,seed:1000+s}); pRuns.push(p.p); flagRuns.push(p.eta.map(e=>e>=TAU));}
    const psd=[]; for(let i=0;i<N;i++) psd.push(sd(pRuns.map(p=>p[i])));
    let flips=0; for(let i=0;i<N;i++){const vv=flagRuns.map(f=>f[i]); if(vv.some(x=>x)&&vv.some(x=>!x)) flips++;}
    console.log("  iterations="+String(iters).padEnd(4)+" | mean per-resp p SD = "+mean(psd).toFixed(4)+
      " | flag-unstable = "+flips+"/"+N+" ("+Math.round(100*flips/N)+"%)"); }
}

console.log("\n== C. Idempotence: clean, reload, re-flag? (absolute eta>=TAU) ==");
{ let keep=[...Array(N).keys()]; console.log("  round 0: n="+keep.length);
  for(let round=1;round<=6;round++){
    if(keep.length<15){console.log("  -> stop: n="+keep.length+" too few to re-standardise");break;}
    const res=ens(keep,{iterations:300,seed:42});
    const flagged=[]; keep.forEach((gi,li)=>{ if(res.eta[li]>=TAU) flagged.push(gi); });
    const nCareful=flagged.filter(gi=>y[gi]===0).length;
    console.log("  round "+round+": on n="+keep.length+" -> flags "+flagged.length+
      " ("+(flagged.length-nCareful)+" true-careless, "+nCareful+" careful) | diag="+res.diagnostic.level);
    if(res.diagnostic.level==="no_structure"){console.log("  -> STOP: diagnostic=no_structure — cleaned set has no careless structure left; the tool warns against a further pass (safeguard, not threshold convergence: within-dataset standardisation re-centres each pass).");break;}
    if(flagged.length===0){console.log("  -> CONVERGED: a cleaned set produces no new flags.");break;}
    keep=keep.filter(gi=>flagged.indexOf(gi)<0);
  }
}

console.log("\n== C2. Reload a set of ONLY verified-careful respondents ==");
{ const carefulOnly=[...Array(N).keys()].filter(i=>y[i]===0);
  const res=ens(carefulOnly,{iterations:300,seed:9});
  const flagged=res.eta.filter(e=>e>=TAU).length;
  console.log("  "+carefulOnly.length+" verified-careful respondents, no careless present.");
  console.log("  ensemble flags "+flagged+" of them ("+Math.round(100*flagged/carefulOnly.length)+"%) at eta>=TAU (absolute)");
  // also prevalence-based (what the deployed tool does): top-15%
  const p15=quantile(res.p,0.85); const f15=res.p.filter(v=>v>=p15).length;
  console.log("  (deployable top-15% rule would flag ~"+f15+" = 15% by construction — use a realistic prevalence)"); }
