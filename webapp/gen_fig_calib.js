/* gen_fig_calib.js — data for the new Fig 3: pi-hat (mixture prevalence, deployed autoFlag semantics)
 * vs true rate, with a 5-95% coverage band. 300 resamples per achievable rate k/(70+k). */
const R=require("./site/reread.js"); const W=R.WEIGHTS; const fs=require("fs"); const rng=R._rng(31459);
const Ms=R.parseCSV(fs.readFileSync("_study_matrix.csv","utf8")), J=Ms.header.length;
const allR=Ms.rows.map(r=>r.map(v=>{const x=Number(v);return Number.isFinite(x)?x:NaN;}));
const y=R.parseCSV(fs.readFileSync("_study_labels.csv","utf8")).rows.map(r=>Number(r[0]));
const cf=y.map((v,i)=>v===0?i:-1).filter(i=>i>=0), cl=y.map((v,i)=>v===1?i:-1).filter(i=>i>=0);
function sampleK(pool,k){const a=pool.slice();for(let i=0;i<k;i++){const j=i+Math.floor(rng()*(a.length-i));const t=a[i];a[i]=a[j];a[j]=t;}return a.slice(0,k);}
function gateOK(gm){const pooled=Math.sqrt((gm.v1+gm.v2)/2)||1,sep=(gm.m2-gm.m1)/pooled;return (gm.bic1-gm.bic2>2)&&gm.pi>0.01&&gm.pi<0.7&&sep>1.0;}
function eta(idx){const n=idx.length,mat=new Float64Array(n*J);for(let i=0;i<n;i++)for(let j=0;j<J;j++)mat[i*J+j]=allR[idx[i]][j];
  const ef=R.ensembleFeatures(mat,n,J,{seed:1,iterations:200}),O=ef.oriented,g=Math.max(0,Math.min(1,(ef.profileInfo-0.01)/0.05));
  const e=new Array(n);for(let i=0;i<n;i++)e[i]=W.b0+W.rr*O.rr[i]+g*(W.longstring*O.longstring[i]+W.person_total*O.person_total[i]);return e;}
const pct=(a,q)=>{const s=[...a].sort((x,y)=>x-y);return s[Math.min(s.length-1,Math.max(0,Math.round(q*(s.length-1))))];};
const mn=a=>a.reduce((s,v)=>s+v,0)/a.length, NC=cf.length, REPS=300;
const KS=[4,6,8,10,12,15,18,21,24,28,33,40,47,58]; // rates ~5% .. ~45%
// Two quantities per rate: (a) pi-hat, the mixture's estimate of how many careless are PRESENT;
// (b) flag_*, the share of the sample the shipped Standard cut actually REMOVES (autoFlag.estRate).
const out=[["true","pi_mean","pi_p5","pi_p95","gate_on","flag_mean","flag_p5","flag_p95"]];
console.log("true |  pi_mean  [p5,p95]   flagged  [p5,p95]   gateOn");
for(const k of KS){ if(k>cl.length) break; const tr=k/(NC+k); const P=[], Fl=[]; let gon=0;
  for(let rep=0;rep<REPS;rep++){const idx=[...cf,...sampleK(cl,k)];const e=eta(idx);
    const af=R.autoFlag(e,"standard");                   // gate + one-sided estimator + shipped cut
    if(af.twoComp)gon++; P.push(af.pi); Fl.push(af.estRate);}
  const m=mn(P),lo=pct(P,0.05),hi=pct(P,0.95);
  const fm=mn(Fl),flo=pct(Fl,0.05),fhi=pct(Fl,0.95);
  out.push([tr.toFixed(4),m.toFixed(4),lo.toFixed(4),hi.toFixed(4),(gon/REPS).toFixed(3),
            fm.toFixed(4),flo.toFixed(4),fhi.toFixed(4)]);
  console.log((100*tr).toFixed(1).padStart(5)+"% |  "+(100*m).toFixed(1).padStart(5)+"%  ["+(100*lo).toFixed(0)+","+(100*hi).toFixed(0)+"]   "+
              (100*fm).toFixed(1).padStart(5)+"%  ["+(100*flo).toFixed(0)+","+(100*fhi).toFixed(0)+"]   "+(100*gon/REPS).toFixed(0)+"%");
}
fs.writeFileSync("fig_calib_pi.csv",out.map(r=>r.join(",")).join("\n"));
console.log("\nwrote fig_calib_pi.csv");
