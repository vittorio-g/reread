/* sweep_injection_conv.js — injection sweep with CONVERGENT-GT careless pools (>=2 independent
 * families) instead of check-only pools. Prediction: certified careless -> slopes closer to 1.
 * Prints band (<=25%) and full OLS slopes on estimable resamples directly. */
const R = require("./site/reread.js"); const W = R.WEIGHTS; const fs = require("fs");
let SEED = 515151; const rng = () => { SEED = (SEED*1103515245+12345)&0x7fffffff; return SEED/0x7fffffff; };
const REPS = 10, ITER = 100, NCAP = 1200, NMIN = 180;
const TARGETS = [0.02,0.05,0.08,0.12,0.16,0.20,0.25,0.30,0.35,0.40];
const D = "../Dataset/gt_benchmark_candidates";
const SETS = [["warning IPIP", D+"/warning_ipipneo300"], ["smarvus", D+"/smarvus"],
              ["Kay S1", D+"/kay_idris_idria/s1"]];
function shuffle(a){for(let i=a.length-1;i>0;i--){const j=Math.floor(rng()*(i+1));const t=a[i];a[i]=a[j];a[j]=t;}return a;}
const out=[["dataset","target","true_prev","rep","piHat","estimable"]];
for (const [nm, dir] of SETS) {
  const P = R.parseCSV(fs.readFileSync(dir+"/_matrix.csv","utf8")), J = P.header.length;
  const rows = P.rows.map(r=>r.map(v=>{const x=Number(v);return Number.isFinite(x)?x:NaN;}));
  const y = R.parseCSV(fs.readFileSync(dir+"/_labels_conv.csv","utf8")).rows.map(r=>Number(r[0]));
  const n0 = Math.min(rows.length,y.length), att=[], car=[];
  for(let i=0;i<n0;i++)(y[i]?car:att).push(rows[i]);
  console.log("\n"+nm+"  J="+J+"  att="+att.length+"  CONV careless="+car.length);
  const pts=[];
  for (const p of TARGETS) {
    let k=Math.floor(Math.min(car.length, att.length*p/(1-p)));
    let m=Math.round(k*(1-p)/p);
    if(m>att.length){m=att.length;k=Math.round(m*p/(1-p));}
    if(k+m>NCAP){const s=NCAP/(k+m);k=Math.max(1,Math.round(k*s));m=Math.round(m*s);}
    if(k<3||k+m<NMIN) continue;
    const truep=k/(k+m); let pb=0,ne=0;
    for(let r=0;r<REPS;r++){
      const A=shuffle(att.slice()).slice(0,m), C=shuffle(car.slice()).slice(0,k);
      const sub=A.concat(C), n=sub.length, mat=new Float64Array(n*J);
      for(let i=0;i<n;i++)for(let j=0;j<J;j++)mat[i*J+j]=sub[i][j];
      const ef=R.ensembleFeatures(mat,n,J,{seed:21+r,iterations:ITER}),O=ef.oriented;
      const g=Math.max(0,Math.min(1,(ef.profileInfo-0.01)/0.05));
      const eta=new Array(n);
      for(let i=0;i<n;i++)eta[i]=W.b0+W.rr*O.rr[i]+g*(W.longstring*O.longstring[i]+W.person_total*O.person_total[i]);
      const af=R.autoFlag(eta,"standard");
      out.push([nm,p.toFixed(2),truep.toFixed(4),r,af.pi.toFixed(4),af.piEstimable?1:0]);
      if(af.piEstimable){pts.push([truep,af.pi]);pb+=af.pi;ne++;}
    }
    console.log("  true="+(100*truep).toFixed(1).padStart(5)+"%  n="+String(k+m).padStart(5)+
      "  piHat="+(ne?(100*pb/ne).toFixed(1):"  na")+"%  estimable "+ne+"/"+REPS);
  }
  const slope=(ps)=>{if(ps.length<8)return null;const mx=ps.reduce((s,p)=>s+p[0],0)/ps.length,my=ps.reduce((s,p)=>s+p[1],0)/ps.length;
    const sxx=ps.reduce((s,p)=>s+(p[0]-mx)**2,0);return sxx>0?ps.reduce((s,p)=>s+(p[0]-mx)*(p[1]-my),0)/sxx:null;};
  const sb=slope(pts.filter(p=>p[0]<=0.25)), sf=slope(pts);
  console.log("  >> SLOPE band(0-25%)="+(sb===null?"na":sb.toFixed(3))+"   full="+(sf===null?"na":sf.toFixed(3)));
}
fs.writeFileSync("injection_curve_conv.csv", out.map(r=>r.join(",")).join("\n"));
console.log("\nwrote injection_curve_conv.csv");
