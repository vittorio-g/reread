/* test_leftfit_fixes.js — the one-sided fit occasionally collapses on small samples: the shipped
 * rule takes the LEFT-MOST kde peak above 15% of the maximum, so a spurious bump in the lower tail
 * wins over the dominant attentive mode; sigma is then estimated from the handful of points below
 * it and the cut collapses. Two candidate fixes, both tested against the shipped behaviour:
 *   A  guard   : keep the fit, but abstain when the implied flag share is absurd (> LIMIT)
 *   B1 global  : take the TALLEST kde peak (capped at the median) instead of the left-most one
 *   B2 hsm     : half-sample mode (classic robust mode estimator)
 * Reported on clean simulated data (false alarms) and on the six real datasets (detection).
 */
const R = require("./site/reread.js"); const W = R.WEIGHTS; const fs = require("fs");
const Z = 2.5, LIMIT = 0.25;
let SEED = 424242; const rng = () => { SEED = (SEED*1103515245+12345)&0x7fffffff; return SEED/0x7fffffff; };
function gauss(){let u=0,v=0;while(!u)u=rng();while(!v)v=rng();return Math.sqrt(-2*Math.log(u))*Math.cos(2*Math.PI*v);}
const K5=5, THR=[-0.84,-0.25,0.25,0.84];
function genClean(F,ipf,n){const J=F*ipf;const load=[],fac=[],tau=[];
  for(let j=0;j<J;j++){let l=0.4+0.4*rng();if(rng()<0.4)l=-l;load.push(l);fac.push(Math.floor(j/ipf));tau.push(0.45*gauss());}
  const rho=0.3,sq=Math.sqrt(rho),sq1=Math.sqrt(1-rho);const rows=[];
  for(let i=0;i<n;i++){const c=gauss(),f=[];for(let k=0;k<F;k++)f.push(sq*c+sq1*gauss());
    const row=new Array(J);
    for(let j=0;j<J;j++){let y=tau[j]+load[j]*f[fac[j]];y+=Math.sqrt(Math.max(0.05,1-load[j]*load[j]))*gauss();
      let cat=1;for(let q=0;q<K5-1;q++)if(y>THR[q])cat=q+2;row[j]=cat;}
    rows.push(row);}
  return {rows,J};}
const toMat=(rows,J)=>{const n=rows.length,m=new Float64Array(n*J);for(let i=0;i<n;i++)for(let j=0;j<J;j++)m[i*J+j]=rows[i][j];return m;};
function etaOf(mat,n,J,it){const ef=R.ensembleFeatures(mat,n,J,{iterations:it,seed:1}),O=ef.oriented;
  const g=Math.max(0,Math.min(1,(ef.profileInfo-0.01)/0.05));
  const e=new Array(n);for(let i=0;i<n;i++)e[i]=W.b0+W.rr*O.rr[i]+g*(W.longstring*O.longstring[i]+W.person_total*O.person_total[i]);
  return e;}
/* ---- kde shared by the mode variants ---- */
function kde(eta){const n=eta.length,s=[...eta].sort((a,b)=>a-b);
  const q1=s[Math.floor(.25*n)],q3=s[Math.floor(.75*n)],iqr=q3-q1;
  let mu=0;for(const v of eta)mu+=v;mu/=n; let va=0;for(const v of eta)va+=(v-mu)*(v-mu);va/=n;
  const sd=Math.sqrt(va); let h=0.9*Math.min(sd||1,(iqr/1.34)||(sd||1))*Math.pow(n,-0.2);
  if(!(h>0))h=(sd||1)*0.3;
  const G=512,lo=s[0]-h,hi=s[n-1]+h,dens=new Float64Array(G);
  for(let g=0;g<G;g++){const xg=lo+(hi-lo)*g/(G-1);let d=0;
    for(let i=0;i<n;i++){const z=(xg-eta[i])/h;d+=Math.exp(-0.5*z*z);}dens[g]=d;}
  return {s,dens,lo,hi,G,sd};}
function sigmaBelow(eta,m,sd){const below=[];for(const v of eta)if(v<m)below.push(m-v);
  below.sort((a,b)=>a-b); let sig=below.length?1.4826*below[Math.floor(below.length/2)]:0;
  if(!(sig>0))sig=(sd||1e-9); return {sigma:sig,nBelow:below.length};}
function fitShipped(eta){const {s,dens,lo,hi,G,sd}=kde(eta);let dmax=0;for(let g=0;g<G;g++)if(dens[g]>dmax)dmax=dens[g];
  let m=s[Math.floor(eta.length/2)];
  for(let g=1;g<G-1;g++) if(dens[g]>=dens[g-1]&&dens[g]>=dens[g+1]&&dens[g]>0.15*dmax){m=lo+(hi-lo)*g/(G-1);break;}
  return {m,...sigmaBelow(eta,m,sd)};}
function fitGlobal(eta){const {s,dens,lo,hi,G,sd}=kde(eta);let best=0,bi=0;
  for(let g=0;g<G;g++)if(dens[g]>best){best=dens[g];bi=g;}
  let m=lo+(hi-lo)*bi/(G-1); const med=s[Math.floor(eta.length/2)];
  if(m>med)m=med;                                    // never above the median
  return {m,...sigmaBelow(eta,m,sd)};}
function fitHSM(eta){let s=[...eta].sort((a,b)=>a-b);
  while(s.length>3){const k=Math.ceil(s.length/2);let bw=Infinity,bi=0;
    for(let i=0;i+k-1<s.length;i++){const w=s[i+k-1]-s[i];if(w<bw){bw=w;bi=i;}}
    s=s.slice(bi,bi+k);}
  const m=s.reduce((a,b)=>a+b,0)/s.length;
  const {sd}=kde(eta); return {m,...sigmaBelow(eta,m,sd)};}
const FITS={shipped:fitShipped,global:fitGlobal,hsm:fitHSM};
function flagFrac(eta,fit){const f=fit(eta);const tau=f.m+Z*f.sigma;
  return eta.reduce((s,v)=>s+(v>tau?1:0),0)/eta.length;}
function mcc(flag,y){let tp=0,tn=0,fp=0,fn=0;
  for(let i=0;i<y.length;i++){if(flag[i]){y[i]?tp++:fp++;}else{y[i]?fn++:tn++;}}
  const d=Math.sqrt((tp+fp)*(tp+fn)*(tn+fp)*(tn+fn));return d?(tp*tn-fp*fn)/d:0;}

console.log("=== CLEAN simulated data: share flagged (mean / worst) and how often it is absurd ===");
console.log("config              | shipped mean/worst  >25%  | global mean/worst  >25%  | hsm mean/worst  >25%");
for(const [F,ipf,N,REPS] of [[10,6,200,40],[30,6,500,40],[10,6,500,30],[20,6,2000,15]]){
  const J=F*ipf; const acc={shipped:[],global:[],hsm:[]};
  for(let r=0;r<REPS;r++){ const {rows}=genClean(F,ipf,N); const eta=etaOf(toMat(rows,J),N,J,100);
    for(const k in FITS) acc[k].push(flagFrac(eta,FITS[k])); }
  const fmt=a=>{const mean=a.reduce((s,v)=>s+v,0)/a.length,worst=Math.max(...a),bad=a.filter(v=>v>LIMIT).length;
    return (100*mean).toFixed(1).padStart(5)+"/"+(100*worst).toFixed(1).padStart(5)+"  "+String(bad).padStart(2)+"/"+a.length;};
  console.log(("J="+J+" n="+N).padEnd(20)+"| "+fmt(acc.shipped)+" | "+fmt(acc.global)+" | "+fmt(acc.hsm));
}
console.log("\n=== REAL datasets: share flagged and MCC (GT >=2 checks; TISP >=1) ===");
const D="../Dataset/gt_benchmark_candidates";
const SETS=[["TISP",D+"/tisp","_labels.csv"],["TUMI",D+"/tumi","_labels2.csv"],
  ["warnIPIP",D+"/warning_ipipneo300","_labels2.csv"],["smarvus",D+"/smarvus","_labels2.csv"],
  ["KayS2",D+"/kay_idris_idria/s2","_labels2.csv"],["KayS1",D+"/kay_idris_idria/s1","_labels2.csv"]];
console.log("dataset     prev  | shipped flag/MCC | global flag/MCC | hsm flag/MCC");
for(const [nm,dir,lf] of SETS){
  const P=R.parseCSV(fs.readFileSync(dir+"/_matrix.csv","utf8")),J=P.header.length;
  const rows=P.rows.map(r=>r.map(v=>{const x=Number(v);return Number.isFinite(x)?x:NaN;}));
  const y=R.parseCSV(fs.readFileSync(dir+"/"+lf,"utf8")).rows.map(r=>Number(r[0]));
  const n=Math.min(rows.length,y.length);
  const eta=etaOf(toMat(rows.slice(0,n),J),n,J,200);
  const cells=[];
  for(const k in FITS){ const f=FITS[k](eta), tau=f.m+Z*f.sigma, fl=eta.map(v=>v>tau);
    cells.push((100*fl.filter(Boolean).length/n).toFixed(1).padStart(5)+"% "+mcc(fl,y).toFixed(3)); }
  console.log(nm.padEnd(11)+(100*y.slice(0,n).reduce((s,v)=>s+v,0)/n).toFixed(1).padStart(5)+"% | "+cells.join(" | "));
}
