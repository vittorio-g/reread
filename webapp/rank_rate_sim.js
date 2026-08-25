/* rank_rate_sim.js — Study 3 (simulation): ranking robustness vs careless share.
 * Cutoff-free (AUC only). Unlike the Study 1 mixture, here N is held CONSTANT (n=500) and
 * only the careless proportion changes, so the rate effect is isolated from sample size.
 * Engine recomputed per dataset (pair selection / reference profile are sample-level, and
 * that is exactly what contamination attacks). Generator ported from sim_rr_scaling.js. */
const fs = require("fs");
const R = require("./site/rerere.js");
const rng = R._rng(20260728);

function gauss(){let u=0,v=0;while(u===0)u=rng();while(v===0)v=rng();return Math.sqrt(-2*Math.log(u))*Math.cos(2*Math.PI*v);}
function qnorm(p){
  const a=[-3.969683028665376e+01,2.209460984245205e+02,-2.759285104469687e+02,1.383577518672690e+02,-3.066479806614716e+01,2.506628277459239e+00];
  const b=[-5.447609879822406e+01,1.615858368580409e+02,-1.556989798598866e+02,6.680131188771972e+01,-1.328068155288572e+01];
  const c=[-7.784894002430293e-03,-3.223964580411365e-01,-2.400758277161838e+00,-2.549732539343734e+00,4.374664141464968e+00,2.938163982698783e+00];
  const d=[7.784695709041462e-03,3.224671290700398e-01,2.445134137142996e+00,3.754408661907416e+00];
  const pl=0.02425;let q,r;
  if(p<pl){q=Math.sqrt(-2*Math.log(p));return (((((c[0]*q+c[1])*q+c[2])*q+c[3])*q+c[4])*q+c[5])/((((d[0]*q+d[1])*q+d[2])*q+d[3])*q+1);}
  if(p<=1-pl){q=p-0.5;r=q*q;return (((((a[0]*r+a[1])*r+a[2])*r+a[3])*r+a[4])*r+a[5])*q/(((((b[0]*r+b[1])*r+b[2])*r+b[3])*r+b[4])*r+1);}
  q=Math.sqrt(-2*Math.log(1-p));return -(((((c[0]*q+c[1])*q+c[2])*q+c[3])*q+c[4])*q+c[5])/((((d[0]*q+d[1])*q+d[2])*q+d[3])*q+1);
}
const K=5, THR=[]; for(let c=1;c<K;c++) THR.push(qnorm(c/K));
function genClean(F,ipf,n){
  const J=F*ipf, load=[],fac=[],tau=[];
  for(let j=0;j<J;j++){let l=0.4+0.4*rng(); if(rng()<0.4) l=-l; load.push(l); fac.push(Math.floor(j/ipf)); tau.push(0.45*gauss());}
  const rho=0.3,sq=Math.sqrt(rho),sq1=Math.sqrt(1-rho); const rows=[];
  for(let i=0;i<n;i++){
    const common=gauss(),f=[];for(let k=0;k<F;k++)f.push(sq*common+sq1*gauss());
    const row=new Array(J);
    for(let j=0;j<J;j++){
      const yv=tau[j]+load[j]*f[fac[j]]+Math.sqrt(1-load[j]*load[j])*gauss();
      let cat=1;for(let c=0;c<K-1;c++)if(yv>THR[c])cat=c+2;row[j]=cat;
    }
    rows.push(row);
  }
  return {rows,J};
}
/* six careless patterns at graded corruption, as in app.js/eval_method.js */
function inject(rows,J,pct){
  const N=rows.length, LV=[1,2,3,4,5];
  const ri=n=>Math.floor(rng()*n), pick=a=>a[ri(a.length)];
  const sk=k=>{const a=[...Array(J).keys()];for(let i=0;i<k;i++){const j=i+ri(J-i);[a[i],a[j]]=[a[j],a[i]];}return a.slice(0,k);};
  const out=rows.map(r=>r.slice()), lab=new Array(N).fill(0);
  const ids=[...Array(N).keys()]; for(let i=N-1;i>0;i--){const j=ri(i+1);[ids[i],ids[j]]=[ids[j],ids[i]];}
  const M=Math.round(N*pct), pats=["random","long","straight","acq","mixed","fatigue"];
  for(let m=0;m<M;m++){
    const i=ids[m], lvl=[.5,.6,.7,.8,.9,1][ri(6)], k=Math.max(1,Math.round(J*lvl)), pat=pats[m%6], row=out[i];
    if(pat==="random") sk(k).forEach(p=>row[p]=pick(LV));
    else if(pat==="straight"){const v=pick(LV); sk(k).forEach(p=>row[p]=v);}
    else if(pat==="acq") sk(k).forEach(p=>row[p]=Math.max(1,Math.min(5,4+pick([-1,0,0,0,1]))));
    else if(pat==="long"){const idx=sk(k);let v=pick(LV);idx.forEach((p,t)=>{if(t%12===0)v=pick(LV);row[p]=v;});}
    else if(pat==="mixed"){const idx=sk(k),h=Math.round(k/2);let v=pick(LV);idx.slice(0,h).forEach((p,t)=>{if(t%10===0)v=pick(LV);row[p]=v;});idx.slice(h).forEach(p=>row[p]=pick(LV));}
    else {const idx=[];for(let t=J-k;t<J;t++)idx.push(t);let v=pick(LV);idx.forEach((p,t)=>{if(t%10===0)v=pick(LV);row[p]=v;});}
    lab[i]=1;
  }
  return {out,lab};
}
function auc(sc,lab){
  const idx=[...sc.keys()].filter(i=>Number.isFinite(sc[i])).sort((a,b)=>sc[a]-sc[b]);
  const rk=[]; idx.forEach((id,r)=>rk[id]=r+1);
  let n1=0,n0=0,s=0; for(const i of idx){ if(lab[i]===1){n1++;s+=rk[i];} else n0++; }
  return n1&&n0?(s-n1*(n1+1)/2)/(n1*n0):NaN;
}
const mean=v=>v.reduce((s,x)=>s+x,0)/v.length;
const sd=v=>{const m=mean(v);return Math.sqrt(mean(v.map(x=>(x-m)**2)));};

const N=500, REPS=10;
const CONFIGS=[{F:20,ipf:6,tag:"120 items (20x6)"},{F:10,ipf:6,tag:"60 items (10x6)"}];
const out=[["items","true_rate","auc_rc","sd_rc","auc_ens","sd_ens","auc_ls","auc_pt"]];
for(const cfg of CONFIGS){
  console.log("\n=== " + cfg.tag + ", n=" + N + ", " + REPS + " reps/rate (AUC only) ===");
  console.log("rate | AUC rc         | AUC ensemble   | AUC LongString | AUC PersonTotal");
  for(const rate of [0.05,0.10,0.20,0.30,0.40,0.50,0.60,0.70,0.80,0.90]){
    const arc=[],aen=[],als=[],apt=[];
    for(let rep=0;rep<REPS;rep++){
      const {rows,J}=genClean(cfg.F,cfg.ipf,N);
      const {out:data,lab}=inject(rows,J,rate);
      const m=new Float64Array(N*J);
      for(let i=0;i<N;i++)for(let j=0;j<J;j++)m[i*J+j]=data[i][j];
      const res=R.ensemble(m,N,J,{iterations:150,seed:500+rep});
      arc.push(auc(res.rr.map(v=>-v),lab));
      aen.push(auc(res.p,lab));
      als.push(auc(res.oriented.longstring,lab));
      apt.push(auc(res.oriented.person_total,lab));
    }
    console.log((100*rate).toFixed(0).padStart(3)+"% | "+mean(arc).toFixed(3)+" ("+sd(arc).toFixed(3)+") | "+
      mean(aen).toFixed(3)+" ("+sd(aen).toFixed(3)+") | "+mean(als).toFixed(3)+"          | "+mean(apt).toFixed(3));
    out.push([cfg.F*cfg.ipf,rate,mean(arc).toFixed(4),sd(arc).toFixed(4),mean(aen).toFixed(4),sd(aen).toFixed(4),mean(als).toFixed(4),mean(apt).toFixed(4)]);
  }
}
fs.writeFileSync("rank_rate_sim.csv",out.map(r=>r.join(",")).join("\n"));
console.log("\nwrote rank_rate_sim.csv");
