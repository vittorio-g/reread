/* eval_onesided.js — SOLUTION 1: one-sided (robust) estimation of the ATTENTIVE component.
 *
 * Rationale: the shipped estimator fits BOTH mixture components on contaminated data, so when
 * careless respondents are numerous the "attentive" component absorbs some of them and pi shrinks
 * (the collapse). Here the attentive component is estimated from the LEFT side only, which is
 * (nearly) careless-free by construction, and the careless mass follows by subtraction.
 *
 *   m        = left-most prominent mode of the score distribution (KDE)
 *   mirror   : attentive are symmetric about m  ->  n_att = 2*#{x < m},        pi = 1 - n_att/n
 *   onesigma : sigma from the LEFT side only (MAD below m); cut c = m + sigma,
 *              expected attentive share below c is Phi(1)=.8413 -> n_att = #{x<c}/.8413
 *
 * Baseline: the shipped two-component mixture pi_eta.
 * Bench: Study 1 (real, prevalence swept) + Study 3 (simulated, prevalence x length).
 * Output: onesided_eval.csv
 */
const R=require("./site/reread.js"); const W=R.WEIGHTS; const fs=require("fs");
let SEED=4242; const rng=()=>{SEED=(SEED*1103515245+12345)&0x7fffffff;return SEED/0x7fffffff;};
function gauss(){let u=0,v=0;while(!u)u=rng();while(!v)v=rng();return Math.sqrt(-2*Math.log(u))*Math.cos(2*Math.PI*v);}
const K=5, THR=[-0.84,-0.25,0.25,0.84];
function genClean(F,ipf,n){const J=F*ipf;const load=[],fac=[],tau=[];
  for(let j=0;j<J;j++){let l=0.4+0.4*rng();if(rng()<0.4)l=-l;load.push(l);fac.push(Math.floor(j/ipf));tau.push(0.45*gauss());}
  const rho=0.3,sq=Math.sqrt(rho),sq1=Math.sqrt(1-rho);const rows=[];
  for(let i=0;i<n;i++){const c=gauss(),f=[];for(let k=0;k<F;k++)f.push(sq*c+sq1*gauss());
    const row=new Array(J);
    for(let j=0;j<J;j++){let y=tau[j]+load[j]*f[fac[j]];y+=Math.sqrt(Math.max(0.05,1-load[j]*load[j]))*gauss();
      let cat=1;for(let q=0;q<K-1;q++)if(y>THR[q])cat=q+2;row[j]=cat;}
    rows.push(row);}
  return {rows,J};}
function inject(rows,J,pct){const N=rows.length,out=rows.map(r=>r.slice()),lab=new Array(N).fill(0);
  const M=Math.round(N*pct),ord=[...Array(N).keys()];
  for(let i=0;i<M;i++){const j=i+Math.floor(rng()*(N-i));const t=ord[i];ord[i]=ord[j];ord[j]=t;}
  const pats=["random","longstring","pure_straight","acquiescent","mixed","fatigue"];
  for(let m=0;m<M;m++){const i=ord[m],pat=pats[m%6],lvl=0.5+0.5*rng(),k=Math.max(1,Math.round(J*lvl)),row=out[i];
    const idx=[...Array(J).keys()];for(let t=0;t<k;t++){const j=t+Math.floor(rng()*(J-t));const tm=idx[t];idx[t]=idx[j];idx[j]=tm;}
    const sel=idx.slice(0,k);
    if(pat==="pure_straight"||pat==="longstring"){const v=1+Math.floor(rng()*K);sel.forEach(p=>row[p]=v);}
    else if(pat==="acquiescent"){sel.forEach(p=>row[p]=Math.max(1,Math.min(K,4+Math.round(rng()))));}
    else if(pat==="fatigue"){for(let t=J-k;t<J;t++)row[t]=1+Math.floor(rng()*K);}
    else sel.forEach(p=>row[p]=1+Math.floor(rng()*K));
    lab[i]=1;}
  return {out,lab};}
const toMat=(rows,J)=>{const n=rows.length,m=new Float64Array(n*J);for(let i=0;i<n;i++)for(let j=0;j<J;j++)m[i*J+j]=rows[i][j];return m;};
function gateOK(gm){const p=Math.sqrt((gm.v1+gm.v2)/2)||1;return (gm.bic1-gm.bic2>2)&&gm.pi>0.01&&gm.pi<0.7&&(gm.m2-gm.m1)/p>1.0;}
const med=a=>{const s=[...a].sort((x,y)=>x-y);const h=s.length>>1;return s.length%2?s[h]:(s[h-1]+s[h])/2;};
/* left-most prominent KDE mode */
function leftMode(x){const n=x.length,s=[...x].sort((a,b)=>a-b);
  const q1=s[Math.floor(.25*n)],q3=s[Math.floor(.75*n)],iqr=q3-q1;
  const mu=x.reduce((a,b)=>a+b,0)/n, sd=Math.sqrt(x.reduce((a,b)=>a+(b-mu)*(b-mu),0)/n);
  let h=0.9*Math.min(sd||1,(iqr/1.34)||(sd||1))*Math.pow(n,-0.2); if(!(h>0))h=(sd||1)*0.3;
  const G=512,lo=s[0]-h,hi=s[n-1]+h,dens=new Float64Array(G);
  for(let g=0;g<G;g++){const xg=lo+(hi-lo)*g/(G-1);let d=0;
    for(let i=0;i<n;i++){const z=(xg-x[i])/h;d+=Math.exp(-0.5*z*z);} dens[g]=d;}
  let dmax=0;for(let g=0;g<G;g++)if(dens[g]>dmax)dmax=dens[g];
  for(let g=1;g<G-1;g++) if(dens[g]>=dens[g-1]&&dens[g]>=dens[g+1]&&dens[g]>0.15*dmax) return lo+(hi-lo)*g/(G-1);
  return med(x);}
const PHI1=0.8413447;
function oneSided(x){const n=x.length,m=leftMode(x);
  const below=x.filter(v=>v<m);
  const piMirror=Math.max(0,Math.min(1,1-2*below.length/n));
  let sd=below.length?1.4826*med(below.map(v=>m-v)):0; if(!(sd>0))sd=1e-9;
  const c=m+sd, nBelowC=x.filter(v=>v<c).length;
  const piSigma=Math.max(0,Math.min(1,1-(nBelowC/PHI1)/n));
  return {piMirror,piSigma};}
function etaOf(mat,n,J,iters){const ef=R.ensembleFeatures(mat,n,J,{iterations:iters,seed:1}),O=ef.oriented;
  const g=Math.max(0,Math.min(1,(ef.profileInfo-0.01)/0.05));
  const e=new Array(n);for(let i=0;i<n;i++)e[i]=W.b0+W.rr*O.rr[i]+g*(W.longstring*O.longstring[i]+W.person_total*O.person_total[i]);
  return e;}
const mn=a=>a.length?a.reduce((s,v)=>s+v,0)/a.length:NaN;
const out=[["study","J","true_rate","reps","err_mirror","err_1sigma","err_eta","bias_mirror","bias_1sigma","bias_eta"]];
function report(study,J,tr,res,reps){
  const eM=mn(res.map(r=>Math.abs(r.piMirror-tr))),eS=mn(res.map(r=>Math.abs(r.piSigma-tr))),eE=mn(res.map(r=>Math.abs(r.piEta-tr)));
  const bM=mn(res.map(r=>r.piMirror-tr)),bS=mn(res.map(r=>r.piSigma-tr)),bE=mn(res.map(r=>r.piEta-tr));
  out.push([study,J,tr.toFixed(4),reps,eM.toFixed(4),eS.toFixed(4),eE.toFixed(4),bM.toFixed(4),bS.toFixed(4),bE.toFixed(4)]);
  console.log(study.padEnd(7)+"J="+String(J).padStart(3)+" true="+(100*tr).toFixed(0).padStart(3)+"% | |err| mirror "+
    (100*eM).toFixed(1).padStart(5)+"  1sig "+(100*eS).toFixed(1).padStart(5)+"  eta "+(100*eE).toFixed(1).padStart(5)+
    "  | bias "+(100*bM).toFixed(1).padStart(6)+" "+(100*bS).toFixed(1).padStart(6)+" "+(100*bE).toFixed(1).padStart(6));}
console.log("=== STUDY 1 (real data, prevalence swept) ===");
{ const Ms=R.parseCSV(fs.readFileSync("_study_matrix.csv","utf8")),J=Ms.header.length;
  const allR=Ms.rows.map(r=>r.map(v=>{const x=Number(v);return Number.isFinite(x)?x:NaN;}));
  const y=R.parseCSV(fs.readFileSync("_study_labels.csv","utf8")).rows.map(r=>Number(r[0]));
  const cf=y.map((v,i)=>v===0?i:-1).filter(i=>i>=0), cl=y.map((v,i)=>v===1?i:-1).filter(i=>i>=0);
  const REPS=60;
  for(const k of [4,8,12,18,24,33,47,58]){ if(k>cl.length) continue; const tr=k/(cf.length+k),res=[];
    for(let rep=0;rep<REPS;rep++){ const a=cl.slice();
      for(let i=0;i<k;i++){const j=i+Math.floor(rng()*(a.length-i));const t=a[i];a[i]=a[j];a[j]=t;}
      const idx=[...cf,...a.slice(0,k)],n=idx.length,mat=new Float64Array(n*J);
      for(let i=0;i<n;i++)for(let j=0;j<J;j++)mat[i*J+j]=allR[idx[i]][j];
      const e=etaOf(mat,n,J,150),os=oneSided(e),gm=R.gaussMix1D(e);
      res.push({...os,piEta:gateOK(gm)?gm.pi:0}); }
    report("Study1",J,tr,res,REPS); } }
console.log("\n=== STUDY 3 (simulated, six patterns) ===");
{ const REPS=12,N=400;
  for(const [F,ipf] of [[10,6],[30,6]]){ const Jn=F*ipf;
    for(const pct of [0.05,0.10,0.15,0.20,0.25,0.30,0.35,0.40,0.45]){ const res=[];
      for(let rep=0;rep<REPS;rep++){ const {rows,J}=genClean(F,ipf,N); const {out:o}=inject(rows,J,pct);
        const e=etaOf(toMat(o,J),N,J,100),os=oneSided(e),gm=R.gaussMix1D(e);
        res.push({...os,piEta:gateOK(gm)?gm.pi:0}); }
      report("Study3",Jn,pct,res,REPS); } } }
fs.writeFileSync("onesided_eval.csv",out.map(r=>r.join(",")).join("\n"));
console.log("\nwrote onesided_eval.csv");
