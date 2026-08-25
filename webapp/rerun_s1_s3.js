/* rerun_s1_s3.js — full re-evaluation of Study 1 (real) and Study 3 (simulated) under the
 * updated engine (HWHM sigma, B1 global mode, B3 tail-OR gate, piEstimable, heavyContamination).
 * Reports at each true rate: auto-flag MCC, pi-hat, estimable count, heavy-warning count. */
const R = require("./site/reread.js"); const W = R.WEIGHTS; const fs = require("fs");
let SEED = 606; const rng = () => { SEED = (SEED*1103515245+12345)&0x7fffffff; return SEED/0x7fffffff; };
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
function inject(rows,J,pct){const N=rows.length,out=rows.map(r=>r.slice()),lab=new Array(N).fill(0);
  const M=Math.round(N*pct),ord=[...Array(N).keys()];
  for(let i=0;i<M;i++){const j=i+Math.floor(rng()*(N-i));const t=ord[i];ord[i]=ord[j];ord[j]=t;}
  const pats=["random","longstring","pure_straight","acquiescent","mixed","fatigue"];
  for(let m=0;m<M;m++){const i=ord[m],pat=pats[m%6],lvl=0.5+0.5*rng(),k=Math.max(1,Math.round(J*lvl)),row=out[i];
    const idx=[...Array(J).keys()];for(let t=0;t<k;t++){const j=t+Math.floor(rng()*(J-t));const tm=idx[t];idx[t]=idx[j];idx[j]=tm;}
    const sel=idx.slice(0,k);
    if(pat==="pure_straight"||pat==="longstring"){const v=1+Math.floor(rng()*K5);sel.forEach(p=>row[p]=v);}
    else if(pat==="acquiescent"){sel.forEach(p=>row[p]=Math.max(1,Math.min(K5,4+Math.round(rng()))));}
    else if(pat==="fatigue"){for(let t=J-k;t<J;t++)row[t]=1+Math.floor(rng()*K5);}
    else sel.forEach(p=>row[p]=1+Math.floor(rng()*K5));
    lab[i]=1;}
  return {out,lab};}
const toMat=(rows,J)=>{const n=rows.length,m=new Float64Array(n*J);for(let i=0;i<n;i++)for(let j=0;j<J;j++)m[i*J+j]=rows[i][j];return m;};
function etaOf(mat,n,J,it){const ef=R.ensembleFeatures(mat,n,J,{iterations:it,seed:1}),O=ef.oriented;
  const g=Math.max(0,Math.min(1,(ef.profileInfo-0.01)/0.05));
  const e=new Array(n);for(let i=0;i<n;i++)e[i]=W.b0+W.rr*O.rr[i]+g*(W.longstring*O.longstring[i]+W.person_total*O.person_total[i]);
  return e;}
function mcc(flag,y){let tp=0,tn=0,fp=0,fn=0;
  for(let i=0;i<y.length;i++){if(flag[i]){y[i]?tp++:fp++;}else{y[i]?fn++:tn++;}}
  const d=Math.sqrt((tp+fp)*(tp+fn)*(tn+fp)*(tn+fn));return d?(tp*tn-fp*fn)/d:0;}
const mn=a=>a.reduce((s,v)=>s+v,0)/a.length;
function block(tag,res,tr){
  console.log(tag+"  true="+(100*tr).toFixed(0).padStart(3)+"% | MCC="+mn(res.map(r=>r.m)).toFixed(3)+
    "  piHat="+(100*mn(res.map(r=>r.pi))).toFixed(1).padStart(5)+"%  estimable "+res.filter(r=>r.e).length+"/"+res.length+
    "  heavy "+res.filter(r=>r.h).length+"/"+res.length);}
console.log("=== STUDY 1 (70 attentive + k real careless, updated engine, 10 reps) ===");
{ const Ms=R.parseCSV(fs.readFileSync("_study_matrix.csv","utf8")),J=Ms.header.length;
  const allR=Ms.rows.map(r=>r.map(v=>{const x=Number(v);return Number.isFinite(x)?x:NaN;}));
  const y=R.parseCSV(fs.readFileSync("_study_labels.csv","utf8")).rows.map(r=>Number(r[0]));
  const cf=y.map((v,i)=>v===0?i:-1).filter(i=>i>=0), cl=y.map((v,i)=>v===1?i:-1).filter(i=>i>=0);
  for(const k of [4,8,12,18,24,33,47,58,87]){ if(k>cl.length) continue; const tr=k/(cf.length+k),res=[];
    const REPS=k===87?1:10;
    for(let rep=0;rep<REPS;rep++){ const a=cl.slice();
      for(let i=0;i<k;i++){const j=i+Math.floor(rng()*(a.length-i));const t=a[i];a[i]=a[j];a[j]=t;}
      const idx=[...cf,...a.slice(0,k)],n=idx.length,mat=new Float64Array(n*J),lab=idx.map(i=>y[i]);
      for(let i=0;i<n;i++)for(let j=0;j<J;j++)mat[i*J+j]=allR[idx[i]][j];
      const eta=etaOf(mat,n,J,150),af=R.autoFlag(eta,"standard");
      res.push({m:mcc(af.flagged,lab),pi:af.pi,e:af.piEstimable,h:af.heavyContamination});}
    block("S1",res,tr);}}
console.log("\n=== STUDY 3 (simulated, 6 patterns, J=120 n=500, updated engine, 10 reps) ===");
{ const REPS=10,N=500,F=20,ipf=6;
  for(const pct of [0.05,0.10,0.20,0.30,0.40]){ const res=[];
    for(let rep=0;rep<REPS;rep++){ const {rows,J}=genClean(F,ipf,N); const {out:o,lab}=inject(rows,J,pct);
      const eta=etaOf(toMat(o,J),N,J,100),af=R.autoFlag(eta,"standard");
      res.push({m:mcc(af.flagged,lab),pi:af.pi,e:af.piEstimable,h:af.heavyContamination});}
    block("S3",res,pct);}}
