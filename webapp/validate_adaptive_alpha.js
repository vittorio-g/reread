/* validate_adaptive_alpha.js — does the n-adaptive tail threshold (1e-3 below n=400) recover
 * availability at low prevalence without paying on clean data?  Two checks:
 *   A. clean simulated data at n=200/350 (the newly-relaxed regime): gate fire rate + flagged share
 *   B. Study 1 at 5/10/15% true rate: availability and MCC (was 5/10 estimable at 5%). */
const R = require("./site/reread.js"); const W = R.WEIGHTS; const fs = require("fs");
let SEED = 2468; const rng = () => { SEED = (SEED*1103515245+12345)&0x7fffffff; return SEED/0x7fffffff; };
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
function mcc(flag,y){let tp=0,tn=0,fp=0,fn=0;
  for(let i=0;i<y.length;i++){if(flag[i]){y[i]?tp++:fp++;}else{y[i]?fn++:tn++;}}
  const d=Math.sqrt((tp+fp)*(tp+fn)*(tn+fp)*(tn+fn));return d?(tp*tn-fp*fn)/d:0;}
console.log("=== A. CLEAN data in the relaxed regime (alpha 1e-3) ===");
for(const [N,REPS] of [[200,25],[350,25]]){
  let fire=0,frac=0,worst=0;
  for(let r=0;r<REPS;r++){const {rows,J}=genClean(10,6,N);
    const af=R.autoFlag(etaOf(toMat(rows,J),N,J,100),"standard");
    if(af.twoComp){fire++;frac+=af.estRate;worst=Math.max(worst,af.estRate);}}
  console.log("clean J=60 n="+N+" | fires "+(100*fire/REPS).toFixed(0)+"%  mean flagged (fired) "+
    (fire?(100*frac/fire).toFixed(1):"0.0")+"%  worst "+(100*worst).toFixed(1)+"%");}
console.log("\n=== B. Study 1 low prevalence (availability + MCC) ===");
{ const Ms=R.parseCSV(fs.readFileSync("_study_matrix.csv","utf8")),J=Ms.header.length;
  const allR=Ms.rows.map(r=>r.map(v=>{const x=Number(v);return Number.isFinite(x)?x:NaN;}));
  const y=R.parseCSV(fs.readFileSync("_study_labels.csv","utf8")).rows.map(r=>Number(r[0]));
  const cf=y.map((v,i)=>v===0?i:-1).filter(i=>i>=0), cl=y.map((v,i)=>v===1?i:-1).filter(i=>i>=0);
  for(const k of [4,8,12]){ const tr=k/(cf.length+k); let est=0,ms=[],pis=[];
    for(let rep=0;rep<10;rep++){ const a=cl.slice();
      for(let i=0;i<k;i++){const j=i+Math.floor(rng()*(a.length-i));const t=a[i];a[i]=a[j];a[j]=t;}
      const idx=[...cf,...a.slice(0,k)],n=idx.length,mat=new Float64Array(n*J),lab=idx.map(i=>y[i]);
      for(let i=0;i<n;i++)for(let j=0;j<J;j++)mat[i*J+j]=allR[idx[i]][j];
      const af=R.autoFlag(etaOf(mat,n,J,150),"standard");
      if(af.piEstimable)est++; ms.push(mcc(af.flagged,lab)); pis.push(af.pi);}
    const mn=a=>a.reduce((s,v)=>s+v,0)/a.length;
    console.log("S1 true="+(100*tr).toFixed(0)+"% | estimable "+est+"/10  MCC="+mn(ms).toFixed(3)+
      "  piHat="+(100*mn(pis)).toFixed(1)+"%");}}
