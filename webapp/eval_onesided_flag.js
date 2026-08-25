/* eval_onesided_flag.js — (c) use the one-sided logic to PLACE THE FLAGGING CUT, not just to
 * estimate the prevalence, then re-score detection quality on Study 1, Study 2 and Study 3.
 *
 *  current rule : gaussMix1D(eta) -> cut = m1 + z*sd1, and NOTHING is flagged if the gate rejects
 *  one-sided    : m = left-most KDE mode; sigma = 1.4826*median(m - x | x<m)   [left side only]
 *                 cut = m + z*sigma, always available (no gate)
 *  Both use the shipped z: Standard 2.5, High 1.5. Same nominal false-positive rate, so the
 *  comparison is like-for-like.
 * Reported: MCC / precision / recall vs the true label, plus the false-positive rate on
 * genuinely CLEAN simulated data (the property the current gate protects).
 * Output: onesided_flag.csv
 */
const R=require("./site/rerere.js"); const W=R.WEIGHTS; const fs=require("fs");
let SEED=1357; const rng=()=>{SEED=(SEED*1103515245+12345)&0x7fffffff;return SEED/0x7fffffff;};
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
const med=a=>{const s=[...a].sort((x,y)=>x-y);const h=s.length>>1;return s.length%2?s[h]:(s[h-1]+s[h])/2;};
function leftMode(x){const n=x.length,s=[...x].sort((a,b)=>a-b);
  const q1=s[Math.floor(.25*n)],q3=s[Math.floor(.75*n)],iqr=q3-q1;
  const mu=x.reduce((a,b)=>a+b,0)/n, sd=Math.sqrt(x.reduce((a,b)=>a+(b-mu)*(b-mu),0)/n);
  let h=0.9*Math.min(sd||1,(iqr/1.34)||(sd||1))*Math.pow(n,-0.2); if(!(h>0))h=(sd||1)*0.3;
  const G=512,lo=s[0]-h,hi=s[n-1]+h,dens=new Float64Array(G);
  for(let g=0;g<G;g++){const xg=lo+(hi-lo)*g/(G-1);let d=0;for(let i=0;i<n;i++){const z=(xg-x[i])/h;d+=Math.exp(-0.5*z*z);}dens[g]=d;}
  let dmax=0;for(let g=0;g<G;g++)if(dens[g]>dmax)dmax=dens[g];
  for(let g=1;g<G-1;g++) if(dens[g]>=dens[g-1]&&dens[g]>=dens[g+1]&&dens[g]>0.15*dmax) return lo+(hi-lo)*g/(G-1);
  return med(x);}
function cutOneSided(e,z){const m=leftMode(e),below=e.filter(v=>v<m);
  let sd=below.length?1.4826*med(below.map(v=>m-v)):0; if(!(sd>0))sd=1e-9; return m+z*sd;}
function gateOK(gm){const p=Math.sqrt((gm.v1+gm.v2)/2)||1;return (gm.bic1-gm.bic2>2)&&gm.pi>0.01&&gm.pi<0.7&&(gm.m2-gm.m1)/p>1.0;}
function mccOf(flag,lab){let tp=0,tn=0,fp=0,fn=0;
  for(let i=0;i<lab.length;i++){if(flag[i]){lab[i]?tp++:fp++;}else{lab[i]?fn++:tn++;}}
  const d=Math.sqrt((tp+fp)*(tp+fn)*(tn+fp)*(tn+fn));
  return {mcc:d?(tp*tn-fp*fn)/d:0, prec:(tp+fp)?tp/(tp+fp):0, rec:(tp+fn)?tp/(tp+fn):0, frac:(tp+fp)/lab.length};}
function etaOf(mat,n,J,iters){const ef=R.ensembleFeatures(mat,n,J,{iterations:iters,seed:1}),O=ef.oriented;
  const g=Math.max(0,Math.min(1,(ef.profileInfo-0.01)/0.05));
  const e=new Array(n);for(let i=0;i<n;i++)e[i]=W.b0+W.rr*O.rr[i]+g*(W.longstring*O.longstring[i]+W.person_total*O.person_total[i]);
  return e;}
function bothRules(e,lab,z){
  const gm=R.gaussMix1D(e);
  const cutCur=gateOK(gm)?(gm.m1+z*Math.sqrt(gm.v1)):Infinity;   // current: nothing flagged if gate rejects
  const cutOS=cutOneSided(e,z);
  return {cur:mccOf(e.map(v=>v>cutCur),lab), os:mccOf(e.map(v=>v>cutOS),lab)};}
const mn=a=>a.length?a.reduce((s,v)=>s+v,0)/a.length:NaN;
const out=[["study","J","true_rate","reps","mcc_cur","mcc_os","prec_cur","prec_os","rec_cur","rec_os"]];
function rep(study,J,tr,res,reps){
  const g=k=>mn(res.map(r=>r[k]));
  out.push([study,J,tr.toFixed(4),reps,g("mc").toFixed(4),g("mo").toFixed(4),g("pc").toFixed(4),g("po").toFixed(4),g("rc").toFixed(4),g("ro").toFixed(4)]);
  console.log(study.padEnd(7)+"J="+String(J).padStart(3)+" true="+(100*tr).toFixed(0).padStart(3)+
    "% | MCC cur "+g("mc").toFixed(3)+" -> OS "+g("mo").toFixed(3)+
    "   | prec "+g("pc").toFixed(2)+"/"+g("po").toFixed(2)+"  rec "+g("rc").toFixed(2)+"/"+g("ro").toFixed(2));}
const pack=o=>({mc:o.cur.mcc,mo:o.os.mcc,pc:o.cur.prec,po:o.os.prec,rc:o.cur.rec,ro:o.os.rec});
const Z=2.5;
console.log("=== STUDY 1 (real, prevalence swept; Standard z=2.5) ===");
{ const Ms=R.parseCSV(fs.readFileSync("_study_matrix.csv","utf8")),J=Ms.header.length;
  const allR=Ms.rows.map(r=>r.map(v=>{const x=Number(v);return Number.isFinite(x)?x:NaN;}));
  const y=R.parseCSV(fs.readFileSync("_study_labels.csv","utf8")).rows.map(r=>Number(r[0]));
  const cf=y.map((v,i)=>v===0?i:-1).filter(i=>i>=0), cl=y.map((v,i)=>v===1?i:-1).filter(i=>i>=0);
  const REPS=40;
  for(const k of [4,8,12,18,24,33,47,58]){ if(k>cl.length) continue; const tr=k/(cf.length+k),res=[];
    for(let r0=0;r0<REPS;r0++){ const a=cl.slice();
      for(let i=0;i<k;i++){const j=i+Math.floor(rng()*(a.length-i));const t=a[i];a[i]=a[j];a[j]=t;}
      const idx=[...cf,...a.slice(0,k)],n=idx.length,mat=new Float64Array(n*J),lab=idx.map(i=>y[i]);
      for(let i=0;i<n;i++)for(let j=0;j<J;j++)mat[i*J+j]=allR[idx[i]][j];
      res.push(pack(bothRules(etaOf(mat,n,J,150),lab,Z))); }
    rep("Study1",J,tr,res,REPS); }
  // full sample
  { const n=allR.length,mat=new Float64Array(n*J);
    for(let i=0;i<n;i++)for(let j=0;j<J;j++)mat[i*J+j]=allR[i][j];
    const o=bothRules(etaOf(mat,n,J,400),y,Z);
    console.log("Study1 FULL n=157 true=55% | MCC cur "+o.cur.mcc.toFixed(3)+" -> OS "+o.os.mcc.toFixed(3)+
      "   | flagged cur "+(100*o.cur.frac).toFixed(1)+"%  OS "+(100*o.os.frac).toFixed(1)+"%");
    out.push(["Study1full",J,"0.5541",1,o.cur.mcc.toFixed(4),o.os.mcc.toFixed(4),o.cur.prec.toFixed(4),o.os.prec.toFixed(4),o.cur.rec.toFixed(4),o.os.rec.toFixed(4)]); } }
console.log("\n=== STUDY 3 (simulated) ===");
{ const REPS=10,N=400;
  for(const [F,ipf] of [[10,6],[30,6]]){ const Jn=F*ipf;
    for(const pct of [0.05,0.10,0.20,0.30,0.40,0.45]){ const res=[];
      for(let r0=0;r0<REPS;r0++){ const {rows,J}=genClean(F,ipf,N); const {out:o,lab}=inject(rows,J,pct);
        res.push(pack(bothRules(etaOf(toMat(o,J),N,J,100),lab,Z))); }
      rep("Study3",Jn,pct,res,REPS); } } }
console.log("\n=== CLEAN DATA (no careless): false-positive rate ===");
{ const REPS=25,N=400;
  for(const [F,ipf] of [[10,6],[30,6]]){ const J=F*ipf; let fc=[],fo=[];
    for(let r0=0;r0<REPS;r0++){ const {rows}=genClean(F,ipf,N); const e=etaOf(toMat(rows,J),N,J,100);
      const gm=R.gaussMix1D(e); const cutCur=gateOK(gm)?(gm.m1+Z*Math.sqrt(gm.v1)):Infinity;
      const cutOS=cutOneSided(e,Z);
      fc.push(e.filter(v=>v>cutCur).length/N); fo.push(e.filter(v=>v>cutOS).length/N); }
    console.log("clean J="+String(J).padStart(3)+" | flagged: current "+(100*mn(fc)).toFixed(2)+"%   one-sided "+(100*mn(fo)).toFixed(2)+"%");
    out.push(["clean",J,"0",REPS,"","","","","",""]); } }
fs.writeFileSync("onesided_flag.csv",out.map(r=>r.join(",")).join("\n"));
console.log("\nwrote onesided_flag.csv");
