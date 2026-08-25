/* sim_shipped.js — controlled simulation using the SHIPPED ReReRe model
 * (logistic triad rr+LongString+Person-Total, two-groups automatic calibration),
 * not the legacy 5-feature random forest. Reports MCC vs questionnaire length under
 * BOTH sensitivity settings (Standard z=2.5, High z=1.5) and the rr ablation. */
const R=require("./site/reread.js");
const W=R.WEIGHTS;
const rng=R._rng(20260709);

function gauss(){let u=0,v=0;while(u===0)u=rng();while(v===0)v=rng();return Math.sqrt(-2*Math.log(u))*Math.cos(2*Math.PI*v);}
// inverse normal CDF (Acklam)
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
// CFA clean data: F factors x ipf items, n respondents, 5-point Likert.
// Realism (so the shipped partner gate + Standard/High cut behave as on real data):
//  - per-item intercepts tau ~ N(0, 0.45) -> item means vary (SD(means)~0.47, matching
//    the real study), which is what makes Person-Total / LongString reliable (profileInfo>0);
//  - ~40% reverse-keyed items (negative loadings), as in real questionnaires (rr aligns signs).
const K=5, THR=[]; for(let c=1;c<K;c++) THR.push(qnorm(c/K));
function genClean(F,ipf,n){
  const J=F*ipf, load=[],fac=[],tau=[];
  for(let j=0;j<J;j++){let l=0.4+0.4*rng(); if(rng()<0.4) l=-l; load.push(l); fac.push(Math.floor(j/ipf)); tau.push(0.45*gauss());}
  const rho=0.3,sq=Math.sqrt(rho),sq1=Math.sqrt(1-rho);
  const rows=[];
  for(let i=0;i<n;i++){
    const common=gauss(),f=[];for(let k=0;k<F;k++)f.push(sq*common+sq1*gauss());
    const row=new Array(J);
    for(let j=0;j<J;j++){
      const y=tau[j]+load[j]*f[fac[j]]+Math.sqrt(1-load[j]*load[j])*gauss();
      let cat=1;for(let c=0;c<K-1;c++)if(y>THR[c])cat=c+2;row[j]=cat;
    }
    rows.push(row);
  }
  return {rows,J};
}
// port of app.js injectCareless (six patterns, corruption levels)
function inject(rows,J,pct){
  const N=rows.length;
  let lo=Infinity,hi=-Infinity;for(let i=0;i<N;i++)for(let j=0;j<J;j++){const v=rows[i][j];if(v<lo)lo=v;if(v>hi)hi=v;}
  const LEVELS=[];for(let v=lo;v<=hi;v++)LEVELS.push(v);
  const highLo=lo+0.65*(hi-lo),highLevels=LEVELS.filter(v=>v>=highLo);
  const ri=n=>Math.floor(rng()*n),pick=a=>a[ri(a.length)];
  const sampleK=(pool,k)=>{const a=pool.slice();for(let i=0;i<k;i++){const j=i+ri(a.length-i);const t=a[i];a[i]=a[j];a[j]=t;}return a.slice(0,k);};
  const allIdx=[...Array(J).keys()],clamp=x=>Math.max(lo,Math.min(hi,x));
  const chunks=(row,idx,k,seq)=>{const nc=k>=4?2+ri(Math.min(5,k)-1):(k>=2?2:1);const asg=[];for(let t=0;t<k;t++)asg.push(t%nc);
    if(!seq)for(let t=k-1;t>0;t--){const j=ri(t+1);const tmp=asg[t];asg[t]=asg[j];asg[j]=tmp;}
    for(let ch=0;ch<nc;ch++){const val=pick(LEVELS);for(let t=0;t<k;t++)if(asg[t]===ch)row[idx[t]]=val;}};
  const patterns=["random","longstring","pure_straight","acquiescent","mixed","fatigue"],levels=[0.5,0.6,0.7,0.8,0.9,1.0];
  const M=Math.round(N*pct),order=sampleK([...Array(N).keys()],N),ids=order.slice(0,M);
  const labels=new Array(N).fill(0),sev=new Array(N).fill(0),out=rows.map(r=>r.slice());
  for(let m=0;m<M;m++){
    const i=ids[m],pat=patterns[m%patterns.length],lvl=pick(levels),k=Math.min(J,Math.max(1,Math.round(J*lvl))),row=out[i];
    if(pat==="random")sampleK(allIdx,k).forEach(p=>row[p]=pick(LEVELS));
    else if(pat==="longstring")chunks(row,sampleK(allIdx,k),k,false);
    else if(pat==="pure_straight"){const v=pick(LEVELS);sampleK(allIdx,k).forEach(p=>row[p]=v);}
    else if(pat==="acquiescent"){const b=pick(highLevels);sampleK(allIdx,k).forEach(p=>row[p]=clamp(b+pick([-1,0,0,0,1])));}
    else if(pat==="mixed"){const kl=Math.max(1,Math.round(k/2)),idx=sampleK(allIdx,k);chunks(row,idx.slice(0,kl),kl,false);idx.slice(kl).forEach(p=>row[p]=pick(LEVELS));}
    else if(pat==="fatigue"){const idx=[];for(let t=J-k;t<J;t++)idx.push(t);chunks(row,idx,k,true);}
    labels[i]=1; sev[i]=lvl;   // severity = corruption fraction (clean=0)
  }
  return {out,labels,sev};
}
function mcc(flag,lab){let tp=0,fp=0,fn=0,tn=0;for(let i=0;i<lab.length;i++){if(lab[i]===1){flag[i]?tp++:fn++;}else{flag[i]?fp++:tn++;}}
  const d=Math.sqrt((tp+fp)*(tp+fn)*(tn+fp)*(tn+fn));return d?(tp*tn-fp*fn)/d:0;}
function auc(s,y){let po=[],ne=[];for(let i=0;i<y.length;i++)(y[i]?po:ne).push(s[i]);let c=0;for(const a of po)for(const b of ne)c+=a>b?1:a===b?0.5:0;return c/(po.length*ne.length);}
// severity weight: heavy for clearly-clean (s=0) and clearly-careless (s=1), light for the ambiguous middle
function uWeight(s){return (2*s-1)*(2*s-1);}
// threshold-free detection ceiling: best (severity-weighted) MCC over thresholds (oracle)
function wmcc(f,lab,w){let tp=0,fp=0,fn=0,tn=0;for(let i=0;i<lab.length;i++){const wi=w[i];if(lab[i]===1){f[i]?tp+=wi:fn+=wi;}else{f[i]?fp+=wi:tn+=wi;}}
  const d=Math.sqrt((tp+fp)*(tp+fn)*(tn+fp)*(tn+fn));return d?(tp*tn-fp*fn)/d:0;}
function oracleMCC(eta,lab,w){const s=[...eta].sort((a,b)=>a-b);let best=0;const step=Math.max(1,Math.floor(s.length/100));
  for(let t=0;t<s.length;t+=step){const thr=s[t];const f=eta.map(v=>v>=thr);const m=w?wmcc(f,lab,w):mcc(f,lab);if(m>best)best=m;}return best;}
function scoreDataset(out,J,labels,sev){
  const n=out.length,mat=new Float64Array(n*J);for(let i=0;i<n;i++)for(let j=0;j<J;j++)mat[i*J+j]=out[i][j];
  const ef=R.ensembleFeatures(mat,n,J,{seed:1,iterations:120}),O=ef.oriented;
  const g=Math.max(0,Math.min(1,(ef.profileInfo-0.01)/0.05));
  const etaF=new Array(n),etaN=new Array(n);
  for(let i=0;i<n;i++){const partners=g*(W.longstring*O.longstring[i]+W.person_total*O.person_total[i]);
    etaF[i]=W.b0+W.rr*O.rr[i]+partners;etaN[i]=W.b0+partners;}   // full (rr+partners) vs no-rr (partners only)
  const w=sev.map(uWeight);
  // rr-alone AUC + z-separation (oriented rr: high=careless)
  const care=[],att=[];for(let i=0;i<n;i++)(labels[i]?care:att).push(O.rr[i]);
  const zsep=(care.reduce((s,v)=>s+v,0)/care.length)-(att.reduce((s,v)=>s+v,0)/att.length);
  return {aucF:auc(etaF,labels),aucN:auc(etaN,labels),aucRR:auc(O.rr,labels),zsep:zsep,omF:oracleMCC(etaF,labels,w),omN:oracleMCC(etaN,labels,w)};
}
const mn=a=>a.reduce((s,v)=>s+v,0)/a.length,sd=a=>{const m=mn(a);return Math.sqrt(mn(a.map(v=>(v-m)**2)));};
function bootCI(deltas,B){const out=[];for(let b=0;b<B;b++){let s=0;for(let i=0;i<deltas.length;i++)s+=deltas[Math.floor(rng()*deltas.length)];out.push(s/deltas.length);}
  out.sort((a,b)=>a-b);return {lo:out[Math.floor(0.025*B)],hi:out[Math.floor(0.975*B)],p:out.filter(v=>v>0).length/B};}

const SIZES=[[5,6],[10,6],[20,6],[30,6],[40,6],[50,6]];  // F x ipf -> 30,60,120,180,240,300
const N=500, RATE=0.20, REPS=20;
const se=a=>sd(a)/Math.sqrt(a.length);
const rows=[["items","auc_rr","auc_rr_se","auc_ens","auc_ens_se","zsep","zsep_se"]];
console.log("items  AUC(rr)   AUC(ens)   z-sep(careless-attentive)");
for(const [F,ipf] of SIZES){
  const J=F*ipf, acc={rr:[],ens:[],zs:[]};
  for(let rep=0;rep<REPS;rep++){
    const {rows:clean}=genClean(F,ipf,N);
    const {out,labels,sev}=inject(clean,J,RATE);
    const r=scoreDataset(out,J,labels,sev);
    acc.rr.push(r.aucRR);acc.ens.push(r.aucF);acc.zs.push(r.zsep);
  }
  rows.push([J,mn(acc.rr).toFixed(4),se(acc.rr).toFixed(4),mn(acc.ens).toFixed(4),se(acc.ens).toFixed(4),
    mn(acc.zs).toFixed(4),se(acc.zs).toFixed(4)]);
  console.log(String(J).padStart(4)+"    "+mn(acc.rr).toFixed(3)+"     "+mn(acc.ens).toFixed(3)+"      "+mn(acc.zs).toFixed(3));
}
require("fs").writeFileSync("sim_rr_scaling.csv",rows.map(r=>r.join(",")).join("\n"));
console.log("\nwrote sim_rr_scaling.csv");
