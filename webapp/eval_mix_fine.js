const fs=require("fs");const R=require("./site/rerere.js");
const Ms=R.parseCSV(fs.readFileSync("_study_matrix.csv","utf8"));
const J=Ms.header.length, allRows=Ms.rows.map(r=>r.map(Number));
const y=R.parseCSV(fs.readFileSync("_study_labels.csv","utf8")).rows.map(r=>Number(r[0]));
const carefulIdx=y.map((v,i)=>v===0?i:-1).filter(i=>i>=0);
const carelessIdx=y.map((v,i)=>v===1?i:-1).filter(i=>i>=0);
function mat(idx){const n=idx.length,m=new Float64Array(n*J);for(let i=0;i<n;i++)for(let j=0;j<J;j++)m[i*J+j]=allRows[idx[i]][j];return {m,n};}
function auc(sc,lab){const idx=[...sc.keys()].sort((a,b)=>sc[a]-sc[b]);const rk=[];idx.forEach((id,r)=>rk[id]=r+1);
  let n1=0,n0=0,s=0;for(let i=0;i<sc.length;i++){if(lab[i]===1){n1++;s+=rk[i];}else n0++;}return n1&&n0?(s-n1*(n1+1)/2)/(n1*n0):NaN;}
const rng=R._rng(4242);
function sampleK(arr,k){const a=arr.slice();for(let i=0;i<k;i++){const j=i+Math.floor(rng()*(a.length-i));const t=a[i];a[i]=a[j];a[j]=t;}return a.slice(0,k);}
const mn=a=>a.reduce((s,v)=>s+v,0)/a.length, sd=a=>{const m=mn(a);return Math.sqrt(mn(a.map(v=>(v-m)**2)));};
function pct(a,q){const s=[...a].sort((x,y)=>x-y);return s[Math.min(s.length-1,Math.floor(q*s.length))];}
const REPS=300, NC=carefulIdx.length;
const out=[["true_rate","N","k","auc",
  "std_mean","std_sd","std_p10","std_p90",
  "high_mean","high_sd","high_p10","high_p90",
  "pi_mean","pi_sd","reps"]];
console.log("target trueRate  AUC   std(z2.5)   high(z1.5)   pi(mixture)   (300 resamples)");
for(let t=10;t<=30;t+=2){
  const rate=t/100, k=Math.round(NC*rate/(1-rate)), trueRate=k/(NC+k);
  const eStd=[],eHigh=[],ePi=[],aucs=[];
  for(let rep=0;rep<REPS;rep++){
    const idx=[...carefulIdx, ...sampleK(carelessIdx,k)]; const lab=idx.map(i=>y[i]); const {m,n}=mat(idx);
    const ra=R.ensemble(m,n,J,{iterations:200,seed:1,skipD2:true});   // default = standard
    const etaArr=Array.from(ra.eta);
    const afS=R.autoFlag(etaArr,"standard"), afH=R.autoFlag(etaArr,"high");
    eStd.push(afS.estRate); eHigh.push(afH.estRate); ePi.push(afS.pi); aucs.push(auc(ra.p,lab));
  }
  out.push([trueRate.toFixed(4),NC+k,k,mn(aucs).toFixed(4),
    mn(eStd).toFixed(4),sd(eStd).toFixed(4),pct(eStd,0.1).toFixed(4),pct(eStd,0.9).toFixed(4),
    mn(eHigh).toFixed(4),sd(eHigh).toFixed(4),pct(eHigh,0.1).toFixed(4),pct(eHigh,0.9).toFixed(4),
    mn(ePi).toFixed(4),sd(ePi).toFixed(4),REPS]);
  console.log(" "+String(t).padStart(3)+"%   "+(100*trueRate).toFixed(1).padStart(4)+"%  "+mn(aucs).toFixed(2)+
    "   "+(100*mn(eStd)).toFixed(1).padStart(4)+"%±"+(100*sd(eStd)).toFixed(1)+
    "   "+(100*mn(eHigh)).toFixed(1).padStart(4)+"%±"+(100*sd(eHigh)).toFixed(1)+
    "   "+(100*mn(ePi)).toFixed(1).padStart(4)+"%±"+(100*sd(ePi)).toFixed(1));
}
fs.writeFileSync("eval_mix_fine.csv", out.map(r=>r.join(",")).join("\n"));
console.log("\nwrote eval_mix_fine.csv");
