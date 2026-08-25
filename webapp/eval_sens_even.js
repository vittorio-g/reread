const fs=require("fs");const R=require("./site/rerere.js");
const Ms=R.parseCSV(fs.readFileSync("_study_matrix.csv","utf8"));
const J=Ms.header.length, allRows=Ms.rows.map(r=>r.map(Number));
const y=R.parseCSV(fs.readFileSync("_study_labels.csv","utf8")).rows.map(r=>Number(r[0]));
const carefulIdx=y.map((v,i)=>v===0?i:-1).filter(i=>i>=0);   // 38
const carelessIdx=y.map((v,i)=>v===1?i:-1).filter(i=>i>=0);  // 46
const NCf=carefulIdx.length, NCl=carelessIdx.length;
// pick (m careful, k careless) hitting target rate t as closely as possible, max sample size
function pickMK(t){let best=null;
  for(let m=18;m<=NCf;m++)for(let k=1;k<=NCl;k++){const r=k/(k+m),err=Math.abs(r-t);
    if(!best||err<best.err-1e-9||(Math.abs(err-best.err)<1e-9&&(m+k)>best.N))best={m,k,err,N:m+k,r};}
  return best;}
function mat(idx){const n=idx.length,mm=new Float64Array(n*J);for(let i=0;i<n;i++)for(let j=0;j<J;j++)mm[i*J+j]=allRows[idx[i]][j];return {m:mm,n};}
function mcc(f,lab){let tp=0,fp=0,fn=0,tn=0;for(let i=0;i<lab.length;i++){if(lab[i]===1){f[i]?tp++:fn++;}else{f[i]?fp++:tn++;}}
  const d=Math.sqrt((tp+fp)*(tp+fn)*(tn+fp)*(tn+fn));return d?(tp*tn-fp*fn)/d:0;}
const rng=R._rng(4242), REPS=120;
function sampleK(a,k){a=a.slice();for(let i=0;i<k;i++){const j=i+Math.floor(rng()*(a.length-i));const t=a[i];a[i]=a[j];a[j]=t;}return a.slice(0,k);}
const mn=a=>a.reduce((s,v)=>s+v,0)/a.length, sd=a=>{const m=mn(a);return Math.sqrt(mn(a.map(v=>(v-m)**2)));};
const out=[["target","true_rate","m_careful","k_careless","low","medium","high","med_sd"]];
console.log("target achieved  m  k   low  med  high");
for(let t=5;t<=50;t++){
  const p=pickMK(t/100), acc={low:[],medium:[],high:[]};
  for(let rep=0;rep<REPS;rep++){
    const idx=[...sampleK(carefulIdx,p.m),...sampleK(carelessIdx,p.k)], lab=idx.map(i=>y[i]), {m,n}=mat(idx);
    const eta=R.ensemble(m,n,J,{iterations:200,seed:1}).eta;
    for(const L of ["low","medium","high"]) acc[L].push(mcc(R.autoFlag(eta,L).flagged,lab));
  }
  out.push([t,p.r.toFixed(4),p.m,p.k,mn(acc.low).toFixed(4),mn(acc.medium).toFixed(4),mn(acc.high).toFixed(4),sd(acc.medium).toFixed(4)]);
  console.log(String(t).padStart(3)+"%   "+(100*p.r).toFixed(1).padStart(4)+"%  "+String(p.m).padStart(2)+" "+String(p.k).padStart(2)+"   "+mn(acc.low).toFixed(2)+" "+mn(acc.medium).toFixed(2)+" "+mn(acc.high).toFixed(2));
}
fs.writeFileSync("eval_sens_even.csv",out.map(r=>r.join(",")).join("\n"));
console.log("wrote eval_sens_even.csv");
