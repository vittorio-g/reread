const fs=require("fs");const R=require("./site/reread.js");
const Ms=R.parseCSV(fs.readFileSync("_study_matrix.csv","utf8"));
const J=Ms.header.length, allRows=Ms.rows.map(r=>r.map(Number));
const y=R.parseCSV(fs.readFileSync("_study_labels.csv","utf8")).rows.map(r=>Number(r[0]));
const carefulIdx=y.map((v,i)=>v===0?i:-1).filter(i=>i>=0);
const carelessIdx=y.map((v,i)=>v===1?i:-1).filter(i=>i>=0);
function mat(idx){const n=idx.length,m=new Float64Array(n*J);for(let i=0;i<n;i++)for(let j=0;j<J;j++)m[i*J+j]=allRows[idx[i]][j];return {m,n};}
function mcc(flag,lab){let tp=0,fp=0,fn=0,tn=0;for(let i=0;i<lab.length;i++){if(lab[i]===1){flag[i]?tp++:fn++;}else{flag[i]?fp++:tn++;}}
  const d=Math.sqrt((tp+fp)*(tp+fn)*(tn+fp)*(tn+fn));return d?(tp*tn-fp*fn)/d:0;}
const rng=R._rng(4242);
function sampleK(a,k){a=a.slice();for(let i=0;i<k;i++){const j=i+Math.floor(rng()*(a.length-i));const t=a[i];a[i]=a[j];a[j]=t;}return a.slice(0,k);}
const mn=a=>a.reduce((s,v)=>s+v,0)/a.length, NC=carefulIdx.length, REPS=30;
const LEVELS=["low","medium","high"];
const out=[["true_rate","low","medium","high"]];
console.log("trueRate   low   medium   high   (auto MCC)");
for(let t=5;t<=50;t+=5){
  const rate=t/100, k=Math.round(NC*rate/(1-rate)), trueRate=k/(NC+k);
  const acc={low:[],medium:[],high:[]};
  for(let rep=0;rep<REPS;rep++){
    const idx=[...carefulIdx,...sampleK(carelessIdx,k)], lab=idx.map(i=>y[i]), {m,n}=mat(idx);
    const eta=R.ensemble(m,n,J,{iterations:250,seed:1}).eta;   // one permutation run
    for(const L of LEVELS) acc[L].push(mcc(R.autoFlag(eta,L).flagged,lab));  // reuse eta
  }
  out.push([trueRate.toFixed(4),mn(acc.low).toFixed(4),mn(acc.medium).toFixed(4),mn(acc.high).toFixed(4)]);
  console.log(" "+(100*trueRate).toFixed(0).padStart(3)+"%   "+mn(acc.low).toFixed(2)+"   "+mn(acc.medium).toFixed(2)+"    "+mn(acc.high).toFixed(2));
}
fs.writeFileSync("eval_sens_sweep.csv",out.map(r=>r.join(",")).join("\n"));
console.log("\nwrote eval_sens_sweep.csv");
