const fs=require("fs");const R=require("./site/rerere.js");
const Ms=R.parseCSV(fs.readFileSync("_study_matrix.csv","utf8"));
const J=Ms.header.length, allRows=Ms.rows.map(r=>r.map(Number));
const y=R.parseCSV(fs.readFileSync("_study_labels.csv","utf8")).rows.map(r=>Number(r[0]));
const carefulIdx=y.map((v,i)=>v===0?i:-1).filter(i=>i>=0);
const carelessIdx=y.map((v,i)=>v===1?i:-1).filter(i=>i>=0);
function mat(idx){const n=idx.length,m=new Float64Array(n*J);for(let i=0;i<n;i++)for(let j=0;j<J;j++)m[i*J+j]=allRows[idx[i]][j];return {m,n};}
function mcc(f,lab){let tp=0,fp=0,fn=0,tn=0;for(let i=0;i<lab.length;i++){if(lab[i]===1){f[i]?tp++:fn++;}else{f[i]?fp++:tn++;}}
  const d=Math.sqrt((tp+fp)*(tp+fn)*(tn+fp)*(tn+fn));return d?(tp*tn-fp*fn)/d:0;}
const rng=R._rng(4242), NC=carefulIdx.length, REPS=200;
function sampleK(a,k){a=a.slice();for(let i=0;i<k;i++){const j=i+Math.floor(rng()*(a.length-i));const t=a[i];a[i]=a[j];a[j]=t;}return a.slice(0,k);}
const mn=a=>a.reduce((s,v)=>s+v,0)/a.length, sd=a=>{const m=mn(a);return Math.sqrt(mn(a.map(v=>(v-m)**2)));};
const out=[["true_rate","k","low","medium","high","med_sd","sep","twoComp_rate","pi"]];
console.log("k   trueRate  low  med  high | medSD  sep  2comp%");
for(let k=2;k<=38;k++){
  const trueRate=k/(NC+k), acc={low:[],medium:[],high:[]}, sep=[],two=[],pi=[];
  for(let rep=0;rep<REPS;rep++){
    const idx=[...carefulIdx,...sampleK(carelessIdx,k)], lab=idx.map(i=>y[i]), {m,n}=mat(idx);
    const eta=R.ensemble(m,n,J,{iterations:200,seed:1}).eta;
    const afM=R.autoFlag(eta,"medium"); sep.push(afM.separation); two.push(afM.twoComp?1:0); pi.push(afM.pi);
    acc.medium.push(mcc(afM.flagged,lab));
    acc.low.push(mcc(R.autoFlag(eta,"low").flagged,lab));
    acc.high.push(mcc(R.autoFlag(eta,"high").flagged,lab));
  }
  out.push([trueRate.toFixed(4),k,mn(acc.low).toFixed(4),mn(acc.medium).toFixed(4),mn(acc.high).toFixed(4),
            sd(acc.medium).toFixed(4),mn(sep).toFixed(3),(mn(two)).toFixed(3),mn(pi).toFixed(3)]);
  console.log(String(k).padStart(2)+"  "+(100*trueRate).toFixed(1).padStart(4)+"%  "+mn(acc.low).toFixed(2)+" "+mn(acc.medium).toFixed(2)+" "+mn(acc.high).toFixed(2)+
    " | "+sd(acc.medium).toFixed(3)+"  "+mn(sep).toFixed(2)+"  "+(100*mn(two)).toFixed(0)+"%");
}
fs.writeFileSync("eval_sens_sweep_fine.csv",out.map(r=>r.join(",")).join("\n"));
console.log("wrote eval_sens_sweep_fine.csv");
