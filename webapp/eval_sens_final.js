const fs=require("fs");const R=require("./site/reread.js");
const Ms=R.parseCSV(fs.readFileSync("_study_matrix.csv","utf8"));
const J=Ms.header.length, allRows=Ms.rows.map(r=>r.map(Number));
const y=R.parseCSV(fs.readFileSync("_study_labels.csv","utf8")).rows.map(r=>Number(r[0]));
const carefulIdx=y.map((v,i)=>v===0?i:-1).filter(i=>i>=0);   // careful base (dynamic, n=67 now)
const carelessIdx=y.map((v,i)=>v===1?i:-1).filter(i=>i>=0);
const NCf=carefulIdx.length;
function mat(idx){const n=idx.length,mm=new Float64Array(n*J);for(let i=0;i<n;i++)for(let j=0;j<J;j++)mm[i*J+j]=allRows[idx[i]][j];return {m:mm,n};}
function mcc(f,lab){let tp=0,fp=0,fn=0,tn=0;for(let i=0;i<lab.length;i++){if(lab[i]===1){f[i]?tp++:fn++;}else{f[i]?fp++:tn++;}}
  const d=Math.sqrt((tp+fp)*(tp+fn)*(tn+fp)*(tn+fn));return d?(tp*tn-fp*fn)/d:0;}
const rng=R._rng(4242), REPS=120;
function sampleK(a,k){a=a.slice();for(let i=0;i<k;i++){const j=i+Math.floor(rng()*(a.length-i));const t=a[i];a[i]=a[j];a[j]=t;}return a.slice(0,k);}
const mn=a=>a.reduce((s,v)=>s+v,0)/a.length;
// independent resamples for one level: all careful + k careless, each rep its own careless draw
function scoreLevel(k,L){const acc=[];
  for(let rep=0;rep<REPS;rep++){
    const idx=[...carefulIdx,...sampleK(carelessIdx,k)], lab=idx.map(i=>y[i]), {m,n}=mat(idx);
    const eta=R.ensemble(m,n,J,{iterations:150,seed:1,skipD2:true}).eta;
    acc.push(mcc(R.autoFlag(eta,L).flagged,lab));
  } return mn(acc);}
const out=[["true_rate","k","low","medium","high"]];
console.log("k  trueRate  low  med  high  (independent per level, careful base)");
for(let k=2;k<=62;k++){
  const tr=k/(NCf+k);
  const low=scoreLevel(k,"low"), medium=scoreLevel(k,"medium"), high=scoreLevel(k,"high");
  out.push([tr.toFixed(4),k,low.toFixed(4),medium.toFixed(4),high.toFixed(4)]);
  console.log(String(k).padStart(2)+"  "+(100*tr).toFixed(1).padStart(4)+"%  "+low.toFixed(2)+" "+medium.toFixed(2)+" "+high.toFixed(2));
}
fs.writeFileSync("eval_sens_final.csv",out.map(r=>r.join(",")).join("\n"));
console.log("wrote eval_sens_final.csv");
