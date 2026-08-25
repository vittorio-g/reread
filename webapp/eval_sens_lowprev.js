const fs=require("fs");const R=require("./site/reread.js");
const Ms=R.parseCSV(fs.readFileSync("_study_matrix.csv","utf8"));
const J=Ms.header.length, allRows=Ms.rows.map(r=>r.map(Number));
const y=R.parseCSV(fs.readFileSync("_study_labels.csv","utf8")).rows.map(r=>Number(r[0]));
const carefulIdx=y.map((v,i)=>v===0?i:-1).filter(i=>i>=0);
const carelessIdx=y.map((v,i)=>v===1?i:-1).filter(i=>i>=0);
function mat(idx){const n=idx.length,m=new Float64Array(n*J);for(let i=0;i<n;i++)for(let j=0;j<J;j++)m[i*J+j]=allRows[idx[i]][j];return {m,n};}
const rng=R._rng(4242), NC=carefulIdx.length, REPS=40;
function sampleK(a,k){a=a.slice();for(let i=0;i<k;i++){const j=i+Math.floor(rng()*(a.length-i));const t=a[i];a[i]=a[j];a[j]=t;}return a.slice(0,k);}
const mn=a=>a.reduce((s,v)=>s+v,0)/a.length;
console.log("rate  sens     detect%  flagged  TP    FP    recall  precision");
for(const t of [5,10,15,20]){
  const rate=t/100, k=Math.round(NC*rate/(1-rate)), trueRate=k/(NC+k), nTrue=k;
  const S={low:[],medium:[],high:[]};
  for(let rep=0;rep<REPS;rep++){
    const idx=[...carefulIdx,...sampleK(carelessIdx,k)], lab=idx.map(i=>y[i]), {m,n}=mat(idx);
    const eta=R.ensemble(m,n,J,{iterations:250,seed:1}).eta;
    for(const L of ["low","medium","high"]){
      const af=R.autoFlag(eta,L), f=af.flagged;
      let tp=0,fp=0;for(let i=0;i<n;i++){if(f[i]){lab[i]?tp++:fp++;}}
      S[L].push({flag:tp+fp,tp,fp,rec:tp/nTrue,prec:(tp+fp)?tp/(tp+fp):0,det:af.twoComp?1:0});
    }
  }
  for(const L of ["low","medium","high"]){
    const a=S[L];
    console.log(" "+String(t).padStart(2)+"%  "+L.padEnd(7)+"  "+(100*mn(a.map(x=>x.det))).toFixed(0).padStart(4)+"%    "+
      mn(a.map(x=>x.flag)).toFixed(1).padStart(4)+"   "+mn(a.map(x=>x.tp)).toFixed(1).padStart(4)+"  "+
      mn(a.map(x=>x.fp)).toFixed(1).padStart(4)+"   "+mn(a.map(x=>x.rec)).toFixed(2)+"    "+mn(a.map(x=>x.prec)).toFixed(2));
  }
  console.log("      (true careless k="+k+" of N="+(NC+k)+")");
}
