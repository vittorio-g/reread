const fs=require("fs");const R=require("./site/reread.js");
const Ms=R.parseCSV(fs.readFileSync("_study_matrix.csv","utf8"));
const J=Ms.header.length, allRows=Ms.rows.map(r=>r.map(Number));
const y=R.parseCSV(fs.readFileSync("_study_labels.csv","utf8")).rows.map(r=>Number(r[0]));
const carefulIdx=y.map((v,i)=>v===0?i:-1).filter(i=>i>=0);
const carelessIdx=y.map((v,i)=>v===1?i:-1).filter(i=>i>=0);
function mat(idx){const n=idx.length,m=new Float64Array(n*J);for(let i=0;i<n;i++)for(let j=0;j<J;j++)m[i*J+j]=allRows[idx[i]][j];return {m,n};}
function auc(sc,lab){const idx=[...sc.keys()].sort((a,b)=>sc[a]-sc[b]);const rk=[];idx.forEach((id,r)=>rk[id]=r+1);
  let n1=0,n0=0,s=0;for(let i=0;i<sc.length;i++){if(lab[i]===1){n1++;s+=rk[i];}else n0++;}return n1&&n0?(s-n1*(n1+1)/2)/(n1*n0):NaN;}
function metrics(flag,lab){let tp=0,fp=0,fn=0,tn=0;for(let i=0;i<lab.length;i++){if(lab[i]===1){flag[i]?tp++:fn++;}else{flag[i]?fp++:tn++;}}
  const prec=tp+fp?tp/(tp+fp):(tp?1:0),rec=tp+fn?tp/(tp+fn):0,spec=tn+fp?tn/(tn+fp):1;
  const d=Math.sqrt((tp+fp)*(tp+fn)*(tn+fp)*(tn+fn)),mcc=d?(tp*tn-fp*fn)/d:0;return {prec,rec,spec,mcc};}
const rng=R._rng(4242);
function sampleK(arr,k){const a=arr.slice();for(let i=0;i<k;i++){const j=i+Math.floor(rng()*(a.length-i));const t=a[i];a[i]=a[j];a[j]=t;}return a.slice(0,k);}
const mn=a=>a.reduce((s,v)=>s+v,0)/a.length, sd=a=>{const m=mn(a);return Math.sqrt(mn(a.map(v=>(v-m)**2)));};
const REPS=30, NC=carefulIdx.length, out=[["true_rate","N","k","auc","detect","est_mean","est_sd","auto_mcc","auto_rec","auto_prec","high_mcc","high_rec","high_prec","oracle_mcc","oracle_rec"]];
console.log("target trueRate AUC detect%  est(mean±sd)   autoMCC  oracleMCC");
for(let t=5;t<=50;t+=5){
  const rate=t/100, k=Math.round(NC*rate/(1-rate)), trueRate=k/(NC+k);
  const A={auc:[],est:[],rec:[],prec:[],mcc:[],det:[]}, H={rec:[],prec:[],mcc:[]}, O={rec:[],mcc:[]};
  for(let rep=0;rep<REPS;rep++){
    const idx=[...carefulIdx, ...sampleK(carelessIdx,k)]; const lab=idx.map(i=>y[i]); const {m,n}=mat(idx);
    const ra=R.ensemble(m,n,J,{iterations:250,seed:1,skipD2:true});     // default = Standard
    A.auc.push(auc(ra.p,lab)); A.est.push(ra.estimatedRate); A.det.push(ra.twoComponent?1:0);
    const ma=metrics(ra.flagged,lab); A.rec.push(ma.rec); A.prec.push(ma.prec); A.mcc.push(ma.mcc);
    const afH=R.autoFlag(Array.from(ra.eta),"high"); const mh=metrics(afH.flagged,lab);   // High setting
    H.rec.push(mh.rec); H.prec.push(mh.prec); H.mcc.push(mh.mcc);
    const ro=R.ensemble(m,n,J,{iterations:250,seed:1,prevalence:trueRate,skipD2:true}); const mo=metrics(ro.flagged,lab);
    O.rec.push(mo.rec); O.mcc.push(mo.mcc);
  }
  out.push([trueRate.toFixed(4),NC+k,k,mn(A.auc).toFixed(4),mn(A.det).toFixed(3),mn(A.est).toFixed(4),sd(A.est).toFixed(4),
            mn(A.mcc).toFixed(4),mn(A.rec).toFixed(4),mn(A.prec).toFixed(4),
            mn(H.mcc).toFixed(4),mn(H.rec).toFixed(4),mn(H.prec).toFixed(4),
            mn(O.mcc).toFixed(4),mn(O.rec).toFixed(4)]);
  console.log(" "+String(t).padStart(3)+"%   "+(100*trueRate).toFixed(0).padStart(3)+"%  "+mn(A.auc).toFixed(2)+"  "+
    (100*mn(A.det)).toFixed(0).padStart(3)+"%   "+(100*mn(A.est)).toFixed(0).padStart(2)+"%±"+(100*sd(A.est)).toFixed(0)+
    "     "+mn(A.mcc).toFixed(2)+"     "+mn(O.mcc).toFixed(2));
}
fs.writeFileSync("eval_mix.csv", out.map(r=>r.join(",")).join("\n"));
console.log("\nwrote eval_mix.csv");
