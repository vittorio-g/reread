/* real_data_test.js — test all 7 flagging methods (incl. E, 3-component) on the REAL collected
 * study data (n=157: 70 careful + 87 real scrambled careless). Prevalence sweep: keep all careful,
 * sample k real careless to hit target rates; average over resamples. Real careless are extreme
 * (scrambled text) => they DO form a cluster, unlike the heterogeneous simulated careless. */
const R=require("./site/rerere.js"); const W=R.WEIGHTS; const fs=require("fs");
const rng=R._rng(20260714);
const src=fs.readFileSync("sim_v2.js","utf8");
eval(src.slice(src.indexOf("function qnorm"), src.indexOf("function scoreDataset")));  // helpers + metrics + allFlags
const Ms=R.parseCSV(fs.readFileSync("_study_matrix.csv","utf8"));
const J=Ms.header.length, allRows=Ms.rows.map(r=>r.map(v=>{const x=Number(v);return Number.isFinite(x)?x:NaN;}));
const y=R.parseCSV(fs.readFileSync("_study_labels.csv","utf8")).rows.map(r=>Number(r[0]));
const carefulIdx=y.map((v,i)=>v===0?i:-1).filter(i=>i>=0);
const carelessIdx=y.map((v,i)=>v===1?i:-1).filter(i=>i>=0);
console.log("REAL study: "+carefulIdx.length+" careful + "+carelessIdx.length+" careless, "+J+" items");
function sampleK(pool,k){const a=pool.slice();for(let i=0;i<k;i++){const j=i+Math.floor(rng()*(a.length-i));const t=a[i];a[i]=a[j];a[j]=t;}return a.slice(0,k);}
function scoreReal(idx){
  const n=idx.length, mat=new Float64Array(n*J);
  for(let i=0;i<n;i++)for(let j=0;j<J;j++)mat[i*J+j]=allRows[idx[i]][j];
  const ef=R.ensembleFeatures(mat,n,J,{seed:1,iterations:200}), O=ef.oriented;
  const gg=Math.max(0,Math.min(1,(ef.profileInfo-0.01)/0.05));
  const eta=[];for(let i=0;i<n;i++)eta.push(W.b0+W.rr*O.rr[i]+gg*(W.longstring*O.longstring[i]+W.person_total*O.person_total[i]));
  const lab=idx.map(i=>y[i]), F=allFlags(eta), r={oracle:oracleMCC(eta,lab),auc:auc(eta,lab)};
  for(const m of ["std","high","A","B","C","D","E"]){r[m+"_mcc"]=mcc(F[m],lab);r[m+"_rec"]=recall(F[m],lab);r[m+"_prec"]=precision(F[m],lab);}
  return r;
}
const mn=a=>a.reduce((s,v)=>s+v,0)/a.length;
const NCf=carefulIdx.length, REPS=30, M=["std","high","A","B","C","D","E"];
console.log("\ntrueRate  AUC oracle | "+M.map(m=>m.padEnd(4)).join("")+"   (MCC per method)");
const out=[["true_rate","auc","oracle",...M.flatMap(m=>[m+"_mcc",m+"_rec",m+"_prec"])]];
for(const p of [0.05,0.10,0.15,0.20,0.25,0.30,0.40,0.50]){
  const k=Math.round(NCf*p/(1-p)); if(k>carelessIdx.length) continue;
  const acc={auc:[],oracle:[]}; for(const m of M)for(const s of ["_mcc","_rec","_prec"])acc[m+s]=[];
  for(let rep=0;rep<REPS;rep++){const idx=[...carefulIdx,...sampleK(carelessIdx,k)];const r=scoreReal(idx);
    acc.auc.push(r.auc);acc.oracle.push(r.oracle);for(const m of M)for(const s of ["_mcc","_rec","_prec"])acc[m+s].push(r[m+s]);}
  const tr=k/(NCf+k);
  out.push([tr.toFixed(4),mn(acc.auc).toFixed(4),mn(acc.oracle).toFixed(4),...M.flatMap(m=>[mn(acc[m+"_mcc"]).toFixed(4),mn(acc[m+"_rec"]).toFixed(4),mn(acc[m+"_prec"]).toFixed(4)])]);
  console.log(String((100*tr).toFixed(0)).padStart(6)+"%  "+mn(acc.auc).toFixed(2)+"  "+mn(acc.oracle).toFixed(2)+"  | "+M.map(m=>mn(acc[m+"_mcc"]).toFixed(2)+" ").join(""));
}
fs.writeFileSync("real_data_test.csv",out.map(r=>r.join(",")).join("\n"));
console.log("\nwrote real_data_test.csv  (columns: mcc/rec/prec per method)");
