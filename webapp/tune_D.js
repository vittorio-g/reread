/* tune_D.js — find the best SINGLE cut rule (to replace the Standard/High double setting).
 * Evaluates fixed-z cuts and adaptive rules on BOTH simulated (weighted-MCC, realistic careless
 * + distractors) and the real study data (MCC, extreme careless). A good single default must do
 * well on both. */
const R=require("./site/reread.js"); const W=R.WEIGHTS; const fs=require("fs"); const rng=R._rng(7788);
const src=fs.readFileSync("sim_v2.js","utf8");
eval(src.slice(src.indexOf("function gauss"), src.indexOf("function scoreDataset")));  // genClean, inject2, uWeight, wmcc, mcc, gateOK, ...
function getEta(out,J){const n=out.length,mat=new Float64Array(n*J);for(let i=0;i<n;i++)for(let j=0;j<J;j++)mat[i*J+j]=out[i][j];
  const ef=R.ensembleFeatures(mat,n,J,{seed:1,iterations:150}),O=ef.oriented,g=Math.max(0,Math.min(1,(ef.profileInfo-0.01)/0.05));
  const eta=new Array(n);for(let i=0;i<n;i++)eta[i]=W.b0+W.rr*O.rr[i]+g*(W.longstring*O.longstring[i]+W.person_total*O.person_total[i]);return eta;}
function cutFlags(eta,rule){const gm=R.gaussMix1D(eta);if(!gateOK(gm))return eta.map(_=>false);const s1=Math.sqrt(gm.v1);
  const pooled=Math.sqrt((gm.v1+gm.v2)/2)||1,sep=(gm.m2-gm.m1)/pooled;let z;
  if(rule.t==="fix")z=rule.z; else if(rule.t==="sep")z=Math.max(rule.lo,Math.min(rule.hi,rule.a-rule.b*sep));
  else z=Math.max(rule.lo,Math.min(rule.hi,rule.a-rule.b*gm.pi));   // prevalence-adaptive
  const cut=gm.m1+z*s1;return eta.map(v=>v>cut);}
const RULES={
  "z1.2":{t:"fix",z:1.2},"z1.5(High)":{t:"fix",z:1.5},"z1.8":{t:"fix",z:1.8},"z2.0":{t:"fix",z:2.0},"z2.5(Std)":{t:"fix",z:2.5},
  "Dcur(sep)":{t:"sep",a:3.0,b:0.45,lo:1.2,hi:2.5},
  "Dsep-lo":{t:"sep",a:2.6,b:0.35,lo:1.0,hi:2.3},
  "Dprev":{t:"prev",a:2.4,b:5.0,lo:1.1,hi:2.5},        // lower z when more careless estimated
  "Dprev2":{t:"prev",a:2.2,b:3.5,lo:1.2,hi:2.4}
};
const RN=Object.keys(RULES);
const mn=a=>a.reduce((s,v)=>s+v,0)/a.length;
// ---- SIM: weighted-MCC, realistic careless + distractors, a few prevalences ----
console.log("=== SIM (weighted-MCC, 120 items, +15% distractors, 25 reps) ===");
console.log("prev  | "+RN.map(r=>r.padEnd(11)).join(""));
const simAgg={}; for(const r of RN)simAgg[r]=[];
for(const p of [0.10,0.15,0.20,0.25]){
  const acc={};for(const r of RN)acc[r]=[];
  for(let rep=0;rep<25;rep++){const {rows}=genClean(20,6,600);const {out,labels,dist,sev}=inject2(rows,120,p,0.15);
    const eta=getEta(out,120),w=sev.map(uWeight);
    for(const r of RN)acc[r].push(wmcc(cutFlags(eta,RULES[r]),labels,w));}
  console.log(String(p.toFixed(2)).padStart(4)+"  | "+RN.map(r=>mn(acc[r]).toFixed(3).padEnd(11)).join(""));
  for(const r of RN)simAgg[r].push(mn(acc[r]));
}
// ---- REAL study data: MCC, prevalence sweep ----
console.log("\n=== REAL study (MCC, prevalence sweep, 30 resamples) ===");
const Ms=R.parseCSV(fs.readFileSync("_study_matrix.csv","utf8")),Jr=Ms.header.length;
const allR=Ms.rows.map(r=>r.map(v=>{const x=Number(v);return Number.isFinite(x)?x:NaN;}));
const yR=R.parseCSV(fs.readFileSync("_study_labels.csv","utf8")).rows.map(r=>Number(r[0]));
const cf=yR.map((v,i)=>v===0?i:-1).filter(i=>i>=0),cl=yR.map((v,i)=>v===1?i:-1).filter(i=>i>=0);
function sampleK(pool,k){const a=pool.slice();for(let i=0;i<k;i++){const j=i+Math.floor(rng()*(a.length-i));const t=a[i];a[i]=a[j];a[j]=t;}return a.slice(0,k);}
function etaReal(idx){const n=idx.length,mat=new Float64Array(n*Jr);for(let i=0;i<n;i++)for(let j=0;j<Jr;j++)mat[i*Jr+j]=allR[idx[i]][j];
  const ef=R.ensembleFeatures(mat,n,Jr,{seed:1,iterations:200}),O=ef.oriented,g=Math.max(0,Math.min(1,(ef.profileInfo-0.01)/0.05));
  const eta=new Array(n);for(let i=0;i<n;i++)eta[i]=W.b0+W.rr*O.rr[i]+g*(W.longstring*O.longstring[i]+W.person_total*O.person_total[i]);return eta;}
console.log("prev  | "+RN.map(r=>r.padEnd(11)).join(""));
const realAgg={};for(const r of RN)realAgg[r]=[];
for(const p of [0.10,0.15,0.20,0.25,0.30,0.40]){
  const k=Math.round(cf.length*p/(1-p));if(k>cl.length)continue;
  const acc={};for(const r of RN)acc[r]=[];
  for(let rep=0;rep<30;rep++){const idx=[...cf,...sampleK(cl,k)],lab=idx.map(i=>yR[i]),eta=etaReal(idx);
    for(const r of RN)acc[r].push(mcc(cutFlags(eta,RULES[r]),lab));}
  console.log(String(p.toFixed(2)).padStart(4)+"  | "+RN.map(r=>mn(acc[r]).toFixed(3).padEnd(11)).join(""));
  for(const r of RN)realAgg[r].push(mn(acc[r]));
}
console.log("\n=== SUMMARY: mean over conditions ===");
console.log("rule          sim(wMCC)  real(MCC)  combined");
const rows=[["rule","sim_wmcc","real_mcc","combined"]];
for(const r of RN){const s=mn(simAgg[r]),re=mn(realAgg[r]),c=(s+re)/2;
  console.log(r.padEnd(13)+"  "+s.toFixed(3)+"      "+re.toFixed(3)+"      "+c.toFixed(3));
  rows.push([r,s.toFixed(4),re.toFixed(4),c.toFixed(4)]);}
fs.writeFileSync("tune_D.csv",rows.map(r=>r.join(",")).join("\n"));
console.log("\nwrote tune_D.csv");
