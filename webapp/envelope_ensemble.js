/* envelope_ensemble.js — operating envelope of the ReReRe ENSEMBLE (triad
 * rr+longstring+person_total, shipped weights) vs the rr index alone.
 * Subsamples J items x n respondents (stratified) from the 84x109 study data,
 * scores each subsample exactly as the deployed tool (features restandardized
 * within the subsample, fixed weights), records AUC(careful vs careless). */
const fs = require("fs");
const R = require("./site/reread.js");
const M = R.parseCSV(fs.readFileSync("_study_matrix.csv", "utf8"));
const L = R.parseCSV(fs.readFileSync("_study_labels.csv", "utf8"));
const rows = M.rows.map(r => r.map(Number)), y = L.rows.map(r => Number(r[0]));
const Jall = M.header.length;
const carefulIdx = y.map((v,i)=>v===0?i:-1).filter(i=>i>=0);
const carelessIdx = y.map((v,i)=>v===1?i:-1).filter(i=>i>=0);
const rng = R._rng(12345);
function sampleK(arr,k){const a=arr.slice();for(let i=0;i<k;i++){const j=i+Math.floor(rng()*(a.length-i));const t=a[i];a[i]=a[j];a[j]=t;}return a.slice(0,k);}
function aucHigh(s,lab){const arr=s.map((v,i)=>({v,y:lab[i]})).sort((a,b)=>a.v-b.v);let r=0,np=0,nn=0;arr.forEach((o,i)=>{if(o.y){r+=i+1;np++;}else nn++;});return (np&&nn)?(r-np*(np+1)/2)/(np*nn):NaN;}
function aucLow(s,lab){return aucHigh(s.map(v=>-v),lab);}

const J_GRID=[10,20,30,45,60,80,109], N_GRID=[20,30,40,50,60,70,76], REPS=30, ITERS=100;
function subrun(J,nTot){
  const half=Math.floor(nTot/2);
  if(half>carefulIdx.length||half>carelessIdx.length) return null;
  const cols=sampleK([...Array(Jall).keys()],J);
  const pick=[...sampleK(carefulIdx,half),...sampleK(carelessIdx,nTot-half)];
  const n=pick.length,sub=new Float64Array(n*J),lab=new Array(n);
  for(let i=0;i<n;i++){lab[i]=y[pick[i]];for(let j=0;j<J;j++)sub[i*J+j]=rows[pick[i]][cols[j]];}
  try{ const res=R.ensemble(sub,n,J,{iterations:ITERS,seed:1+Math.floor(rng()*1e6),skipD2:true});
       return {ens:aucHigh(res.p,lab), rr:aucLow(res.rr,lab), level:res.diagnostic.level}; }
  catch(e){ return {ens:NaN,rr:NaN,level:"err"}; }
}
console.log("ReReRe ENSEMBLE operating envelope — mean AUC over "+REPS+" stratified subsamples");
console.log("dataset: "+rows.length+" resp ("+carefulIdx.length+" careful/"+carelessIdx.length+" careless), "+Jall+" items\n");
const grid={};
for(const nTot of N_GRID) for(const J of J_GRID){
  const ea=[],ra=[]; for(let r=0;r<REPS;r++){const s=subrun(J,nTot); if(s&&Number.isFinite(s.ens)){ea.push(s.ens);ra.push(s.rr);}}
  ea.sort((a,b)=>a-b);
  grid[nTot+"_"+J]={ens:ea.reduce((s,v)=>s+v,0)/ea.length, rr:ra.reduce((s,v)=>s+v,0)/ra.length,
    pGood:ea.filter(a=>a>=0.8).length/ea.length};
}
function tbl(field,label){console.log("\n=== "+label+" ===");let h="n \ J".padEnd(8);for(const J of J_GRID)h+=("J="+J).padStart(8);console.log(h);
  for(const nTot of N_GRID){let ln=("n="+nTot).padEnd(8);for(const J of J_GRID){const v=grid[nTot+"_"+J][field];ln+=(Number.isFinite(v)?v.toFixed(2):"  -").padStart(8);}console.log(ln);}}
tbl("ens","Mean AUC — ENSEMBLE (triad)");
tbl("rr","Mean AUC — rr index alone");
tbl("pGood","Fraction of ensemble subsamples with AUC >= 0.80");
console.log("\n=== Frontier: smallest J reaching mean AUC threshold (ENSEMBLE) ===");
for(const thr of [0.9,0.8,0.7]){console.log("  AUC >= "+thr+":");
  for(const nTot of N_GRID){const js=J_GRID.filter(J=>grid[nTot+"_"+J].ens>=thr);
    console.log("    n="+String(nTot).padEnd(3)+" -> "+(js.length?"J >= "+js[0]:"not reached"));}}
// headline: ensemble vs rr improvement, averaged
let de=0,cnt=0; for(const k in grid){if(Number.isFinite(grid[k].ens)&&Number.isFinite(grid[k].rr)){de+=grid[k].ens-grid[k].rr;cnt++;}}
console.log("\nMean AUC gain ENSEMBLE - rr across all "+cnt+" cells: +"+(de/cnt).toFixed(3));
