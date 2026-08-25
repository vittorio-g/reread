/* test_gate_tailtest.js — evaluate the excess-tail gate against the shipped mixture gate.
 *  PART A: the 6 external datasets with bogus/attention-check ground truth - does the gate open
 *          where careless respondents demonstrably exist, and what MCC follows?
 *  PART B: genuinely CLEAN simulated data - false-alarm rate (the property the gate protects).
 */
const R = require("./site/rerere.js"); const W = R.WEIGHTS; const fs = require("fs");
const { tailGate } = require("./gate_tailtest.js");
const ALPHA = Number(process.argv[2] || 0.01);
let SEED = 987; const rng = () => { SEED = (SEED*1103515245+12345)&0x7fffffff; return SEED/0x7fffffff; };
function gauss(){let u=0,v=0;while(!u)u=rng();while(!v)v=rng();return Math.sqrt(-2*Math.log(u))*Math.cos(2*Math.PI*v);}
const K5=5, THR=[-0.84,-0.25,0.25,0.84];
function genClean(F,ipf,n){const J=F*ipf;const load=[],fac=[],tau=[];
  for(let j=0;j<J;j++){let l=0.4+0.4*rng();if(rng()<0.4)l=-l;load.push(l);fac.push(Math.floor(j/ipf));tau.push(0.45*gauss());}
  const rho=0.3,sq=Math.sqrt(rho),sq1=Math.sqrt(1-rho);const rows=[];
  for(let i=0;i<n;i++){const c=gauss(),f=[];for(let k=0;k<F;k++)f.push(sq*c+sq1*gauss());
    const row=new Array(J);
    for(let j=0;j<J;j++){let y=tau[j]+load[j]*f[fac[j]];y+=Math.sqrt(Math.max(0.05,1-load[j]*load[j]))*gauss();
      let cat=1;for(let q=0;q<K5-1;q++)if(y>THR[q])cat=q+2;row[j]=cat;}
    rows.push(row);}
  return {rows,J};}
const toMat=(rows,J)=>{const n=rows.length,m=new Float64Array(n*J);for(let i=0;i<n;i++)for(let j=0;j<J;j++)m[i*J+j]=rows[i][j];return m;};
function etaOf(mat,n,J,iters){const ef=R.ensembleFeatures(mat,n,J,{iterations:iters,seed:1}),O=ef.oriented;
  const g=Math.max(0,Math.min(1,(ef.profileInfo-0.01)/0.05));
  const e=new Array(n);for(let i=0;i<n;i++)e[i]=W.b0+W.rr*O.rr[i]+g*(W.longstring*O.longstring[i]+W.person_total*O.person_total[i]);
  return e;}
function mcc(flag,y){let tp=0,tn=0,fp=0,fn=0;
  for(let i=0;i<y.length;i++){if(flag[i]){y[i]?tp++:fp++;}else{y[i]?fn++:tn++;}}
  const d=Math.sqrt((tp+fp)*(tp+fn)*(tn+fp)*(tn+fn));
  return {mcc:d?(tp*tn-fp*fn)/d:0,prec:(tp+fp)?tp/(tp+fp):0,rec:(tp+fn)?tp/(tp+fn):0};}

const D="../Dataset/gt_benchmark_candidates";
const SETS=[["TISP",D+"/tisp","_labels.csv"],["TUMI",D+"/tumi","_labels2.csv"],
  ["warnIPIP",D+"/warning_ipipneo300","_labels2.csv"],["smarvus",D+"/smarvus","_labels2.csv"],
  ["KayS2",D+"/kay_idris_idria/s2","_labels2.csv"],["KayS1",D+"/kay_idris_idria/s1","_labels2.csv"]];
console.log("=== PART A: external datasets (GT = failed >=2 checks, TISP >=1)  alpha=" + ALPHA + " ===");
console.log("dataset     n     prev  | OLD gate  MCC   | NEW tail gate: k     exp    p-value    open  MCC");
for(const [nm,dir,lf] of SETS){
  const P=R.parseCSV(fs.readFileSync(dir+"/_matrix.csv","utf8")),J=P.header.length;
  const rows=P.rows.map(r=>r.map(v=>{const x=Number(v);return Number.isFinite(x)?x:NaN;}));
  const y=R.parseCSV(fs.readFileSync(dir+"/"+lf,"utf8")).rows.map(r=>Number(r[0]));
  const n=Math.min(rows.length,y.length);
  const eta=etaOf(toMat(rows.slice(0,n),J),n,J,200);
  const af=R.autoFlag(eta,"standard");
  const mOld=mcc(af.flagged,y);
  const tg=tailGate(eta,ALPHA);
  const flagsNew=tg.open?eta.map(v=>v>tg.tau):new Array(n).fill(false);
  const mNew=mcc(flagsNew,y);
  console.log([nm.padEnd(10),String(n).padStart(5),(100*y.slice(0,n).reduce((s,v)=>s+v,0)/n).toFixed(1).padStart(5)+"%","|",
    (af.twoComp?"ON ":"OFF"),mOld.mcc.toFixed(3).padStart(6),"|",
    String(tg.k).padStart(5),tg.expected.toFixed(1).padStart(7),
    tg.p.toExponential(1).padStart(9),(tg.open?"OPEN":"shut").padStart(6),mNew.mcc.toFixed(3).padStart(6)].join(" "));
}
console.log("\n=== PART B: genuinely CLEAN simulated data - false alarm rate ===");
for(const [F,ipf,N] of [[10,6,222],[10,6,800],[20,6,2000],[30,6,4000]]){
  const J=F*ipf; const REPS=25; let oldFire=0,newFire=0,oldFrac=0,newFrac=0;
  for(let r=0;r<REPS;r++){
    const {rows}=genClean(F,ipf,N); const eta=etaOf(toMat(rows,J),N,J,100);
    const af=R.autoFlag(eta,"standard"); if(af.twoComp){oldFire++; oldFrac+=af.estRate;}
    const tg=tailGate(eta,ALPHA); if(tg.open){newFire++; newFrac+=tg.flaggedFrac;}
  }
  console.log("clean J="+String(J).padStart(3)+" n="+String(N).padStart(5)+
    " | OLD fires "+(100*oldFire/REPS).toFixed(0).padStart(3)+"% (flags "+(oldFire?(100*oldFrac/oldFire).toFixed(1):"0.0")+"%)"+
    " | NEW fires "+(100*newFire/REPS).toFixed(0).padStart(3)+"% (flags "+(newFire?(100*newFrac/newFire).toFixed(1):"0.0")+"%)");
}
