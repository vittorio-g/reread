/* calibrate_gate.js — PHASE 1 of calibrating the excess-tail gate.
 * The gate assumed a Gaussian tail p0 = 1 - Phi(2.5) = .0062 above the one-sided cut. On clean data
 * eta actually has heavier tails than that, which is why the test false-alarmed. Here we MEASURE
 * the clean tail rate:
 *    (a) on simulated clean data across battery length and sample size
 *    (b) on the REAL clean sample: the 70 attentive respondents of Study 1
 * so that p0 can be set empirically instead of assumed.
 */
const R = require("./site/rerere.js"); const W = R.WEIGHTS; const fs = require("fs");
const Z = 2.5;
let SEED = 20260805; const rng = () => { SEED = (SEED*1103515245+12345)&0x7fffffff; return SEED/0x7fffffff; };
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
function etaOf(mat,n,J,it){const ef=R.ensembleFeatures(mat,n,J,{iterations:it,seed:1}),O=ef.oriented;
  const g=Math.max(0,Math.min(1,(ef.profileInfo-0.01)/0.05));
  const e=new Array(n);for(let i=0;i<n;i++)e[i]=W.b0+W.rr*O.rr[i]+g*(W.longstring*O.longstring[i]+W.person_total*O.person_total[i]);
  return e;}
function tailFrac(eta){ const af=R.autoFlag(eta,"standard");
  const tau=af.attentiveMode+Z*af.attentiveSigma;
  return eta.reduce((s,v)=>s+(v>tau?1:0),0)/eta.length; }
const q=(a,p)=>{const s=[...a].sort((x,y)=>x-y);return s[Math.min(s.length-1,Math.max(0,Math.round(p*(s.length-1))))];};
const mn=a=>a.reduce((s,v)=>s+v,0)/a.length;

console.log("=== clean tail rate above m + 2.5*sigma  (Gaussian assumption says .0062) ===");
console.log("source                   reps   mean    p50    p90    p95    max");
const all=[];
for(const [F,ipf,N,REPS] of [[10,6,200,30],[10,6,500,30],[10,6,2000,20],
                             [20,6,200,30],[20,6,500,30],[20,6,2000,20],
                             [30,6,500,20],[40,6,500,20]]){
  const J=F*ipf, fr=[];
  for(let r=0;r<REPS;r++){ const {rows}=genClean(F,ipf,N); fr.push(tailFrac(etaOf(toMat(rows,J),N,J,100))); }
  all.push(...fr);
  console.log(("sim J="+J+" n="+N).padEnd(24)+String(REPS).padStart(5)+
    ["",mn(fr),q(fr,.5),q(fr,.9),q(fr,.95),Math.max(...fr)].slice(1).map(v=>v.toFixed(4).padStart(7)).join(""));
}
/* real clean sample: the attentive respondents of Study 1 */
{ const Ms=R.parseCSV(fs.readFileSync("_study_matrix.csv","utf8")),J=Ms.header.length;
  const rows=Ms.rows.map(r=>r.map(v=>{const x=Number(v);return Number.isFinite(x)?x:NaN;}));
  const y=R.parseCSV(fs.readFileSync("_study_labels.csv","utf8")).rows.map(r=>Number(r[0]));
  const att=rows.filter((_,i)=>y[i]===0);
  const n=att.length, mat=new Float64Array(n*J);
  for(let i=0;i<n;i++)for(let j=0;j<J;j++)mat[i*J+j]=att[i][j];
  const f=tailFrac(etaOf(mat,n,J,400));
  console.log(("STUDY1 attentive-only J="+J+" n="+n).padEnd(24)+"    1"+f.toFixed(4).padStart(7).repeat(1));
  console.log("   (real clean sample: observed tail = "+(100*f).toFixed(2)+"%)");
}
console.log("\npooled simulated clean: mean="+mn(all).toFixed(4)+"  p90="+q(all,.9).toFixed(4)+
            "  p95="+q(all,.95).toFixed(4)+"  p99="+q(all,.99).toFixed(4)+"  max="+Math.max(...all).toFixed(4));
