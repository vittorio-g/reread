/* sweep_gate_prevalence.js — does the gate close because the careless rate is genuinely LOW,
 * or does it just fail?  Controlled test on REAL respondents: hold the attentive pool fixed and
 * mix in k real careless respondents to hit a target prevalence, then record, at each rate:
 *   - does the gate open (shipped BIC-conjunction, and the B1+B3 variant)
 *   - the reported pi-hat, the share actually flagged, and MCC
 * so we can see the operating characteristic of the gate against a KNOWN prevalence.
 */
const R = require("./site/reread.js"); const W = R.WEIGHTS; const fs = require("fs");
const { binomUpper } = require("./gate_tailtest.js");
const Z = 2.5, P0 = 0.006209665, TAILALPHA = 1e-6;
let SEED = 24680; const rng = () => { SEED = (SEED*1103515245+12345)&0x7fffffff; return SEED/0x7fffffff; };
function kde(eta){const n=eta.length,s=[...eta].sort((a,b)=>a-b);
  const q1=s[Math.floor(.25*n)],q3=s[Math.floor(.75*n)],iqr=q3-q1;
  let mu=0;for(const v of eta)mu+=v;mu/=n;let va=0;for(const v of eta)va+=(v-mu)*(v-mu);va/=n;
  const sd=Math.sqrt(va);let h=0.9*Math.min(sd||1,(iqr/1.34)||(sd||1))*Math.pow(n,-0.2);
  if(!(h>0))h=(sd||1)*0.3;
  const G=512,lo=s[0]-h,hi=s[n-1]+h,dens=new Float64Array(G);
  for(let g=0;g<G;g++){const xg=lo+(hi-lo)*g/(G-1);let d=0;
    for(let i=0;i<n;i++){const z=(xg-eta[i])/h;d+=Math.exp(-0.5*z*z);}dens[g]=d;}
  return {s,dens,lo,hi,G,sd};}
function sigmaBelow(eta,m,sd){const b=[];for(const v of eta)if(v<m)b.push(m-v);b.sort((a,x)=>a-x);
  let s=b.length?1.4826*b[Math.floor(b.length/2)]:0;if(!(s>0))s=(sd||1e-9);return s;}
function modeGlobal(eta){const {s,dens,lo,hi,G}=kde(eta);let best=0,bi=0;
  for(let g=0;g<G;g++)if(dens[g]>best){best=dens[g];bi=g;}
  let m=lo+(hi-lo)*bi/(G-1);const med=s[Math.floor(eta.length/2)];return m>med?med:m;}
function bicGate(eta){const gm=R.gaussMix1D(eta);const p=Math.sqrt((gm.v1+gm.v2)/2)||1;
  return {open:(gm.bic1-gm.bic2>2)&&gm.pi>0.01&&gm.pi<0.7&&(gm.m2-gm.m1)/p>1.0,
          dbic:gm.bic1-gm.bic2, sep:(gm.m2-gm.m1)/p, w:gm.pi};}
const PHI1=0.8413447;
function evaluate(eta,y){
  const bg=bicGate(eta);
  const m=modeGlobal(eta),{sd}=kde(eta),sig=sigmaBelow(eta,m,sd),tau=m+Z*sig;
  const k=eta.reduce((s,v)=>s+(v>tau?1:0),0);
  const ptail=binomUpper(k,eta.length,P0);
  const openNew = bg.open || ptail<TAILALPHA;
  const nBelow=eta.reduce((s,v)=>s+(v<m+sig?1:0),0);
  const piHat=Math.max(0,Math.min(1,1-(nBelow/PHI1)/eta.length));
  const flags=openNew?eta.map(v=>v>tau):new Array(eta.length).fill(false);
  let tp=0,tn=0,fp=0,fn=0;
  for(let i=0;i<y.length;i++){if(flags[i]){y[i]?tp++:fp++;}else{y[i]?fn++:tn++;}}
  const d=Math.sqrt((tp+fp)*(tp+fn)*(tn+fp)*(tn+fn));
  return {bicOpen:bg.open,dbic:bg.dbic,sep:bg.sep,newOpen:openNew,ptail,
          piHat, flagged:k/eta.length, mcc:d?(tp*tn-fp*fn)/d:0};
}
const D="../Dataset/gt_benchmark_candidates";
const SETS=[["TISP DEU",D+"/tisp","_labels_DEU.csv","_matrix_DEU.csv"],
            ["smarvus",D+"/smarvus","_labels2.csv","_matrix.csv"]];
for(const [nm,dir,lf,mf] of SETS){
  const P=R.parseCSV(fs.readFileSync(dir+"/"+mf,"utf8")),J=P.header.length;
  const rows=P.rows.map(r=>r.map(v=>{const x=Number(v);return Number.isFinite(x)?x:NaN;}));
  const y=R.parseCSV(fs.readFileSync(dir+"/"+lf,"utf8")).rows.map(r=>Number(r[0]));
  const n0=Math.min(rows.length,y.length);
  const att=[],car=[];
  for(let i=0;i<n0;i++)(y[i]?car:att).push(rows[i]);
  const NATT=Math.min(att.length,1500);
  console.log("\n=== "+nm+" (J="+J+", attentive pool "+NATT+", careless pool "+car.length+") ===");
  console.log(" true%   n   | BIC gate  dBIC     sep   | NEW gate  p(tail)  | piHat  flagged   MCC");
  for(const target of [0,0.02,0.05,0.08,0.12,0.18,0.25,0.35]){
    const k=Math.round(NATT*target/(1-target));
    if(k>car.length) continue;
    const ai=att.slice(0,NATT);
    const ci=car.slice(); for(let i=0;i<k;i++){const j=i+Math.floor(rng()*(ci.length-i));const t=ci[i];ci[i]=ci[j];ci[j]=t;}
    const sub=[...ai,...ci.slice(0,k)], lab=[...new Array(NATT).fill(0),...new Array(k).fill(1)];
    const n=sub.length,mat=new Float64Array(n*J);
    for(let i=0;i<n;i++)for(let j=0;j<J;j++)mat[i*J+j]=sub[i][j];
    const ef=R.ensembleFeatures(mat,n,J,{seed:1,iterations:200}),O=ef.oriented;
    const g=Math.max(0,Math.min(1,(ef.profileInfo-0.01)/0.05));
    const eta=new Array(n);
    for(let i=0;i<n;i++)eta[i]=W.b0+W.rr*O.rr[i]+g*(W.longstring*O.longstring[i]+W.person_total*O.person_total[i]);
    const r=evaluate(eta,lab);
    console.log([(100*k/n).toFixed(1).padStart(5)+"%",String(n).padStart(5),"|",
      (r.bicOpen?"OPEN":"shut").padStart(5),r.dbic.toFixed(0).padStart(7),r.sep.toFixed(2).padStart(7),"|",
      (r.newOpen?"OPEN":"shut").padStart(5),r.ptail.toExponential(1).padStart(9),"|",
      (100*r.piHat).toFixed(1).padStart(5)+"%",(100*r.flagged).toFixed(1).padStart(6)+"%",
      r.mcc.toFixed(3).padStart(6)].join(" "));
  }
}
