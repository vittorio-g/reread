/* test_fixes_matrix.js — evaluate each candidate algorithm fix SEPARATELY, on the cleaned subsets
 * (single country / single experimental arm / n>J halves) and on clean simulated data.
 *
 *   shipped : gate = BIC-conjunction        , mode = left-most kde peak >15% of max
 *   B1      : gate = BIC-conjunction        , mode = TALLEST kde peak (capped at the median)
 *   B3      : gate = BIC-conj OR extreme tail (p < 1e-6), mode = left-most
 *   B1+B3   : both fixes together
 * The cut is always m + 2.5*sigma with sigma from the left side only.
 */
const R = require("./site/rerere.js"); const W = R.WEIGHTS; const fs = require("fs");
const { binomUpper } = require("./gate_tailtest.js");
const Z = 2.5, P0 = 0.006209665, TAILALPHA = 1e-6;
let SEED=13579; const rng=()=>{SEED=(SEED*1103515245+12345)&0x7fffffff;return SEED/0x7fffffff;};
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
function modeLeft(eta){const {s,dens,lo,hi,G}=kde(eta);let dmax=0;for(let g=0;g<G;g++)if(dens[g]>dmax)dmax=dens[g];
  for(let g=1;g<G-1;g++)if(dens[g]>=dens[g-1]&&dens[g]>=dens[g+1]&&dens[g]>0.15*dmax)return lo+(hi-lo)*g/(G-1);
  return s[Math.floor(eta.length/2)];}
function modeGlobal(eta){const {s,dens,lo,hi,G}=kde(eta);let best=0,bi=0;
  for(let g=0;g<G;g++)if(dens[g]>best){best=dens[g];bi=g;}
  let m=lo+(hi-lo)*bi/(G-1);const med=s[Math.floor(eta.length/2)];return m>med?med:m;}
function bicGate(eta){const gm=R.gaussMix1D(eta);const p=Math.sqrt((gm.v1+gm.v2)/2)||1;
  return (gm.bic1-gm.bic2>2)&&gm.pi>0.01&&gm.pi<0.7&&(gm.m2-gm.m1)/p>1.0;}
function run(eta,useGlobal,useTailOr){
  const m=useGlobal?modeGlobal(eta):modeLeft(eta);
  const {sd}=kde(eta); const sig=sigmaBelow(eta,m,sd); const tau=m+Z*sig;
  const k=eta.reduce((s,v)=>s+(v>tau?1:0),0);
  let open=bicGate(eta);
  if(!open&&useTailOr) open = binomUpper(k,eta.length,P0) < TAILALPHA;
  return {open, flags: open?eta.map(v=>v>tau):new Array(eta.length).fill(false), frac:k/eta.length};
}
function mcc(flag,y){let tp=0,tn=0,fp=0,fn=0;
  for(let i=0;i<y.length;i++){if(flag[i]){y[i]?tp++:fp++;}else{y[i]?fn++:tn++;}}
  const d=Math.sqrt((tp+fp)*(tp+fn)*(tn+fp)*(tn+fn));return d?(tp*tn-fp*fn)/d:0;}
const VAR=[["shipped",false,false],["B1 global",true,false],["B3 tailOR",false,true],["B1+B3",true,true]];
const D="../Dataset/gt_benchmark_candidates";
const SETS=[
 ["TISP DEU",D+"/tisp","_labels_DEU.csv","_matrix_DEU.csv"],
 ["TISP POL",D+"/tisp","_labels_POL.csv","_matrix_POL.csv"],
 ["TISP AUS",D+"/tisp","_labels_AUS.csv","_matrix_AUS.csv"],
 ["smarvus Eng",D+"/smarvus_England","_labels2.csv","_matrix.csv"],
 ["smarvus Egy",D+"/smarvus_Egypt","_labels2.csv","_matrix.csv"],
 ["smarvus Can",D+"/smarvus_Canada","_labels2.csv","_matrix.csv"],
 ["warn control",D+"/warning_control","_labels2.csv","_matrix.csv"],
 ["warn act",D+"/warning_warn_act","_labels2.csv","_matrix.csv"],
 ["warn pas",D+"/warning_warn_pas","_labels2.csv","_matrix.csv"],
 ["TUMI A",D+"/tumi_A","_labels2.csv","_matrix.csv"],
 ["TUMI B",D+"/tumi_B","_labels2.csv","_matrix.csv"],
 ["Kay S2",D+"/kay_idris_idria/s2","_labels2.csv","_matrix.csv"],
 ["Kay S1",D+"/kay_idris_idria/s1","_labels2.csv","_matrix.csv"]];
console.log("=== MCC by variant (gate closed => MCC 0) ===");
console.log("dataset       prev  |  shipped  B1 global  B3 tailOR   B1+B3   | gates open");
const tot={};VAR.forEach(v=>tot[v[0]]=0);
for(const [nm,dir,lf,mf] of SETS){
  const P=R.parseCSV(fs.readFileSync(dir+"/"+mf,"utf8")),J=P.header.length;
  const rows=P.rows.map(r=>r.map(v=>{const x=Number(v);return Number.isFinite(x)?x:NaN;}));
  const y=R.parseCSV(fs.readFileSync(dir+"/"+lf,"utf8")).rows.map(r=>Number(r[0]));
  const n=Math.min(rows.length,y.length);
  const eta=etaOf(toMat(rows.slice(0,n),J),n,J,200);
  const cells=[],gates=[];
  for(const [vn,gl,to] of VAR){const r=run(eta,gl,to);const m=mcc(r.flags,y.slice(0,n));
    cells.push(m.toFixed(3).padStart(8));gates.push(r.open?"O":"-");tot[vn]+=m;}
  console.log(nm.padEnd(13)+(100*y.slice(0,n).reduce((a,b)=>a+b,0)/n).toFixed(1).padStart(5)+"% |"+cells.join(" ")+"   | "+gates.join(" "));
}
console.log("MEAN".padEnd(13)+"      |"+VAR.map(v=>(tot[v[0]]/SETS.length).toFixed(3).padStart(8)).join(" "));
console.log("\n=== CLEAN simulated data: how often does each variant flag anyone? ===");
for(const [F,ipf,N,REPS] of [[10,6,200,30],[10,6,500,25],[20,6,2000,12]]){
  const J=F*ipf;const fire={};VAR.forEach(v=>fire[v[0]]=0);const frac={};VAR.forEach(v=>frac[v[0]]=0);
  for(let r=0;r<REPS;r++){const {rows}=genClean(F,ipf,N);const eta=etaOf(toMat(rows,J),N,J,100);
    for(const [vn,gl,to] of VAR){const res=run(eta,gl,to);
      if(res.open){fire[vn]++;frac[vn]+=res.frac;}}}
  console.log(("clean J="+J+" n="+N).padEnd(18)+VAR.map(v=>{
    const f=fire[v[0]];return (v[0]+": "+(100*f/REPS).toFixed(0)+"% fire, "+(f?(100*frac[v[0]]/f).toFixed(1):"0.0")+"% flagged");}).join(" | "));
}
