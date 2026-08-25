/* diag_bimodality.js — is eta a CONTINUUM (unimodal, wide) or a GROUP structure (bimodal)?
 * Candidate applicability check for the engine. Per dataset:
 *   valley depth = 1 - dens(valley)/min(dens(peak1),dens(peak2))  on the eta KDE
 *     (~0 = unimodal continuum -> counting machinery inapplicable; ->1 = separable group)
 *   mixture sep (gaussMix1D) and GT separation d as references.
 * Prediction: warning IPIP valley ~0; Kay S2 deep valley. */
const R = require("./site/rerere.js"); const W = R.WEIGHTS; const fs = require("fs");
const D = "../Dataset/gt_benchmark_candidates";
const SETS = [
  ["warning pooled", D+"/warning_ipipneo300", "_labels2.csv", "_matrix.csv"],
  ["Kay S2",  D+"/kay_idris_idria/s2", "_labels2.csv", "_matrix.csv"],
  ["Kay S1",  D+"/kay_idris_idria/s1", "_labels2.csv", "_matrix.csv"],
  ["smarvus", D+"/smarvus", "_labels2.csv", "_matrix.csv"],
  ["TISP DEU", D+"/tisp", "_labels_DEU.csv", "_matrix_DEU.csv"],
  ["TUMI",    D+"/tumi", "_labels2.csv", "_matrix.csv"],
  ["COVID-Dyn", D+"/covid_dynamic", "_labels_strict.csv", "_matrix.csv"],
];
function kdeGrid(x){const n=x.length,s=[...x].sort((a,b)=>a-b);
  const q1=s[Math.floor(.25*n)],q3=s[Math.floor(.75*n)],iqr=q3-q1;
  let mu=0;for(const v of x)mu+=v;mu/=n;let va=0;for(const v of x)va+=(v-mu)*(v-mu);va/=n;
  const sd=Math.sqrt(va);let h=0.9*Math.min(sd||1,(iqr/1.34)||(sd||1))*Math.pow(n,-0.2);
  if(!(h>0))h=(sd||1)*0.3;
  const G=512,lo=s[0]-h,hi=s[n-1]+h,dens=new Float64Array(G);
  for(let g=0;g<G;g++){const xg=lo+(hi-lo)*g/(G-1);let d=0;
    for(let i=0;i<n;i++){const z=(xg-x[i])/h;d+=Math.exp(-0.5*z*z);}dens[g]=d;}
  return dens;}
function valleyDepth(dens){const G=dens.length;
  const peaks=[];
  for(let g=1;g<G-1;g++)if(dens[g]>=dens[g-1]&&dens[g]>=dens[g+1])peaks.push([g,dens[g]]);
  peaks.sort((a,b)=>b[1]-a[1]);
  const dmax=peaks.length?peaks[0][1]:1;
  const main=peaks.filter(p=>p[1]>0.05*dmax);
  if(main.length<2)return {vd:0,npk:main.length};
  const [p1,p2]=[main[0],main[1]].sort((a,b)=>a[0]-b[0]);
  let vmin=Infinity;for(let g=p1[0];g<=p2[0];g++)if(dens[g]<vmin)vmin=dens[g];
  return {vd:1-vmin/Math.min(p1[1],p2[1]),npk:main.length};}
console.log("dataset        n     J   | valleyDepth #peaks | mixSep | GT d'");
for(const [nm,dir,lf,mf] of SETS){
  const P=R.parseCSV(fs.readFileSync(dir+"/"+mf,"utf8")),J=P.header.length;
  const rows=P.rows.map(r=>r.map(v=>{const x=Number(v);return Number.isFinite(x)?x:NaN;}));
  const y=R.parseCSV(fs.readFileSync(dir+"/"+lf,"utf8")).rows.map(r=>Number(r[0]));
  const n=Math.min(rows.length,y.length);
  const mat=new Float64Array(n*J);
  for(let i=0;i<n;i++)for(let j=0;j<J;j++)mat[i*J+j]=rows[i][j];
  const ef=R.ensembleFeatures(mat,n,J,{seed:1,iterations:150}),O=ef.oriented;
  const g=Math.max(0,Math.min(1,(ef.profileInfo-0.01)/0.05));
  const eta=new Array(n);
  for(let i=0;i<n;i++)eta[i]=W.b0+W.rr*O.rr[i]+g*(W.longstring*O.longstring[i]+W.person_total*O.person_total[i]);
  const {vd,npk}=valleyDepth(kdeGrid(eta));
  const gm=R.gaussMix1D(eta);const sep=(gm.m2-gm.m1)/(Math.sqrt((gm.v1+gm.v2)/2)||1);
  const e0=eta.filter((_,i)=>y[i]===0),e1=eta.filter((_,i)=>y[i]===1);
  const mu=a=>a.reduce((s,v)=>s+v,0)/a.length;
  const sd2=a=>{const m=mu(a);return a.reduce((s,v)=>s+(v-m)*(v-m),0)/a.length;};
  const dPr=e1.length?(mu(e1)-mu(e0))/Math.sqrt((sd2(e0)+sd2(e1))/2):NaN;
  console.log(nm.padEnd(14)+String(n).padStart(5)+" "+String(J).padStart(5)+" |   "+
    vd.toFixed(3)+"    "+String(npk).padStart(2)+"    |  "+sep.toFixed(2)+"  |  "+dPr.toFixed(2));
}
