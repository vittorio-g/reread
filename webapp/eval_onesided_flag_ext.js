/* eval_onesided_flag_ext.js — the same current-vs-one-sided flagging comparison on the 16
 * external Study-2 datasets. NOTE: their "ground truth" is a proxy criterion (attention checks,
 * self-report, speeding), not a true careless label, so MCC here is a RELATIVE comparison of two
 * thresholding rules on identical scores and identical labels, not an absolute quality figure. */
const R=require("./site/rerere.js"); const W=R.WEIGHTS; const fs=require("fs"); const rng=R._rng(20260717);
const N_CAP=4000, ITER=200, Z=2.5, D="../Dataset/";
const SETS=[
 ["Kay S2","core",D+"gt_benchmark_candidates/kay_idris_idria/kay_matrix_s2.csv",D+"gt_benchmark_candidates/kay_idris_idria/kay_labels_s2.csv"],
 ["warning IPIP","core",D+"gt_benchmark_candidates/warning_ipipneo300/_matrix.csv",D+"gt_benchmark_candidates/warning_ipipneo300/_labels.csv"],
 ["Kay S1","core",D+"gt_benchmark_candidates/kay_idris_idria/kay_matrix_s1.csv",D+"gt_benchmark_candidates/kay_idris_idria/kay_labels_s1.csv"],
 ["opsy 16PF","core",D+"gt_benchmark_candidates/opsy_16pf_matrix.csv",D+"gt_benchmark_candidates/opsy_16pf_labels.csv"],
 ["smarvus","core",D+"gt_benchmark_candidates/smarvus/smarvus_matrix.csv",D+"gt_benchmark_candidates/smarvus/smarvus_labels.csv"],
 ["Kay S6","core",D+"gt_benchmark_candidates/kay_idris_idria/kay_matrix_s6.csv",D+"gt_benchmark_candidates/kay_idris_idria/kay_labels_s6.csv"],
 ["Kay S5","core",D+"gt_benchmark_candidates/kay_idris_idria/kay_matrix_s5.csv",D+"gt_benchmark_candidates/kay_idris_idria/kay_labels_s5.csv"],
 ["Duckworth VCL","core",D+"gt_benchmark_candidates/duckworth_grit_vcl/duckworth_matrix.csv",D+"gt_benchmark_candidates/duckworth_grit_vcl/duckworth_labels.csv"],
 ["Douglas 2023","core",D+"gt_benchmark_candidates/douglas2023_dataquality/douglas_matrix.csv",D+"gt_benchmark_candidates/douglas2023_dataquality/douglas_labels.csv"],
 ["Krause youth","core",D+"gt_benchmark_candidates/krause_ier_youth/krause_matrix.csv",D+"gt_benchmark_candidates/krause_ier_youth/krause_labels.csv"],
 ["Mastroianni","bound",D+"learning_to_pay_attention/mastroianni2022_attitude/_matrix.csv",D+"learning_to_pay_attention/mastroianni2022_attitude/_labels.csv"],
 ["Ivanov 2021","bound",D+"learning_to_pay_attention/ivanov2021_racialresentment/_matrix.csv",D+"learning_to_pay_attention/ivanov2021_racialresentment/_labels.csv"],
 ["O'Grady 2019","bound",D+"learning_to_pay_attention/ogrady2019_moralfoundations/_matrix.csv",D+"learning_to_pay_attention/ogrady2019_moralfoundations/_labels.csv"],
 ["Pennycook 2020","bound",D+"learning_to_pay_attention/pennycook2020_covid/_matrix.csv",D+"learning_to_pay_attention/pennycook2020_covid/_labels.csv"],
 ["Moss 2023","bound",D+"learning_to_pay_attention/moss2023_ethical/_matrix.csv",D+"learning_to_pay_attention/moss2023_ethical/_labels.csv"],
 ["Alvarez 2019","bound",D+"learning_to_pay_attention/alvarez2019_inattentive/_matrix.csv",D+"learning_to_pay_attention/alvarez2019_inattentive/_labels.csv"],
];
const med=a=>{const s=[...a].sort((x,y)=>x-y);const h=s.length>>1;return s.length%2?s[h]:(s[h-1]+s[h])/2;};
function leftMode(x){const n=x.length,s=[...x].sort((a,b)=>a-b);
  const q1=s[Math.floor(.25*n)],q3=s[Math.floor(.75*n)],iqr=q3-q1;
  const mu=x.reduce((a,b)=>a+b,0)/n, sd=Math.sqrt(x.reduce((a,b)=>a+(b-mu)*(b-mu),0)/n);
  let h=0.9*Math.min(sd||1,(iqr/1.34)||(sd||1))*Math.pow(n,-0.2); if(!(h>0))h=(sd||1)*0.3;
  const G=512,lo=s[0]-h,hi=s[n-1]+h,dens=new Float64Array(G);
  for(let g=0;g<G;g++){const xg=lo+(hi-lo)*g/(G-1);let d=0;for(let i=0;i<n;i++){const z=(xg-x[i])/h;d+=Math.exp(-0.5*z*z);}dens[g]=d;}
  let dmax=0;for(let g=0;g<G;g++)if(dens[g]>dmax)dmax=dens[g];
  for(let g=1;g<G-1;g++) if(dens[g]>=dens[g-1]&&dens[g]>=dens[g+1]&&dens[g]>0.15*dmax) return lo+(hi-lo)*g/(G-1);
  return med(x);}
function cutOneSided(e,z){const m=leftMode(e),below=e.filter(v=>v<m);
  let sd=below.length?1.4826*med(below.map(v=>m-v)):0; if(!(sd>0))sd=1e-9; return m+z*sd;}
function gateOK(gm){const p=Math.sqrt((gm.v1+gm.v2)/2)||1;return (gm.bic1-gm.bic2>2)&&gm.pi>0.01&&gm.pi<0.7&&(gm.m2-gm.m1)/p>1.0;}
function mccOf(flag,lab){let tp=0,tn=0,fp=0,fn=0;
  for(let i=0;i<lab.length;i++){if(flag[i]){lab[i]?tp++:fp++;}else{lab[i]?fn++:tn++;}}
  const d=Math.sqrt((tp+fp)*(tp+fn)*(tn+fp)*(tn+fn));
  return {mcc:d?(tp*tn-fp*fn)/d:0,prec:(tp+fp)?tp/(tp+fp):0,rec:(tp+fn)?tp/(tp+fn):0,frac:(tp+fp)/lab.length};}
const out=[["dataset","group","n","J","gt_rate","mcc_cur","mcc_os","prec_cur","prec_os","rec_cur","rec_os","flag_cur","flag_os"]];
let sc=0,so=0,cnt=0;
console.log("dataset            n     J   GT    | MCC cur -> OS   | prec cur/OS | rec cur/OS | flagged cur/OS");
for(const [name,tag,mf,lf] of SETS){
  let P,lab; try{P=R.parseCSV(fs.readFileSync(mf,"utf8"));lab=R.parseCSV(fs.readFileSync(lf,"utf8")).rows.map(r=>Number(r[0]));}
  catch(e){console.log(name+" SKIP");continue;}
  const J=P.header.length;
  let rows=P.rows.map(r=>r.map(v=>{const x=Number(v);return Number.isFinite(x)?x:NaN;}));
  const m0=Math.min(rows.length,lab.length); rows=rows.slice(0,m0); lab=lab.slice(0,m0);
  if(rows.length>N_CAP){const o=rows.map((_,i)=>i);for(let i=0;i<N_CAP;i++){const j=i+Math.floor(rng()*(o.length-i));const t=o[i];o[i]=o[j];o[j]=t;}
    const s=o.slice(0,N_CAP);rows=s.map(i=>rows[i]);lab=s.map(i=>lab[i]);}
  const n=rows.length,mat=new Float64Array(n*J);
  for(let i=0;i<n;i++)for(let j=0;j<J;j++)mat[i*J+j]=rows[i][j];
  const ef=R.ensembleFeatures(mat,n,J,{iterations:ITER,seed:1}),O=ef.oriented;
  const g=Math.max(0,Math.min(1,(ef.profileInfo-0.01)/0.05));
  const e=new Array(n);for(let i=0;i<n;i++)e[i]=W.b0+W.rr*O.rr[i]+g*(W.longstring*O.longstring[i]+W.person_total*O.person_total[i]);
  const gm=R.gaussMix1D(e), cutCur=gateOK(gm)?(gm.m1+Z*Math.sqrt(gm.v1)):Infinity, cutOS=cutOneSided(e,Z);
  const A=mccOf(e.map(v=>v>cutCur),lab), B=mccOf(e.map(v=>v>cutOS),lab);
  const gt=lab.reduce((s,v)=>s+v,0)/n; sc+=A.mcc; so+=B.mcc; cnt++;
  out.push([name,tag,n,J,gt.toFixed(4),A.mcc.toFixed(4),B.mcc.toFixed(4),A.prec.toFixed(4),B.prec.toFixed(4),A.rec.toFixed(4),B.rec.toFixed(4),A.frac.toFixed(4),B.frac.toFixed(4)]);
  console.log(name.padEnd(16)+String(n).padStart(5)+String(J).padStart(5)+(100*gt).toFixed(0).padStart(5)+"% |  "+
    A.mcc.toFixed(3).padStart(6)+" -> "+B.mcc.toFixed(3).padStart(6)+"  | "+A.prec.toFixed(2)+"/"+B.prec.toFixed(2)+
    "   | "+A.rec.toFixed(2)+"/"+B.rec.toFixed(2)+"  | "+(100*A.frac).toFixed(1)+"%/"+(100*B.frac).toFixed(1)+"%");
}
console.log("\nmean MCC: current "+(sc/cnt).toFixed(3)+"   one-sided "+(so/cnt).toFixed(3));
fs.writeFileSync("onesided_flag_ext.csv",out.map(r=>r.join(",")).join("\n"));
console.log("wrote onesided_flag_ext.csv");
