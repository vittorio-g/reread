/* test_cut_ext.js — same cut-strategy comparison on the 12 external GT datasets (natural prevalence).
 * MCC vs each dataset's GT label. Big datasets subsampled to N_CAP (prevalence preserved). */
const R=require("./site/reread.js"); const W=R.WEIGHTS; const fs=require("fs"); const rng=R._rng(20260714);
const N_CAP=5000, ITER=150;
const ncdf=x=>{const t=1/(1+0.2316419*Math.abs(x)),d=0.3989423*Math.exp(-x*x/2);const p=d*t*(0.3193815+t*(-0.3565638+t*(1.781478+t*(-1.821256+t*1.330274))));return x>0?1-p:p;};
function gateOK(gm){const pooled=Math.sqrt((gm.v1+gm.v2)/2)||1,sep=(gm.m2-gm.m1)/pooled;return (gm.bic1-gm.bic2>2)&&gm.pi>0.01&&gm.pi<0.7&&sep>1.0;}
function mcc(flag,lab){let tp=0,fp=0,tn=0,fn=0;for(let i=0;i<lab.length;i++){if(flag[i]){lab[i]?tp++:fp++;}else{lab[i]?fn++:tn++;}}
  const d=Math.sqrt((tp+fp)*(tp+fn)*(tn+fp)*(tn+fn));return d===0?0:(tp*tn-fp*fn)/d;}
function topK(e,k){const idx=e.map((v,i)=>[v,i]).sort((a,b)=>b[0]-a[0]);const f=new Array(e.length).fill(false);for(let i=0;i<k;i++)f[idx[i][1]]=true;return f;}
function oracleMCC(e,lab){const th=[...new Set(e)].sort((a,b)=>a-b);let best=-1;for(const t of th){const f=e.map(v=>v>=t);const m=mcc(f,lab);if(m>best)best=m;}return best;}
function etaOf(mat,n,J){const ef=R.ensembleFeatures(mat,n,J,{seed:1,iterations:ITER}),O=ef.oriented,g=Math.max(0,Math.min(1,(ef.profileInfo-0.01)/0.05));
  const e=new Array(n);for(let i=0;i<n;i++)e[i]=W.b0+W.rr*O.rr[i]+g*(W.longstring*O.longstring[i]+W.person_total*O.person_total[i]);return e;}
function strat(e,lab){
  const gm=R.gaussMix1D(e),ok=gateOK(gm),n=e.length,s1=Math.sqrt(gm.v1),s2=Math.sqrt(gm.v2);
  const std=ok?e.map(v=>v>gm.m1+2.5*s1):e.map(()=>false), high=ok?e.map(v=>v>gm.m1+1.5*s1):e.map(()=>false);
  const pi=ok?gm.pi:0; let leak=pi;
  if(ok){const cut=gm.m1+2.5*s1,Satt=1-ncdf((cut-gm.m1)/s1),Scare=1-ncdf((cut-gm.m2)/s2),den=Scare-Satt,obs=std.filter(Boolean).length/n;
    leak=Math.abs(den)>0.05?Math.max(0,Math.min(1,(obs-Satt)/den)):pi;}
  return {std:mcc(std,lab),high:mcc(high,lab),picut:mcc(topK(e,Math.round(pi*n)),lab),leakcut:mcc(topK(e,Math.round(leak*n)),lab),
    oracle:oracleMCC(e,lab),piHat:pi,gate:ok?1:0};
}
const D="../Dataset/", SETS=[
 ["warning300",D+"gt_benchmark_candidates/warning_ipipneo300/_matrix.csv",D+"gt_benchmark_candidates/warning_ipipneo300/_labels.csv"],
 ["smarvus",D+"gt_benchmark_candidates/smarvus/smarvus_matrix.csv",D+"gt_benchmark_candidates/smarvus/smarvus_labels.csv"],
 ["opsy_16pf",D+"gt_benchmark_candidates/opsy_16pf_matrix.csv",D+"gt_benchmark_candidates/opsy_16pf_labels.csv"],
 ["douglas",D+"gt_benchmark_candidates/douglas2023_dataquality/douglas_matrix.csv",D+"gt_benchmark_candidates/douglas2023_dataquality/douglas_labels.csv"],
 ["duckworth",D+"gt_benchmark_candidates/duckworth_grit_vcl/duckworth_matrix.csv",D+"gt_benchmark_candidates/duckworth_grit_vcl/duckworth_labels.csv"],
 ["krause",D+"gt_benchmark_candidates/krause_ier_youth/krause_matrix.csv",D+"gt_benchmark_candidates/krause_ier_youth/krause_labels.csv"],
 ["alvarez19",D+"learning_to_pay_attention/alvarez2019_inattentive/_matrix.csv",D+"learning_to_pay_attention/alvarez2019_inattentive/_labels.csv"],
 ["ivanov21",D+"learning_to_pay_attention/ivanov2021_racialresentment/_matrix.csv",D+"learning_to_pay_attention/ivanov2021_racialresentment/_labels.csv"],
 ["mastroianni22",D+"learning_to_pay_attention/mastroianni2022_attitude/_matrix.csv",D+"learning_to_pay_attention/mastroianni2022_attitude/_labels.csv"],
 ["moss23",D+"learning_to_pay_attention/moss2023_ethical/_matrix.csv",D+"learning_to_pay_attention/moss2023_ethical/_labels.csv"],
 ["ogrady19",D+"learning_to_pay_attention/ogrady2019_moralfoundations/_matrix.csv",D+"learning_to_pay_attention/ogrady2019_moralfoundations/_labels.csv"],
 ["pennycook20",D+"learning_to_pay_attention/pennycook2020_covid/_matrix.csv",D+"learning_to_pay_attention/pennycook2020_covid/_labels.csv"],
];
console.log("EXTERNAL GT datasets — MCC by cut strategy (cap n="+N_CAP+")\n");
console.log("dataset          n     J   rate  piHat | Standard  High   pi-cut leak-cut | oracle | best-auto");
const out=[["dataset","n","J","rate","piHat","std","high","picut","leakcut","oracle"]];
for(const [name,mf,lf] of SETS){
  try{
  const P=R.parseCSV(fs.readFileSync(mf,"utf8")),J=P.header.length;
  let rows=P.rows.map(r=>r.map(v=>{const x=Number(v);return Number.isFinite(x)?x:NaN;}));
  let lab=R.parseCSV(fs.readFileSync(lf,"utf8")).rows.map(r=>Number(r[0]));
  const m=Math.min(rows.length,lab.length); rows=rows.slice(0,m); lab=lab.slice(0,m);
  if(rows.length>N_CAP){const ord=rows.map((_,i)=>i);for(let i=0;i<N_CAP;i++){const j=i+Math.floor(rng()*(ord.length-i));const t=ord[i];ord[i]=ord[j];ord[j]=t;}
    const sel=ord.slice(0,N_CAP);rows=sel.map(i=>rows[i]);lab=sel.map(i=>lab[i]);}
  const n=rows.length,mat=new Float64Array(n*J);for(let i=0;i<n;i++)for(let j=0;j<J;j++)mat[i*J+j]=rows[i][j];
  const rate=lab.reduce((s,v)=>s+v,0)/n, s=strat(etaOf(mat,n,J),lab);
  out.push([name,n,J,rate.toFixed(3),s.piHat.toFixed(3),s.std.toFixed(3),s.high.toFixed(3),s.picut.toFixed(3),s.leakcut.toFixed(3),s.oracle.toFixed(3)]);
  const autos=["std","high","picut","leakcut"],ba=autos.reduce((a,b)=>s[b]>s[a]?b:a,"std");
  console.log(name.padEnd(15)+String(n).padStart(6)+String(J).padStart(5)+(100*rate).toFixed(0).padStart(5)+"%"+s.piHat.toFixed(2).padStart(7)+" | "+
    s.std.toFixed(3)+"   "+s.high.toFixed(3)+"  "+s.picut.toFixed(3)+"  "+s.leakcut.toFixed(3)+"  |  "+s.oracle.toFixed(3)+" |  "+ba+"("+s[ba].toFixed(3)+")");
  }catch(err){console.log(name.padEnd(15)+" ERROR "+err.message);}
}
fs.writeFileSync("test_cut_ext.csv",out.map(r=>r.join(",")).join("\n"));
console.log("\nwrote test_cut_ext.csv");
