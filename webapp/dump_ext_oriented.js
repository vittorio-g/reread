/* dump_ext_oriented.js — like dump_ext.js but also dumps the oriented auxiliary
 * features (O.longstring, O.person_total) and the per-dataset gate g, so the
 * weight-robustness analysis (rounding / perturbation / fixed-vs-refit) can be
 * repeated on the EXTERNAL datasets. Output: ext_bench_oriented/. */
const R=require("./site/reread.js"); const W=R.WEIGHTS; const fs=require("fs"); const rng=R._rng(20260714);
const N_CAP=4000, ITER=200, OUT="ext_bench_oriented";
if(!fs.existsSync(OUT)) fs.mkdirSync(OUT);
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
const idx=[["name","n","J","rate","gate"]];
for(const [name,mf,lf] of SETS){
  const P=R.parseCSV(fs.readFileSync(mf,"utf8")),J=P.header.length;
  let rows=P.rows.map(r=>r.map(v=>{const x=Number(v);return Number.isFinite(x)?x:NaN;}));
  let lab=R.parseCSV(fs.readFileSync(lf,"utf8")).rows.map(r=>Number(r[0]));
  const m=Math.min(rows.length,lab.length); rows=rows.slice(0,m); lab=lab.slice(0,m);
  if(rows.length>N_CAP){const o=rows.map((_,i)=>i);for(let i=0;i<N_CAP;i++){const j=i+Math.floor(rng()*(o.length-i));const t=o[i];o[i]=o[j];o[j]=t;}
    const s=o.slice(0,N_CAP);rows=s.map(i=>rows[i]);lab=s.map(i=>lab[i]);}
  const n=rows.length,mat=new Float64Array(n*J);for(let i=0;i<n;i++)for(let j=0;j<J;j++)mat[i*J+j]=rows[i][j];
  const ef=R.ensembleFeatures(mat,n,J,{seed:1,iterations:ITER}),O=ef.oriented,g=Math.max(0,Math.min(1,(ef.profileInfo-0.01)/0.05));
  const Yo=["y,rr,longstring,person_total,eta"];
  for(let i=0;i<n;i++){const eta=W.b0+W.rr*O.rr[i]+g*(W.longstring*O.longstring[i]+W.person_total*O.person_total[i]);
    Yo.push([lab[i],O.rr[i].toFixed(6),O.longstring[i].toFixed(6),O.person_total[i].toFixed(6),eta.toFixed(6)].join(","));}
  fs.writeFileSync(OUT+"/"+name+"_y.csv",Yo.join("\n"));
  idx.push([name,n,J,(lab.reduce((s,v)=>s+v,0)/n).toFixed(3),g.toFixed(4)]);
  console.log(name.padEnd(14)+" n="+n+" J="+J+" rate="+(100*lab.reduce((s,v)=>s+v,0)/n).toFixed(0)+"% gate="+g.toFixed(3));
}
fs.writeFileSync(OUT+"/index.csv",idx.map(r=>r.join(",")).join("\n"));
console.log("\nwrote "+OUT+"/*");
