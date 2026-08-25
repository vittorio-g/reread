/* dump_ext2.js — external competitive benchmark on the SAME datasets as Study 2 tables
 * (tab:ext core + tab:ltpa boundary) that have archived matrix+labels. Full n (cap 6000,
 * prevalence preserved). Dumps matrix + labels + shipped rr/eta on identical rows for R. */
const R=require("./site/reread.js"); const W=R.WEIGHTS; const fs=require("fs"); const rng=R._rng(20260715);
const N_CAP=Infinity, ITER=200, OUT="ext_bench2"; if(!fs.existsSync(OUT)) fs.mkdirSync(OUT);  // full n for reproducible Table 3
const G="../Dataset/gt_benchmark_candidates/", L="../Dataset/learning_to_pay_attention/";
const SETS=[
 ["kay_s2",G+"kay_idris_idria/kay_matrix_s2.csv",G+"kay_idris_idria/kay_labels_s2.csv","core"],
 ["warning300",G+"warning_ipipneo300/_matrix.csv",G+"warning_ipipneo300/_labels.csv","core"],
 ["kay_s1",G+"kay_idris_idria/kay_matrix_s1.csv",G+"kay_idris_idria/kay_labels_s1.csv","core"],
 ["opsy_16pf",G+"opsy_16pf_matrix.csv",G+"opsy_16pf_labels.csv","core"],
 ["smarvus",G+"smarvus/smarvus_matrix.csv",G+"smarvus/smarvus_labels.csv","core"],
 ["kay_s6",G+"kay_idris_idria/kay_matrix_s6.csv",G+"kay_idris_idria/kay_labels_s6.csv","core"],
 ["kay_s5",G+"kay_idris_idria/kay_matrix_s5.csv",G+"kay_idris_idria/kay_labels_s5.csv","core"],
 ["duckworth",G+"duckworth_grit_vcl/duckworth_matrix.csv",G+"duckworth_grit_vcl/duckworth_labels.csv","core"],
 ["douglas",G+"douglas2023_dataquality/douglas_matrix.csv",G+"douglas2023_dataquality/douglas_labels.csv","core"],
 ["krause",G+"krause_ier_youth/krause_matrix.csv",G+"krause_ier_youth/krause_labels.csv","core"],
 ["ivanov21",L+"ivanov2021_racialresentment/_matrix.csv",L+"ivanov2021_racialresentment/_labels.csv","bound"],
 ["moss23",L+"moss2023_ethical/_matrix.csv",L+"moss2023_ethical/_labels.csv","bound"],
 ["alvarez19",L+"alvarez2019_inattentive/_matrix.csv",L+"alvarez2019_inattentive/_labels.csv","bound"],
 ["ogrady19",L+"ogrady2019_moralfoundations/_matrix.csv",L+"ogrady2019_moralfoundations/_labels.csv","bound"],
 ["mastroianni22",L+"mastroianni2022_attitude/_matrix.csv",L+"mastroianni2022_attitude/_labels.csv","bound"],
 ["pennycook20",L+"pennycook2020_covid/_matrix.csv",L+"pennycook2020_covid/_labels.csv","bound"],
 ["buchanan18",L+"buchanan2018_lowquality/_matrix.csv",L+"buchanan2018_lowquality/_labels.csv","bound"],
];
const idx=[["name","n","J","rate","kind"]];
for(const [name,mf,lf,kind] of SETS){
  if(!fs.existsSync(mf)||!fs.existsSync(lf)){console.log(name+" MISSING file, skip");continue;}
  const P=R.parseCSV(fs.readFileSync(mf,"utf8")),J=P.header.length;
  let rows=P.rows.map(r=>r.map(v=>{const x=Number(v);return Number.isFinite(x)?x:NaN;}));
  let lab=R.parseCSV(fs.readFileSync(lf,"utf8")).rows.map(r=>Number(r[0]));
  const m=Math.min(rows.length,lab.length); rows=rows.slice(0,m); lab=lab.slice(0,m);
  if(rows.length>N_CAP){const o=rows.map((_,i)=>i);for(let i=0;i<N_CAP;i++){const j=i+Math.floor(rng()*(o.length-i));const t=o[i];o[i]=o[j];o[j]=t;}
    const s=o.slice(0,N_CAP);rows=s.map(i=>rows[i]);lab=s.map(i=>lab[i]);}
  const n=rows.length,mat=new Float64Array(n*J);for(let i=0;i<n;i++)for(let j=0;j<J;j++)mat[i*J+j]=rows[i][j];
  const ef=R.ensembleFeatures(mat,n,J,{seed:1,iterations:ITER}),O=ef.oriented,g=Math.max(0,Math.min(1,(ef.profileInfo-0.01)/0.05));
  const Mo=[];for(let i=0;i<n;i++)Mo.push(rows[i].map(v=>Number.isFinite(v)?v:"").join(","));
  const Yo=["y,rr,eta,pt"];for(let i=0;i<n;i++){const eta=W.b0+W.rr*O.rr[i]+g*(W.longstring*O.longstring[i]+W.person_total*O.person_total[i]);
    Yo.push(lab[i]+","+O.rr[i].toFixed(6)+","+eta.toFixed(6)+","+O.person_total[i].toFixed(6));}
  fs.writeFileSync(OUT+"/"+name+"_M.csv",Mo.join("\n"));
  fs.writeFileSync(OUT+"/"+name+"_y.csv",Yo.join("\n"));
  idx.push([name,n,J,(lab.reduce((s,v)=>s+v,0)/n).toFixed(3),kind]);
  console.log(name.padEnd(14)+" n="+n+" J="+J+" rate="+(100*lab.reduce((s,v)=>s+v,0)/n).toFixed(0)+"% gate="+g.toFixed(2));
}
fs.writeFileSync(OUT+"/index.csv",idx.map(r=>r.join(",")).join("\n"));
console.log("\nwrote "+OUT+"/* ("+(idx.length-1)+" datasets)");
