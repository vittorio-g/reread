/* exp_combine.js — test alternative ways to COMBINE the oriented detectors
 * (rc, PT, LongString), on Study 1 (induced) and the 6 bogus external datasets.
 * All scores are oriented robust-z (high = careless); AUC vs the GT. */
const R=require("./site/reread.js"); const W=R.WEIGHTS; const fs=require("fs");
const auc=(s,y)=>{let po=[],ne=[];for(let i=0;i<y.length;i++)(y[i]?po:ne).push(s[i]);
  if(!po.length||!ne.length)return NaN;let c=0;for(const a of po)for(const b of ne)c+=a>b?1:a===b?0.5:0;return c/(po.length*ne.length);};
const D="../Dataset/gt_benchmark_candidates";
const SETS=[
  ["Study 1 (induced)", "_study_matrix.csv", "_study_labels.csv", "careless"],
  ["Kay S1", D+"/kay_idris_idria/kay_matrix_s1.csv", "bogus_bench/kay_s1_gt.csv", null],
  ["Kay S2", D+"/kay_idris_idria/kay_matrix_s2.csv", "bogus_bench/kay_s2_gt.csv", null],
  ["Kay S5", D+"/kay_idris_idria/kay_matrix_s5.csv", "bogus_bench/kay_s5_gt.csv", null],
  ["warning", D+"/warning_ipipneo300/_matrix.csv", "bogus_bench/warning_gt.csv", null],
  ["smarvus", D+"/smarvus/smarvus_matrix.csv", "bogus_bench/smarvus_gt.csv", null],
  ["Krause", D+"/krause_ier_youth/krause_matrix.csv", "bogus_bench/krause_gt.csv", null],
];
// combination rules on oriented z (rc, pt, ls) and gate g
const RULES={
  "shipped (2.68rc+g(1.24ls+1.44pt))": (rc,pt,ls,g)=>W.b0+W.rr*rc+g*(W.longstring*ls+W.person_total*pt),
  "linear eq (rc+pt)":                 (rc,pt,ls,g)=>rc+pt,
  "linear eq (rc+pt+ls)":              (rc,pt,ls,g)=>rc+pt+ls,
  "MAX(rc,pt,ls)":                     (rc,pt,ls,g)=>Math.max(rc,pt,ls),
  "MAX(rc,pt)":                        (rc,pt,ls,g)=>Math.max(rc,pt),
  "user: MAX((rc+pt)/2, ls)":          (rc,pt,ls,g)=>Math.max((rc+pt)/2,ls),
  "MAX(rc,pt,gated ls)":               (rc,pt,ls,g)=>Math.max(rc,pt,g*ls),
};
const names=Object.keys(RULES);
const rows=[["dataset",...names]];
const acc={}; names.forEach(k=>acc[k]=[]);
console.log("dataset      "+names.map(n=>n.slice(0,10).padEnd(11)).join(""));
for(const [name,mp,gp,ycol] of SETS){
  const P=R.parseCSV(fs.readFileSync(mp,"utf8")),J=P.header.length;
  const mat=P.rows.map(r=>r.map(v=>{const x=Number(v);return Number.isFinite(x)?x:NaN;}));
  let gt;
  if(ycol){const L=R.parseCSV(fs.readFileSync(gp,"utf8"));const yi=L.header.indexOf(ycol);gt=L.rows.map(r=>Number(r[yi]));}
  else gt=R.parseCSV(fs.readFileSync(gp,"utf8")).rows.map(r=>Number(r[0]));
  const n=Math.min(mat.length,gt.length);gt=gt.slice(0,n);
  const M=new Float64Array(n*J);for(let i=0;i<n;i++)for(let j=0;j<J;j++)M[i*J+j]=mat[i][j];
  const ef=R.ensembleFeatures(M,n,J,{iterations:300,seed:1}),O=ef.oriented,g=Math.max(0,Math.min(1,(ef.profileInfo-0.01)/0.05));
  const line=[name];
  for(const k of names){const s=new Array(n);for(let i=0;i<n;i++)s[i]=RULES[k](O.rr[i],O.person_total[i],O.longstring[i],g);
    const a=auc(s,gt);acc[k].push(a);line.push(a.toFixed(3));}
  rows.push(line);
  console.log(name.padEnd(12)+line.slice(1).map(v=>v.padEnd(11)).join(""));
}
const mean=a=>a.reduce((s,v)=>s+v,0)/a.length;
console.log("\n--- means ---");
console.log("Study1       "+names.map(k=>acc[k][0].toFixed(3).padEnd(11)).join(""));
console.log("6 bogus mean "+names.map(k=>mean(acc[k].slice(1)).toFixed(3).padEnd(11)).join(""));
fs.writeFileSync("exp_combine.csv",rows.map(r=>r.join(",")).join("\n"));
