/* bogus_score.js — score each >=5-bogus dataset with the shipped engine and
 * compute AUC of rc / ensemble / Person-Total against the bogus/instructed GT.
 * Bootstrap 95% CI on the rc AUC. */
const R=require("./site/rerere.js"); const W=R.WEIGHTS; const fs=require("fs"); const rng=R._rng(7);
const D="../Dataset/gt_benchmark_candidates";
const SETS=[
  ["Kay S1 (IDRIS)", D+"/kay_idris_idria/kay_matrix_s1.csv", "bogus_bench/kay_s1_gt.csv", 250],
  ["Kay S2 (IDRIS)", D+"/kay_idris_idria/kay_matrix_s2.csv", "bogus_bench/kay_s2_gt.csv", 363],
  ["Kay S5 (IDRIS)", D+"/kay_idris_idria/kay_matrix_s5.csv", "bogus_bench/kay_s5_gt.csv", 62],
  ["warning (6 bogus)", D+"/warning_ipipneo300/_matrix.csv", "bogus_bench/warning_gt.csv", 301],
  ["smarvus (6 instr.)", D+"/smarvus/smarvus_matrix.csv", "bogus_bench/smarvus_gt.csv", 141],
  ["Krause (6 checks)", D+"/krause_ier_youth/krause_matrix.csv", "bogus_bench/krause_gt.csv", 49],
];
const auc=(s,y)=>{let po=[],ne=[];for(let i=0;i<y.length;i++)(y[i]?po:ne).push(s[i]);
  if(!po.length||!ne.length)return NaN;let c=0;for(const a of po)for(const b of ne)c+=a>b?1:a===b?0.5:0;return c/(po.length*ne.length);};
function bootCI(s,y,B){const n=y.length,out=[];for(let b=0;b<B;b++){const ss=[],yy=[];for(let i=0;i<n;i++){const j=Math.floor(rng()*n);ss.push(s[j]);yy.push(y[j]);}out.push(auc(ss,yy));}
  out.sort((a,b)=>a-b);return [out[Math.floor(0.025*B)],out[Math.floor(0.975*B)]];}
console.log("dataset             items    N   rate   rc  [95% CI]        ens    PT");
const rows=[["dataset","items","N","rate","rc","rc_lo","rc_hi","ens","pt"]];
for(const [name,mp,gp,items] of SETS){
  const P=R.parseCSV(fs.readFileSync(mp,"utf8")),J=P.header.length;
  let mat=P.rows.map(r=>r.map(v=>{const x=Number(v);return Number.isFinite(x)?x:NaN;}));
  let gt=R.parseCSV(fs.readFileSync(gp,"utf8")).rows.map(r=>Number(r[0]));
  const n=Math.min(mat.length,gt.length); mat=mat.slice(0,n); gt=gt.slice(0,n);
  const M=new Float64Array(n*J);for(let i=0;i<n;i++)for(let j=0;j<J;j++)M[i*J+j]=mat[i][j];
  const ef=R.ensembleFeatures(M,n,J,{iterations:300,seed:1}),O=ef.oriented,g=Math.max(0,Math.min(1,(ef.profileInfo-0.01)/0.05));
  const eta=new Array(n);for(let i=0;i<n;i++)eta[i]=W.b0+W.rr*O.rr[i]+g*(W.longstring*O.longstring[i]+W.person_total*O.person_total[i]);
  const arc=auc(O.rr,gt), aens=auc(eta,gt), apt=auc(O.person_total,gt), ci=bootCI(O.rr,gt,1000);
  const rate=gt.reduce((s,v)=>s+v,0)/n;
  rows.push([name,items,n,(rate*100).toFixed(0)+"%",arc.toFixed(3),ci[0].toFixed(3),ci[1].toFixed(3),aens.toFixed(3),apt.toFixed(3)]);
  console.log(name.padEnd(19)+String(items).padStart(4)+String(n).padStart(6)+(rate*100).toFixed(0).padStart(5)+"%  "+
    arc.toFixed(3)+" ["+ci[0].toFixed(3)+","+ci[1].toFixed(3)+"]   "+aens.toFixed(3)+"  "+apt.toFixed(3));
}
fs.writeFileSync("bogus_bench/results.csv",rows.map(r=>r.join(",")).join("\n"));
// mean rc weighted equally
const rc=rows.slice(1).map(r=>Number(r[4]));
console.log("\nmean rc AUC (6 datasets):", (rc.reduce((s,v)=>s+v,0)/rc.length).toFixed(3));
