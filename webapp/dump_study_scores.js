/* dump_study_scores.js — per-respondent shipped rr (oriented high=careless) and ensemble log-odds
 * (high=careless) for the Study-1 matrix, for the competitive benchmark. */
const R=require("./site/reread.js"); const W=R.WEIGHTS; const fs=require("fs");
const Ms=R.parseCSV(fs.readFileSync("_study_matrix.csv","utf8")), J=Ms.header.length;
const rows=Ms.rows.map(r=>r.map(v=>{const x=Number(v);return Number.isFinite(x)?x:NaN;}));
const n=rows.length, mat=new Float64Array(n*J);
for(let i=0;i<n;i++)for(let j=0;j<J;j++)mat[i*J+j]=rows[i][j];
const ef=R.ensembleFeatures(mat,n,J,{seed:1,iterations:400}), O=ef.oriented;
const g=Math.max(0,Math.min(1,(ef.profileInfo-0.01)/0.05));
const out=[["rr","longstring","person_total","irv","d2","eta"]];
for(let i=0;i<n;i++){
  const eta=W.b0+W.rr*O.rr[i]+g*(W.longstring*O.longstring[i]+W.person_total*O.person_total[i]);
  out.push([O.rr[i],O.longstring[i],O.person_total[i],O.irv[i],O.d2[i],eta].map(x=>Number(x).toFixed(6)));
}
fs.writeFileSync("_study_scores.csv",out.map(r=>r.join(",")).join("\n"));
console.log("wrote _study_scores.csv  (n="+n+", J="+J+", gate g="+g.toFixed(3)+", profileInfo="+ef.profileInfo.toFixed(4)+")");
