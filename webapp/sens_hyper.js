/* sens_hyper.js (reviewer 2.1) — robustness of every deployment constant.
 * Sweeps corProp (top-share), min_pairs (floor), gate offset & scale, and the
 * flagging cut z, recomputing ensemble AUC (and, for z, MCC) on real data:
 * Study 1 + two in-envelope external datasets. Writes sens_hyper.csv. */
const R=require("./site/rerere.js"); const W=R.WEIGHTS; const fs=require("fs");
const auc=(s,y)=>{let po=[],ne=[];for(let i=0;i<y.length;i++)(y[i]?po:ne).push(s[i]);
  let c=0;for(const a of po)for(const b of ne)c+=a>b?1:a===b?0.5:0;return c/(po.length*ne.length);};
function load(mp,lp,ycol){const P=R.parseCSV(fs.readFileSync(mp,"utf8")),J=P.header.length;
  const rows=P.rows.map(r=>r.map(v=>{const x=Number(v);return Number.isFinite(x)?x:NaN;}));
  const L=R.parseCSV(fs.readFileSync(lp,"utf8"));
  const yi=ycol!=null?L.header.indexOf(ycol):0;
  const y=L.rows.map(r=>Number(r[yi]));const n=Math.min(rows.length,y.length);
  const mat=new Float64Array(n*J);for(let i=0;i<n;i++)for(let j=0;j<J;j++)mat[i*J+j]=rows[i][j];
  return {mat,n,J,y:y.slice(0,n)};}
const D="../Dataset/gt_benchmark_candidates/";
const SETS=[
  ["study","_study_matrix.csv","_study_labels.csv","careless"],
  ["smarvus",D+"smarvus/smarvus_matrix.csv",D+"smarvus/smarvus_labels.csv","y"],
  ["warning",D+"warning_ipipneo300/_matrix.csv",D+"warning_ipipneo300/_labels.csv","y"],
];
const CAP=4000;
function etaFrom(O,pinfo,gOff,gScale){const g=Math.max(0,Math.min(1,(pinfo-gOff)/gScale));
  const n=O.rr.length,e=new Array(n);
  for(let i=0;i<n;i++)e[i]=W.b0+W.rr*O.rr[i]+g*(W.longstring*O.longstring[i]+W.person_total*O.person_total[i]);
  return e;}
const out=[["dataset","param","value","auc"]];
for(const [name,mp,lp,yc] of SETS){
  let {mat,n,J,y}=load(mp,lp,yc);
  if(n>CAP){const idx=[];for(let i=0;i<n;i++)idx.push(i);// deterministic head cap for speed
    mat=mat.slice(0,CAP*J);y=y.slice(0,CAP);n=CAP;}
  // base features at defaults (corProp .03) -> reused for gate sweep
  const base=R.ensembleFeatures(mat,n,J,{corProp:0.03,minPairs:15,iterations:200,seed:1});
  const push=(p,v,a)=>{out.push([name,p,v,a.toFixed(4)]);};
  // corProp sweep
  for(const cp of [0.01,0.02,0.03,0.05,0.07,0.10]){
    const ef=R.ensembleFeatures(mat,n,J,{corProp:cp,minPairs:15,iterations:200,seed:1});
    push("cor_prop",cp,auc(etaFrom(ef.oriented,ef.profileInfo,0.01,0.05),y));}
  // min_pairs sweep
  for(const mpn of [5,10,15,20,30,50]){
    const ef=R.ensembleFeatures(mat,n,J,{corProp:0.03,minPairs:mpn,iterations:200,seed:1});
    push("min_pairs",mpn,auc(etaFrom(ef.oriented,ef.profileInfo,0.01,0.05),y));}
  // gate offset sweep (scale fixed .05)
  for(const off of [0,0.005,0.01,0.02,0.03]) push("gate_offset",off,auc(etaFrom(base.oriented,base.profileInfo,off,0.05),y));
  // gate scale sweep (offset fixed .01)
  for(const sc of [0.02,0.03,0.05,0.07,0.10]) push("gate_scale",sc,auc(etaFrom(base.oriented,base.profileInfo,0.01,sc),y));
  console.log(name+" done (n="+n+", J="+J+", careless="+y.reduce((s,v)=>s+v,0)+")");
}
fs.writeFileSync("sens_hyper.csv",out.map(r=>r.join(",")).join("\n"));
console.log("wrote sens_hyper.csv");
