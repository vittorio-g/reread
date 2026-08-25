const fs=require("fs");const R=require("./site/rerere.js");
const M=R.parseCSV(fs.readFileSync("_study_matrix.csv","utf8"));
const n=M.rows.length,J=M.header.length,mat=new Float64Array(n*J);
for(let i=0;i<n;i++)for(let j=0;j<J;j++)mat[i*J+j]=R.toNum(M.rows[i][j]);
const y=R.parseCSV(fs.readFileSync("_study_labels.csv","utf8")).rows.map(r=>Number(r[0]));
console.log("n="+n+"  careless="+y.reduce((s,v)=>s+v,0)+"  careful="+y.filter(v=>v===0).length);
const ef=R.ensembleFeatures(mat,n,J,{seed:1,iterations:400}); const O=ef.oriented;
const raw={rr:ef.base.z,irv:ef.base.irv,longstring:ef.base.longstring,d2:ef.d2,person_total:ef.base.personTotal};
function auc(sc,hiCare){const s=hiCare?sc:sc.map(v=>-v);const idx=[...s.keys()].sort((a,b)=>s[a]-s[b]);
  const rk=[];idx.forEach((id,r)=>rk[id]=r+1);let n1=0,n0=0,sr=0;
  for(let i=0;i<s.length;i++){if(y[i]===1){n1++;sr+=rk[i];}else n0++;}return (sr-n1*(n1+1)/2)/(n1*n0);}
console.log("\n=== per-feature standalone AUC (higher=careless) ===");
for(const[k,d] of [["person_total",false],["rr",false],["longstring",true],["d2",true],["irv",true]])
  console.log("  "+k.padEnd(13),auc(raw[k],d).toFixed(3));
// logistic on oriented features
function solve(A,b){const m=b.length,Au=A.map((r,i)=>r.concat(b[i]));
 for(let c=0;c<m;c++){let p=c;for(let r=c+1;r<m;r++)if(Math.abs(Au[r][c])>Math.abs(Au[p][c]))p=r;[Au[c],Au[p]]=[Au[p],Au[c]];const d=Au[c][c]||1e-12;
   for(let r=0;r<m;r++)if(r!==c){const f=Au[r][c]/d;for(let k=c;k<=m;k++)Au[r][k]-=f*Au[c][k];}}return Au.map((r,i)=>r[m]/(r[i]||1e-12));}
function fit(rows,yy,P){let b=new Array(P).fill(0);for(let it=0;it<60;it++){const g=new Array(P).fill(0),H=Array.from({length:P},()=>new Array(P).fill(0));
  for(let i=0;i<rows.length;i++){let e=0;for(let k=0;k<P;k++)e+=b[k]*rows[i][k];const p=1/(1+Math.exp(-e)),w=Math.max(p*(1-p),1e-6);
    for(let k=0;k<P;k++){g[k]+=(yy[i]-p)*rows[i][k];for(let l=0;l<P;l++)H[k][l]+=w*rows[i][k]*rows[i][l];}}
  for(let k=1;k<P;k++){g[k]-=1.0*b[k];H[k][k]+=1.0;}const s=solve(H,g);let mx=0;for(let k=0;k<P;k++){b[k]+=s[k];mx=Math.max(mx,Math.abs(s[k]));}if(mx<1e-9)break;}return b;}
const pred=(b,row)=>{let e=0;for(let k=0;k<b.length;k++)e+=b[k]*row[k];return 1/(1+Math.exp(-e));};
function mcc(yy,yh){let tp=0,tn=0,fp=0,fn=0;for(let i=0;i<yy.length;i++){if(yy[i]&&yh[i])tp++;else if(!yy[i]&&!yh[i])tn++;else if(!yy[i]&&yh[i])fp++;else fn++;}
  const d=Math.sqrt((tp+fp)*(tp+fn)*(tn+fp)*(tn+fn));return d?(tp*tn-fp*fn)/d:0;}
function looEval(feats){const cols=feats.map(f=>O[f]),P=feats.length+1;
  const X=Array.from({length:n},(_,i)=>[1,...cols.map(c=>c[i])]);const p=new Array(n);
  for(let i=0;i<n;i++){const rows=[],yy=[];for(let j=0;j<n;j++)if(j!==i){rows.push(X[j]);yy.push(y[j]);}p[i]=pred(fit(rows,yy,P),X[i]);}
  const rate=y.reduce((s,v)=>s+v,0)/n,srt=[...p].sort((a,b)=>a-b),thr=srt[Math.floor((1-rate)*n)];
  return {auc:auc(p,true),mcc:mcc(y,p.map(v=>v>=thr?1:0))};}
console.log("\n=== ablation LOO (new data) ===");
const sets={"full5":["rr","irv","longstring","d2","person_total"],"triad (rr+ls+pt)":["rr","longstring","person_total"],
  "no-rr":["irv","longstring","d2","person_total"],"rr only":["rr"]};
for(const[nm,f] of Object.entries(sets)){const m=looEval(f);console.log("  "+nm.padEnd(20),"AUC",m.auc.toFixed(3),"MCC",m.mcc.toFixed(3));}
// fit shipped triad on full data, write weights
const feats=["rr","longstring","person_total"],P=4,b=fit(Array.from({length:n},(_,i)=>[1,O.rr[i],O.longstring[i],O.person_total[i]]),y,P);
const W={b0:+b[0].toFixed(4),rr:+b[1].toFixed(4),irv:0,longstring:+b[2].toFixed(4),d2:0,person_total:+b[3].toFixed(4)};
fs.writeFileSync("ensemble_weights.json",JSON.stringify(W,null,2));
console.log("\n=== NEW shipped triad weights ==="); console.log(JSON.stringify(W));
