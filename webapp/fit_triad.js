const fs=require("fs");const R=require("./site/rerere.js");
const M=R.parseCSV(fs.readFileSync("_study_matrix.csv","utf8"));
const n=M.rows.length,J=M.header.length,mat=new Float64Array(n*J);
for(let i=0;i<n;i++)for(let j=0;j<J;j++)mat[i*J+j]=R.toNum(M.rows[i][j]);
const y=R.parseCSV(fs.readFileSync("_study_labels.csv","utf8")).rows.map(r=>Number(r[0]));
const O=R.ensembleFeatures(mat,n,J,{seed:1,iterations:300}).oriented;
const feats=["rr","longstring","person_total"];
const X=Array.from({length:n},(_,i)=>[1,...feats.map(f=>O[f][i])]),P=feats.length+1;
function solve(A,b){const m=b.length,Aug=A.map((r,i)=>r.concat(b[i]));
 for(let c=0;c<m;c++){let pv=c;for(let r=c+1;r<m;r++)if(Math.abs(Aug[r][c])>Math.abs(Aug[pv][c]))pv=r;[Aug[c],Aug[pv]]=[Aug[pv],Aug[c]];const d=Aug[c][c]||1e-12;
   for(let r=0;r<m;r++)if(r!==c){const f=Aug[r][c]/d;for(let k=c;k<=m;k++)Aug[r][k]-=f*Aug[c][k];}}return Aug.map((r,i)=>r[m]/(r[i]||1e-12));}
function fit(rows,yy){let b=new Array(P).fill(0);for(let it=0;it<50;it++){const g=new Array(P).fill(0),H=Array.from({length:P},()=>new Array(P).fill(0));
 for(let i=0;i<rows.length;i++){let e=0;for(let k=0;k<P;k++)e+=b[k]*rows[i][k];const p=1/(1+Math.exp(-e)),w=Math.max(p*(1-p),1e-6);
   for(let k=0;k<P;k++){g[k]+=(yy[i]-p)*rows[i][k];for(let l=0;l<P;l++)H[k][l]+=w*rows[i][k]*rows[i][l];}}
 for(let k=1;k<P;k++){g[k]-=1.0*b[k];H[k][k]+=1.0;}const s=solve(H,g);let mx=0;for(let k=0;k<P;k++){b[k]+=s[k];mx=Math.max(mx,Math.abs(s[k]));}if(mx<1e-8)break;}return b;}
const b=fit(X,y);
const W={b0:+b[0].toFixed(4),rr:+b[1].toFixed(4),irv:0,longstring:+b[2].toFixed(4),d2:0,person_total:+b[3].toFixed(4)};
fs.writeFileSync("ensemble_weights.json",JSON.stringify(W,null,2));
console.log(JSON.stringify(W));
