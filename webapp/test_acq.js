/* test_acq.js — does adding an ACQUIESCENCE feature help separate response-style distractors
 * (acquiescent yea-sayers we want to KEEP) from real junk careless? Acquiescence index = the
 * person's raw mean response (on a mixed-keyed scale a yea-sayer scores high on both regular and
 * reverse items -> high raw mean). We refit the logistic ensemble WITH vs WITHOUT acq (5-fold CV)
 * on simulated data (with acquiescent distractors) and report weighted-MCC + distractor handling. */
const R=require("./site/reread.js"); const W=R.WEIGHTS; const fs=require("fs"); const rng=R._rng(9911);
const src=fs.readFileSync("sim_v2.js","utf8");
eval(src.slice(src.indexOf("function gauss"), src.indexOf("function scoreDataset")));  // helpers
// ---- assemble one big labelled sim sample (120 items, 15% careless + 15% distractors) ----
function features(out,J){
  const n=out.length,mat=new Float64Array(n*J);for(let i=0;i<n;i++)for(let j=0;j<J;j++)mat[i*J+j]=out[i][j];
  const ef=R.ensembleFeatures(mat,n,J,{seed:1,iterations:150}),O=ef.oriented;
  const g=Math.max(0,Math.min(1,(ef.profileInfo-0.01)/0.05));
  // acquiescence = robust-z of the person's raw mean response (high = yea-sayer)
  const rawMean=[];for(let i=0;i<n;i++){let s=0;for(let j=0;j<J;j++)s+=out[i][j];rawMean.push(s/J);}
  const acq=R.robustZ(rawMean,n);
  const X=[];for(let i=0;i<n;i++)X.push([1,O.rr[i],g*O.longstring[i],g*O.person_total[i],acq[i]]);  // last col = acq
  return X;
}
function auc(s,y){let po=[],ne=[];for(let i=0;i<y.length;i++)(y[i]?po:ne).push(s[i]);let c=0;for(const a of po)for(const b of ne)c+=a>b?1:a===b?0.5:0;return c/(po.length*ne.length);}
function oracleWMCC(p,y,w){const s=[...p].sort((a,b)=>a-b);let best=0,bt=0;const step=Math.max(1,Math.floor(s.length/100));
  for(let t=0;t<s.length;t+=step){const thr=s[t];const f=p.map(v=>v>=thr);const m=wmcc(f,y,w);if(m>best){best=m;bt=thr;}}return {mcc:best,thr:bt};}
function solve(A,b){const m=b.length,Au=A.map((r,i)=>r.concat(b[i]));for(let c=0;c<m;c++){let p=c;for(let r=c+1;r<m;r++)if(Math.abs(Au[r][c])>Math.abs(Au[p][c]))p=r;[Au[c],Au[p]]=[Au[p],Au[c]];const d=Au[c][c]||1e-9;for(let r=0;r<m;r++)if(r!==c){const f=Au[r][c]/d;for(let k=c;k<=m;k++)Au[r][k]-=f*Au[c][k];}}return Au.map((r,i)=>r[m]/(r[i]||1e-9));}
function fit(rows,y,P){let b=new Array(P).fill(0);for(let it=0;it<80;it++){const g=new Array(P).fill(0),H=Array.from({length:P},()=>new Array(P).fill(0));
  for(let i=0;i<rows.length;i++){let e=0;for(let k=0;k<P;k++)e+=b[k]*rows[i][k];const p=1/(1+Math.exp(-e)),w=Math.max(p*(1-p),1e-6);
    for(let k=0;k<P;k++){g[k]+=(y[i]-p)*rows[i][k];for(let l=0;l<P;l++)H[k][l]+=w*rows[i][k]*rows[i][l];}}
  for(let k=1;k<P;k++){g[k]-=0.5*b[k];H[k][k]+=0.5;}const s=solve(H,g);let mx=0;for(let k=0;k<P;k++){b[k]+=s[k];mx=Math.max(mx,Math.abs(s[k]));}if(mx<1e-8)break;}return b;}
const pred=(b,row)=>{let e=0;for(let k=0;k<b.length;k++)e+=b[k]*row[k];return 1/(1+Math.exp(-e));};

// build a pooled dataset (several reps -> one big matrix scored per-rep, features stacked)
let X=[],y=[],dist=[],sev=[];
for(let rep=0;rep<12;rep++){
  const {rows}=genClean(20,6,600);const {out,labels,dist:d,sev:s}=inject2(rows,120,0.15,0.15);
  const F=features(out,120);
  for(let i=0;i<out.length;i++){X.push(F[i]);y.push(labels[i]);dist.push(d[i]);sev.push(s[i]);}
}
const n=X.length,w=sev.map(uWeight);
console.log("pooled sim: "+n+" respondents ("+y.reduce((a,b)=>a+b,0)+" careless, "+dist.filter(Boolean).length+" distractors)");
// 5-fold CV, base (rr+ls+pt) vs +acq
function cv(cols){const P=cols.length+1,oof=new Array(n);
  for(let f=0;f<5;f++){const tr=[],ytr=[],te=[];for(let i=0;i<n;i++)(i%5===f?te:tr).push(i);
    const Xtr=tr.map(i=>[1,...cols.map(c=>X[i][c])]),Ytr=tr.map(i=>y[i]);const b=fit(Xtr,Ytr,P);
    for(const i of te)oof[i]=pred(b,[1,...cols.map(c=>X[i][c])]);}
  return oof;}
function report(name,cols){const p=cv(cols);const A=auc(p,y);const ow=oracleWMCC(p,y,w);
  const flag=p.map(v=>v>=ow.thr);
  let dtot=0,dflag=0;for(let i=0;i<n;i++)if(dist[i]){dtot++;if(flag[i])dflag++;}
  const meanCare=y.reduce((a,b,i)=>a+(b?p[i]:0),0)/y.reduce((a,b)=>a+b,0);
  const meanDist=dist.reduce((a,b,i)=>a+(b?p[i]:0),0)/dist.filter(Boolean).length;
  console.log(name.padEnd(16)+" AUC "+A.toFixed(3)+"  wMCC(oracle) "+ow.mcc.toFixed(3)+"  distractor-FPR "+(dflag/dtot).toFixed(3)+
    "   | mean score: careless "+meanCare.toFixed(2)+" vs distractor "+meanDist.toFixed(2));
  return {A,wmcc:ow.mcc,dfpr:dflag/dtot};}
console.log("\n5-fold CV, weighted-MCC (extremes heavy):");
const base=report("rr+LS+PT",[1,2,3]);          // cols 1,2,3 = rr, g*ls, g*pt
const wacq=report("rr+LS+PT+ACQ",[1,2,3,4]);     // + col 4 = acq
console.log("\nΔ from acquiescence feature: wMCC "+(wacq.wmcc-base.wmcc>=0?"+":"")+(wacq.wmcc-base.wmcc).toFixed(3)+
  ",  distractor-FPR "+(wacq.dfpr-base.dfpr>=0?"+":"")+(wacq.dfpr-base.dfpr).toFixed(3)+"  (lower FPR = better)");
fs.writeFileSync("test_acq.csv","model,auc,wmcc,dfpr\nbase,"+base.A.toFixed(4)+","+base.wmcc.toFixed(4)+","+base.dfpr.toFixed(4)+
  "\nacq,"+wacq.A.toFixed(4)+","+wacq.wmcc.toFixed(4)+","+wacq.dfpr.toFixed(4)+"\n");
console.log("wrote test_acq.csv");
