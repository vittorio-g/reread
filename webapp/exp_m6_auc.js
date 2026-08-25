/* exp_m6_auc.js — broader AUC parity check: shipped |cor(A,B)| vs swap-invariant ICC,
 * baseline z (fixed null), across the GT benchmark datasets. Confirms the ICC preserves
 * detection everywhere (the swap-invariance was already shown in exp_m6_swap.js). */
const R = require("./site/rerere.js"); const fs = require("fs");
const auc = (s, y) => { let po = [], ne = []; for (let i = 0; i < y.length; i++)(y[i] ? po : ne).push(s[i]);
  if (!po.length || !ne.length) return NaN; let c = 0; for (const a of po) for (const b of ne) c += a > b ? 1 : a === b ? 0.5 : 0; return c / (po.length * ne.length); };
function rank(a){const idx=a.map((v,i)=>[v,i]).sort((p,q)=>p[0]-q[0]);const r=new Array(a.length);for(let k=0;k<idx.length;){let j=k;while(j<idx.length&&idx[j][0]===idx[k][0])j++;const av=(k+j-1)/2+1;for(let t=k;t<j;t++)r[idx[t][1]]=av;k=j;}return r;}
const spearman=(x,y)=>{const n=x.length,rx=rank(x),ry=rank(y);let sx=0,sy=0,sxx=0,syy=0,sxy=0;for(let i=0;i<n;i++){sx+=rx[i];sy+=ry[i];}sx/=n;sy/=n;for(let i=0;i<n;i++){const dx=rx[i]-sx,dy=ry[i]-sy;sxx+=dx*dx;syy+=dy*dy;sxy+=dx*dy;}return sxy/Math.sqrt(sxx*syy);};
function prepare(M,n,J){const max=new Float64Array(J),min=new Float64Array(J),med=new Float64Array(J);
  for(let j=0;j<J;j++){const c=[];for(let i=0;i<n;i++){const v=M[i*J+j];if(Number.isFinite(v))c.push(v);}c.sort((a,b)=>a-b);max[j]=c[c.length-1];min[j]=c[0];med[j]=c[Math.floor(c.length/2)];}
  const P=new Float64Array(n*J);for(let i=0;i<n;i++)for(let j=0;j<J;j++){let v=M[i*J+j];if(!Number.isFinite(v))v=med[j];P[i*J+j]=v/max[j];}
  const minP=new Float64Array(J);for(let j=0;j<J;j++)minP[j]=min[j]/max[j];return{P,minP};}
function corPairs(P,n,J){const mean=new Float64Array(J),sdv=new Float64Array(J);
  for(let j=0;j<J;j++){let s=0;for(let i=0;i<n;i++)s+=P[i*J+j];mean[j]=s/n;let q=0;for(let i=0;i<n;i++){const d=P[i*J+j]-mean[j];q+=d*d;}sdv[j]=Math.sqrt(q/(n-1))||0;}
  const nP=J*(J-1)/2;const pa=new Int32Array(nP),pb=new Int32Array(nP),pr=new Float64Array(nP);let p=0;
  const Z=new Float64Array(n*J);for(let i=0;i<n;i++)for(let j=0;j<J;j++)Z[i*J+j]=sdv[j]>0?(P[i*J+j]-mean[j])/sdv[j]:0;
  for(let a=0;a<J;a++)for(let b=a+1;b<J;b++,p++){let s=0;for(let i=0;i<n;i++)s+=Z[i*J+a]*Z[i*J+b];pa[p]=a;pb[p]=b;pr[p]=(sdv[a]>0&&sdv[b]>0)?s/(n-1):0;}
  return{pa,pb,pr,nPairs:nP};}
function cohCor(P,minP,n,J,pa,pb,pr,idx){const k=idx.length,out=new Float64Array(n),A=new Float64Array(k),B=new Float64Array(k);
  for(let i=0;i<n;i++){for(let t=0;t<k;t++){const p=idx[t];A[t]=P[i*J+pa[p]];const bv=P[i*J+pb[p]];B[t]=pr[p]<0?(1+minP[pb[p]])-bv:bv;}
    let ma=0,mb=0;for(let t=0;t<k;t++){ma+=A[t];mb+=B[t];}ma/=k;mb/=k;let sab=0,sa=0,sb=0;for(let t=0;t<k;t++){const da=A[t]-ma,db=B[t]-mb;sab+=da*db;sa+=da*da;sb+=db*db;}
    out[i]=(sa>0&&sb>0)?Math.abs(sab/Math.sqrt(sa*sb)):0;}return out;}
function cohICC(P,minP,n,J,pa,pb,pr,idx){const k=idx.length,out=new Float64Array(n);
  for(let i=0;i<n;i++){let sAB=0,sSq=0,sSum=0;for(let t=0;t<k;t++){const p=idx[t];const a=P[i*J+pa[p]];const bv=P[i*J+pb[p]];const b=pr[p]<0?(1+minP[pb[p]])-bv:bv;sAB+=a*b;sSq+=a*a+b*b;sSum+=a+b;}
    const m=sSum/(2*k);const num=sAB/k-m*m;const den=sSq/(2*k)-m*m;out[i]=den>1e-12?Math.abs(num/den):0;}return out;}
function topIdx(pr,nP,k){const o=Array.from({length:nP},(_,i)=>i).sort((x,y)=>Math.abs(pr[y])-Math.abs(pr[x]));return Int32Array.from(o.slice(0,k));}
function nullZ(fn,n,k,nP,iters,seed){const rng=R._rng(seed);const sum=new Float64Array(n),sq=new Float64Array(n);
  for(let b=0;b<iters;b++){const pool=Int32Array.from({length:nP},(_,i)=>i);const s=new Int32Array(k);for(let t=0;t<k;t++){const r=t+Math.floor(rng()*(nP-t));const tm=pool[t];pool[t]=pool[r];pool[r]=tm;s[t]=pool[t];}const c=fn(s);for(let i=0;i<n;i++){sum[i]+=c[i];sq[i]+=c[i]*c[i];}}
  const mu=new Float64Array(n),sd=new Float64Array(n);for(let i=0;i<n;i++){mu[i]=sum[i]/iters;sd[i]=Math.sqrt(Math.max(0,sq[i]/iters-mu[i]*mu[i]));}return{mu,sd};}
const readM=p=>{const P=R.parseCSV(fs.readFileSync(p,"utf8"));return{rows:P.rows.map(r=>r.map(Number)),J:P.header.length};};
const CAP=2500, ITERS=500;
const D="../Dataset/gt_benchmark_candidates";
const SETS=[
  ["Study1","_study_matrix.csv","_study_labels.csv","careless"],
  ["Kay_S1",D+"/kay_idris_idria/kay_matrix_s1.csv","bogus_bench/kay_s1_gt.csv","gt"],
  ["Kay_S2",D+"/kay_idris_idria/kay_matrix_s2.csv","bogus_bench/kay_s2_gt.csv","gt"],
  ["Kay_S5",D+"/kay_idris_idria/kay_matrix_s5.csv","bogus_bench/kay_s5_gt.csv","gt"],
  ["warning",D+"/warning_ipipneo300/_matrix.csv","bogus_bench/warning_gt.csv","gt"],
  ["smarvus",D+"/smarvus/smarvus_matrix.csv","bogus_bench/smarvus_gt.csv","gt"],
  ["krause",D+"/krause_ier_youth/krause_matrix.csv","bogus_bench/krause_gt.csv","gt"],
  ["duckworth",D+"/duckworth_grit_vcl/duckworth_gt.csv","",""], // placeholder skip (overclaiming control)
];
console.log("dataset      n     J    AUC_cor  AUC_ICC   dAUC   rank_agree");
let sc=0,si=0,cnt=0;
for(const [name,mp,gp,yc] of SETS){
  if(!gp) continue;
  const M=readM(mp); let n=M.rows.length, J=M.J;
  const L=R.parseCSV(fs.readFileSync(gp,"utf8")); const yi=L.header.indexOf(yc); let y=L.rows.slice(0,n).map(r=>Number(r[yi]));
  let rows=M.rows;
  if(n>CAP){ const rng=R._rng(42); const id=rows.map((_,i)=>i); for(let i=id.length-1;i>0;i--){const j=Math.floor(rng()*(i+1));const t=id[i];id[i]=id[j];id[j]=t;} const keep=id.slice(0,CAP); rows=keep.map(i=>rows[i]); y=keep.map(i=>y[i]); n=CAP; }
  const flat=new Float64Array(n*J);for(let i=0;i<n;i++)for(let j=0;j<J;j++)flat[i*J+j]=rows[i][j];
  const {P,minP}=prepare(flat,n,J); const {pa,pb,pr,nPairs}=corPairs(P,n,J);
  const k=Math.min(Math.max(15,Math.round(0.03*nPairs)),nPairs); const idx=topIdx(pr,nPairs,k);
  const nc=nullZ(s=>cohCor(P,minP,n,J,pa,pb,pr,s),n,k,nPairs,ITERS,7); const ni=nullZ(s=>cohICC(P,minP,n,J,pa,pb,pr,s),n,k,nPairs,ITERS,7);
  const cc=cohCor(P,minP,n,J,pa,pb,pr,idx), ci=cohICC(P,minP,n,J,pa,pb,pr,idx);
  const zc=Array.from({length:n},(_,i)=>nc.sd[i]>1e-12?(cc[i]-nc.mu[i])/nc.sd[i]:0);
  const zi=Array.from({length:n},(_,i)=>ni.sd[i]>1e-12?(ci[i]-ni.mu[i])/ni.sd[i]:0);
  const aC=auc(zc.map(v=>-v),y), aI=auc(zi.map(v=>-v),y);
  console.log(`${name.padEnd(11)}${String(n).padStart(5)}${String(J).padStart(5)}   ${aC.toFixed(3)}    ${aI.toFixed(3)}   ${(aI-aC>=0?"+":"")}${(aI-aC).toFixed(3)}   ${spearman(zc,zi).toFixed(3)}`);
  sc+=aC; si+=aI; cnt++;
}
console.log(`\nmean AUC   cor ${(sc/cnt).toFixed(3)}   ICC ${(si/cnt).toFixed(3)}   dAUC ${((si-sc)/cnt>=0?"+":"")}${((si-sc)/cnt).toFixed(3)}`);
