/* eval_broad.js — HEAD-TO-HEAD across every dataset we have, at KNOWN prevalence.
 *
 * External datasets do not come with a known careless rate, so we rebuild the Study-1 design on
 * each of them: GT-negative rows form the attentive base, GT-positive rows are mixed in at
 * controlled rates. The prevalence is then known BY CONSTRUCTION on 16 real data-generating
 * processes (plus Study 1), which is what a fair comparison of prevalence estimators needs.
 *
 * Compared (all computed on the SAME resamples):
 *   solution 1  : one-sided mirror, one-sided 1sigma      (attentive mode estimated from the left)
 *   solution 2  : colperm+corr, rowshufLo ratio, unif+corr (synthetic careless as calibration curve)
 *   baseline    : the shipped two-component mixture pi_eta
 *
 * Caveat kept in view: the GT label is a proxy criterion (bogus items, attention checks,
 * self-report), so "careless" here means "flagged by that criterion" - the mixing is exact, the
 * construct is the criterion's. Output: broad_eval.csv
 */
const R=require("./site/rerere.js"); const W=R.WEIGHTS; const fs=require("fs");
let SEED=20260804; const rng=()=>{SEED=(SEED*1103515245+12345)&0x7fffffff;return SEED/0x7fffffff;};
const mn=a=>a.length?a.reduce((s,v)=>s+v,0)/a.length:NaN;
const med=a=>{const s=[...a].sort((x,y)=>x-y);const h=s.length>>1;return s.length%2?s[h]:(s[h-1]+s[h])/2;};
const clip=p=>Math.max(0,Math.min(1,p));
const phi=z=>Math.exp(-0.5*z*z)/Math.sqrt(2*Math.PI);
function Phi(z){const t=1/(1+0.2316419*Math.abs(z));
  const d=phi(z)*t*(0.319381530+t*(-0.356563782+t*(1.781477937+t*(-1.821255978+t*1.330274429))));
  return z>=0?1-d:d;}
/* ---------- frozen-frame scoring (from nested_loo.js / eval_synthetic.js) ---------- */
function colStats(M,n,J){const max=new Float64Array(J).fill(-Infinity),min=new Float64Array(J).fill(Infinity),median=new Float64Array(J),buf=[];
  for(let j=0;j<J;j++){buf.length=0;
    for(let i=0;i<n;i++){const v=M[i*J+j];if(Number.isFinite(v)){buf.push(v);if(v>max[j])max[j]=v;if(v<min[j])min[j]=v;}}
    buf.sort((a,b)=>a-b); median[j]=buf.length?buf[Math.floor(buf.length/2)]:0;
    if(!Number.isFinite(max[j])){max[j]=1;min[j]=0;} if(max[j]===0)max[j]=1;}
  return {max,min,median};}
function prepareWith(M,J,rowsIdx,cs){const n=rowsIdx.length,P=new Float64Array(n*J);
  for(let a=0;a<n;a++){const i=rowsIdx[a];
    for(let j=0;j<J;j++){let v=M[i*J+j];if(!Number.isFinite(v))v=cs.median[j];P[a*J+j]=v/cs.max[j];}}
  const minP=new Float64Array(J);for(let j=0;j<J;j++)minP[j]=cs.min[j]/cs.max[j];
  return {P,minP};}
function correlationPairs(P,n,J){const mean=new Float64Array(J),sd=new Float64Array(J);
  for(let j=0;j<J;j++){let s=0;for(let i=0;i<n;i++)s+=P[i*J+j];mean[j]=s/n;
    let q=0;for(let i=0;i<n;i++){const d=P[i*J+j]-mean[j];q+=d*d;}sd[j]=Math.sqrt(q/(n-1))||0;}
  const Z=new Float64Array(n*J);
  for(let i=0;i<n;i++)for(let j=0;j<J;j++)Z[i*J+j]=sd[j]>0?(P[i*J+j]-mean[j])/sd[j]:0;
  const nPairs=(J*(J-1))/2,pa=new Int32Array(nPairs),pb=new Int32Array(nPairs),pr=new Float64Array(nPairs);let p=0;
  for(let a=0;a<J;a++)for(let b=a+1;b<J;b++,p++){let s=0;for(let i=0;i<n;i++)s+=Z[i*J+a]*Z[i*J+b];
    pa[p]=a;pb[p]=b;pr[p]=(sd[a]>0&&sd[b]>0)?s/(n-1):0;}
  return {pa,pb,pr,nPairs};}
function indCors(P,minP,n,J,pa,pb,pr,idx){const k=idx.length,out=new Float64Array(n);
  for(let i=0;i<n;i++){let sAB=0,sSq=0,sSum=0;
    for(let t=0;t<k;t++){const p=idx[t];const a=P[i*J+pa[p]];const bj=pb[p],bv=P[i*J+bj];
      const b=pr[p]<0?(1+minP[bj])-bv:bv;sAB+=a*b;sSq+=a*a+b*b;sSum+=a+b;}
    const m=sSum/(2*k);const num=sAB/k-m*m;const den=sSq/(2*k)-m*m;
    out[i]=den>1e-12?Math.abs(num/den):0;}
  return out;}
function sampleIdx(nPairs,k,rand){const pool=new Int32Array(nPairs);for(let i=0;i<nPairs;i++)pool[i]=i;
  const out=new Int32Array(k);
  for(let t=0;t<k;t++){const r=t+Math.floor(rand()*(nPairs-t));const tmp=pool[t];pool[t]=pool[r];pool[r]=tmp;out[t]=pool[t];}
  return out;}
function auxWithRef(P,n,J,cm){const personTotal=new Float64Array(n),longstring=new Int32Array(n);
  let mc=0;for(let j=0;j<J;j++)mc+=cm[j];mc/=J;
  for(let i=0;i<n;i++){let m=0;for(let j=0;j<J;j++)m+=P[i*J+j];m/=J;
    let sab=0,sa=0,sb=0;
    for(let j=0;j<J;j++){const da=P[i*J+j]-m,db=cm[j]-mc;sab+=da*db;sa+=da*da;sb+=db*db;}
    personTotal[i]=(sa>0&&sb>0)?sab/Math.sqrt(sa*sb):0;
    let best=1,run=1;for(let j=1;j<J;j++){if(P[i*J+j]===P[i*J+j-1]){run++;if(run>best)best=run;}else run=1;}
    longstring[i]=best;}
  return {personTotal,longstring};}
function robustParams(x,trainIdx){const v=trainIdx.map(i=>x[i]).sort((a,b)=>a-b),m=v[Math.floor(v.length/2)];
  const dev=trainIdx.map(i=>Math.abs(x[i]-m)).sort((a,b)=>a-b);
  let mad=dev[Math.floor(dev.length/2)]*1.4826;
  if(!(mad>1e-9)){let mm=0;for(const i of trainIdx)mm+=x[i];mm/=trainIdx.length;
    let q=0;for(const i of trainIdx){const d=x[i]-mm;q+=d*d;}
    const sd=Math.sqrt(q/Math.max(trainIdx.length-1,1))||1;return {c:mm,s:sd};}
  return {c:m,s:mad};}
function subMatrix(M,J,idx){const m=new Float64Array(idx.length*J);for(let a=0;a<idx.length;a++)for(let j=0;j<J;j++)m[a*J+j]=M[idx[a]*J+j];return m;}
function subRows(P,J,idx){const m=new Float64Array(idx.length*J);for(let a=0;a<idx.length;a++)for(let j=0;j<J;j++)m[a*J+j]=P[idx[a]*J+j];return m;}
function etaFrozen(M,n,J,trainIdx,iters,seed){
  const allIdx=[...Array(n).keys()];
  const cs=colStats(subMatrix(M,J,trainIdx),trainIdx.length,J);
  const {P:Pall,minP}=prepareWith(M,J,allIdx,cs);
  const Ptr=subRows(Pall,J,trainIdx);
  const {pa,pb,pr,nPairs}=correlationPairs(Ptr,trainIdx.length,J);
  const k=Math.min(Math.max(15,Math.round(0.03*nPairs)),nPairs);
  const order=Array.from({length:nPairs},(_,i)=>i).sort((x,y)=>Math.abs(pr[y])-Math.abs(pr[x]));
  const topIdx=Int32Array.from(order.slice(0,k));
  const coupled=indCors(Pall,minP,n,J,pa,pb,pr,topIdx);
  const rand=R._rng(seed),sum=new Float64Array(n),sq=new Float64Array(n);
  for(let b=0;b<iters;b++){const rc=indCors(Pall,minP,n,J,pa,pb,pr,sampleIdx(nPairs,k,rand));
    for(let i=0;i<n;i++){sum[i]+=rc[i];sq[i]+=rc[i]*rc[i];}}
  const z=new Float64Array(n);let mv=Infinity;
  for(let i=0;i<n;i++){const mu=sum[i]/iters,va=sq[i]/iters-mu*mu,sd=Math.sqrt(Math.max(va,0));
    z[i]=sd>1e-12?(coupled[i]-mu)/sd:NaN; if(Number.isFinite(z[i])&&z[i]<mv)mv=z[i];}
  for(let i=0;i<n;i++)if(!Number.isFinite(z[i]))z[i]=mv;
  const cm=new Float64Array(J);for(const i of trainIdx)for(let j=0;j<J;j++)cm[j]+=Pall[i*J+j];
  for(let j=0;j<J;j++)cm[j]/=trainIdx.length;
  const aux=auxWithRef(Pall,n,J,cm);
  const cmv=new Float64Array(J);
  for(let j=0;j<J;j++){let q=0;for(const i of trainIdx){const d=Pall[i*J+j]-cm[j];q+=d*d;}cmv[j]=q/(trainIdx.length-1);}
  let gm2=0;for(let j=0;j<J;j++)gm2+=cm[j];gm2/=J;
  let vm=0;for(let j=0;j<J;j++)vm+=(cm[j]-gm2)*(cm[j]-gm2);vm/=J;
  let mv2=0;for(let j=0;j<J;j++)mv2+=cmv[j];mv2/=J;
  const g=Math.max(0,Math.min(1,((mv2>0?vm/mv2:0)-0.01)/0.05));
  const pRR=robustParams(z,trainIdx),pLS=robustParams(aux.longstring,trainIdx),pPT=robustParams(aux.personTotal,trainIdx);
  const eta=new Array(n);
  for(let i=0;i<n;i++){const oRR=-((z[i]-pRR.c)/pRR.s),oLS=(aux.longstring[i]-pLS.c)/pLS.s,oPT=-((aux.personTotal[i]-pPT.c)/pPT.s);
    eta[i]=W.b0+W.rr*oRR+g*(W.longstring*oLS+W.person_total*oPT);}
  return eta;}
/* ---------- synthetic generators ---------- */
function genSynthetic(rows,J,kind,nSyn,loHalfIdx){const N=rows.length,out=[];
  const mx=new Array(J).fill(-Infinity),mi=new Array(J).fill(Infinity);
  for(const r of rows)for(let j=0;j<J;j++){const v=r[j];if(Number.isFinite(v)){if(v>mx[j])mx[j]=v;if(v<mi[j])mi[j]=v;}}
  for(let s=0;s<nSyn;s++){
    if(kind==="unif"){const row=new Array(J);
      for(let j=0;j<J;j++)row[j]=mi[j]+Math.floor(rng()*(mx[j]-mi[j]+1));out.push(row);}
    else if(kind==="colperm"){const row=new Array(J);
      for(let j=0;j<J;j++)row[j]=rows[Math.floor(rng()*N)][j];out.push(row);}
    else {const pool=(loHalfIdx&&loHalfIdx.length)?loHalfIdx:[...Array(N).keys()];
      const src=rows[pool[Math.floor(rng()*pool.length)]].slice();
      for(let t=src.length-1;t>0;t--){const q=Math.floor(rng()*(t+1));const tmp=src[t];src[t]=src[q];src[q]=tmp;}
      out.push(src);}}
  return out;}
/* ---------- estimators ---------- */
function survival(x,t){let c=0;for(const v of x)if(v>t)c++;return c/x.length;}
function leftFit(x){const n=x.length,s=[...x].sort((a,b)=>a-b);
  const q1=s[Math.floor(.25*n)],q3=s[Math.floor(.75*n)],iqr=q3-q1;
  const mu=mn(x),sd=Math.sqrt(mn(x.map(v=>(v-mu)*(v-mu))));
  let h=0.9*Math.min(sd||1,(iqr/1.34)||(sd||1))*Math.pow(n,-0.2);if(!(h>0))h=(sd||1)*0.3;
  const G=512,lo=s[0]-h,hi=s[n-1]+h,d=new Float64Array(G);
  for(let g=0;g<G;g++){const xg=lo+(hi-lo)*g/(G-1);let v=0;for(let i=0;i<n;i++){const zz=(xg-x[i])/h;v+=Math.exp(-0.5*zz*zz);}d[g]=v;}
  let dm=0;for(let g=0;g<G;g++)if(d[g]>dm)dm=d[g];
  let m=med(x),found=false;
  for(let g=1;g<G-1;g++)if(d[g]>=d[g-1]&&d[g]>=d[g+1]&&d[g]>0.15*dm){m=lo+(hi-lo)*g/(G-1);found=true;break;}
  const below=x.filter(v=>v<m);let sg=below.length?1.4826*med(below.map(v=>m-v)):0;if(!(sg>0))sg=1e-9;
  return {m,sd:sg,nBelow:below.length};}
function sol1(x){const {m,sd,nBelow}=leftFit(x),n=x.length;
  const c=m+sd;
  return {mirror:clip(1-2*nBelow/n), sigma:clip(1-(x.filter(v=>v<c).length/0.8413447)/n)};}
function synEst(obs,syn){const {m,sd}=leftFit(obs);
  const all=[...obs].sort((a,b)=>a-b),grid=[];
  for(let q=0.50;q<=0.97;q+=0.01)grid.push(all[Math.min(all.length-1,Math.floor(q*all.length))]);
  const rat=[],cor=[];
  for(const t of grid){const So=survival(obs,t),Ss=survival(syn,t),Sc=1-Phi((t-m)/sd);
    if(Ss>0.10&&Ss<0.98)rat.push(So/Ss);
    if((Ss-Sc)>0.08&&Ss<0.98)cor.push((So-Sc)/(Ss-Sc));}
  return {ratio:rat.length?clip(med(rat)):NaN, corr:cor.length?clip(med(cor)):NaN};}
function gateOK(gm){const p=Math.sqrt((gm.v1+gm.v2)/2)||1;return (gm.bic1-gm.bic2>2)&&gm.pi>0.01&&gm.pi<0.7&&(gm.m2-gm.m1)/p>1.0;}
/* ---------- data ---------- */
const DIR="ext_bench2";
function loadCSVnoHeader(p){return fs.readFileSync(p,"utf8").trim().split(/\r?\n/).map(l=>l.split(",").map(Number));}
function loadY(p){const L=fs.readFileSync(p,"utf8").trim().split(/\r?\n/);return L.slice(1).map(l=>Number(l.split(",")[0]));}
const idx=fs.readFileSync(DIR+"/index.csv","utf8").trim().split(/\r?\n/).slice(1)
  .map(l=>{const c=l.split(",");return {name:c[0],J:+c[2],kind:c[4]};});
/* Study 1 as a 17th dataset (induced labels) */
function loadStudy1(){const Ms=R.parseCSV(fs.readFileSync("_study_matrix.csv","utf8"));
  const rows=Ms.rows.map(r=>r.map(Number));
  const y=R.parseCSV(fs.readFileSync("_study_labels.csv","utf8")).rows.map(r=>Number(r[0]));
  return {name:"Study1_induced",J:Ms.header.length,kind:"induced",rows,y};}
const RATES=[0.05,0.15,0.25,0.35], NTARGET=300, REPS=6, ITERS=80;
const COLS=["mirror","sigma","eta","colperm_corr","rowLo_ratio","unif_corr"];
const out=[["dataset","kind","J","target_rate","actual_rate","n","reps"].concat(COLS.flatMap(c=>["err_"+c,"bias_"+c]))];
function runRep(rows,J){
  const nR=rows.length;
  const M0=new Float64Array(nR*J);for(let i=0;i<nR;i++)for(let j=0;j<J;j++)M0[i*J+j]=rows[i][j];
  const obs=etaFrozen(M0,nR,J,[...Array(nR).keys()],ITERS,1);
  const lo=[...Array(nR).keys()].sort((a,b)=>obs[a]-obs[b]).slice(0,Math.floor(nR/2));
  const s1=sol1(obs),gm=R.gaussMix1D(obs);
  const res={mirror:s1.mirror,sigma:s1.sigma,eta:gateOK(gm)?gm.pi:0};
  for(const [kind,tag,which] of [["colperm","colperm_corr","corr"],["rowshufLo","rowLo_ratio","ratio"],["unif","unif_corr","corr"]]){
    const syn=genSynthetic(rows,J,kind,nR,lo);
    const nT=nR+syn.length,MT=new Float64Array(nT*J);
    for(let i=0;i<nR;i++)for(let j=0;j<J;j++)MT[i*J+j]=rows[i][j];
    for(let s=0;s<syn.length;s++)for(let j=0;j<J;j++)MT[(nR+s)*J+j]=syn[s][j];
    const eta=etaFrozen(MT,nT,J,[...Array(nR).keys()],ITERS,1);
    const e=synEst(eta.slice(0,nR),eta.slice(nR));
    res[tag]=which==="corr"?e.corr:e.ratio;}
  return res;}
console.log("dataset          kind      J   target  actual   n  | "+COLS.map(c=>c.slice(0,9).padStart(10)).join(""));
const DATASETS=[loadStudy1()].concat(idx.map(d=>{
  try{return {name:d.name,J:d.J,kind:d.kind,rows:loadCSVnoHeader(`${DIR}/${d.name}_M.csv`),y:loadY(`${DIR}/${d.name}_y.csv`)};}
  catch(e){return null;}}).filter(Boolean));
for(const D of DATASETS){
  const J=D.J, rows=D.rows, y=D.y;
  const neg=y.map((v,i)=>v===0?i:-1).filter(i=>i>=0), pos=y.map((v,i)=>v===1?i:-1).filter(i=>i>=0);
  for(const tr of RATES){
    let k=Math.round(tr*NTARGET); if(k>pos.length)k=pos.length;
    let nCf=Math.round(k*(1-tr)/tr); if(nCf>neg.length){nCf=neg.length;k=Math.round(tr*nCf/(1-tr));}
    if(k<5||nCf<40||k>pos.length) continue;
    const actual=k/(k+nCf), res=[];
    for(let rep=0;rep<REPS;rep++){
      const pk=pos.slice(),nk=neg.slice();
      for(let i=0;i<k;i++){const j=i+Math.floor(rng()*(pk.length-i));const t=pk[i];pk[i]=pk[j];pk[j]=t;}
      for(let i=0;i<nCf;i++){const j=i+Math.floor(rng()*(nk.length-i));const t=nk[i];nk[i]=nk[j];nk[j]=t;}
      const sel=[...nk.slice(0,nCf),...pk.slice(0,k)].map(i=>rows[i]);
      res.push(runRep(sel,J));}
    const row=[D.name,D.kind,J,tr.toFixed(2),actual.toFixed(4),k+nCf,REPS];
    for(const c of COLS){const v=res.map(r=>r[c]).filter(Number.isFinite);
      row.push(mn(v.map(x=>Math.abs(x-actual))).toFixed(4),mn(v.map(x=>x-actual)).toFixed(4));}
    out.push(row);
    const f=c=>{const v=res.map(r=>r[c]).filter(Number.isFinite);
      return v.length?(100*mn(v.map(x=>Math.abs(x-actual)))).toFixed(1).padStart(10):"        NA";};
    console.log(D.name.padEnd(17)+D.kind.padEnd(9)+String(J).padStart(4)+"  "+(100*tr).toFixed(0).padStart(5)+"%  "+
      (100*actual).toFixed(1).padStart(5)+"% "+String(k+nCf).padStart(4)+" |"+COLS.map(f).join(""));}
  fs.writeFileSync("broad_eval.csv",out.map(r=>r.join(",")).join("\n"));}
console.log("\nwrote broad_eval.csv");
