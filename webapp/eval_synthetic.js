/* eval_synthetic.js — SOLUTION 2: estimate the prevalence with a SYNTHETIC careless reference.
 *
 * The previous two experiments established the asymmetry: the attentive form an estimable MODE,
 * the careless form a heavy tail with no mode of its own, so the careless component cannot be
 * read off the data. This script does not try to: it MANUFACTURES the careless reference by
 * simulating content-non-responsive answering, scores those synthetic rows in the FROZEN frame of
 * the real sample (they never touch pair selection, sign alignment, the Person-Total profile, the
 * median/MAD standardization or the gate), and uses their score distribution as a calibration
 * curve: knowing what fraction of a careless respondent's scores exceed t converts "how many
 * observations exceed t" into "how many careless are present".
 *
 * Synthetic careless generators (all derived from the real matrix, as deployment would):
 *   unif      : uniform random category per item over the observed item range
 *   colperm   : each item column permuted across respondents (keeps item marginals, kills profile)
 *   rowshuf   : a real respondent's own answers permuted across items (keeps their response STYLE
 *               - acquiescence, variance, straight-lining - and kills only content responsiveness)
 *   rowshufLo : same, but seeded only from the lower half of the score distribution (the "careless
 *               version of a typical attentive person")
 *
 * Estimators, with S(t) the survival function on a threshold grid:
 *   synRatio  : pi = median_t S_obs(t)/S_syn(t)                     (no attentive model at all)
 *   synCorr   : pi = median_t [S_obs - S_cf]/[S_syn - S_cf]         (S_cf = normal fit of the
 *               attentive mode from the LEFT side, i.e. solution 1 - uses only what IS estimable)
 *   synLS     : grid pi minimizing || S_obs - ((1-pi) S_cf + pi S_syn) ||^2
 *
 * Comparators: one-sided 1sigma (solution 1) and the shipped mixture. Output: synthetic_eval.csv
 */
const R=require("./site/reread.js"); const W=R.WEIGHTS; const fs=require("fs");
let SEED=7717; const rng=()=>{SEED=(SEED*1103515245+12345)&0x7fffffff;return SEED/0x7fffffff;};
function gaussR(){let u=0,v=0;while(!u)u=rng();while(!v)v=rng();return Math.sqrt(-2*Math.log(u))*Math.cos(2*Math.PI*v);}
const mn=a=>a.length?a.reduce((s,v)=>s+v,0)/a.length:NaN;
const med=a=>{const s=[...a].sort((x,y)=>x-y);const h=s.length>>1;return s.length%2?s[h]:(s[h-1]+s[h])/2;};
const clip=p=>Math.max(0,Math.min(1,p));
const phi=z=>Math.exp(-0.5*z*z)/Math.sqrt(2*Math.PI);
function Phi(z){const t=1/(1+0.2316419*Math.abs(z));
  const d=phi(z)*t*(0.319381530+t*(-0.356563782+t*(1.781477937+t*(-1.821255978+t*1.330274429))));
  return z>=0?1-d:d;}
/* ---------- frozen-frame scoring machinery (verbatim from nested_loo.js) ---------- */
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
/* eta for ALL rows, every sample-level quantity taken from trainIdx (the real sample) */
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
  // gate from the TRAIN slice only
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
function subMatrix(M,J,idx){const m=new Float64Array(idx.length*J);for(let a=0;a<idx.length;a++)for(let j=0;j<J;j++)m[a*J+j]=M[idx[a]*J+j];return m;}
function subRows(P,J,idx){const m=new Float64Array(idx.length*J);for(let a=0;a<idx.length;a++)for(let j=0;j<J;j++)m[a*J+j]=P[idx[a]*J+j];return m;}
/* ---------- synthetic careless generators (from the real rows) ---------- */
function genSynthetic(rows,J,kind,nSyn,loHalfIdx){const N=rows.length,out=[];
  const mx=new Array(J).fill(-Infinity),mi=new Array(J).fill(Infinity);
  for(const r of rows)for(let j=0;j<J;j++){const v=r[j];if(Number.isFinite(v)){if(v>mx[j])mx[j]=v;if(v<mi[j])mi[j]=v;}}
  for(let s=0;s<nSyn;s++){
    if(kind==="unif"){const row=new Array(J);
      for(let j=0;j<J;j++)row[j]=mi[j]+Math.floor(rng()*(mx[j]-mi[j]+1));out.push(row);}
    else if(kind==="colperm"){const row=new Array(J);
      for(let j=0;j<J;j++)row[j]=rows[Math.floor(rng()*N)][j];out.push(row);}   // independent draw per column
    else { const pool=(kind==="rowshufLo"&&loHalfIdx&&loHalfIdx.length)?loHalfIdx:[...Array(N).keys()];
      const src=rows[pool[Math.floor(rng()*pool.length)]].slice();
      for(let t=src.length-1;t>0;t--){const q=Math.floor(rng()*(t+1));const tmp=src[t];src[t]=src[q];src[q]=tmp;}
      out.push(src);} }
  return out;}
/* ---------- estimators ---------- */
function survival(x,t){let c=0;for(const v of x)if(v>t)c++;return c/x.length;}
function leftFit(x){ // solution-1 attentive fit: left-most KDE mode + left-side MAD
  const n=x.length,s=[...x].sort((a,b)=>a-b);
  const q1=s[Math.floor(.25*n)],q3=s[Math.floor(.75*n)],iqr=q3-q1;
  const mu=mn(x),sd=Math.sqrt(mn(x.map(v=>(v-mu)*(v-mu))));
  let h=0.9*Math.min(sd||1,(iqr/1.34)||(sd||1))*Math.pow(n,-0.2);if(!(h>0))h=(sd||1)*0.3;
  const G=512,lo=s[0]-h,hi=s[n-1]+h,d=new Float64Array(G);
  for(let g=0;g<G;g++){const xg=lo+(hi-lo)*g/(G-1);let v=0;for(let i=0;i<n;i++){const zz=(xg-x[i])/h;v+=Math.exp(-0.5*zz*zz);}d[g]=v;}
  let dm=0;for(let g=0;g<G;g++)if(d[g]>dm)dm=d[g];
  let m=med(x);for(let g=1;g<G-1;g++)if(d[g]>=d[g-1]&&d[g]>=d[g+1]&&d[g]>0.15*dm){m=lo+(hi-lo)*g/(G-1);break;}
  const below=x.filter(v=>v<m);let sg=below.length?1.4826*med(below.map(v=>m-v)):0;if(!(sg>0))sg=1e-9;
  return {m,sd:sg};}
function estimators(obs,syn){
  const {m,sd}=leftFit(obs);
  const all=[...obs].sort((a,b)=>a-b);
  const grid=[];for(let q=0.50;q<=0.97;q+=0.01)grid.push(all[Math.min(all.length-1,Math.floor(q*all.length))]);
  const rat=[],cor=[];
  for(const t of grid){const So=survival(obs,t),Ss=survival(syn,t),Sc=1-Phi((t-m)/sd);
    if(Ss>0.10&&Ss<0.98) rat.push(So/Ss);
    if((Ss-Sc)>0.08&&Ss<0.98) cor.push((So-Sc)/(Ss-Sc));}
  let best=0,bd=Infinity;
  for(let p=0;p<=1.0001;p+=0.005){let d=0,cnt=0;
    for(const t of grid){const So=survival(obs,t),Ss=survival(syn,t),Sc=1-Phi((t-m)/sd);
      const e=So-((1-p)*Sc+p*Ss);d+=e*e;cnt++;}
    if(d<bd){bd=d;best=p;}}
  return {ratio:rat.length?clip(med(rat)):NaN, corr:cor.length?clip(med(cor)):NaN, ls:clip(best)};}
function oneSidedSigma(x){const {m,sd}=leftFit(x);const c=m+sd;
  return clip(1-(x.filter(v=>v<c).length/0.8413447)/x.length);}
function gateOK(gm){const p=Math.sqrt((gm.v1+gm.v2)/2)||1;return (gm.bic1-gm.bic2>2)&&gm.pi>0.01&&gm.pi<0.7&&(gm.m2-gm.m1)/p>1.0;}
/* ---------- simulated data (Study 3), identical to the other benches ---------- */
const K=5,THR=[-0.84,-0.25,0.25,0.84];
function genClean(F,ipf,n){const J=F*ipf;const load=[],fac=[],tau=[];
  for(let j=0;j<J;j++){let l=0.4+0.4*rng();if(rng()<0.4)l=-l;load.push(l);fac.push(Math.floor(j/ipf));tau.push(0.45*gaussR());}
  const rho=0.3,sq=Math.sqrt(rho),sq1=Math.sqrt(1-rho);const rows=[];
  for(let i=0;i<n;i++){const c=gaussR(),f=[];for(let k=0;k<F;k++)f.push(sq*c+sq1*gaussR());
    const row=new Array(J);
    for(let j=0;j<J;j++){let y=tau[j]+load[j]*f[fac[j]];y+=Math.sqrt(Math.max(0.05,1-load[j]*load[j]))*gaussR();
      let cat=1;for(let q=0;q<K-1;q++)if(y>THR[q])cat=q+2;row[j]=cat;}
    rows.push(row);}
  return {rows,J};}
function inject(rows,J,pct){const N=rows.length,out=rows.map(r=>r.slice()),lab=new Array(N).fill(0);
  const M=Math.round(N*pct),ord=[...Array(N).keys()];
  for(let i=0;i<M;i++){const j=i+Math.floor(rng()*(N-i));const t=ord[i];ord[i]=ord[j];ord[j]=t;}
  const pats=["random","longstring","pure_straight","acquiescent","mixed","fatigue"];
  for(let m=0;m<M;m++){const i=ord[m],pat=pats[m%6],lvl=0.5+0.5*rng(),k=Math.max(1,Math.round(J*lvl)),row=out[i];
    const idx=[...Array(J).keys()];for(let t=0;t<k;t++){const j=t+Math.floor(rng()*(J-t));const tm=idx[t];idx[t]=idx[j];idx[j]=tm;}
    const sel=idx.slice(0,k);
    if(pat==="pure_straight"||pat==="longstring"){const v=1+Math.floor(rng()*K);sel.forEach(p=>row[p]=v);}
    else if(pat==="acquiescent"){sel.forEach(p=>row[p]=Math.max(1,Math.min(K,4+Math.round(rng()))));}
    else if(pat==="fatigue"){for(let t=J-k;t<J;t++)row[t]=1+Math.floor(rng()*K);}
    else sel.forEach(p=>row[p]=1+Math.floor(rng()*K));
    lab[i]=1;}
  return {out,lab};}
const KINDS=["unif","colperm","rowshuf","rowshufLo"];
const COLS=["piSigma","piEta"].concat(KINDS.flatMap(k=>[k+"_ratio",k+"_corr",k+"_ls"]));
const out=[["study","J","true_rate","reps"].concat(COLS.flatMap(c=>["err_"+c,"bias_"+c]))];
function runOne(rows,J,iters,seed){
  const nR=rows.length;
  // first pass: score the real sample alone to get eta_obs (deployment view) and the low half
  const M0=new Float64Array(nR*J);for(let i=0;i<nR;i++)for(let j=0;j<J;j++)M0[i*J+j]=rows[i][j];
  const obs=etaFrozen(M0,nR,J,[...Array(nR).keys()],iters,seed);
  const ordLo=[...Array(nR).keys()].sort((a,b)=>obs[a]-obs[b]).slice(0,Math.floor(nR/2));
  const res={};
  for(const kind of KINDS){
    const syn=genSynthetic(rows,J,kind,nR,ordLo);
    const nT=nR+syn.length,MT=new Float64Array(nT*J);
    for(let i=0;i<nR;i++)for(let j=0;j<J;j++)MT[i*J+j]=rows[i][j];
    for(let s=0;s<syn.length;s++)for(let j=0;j<J;j++)MT[(nR+s)*J+j]=syn[s][j];
    const eta=etaFrozen(MT,nT,J,[...Array(nR).keys()],iters,seed);
    const e=estimators(eta.slice(0,nR),eta.slice(nR));
    res[kind+"_ratio"]=e.ratio;res[kind+"_corr"]=e.corr;res[kind+"_ls"]=e.ls;}
  res.piSigma=oneSidedSigma(obs);
  const gm=R.gaussMix1D(obs);res.piEta=gateOK(gm)?gm.pi:0;
  return res;}
function report(study,J,tr,res,reps){
  const row=[study,J,tr.toFixed(4),reps];
  for(const c of COLS){const v=res.map(r=>r[c]).filter(Number.isFinite);
    row.push(mn(v.map(x=>Math.abs(x-tr))).toFixed(4),mn(v.map(x=>x-tr)).toFixed(4));}
  out.push(row);
  const f=c=>{const v=res.map(r=>r[c]).filter(Number.isFinite);return v.length?(100*mn(v.map(x=>Math.abs(x-tr)))).toFixed(1).padStart(5):"   NA";};
  console.log(study.padEnd(7)+"J="+String(J).padStart(3)+" true="+(100*tr).toFixed(0).padStart(3)+"% | 1sig"+f("piSigma")+" eta"+f("piEta")+
    " || unif"+f("unif_corr")+" colp"+f("colperm_corr")+" rowsh"+f("rowshuf_corr")+" rowLo"+f("rowshufLo_corr")+
    " | ratio(rowLo)"+f("rowshufLo_ratio")+" ls(rowLo)"+f("rowshufLo_ls"));}
console.log("=== STUDY 1 (real data, prevalence swept) — corr = synthetic + left-side attentive fit ===");
{ const Ms=R.parseCSV(fs.readFileSync("_study_matrix.csv","utf8")),J=Ms.header.length;
  const allR=Ms.rows.map(r=>r.map(v=>{const x=Number(v);return Number.isFinite(x)?x:NaN;}));
  const y=R.parseCSV(fs.readFileSync("_study_labels.csv","utf8")).rows.map(r=>Number(r[0]));
  const cf=y.map((v,i)=>v===0?i:-1).filter(i=>i>=0), cl=y.map((v,i)=>v===1?i:-1).filter(i=>i>=0);
  const REPS=20;
  for(const k of [4,8,12,18,24,33,47,58]){ if(k>cl.length) continue; const tr=k/(cf.length+k),res=[];
    for(let rep=0;rep<REPS;rep++){ const a=cl.slice();
      for(let i=0;i<k;i++){const j=i+Math.floor(rng()*(a.length-i));const t=a[i];a[i]=a[j];a[j]=t;}
      const idx=[...cf,...a.slice(0,k)],rows=idx.map(i=>allR[i]);
      res.push(runOne(rows,J,120,1)); }
    report("Study1",J,tr,res,REPS); } }
console.log("\n=== STUDY 3 (simulated, six patterns) ===");
{ const REPS=6,N=300;
  for(const [F,ipf] of [[10,6],[30,6]]){ const Jn=F*ipf;
    for(const pct of [0.05,0.15,0.25,0.35,0.45]){ const res=[];
      for(let rep=0;rep<REPS;rep++){ const {rows,J}=genClean(F,ipf,N); const {out:o}=inject(rows,J,pct);
        res.push(runOne(o,J,80,1)); }
      report("Study3",Jn,pct,res,REPS); } } }
fs.writeFileSync("synthetic_eval.csv",out.map(r=>r.join(",")).join("\n"));
console.log("\nwrote synthetic_eval.csv");
