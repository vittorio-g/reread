/* eval_rightside.js — MIRROR OF SOLUTION 1: estimate the CARELESS component from the RIGHT side
 * only (the top of the score distribution), instead of the attentive component from the left.
 *
 * Solution 1 works because the left tail is careless-free by construction. The symmetric bet is
 * that the right tail is careful-free, so the careless component can be estimated there and its
 * mass gives pi directly (no subtraction).
 *
 *   M        = right-most prominent mode of the score distribution (KDE)
 *   rMirror  : careless symmetric about M      -> pi = 2*#{x > M}/n
 *   r1sigma  : sigma from the RIGHT side only (MAD above M); cut c = M - sigma,
 *              expected careless share above c is Phi(1)=.8413 -> pi = (#{x>c}/.8413)/n
 *   rTail(q) : uses ONLY the top q% (Vittorio's "top 5%"). The top q% is the careless component
 *              truncated at t=quantile(1-q); a truncated-normal moment fit recovers (mu1,sigma1),
 *              hence lambda=(t-mu1)/sigma1, and since the observed mass above t is exactly q,
 *              pi = q / (1 - Phi(lambda)).
 *
 * Predictions to falsify: right-side estimation should be strongest where solution 1 is weakest
 * (low prevalence is the risky regime here: with few careless the top q% gets contaminated by
 * careful, and the heterogeneous careless subtypes make the tail unrepresentative).
 *
 * Comparators: one-sided LEFT (mirror, 1sigma) and the shipped mixture. Same bench as
 * eval_onesided.js / eval_anchor.js. Output: rightside_eval.csv (+ rightside_diag.csv)
 */
const R=require("./site/reread.js"); const W=R.WEIGHTS; const fs=require("fs");
let SEED=4242; const rng=()=>{SEED=(SEED*1103515245+12345)&0x7fffffff;return SEED/0x7fffffff;};
function gauss(){let u=0,v=0;while(!u)u=rng();while(!v)v=rng();return Math.sqrt(-2*Math.log(u))*Math.cos(2*Math.PI*v);}
const K=5, THR=[-0.84,-0.25,0.25,0.84];
function genClean(F,ipf,n){const J=F*ipf;const load=[],fac=[],tau=[];
  for(let j=0;j<J;j++){let l=0.4+0.4*rng();if(rng()<0.4)l=-l;load.push(l);fac.push(Math.floor(j/ipf));tau.push(0.45*gauss());}
  const rho=0.3,sq=Math.sqrt(rho),sq1=Math.sqrt(1-rho);const rows=[];
  for(let i=0;i<n;i++){const c=gauss(),f=[];for(let k=0;k<F;k++)f.push(sq*c+sq1*gauss());
    const row=new Array(J);
    for(let j=0;j<J;j++){let y=tau[j]+load[j]*f[fac[j]];y+=Math.sqrt(Math.max(0.05,1-load[j]*load[j]))*gauss();
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
const toMat=(rows,J)=>{const n=rows.length,m=new Float64Array(n*J);for(let i=0;i<n;i++)for(let j=0;j<J;j++)m[i*J+j]=rows[i][j];return m;};
function gateOK(gm){const p=Math.sqrt((gm.v1+gm.v2)/2)||1;return (gm.bic1-gm.bic2>2)&&gm.pi>0.01&&gm.pi<0.7&&(gm.m2-gm.m1)/p>1.0;}
const med=a=>{const s=[...a].sort((x,y)=>x-y);const h=s.length>>1;return s.length%2?s[h]:(s[h-1]+s[h])/2;};
const mn=a=>a.length?a.reduce((s,v)=>s+v,0)/a.length:NaN;
const clip=p=>Math.max(0,Math.min(1,p));
/* normal pdf/cdf */
const phi=z=>Math.exp(-0.5*z*z)/Math.sqrt(2*Math.PI);
function Phi(z){const t=1/(1+0.2316419*Math.abs(z));
  const d=phi(z)*t*(0.319381530+t*(-0.356563782+t*(1.781477937+t*(-1.821255978+t*1.330274429))));
  return z>=0?1-d:d;}
const PHI1=0.8413447;
/* KDE modes: leftmost (for the solution-1 comparators) and rightmost (for the new estimators) */
function kdeModes(x){const n=x.length,s=[...x].sort((a,b)=>a-b);
  const q1=s[Math.floor(.25*n)],q3=s[Math.floor(.75*n)],iqr=q3-q1;
  const mu=mn(x), sd=Math.sqrt(x.reduce((a,b)=>a+(b-mu)*(b-mu),0)/n);
  let h=0.9*Math.min(sd||1,(iqr/1.34)||(sd||1))*Math.pow(n,-0.2); if(!(h>0))h=(sd||1)*0.3;
  const G=512,lo=s[0]-h,hi=s[n-1]+h,dens=new Float64Array(G);
  for(let g=0;g<G;g++){const xg=lo+(hi-lo)*g/(G-1);let d=0;
    for(let i=0;i<n;i++){const z=(xg-x[i])/h;d+=Math.exp(-0.5*z*z);} dens[g]=d;}
  let dmax=0;for(let g=0;g<G;g++)if(dens[g]>dmax)dmax=dens[g];
  const peaks=[];
  for(let g=1;g<G-1;g++) if(dens[g]>=dens[g-1]&&dens[g]>=dens[g+1]&&dens[g]>0.15*dmax) peaks.push(lo+(hi-lo)*g/(G-1));
  if(!peaks.length){const m=med(x);return {left:m,right:m};}
  return {left:peaks[0], right:peaks[peaks.length-1]};}
/* LEFT-side comparators (solution 1) */
function leftSide(x,m){const n=x.length,below=x.filter(v=>v<m);
  const piMirror=clip(1-2*below.length/n);
  let sd=below.length?1.4826*med(below.map(v=>m-v)):0; if(!(sd>0))sd=1e-9;
  const c=m+sd;
  return {piMirror, piSigma:clip(1-(x.filter(v=>v<c).length/PHI1)/n)};}
/* RIGHT-side estimators (new) */
function rightSide(x,M){const n=x.length,above=x.filter(v=>v>M);
  const rMirror=clip(2*above.length/n);
  let sd=above.length?1.4826*med(above.map(v=>v-M)):0; if(!(sd>0))sd=1e-9;
  const c=M-sd;
  return {rMirror, rSigma:clip((x.filter(v=>v>c).length/PHI1)/n)};}
/* "top q% ONLY": truncated-normal moment fit on the top q% -> pi = q/(1-Phi(lambda)) */
function rTail(x,q){const n=x.length,s=[...x].sort((a,b)=>a-b);
  const k=Math.max(4,Math.round(q*n)); const tail=s.slice(n-k); const t=s[n-k-1>=0?n-k-1:0];
  const mt=mn(tail); const st=Math.sqrt(mn(tail.map(v=>(v-mt)*(v-mt))));
  if(!(st>0)) return NaN;
  const h=l=>{const d=1-Phi(l); return d>1e-12? phi(l)/d : (l+Math.sqrt(l*l+4))/2;};   // inverse Mills
  // solve f(lambda) = lambda - h(lambda) - (t-mt)/sigma(lambda) = 0, sigma from the variance eq.
  const sig=l=>{const hl=h(l), v=1+l*hl-hl*hl; return v>1e-9? st/Math.sqrt(v) : NaN;};
  const f=l=>{const sg=sig(l); if(!Number.isFinite(sg)||sg<=0) return NaN; return l-h(l)-(t-mt)/sg;};
  let lo=-6,hi=6,flo=f(lo),fhi=f(hi);
  if(!Number.isFinite(flo)||!Number.isFinite(fhi)||flo*fhi>0) return NaN;
  for(let it=0;it<80;it++){const mid=(lo+hi)/2,fm=f(mid);
    if(!Number.isFinite(fm)) return NaN;
    if(flo*fm<=0){hi=mid;fhi=fm;} else {lo=mid;flo=fm;} }
  const lam=(lo+hi)/2, tailProb=1-Phi(lam);
  return tailProb>1e-6 ? clip(q/tailProb) : NaN; }
function etaOf(mat,n,J,iters){const ef=R.ensembleFeatures(mat,n,J,{iterations:iters,seed:1}),O=ef.oriented;
  const g=Math.max(0,Math.min(1,(ef.profileInfo-0.01)/0.05));
  const e=new Array(n);for(let i=0;i<n;i++)e[i]=W.b0+W.rr*O.rr[i]+g*(W.longstring*O.longstring[i]+W.person_total*O.person_total[i]);
  return e;}
function estimate(e,gm){const {left,right}=kdeModes(e);
  return {...leftSide(e,left), ...rightSide(e,right),
          rTail5:rTail(e,0.05), rTail10:rTail(e,0.10), rTail20:rTail(e,0.20),
          piEta:gateOK(gm)?gm.pi:0};}
const COLS=["piSigma","piEta","rMirror","rSigma","rTail5","rTail10","rTail20"];
const out=[["study","J","true_rate","reps"].concat(COLS.flatMap(c=>["err_"+c,"bias_"+c,"na_"+c]))];
function report(study,J,tr,res,reps){
  const row=[study,J,tr.toFixed(4),reps];
  for(const c of COLS){const v=res.map(r=>r[c]).filter(Number.isFinite);
    row.push(mn(v.map(x=>Math.abs(x-tr))).toFixed(4), mn(v.map(x=>x-tr)).toFixed(4), (1-v.length/res.length).toFixed(2));}
  out.push(row);
  const f=c=>{const v=res.map(r=>r[c]).filter(Number.isFinite);
    return v.length? (100*mn(v.map(x=>Math.abs(x-tr)))).toFixed(1).padStart(5) : "   NA";};
  const b=c=>{const v=res.map(r=>r[c]).filter(Number.isFinite);
    return v.length? (100*mn(v.map(x=>x-tr))).toFixed(1).padStart(6) : "    NA";};
  console.log(study.padEnd(7)+"J="+String(J).padStart(3)+" true="+(100*tr).toFixed(0).padStart(3)+"% | |err| L1sig"+f("piSigma")+
    " eta"+f("piEta")+" || rMir"+f("rMirror")+" rSig"+f("rSigma")+" rT5"+f("rTail5")+" rT10"+f("rTail10")+" rT20"+f("rTail20")+
    " | bias rT5"+b("rTail5")+" rSig"+b("rSigma"));}
/* diagnostic: purity of the top q% (labels known) */
const diag=[["study","true_rate","q","purity_top","n_top","mu1_true","tailmean"]];
function diagnose(study,tr,e,lab){const n=e.length,idx=[...Array(n).keys()].sort((a,b)=>e[a]-e[b]);
  const mu1=mn(e.filter((_,i)=>lab[i]===1));
  for(const q of [0.05,0.10,0.20]){const k=Math.max(4,Math.round(q*n)),hiI=idx.slice(n-k);
    diag.push([study,tr.toFixed(4),q,(hiI.filter(i=>lab[i]===1).length/k).toFixed(3),k,
               Number.isFinite(mu1)?mu1.toFixed(3):"NA",mn(hiI.map(i=>e[i])).toFixed(3)]);}}
console.log("=== STUDY 1 (real data, prevalence swept) ===");
{ const Ms=R.parseCSV(fs.readFileSync("_study_matrix.csv","utf8")),J=Ms.header.length;
  const allR=Ms.rows.map(r=>r.map(v=>{const x=Number(v);return Number.isFinite(x)?x:NaN;}));
  const y=R.parseCSV(fs.readFileSync("_study_labels.csv","utf8")).rows.map(r=>Number(r[0]));
  const cf=y.map((v,i)=>v===0?i:-1).filter(i=>i>=0), cl=y.map((v,i)=>v===1?i:-1).filter(i=>i>=0);
  const REPS=40;
  for(const k of [4,8,12,18,24,33,47,58]){ if(k>cl.length) continue; const tr=k/(cf.length+k),res=[];
    for(let rep=0;rep<REPS;rep++){ const a=cl.slice();
      for(let i=0;i<k;i++){const j=i+Math.floor(rng()*(a.length-i));const t=a[i];a[i]=a[j];a[j]=t;}
      const idx=[...cf,...a.slice(0,k)],n=idx.length,mat=new Float64Array(n*J);
      for(let i=0;i<n;i++)for(let j=0;j<J;j++)mat[i*J+j]=allR[idx[i]][j];
      const e=etaOf(mat,n,J,150), gm=R.gaussMix1D(e);
      res.push(estimate(e,gm));
      if(rep===0) diagnose("Study1",tr,e,idx.map((_,i)=>i<cf.length?0:1)); }
    report("Study1",J,tr,res,REPS); } }
console.log("\n=== STUDY 3 (simulated, six patterns) ===");
{ const REPS=10,N=400;
  for(const [F,ipf] of [[10,6],[30,6]]){ const Jn=F*ipf;
    for(const pct of [0.05,0.10,0.15,0.20,0.25,0.30,0.35,0.40,0.45]){ const res=[];
      for(let rep=0;rep<REPS;rep++){ const {rows,J}=genClean(F,ipf,N); const {out:o,lab}=inject(rows,J,pct);
        const e=etaOf(toMat(o,J),N,J,100), gm=R.gaussMix1D(e);
        res.push(estimate(e,gm));
        if(rep===0) diagnose("Study3_J"+Jn,pct,e,lab); }
      report("Study3",Jn,pct,res,REPS); } } }
fs.writeFileSync("rightside_eval.csv",out.map(r=>r.join(",")).join("\n"));
fs.writeFileSync("rightside_diag.csv",diag.map(r=>r.join(",")).join("\n"));
console.log("\nwrote rightside_eval.csv + rightside_diag.csv");
