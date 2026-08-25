/* sens_sim.js — reviewer 2.2/2.3/2.4 on simulated data with KNOWN structure.
 *  A) gate false-alarm: many genuinely clean datasets -> how often does the BIC
 *     gate fire / flag anyone?  (quantifies "on clean data it does not").
 *  B) pair-selection contamination vs prevalence: fraction of the selected
 *     top-3% pairs that are genuinely same-factor, as prevalence 0->40%.
 *  C) adverse conditions: low loadings / bifactor / high prevalence -> rr AUC.
 *  D) pi coverage vs overlap: pi-hat vs true prevalence at weak vs strong sep. */
const R=require("./site/reread.js"); const W=R.WEIGHTS; const fs=require("fs");
let SEED=12345; const rng=()=>{SEED=(SEED*1103515245+12345)&0x7fffffff;return SEED/0x7fffffff;};
function gauss(){let u=0,v=0;while(!u)u=rng();while(!v)v=rng();return Math.sqrt(-2*Math.log(u))*Math.cos(2*Math.PI*v);}
const K=5, THR=[-0.84,-0.25,0.25,0.84];
// genClean with options: loRange=[lo,hi] loadings, bifactor adds a general factor g.
function genClean(F,ipf,n,opts){opts=opts||{};const J=F*ipf,lo=opts.loRange?opts.loRange[0]:0.4,hi=opts.loRange?opts.loRange[1]:0.8;
  const gLoad=opts.bifactor?0.4:0; const load=[],fac=[],tau=[];
  for(let j=0;j<J;j++){let l=lo+(hi-lo)*rng();if(rng()<0.4)l=-l;load.push(l);fac.push(Math.floor(j/ipf));tau.push(0.45*gauss());}
  const rho=0.3,sq=Math.sqrt(rho),sq1=Math.sqrt(1-rho);const rows=[];
  for(let i=0;i<n;i++){const common=gauss(),gen=gauss(),f=[];for(let k=0;k<F;k++)f.push(sq*common+sq1*gauss());
    const row=new Array(J);
    for(let j=0;j<J;j++){let y=tau[j]+load[j]*f[fac[j]]+ (opts.bifactor?gLoad*gen:0);
      const uniq=Math.sqrt(Math.max(0.05,1-load[j]*load[j]-(opts.bifactor?gLoad*gLoad:0)));y+=uniq*gauss();
      let cat=1;for(let c=0;c<K-1;c++)if(y>THR[c])cat=c+2;row[j]=cat;}
    rows.push(row);}
  return {rows,J,load,fac};}
function inject(rows,J,pct){const N=rows.length,out=rows.map(r=>r.slice()),lab=new Array(N).fill(0);
  const M=Math.round(N*pct);const ord=[...Array(N).keys()];for(let i=0;i<M;i++){const j=i+Math.floor(rng()*(N-i));const t=ord[i];ord[i]=ord[j];ord[j]=t;}
  const pats=["random","longstring","pure_straight","acquiescent","mixed","fatigue"];
  for(let m=0;m<M;m++){const i=ord[m],pat=pats[m%6],lvl=0.5+0.5*rng(),k=Math.max(1,Math.round(J*lvl)),row=out[i];
    const idx=[...Array(J).keys()];for(let t=0;t<k;t++){const j=t+Math.floor(rng()*(J-t));const tm=idx[t];idx[t]=idx[j];idx[j]=tm;}
    const sel=idx.slice(0,k);
    if(pat==="pure_straight"||pat==="longstring"){const v=1+Math.floor(rng()*K);sel.forEach(p=>row[p]=v);}
    else if(pat==="acquiescent"){sel.forEach(p=>row[p]=Math.max(1,Math.min(K,4+Math.round(rng()))));}
    else if(pat==="fatigue"){for(let t=J-k;t<J;t++)row[t]=1+Math.floor(rng()*K);}
    else sel.forEach(p=>row[p]=1+Math.floor(rng()*K)); // random, mixed ~ random
    lab[i]=1;}
  return {out,lab};}
const toMat=(rows,J)=>{const n=rows.length,m=new Float64Array(n*J);for(let i=0;i<n;i++)for(let j=0;j<J;j++)m[i*J+j]=rows[i][j];return m;};
const auc=(s,y)=>{let po=[],ne=[];for(let i=0;i<y.length;i++)(y[i]?po:ne).push(s[i]);let c=0;for(const a of po)for(const b of ne)c+=a>b?1:a===b?0.5:0;return c/(po.length*ne.length||1);};
function gateOK(gm){const p=Math.sqrt((gm.v1+gm.v2)/2)||1;return (gm.bic1-gm.bic2>2)&&gm.pi>0.01&&gm.pi<0.7&&(gm.m2-gm.m1)/p>1.0;}
function etaOf(mat,n,J){const ef=R.ensembleFeatures(mat,n,J,{iterations:150,seed:1}),O=ef.oriented,g=Math.max(0,Math.min(1,(ef.profileInfo-0.01)/0.05));
  const e=new Array(n);for(let i=0;i<n;i++)e[i]=W.b0+W.rr*O.rr[i]+g*(W.longstring*O.longstring[i]+W.person_total*O.person_total[i]);return e;}

const outC=[["condition","rep","rr_auc","ens_auc"]];
const conds={normal:{},low_load:{loRange:[0.2,0.5]},bifactor:{bifactor:true}};
for(const [cn,opts] of Object.entries(conds)){for(let rep=0;rep<8;rep++){
  const {rows,J}=genClean(15,6,500,opts);const {out,lab}=inject(rows,J,0.2);const mat=toMat(out,J);
  const ef=R.ensembleFeatures(mat,500,J,{iterations:100,seed:1});const e=etaOf(mat,500,J);
  outC.push([cn,rep,auc(ef.oriented.rr,lab).toFixed(3),auc(e,lab).toFixed(3)]);}}
fs.writeFileSync("sens_adverse.csv",outC.map(r=>r.join(",")).join("\n"));console.log("C done");
const outD=[["true_prev","rep","pi_hat","gate_fired"]];
for(const pct of [0.05,0.1,0.15,0.2,0.25,0.3]){for(let rep=0;rep<20;rep++){
  const {rows,J}=genClean(15,6,500);const {out}=inject(rows,J,pct);const e=etaOf(toMat(out,J),500,J);
  const gm=R.gaussMix1D(e),ok=gateOK(gm);outD.push([pct,rep,(ok?gm.pi:0).toFixed(4),ok?1:0]);}}
fs.writeFileSync("sens_pi.csv",outD.map(r=>r.join(",")).join("\n"));console.log("D done");
