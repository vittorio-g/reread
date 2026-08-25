/* eval_A_B_study13.js — evaluate, PER REPLICATE, the two proposed fixes:
 *   A  regime test : verdict "collapsed" iff pi_rc > pi_eta. Truth: is the true rate > 25%?
 *                    -> sensitivity / specificity, not just a mean curve.
 *   B' midpoint    : prevalence estimate = mean of the AVAILABLE component estimates
 *                    (pi_rc, pi_PT, counting only those whose mixture gate fired),
 *                    compared against the shipped pi_eta.
 * Run on:  (1) STUDY 1  — 70 real attentive + k real careless, prevalence swept
 *          (2) STUDY 3  — simulated multi-factor Likert, prevalence x length swept,
 *                         all six careless patterns (as in the paper's simulation).
 * Output: eval_A_B.csv
 */
const R=require("./site/rerere.js"); const W=R.WEIGHTS; const fs=require("fs");
let SEED=777; const rng=()=>{SEED=(SEED*1103515245+12345)&0x7fffffff;return SEED/0x7fffffff;};
function gauss(){let u=0,v=0;while(!u)u=rng();while(!v)v=rng();return Math.sqrt(-2*Math.log(u))*Math.cos(2*Math.PI*v);}
const K=5, THR=[-0.84,-0.25,0.25,0.84];
/* ---- Study 3 generators (same design as the paper's simulation) ---- */
function genClean(F,ipf,n){const J=F*ipf,lo=0.4,hi=0.8;const load=[],fac=[],tau=[];
  for(let j=0;j<J;j++){let l=lo+(hi-lo)*rng();if(rng()<0.4)l=-l;load.push(l);fac.push(Math.floor(j/ipf));tau.push(0.45*gauss());}
  const rho=0.3,sq=Math.sqrt(rho),sq1=Math.sqrt(1-rho);const rows=[];
  for(let i=0;i<n;i++){const common=gauss(),f=[];for(let k=0;k<F;k++)f.push(sq*common+sq1*gauss());
    const row=new Array(J);
    for(let j=0;j<J;j++){let y=tau[j]+load[j]*f[fac[j]];const uniq=Math.sqrt(Math.max(0.05,1-load[j]*load[j]));y+=uniq*gauss();
      let cat=1;for(let c=0;c<K-1;c++)if(y>THR[c])cat=c+2;row[j]=cat;}
    rows.push(row);}
  return {rows,J};}
function inject(rows,J,pct){const N=rows.length,out=rows.map(r=>r.slice()),lab=new Array(N).fill(0);
  const M=Math.round(N*pct);const ord=[...Array(N).keys()];for(let i=0;i<M;i++){const j=i+Math.floor(rng()*(N-i));const t=ord[i];ord[i]=ord[j];ord[j]=t;}
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
function piOf(a){const gm=R.gaussMix1D(a);const ok=gateOK(gm);return {pi:ok?gm.pi:0,ok};}
/* one replicate -> {verdictCollapsed, midpoint (or null), piEta} */
function oneRep(mat,n,J,iters){
  const ef=R.ensembleFeatures(mat,n,J,{iterations:iters,seed:1}),O=ef.oriented;
  const g=Math.max(0,Math.min(1,(ef.profileInfo-0.01)/0.05));
  const eta=new Array(n);
  for(let i=0;i<n;i++)eta[i]=W.b0+W.rr*O.rr[i]+g*(W.longstring*O.longstring[i]+W.person_total*O.person_total[i]);
  const E=piOf(eta),Rc=piOf(Array.from(O.rr)),Pt=piOf(Array.from(O.person_total));
  const avail=[]; if(Rc.ok)avail.push(Rc.pi); if(Pt.ok)avail.push(Pt.pi);
  const mid=avail.length?avail.reduce((s,v)=>s+v,0)/avail.length:null;
  return {collapsed:(Rc.pi>E.pi), mid, piEta:E.pi};
}
const mn=a=>a.length?a.reduce((s,v)=>s+v,0)/a.length:NaN;
const out=[["study","J","true_rate","reps","A_frac_collapsed","A_correct","B_mid_avail","B_abs_err_mid","abs_err_eta"]];
function report(study,J,tr,res,reps){
  const truthHigh = tr>0.25;
  const Acorrect = mn(res.map(r=>(r.collapsed===truthHigh)?1:0));
  const mids=res.filter(r=>r.mid!==null);
  const eMid=mn(mids.map(r=>Math.abs(r.mid-tr))), eEta=mn(res.map(r=>Math.abs(r.piEta-tr)));
  out.push([study,J,tr.toFixed(4),reps,mn(res.map(r=>r.collapsed?1:0)).toFixed(3),Acorrect.toFixed(3),
            (mids.length/reps).toFixed(3),eMid.toFixed(4),eEta.toFixed(4)]);
  console.log(study.padEnd(9)+"J="+String(J).padStart(3)+" true="+(100*tr).toFixed(0).padStart(3)+"% | A:collapsed "+
    (100*mn(res.map(r=>r.collapsed?1:0))).toFixed(0).padStart(3)+"%  correct "+(100*Acorrect).toFixed(0).padStart(3)+
    "% | B': avail "+(100*mids.length/reps).toFixed(0).padStart(3)+"%  |err| "+(100*eMid).toFixed(1).padStart(5)+
    "   (eta "+(100*eEta).toFixed(1).padStart(5)+")");
}
/* ================= STUDY 1 ================= */
console.log("=== STUDY 1 (real attentive + real careless, prevalence swept) ===");
{ const Ms=R.parseCSV(fs.readFileSync("_study_matrix.csv","utf8")),J=Ms.header.length;
  const allR=Ms.rows.map(r=>r.map(v=>{const x=Number(v);return Number.isFinite(x)?x:NaN;}));
  const y=R.parseCSV(fs.readFileSync("_study_labels.csv","utf8")).rows.map(r=>Number(r[0]));
  const cf=y.map((v,i)=>v===0?i:-1).filter(i=>i>=0), cl=y.map((v,i)=>v===1?i:-1).filter(i=>i>=0);
  const REPS=60;
  for(const k of [4,8,12,18,24,33,47,58]){ if(k>cl.length) continue; const tr=k/(cf.length+k); const res=[];
    for(let rep=0;rep<REPS;rep++){ const a=cl.slice();
      for(let i=0;i<k;i++){const j=i+Math.floor(rng()*(a.length-i));const t=a[i];a[i]=a[j];a[j]=t;}
      const idx=[...cf,...a.slice(0,k)],n=idx.length,mat=new Float64Array(n*J);
      for(let i=0;i<n;i++)for(let j=0;j<J;j++)mat[i*J+j]=allR[idx[i]][j];
      res.push(oneRep(mat,n,J,150)); }
    report("Study1",J,tr,res,REPS); } }
/* ================= STUDY 3 ================= */
console.log("\n=== STUDY 3 (simulated, six careless patterns, prevalence x length) ===");
{ const REPS=12, NSUB=400;
  for(const [F,ipf] of [[10,6],[30,6]]){ const Jn=F*ipf;
    for(const pct of [0.05,0.10,0.15,0.20,0.25,0.30,0.35,0.40,0.45]){ const res=[];
      for(let rep=0;rep<REPS;rep++){ const {rows,J}=genClean(F,ipf,NSUB); const {out:o,lab}=inject(rows,J,pct);
        res.push(oneRep(toMat(o,J),NSUB,J,100)); }
      report("Study3",Jn,pct,res,REPS); } } }
fs.writeFileSync("eval_A_B.csv",out.map(r=>r.join(",")).join("\n"));
console.log("\nwrote eval_A_B.csv");
