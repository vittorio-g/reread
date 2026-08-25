/* sim_v2.js — realistic careless definition + flexible flagging.
 * GT: careless(label 1) = junk patterns {random,mixed,longstring,pure_straight} with >50% corruption.
 * Distractors kept IN the data with label 0: acquiescent (response style) + fatigue (partial, <=50%).
 * Flagging methods compared: Standard(z2.5), High(z1.5), A(pi top-share), B(Bayes posterior>0.5).
 * References: AUC and oracle-MCC (ranking / ceiling, method-independent). */
const R=require("./site/rerere.js");const W=R.WEIGHTS;const fs=require("fs");
const rng=R._rng(20260713);
function gauss(){let u=0,v=0;while(u===0)u=rng();while(v===0)v=rng();return Math.sqrt(-2*Math.log(u))*Math.cos(2*Math.PI*v);}
function qnorm(p){const a=[-39.6968302,220.9460984,-275.9285104,138.3577518,-30.66479806,2.506628277];const b=[-54.47609879,161.5858368,-155.6989798,66.80131188,-13.28068155];const c=[-0.007784894002,-0.3223964580,-2.400758277,-2.549732539,4.374664141,2.938163982];const d=[0.007784695709,0.3224671290,2.445134137,3.754408661];const pl=0.02425;let q,r;if(p<pl){q=Math.sqrt(-2*Math.log(p));return(((((c[0]*q+c[1])*q+c[2])*q+c[3])*q+c[4])*q+c[5])/((((d[0]*q+d[1])*q+d[2])*q+d[3])*q+1);}if(p<=1-pl){q=p-0.5;r=q*q;return(((((a[0]*r+a[1])*r+a[2])*r+a[3])*r+a[4])*r+a[5])*q/(((((b[0]*r+b[1])*r+b[2])*r+b[3])*r+b[4])*r+1);}q=Math.sqrt(-2*Math.log(1-p));return-(((((c[0]*q+c[1])*q+c[2])*q+c[3])*q+c[4])*q+c[5])/((((d[0]*q+d[1])*q+d[2])*q+d[3])*q+1);}
function ncdf(x){const t=1/(1+0.2316419*Math.abs(x)),d=0.3989423*Math.exp(-x*x/2);const p=d*t*(0.3193815+t*(-0.3565638+t*(1.781478+t*(-1.821256+t*1.330274))));return x>0?1-p:p;}
const K=5, THR=[]; for(let c=1;c<K;c++) THR.push(qnorm(c/K));
function genClean(F,ipf,n){
  const J=F*ipf, load=[],fac=[],tau=[];
  for(let j=0;j<J;j++){let l=0.4+0.4*rng(); if(rng()<0.4) l=-l; load.push(l); fac.push(Math.floor(j/ipf)); tau.push(0.45*gauss());}
  const rho=0.3,sq=Math.sqrt(rho),sq1=Math.sqrt(1-rho);
  const rows=[];
  for(let i=0;i<n;i++){
    const common=gauss(),f=[];for(let k=0;k<F;k++)f.push(sq*common+sq1*gauss());
    const row=new Array(J);
    for(let j=0;j<J;j++){const y=tau[j]+load[j]*f[fac[j]]+Math.sqrt(1-load[j]*load[j])*gauss();let cat=1;for(let c=0;c<K-1;c++)if(y>THR[c])cat=c+2;row[j]=cat;}
    rows.push(row);
  }
  return {rows,J};
}
// apply one careless pattern to a row over k items
const LEVELS=[]; for(let v=1;v<=K;v++) LEVELS.push(v);
const highLevels=LEVELS.filter(v=>v>=1+0.65*(K-1));
const ri=n=>Math.floor(rng()*n), pick=a=>a[ri(a.length)], clamp=x=>Math.max(1,Math.min(K,x));
function sampleK(pool,k){const a=pool.slice();for(let i=0;i<k;i++){const j=i+ri(a.length-i);const t=a[i];a[i]=a[j];a[j]=t;}return a.slice(0,k);}
function chunks(row,idx,k,seq){const nc=k>=4?2+ri(Math.min(5,k)-1):(k>=2?2:1);const asg=[];for(let t=0;t<k;t++)asg.push(t%nc);
  if(!seq)for(let t=k-1;t>0;t--){const j=ri(t+1);const tmp=asg[t];asg[t]=asg[j];asg[j]=tmp;}
  for(let ch=0;ch<nc;ch++){const val=pick(LEVELS);for(let t=0;t<k;t++)if(asg[t]===ch)row[idx[t]]=val;}}
function applyPattern(row,J,pat,lvl){
  const k=Math.min(J,Math.max(1,Math.round(J*lvl))), allIdx=[...Array(J).keys()];
  if(pat==="random")sampleK(allIdx,k).forEach(p=>row[p]=pick(LEVELS));
  else if(pat==="longstring")chunks(row,sampleK(allIdx,k),k,false);
  else if(pat==="pure_straight"){const v=pick(LEVELS);sampleK(allIdx,k).forEach(p=>row[p]=v);}
  else if(pat==="acquiescent"){const b=pick(highLevels);sampleK(allIdx,k).forEach(p=>row[p]=clamp(b+pick([-1,0,0,0,1])));}
  else if(pat==="mixed"){const kl=Math.max(1,Math.round(k/2)),idx=sampleK(allIdx,k);chunks(row,idx.slice(0,kl),kl,false);idx.slice(kl).forEach(p=>row[p]=pick(LEVELS));}
  else if(pat==="fatigue"){const idx=[];for(let t=J-k;t<J;t++)idx.push(t);chunks(row,idx,k,true);}
}
const JUNK=["random","mixed","longstring","pure_straight"], JUNK_CORR=[0.6,0.7,0.8,0.9,1.0];
const ACQ_CORR=[0.5,0.6,0.7,0.8,0.9,1.0], FAT_CORR=[0.3,0.4,0.5];
// inject: pCare true careless (junk >50% -> label 1); pDist distractors (acquiescent/fatigue -> label 0, kept)
function inject2(rows,J,pCare,pDist){
  const n=rows.length, out=rows.map(r=>r.slice()), labels=new Array(n).fill(0), dist=new Array(n).fill(false);
  const sev=new Array(n).fill(0);   // continuous "how clearly this should be excluded": clean 0 .. extreme junk 1
  const order=sampleK([...Array(n).keys()],n);
  const nCare=Math.round(n*pCare), nDist=Math.round(n*pDist);
  for(let m=0;m<nCare;m++){const i=order[m];const c=pick(JUNK_CORR);applyPattern(out[i],J,pick(JUNK),c);labels[i]=1;sev[i]=c;/* junk severity = corruption 0.6-1.0 */}
  for(let m=nCare;m<nCare+nDist;m++){const i=order[m];
    if(rng()<0.5){applyPattern(out[i],J,"acquiescent",pick(ACQ_CORR));sev[i]=0.30;/* response style: low junk severity */}
    else{const c=pick(FAT_CORR);applyPattern(out[i],J,"fatigue",c);sev[i]=c;/* partial junk 0.3-0.5 */}
    dist[i]=true; /* label stays 0 */}
  return {out,labels,dist,sev};
}
// ---------- metrics + flagging ----------
function mcc(f,y){let tp=0,fp=0,fn=0,tn=0;for(let i=0;i<y.length;i++){if(y[i]){f[i]?tp++:fn++;}else{f[i]?fp++:tn++;}}const d=Math.sqrt((tp+fp)*(tp+fn)*(tn+fp)*(tn+fn));return d?(tp*tn-fp*fn)/d:0;}
function recall(f,y){let tp=0,fn=0;for(let i=0;i<y.length;i++)if(y[i]){f[i]?tp++:fn++;}return tp+fn?tp/(tp+fn):0;}
function precision(f,y){let tp=0,fp=0;for(let i=0;i<y.length;i++){if(f[i]){y[i]?tp++:fp++;}}return tp+fp?tp/(tp+fp):(0);}
function fprOn(f,mask){let flagged=0,tot=0;for(let i=0;i<mask.length;i++)if(mask[i]){tot++;if(f[i])flagged++;}return tot?flagged/tot:0;}
function auc(s,y){let po=[],ne=[];for(let i=0;i<y.length;i++)(y[i]?po:ne).push(s[i]);let c=0;for(const a of po)for(const b of ne)c+=a>b?1:a===b?0.5:0;return c/(po.length*ne.length);}
function oracleMCC(eta,y){const s=[...eta].sort((a,b)=>a-b);let best=0;const step=Math.max(1,Math.floor(s.length/100));for(let t=0;t<s.length;t+=step){const thr=s[t];const f=eta.map(v=>v>=thr);const m=mcc(f,y);if(m>best)best=m;}return best;}
// U-shaped severity weight: heavy for clearly-careful (s~0) and clearly-careless (s~1), light for the ambiguous middle (s~0.5)
function uWeight(s){return (2*s-1)*(2*s-1);}
// weighted MCC: same as MCC but each respondent contributes its U-weight to the confusion cells
function wmcc(f,y,w){let tp=0,fp=0,fn=0,tn=0;for(let i=0;i<y.length;i++){const wi=w[i];if(y[i]){f[i]?tp+=wi:fn+=wi;}else{f[i]?fp+=wi:tn+=wi;}}
  const d=Math.sqrt((tp+fp)*(tp+fn)*(tn+fp)*(tn+fn));return d?(tp*tn-fp*fn)/d:0;}
// Spearman rank correlation between the detector score and the graded severity (index-quality reference)
function spearman(x,s){const n=x.length,rk=a=>{const idx=[...a.keys()].sort((p,q)=>a[p]-a[q]);const r=new Array(n);for(let i=0;i<n;i++)r[idx[i]]=i;return r;};
  const rx=rk(x),rs=rk(s);const mx=(n-1)/2;let num=0,dx=0,ds=0;for(let i=0;i<n;i++){const a=rx[i]-mx,b=rs[i]-mx;num+=a*b;dx+=a*a;ds+=b*b;}return num/Math.sqrt(dx*ds+1e-9);}
function gateOK(gm){const pooled=Math.sqrt((gm.v1+gm.v2)/2)||1,sep=(gm.m2-gm.m1)/pooled;return (gm.bic1-gm.bic2>2)&&gm.pi>0.01&&gm.pi<0.7&&sep>1.0;}
// 3-component 1-D Gaussian mixture (clean / response-style / careless), EM, sorted by mean asc
function gaussMix3(x){
  const n=x.length, s=[...x].sort((a,b)=>a-b), mean=x.reduce((a,b)=>a+b,0)/n;
  let varr=0;for(let i=0;i<n;i++)varr+=(x[i]-mean)**2;varr/=n;const vfloor=Math.max(1e-4,varr*0.02);
  let m=[s[Math.floor(0.15*n)],s[Math.floor(0.5*n)],s[Math.floor(0.88*n)]], v=[varr/3,varr/3,varr/3], p=[0.6,0.25,0.15];
  const norm=(x,mu,s2)=>Math.exp(-((x-mu)*(x-mu))/(2*s2))/Math.sqrt(2*Math.PI*s2);
  const r=[new Array(n),new Array(n),new Array(n)];
  for(let it=0;it<300;it++){
    const sp=[0,0,0];
    for(let i=0;i<n;i++){const a=[p[0]*norm(x[i],m[0],v[0]),p[1]*norm(x[i],m[1],v[1]),p[2]*norm(x[i],m[2],v[2])],z=a[0]+a[1]+a[2]+1e-300;for(let g=0;g<3;g++){r[g][i]=a[g]/z;sp[g]+=r[g][i];}}
    for(let g=0;g<3;g++)p[g]=Math.max(1e-3,sp[g]/n);
    for(let g=0;g<3;g++){let mu=0;for(let i=0;i<n;i++)mu+=r[g][i]*x[i];m[g]=mu/(sp[g]+1e-9);}
    for(let g=0;g<3;g++){let q=0;for(let i=0;i<n;i++)q+=r[g][i]*(x[i]-m[g])**2;v[g]=Math.max(vfloor,q/(sp[g]+1e-9));}
  }
  const idx=[0,1,2].sort((a,b)=>m[a]-m[b]);
  const comps=idx.map(g=>({m:m[g],v:v[g],p:p[g],r:r[g]}));
  let ll=0;for(let i=0;i<n;i++)ll+=Math.log(comps[0].p*norm(x[i],comps[0].m,comps[0].v)+comps[1].p*norm(x[i],comps[1].m,comps[1].v)+comps[2].p*norm(x[i],comps[2].m,comps[2].v)+1e-300);
  return {comps, bic3:-2*ll+8*Math.log(n), post3:comps[2].r};
}
function allFlags(eta){
  const gm=R.gaussMix1D(eta),n=eta.length,ok=gateOK(gm),s1=Math.sqrt(gm.v1);
  const off=eta.map(_=>false);
  if(!ok) return {std:off,high:off,A:off,B:off,C:off,D:off,E:off};
  const cutStd=gm.m1+2.5*s1, cutHigh=gm.m1+1.5*s1;
  const srt=[...eta].sort((a,b)=>a-b), cutA=srt[Math.min(n-1,Math.floor((1-gm.pi)*n))];
  // C: Benjamini-Hochberg FDR (q=0.10) on upper-tail p-values vs the attentive component
  const pv=eta.map(v=>1-ncdf((v-gm.m1)/s1)), idx=[...pv.keys()].sort((a,b)=>pv[a]-pv[b]);
  let kmax=-1; for(let k=0;k<n;k++) if(pv[idx[k]]<=((k+1)/n)*0.10) kmax=k;
  const C=eta.map(_=>false); for(let k=0;k<=kmax;k++) C[idx[k]]=true;
  // D: adaptive z on the attentive SD — lower z (more aggressive) when the two clusters separate,
  //    which happens as questionnaire length grows -> the cut "scales" with items.
  const pooled=Math.sqrt((gm.v1+gm.v2)/2)||1, sep=(gm.m2-gm.m1)/pooled;
  const zD=Math.max(1.2,Math.min(2.5,3.0-0.45*sep)), cutD=gm.m1+zD*s1;
  // E: 3-component model — flag only the TOP (careless) component; the middle (response-style)
  //    component is modelled explicitly so it no longer drags the cut. Fall back to Standard's
  //    cut if the top component is not distinct from the middle (e.g. no real 3rd group).
  const g3=gaussMix3(eta), c=g3.comps, sep32=(c[2].m-c[1].m)/Math.sqrt((c[1].v+c[2].v)/2);
  // flag if the TOP (careless) component is the most likely one (argmax over the 3 posteriors)
  const Eflag=eta.map((_,i)=> c[2].r[i]>=c[1].r[i] && c[2].r[i]>=c[0].r[i]);
  const E=(sep32>=0.8 && c[2].p>=0.01 && c[2].p<=0.6) ? Eflag : eta.map(v=>v>cutStd);
  return {
    std: eta.map(v=>v>cutStd),
    high: eta.map(v=>v>cutHigh),
    A: eta.map(v=>v>=cutA),
    B: gm.post.map(p=>p>0.5),
    C: C,
    D: eta.map(v=>v>cutD),
    E: E
  };
}
function scoreDataset(out,J,labels,dist,sev){
  const n=out.length,mat=new Float64Array(n*J);for(let i=0;i<n;i++)for(let j=0;j<J;j++)mat[i*J+j]=out[i][j];
  const ef=R.ensembleFeatures(mat,n,J,{seed:1,iterations:150}),O=ef.oriented;
  const g=Math.max(0,Math.min(1,(ef.profileInfo-0.01)/0.05));
  const eta=new Array(n);for(let i=0;i<n;i++)eta[i]=W.b0+W.rr*O.rr[i]+g*(W.longstring*O.longstring[i]+W.person_total*O.person_total[i]);
  const w=sev.map(uWeight);
  const F=allFlags(eta), r={auc:auc(eta,labels),oracle:oracleMCC(eta,labels),spearman:spearman(eta,sev)};
  for(const m of ["std","high","A","B","C","D","E"]){
    r[m+"_mcc"]=mcc(F[m],labels); r[m+"_wmcc"]=wmcc(F[m],labels,w); r[m+"_rec"]=recall(F[m],labels);
    r[m+"_prec"]=precision(F[m],labels); r[m+"_est"]=F[m].filter(Boolean).length/n; r[m+"_fprD"]=fprOn(F[m],dist);
  }
  return r;
}
const mn=a=>a.reduce((s,v)=>s+v,0)/a.length, se=a=>{const m=mn(a);return Math.sqrt(a.reduce((s,v)=>s+(v-m)**2,0)/(a.length-1))/Math.sqrt(a.length);};
const METHODS=["std","high","A","B","C","D","E"];
const STATS=["mcc","wmcc","rec","prec","est","fprD"];
function aggRow(key,val,accs){
  const row={[key]:val,auc:mn(accs.auc).toFixed(4),oracle:mn(accs.oracle).toFixed(4),spearman:mn(accs.spearman).toFixed(4)};
  for(const m of METHODS){for(const s of STATS){row[m+"_"+s]=mn(accs[m+"_"+s]).toFixed(4);}row[m+"_mcc_se"]=se(accs[m+"_mcc"]).toFixed(4);row[m+"_wmcc_se"]=se(accs[m+"_wmcc"]).toFixed(4);}
  return row;
}
function newAcc(){const a={auc:[],oracle:[],spearman:[]};for(const m of METHODS)for(const s of STATS)a[m+"_"+s]=[];return a;}
function pushAcc(a,r){a.auc.push(r.auc);a.oracle.push(r.oracle);a.spearman.push(r.spearman);for(const m of METHODS)for(const s of STATS)a[m+"_"+s].push(r[m+"_"+s]);}
function writeCSV(fn,rows){const cols=Object.keys(rows[0]);fs.writeFileSync(fn,cols.join(",")+"\n"+rows.map(r=>cols.map(c=>r[c]).join(",")).join("\n"));}

// ================= SWEEP 1: length =================
const N=600, PDIST=0.15, REPS=40;
console.log("=== SCALING: rate=0.15 careless + 0.15 distractors, "+REPS+" reps, n="+N+" ===");
console.log("items  spearman | weighted-MCC per method: "+METHODS.join(" "));
const scaleRows=[];
for(const [F,ipf] of [[5,6],[10,6],[20,6],[30,6],[40,6],[50,6]]){  // 30,60,120,180,240,300
  const J=F*ipf, acc=newAcc();
  for(let rep=0;rep<REPS;rep++){const {rows}=genClean(F,ipf,N);const {out,labels,dist,sev}=inject2(rows,J,0.15,PDIST);pushAcc(acc,scoreDataset(out,J,labels,dist,sev));}
  const row=aggRow("items",J,acc); scaleRows.push(row);
  console.log(String(J).padStart(4)+"    "+mn(acc.spearman).toFixed(2)+"   | "+METHODS.map(m=>mn(acc[m+"_wmcc"]).toFixed(2)).join("  "));
}
writeCSV("sim_v2_scale.csv",scaleRows);

// ================= SWEEP 2: prevalence (fixed length 120) =================
console.log("\n=== PREVALENCE: 120 items, vary true careless rate, +0.15 distractors, "+REPS+" reps ===");
console.log("trueRate | weighted-MCC per method: "+METHODS.join(" "));
const prevRows=[];
for(const p of [0.05,0.10,0.15,0.20,0.25,0.30]){
  const F=20,ipf=6,J=120, acc=newAcc();
  for(let rep=0;rep<REPS;rep++){const {rows}=genClean(F,ipf,N);const {out,labels,dist,sev}=inject2(rows,J,p,PDIST);pushAcc(acc,scoreDataset(out,J,labels,dist,sev));}
  const row=aggRow("true_rate",p,acc); prevRows.push(row);
  console.log(String(p.toFixed(2)).padStart(7)+"  | "+METHODS.map(m=>mn(acc[m+"_wmcc"]).toFixed(2)).join("  "));
}
writeCSV("sim_v2_prev.csv",prevRows);
console.log("\nwrote sim_v2_scale.csv, sim_v2_prev.csv");
