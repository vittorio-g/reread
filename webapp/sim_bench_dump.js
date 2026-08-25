/* sim_bench_dump.js — dump simulated datasets (matrix + labels + factor/reverse structure +
 * shipped rr/eta scores) across questionnaire lengths, for the Study-3 competitive benchmark in R.
 * Reuses the exact generative model of sim_shipped.js. */
const R=require("./site/reread.js"); const W=R.WEIGHTS; const fs=require("fs"); const rng=R._rng(20260709);
const OUT="sim_bench"; if(!fs.existsSync(OUT)) fs.mkdirSync(OUT);
function gauss(){let u=0,v=0;while(u===0)u=rng();while(v===0)v=rng();return Math.sqrt(-2*Math.log(u))*Math.cos(2*Math.PI*v);}
function qnorm(p){const a=[-3.969683028665376e+01,2.209460984245205e+02,-2.759285104469687e+02,1.383577518672690e+02,-3.066479806614716e+01,2.506628277459239e+00];
  const b=[-5.447609879822406e+01,1.615858368580409e+02,-1.556989798598866e+02,6.680131188771972e+01,-1.328068155288572e+01];
  const c=[-7.784894002430293e-03,-3.223964580411365e-01,-2.400758277161838e+00,-2.549732539343734e+00,4.374664141464968e+00,2.938163982698783e+00];
  const d=[7.784695709041462e-03,3.224671290700398e-01,2.445134137142996e+00,3.754408661907416e+00];const pl=0.02425;let q,r;
  if(p<pl){q=Math.sqrt(-2*Math.log(p));return (((((c[0]*q+c[1])*q+c[2])*q+c[3])*q+c[4])*q+c[5])/((((d[0]*q+d[1])*q+d[2])*q+d[3])*q+1);}
  if(p<=1-pl){q=p-0.5;r=q*q;return (((((a[0]*r+a[1])*r+a[2])*r+a[3])*r+a[4])*r+a[5])*q/(((((b[0]*r+b[1])*r+b[2])*r+b[3])*r+b[4])*r+1);}
  q=Math.sqrt(-2*Math.log(1-p));return -(((((c[0]*q+c[1])*q+c[2])*q+c[3])*q+c[4])*q+c[5])/((((d[0]*q+d[1])*q+d[2])*q+d[3])*q+1);}
const K=5, THR=[]; for(let c=1;c<K;c++) THR.push(qnorm(c/K));
function genClean(F,ipf,n){const J=F*ipf,load=[],fac=[],tau=[],rev=[];
  for(let j=0;j<J;j++){let l=0.4+0.4*rng();if(rng()<0.4)l=-l;load.push(l);fac.push(Math.floor(j/ipf));tau.push(0.45*gauss());rev.push(l<0?1:0);}
  const rho=0.3,sq=Math.sqrt(rho),sq1=Math.sqrt(1-rho),rows=[];
  for(let i=0;i<n;i++){const common=gauss(),f=[];for(let k=0;k<F;k++)f.push(sq*common+sq1*gauss());const row=new Array(J);
    for(let j=0;j<J;j++){const y=tau[j]+load[j]*f[fac[j]]+Math.sqrt(1-load[j]*load[j])*gauss();let cat=1;for(let c=0;c<K-1;c++)if(y>THR[c])cat=c+2;row[j]=cat;}
    rows.push(row);}
  return {rows,J,fac,rev};}
function inject(rows,J,pct){const N=rows.length;let lo=Infinity,hi=-Infinity;for(let i=0;i<N;i++)for(let j=0;j<J;j++){const v=rows[i][j];if(v<lo)lo=v;if(v>hi)hi=v;}
  const LEVELS=[];for(let v=lo;v<=hi;v++)LEVELS.push(v);const highLo=lo+0.65*(hi-lo),highLevels=LEVELS.filter(v=>v>=highLo);
  const ri=n=>Math.floor(rng()*n),pick=a=>a[ri(a.length)];
  const sampleK=(pool,k)=>{const a=pool.slice();for(let i=0;i<k;i++){const j=i+ri(a.length-i);const t=a[i];a[i]=a[j];a[j]=t;}return a.slice(0,k);};
  const allIdx=[...Array(J).keys()],clamp=x=>Math.max(lo,Math.min(hi,x));
  const chunks=(row,idx,k,seq)=>{const nc=k>=4?2+ri(Math.min(5,k)-1):(k>=2?2:1);const asg=[];for(let t=0;t<k;t++)asg.push(t%nc);
    if(!seq)for(let t=k-1;t>0;t--){const j=ri(t+1);const tmp=asg[t];asg[t]=asg[j];asg[j]=tmp;}
    for(let ch=0;ch<nc;ch++){const val=pick(LEVELS);for(let t=0;t<k;t++)if(asg[t]===ch)row[idx[t]]=val;}};
  const patterns=["random","longstring","pure_straight","acquiescent","mixed","fatigue"],levels=[0.5,0.6,0.7,0.8,0.9,1.0];
  const M=Math.round(N*pct),order=sampleK([...Array(N).keys()],N),ids=order.slice(0,M);
  const labels=new Array(N).fill(0),sev=new Array(N).fill(0),out=rows.map(r=>r.slice());
  for(let m=0;m<M;m++){const i=ids[m],pat=patterns[m%patterns.length],lvl=pick(levels),k=Math.min(J,Math.max(1,Math.round(J*lvl))),row=out[i];
    if(pat==="random")sampleK(allIdx,k).forEach(p=>row[p]=pick(LEVELS));
    else if(pat==="longstring")chunks(row,sampleK(allIdx,k),k,false);
    else if(pat==="pure_straight"){const v=pick(LEVELS);sampleK(allIdx,k).forEach(p=>row[p]=v);}
    else if(pat==="acquiescent"){const b=pick(highLevels);sampleK(allIdx,k).forEach(p=>row[p]=clamp(b+pick([-1,0,0,0,1])));}
    else if(pat==="mixed"){const kl=Math.max(1,Math.round(k/2)),idx=sampleK(allIdx,k);chunks(row,idx.slice(0,kl),kl,false);idx.slice(kl).forEach(p=>row[p]=pick(LEVELS));}
    else if(pat==="fatigue"){const idx=[];for(let t=J-k;t<J;t++)idx.push(t);chunks(row,idx,k,true);}
    labels[i]=1;sev[i]=lvl;}
  return {out,labels,sev};}

const SIZES=[[5,6],[10,6],[20,6],[30,6],[40,6],[50,6]];  // 30,60,120,180,240,300
const N=500, RATE=0.20, REPS=10;
const idx=[["file","items","rep","n","rate"]];
for(const [F,ipf] of SIZES){const J=F*ipf;
  for(let rep=0;rep<REPS;rep++){
    const {rows:clean,fac,rev}=genClean(F,ipf,N);
    const {out,labels,sev}=inject(clean,J,RATE);
    const n=out.length,mat=new Float64Array(n*J);for(let i=0;i<n;i++)for(let j=0;j<J;j++)mat[i*J+j]=out[i][j];
    const ef=R.ensembleFeatures(mat,n,J,{seed:1,iterations:120}),O=ef.oriented,g=Math.max(0,Math.min(1,(ef.profileInfo-0.01)/0.05));
    const tag=J+"_"+rep;
    const Mo=[];for(let i=0;i<n;i++)Mo.push(out[i].join(","));
    const Yo=["y,rr,eta,sev"];for(let i=0;i<n;i++){const eta=W.b0+W.rr*O.rr[i]+g*(W.longstring*O.longstring[i]+W.person_total*O.person_total[i]);
      Yo.push(labels[i]+","+O.rr[i].toFixed(6)+","+eta.toFixed(6)+","+sev[i].toFixed(3));}
    const Me=["factor,reverse"];for(let j=0;j<J;j++)Me.push(fac[j]+","+rev[j]);
    fs.writeFileSync(OUT+"/"+tag+"_M.csv",Mo.join("\n"));
    fs.writeFileSync(OUT+"/"+tag+"_y.csv",Yo.join("\n"));
    fs.writeFileSync(OUT+"/"+tag+"_meta.csv",Me.join("\n"));
    idx.push([tag,J,rep,n,(labels.reduce((s,v)=>s+v,0)/n).toFixed(3)]);
  }
  console.log("dumped J="+J+" ("+REPS+" reps)");
}
fs.writeFileSync(OUT+"/index.csv",idx.map(r=>r.join(",")).join("\n"));
console.log("wrote "+OUT+"/* ("+(idx.length-1)+" datasets)");
