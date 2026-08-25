/* test_cut.js — does using pi (mixture prevalence) to PLACE the cut beat the fixed z-cut?
 * Metric = MCC vs true labels (that's what matters: where to cut). Real study sweep, 200 resamples/rate.
 * Strategies: Standard(z2.5), High(z1.5), pi-cut(top-pi by eta), leak-cut(top-leak), oracle(best thr). */
const R=require("./site/rerere.js"); const W=R.WEIGHTS; const fs=require("fs"); const rng=R._rng(31459);
const Ms=R.parseCSV(fs.readFileSync("_study_matrix.csv","utf8")), J=Ms.header.length;
const allR=Ms.rows.map(r=>r.map(v=>{const x=Number(v);return Number.isFinite(x)?x:NaN;}));
const y=R.parseCSV(fs.readFileSync("_study_labels.csv","utf8")).rows.map(r=>Number(r[0]));
const cf=y.map((v,i)=>v===0?i:-1).filter(i=>i>=0), cl=y.map((v,i)=>v===1?i:-1).filter(i=>i>=0);
function sampleK(pool,k){const a=pool.slice();for(let i=0;i<k;i++){const j=i+Math.floor(rng()*(a.length-i));const t=a[i];a[i]=a[j];a[j]=t;}return a.slice(0,k);}
function gateOK(gm){const pooled=Math.sqrt((gm.v1+gm.v2)/2)||1,sep=(gm.m2-gm.m1)/pooled;return (gm.bic1-gm.bic2>2)&&gm.pi>0.01&&gm.pi<0.7&&sep>1.0;}
const ncdf=x=>{const t=1/(1+0.2316419*Math.abs(x)),d=0.3989423*Math.exp(-x*x/2);const p=d*t*(0.3193815+t*(-0.3565638+t*(1.781478+t*(-1.821256+t*1.330274))));return x>0?1-p:p;};
function eta(idx){const n=idx.length,mat=new Float64Array(n*J);for(let i=0;i<n;i++)for(let j=0;j<J;j++)mat[i*J+j]=allR[idx[i]][j];
  const ef=R.ensembleFeatures(mat,n,J,{seed:1,iterations:200}),O=ef.oriented,g=Math.max(0,Math.min(1,(ef.profileInfo-0.01)/0.05));
  const e=new Array(n);for(let i=0;i<n;i++)e[i]=W.b0+W.rr*O.rr[i]+g*(W.longstring*O.longstring[i]+W.person_total*O.person_total[i]);return e;}
function mcc(flag,lab){let tp=0,fp=0,tn=0,fn=0;for(let i=0;i<lab.length;i++){if(flag[i]){lab[i]?tp++:fp++;}else{lab[i]?fn++:tn++;}}
  const d=Math.sqrt((tp+fp)*(tp+fn)*(tn+fp)*(tn+fn));return d===0?0:(tp*tn-fp*fn)/d;}
function topK(e,k){const idx=e.map((v,i)=>[v,i]).sort((a,b)=>b[0]-a[0]);const f=new Array(e.length).fill(false);for(let i=0;i<k;i++)f[idx[i][1]]=true;return f;}
function oracleMCC(e,lab){const th=[...e].sort((a,b)=>a-b);let best=-1;for(let t=0;t<th.length;t++){const f=e.map(v=>v>=th[t]);const m=mcc(f,lab);if(m>best)best=m;}return best;}
function strat(e,lab){
  const gm=R.gaussMix1D(e),ok=gateOK(gm),n=e.length,s1=Math.sqrt(gm.v1),s2=Math.sqrt(gm.v2);
  const std=ok?e.map(v=>v>gm.m1+2.5*s1):e.map(()=>false);
  const high=ok?e.map(v=>v>gm.m1+1.5*s1):e.map(()=>false);
  const pi=ok?gm.pi:0; let leak=pi;
  if(ok){const cut=gm.m1+2.5*s1,Satt=1-ncdf((cut-gm.m1)/s1),Scare=1-ncdf((cut-gm.m2)/s2),den=Scare-Satt,obs=std.filter(Boolean).length/n;
    leak=Math.abs(den)>0.05?Math.max(0,Math.min(1,(obs-Satt)/den)):pi;}
  const picut=ok?topK(e,Math.round(pi*n)):e.map(()=>false);
  const leakcut=ok?topK(e,Math.round(leak*n)):e.map(()=>false);
  return {std:mcc(std,lab),high:mcc(high,lab),picut:mcc(picut,lab),leakcut:mcc(leakcut,lab),oracle:oracleMCC(e,lab)};
}
const mn=a=>a.reduce((s,v)=>s+v,0)/a.length, NAMES=["std","high","picut","leakcut","oracle"], NC=cf.length, REPS=200;
console.log("REAL study sweep — MCC by cut strategy ("+NC+" careful + up to "+cl.length+" careless, "+REPS+" resamples/rate)\n");
console.log("true |  Standard   High     pi-cut   leak-cut | oracle | best-auto");
const out=[["true",...NAMES]];
for(const p of [0.05,0.10,0.15,0.20,0.25,0.30,0.40]){
  const k=Math.round(NC*p/(1-p)); if(k>cl.length) continue; const tr=k/(NC+k);
  const acc={}; for(const m of NAMES) acc[m]=[];
  for(let rep=0;rep<REPS;rep++){const idx=[...cf,...sampleK(cl,k)];const lab=idx.map(i=>y[i]);const s=strat(eta(idx),lab);for(const m of NAMES)acc[m].push(s[m]);}
  const M={}; for(const m of NAMES) M[m]=mn(acc[m]);
  const autos=["std","high","picut","leakcut"], bestAuto=autos.reduce((a,b)=>M[b]>M[a]?b:a,"std");
  out.push([tr.toFixed(3),...NAMES.map(m=>M[m].toFixed(3))]);
  console.log((100*tr).toFixed(0).padStart(4)+"% |  "+["std","high","picut","leakcut"].map(m=>M[m].toFixed(3)).join("    ")+"  |  "+M.oracle.toFixed(3)+" |  "+bestAuto+" ("+M[bestAuto].toFixed(3)+")");
}
fs.writeFileSync("test_cut_study.csv",out.map(r=>r.join(",")).join("\n"));
console.log("\nwrote test_cut_study.csv");
