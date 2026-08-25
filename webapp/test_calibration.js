/* test_calibration.js — is there a GENUINELY calibrated rate estimator (low bias + good coverage),
 * or are we stuck with two biased estimators? Real study data (careful base + sampled real careless),
 * matching Figure 3. For each true rate, 200 resamples; report mean, bias, and 90% coverage for:
 *   - Standard flagged fraction (current)   - High flagged fraction (current)
 *   - avg(Std,High)                         - mixture pi (EM prevalence MLE)
 *   - leakage-corrected flagged fraction (invert the mixture)                                     */
const R=require("./site/reread.js"); const W=R.WEIGHTS; const fs=require("fs"); const rng=R._rng(31459);
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
function estimators(e){
  const gm=R.gaussMix1D(e),ok=gateOK(gm),n=e.length;
  const std=ok?R.autoFlag(e,"standard").estRate:0, high=ok?R.autoFlag(e,"high").estRate:0;
  const pi=ok?gm.pi:0;
  // leakage-corrected: observed flagged (Standard) = pi*S_care(c)+(1-pi)*S_att(c); solve pi given cut & components
  let leak=std;
  if(ok){const s1=Math.sqrt(gm.v1),s2=Math.sqrt(gm.v2),cut=gm.m1+2.5*s1;
    const Satt=1-ncdf((cut-gm.m1)/s1),Scare=1-ncdf((cut-gm.m2)/s2),denom=Scare-Satt;
    leak=Math.abs(denom)>0.05?Math.max(0,Math.min(1,(std-Satt)/denom)):pi;}
  return {std,high,avg:(std+high)/2,pi,leak};
}
const mn=a=>a.reduce((s,v)=>s+v,0)/a.length, pct=(a,q)=>{const s=[...a].sort((x,y)=>x-y);return s[Math.min(s.length-1,Math.floor(q*s.length))];};
const NAMES=["std","high","avg","pi","leak"], NC=cf.length, REPS=200;
console.log("REAL study: "+NC+" careful + "+cl.length+" careless. 90% coverage in [] (true inside estimator's 5-95 pct).\n");
console.log("true |            std                high               avg                pi                 leak");
const out=[["true",...NAMES.flatMap(m=>[m+"_mean",m+"_bias",m+"_cov"])]];
for(const p of [0.05,0.10,0.15,0.20,0.25,0.30,0.40]){
  const k=Math.round(NC*p/(1-p)); if(k>cl.length) continue; const tr=k/(NC+k);
  const acc={}; for(const m of NAMES) acc[m]=[];
  for(let rep=0;rep<REPS;rep++){const idx=[...cf,...sampleK(cl,k)];const est=estimators(eta(idx));for(const m of NAMES)acc[m].push(est[m]);}
  const row=[tr.toFixed(3)]; let line=(100*tr).toFixed(0).padStart(4)+"% | ";
  for(const m of NAMES){const me=mn(acc[m]),bias=me-tr,lo=pct(acc[m],0.05),hi=pct(acc[m],0.95),cov=(tr>=lo&&tr<=hi)?"Y":"n";
    row.push(me.toFixed(4),bias.toFixed(4),cov);
    line+=(100*me).toFixed(0).padStart(3)+"%("+(bias>=0?"+":"")+(100*bias).toFixed(0)+",cov"+cov+")  ";}
  out.push(row); console.log(line);
}
fs.writeFileSync("test_calibration.csv",out.map(r=>r.join(",")).join("\n"));
console.log("\n(mean%, bias in pp, coverage Y/n)   wrote test_calibration.csv");
