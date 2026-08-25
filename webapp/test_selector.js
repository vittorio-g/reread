/* test_selector.js — use calibrated pi-hat to AUTO-SELECT Standard vs High (not to place the cut).
 * Per resample: if piHat < theta -> Standard cut, else High cut. Compare to always-Std, always-High,
 * and the per-resample envelope max(Std,High) (=perfect selector upper bound). Real study sweep. */
const R=require("./site/rerere.js"); const W=R.WEIGHTS; const fs=require("fs"); const rng=R._rng(31459);
const Ms=R.parseCSV(fs.readFileSync("_study_matrix.csv","utf8")), J=Ms.header.length;
const allR=Ms.rows.map(r=>r.map(v=>{const x=Number(v);return Number.isFinite(x)?x:NaN;}));
const y=R.parseCSV(fs.readFileSync("_study_labels.csv","utf8")).rows.map(r=>Number(r[0]));
const cf=y.map((v,i)=>v===0?i:-1).filter(i=>i>=0), cl=y.map((v,i)=>v===1?i:-1).filter(i=>i>=0);
function sampleK(pool,k){const a=pool.slice();for(let i=0;i<k;i++){const j=i+Math.floor(rng()*(a.length-i));const t=a[i];a[i]=a[j];a[j]=t;}return a.slice(0,k);}
function gateOK(gm){const pooled=Math.sqrt((gm.v1+gm.v2)/2)||1,sep=(gm.m2-gm.m1)/pooled;return (gm.bic1-gm.bic2>2)&&gm.pi>0.01&&gm.pi<0.7&&sep>1.0;}
function eta(idx){const n=idx.length,mat=new Float64Array(n*J);for(let i=0;i<n;i++)for(let j=0;j<J;j++)mat[i*J+j]=allR[idx[i]][j];
  const ef=R.ensembleFeatures(mat,n,J,{seed:1,iterations:200}),O=ef.oriented,g=Math.max(0,Math.min(1,(ef.profileInfo-0.01)/0.05));
  const e=new Array(n);for(let i=0;i<n;i++)e[i]=W.b0+W.rr*O.rr[i]+g*(W.longstring*O.longstring[i]+W.person_total*O.person_total[i]);return e;}
function mcc(flag,lab){let tp=0,fp=0,tn=0,fn=0;for(let i=0;i<lab.length;i++){if(flag[i]){lab[i]?tp++:fp++;}else{lab[i]?fn++:tn++;}}
  const d=Math.sqrt((tp+fp)*(tp+fn)*(tn+fp)*(tn+fn));return d===0?0:(tp*tn-fp*fn)/d;}
const THETAS=[0.13,0.15,0.17,0.19,0.22], mn=a=>a.reduce((s,v)=>s+v,0)/a.length, NC=cf.length, REPS=200;
console.log("AUTO-SELECTOR: piHat<theta -> Standard, else High. Real study sweep, "+REPS+" resamples/rate.\n");
console.log("true | always-Std  always-High | "+THETAS.map(t=>"sel@"+t).join("  ")+" | envelope(max)");
const out=[["true","std","high",...THETAS.map(t=>"sel_"+t),"envelope"]];
for(const p of [0.05,0.10,0.15,0.20,0.25,0.30]){
  const k=Math.round(NC*p/(1-p)); if(k>cl.length) continue; const tr=k/(NC+k);
  const S=[],H=[],ENV=[],SEL={}; for(const t of THETAS) SEL[t]=[];
  for(let rep=0;rep<REPS;rep++){
    const idx=[...cf,...sampleK(cl,k)],lab=idx.map(i=>y[i]),e=eta(idx),gm=R.gaussMix1D(e),ok=gateOK(gm),s1=Math.sqrt(gm.v1);
    const stdF=ok?e.map(v=>v>gm.m1+2.5*s1):e.map(()=>false), highF=ok?e.map(v=>v>gm.m1+1.5*s1):e.map(()=>false), pi=ok?gm.pi:0;
    const ms=mcc(stdF,lab),mh=mcc(highF,lab); S.push(ms);H.push(mh);ENV.push(Math.max(ms,mh));
    for(const t of THETAS) SEL[t].push(pi<t?ms:mh);
  }
  const row=[tr.toFixed(3),mn(S).toFixed(3),mn(H).toFixed(3),...THETAS.map(t=>mn(SEL[t]).toFixed(3)),mn(ENV).toFixed(3)];
  out.push(row);
  console.log((100*tr).toFixed(0).padStart(4)+"% |   "+mn(S).toFixed(3)+"      "+mn(H).toFixed(3)+"    | "+
    THETAS.map(t=>mn(SEL[t]).toFixed(3)).join("  ")+" |   "+mn(ENV).toFixed(3));
}
fs.writeFileSync("test_selector.csv",out.map(r=>r.join(",")).join("\n"));
// pooled means over the realistic 5-30% band (equal weight per rate)
const rows=out.slice(1), cols=out[0];
console.log("\nPOOLED (5-30%, equal weight):");
for(let c=1;c<cols.length;c++){const m=mn(rows.map(r=>Number(r[c])));console.log("  "+cols[c].padEnd(10)+m.toFixed(3));}
console.log("\nwrote test_selector.csv");
