/* gen_calib_decomp.js — WHY does pi-hat collapse at high contamination, and who is to blame?
 * Same design as gen_fig_calib.js (70 attentive + k real careless from Study 1), but at each
 * true rate we fit the SAME two-component mixture separately on:
 *    eta (shipped ensemble), O.rr (rc), O.longstring, O.person_total
 * and we also record, per feature, the AUC against the true label. That separates
 *   "the index lost the signal"  (AUC falls)   from
 *   "the estimator broke"        (AUC holds but pi-hat collapses).
 * Plus two candidate absolute anchors that do NOT depend on sample composition:
 *    frac_rc_low = share of respondents whose RAW rc z-score <= 1.5 (and <= 0.0)
 * Output: calib_decomp.csv
 */
const R=require("./site/rerere.js"); const W=R.WEIGHTS; const fs=require("fs"); const rng=R._rng(31459);
const Ms=R.parseCSV(fs.readFileSync("_study_matrix.csv","utf8")), J=Ms.header.length;
const allR=Ms.rows.map(r=>r.map(v=>{const x=Number(v);return Number.isFinite(x)?x:NaN;}));
const y=R.parseCSV(fs.readFileSync("_study_labels.csv","utf8")).rows.map(r=>Number(r[0]));
const cf=y.map((v,i)=>v===0?i:-1).filter(i=>i>=0), cl=y.map((v,i)=>v===1?i:-1).filter(i=>i>=0);
function sampleK(pool,k){const a=pool.slice();for(let i=0;i<k;i++){const j=i+Math.floor(rng()*(a.length-i));const t=a[i];a[i]=a[j];a[j]=t;}return a.slice(0,k);}
function gateOK(gm){const pooled=Math.sqrt((gm.v1+gm.v2)/2)||1,sep=(gm.m2-gm.m1)/pooled;return (gm.bic1-gm.bic2>2)&&gm.pi>0.01&&gm.pi<0.7&&sep>1.0;}
/* AUC, higher score = more careless (lab 1 = careless) */
function auc(score,lab){const idx=score.map((s,i)=>[s,lab[i]]).sort((a,b)=>a[0]-b[0]);
  let r=0,n1=0,n0=0,i=0;
  while(i<idx.length){let j=i;while(j+1<idx.length&&idx[j+1][0]===idx[i][0])j++;
    const rank=(i+j)/2+1;for(let t=i;t<=j;t++){if(idx[t][1]===1){r+=rank;n1++;}else n0++;}i=j+1;}
  return n1&&n0?(r-n1*(n1+1)/2)/(n1*n0):NaN;}
const pct=(a,q)=>{const s=[...a].sort((x,y)=>x-y);return s[Math.min(s.length-1,Math.max(0,Math.round(q*(s.length-1))))];};
const mn=a=>a.reduce((s,v)=>s+v,0)/a.length;

const REPS=150, KS=[4,6,8,10,12,15,18,21,24,28,33,40,47,58];
const FEATS=["eta","rc","longstring","person_total"];
const hdr=["true"];
for(const f of FEATS) hdr.push("pi_"+f,"pi_"+f+"_p5","pi_"+f+"_p95","gate_"+f,"auc_"+f);
hdr.push("frac_rcraw_le1_5","frac_rcraw_le0","sep_eta");
const out=[hdr];
console.log("true | "+FEATS.map(f=>("pi_"+f).padEnd(9)+"auc").join(" | ")+" | absAnchor");
for(const k of KS){ if(k>cl.length) break; const tr=k/(cf.length+k);
  const acc={}; FEATS.forEach(f=>acc[f]={pi:[],gate:0,auc:[]});
  const anch15=[],anch0=[],seps=[];
  for(let rep=0;rep<REPS;rep++){
    const idx=[...cf,...sampleK(cl,k)], n=idx.length, lab=idx.map(i=>y[i]);
    const mat=new Float64Array(n*J);
    for(let i=0;i<n;i++)for(let j=0;j<J;j++)mat[i*J+j]=allR[idx[i]][j];
    const ef=R.ensembleFeatures(mat,n,J,{seed:1,iterations:200}),O=ef.oriented;
    const g=Math.max(0,Math.min(1,(ef.profileInfo-0.01)/0.05));
    const series={ eta:new Array(n), rc:Array.from(O.rr), longstring:Array.from(O.longstring), person_total:Array.from(O.person_total) };
    for(let i=0;i<n;i++) series.eta[i]=W.b0+W.rr*O.rr[i]+g*(W.longstring*O.longstring[i]+W.person_total*O.person_total[i]);
    for(const f of FEATS){ const gm=R.gaussMix1D(series[f]), ok=gateOK(gm);
      if(ok)acc[f].gate++; acc[f].pi.push(ok?gm.pi:0); acc[f].auc.push(auc(series[f],lab));
      if(f==="eta"){const pooled=Math.sqrt((gm.v1+gm.v2)/2)||1;seps.push((gm.m2-gm.m1)/pooled);} }
    /* absolute anchors on the RAW rc z (z=0 means "no more coherent than own chance") */
    const zraw=ef.base.z;
    anch15.push(zraw.filter(v=>Number.isFinite(v)&&v<=1.5).length/n);
    anch0.push(zraw.filter(v=>Number.isFinite(v)&&v<=0).length/n);
  }
  const row=[tr.toFixed(4)];
  for(const f of FEATS) row.push(mn(acc[f].pi).toFixed(4),pct(acc[f].pi,.05).toFixed(4),pct(acc[f].pi,.95).toFixed(4),
                                (acc[f].gate/REPS).toFixed(3),mn(acc[f].auc).toFixed(4));
  row.push(mn(anch15).toFixed(4),mn(anch0).toFixed(4),mn(seps).toFixed(3));
  out.push(row);
  console.log((100*tr).toFixed(1).padStart(5)+"% | "+FEATS.map(f=>(100*mn(acc[f].pi)).toFixed(1).padStart(5)+"%   "+mn(acc[f].auc).toFixed(3)).join(" | ")
    +" | rc<=1.5: "+(100*mn(anch15)).toFixed(1)+"%");
}
fs.writeFileSync("calib_decomp.csv",out.map(r=>r.join(",")).join("\n"));
console.log("\nwrote calib_decomp.csv");
