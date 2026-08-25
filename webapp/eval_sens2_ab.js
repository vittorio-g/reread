/* eval_sens2.js — Standard vs High automatic-flagging sweep on the real study, with the
 * DEPLOYED calibration (BIC gate + one-sided attentive-mode cut, R.autoFlag). Replaces the
 * data behind fig_sensitivity and the calibration-section prose numbers.
 * Design: all 70 attentive respondents + k real careless resampled, achievable rates
 * k/(70+k) from ~5% to 50%; 200 resamples per rate; MCC/precision/recall of the automatic
 * flags at both settings, plus the oracle-threshold MCC for reference. */
const R=require("./site/reread.js"); const W=R.WEIGHTS; const fs=require("fs"); const rng=R._rng(31459);
const Ms=R.parseCSV(fs.readFileSync("_study_matrix.csv","utf8")), J=Ms.header.length;
const allR=Ms.rows.map(r=>r.map(v=>{const x=Number(v);return Number.isFinite(x)?x:NaN;}));
const y=R.parseCSV(fs.readFileSync("_study_labels.csv","utf8")).rows.map(r=>Number(r[0]));
const cf=y.map((v,i)=>v===0?i:-1).filter(i=>i>=0), cl=y.map((v,i)=>v===1?i:-1).filter(i=>i>=0);
function sampleK(pool,k){const a=pool.slice();for(let i=0;i<k;i++){const j=i+Math.floor(rng()*(a.length-i));const t=a[i];a[i]=a[j];a[j]=t;}return a.slice(0,k);}
function eta(idx){const n=idx.length,mat=new Float64Array(n*J);for(let i=0;i<n;i++)for(let j=0;j<J;j++)mat[i*J+j]=allR[idx[i]][j];
  const ef=R.ensembleFeatures(mat,n,J,{seed:1,iterations:200}),O=ef.oriented,g=Math.max(0,Math.min(1,(ef.profileInfo-0.01)/0.05));
  const e=new Array(n);for(let i=0;i<n;i++)e[i]=W.b0+W.rr*O.rr[i]+g*(W.longstring*O.longstring[i]+W.person_total*O.person_total[i]);return e;}
function conf(flag,lab){let tp=0,fp=0,tn=0,fn=0;for(let i=0;i<lab.length;i++){if(flag[i]){lab[i]?tp++:fp++;}else{lab[i]?fn++:tn++;}}
  const d=Math.sqrt((tp+fp)*(tp+fn)*(tn+fp)*(tn+fn));
  return {mcc:d===0?0:(tp*tn-fp*fn)/d, prec:(tp+fp)?tp/(tp+fp):1, rec:(tp+fn)?tp/(tp+fn):0};}
function oracleMCC(e,lab){const ts=[...new Set(e)].sort((a,b)=>a-b);let best=-2;
  for(const t of ts){const {mcc}=conf(e.map(v=>v>t),lab);if(mcc>best)best=mcc;}return best;}
const mn=a=>a.reduce((s,v)=>s+v,0)/a.length, NC=cf.length, REPS=200;
const KS=[4,6,8,10,12,15,18,21,24,28,33,40,47,58,70];
const out=[["true","k","reps","std_mcc","std_prec","std_rec","high_mcc","high_prec","high_rec","oracle_mcc","gate_on","old_std_mcc","old_high_mcc","old_pi","new_pi"]];
console.log("true |  std MCC (prec/rec)   high MCC (prec/rec)   oracle  gateOn");
for(const k of KS){ if(k>cl.length) break; const tr=k/(NC+k);
  const S={mcc:[],prec:[],rec:[]},H={mcc:[],prec:[],rec:[]},O=[],OS=[],OH=[],OPI=[],NPI=[];let gon=0;
  for(let rep=0;rep<REPS;rep++){
    const idx=[...cf,...sampleK(cl,k)],lab=idx.map(i=>y[i]),e=eta(idx);
    const af=R.autoFlag(e,"standard"),ah=R.autoFlag(e,1.5);
    if(af.twoComp)gon++;
    const cs=conf(af.flagged,lab),ch=conf(ah.flagged,lab);
    S.mcc.push(cs.mcc);S.prec.push(cs.prec);S.rec.push(cs.rec);
    H.mcc.push(ch.mcc);H.prec.push(ch.prec);H.rec.push(ch.rec);
    // OLD calibration on the same replicate: cut anchored on the MIXTURE's attentive component
    const gm=R.gaussMix1D(e),pooled=Math.sqrt((gm.v1+gm.v2)/2)||1,sepq=(gm.m2-gm.m1)/pooled;
    const okOld=(gm.bic1-gm.bic2>2)&&gm.pi>0.01&&gm.pi<0.7&&sepq>1.0,s1=Math.sqrt(gm.v1);
    const oldS=okOld?e.map(v=>v>gm.m1+2.5*s1):e.map(()=>false);
    const oldH=okOld?e.map(v=>v>gm.m1+1.5*s1):e.map(()=>false);
    OS.push(conf(oldS,lab).mcc);OH.push(conf(oldH,lab).mcc);OPI.push(okOld?gm.pi:0);NPI.push(af.pi);
    O.push(oracleMCC(e,lab));}
  out.push([tr.toFixed(4),k,REPS,mn(S.mcc).toFixed(4),mn(S.prec).toFixed(4),mn(S.rec).toFixed(4),
            mn(H.mcc).toFixed(4),mn(H.prec).toFixed(4),mn(H.rec).toFixed(4),mn(O).toFixed(4),(gon/REPS).toFixed(3),
            mn(OS).toFixed(4),mn(OH).toFixed(4),mn(OPI).toFixed(4),mn(NPI).toFixed(4)]);
  console.log((100*tr).toFixed(1).padStart(5)+"% |  "+mn(S.mcc).toFixed(3)+" ("+mn(S.prec).toFixed(2)+"/"+mn(S.rec).toFixed(2)+")      "+
    mn(H.mcc).toFixed(3)+" ("+mn(H.prec).toFixed(2)+"/"+mn(H.rec).toFixed(2)+")      "+mn(O).toFixed(3)+"   "+(100*gon/REPS).toFixed(0)+
    "% | OLD std "+mn(OS).toFixed(3)+" high "+mn(OH).toFixed(3)+" | pi old "+(100*mn(OPI)).toFixed(1)+" new "+(100*mn(NPI)).toFixed(1));}
fs.writeFileSync("eval_sens2_ab.csv",out.map(r=>r.join(",")).join("\n"));
console.log("\nwrote eval_sens2.csv");
