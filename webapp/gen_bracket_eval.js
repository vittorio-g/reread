/* gen_bracket_eval.js — SOLUTION B evaluated honestly.
 * Proposal: replace the single ensemble estimate pi_eta with the BRACKET [pi_rc, pi_PT].
 * Evaluated PER REPLICATE (not on across-replicate means, which would fake the coverage):
 *   cover      = does [min,max] contain the true rate in THIS sample?
 *   width      = how wide is the interval (a vacuous interval is not useful)
 *   mid        = bracket midpoint as a point estimate
 * A pi of 0 means the mixture gate REJECTED (no estimate). A bracket whose lower end is a
 * rejected gate trivially covers small truths, so we report coverage twice:
 *   cover_all   — over all replicates (0 treated as an estimate)
 *   cover_valid — restricted to replicates where BOTH gates fired (an honest bracket)
 * Baseline for comparison: the shipped pi_eta.
 * Output: bracket_eval.csv
 */
const R=require("./site/reread.js"); const W=R.WEIGHTS; const fs=require("fs"); const rng=R._rng(90210);
const Ms=R.parseCSV(fs.readFileSync("_study_matrix.csv","utf8")), J=Ms.header.length;
const allR=Ms.rows.map(r=>r.map(v=>{const x=Number(v);return Number.isFinite(x)?x:NaN;}));
const y=R.parseCSV(fs.readFileSync("_study_labels.csv","utf8")).rows.map(r=>Number(r[0]));
const cf=y.map((v,i)=>v===0?i:-1).filter(i=>i>=0), cl=y.map((v,i)=>v===1?i:-1).filter(i=>i>=0);
function sampleK(pool,k){const a=pool.slice();for(let i=0;i<k;i++){const j=i+Math.floor(rng()*(a.length-i));const t=a[i];a[i]=a[j];a[j]=t;}return a.slice(0,k);}
function gateOK(gm){const p=Math.sqrt((gm.v1+gm.v2)/2)||1,sep=(gm.m2-gm.m1)/p;return (gm.bic1-gm.bic2>2)&&gm.pi>0.01&&gm.pi<0.7&&sep>1.0;}
function piOf(arr){const gm=R.gaussMix1D(arr);const ok=gateOK(gm);return {pi:ok?gm.pi:0, ok};}
const mn=a=>a.length?a.reduce((s,v)=>s+v,0)/a.length:NaN;

const REPS=150, KS=[4,6,8,10,12,15,18,21,24,28,33,40,47,58];
const out=[["true","cover_all","cover_valid","n_valid","width_mean","abs_err_mid","abs_err_eta","pi_rc_mean","pi_pt_mean","pi_eta_mean"]];
console.log("true  | coverAll coverValid (nValid) | width | |err|mid  |err|eta");
for(const k of KS){ if(k>cl.length) break; const tr=k/(cf.length+k);
  const cAll=[],cVal=[],wid=[],eMid=[],eEta=[],pr=[],pp=[],pe=[];
  for(let rep=0;rep<REPS;rep++){
    const idx=[...cf,...sampleK(cl,k)], n=idx.length;
    const mat=new Float64Array(n*J);
    for(let i=0;i<n;i++)for(let j=0;j<J;j++)mat[i*J+j]=allR[idx[i]][j];
    const ef=R.ensembleFeatures(mat,n,J,{seed:1,iterations:200}),O=ef.oriented;
    const g=Math.max(0,Math.min(1,(ef.profileInfo-0.01)/0.05));
    const eta=new Array(n);
    for(let i=0;i<n;i++)eta[i]=W.b0+W.rr*O.rr[i]+g*(W.longstring*O.longstring[i]+W.person_total*O.person_total[i]);
    const E=piOf(eta), Rc=piOf(Array.from(O.rr)), Pt=piOf(Array.from(O.person_total));
    const lo=Math.min(Rc.pi,Pt.pi), hi=Math.max(Rc.pi,Pt.pi);
    const cov=(tr>=lo&&tr<=hi);
    cAll.push(cov?1:0); if(Rc.ok&&Pt.ok){cVal.push(cov?1:0);}
    wid.push(hi-lo); eMid.push(Math.abs((lo+hi)/2-tr)); eEta.push(Math.abs(E.pi-tr));
    pr.push(Rc.pi); pp.push(Pt.pi); pe.push(E.pi);
  }
  out.push([tr.toFixed(4),mn(cAll).toFixed(3),(cVal.length?mn(cVal):NaN).toFixed(3),cVal.length,
            mn(wid).toFixed(4),mn(eMid).toFixed(4),mn(eEta).toFixed(4),
            mn(pr).toFixed(4),mn(pp).toFixed(4),mn(pe).toFixed(4)]);
  console.log((100*tr).toFixed(1).padStart(5)+"% |  "+(100*mn(cAll)).toFixed(0).padStart(3)+"%     "+
    (cVal.length?(100*mn(cVal)).toFixed(0):" na").padStart(3)+"%  ("+String(cVal.length).padStart(3)+")  | "+
    (100*mn(wid)).toFixed(1).padStart(5)+" | "+(100*mn(eMid)).toFixed(1).padStart(5)+"    "+(100*mn(eEta)).toFixed(1).padStart(5));
}
fs.writeFileSync("bracket_eval.csv",out.map(r=>r.join(",")).join("\n"));
console.log("\nwrote bracket_eval.csv");
