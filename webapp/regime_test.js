/* regime_test.js — validate the label-free REGIME TEST on Study 1 and the 16 Study-2 datasets.
 *
 * Test: fit the same 2-component mixture separately on the ensemble log-odds (eta) and on the
 * oriented rc feature. If pi_rc > pi_eta -> "collapsed" regime (true careless rate is high and
 * pi_eta is not identified). If pi_rc < pi_eta -> "identified" regime.
 * We also report pi_PT and the bracket [pi_rc, pi_PT], and compare against each dataset's
 * ground-truth positive rate (a PROXY for the true careless rate, not the truth itself).
 */
const R=require("./site/rerere.js"); const W=R.WEIGHTS; const fs=require("fs"); const rng=R._rng(20260717);
const N_CAP=4000, ITER=200;
const D="../Dataset/";
const SETS=[
 // --- Study 2, Table 3 (core, 10) ---
 ["Kay S2","core",D+"gt_benchmark_candidates/kay_idris_idria/kay_matrix_s2.csv",D+"gt_benchmark_candidates/kay_idris_idria/kay_labels_s2.csv"],
 ["warning IPIP-NEO","core",D+"gt_benchmark_candidates/warning_ipipneo300/_matrix.csv",D+"gt_benchmark_candidates/warning_ipipneo300/_labels.csv"],
 ["Kay S1","core",D+"gt_benchmark_candidates/kay_idris_idria/kay_matrix_s1.csv",D+"gt_benchmark_candidates/kay_idris_idria/kay_labels_s1.csv"],
 ["opsy 16PF","core",D+"gt_benchmark_candidates/opsy_16pf_matrix.csv",D+"gt_benchmark_candidates/opsy_16pf_labels.csv"],
 ["smarvus","core",D+"gt_benchmark_candidates/smarvus/smarvus_matrix.csv",D+"gt_benchmark_candidates/smarvus/smarvus_labels.csv"],
 ["Kay S6","core",D+"gt_benchmark_candidates/kay_idris_idria/kay_matrix_s6.csv",D+"gt_benchmark_candidates/kay_idris_idria/kay_labels_s6.csv"],
 ["Kay S5","core",D+"gt_benchmark_candidates/kay_idris_idria/kay_matrix_s5.csv",D+"gt_benchmark_candidates/kay_idris_idria/kay_labels_s5.csv"],
 ["Duckworth VCL","core",D+"gt_benchmark_candidates/duckworth_grit_vcl/duckworth_matrix.csv",D+"gt_benchmark_candidates/duckworth_grit_vcl/duckworth_labels.csv"],
 ["Douglas 2023","core",D+"gt_benchmark_candidates/douglas2023_dataquality/douglas_matrix.csv",D+"gt_benchmark_candidates/douglas2023_dataquality/douglas_labels.csv"],
 ["Krause youth","core",D+"gt_benchmark_candidates/krause_ier_youth/krause_matrix.csv",D+"gt_benchmark_candidates/krause_ier_youth/krause_labels.csv"],
 // --- Study 2, Table 4 (boundary, 6) ---
 ["Mastroianni 2022","bound",D+"learning_to_pay_attention/mastroianni2022_attitude/_matrix.csv",D+"learning_to_pay_attention/mastroianni2022_attitude/_labels.csv"],
 ["Ivanov 2021","bound",D+"learning_to_pay_attention/ivanov2021_racialresentment/_matrix.csv",D+"learning_to_pay_attention/ivanov2021_racialresentment/_labels.csv"],
 ["O'Grady 2019","bound",D+"learning_to_pay_attention/ogrady2019_moralfoundations/_matrix.csv",D+"learning_to_pay_attention/ogrady2019_moralfoundations/_labels.csv"],
 ["Pennycook 2020","bound",D+"learning_to_pay_attention/pennycook2020_covid/_matrix.csv",D+"learning_to_pay_attention/pennycook2020_covid/_labels.csv"],
 ["Moss 2023","bound",D+"learning_to_pay_attention/moss2023_ethical/_matrix.csv",D+"learning_to_pay_attention/moss2023_ethical/_labels.csv"],
 ["Alvarez 2019","bound",D+"learning_to_pay_attention/alvarez2019_inattentive/_matrix.csv",D+"learning_to_pay_attention/alvarez2019_inattentive/_labels.csv"],
];
function gateOK(gm){const p=Math.sqrt((gm.v1+gm.v2)/2)||1,sep=(gm.m2-gm.m1)/p;return (gm.bic1-gm.bic2>2)&&gm.pi>0.01&&gm.pi<0.7&&sep>1.0;}
function piOf(arr){const gm=R.gaussMix1D(arr);return gateOK(gm)?gm.pi:0;}
function scoreOne(name,tag,rows,lab,J){
  const n=rows.length,mat=new Float64Array(n*J);
  for(let i=0;i<n;i++)for(let j=0;j<J;j++)mat[i*J+j]=rows[i][j];
  const ef=R.ensembleFeatures(mat,n,J,{seed:1,iterations:ITER}),O=ef.oriented;
  const g=Math.max(0,Math.min(1,(ef.profileInfo-0.01)/0.05));
  const eta=new Array(n);
  for(let i=0;i<n;i++)eta[i]=W.b0+W.rr*O.rr[i]+g*(W.longstring*O.longstring[i]+W.person_total*O.person_total[i]);
  const piE=piOf(eta), piR=piOf(Array.from(O.rr)), piP=piOf(Array.from(O.person_total));
  const gt=lab.reduce((s,v)=>s+v,0)/n;
  const verdict=(piR>piE)?"COLLAPSED":"identified";
  const inBracket=(gt>=Math.min(piR,piP)&&gt<=Math.max(piR,piP));
  console.log(name.padEnd(18)+tag.padEnd(7)+"n="+String(n).padStart(5)+" J="+String(J).padStart(3)+
    " | GT="+(100*gt).toFixed(1).padStart(5)+"%  piEta="+(100*piE).toFixed(1).padStart(5)+
    "%  piRc="+(100*piR).toFixed(1).padStart(5)+"%  piPT="+(100*piP).toFixed(1).padStart(5)+
    "% | diff="+(100*(piR-piE)).toFixed(1).padStart(6)+"  -> "+verdict.padEnd(10)+" bracket:"+(inBracket?"yes":"no"));
  return [name,tag,n,J,gt.toFixed(4),piE.toFixed(4),piR.toFixed(4),piP.toFixed(4),(piR-piE).toFixed(4),verdict,inBracket?1:0];
}
const out=[["dataset","group","n","J","gt_rate","pi_eta","pi_rc","pi_pt","pi_rc_minus_eta","verdict","gt_in_bracket"]];
console.log("=== STUDY 1 (full sample; true careless rate is KNOWN = 87/157 = 55.4%) ===");
{ const Ms=R.parseCSV(fs.readFileSync("_study_matrix.csv","utf8")),J=Ms.header.length;
  const rows=Ms.rows.map(r=>r.map(v=>{const x=Number(v);return Number.isFinite(x)?x:NaN;}));
  const lab=R.parseCSV(fs.readFileSync("_study_labels.csv","utf8")).rows.map(r=>Number(r[0]));
  out.push(scoreOne("Study 1 (full)","truth",rows,lab,J)); }
console.log("\n=== STUDY 2 (16 datasets; GT rate is a PROXY criterion, not the true careless rate) ===");
for(const [name,tag,mf,lf] of SETS){
  let P,lab;
  try{ P=R.parseCSV(fs.readFileSync(mf,"utf8")); lab=R.parseCSV(fs.readFileSync(lf,"utf8")).rows.map(r=>Number(r[0])); }
  catch(e){ console.log(name.padEnd(18)+"SKIP (unreadable: "+e.message.slice(0,40)+")"); continue; }
  const J=P.header.length;
  let rows=P.rows.map(r=>r.map(v=>{const x=Number(v);return Number.isFinite(x)?x:NaN;}));
  const m=Math.min(rows.length,lab.length); rows=rows.slice(0,m); lab=lab.slice(0,m);
  if(rows.length>N_CAP){const o=rows.map((_,i)=>i);for(let i=0;i<N_CAP;i++){const j=i+Math.floor(rng()*(o.length-i));const t=o[i];o[i]=o[j];o[j]=t;}
    const s=o.slice(0,N_CAP);rows=s.map(i=>rows[i]);lab=s.map(i=>lab[i]);}
  out.push(scoreOne(name,tag,rows,lab,J));
}
fs.writeFileSync("regime_test.csv",out.map(r=>r.join(",")).join("\n"));
console.log("\nwrote regime_test.csv");
