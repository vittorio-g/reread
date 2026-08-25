/* diag_flatzero.js — why are warning-control and TISP-Australia FLAT AT ZERO in the injection
 * sweep even at 25-40% true careless? Bug or mechanism? Open the box: for one resample at
 * several rates, print the eta distribution, the one-sided fit (m, sigma, cut), the tail count,
 * every gate component (dBIC, w, sep, tailP), and the raw piHat BEFORE clamping.
 */
const R = require("./site/reread.js"); const W = R.WEIGHTS; const fs = require("fs");
let SEED = 555; const rng = () => { SEED = (SEED*1103515245+12345)&0x7fffffff; return SEED/0x7fffffff; };
const Z = 2.5, PHI1 = 0.8413447;
const D = "../Dataset/gt_benchmark_candidates";
const SETS = [
  ["warning control", D+"/warning_control", "_labels2.csv", "_matrix.csv"],
  ["TISP Australia",  D+"/tisp",            "_labels_AUS.csv", "_matrix_AUS.csv"],
];
function shuffle(a){for(let i=a.length-1;i>0;i--){const j=Math.floor(rng()*(i+1));const t=a[i];a[i]=a[j];a[j]=t;}return a;}
const q=(a,p)=>{const s=[...a].sort((x,y)=>x-y);return s[Math.min(s.length-1,Math.max(0,Math.round(p*(s.length-1))))];};
for (const [nm, dir, lf, mf] of SETS) {
  const P = R.parseCSV(fs.readFileSync(dir+"/"+mf,"utf8")), J = P.header.length;
  const rows = P.rows.map(r=>r.map(v=>{const x=Number(v);return Number.isFinite(x)?x:NaN;}));
  const y = R.parseCSV(fs.readFileSync(dir+"/"+lf,"utf8")).rows.map(r=>Number(r[0]));
  const n0 = Math.min(rows.length,y.length);
  const att=[],car=[];
  for(let i=0;i<n0;i++)(y[i]?car:att).push(rows[i]);
  console.log("\n================ "+nm+"  (J="+J+", att="+att.length+", car="+car.length+") ================");
  for (const p of [0.05,0.25,0.40]) {
    let k=Math.floor(Math.min(car.length, att.length*p/(1-p)));
    let m=Math.round(k*(1-p)/p);
    if(m>att.length){m=att.length;k=Math.round(m*p/(1-p));}
    const cap=1400; if(k+m>cap){const s=cap/(k+m);k=Math.round(k*s);m=Math.round(m*s);}
    if(k<3||k+m<150){console.log("p="+(100*p)+"% not constructible");continue;}
    const A=shuffle(att.slice()).slice(0,m), C=shuffle(car.slice()).slice(0,k);
    const sub=A.concat(C), n=sub.length, lab=[...new Array(m).fill(0),...new Array(k).fill(1)];
    const mat=new Float64Array(n*J);
    for(let i=0;i<n;i++)for(let j=0;j<J;j++)mat[i*J+j]=sub[i][j];
    const ef=R.ensembleFeatures(mat,n,J,{seed:7,iterations:150}),O=ef.oriented;
    const g=Math.max(0,Math.min(1,(ef.profileInfo-0.01)/0.05));
    const eta=new Array(n);
    for(let i=0;i<n;i++)eta[i]=W.b0+W.rr*O.rr[i]+g*(W.longstring*O.longstring[i]+W.person_total*O.person_total[i]);
    const af=R.autoFlag(eta,"standard");
    const gm=R.gaussMix1D(eta);
    const pooled=Math.sqrt((gm.v1+gm.v2)/2)||1;
    const cut=af.attentiveMode+Z*af.attentiveSigma;
    const kTail=eta.filter(v=>v>cut).length;
    const c1=af.attentiveMode+af.attentiveSigma;
    const nBelow=eta.filter(v=>v<c1).length;
    const piRaw=1-(nBelow/PHI1)/n;
    /* separation of the TRUE groups in eta, to see whether the signal exists at all */
    const eAtt=eta.slice(0,m), eCar=eta.slice(m);
    const mu=a=>a.reduce((s,v)=>s+v,0)/a.length;
    const sd=a=>{const u=mu(a);return Math.sqrt(mu(a.map(v=>(v-u)*(v-u))));};
    const dTrue=(mu(eCar)-mu(eAtt))/Math.sqrt((sd(eAtt)**2+sd(eCar)**2)/2);
    console.log("--- true "+(100*k/n).toFixed(1)+"%  n="+n+" (att "+m+" + car "+k+")  gate g="+g.toFixed(2));
    console.log("    eta: p5="+q(eta,.05).toFixed(2)+" p50="+q(eta,.5).toFixed(2)+" p95="+q(eta,.95).toFixed(2)+
                "   TRUE group separation d="+dTrue.toFixed(2));
    console.log("    leftFit: m="+af.attentiveMode.toFixed(2)+" sigma="+af.attentiveSigma.toFixed(2)+
                " -> cut="+cut.toFixed(2)+"   kTail="+kTail+"/"+n+" ("+(100*kTail/n).toFixed(1)+"%)");
    console.log("    BIC gate: dBIC="+(gm.bic1-gm.bic2).toFixed(1)+"  w="+gm.pi.toFixed(3)+
                "  sep="+((gm.m2-gm.m1)/pooled).toFixed(2)+"  -> bicGate="+af.bicGate);
    console.log("    tail branch: tailP="+(af.tailP!==undefined?af.tailP.toExponential(2):"n/a")+
                "  -> twoComp="+af.twoComp);
    console.log("    piHat raw="+(100*piRaw).toFixed(1)+"%  (nBelow(m+sigma)="+nBelow+"/"+n+")   reported pi="+
                (100*af.pi).toFixed(1)+"%   flagged="+(100*af.estRate).toFixed(1)+"%");
  }
}
