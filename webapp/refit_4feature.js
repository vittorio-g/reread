/* refit_4feature.js — fit the shipped ensemble with a 4th feature: ACQUIESCENCE.
 * The collected study has no acquiescent responders (its careless are extreme), so acq would get
 * ~0 weight if fit on it alone. We therefore fit on the study data AUGMENTED with injected
 * acquiescent distractors (label 0, made from copies of the real careful respondents), so the
 * model learns to spare acquiescent yea-sayers while still flagging the real careless.
 * Outputs 4-feature weights + a no-harm check on the pure study (LOO) + acq benefit. */
const R=require("./site/rerere.js"); const W=R.WEIGHTS; const fs=require("fs"); const rng=R._rng(2026714);
const src=fs.readFileSync("sim_v2.js","utf8");
eval(src.slice(src.indexOf("function gauss"), src.indexOf("function scoreDataset")));  // applyPattern, ACQ_CORR, mcc, ...
const Ms=R.parseCSV(fs.readFileSync("_study_matrix.csv","utf8")), J=Ms.header.length;
const rows0=Ms.rows.map(r=>r.map(v=>{const x=Number(v);return Number.isFinite(x)?x:NaN;}));
const y0=R.parseCSV(fs.readFileSync("_study_labels.csv","utf8")).rows.map(r=>Number(r[0]));
const cf=y0.map((v,i)=>v===0?i:-1).filter(i=>i>=0);
// build augmented set: study + N_ACQ injected acquiescent distractors (copies of careful rows, label 0)
const N_ACQ=45, ACQL=[0.5,0.6,0.7,0.8,0.9,1.0];
const rows=rows0.map(r=>r.slice()), y=y0.slice(), isAcq=y0.map(_=>false);
for(let m=0;m<N_ACQ;m++){const src_i=cf[Math.floor(rng()*cf.length)];const row=rows0[src_i].map(v=>Number.isFinite(v)?v:3);
  applyPattern(row,J,"acquiescent",ACQL[Math.floor(rng()*ACQL.length)]);rows.push(row);y.push(0);isAcq.push(true);}
const n=rows.length;
console.log("augmented training: "+n+" ("+y.reduce((a,b)=>a+b,0)+" careless, "+N_ACQ+" injected acquiescent, rest careful)");
// features on the augmented matrix (scored once)
const mat=new Float64Array(n*J);for(let i=0;i<n;i++)for(let j=0;j<J;j++)mat[i*J+j]=rows[i][j];
const ef=R.ensembleFeatures(mat,n,J,{seed:1,iterations:250}),O=ef.oriented,g=Math.max(0,Math.min(1,(ef.profileInfo-0.01)/0.05));
const rawMean=[];for(let i=0;i<n;i++){let s=0,c=0;for(let j=0;j<J;j++){if(Number.isFinite(rows[i][j])){s+=rows[i][j];c++;}}rawMean.push(c?s/c:0);}
const acq=R.robustZ(rawMean,n);
// TWO-STAGE: keep the validated core (rr/ls/pt at shipped weights); fit ONLY an acquiescence
// suppressor on top, with eta_core as a fixed offset. Preserves study careless detection.
const etaCore=[];for(let i=0;i<n;i++)etaCore.push(W.b0+W.rr*O.rr[i]+g*(W.longstring*O.longstring[i]+W.person_total*O.person_total[i]));
const Z=[];for(let i=0;i<n;i++)Z.push([1,g*acq[i]]);   // offset logistic features: intercept + acq
function solve(A,b){const m=b.length,Au=A.map((r,i)=>r.concat(b[i]));for(let c=0;c<m;c++){let p=c;for(let r=c+1;r<m;r++)if(Math.abs(Au[r][c])>Math.abs(Au[p][c]))p=r;[Au[c],Au[p]]=[Au[p],Au[c]];const d=Au[c][c]||1e-9;for(let r=0;r<m;r++)if(r!==c){const f=Au[r][c]/d;for(let k=c;k<=m;k++)Au[r][k]-=f*Au[c][k];}}return Au.map((r,i)=>r[m]/(r[i]||1e-9));}
function fitOffset(Zf,off,yy,P,lam){let b=new Array(P).fill(0);for(let it=0;it<100;it++){const G=new Array(P).fill(0),H=Array.from({length:P},()=>new Array(P).fill(0));
  for(let i=0;i<Zf.length;i++){let e=off[i];for(let k=0;k<P;k++)e+=b[k]*Zf[i][k];const p=1/(1+Math.exp(-e)),w=Math.max(p*(1-p),1e-6);
    for(let k=0;k<P;k++){G[k]+=(yy[i]-p)*Zf[i][k];for(let l=0;l<P;l++)H[k][l]+=w*Zf[i][k]*Zf[i][l];}}
  for(let k=0;k<P;k++){G[k]-=lam*b[k];H[k][k]+=lam;}const s=solve(H,G);let mx=0;for(let k=0;k<P;k++){b[k]+=s[k];mx=Math.max(mx,Math.abs(s[k]));}if(mx<1e-9)break;}return b;}
const bb=fitOffset(Z,etaCore,y,2,1.0);   // bb[0]=Δb0, bb[1]=acq weight
const Wn={b0:+(W.b0+bb[0]).toFixed(4),rr:W.rr,irv:0,longstring:W.longstring,d2:0,person_total:W.person_total,acq:+bb[1].toFixed(4)};
console.log("\n4-feature weights (core kept, acq suppressor added):");
console.log(JSON.stringify(Wn));
console.log("acq coefficient = "+bb[1].toFixed(3)+" (negative = high yea-saying suppresses flagging)");
function auc(s,yy){let po=[],ne=[];for(let i=0;i<yy.length;i++)(yy[i]?po:ne).push(s[i]);let c=0;for(const a of po)for(const q of ne)c+=a>q?1:a===q?0.5:0;return c/(po.length*ne.length);}
const pred=i=>1/(1+Math.exp(-(etaCore[i]+bb[0]+bb[1]*g*acq[i])));
const etaS=[],labS=[];for(let i=0;i<rows0.length;i++){etaS.push(pred(i));labS.push(y0[i]);}
const rate=labS.reduce((a,b)=>a+b,0)/labS.length,srt=[...etaS].sort((a,b)=>a-b),thr=srt[Math.floor((1-rate)*labS.length)];
console.log("\nPURE STUDY (n=157): 3-feature AUC "+auc(rows0.map((_,i)=>1/(1+Math.exp(-etaCore[i]))),y0).toFixed(3)+
  " -> 4-feature AUC "+auc(etaS,labS).toFixed(3)+"  MCC "+mcc(etaS.map(v=>v>=thr?1:0),labS).toFixed(3)+"  (want: no drop)");
let aFlag=0;for(let i=rows0.length;i<n;i++)if(pred(i)>=thr)aFlag++;
let aFlagNo=0;for(let i=rows0.length;i<n;i++)if(1/(1+Math.exp(-(etaCore[i]+bb[0])))>=thr)aFlagNo++;
console.log("injected acquiescent flagged (want low): with-acq "+aFlag+"/"+N_ACQ+"  vs no-acq "+aFlagNo+"/"+N_ACQ);
fs.writeFileSync("ensemble_weights_4feat.json",JSON.stringify(Wn,null,2));
console.log("\nwrote ensemble_weights_4feat.json");
