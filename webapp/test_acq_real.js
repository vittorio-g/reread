/* test_acq_real.js — does the acquiescence feature work on REAL data, or only in simulation?
 * Part A: use REAL clean respondents (Johnson IPIP, 300 people x 100 items) as the base, inject the
 *          same careless + response-style distractors, refit with/without acq (per-rep 5-fold CV).
 *          If acq still halves distractor-FPR on REAL profiles, it is not a synthetic artifact.
 * Part B: the collected study (n=157), LOO with/without acq — check acq does NOT harm real careless
 *          detection (the study has no acquiescent distractors, so this is a no-regression check). */
global.self=global; require("./site/demo_data.js"); const csv0=global.RERERE_DEMO_BASE;
const R=require("./site/rerere.js"); const W=R.WEIGHTS; const fs=require("fs"); const rng=R._rng(4242);
const src=fs.readFileSync("sim_v2.js","utf8");
eval(src.slice(src.indexOf("function gauss"), src.indexOf("function scoreDataset")));  // inject2, uWeight, wmcc, mcc, ...
function featRows(out,J){
  const n=out.length,mat=new Float64Array(n*J);for(let i=0;i<n;i++)for(let j=0;j<J;j++)mat[i*J+j]=out[i][j];
  const ef=R.ensembleFeatures(mat,n,J,{seed:1,iterations:150}),O=ef.oriented,g=Math.max(0,Math.min(1,(ef.profileInfo-0.01)/0.05));
  const rawMean=[];for(let i=0;i<n;i++){let s=0,c=0;for(let j=0;j<J;j++){if(Number.isFinite(out[i][j])){s+=out[i][j];c++;}}rawMean.push(c?s/c:0);}
  const acq=R.robustZ(rawMean,n),X=[];
  for(let i=0;i<n;i++)X.push([1,O.rr[i],g*O.longstring[i],g*O.person_total[i],acq[i]]);
  return X;
}
function auc(s,y){let po=[],ne=[];for(let i=0;i<y.length;i++)(y[i]?po:ne).push(s[i]);let c=0;for(const a of po)for(const b of ne)c+=a>b?1:a===b?0.5:0;return po.length&&ne.length?c/(po.length*ne.length):0;}
function solve(A,b){const m=b.length,Au=A.map((r,i)=>r.concat(b[i]));for(let c=0;c<m;c++){let p=c;for(let r=c+1;r<m;r++)if(Math.abs(Au[r][c])>Math.abs(Au[p][c]))p=r;[Au[c],Au[p]]=[Au[p],Au[c]];const d=Au[c][c]||1e-9;for(let r=0;r<m;r++)if(r!==c){const f=Au[r][c]/d;for(let k=c;k<=m;k++)Au[r][k]-=f*Au[c][k];}}return Au.map((r,i)=>r[m]/(r[i]||1e-9));}
function fit(rows,y,P){let b=new Array(P).fill(0);for(let it=0;it<80;it++){const g=new Array(P).fill(0),H=Array.from({length:P},()=>new Array(P).fill(0));
  for(let i=0;i<rows.length;i++){let e=0;for(let k=0;k<P;k++)e+=b[k]*rows[i][k];const p=1/(1+Math.exp(-e)),w=Math.max(p*(1-p),1e-6);
    for(let k=0;k<P;k++){g[k]+=(y[i]-p)*rows[i][k];for(let l=0;l<P;l++)H[k][l]+=w*rows[i][k]*rows[i][l];}}
  for(let k=1;k<P;k++){g[k]-=0.5*b[k];H[k][k]+=0.5;}const s=solve(H,g);let mx=0;for(let k=0;k<P;k++){b[k]+=s[k];mx=Math.max(mx,Math.abs(s[k]));}if(mx<1e-8)break;}return b;}
const pred=(b,row)=>{let e=0;for(let k=0;k<b.length;k++)e+=b[k]*row[k];return 1/(1+Math.exp(-e));};
function oracleWMCC(p,y,w){const s=[...p].sort((a,b)=>a-b);let best=0,bt=0;const step=Math.max(1,Math.floor(s.length/100));for(let t=0;t<s.length;t+=step){const thr=s[t];const f=p.map(v=>v>=thr);const m=wmcc(f,y,w);if(m>best){best=m;bt=thr;}}return {mcc:best,thr:bt};}
function cvOOF(X,y,cols){const n=X.length,P=cols.length+1,oof=new Array(n);
  for(let f=0;f<5;f++){const tr=[],te=[];for(let i=0;i<n;i++)(i%5===f?te:tr).push(i);
    const b=fit(tr.map(i=>[1,...cols.map(c=>X[i][c])]),tr.map(i=>y[i]),P);
    for(const i of te)oof[i]=pred(b,[1,...cols.map(c=>X[i][c])]);}return oof;}
const mn=a=>a.reduce((s,v)=>s+v,0)/a.length;

// ---------- Part A: REAL Johnson base + injection ----------
const pj=R.parseCSV(csv0), baseAll=pj.rows.map(r=>r.slice(1).map(Number)), Jj=pj.header.length-1;
console.log("Part A: REAL Johnson base ("+baseAll.length+" real respondents x "+Jj+" items) + injected careless/distractors");
const RES={base:{wmcc:[],dfpr:[],auc:[]},acq:{wmcc:[],dfpr:[],auc:[]}};
for(let rep=0;rep<15;rep++){
  const idx=sampleK([...Array(baseAll.length).keys()],250), base=idx.map(i=>baseAll[i].slice());
  const {out,labels,dist,sev}=inject2(base,Jj,0.15,0.15);
  const X=featRows(out,Jj), w=sev.map(uWeight);
  for(const [name,cols] of [["base",[1,2,3]],["acq",[1,2,3,4]]]){
    const p=cvOOF(X,labels,cols), ow=oracleWMCC(p,labels,w), flag=p.map(v=>v>=ow.thr);
    let dt=0,df=0;for(let i=0;i<out.length;i++)if(dist[i]){dt++;if(flag[i])df++;}
    RES[name].wmcc.push(ow.mcc);RES[name].dfpr.push(dt?df/dt:0);RES[name].auc.push(auc(p,labels));
  }
}
console.log("  model         AUC    weighted-MCC   distractor-FPR");
console.log("  rr+LS+PT      "+mn(RES.base.auc).toFixed(3)+"  "+mn(RES.base.wmcc).toFixed(3)+"        "+mn(RES.base.dfpr).toFixed(3));
console.log("  +ACQUIESCENCE "+mn(RES.acq.auc).toFixed(3)+"  "+mn(RES.acq.wmcc).toFixed(3)+"        "+mn(RES.acq.dfpr).toFixed(3));
console.log("  Δ acq: AUC "+(mn(RES.acq.auc)-mn(RES.base.auc)>=0?"+":"")+(mn(RES.acq.auc)-mn(RES.base.auc)).toFixed(3)+
  ", distractor-FPR "+(mn(RES.acq.dfpr)-mn(RES.base.dfpr)).toFixed(3)+" (want negative)");

// ---------- Part B: collected study, no-regression ----------
console.log("\nPart B: collected study (n=157), LOO with/without acq — does acq HARM real careless detection?");
const Ms=R.parseCSV(fs.readFileSync("_study_matrix.csv","utf8")), Js=Ms.header.length;
const rowsS=Ms.rows.map(r=>r.map(v=>{const x=Number(v);return Number.isFinite(x)?x:NaN;}));
const yS=R.parseCSV(fs.readFileSync("_study_labels.csv","utf8")).rows.map(r=>Number(r[0]));
const XS=featRows(rowsS,Js), nS=XS.length;
function loo(cols){const P=cols.length+1,p=new Array(nS);
  for(let i=0;i<nS;i++){const tr=[];for(let j=0;j<nS;j++)if(j!==i)tr.push(j);
    const b=fit(tr.map(j=>[1,...cols.map(c=>XS[j][c])]),tr.map(j=>yS[j]),P);p[i]=pred(b,[1,...cols.map(c=>XS[i][c])]);}
  const rate=yS.reduce((a,b)=>a+b,0)/nS,srt=[...p].sort((a,b)=>a-b),thr=srt[Math.floor((1-rate)*nS)];
  return {auc:auc(p,yS),mcc:mcc(p.map(v=>v>=thr?1:0),yS)};}
const bB=loo([1,2,3]), aB=loo([1,2,3,4]);
console.log("  rr+LS+PT       LOO AUC "+bB.auc.toFixed(3)+"  MCC "+bB.mcc.toFixed(3));
console.log("  +ACQUIESCENCE  LOO AUC "+aB.auc.toFixed(3)+"  MCC "+aB.mcc.toFixed(3)+"   (Δ AUC "+(aB.auc-bB.auc>=0?"+":"")+(aB.auc-bB.auc).toFixed(3)+")");
fs.writeFileSync("test_acq_real.csv","which,model,auc,wmcc,dfpr\n"+
  "johnson,base,"+mn(RES.base.auc).toFixed(4)+","+mn(RES.base.wmcc).toFixed(4)+","+mn(RES.base.dfpr).toFixed(4)+"\n"+
  "johnson,acq,"+mn(RES.acq.auc).toFixed(4)+","+mn(RES.acq.wmcc).toFixed(4)+","+mn(RES.acq.dfpr).toFixed(4)+"\n"+
  "study,base,"+bB.auc.toFixed(4)+",,"+"\n"+"study,acq,"+aB.auc.toFixed(4)+",,\n");
console.log("\nwrote test_acq_real.csv");
