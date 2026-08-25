const fs=require("fs");const R=require("./site/reread.js");
function auc(sc,y){const idx=[...sc.keys()].sort((a,b)=>sc[a]-sc[b]);const rk=[];idx.forEach((id,r)=>rk[id]=r+1);
  let n1=0,n0=0,s=0;for(let i=0;i<sc.length;i++){if(y[i]===1){n1++;s+=rk[i];}else n0++;}return n1&&n0?(s-n1*(n1+1)/2)/(n1*n0):NaN;}
function metrics(flag,y){let tp=0,fp=0,fn=0,tn=0;for(let i=0;i<y.length;i++){if(y[i]===1){flag[i]?tp++:fn++;}else{flag[i]?fp++:tn++;}}
  const prec=tp+fp?tp/(tp+fp):0,rec=tp+fn?tp/(tp+fn):0,spec=tn+fp?tn/(tn+fp):1;
  const d=Math.sqrt((tp+fp)*(tp+fn)*(tn+fp)*(tn+fn)),mcc=d?(tp*tn-fp*fn)/d:0,f1=prec+rec?2*prec*rec/(prec+rec):0;
  return {tp,fp,fn,tn,prec,rec,spec,mcc,f1,fpr:fp/(tn+fp||1)};}
const pc=x=>(100*x).toFixed(0)+"%", f2=x=>x.toFixed(2);

// ---------- STUDY (real ground truth) ----------
const Ms=R.parseCSV(fs.readFileSync("_study_matrix.csv","utf8"));
const ns=Ms.rows.length,Js=Ms.header.length,ms=new Float64Array(ns*Js);
for(let i=0;i<ns;i++)for(let j=0;j<Js;j++)ms[i*Js+j]=R.toNum(Ms.rows[i][j]);
const ys=R.parseCSV(fs.readFileSync("_study_labels.csv","utf8")).rows.map(r=>Number(r[0]));
const trueRate=ys.reduce((s,v)=>s+v,0)/ns;
const rA=R.ensemble(ms,ns,Js,{iterations:400,seed:1});
const rO=R.ensemble(ms,ns,Js,{iterations:400,seed:1,prevalence:trueRate});
const mA=metrics(rA.flagged,ys),mO=metrics(rO.flagged,ys);
console.log("=== STUDY (real, n="+ns+", true careless="+ys.reduce((s,v)=>s+v,0)+"/"+ns+" = "+pc(trueRate)+") ===");
console.log("  AUC(ensemble) = "+f2(auc(rA.p,ys))+"   [rr alone AUC = "+f2(auc(rA.rr.map(v=>-v),ys))+"]");
console.log("  AUTO  : 2comp="+rA.twoComponent+" est="+pc(rA.estimatedRate)+" flagged="+rA.nFlagged+
            " | prec="+f2(mA.prec)+" rec="+f2(mA.rec)+" spec="+f2(mA.spec)+" MCC="+f2(mA.mcc)+" F1="+f2(mA.f1));
console.log("  ORACLE(@"+pc(trueRate)+"): flagged="+rO.nFlagged+
            " | prec="+f2(mO.prec)+" rec="+f2(mO.rec)+" spec="+f2(mO.spec)+" MCC="+f2(mO.mcc)+" F1="+f2(mO.f1));

// ---------- SIMULATED: inject careless patterns on real IPIP base ----------
const base=R.parseCSV(fs.readFileSync("demo_real_base.csv","utf8"),",");
const J=base.header.length-1, rows0=base.rows.map(r=>r.slice(1).map(Number)), N=rows0.length;
function mat(rr){const m=new Float64Array(rr.length*J);for(let i=0;i<rr.length;i++)for(let j=0;j<J;j++)m[i*J+j]=rr[i][j];return m;}
function rng(x){let a=x>>>0;return()=>{a|=0;a=(a+0x6D2B79F5)|0;let t=Math.imul(a^(a>>>15),1|a);t=(t+Math.imul(t^(t>>>7),61|t))^t;return((t^(t>>>14))>>>0)/4294967296;};}
function inject(seed,rate){const rand=rng(seed);const out=rows0.map(r=>r.slice());const y=new Array(N).fill(0);
  const ids=[...Array(N).keys()];for(let i=N-1;i>0;i--){const j=Math.floor(rand()*(i+1));[ids[i],ids[j]]=[ids[j],ids[i]];}
  const M=Math.round(N*rate),pats=["random","long","straight","acq","mixed","fatigue"],LV=[1,2,3,4,5];
  const sk=k=>{const a=[...Array(J).keys()];for(let i=0;i<k;i++){const j=i+Math.floor(rand()*(J-i));[a[i],a[j]]=[a[j],a[i]];}return a.slice(0,k);},pick=a=>a[Math.floor(rand()*a.length)];
  for(let m=0;m<M;m++){const i=ids[m],lvl=[.5,.6,.7,.8,.9,1][Math.floor(rand()*6)],k=Math.max(1,Math.round(J*lvl)),pat=pats[m%6],row=out[i];
    if(pat==="random")sk(k).forEach(p=>row[p]=pick(LV));
    else if(pat==="straight"){const v=pick(LV);sk(k).forEach(p=>row[p]=v);}
    else if(pat==="acq")sk(k).forEach(p=>row[p]=Math.max(1,Math.min(5,4+pick([-1,0,0,0,1]))));
    else if(pat==="long"){const idx=sk(k);let v=pick(LV);idx.forEach((p,t)=>{if(t%12===0)v=pick(LV);row[p]=v;});}
    else if(pat==="mixed"){const idx=sk(k),h=Math.round(k/2);let v=pick(LV);idx.slice(0,h).forEach((p,t)=>{if(t%10===0)v=pick(LV);row[p]=v;});idx.slice(h).forEach(p=>row[p]=pick(LV));}
    else{const idx=[];for(let t=J-k;t<J;t++)idx.push(t);let v=pick(LV);idx.forEach((p,t)=>{if(t%10===0)v=pick(LV);row[p]=v;});}
    y[i]=1;}
  return {out,y};}
console.log("\n=== SIMULATED: inject-and-detect on real IPIP base (N="+N+"×"+J+" items, 10 reps/rate) ===");
console.log(" true    AUC    AUTO: est / flagged / prec / rec / MCC        ORACLE@true: prec / rec / MCC");
for(const rate of [0,0.05,0.10,0.20,0.40]){
  const A={est:[],prec:[],rec:[],mcc:[],auc:[],fpr:[]},O={prec:[],rec:[],mcc:[]};
  for(let rep=0;rep<10;rep++){const {out,y}=inject(1000+rep*7+Math.round(rate*100),rate);
    const ra=R.ensemble(mat(out),N,J,{iterations:200,seed:1});const ma=metrics(ra.flagged,y);
    A.est.push(ra.estimatedRate);A.auc.push(rate>0?auc(ra.p,y):NaN);A.fpr.push(ma.fpr);
    if(rate>0){A.prec.push(ma.prec);A.rec.push(ma.rec);A.mcc.push(ma.mcc);
      const ro=R.ensemble(mat(out),N,J,{iterations:200,seed:1,prevalence:rate});const mo=metrics(ro.flagged,y);
      O.prec.push(mo.prec);O.rec.push(mo.rec);O.mcc.push(mo.mcc);}}
  const mn=a=>a.length?a.reduce((s,v)=>s+v,0)/a.length:NaN;
  if(rate===0)console.log(" "+pc(rate).padStart(4)+"    -     est="+pc(mn(A.est))+"  flagged FPR="+f2(mn(A.fpr))+" (clean → should flag ~0)");
  else console.log(" "+pc(rate).padStart(4)+"   "+f2(mn(A.auc))+"   est="+pc(mn(A.est))+" / rec="+f2(mn(A.rec))+" / prec="+f2(mn(A.prec))+" / MCC="+f2(mn(A.mcc))+
                   "      "+f2(mn(O.prec))+" / "+f2(mn(O.rec))+" / "+f2(mn(O.mcc)));
}
