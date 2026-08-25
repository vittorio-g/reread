const R=require("./site/reread.js");const fs=require("fs");
const rng=R._rng(20260713);
const src=fs.readFileSync("sim_shipped.js","utf8");
eval(src.slice(src.indexOf("function gauss"), src.indexOf("function mcc")));  // genClean, inject, qnorm, gauss, K, THR
function met(f,y){let tp=0,fp=0,fn=0,tn=0;for(let i=0;i<y.length;i++){if(y[i]){f[i]?tp++:fn++;}else{f[i]?fp++:tn++;}}
  const rec=tp+fn?tp/(tp+fn):0,prec=tp+fp?tp/(tp+fp):0,d=Math.sqrt((tp+fp)*(tp+fn)*(tn+fp)*(tn+fn)),mcc=d?(tp*tn-fp*fn)/d:0;return{rec,prec,mcc};}
const mn=a=>a.reduce((x,y)=>x+y,0)/a.length;
console.log("EXACT website pipeline (default opts: corProp .03, iter 200, Standard), 12 datasets/length, 20% careless\n");
for(const [F,ipf] of [[10,6],[20,6],[30,6],[50,6]]){  // 60,120,180,300
  const J=F*ipf, recs=[],precs=[],mccs=[],ests=[],fails=[];
  for(let rep=0;rep<12;rep++){
    const {rows}=genClean(F,ipf,500);const {out,labels}=inject(rows,J,0.20);
    // write CSV exactly like a user file: id + item cols
    let csv="id,"+Array.from({length:J},(_,j)=>"item_"+(j+1)).join(",")+"\n";
    csv+=out.map((r,i)=>"R"+(i+1)+","+r.join(",")).join("\n");
    // === identical to worker.js ===
    const parsed=R.parseCSV(csv);
    const det=R.detectColumns(parsed.header,parsed.rows);
    const {M,n,J:JJ}=R.buildMatrix(parsed.rows,det.itemCols);
    const res=R.ensemble(M,n,JJ,{corProp:0.03,iterations:200,minPairs:15,sensitivity:"standard"});
    const m=met(res.flagged,labels);
    recs.push(m.rec);precs.push(m.prec);mccs.push(m.mcc);ests.push(res.estimatedRate);
    if(m.mcc<0.3)fails.push(rep);
  }
  console.log("J="+String(J).padStart(3)+"  MCC: mean="+mn(mccs).toFixed(2)+" min="+Math.min(...mccs).toFixed(2)+" max="+Math.max(...mccs).toFixed(2)+
    "  | est-rate mean="+(100*mn(ests)).toFixed(0)+"% (true 20%)  | recall="+mn(recs).toFixed(2)+" prec="+mn(precs).toFixed(2)+
    "  | reps with MCC<0.3: "+fails.length+"/12");
  console.log("     per-rep MCC: ["+mccs.map(v=>v.toFixed(2)).join(", ")+"]");
}
