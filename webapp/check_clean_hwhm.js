/* check_clean_hwhm.js — false-alarm check for the HWHM sigma: a tighter scale estimate could
 * inflate flagging on genuinely clean data. Same generator as the earlier gate tests. */
const R = require("./site/reread.js"); const W = R.WEIGHTS;
let SEED = 31337; const rng = () => { SEED = (SEED*1103515245+12345)&0x7fffffff; return SEED/0x7fffffff; };
function gauss(){let u=0,v=0;while(!u)u=rng();while(!v)v=rng();return Math.sqrt(-2*Math.log(u))*Math.cos(2*Math.PI*v);}
const K5=5, THR=[-0.84,-0.25,0.25,0.84];
function genClean(F,ipf,n){const J=F*ipf;const load=[],fac=[],tau=[];
  for(let j=0;j<J;j++){let l=0.4+0.4*rng();if(rng()<0.4)l=-l;load.push(l);fac.push(Math.floor(j/ipf));tau.push(0.45*gauss());}
  const rho=0.3,sq=Math.sqrt(rho),sq1=Math.sqrt(1-rho);const rows=[];
  for(let i=0;i<n;i++){const c=gauss(),f=[];for(let k=0;k<F;k++)f.push(sq*c+sq1*gauss());
    const row=new Array(J);
    for(let j=0;j<J;j++){let y=tau[j]+load[j]*f[fac[j]];y+=Math.sqrt(Math.max(0.05,1-load[j]*load[j]))*gauss();
      let cat=1;for(let q=0;q<K5-1;q++)if(y>THR[q])cat=q+2;row[j]=cat;}
    rows.push(row);}
  return {rows,J};}
const toMat=(rows,J)=>{const n=rows.length,m=new Float64Array(n*J);for(let i=0;i<n;i++)for(let j=0;j<J;j++)m[i*J+j]=rows[i][j];return m;};
for (const [F,ipf,N,REPS] of [[10,6,200,30],[10,6,500,25],[20,6,2000,12],[30,6,500,15]]) {
  const J=F*ipf; let fire=0, fracSum=0, worst=0, estim=0;
  for (let r=0;r<REPS;r++){
    const {rows}=genClean(F,ipf,N);
    const mat=toMat(rows,J);
    const ef=R.ensembleFeatures(mat,N,J,{iterations:100,seed:1}),O=ef.oriented;
    const g=Math.max(0,Math.min(1,(ef.profileInfo-0.01)/0.05));
    const eta=new Array(N);
    for(let i=0;i<N;i++)eta[i]=W.b0+W.rr*O.rr[i]+g*(W.longstring*O.longstring[i]+W.person_total*O.person_total[i]);
    const af=R.autoFlag(eta,"standard");
    if(af.twoComp){fire++;fracSum+=af.estRate;if(af.estRate>worst)worst=af.estRate;}
    if(af.piEstimable)estim++;
  }
  console.log("clean J="+String(J).padStart(3)+" n="+String(N).padStart(5)+
    " | gate fires "+(100*fire/REPS).toFixed(0).padStart(3)+"%"+
    "  mean flagged (when fired) "+(fire?(100*fracSum/fire).toFixed(1):"0.0")+"%"+
    "  worst "+(100*worst).toFixed(1)+"%");
}
