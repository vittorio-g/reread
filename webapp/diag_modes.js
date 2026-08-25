/* Does the careless component form its own KDE mode? Counts prominent peaks of the eta
 * distribution and checks whether the right-most peak coincides with the attentive one. */
const R=require("./site/reread.js"); const W=R.WEIGHTS; const fs=require("fs");
let SEED=99; const rng=()=>{SEED=(SEED*1103515245+12345)&0x7fffffff;return SEED/0x7fffffff;};
const mn=a=>a.reduce((s,v)=>s+v,0)/a.length;
const med=a=>{const s=[...a].sort((x,y)=>x-y);const h=s.length>>1;return s.length%2?s[h]:(s[h-1]+s[h])/2;};
function peaks(x){const n=x.length,s=[...x].sort((a,b)=>a-b);
  const q1=s[Math.floor(.25*n)],q3=s[Math.floor(.75*n)],iqr=q3-q1;
  const mu=mn(x), sd=Math.sqrt(x.reduce((a,b)=>a+(b-mu)*(b-mu),0)/n);
  let h=0.9*Math.min(sd||1,(iqr/1.34)||(sd||1))*Math.pow(n,-0.2); if(!(h>0))h=(sd||1)*0.3;
  const G=512,lo=s[0]-h,hi=s[n-1]+h,d=new Float64Array(G);
  for(let g=0;g<G;g++){const xg=lo+(hi-lo)*g/(G-1);let v=0;for(let i=0;i<n;i++){const z=(xg-x[i])/h;v+=Math.exp(-0.5*z*z);}d[g]=v;}
  let dm=0;for(let g=0;g<G;g++)if(d[g]>dm)dm=d[g];
  const P=[];for(let g=1;g<G-1;g++) if(d[g]>=d[g-1]&&d[g]>=d[g+1]&&d[g]>0.15*dm) P.push(lo+(hi-lo)*g/(G-1));
  return P;}
const Ms=R.parseCSV(fs.readFileSync("_study_matrix.csv","utf8")),J=Ms.header.length;
const allR=Ms.rows.map(r=>r.map(Number));
const y=R.parseCSV(fs.readFileSync("_study_labels.csv","utf8")).rows.map(r=>Number(r[0]));
const cf=y.map((v,i)=>v===0?i:-1).filter(i=>i>=0), cl=y.map((v,i)=>v===1?i:-1).filter(i=>i>=0);
console.log("true |  mean#peaks  single-peak%  | rightPeak vs careful-mean (SD units)");
for(const k of [4,12,24,47,58]){ const tr=k/(cf.length+k); const np=[],sp=[],dist=[];
  for(let rep=0;rep<8;rep++){ const a=cl.slice();
    for(let i=0;i<k;i++){const j=i+Math.floor(rng()*(a.length-i));const t=a[i];a[i]=a[j];a[j]=t;}
    const idx=[...cf,...a.slice(0,k)],n=idx.length,mat=new Float64Array(n*J);
    for(let i=0;i<n;i++)for(let j=0;j<J;j++)mat[i*J+j]=allR[idx[i]][j];
    const ef=R.ensembleFeatures(mat,n,J,{iterations:120,seed:1}),O=ef.oriented;
    const g=Math.max(0,Math.min(1,(ef.profileInfo-0.01)/0.05));
    const e=[];for(let i=0;i<n;i++)e.push(W.b0+W.rr*O.rr[i]+g*(W.longstring*O.longstring[i]+W.person_total*O.person_total[i]));
    const P=peaks(e); np.push(P.length); sp.push(P.length<=1?1:0);
    const eCf=e.slice(0,cf.length), mCf=mn(eCf), sCf=Math.sqrt(mn(eCf.map(v=>(v-mCf)*(v-mCf))))||1;
    dist.push((P[P.length-1]-mCf)/sCf); }
  console.log((100*tr).toFixed(0).padStart(4)+"% |   "+mn(np).toFixed(2)+"        "+(100*mn(sp)).toFixed(0).padStart(3)+"%       |  "+mn(dist).toFixed(2));}
