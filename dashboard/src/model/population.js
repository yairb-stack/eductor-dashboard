const ExposurePopulation=(()=>{
 const sum=a=>a.reduce((s,v)=>s+v,0);
 function create(diameters,shares,density,{bins=36,ratio=2}={}){
  const mass=d=>density*Math.PI*d**3/6,min=Math.min(...diameters.map(mass)),max=Math.max(...diameters.map(mass))*128,r=Math.max(ratio,(max/min)**(1/(bins-1))),m=Array.from({length:bins},(_,i)=>min*r**i),q=new Float64Array(bins);
  diameters.forEach((d,i)=>deposit(q,m,mass(d),shares[i]/mass(d)));
  return {m,q,density,initialNumber:sum(q),initialMass:1,overflowNumber:0,overflowMass:0,steps:0,maxMassError:0};
 }
 function deposit(q,m,v,n){if(v>m.at(-1))return {number:n,mass:n*v};let k=0;while(k<m.length-2&&m[k+1]<v)k++;const w=(m[k+1]-v)/(m[k+1]-m[k]);q[k]+=n*w;q[k+1]+=n*(1-w);return {number:0,mass:0};}
 function advance(s,C,dt,kernel,alpha=1){
  if(!(C>0&&dt>0&&alpha>0))return s;
  const n=s.m.length,K=Array.from({length:n},()=>new Map()),getK=(i,j)=>{if(!K[i].has(j)){const v=Math.max(0,alpha*kernel(i,j));K[i].set(j,v);K[j].set(i,v);}return K[i].get(j);};
  const births=Array.from({length:n},()=>new Map());
  function birth(i,j){if(births[i].has(j))return births[i].get(j);const v=s.m[i]+s.m[j];let out;
   if(v>s.m.at(-1))out={bins:[],mass:v};else{let l=0;while(l<n-2&&s.m[l+1]<v)l++;const w=(s.m[l+1]-v)/(s.m[l+1]-s.m[l]);out={bins:[[l,w],[l+1,1-w]],mass:0};}
   out.lossI=(i===j?2:1)-(out.bins.find(b=>b[0]===i)?.[1]||0);out.lossJ=i===j?out.lossI:1-(out.bins.find(b=>b[0]===j)?.[1]||0);births[i].set(j,out);return out;
  }
  function react(i,j,h){const k=C*getK(i,j),a=s.q[i],b=s.q[j];if(!(a>0&&b>0&&k>0))return;const out=birth(i,j),si=out.lossI,sj=out.lossJ;let event;
   if(i===j){const z=.5*k*si*a*h;event=a/si*(z/(1+z));}
   else if(sj<1e-15)event=a/si*(-Math.expm1(-si*k*b*h));
   else {const aa=a/si,bb=b/sj,small=Math.min(aa,bb),large=Math.max(aa,bb),delta=large-small,kk=k*si*sj,z=kk*delta*h;event=delta<large*1e-10?small*(kk*large*h/(1+kk*large*h)):small*large*(-Math.expm1(-z))/(delta+small*(-Math.expm1(-z)));}
   s.q[i]-=si*event;if(j!==i)s.q[j]-=sj*event;
   for(const [bin,w] of out.bins)if(bin!==i&&bin!==j)s.q[bin]+=w*event;
   if(out.mass){s.overflowNumber+=event;s.overflowMass+=event*out.mass;}
  }
  let left=dt;while(left>dt*1e-12){
   const active=[];for(let i=0;i<n;i++)if(s.q[i]*s.m[i]>1e-14*s.initialMass)active.push(i);
   let rate=0,N=sum(s.q),resolvedLoss=0;
   for(const i of active){let loss=0;for(const j of active){const pair=i<=j?birth(i,j):birth(j,i),encounter=C*getK(i,j)*s.q[j];rate+=s.q[i]*encounter;loss+=encounter*(i===j?.5*pair.lossI:i<j?pair.lossI:pair.lossJ);}if(s.q[i]*s.m[i]>1e-8*s.initialMass)resolvedLoss=Math.max(resolvedLoss,loss);}
   const h=Math.min(left,rate>0?.02*N/rate:left,resolvedLoss>0?.05/resolvedLoss:left),pairs=[];
   for(let a=0;a<active.length;a++)for(let b=a;b<active.length;b++)pairs.push([active[a],active[b]]);
   for(const [i,j] of pairs)react(i,j,h/2);for(let k=pairs.length-1;k>=0;k--)react(...pairs[k],h/2);
   s.steps++;left-=h;if(s.steps>1000000)throw Error('Population balance exceeded its substep budget');
  }
  const total=s.q.reduce((v,q,i)=>v+q*s.m[i],s.overflowMass);s.maxMassError=Math.max(s.maxMassError,Math.abs(total-s.initialMass)/s.initialMass);return s;
 }
 function summary(s){const number=sum(s.q)+s.overflowNumber,mass=s.q.reduce((v,q,i)=>v+q*s.m[i],s.overflowMass);return {numberPerKg:number,numberFraction:number/s.initialNumber,numberLoss:1-number/s.initialNumber,mass,massError:s.maxMassError,meanDiameter:Math.cbrt(6*mass/(s.density*Math.PI*number)),overflowMassFraction:s.overflowMass/mass,distribution:s.m.map((m,i)=>({diameter:Math.cbrt(6*m/(s.density*Math.PI)),numberPerKg:s.q[i],massFraction:m*s.q[i]})),steps:s.steps};}
 return {create,advance,summary};
})();
if(typeof module!=='undefined')module.exports=ExposurePopulation;

'use strict';
