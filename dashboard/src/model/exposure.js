/* Reduced particle-exposure model. All SI internally. No CFD fields or fitted slide results. */
const Exposure=(()=>{
 const dispersion=typeof ExposureDispersion!=='undefined'?ExposureDispersion:require('./dispersion.js');
 const numerics=typeof ExposureNumerics!=='undefined'?ExposureNumerics:require('./numerics.js');
 const network=typeof EductorNetwork!=='undefined'?EductorNetwork:typeof require==='function'?require('./network.js'):null;
 const upgrades=typeof EductorUpgrades!=='undefined'?EductorUpgrades:typeof require==='function'?require('../upgrades/physics.js'):null;
 const spatial=typeof EductorSpatial!=='undefined'?EductorSpatial:typeof require==='function'?require('../upgrades/spatial.js'):null;
 const collisions=typeof ExposureCollisions!=='undefined'?ExposureCollisions:typeof require==='function'?require('./collisions.js'):null;
 const population=typeof ExposurePopulation!=='undefined'?ExposurePopulation:require('./population.js');
 const doseKeys=collisions.KEYS,cumulativeKeys=['hazard','probability','screened','brownian','turbulent','unsupported'];
 const Rg=287.05,gamma=1.4,cp=gamma*Rg/(gamma-1),ATM=101325;
 const clamp=(x,a,b)=>Math.max(a,Math.min(b,x)),smooth=x=>(x=clamp(x,0,1),x*x*(3-2*x)),area=d=>Math.PI*d*d/4;
 const mu=T=>1.716e-5*(T/273.15)**1.5*(273.15+110.4)/(T+110.4);
 function particle(d,rhop,p,T,slip=0){
  const {mu:viscosity,rho,lambda}=numerics.gas(p,T),Kn=2*lambda/d;
  const Cc=1+Kn*(1.257+.4*Math.exp(-1.1/Kn)),Re=rho*d*Math.abs(slip)/viscosity,f=Re<1000?1+.15*Re**.687:.44*Re/24;
  return {tau:rhop*d*d*Cc/(18*viscosity*f),lambda,Cc,Re,Kn,rho,mu:viscosity,nu:viscosity/rho};
 }
 function nozzle(p0,p,T,d,Cd){
  if(!(p0>p&&p>0&&T>0))return null;
  const crit=(2/(gamma+1))**(gamma/(gamma-1)),choked=p/p0<=crit,pe=choked?p0*crit:p;
  const Te=T*(pe/p0)**((gamma-1)/gamma),ue=Math.sqrt(2*cp*(T-Te)),mp=Cd*area(d)*pe/(Rg*Te)*ue;
  return {mp,choked,uJet:Math.sqrt(2*cp*T*(1-(p/p0)**((gamma-1)/gamma))),power:mp*cp*T*(1-(p/p0)**((gamma-1)/gamma))};
 }
 function flow(g,p,legacy){
  p={...p,po:p.ps}; // One external reservoir pressure; internal stage pressures remain solved.
  if(![p.p0,p.ps,p.po,p.temp,p.feed].every(Number.isFinite)||p.ps<=0||p.po<=0||p.temp<=-273.15||p.feed<0)return {ok:false,reason:'Positive absolute pressures and temperature, and nonnegative feed, are required.'};
  const po=p.po*1e5,ps=p.ps*1e5,T=p.temp+273.15,p0=p.p0*1e5;
  if(p.mode==='prescribed'){
   const chamber=p.chamber*1e5,n=nozzle(p0,chamber,T,g.dn/1000,g.family==='british'?.94:.82);
   if(!n||!(p.air>0)||p.bleed<0)return {ok:false,reason:'Prescribed mode needs motive pressure above chamber pressure and positive feeder air.'};
   const ms=p.air/1000,mh=g.family==='british'?p.bleed/1000:0;
   return {ok:true,manual:true,caution:true,reason:'Prescribed air flows: the pressure / entrainment balance has not been solved.',...n,ms,mh,mg:n.mp+ms+mh,md:p.feed/3600,p0,ps,po,Tin:293.15,inletP:ps,portP:chamber,throatP:chamber,T0:T,T:(n.mp*T+(ms+mh)*293.15)/(n.mp+ms+mh)};
  }
  if(p.modelVersion===2){const f=upgrades.solve(g,p),bad=numerics.finiteScalars(f);return bad.length?{ok:false,constraints:['non-finite-output'],reason:'Non-finite flow output: '+bad.join(', ')}:f;}
  if(!(p0>po))return {ok:false,reason:'The existing solver requires motive stagnation pressure above outlet pressure.'};
  legacy.setBoundary(ps/ATM,po/ATM);let r,inletP=ps,trainNeeds=false;
  if(g.key.startsWith('nov'))r=legacy.solveCase(g,{suction:ps/ATM,outlet:po/ATM,pg:(p0-po)/1e5,feed:p.feed,temp:p.temp,loss:p.loss});
  else{
   const a=legacy.solveA(p0-po,p.feed,0,0,T),b=legacy.solveB(p0-po,p0-po,p.feed,0,0,T,0),c=legacy.solveC(p0-po,p.feed,0,T);
   const all=legacy.prepareResults(a,b,c);r=all[{revA:0,revB1:1,revB2:2,revC:3}[g.key]];if(g.key==='revB2')inletP=po+all[1].pt;
   trainNeeds=g.key.startsWith('revB')&&all.slice(1,3).some(v=>network.assess(v,{...g,db:v.st.db*1000,nports:1,df:v.st.db*1000},p).needed);
  }
  const gate=network.assess(r,g,p);
  if(gate.needed||trainNeeds){const f=network.solve(g,p);return {...f,gate:{needed:true,reason:trainNeeds?'At least one train stage needs a compressible estimate':gate.reason},raw:f.ok?network.legacyResult(f,g,r,po):r};}
  if(!r.supported)return {ok:false,reason:r.reason,raw:r};
  return {ok:true,manual:false,caution:r.caution,reason:r.caution?'Existing flow estimate is compressibility-sensitive.':'Existing reduced flow estimate passes its screening checks.',mp:r.mp,ms:r.ms,mh:r.mh||0,mg:r.mg,md:r.md,p0,ps,po:g.key==='revB1'?po+r.pe:po,systemOutletP:po,Tin:r.st?.Ts||293.15,inletP,portP:po+(r.pCh??r.p1),throatP:po+(r.pThroat??r.pMix??r.pe),T0:T,T:r.Tmix,uJet:r.n.uj,power:r.n.power,choked:r.n.choked,raw:r};
 }
 function geometry(g,release){
  const pts=g.profile.map(([x,r])=>[x/1000,r/1000]),end=g.end/1000,db=g.db/1000,dn=g.dn/1000;
  const angle=g.angle*Math.PI/180,sn=Math.sin(angle),cs=Math.cos(angle),R0=g.chamber/2000;
  function radius(x){x=clamp(x,0,end);let lo=0,hi=pts.length-1;while(hi-lo>1){const k=(lo+hi)>>1;if(pts[k][0]<x)lo=k;else hi=k;}const a=pts[lo],b=pts[hi],f=b[0]>a[0]?clamp((x-a[0])/(b[0]-a[0]),0,1):0;return a[1]+f*(b[1]-a[1]);}
  const cap=release/1000+.008;
  function clearance(x,y,dp=0){
   const h=1e-5,rr=radius(x),dr=(radius(x+h)-radius(x-h))/(2*h);
   const main=Math.min(x,(rr-Math.abs(y))/Math.sqrt(1+dr*dr),end-x)-dp/2;
   const s=-cs*x+sn*y,c=sn*x+cs*y;
   const branch=Math.min(s,cap-s,db/2-Math.abs(c))-dp/2;
   return Math.max(main,branch);
  }
  return {radius,end,db,dn,sn,cs,R0,cap,clearance,release:(off)=>[-cs*release/1000+sn*off,sn*release/1000+cs*off]};
 }
 function makeFieldCore(g,p,f){
  if(p.modelVersion===2)return spatial.makeField(upgrades.calibratedGeometry(g,p),p,f,geometry(g,p.release));
  const geo=geometry(g,p.release),b0=geo.dn/Math.sqrt(2),h=geo.dn*.035,bleedX=g.xh!=null?(g.xh-(g.xn||0))/1000:geo.end*.5;
  let noiseScale=1;
  function axial(x){
   x=clamp(x,0,geo.end);const r=geo.radius(x),b=b0+p.spread*x;
   const mixH=smooth((x-bleedX+geo.dn)/(2*geo.dn)),mass=f.mp+f.ms+f.mh*mixH;
   const totalT=(f.mp*f.T0+f.ms*f.Tin+f.mh*mixH*293.15)/mass;
   // Port -> mixing bore -> outlet interpolation; a reduced pressure field, not a nozzle/shock solution.
   const split=Math.min(g.L/1000,geo.end*.8),pressure=x<split?f.portP+(f.throatP-f.portP)*x/split:f.throatP+(f.po-f.throatP)*(x-split)/(geo.end-split);
   // Static temperature from adiabatic bulk energy including powder inertia.
   const energyC=(1+f.md/mass)/(2*cp)*(mass*Rg/(pressure*Math.PI*r*r))**2;
   const temp=f.compressible?2*totalT/(1+Math.sqrt(1+4*energyC*totalT)):totalT;
   const rho=pressure/(Rg*temp),mean=mass/(rho*Math.PI*r*r),fb=b*b/(r*r)*(-Math.expm1(-r*r/(b*b)));
   const amp=Math.min(f.uJet*b0/b,mean/Math.max(fb,1e-9));
   return {x,r,b,mass,temp,pressure,rho,mean,fb,amp};
  }
  const integral=(a,y)=>a.rho*(a.mean*y*y/2+a.amp*(a.b*a.b/2*(-Math.expm1(-y*y/(a.b*a.b)))-a.fb*y*y/2));
  function sample(x,y){
   const a=axial(x),eta=Math.abs(y)/a.b,uf=a.mean+a.amp*(Math.exp(-eta*eta)-a.fb);
   const am=axial(x-h),ap=axial(x+h),dx=Math.max(ap.x-am.x,1e-12);
   const radial=Math.abs(y)>1e-8?-(integral(ap,Math.abs(y))-integral(am,Math.abs(y)))/(dx*a.rho*Math.abs(y))*Math.sign(y):0;
   const s=-geo.cs*x+geo.sn*y,c=geo.sn*x+geo.cs*y,w=Math.abs(c)<geo.db/2&&s>0?smooth((s/geo.R0-.08)/.95):0;
   const feederRho=f.inletP/(Rg*f.Tin),us=f.ms/(g.nports*feederRho*area(geo.db));
   const ux=(1-w)*uf+w*geo.cs*us,uy=(1-w)*radial-w*geo.sn*us;
   const dev=1-Math.exp(-((Math.max(0,x)/(2*geo.dn))**2)),sigmaMain=p.intensity*a.amp*Math.sqrt(dev)*Math.exp(-(((eta-.7)/.65)**2));
   const sigma=noiseScale*((1-w)*sigmaMain+w*.15*p.intensity*us),ell=Math.max(geo.dn*.015,(1-w)*p.lengthScale*Math.min(a.b,a.r)+w*.1*geo.db);
   const k=1.5*sigma*sigma,eps=.9*k**1.5/ell;
   const pressure=(1-w)*a.pressure+w*f.inletP,temp=(1-w)*a.temp+w*f.Tin,nu=mu(temp)*Rg*temp/pressure;
   return {ux,uy,k,eps,G:Math.sqrt(eps/nu),sigma,ell,nu,p:pressure,T:temp,rho:pressure/(Rg*temp),powderConcentration:pressure/(Rg*temp)*f.md*((1-w)/a.mass+w/f.ms),b:a.b,r:a.r,TL:numerics.lagrangianTime(sigma,eps)};
  }
  // Bound the modeled turbulent dissipation by available motive mechanical power.
  let lossPower=0;const nx=60,nr=18,dx=geo.end/nx;
  for(let i=0;i<nx;i++){const x=(i+.5)*dx,r=geo.radius(x),dr=r/nr;for(let j=0;j<nr;j++){const y=(j+.5)*dr,v=sample(x,-y);lossPower+=v.rho*v.eps*2*Math.PI*y*dr*dx;}}
  const out=axial(geo.end),kinetic=.5*(f.mg+f.md)*out.mean*out.mean,available=Math.max(0,f.power-kinetic);
  if(lossPower>available&&lossPower>0)noiseScale=(available/lossPower)**(1/3);
  let peak=0,peakX=0,peakY=0,peakNu=mu(f.T)*Rg*f.T/f.po;
  for(let i=0;i<=180;i++){const x=geo.dn*1.001+(geo.end-geo.dn*1.001)*i/180;if(x<0)continue;const r=geo.radius(x);for(let j=0;j<=48;j++){const y=-r+(2*r)*j/48,v=sample(x,y);if(v.G>peak){peak=v.G;peakX=x;peakY=y;peakNu=v.nu;}}}
  return {geo,sample,axial,peak,peakX,peakY,peakNu,noiseScale,lossPower:lossPower*noiseScale**3,available,kinetic};
 }
 function makeField(g,p,f){const field=makeFieldCore(g,p,f),bad=numerics.finiteScalars(field);if(bad.length)throw Error('Non-finite field output: '+bad.join(', '));return field;}
 function random(seed){let x=seed>>>0;const u=()=>{x+=0x6D2B79F5;let t=x;t=Math.imul(t^t>>>15,t|1);t^=t+Math.imul(t^t>>>7,t|61);return ((t^t>>>14)>>>0)/4294967296;};return ()=>Math.sqrt(-2*Math.log(Math.max(u(),1e-12)))*Math.cos(2*Math.PI*u());}
 // Shared physical sections for dilution and particle shear. The first is the
 // release section in the suction branch; downstream sections are axial planes.
 function stationSections(field){
  const peak=field.peakX,end=field.geo.end;
  return [['Suction inlet',null],['Jet entry',0],['Peak shear',peak],['Mixing 7.5%',peak+.075*(end-peak),true],['Mixing 15%',peak+.15*(end-peak),true],['Mixing 22.5%',peak+.225*(end-peak),true],['Mixing 30%',peak+.3*(end-peak),true],['Midway',peak+.5*(end-peak)],['Mixing 70%',peak+.7*(end-peak),true],['Pre-outlet',peak+.85*(end-peak)],['Mixing 95%',peak+.95*(end-peak),true],['Outlet',end]].map(([label,x,diagnosticOnly],i)=>({id:i+1,label,x,...(diagnosticOnly?{diagnosticOnly:true}:{})}));
 }
 function sampleStationCrossings(field,d,sections,samples,a,b,dose=null){
  if(!(b.x>a.x))return;
  for(let i=1;i<sections.length-1;i++){
   const x=sections[i].x;if(samples[i]||x<a.x||x>b.x)continue;
   const w=(x-a.x)/(b.x-a.x),y=a.y+w*(b.y-a.y),z=(a.z||0)+w*((b.z||0)-(a.z||0));
   // Do not count a crossing of the same x-plane high in the suction branch.
   if(Math.hypot(y,z)>field.geo.radius(x)-d/2)continue;
   samples[i]={t:a.t+w*(b.t-a.t),x,y,z,G:field.sample(x,y,z).G,...(sections[i].diagnosticOnly?{diagnosticOnly:true}:{})};
   if(dose){
    const hazard=dose.hazard+integrateRate(dose.rate0,dose.rate1,b.t-a.t,w),screened=dose.screened+integrateRate(dose.screened0,dose.screened1,b.t-a.t,w);
    Object.assign(samples[i],{hazard,screened,probability:coagProbability(hazard,dose.alpha)});
    if(dose.components)for(const key of doseKeys)samples[i][key]=dose.components[key]+integrateRate(dose.start[key],dose.end[key],b.t-a.t,w);
   }
  }
 }
 function crossing(g0,g1,dt,threshold){if(Math.max(g0,g1)<threshold)return 0;if(Math.min(g0,g1)>=threshold)return dt;const f=(threshold-g0)/(g1-g0);return dt*(g1>g0?1-f:f);}
 function turnover(t,q,Gref,nu,d){const eta=Math.sqrt(nu/(q*Gref));return q*Gref*t*Math.min(1,(eta/d)**(2/3));}
 function track(field,p,d,offset,Gref,seed,noisy,keep=false,branchIndex=0){
  const geo=field.geo,rand=random(seed);let [x,y,z=0]=geo.release(offset,branchIndex),a=field.sample(x,y,z,branchIndex),vx=a.ux,vy=a.uy,vz=a.uz||0,t=0,wx=0,wy=0,wz=0;
  if(geo.clearance(x,y,d,z)<=0)throw Error('Release intersects the reduced passage wall; reduce the release offset.');
  const bins=241,horizon=p.horizon/1000,series=new Float64Array(bins),points=keep?[geo.dimension===3?[0,x,y,a.G,z]:[0,x,y,a.G]]:null;
  const sections=stationSections(field),stationSamples=sections.map(()=>null);stationSamples[0]={t:0,x,y,z,G:a.G,hazard:0,screened:0,probability:0,...Object.fromEntries(doseKeys.map(k=>[k,0]))};
  const hazardSeries=new Float64Array(bins),screenedSeries=new Float64Array(bins);let coag=collisionState(a,p,d),hazard=0,screened=0,shearDose=0,epsilonDose=0,maxSt=coag.maxSt,maxSize=coag.maxSize,maxMt=coag.Mt,maxConcentration=coag.concentration,maxVolume=coag.volumeFraction;
  const components=Object.fromEntries(doseKeys.map(k=>[k,0])),componentSeries=Object.fromEntries(doseKeys.map(k=>[k,new Float64Array(bins)]));let maxCollisionRe=coag.maxCollisionRe;
  const populationPath=p.populationBalance?[{t:0,a}]:null;let maxMach=0;
  let nextBin=1,half=0,quarter=0,maxRe=0,peak=0,steps=0,fate='censored',numericalReason=null,tSample=0,maxDt=0,minDt=Infinity;
  let normalized=[0,0,0];
  if(noisy&&p.dispersion==='well-mixed'&&p.initialFluctuation!=='zero'){normalized=normalized.map((_,i)=>numerics.normal(seed,i,0)*(p.antitheticSign??1));[wx,wy,wz]=normalized.map(q=>q*a.sigma);}
  const fluctuationScale=Math.max(a.sigma,Number.isFinite(field.peakX)?field.sample(field.peakX,field.peakY).sigma:0,1e-12);
  const bridge=p.commonNoise?numerics.brownian(seed,horizon,p.antitheticSign??1):null;
  series[0]=a.G;
  while(t<horizon&&steps<(p.maxSteps??2000000)){
   const speed=Math.max(Math.hypot(vx,vy,vz),Math.hypot(a.ux,a.uy,a.uz||0),1);
   const clearance=geo.clearance(x,y,d,z);
   if(clearance<=1e-9){const q=geo.dimension===3?geo.closest(x,y,z):null;fate=x>=geo.end-d/2-2e-6?'outlet':q&&q.s<=d/2+2e-6?'inlet':'wall';break;}
   let dt=numerics.timestep({maximum:p.dt*1e-6,Gref,G:a.G,b:a.b,r:a.r,db:geo.db,speed,remaining:horizon-t,TL:a.TL,noisy,clearance,adaptive:p.adaptiveDt===true});
   const grad=noisy&&p.dispersion==='well-mixed'?dispersion.gradients(field,x,y,z,a):null;
   if(grad&&p.inhomogeneousDt===true)dt=Math.min(dt,dispersion.stepLimit([wx,wy,wz],a,[vx-a.ux,vy-a.uy,vz-(a.uz||0)],grad));
   if(!(dt>1e-15)){fate='numerical';break;}
   maxDt=Math.max(maxDt,dt);minDt=Math.min(minDt,dt);
   const gaussian=axis=>bridge?bridge.increment(t,t+dt,axis)/Math.sqrt(dt):(p.antitheticSign??1)*rand();
   const ex=noisy?Math.exp(-dt/Math.max(a.TL,1e-12)):0,ns=noisy?a.sigma*Math.sqrt(-Math.expm1(-2*dt/Math.max(a.TL,1e-12))):0;
   const normalizedNext=grad?dispersion.advanceNormalized(normalized,a,[vx-a.ux,vy-a.uy,vz-(a.uz||0)],dt,gaussian,grad,{isotropic:p.referenceProtocol===true}):null;
   const mixed=normalizedNext?normalizedNext.map(q=>q*a.sigma):null;
   const nxw=mixed?mixed[0]:noisy?wx*ex+ns*gaussian(0):0,nyw=mixed?mixed[1]:noisy?wy*ex+ns*gaussian(1):0,nzw=mixed?mixed[2]:noisy&&geo.dimension===3?wz*ex+ns*gaussian(2):0,uz=(a.uz||0)+(wz+nzw)/2,ux=a.ux+(wx+nxw)/2,uy=a.uy+(wy+nyw)/2;
   // Reject catastrophic stochastic updates; never count their artificial
   // endpoint as a physical wall contact or exit. No velocity clipping.
   if(mixed&&(!normalizedNext.every(Number.isFinite)||Math.hypot(...mixed)>50*fluctuationScale)){fate='numerical';numericalReason='unresolved-dispersion-gradient';break;}
   const slip=Math.hypot(ux-vx,uy-vy,uz-vz),props=particle(d,p.density,a.p,a.T,slip),tau=props.tau;maxMach=Math.max(maxMach,slip/Math.sqrt(1.4*Rg*a.T));
   maxRe=Math.max(maxRe,props.Re);const ey=uy-9.81*tau,e=Math.exp(-dt/tau),om=-Math.expm1(-dt/tau);
   let xn=x+ux*dt+(vx-ux)*tau*om,yn=y+ey*dt+(vy-ey)*tau*om,zn=z+uz*dt+(vz-uz)*tau*om;
   let frac=1,hit=false;
   const contact=numerics.firstContact((xx,yy,zz)=>geo.clearance(xx,yy,d,zz),[x,y,z],[xn,yn,zn]);if(contact!==null){frac=contact;hit=true;}
   xn=x+frac*(xn-x);yn=y+frac*(yn-y);zn=z+frac*(zn-z);const tn=t+frac*dt,b=field.sample(xn,yn,zn,branchIndex),cb=collisionState(b,p,d,Math.hypot(b.ux-(ux+(vx-ux)*e),b.uy-(ey+(vy-ey)*e),(b.uz||0)-(uz+(vz-uz)*e))),elapsed=tn-t;
   sampleStationCrossings(field,d,sections,stationSamples,{t,x,y,z},{t:tn,x:xn,y:yn,z:zn},{hazard,screened,rate0:coag.rate,rate1:cb.rate,screened0:coag.screened,screened1:cb.screened,alpha:p.alpha??1,components,start:coag,end:cb});
   const nextHazard=hazard+integrateRate(coag.rate,cb.rate,elapsed),nextScreened=screened+integrateRate(coag.screened,cb.screened,elapsed);shearDose+=.5*(a.G+b.G)*elapsed;epsilonDose+=.5*(a.eps+b.eps)*elapsed;maxSt=Math.max(maxSt,cb.maxSt);maxSize=Math.max(maxSize,cb.maxSize);maxMt=Math.max(maxMt,cb.Mt);maxConcentration=Math.max(maxConcentration,cb.concentration);maxVolume=Math.max(maxVolume,cb.volumeFraction);maxCollisionRe=Math.max(maxCollisionRe,cb.maxCollisionRe);
   half+=crossing(a.G,b.G,tn-t,Gref/2);quarter+=crossing(a.G,b.G,tn-t,Gref/4);peak=Math.max(peak,a.G,b.G);
   while(nextBin<bins&&nextBin*horizon/(bins-1)<=tn){const tb=nextBin*horizon/(bins-1);const fracBin=(tb-t)/Math.max(tn-t,1e-15);series[nextBin]=a.G+(b.G-a.G)*fracBin;hazardSeries[nextBin]=hazard+integrateRate(coag.rate,cb.rate,elapsed,fracBin);screenedSeries[nextBin]=screened+integrateRate(coag.screened,cb.screened,elapsed,fracBin);for(const key of doseKeys)componentSeries[key][nextBin]=components[key]+integrateRate(coag[key],cb[key],elapsed,fracBin);nextBin++;}
   if(keep&&(tn>=tSample||hit)){points.push(geo.dimension===3?[tn,xn,yn,b.G,zn]:[tn,xn,yn,b.G]);tSample=tn+Math.max(horizon/1100,2e-6);}
   for(const key of doseKeys)components[key]+=integrateRate(coag[key],cb[key],elapsed);
   const populationPrevious=populationPath?{t,a}:null;
   vx=ux+(vx-ux)*Math.exp(-elapsed/tau);vy=ey+(vy-ey)*Math.exp(-elapsed/tau);vz=uz+(vz-uz)*Math.exp(-elapsed/tau);if(normalizedNext)normalized=normalizedNext;wx=normalizedNext?normalized[0]*b.sigma:nxw;wy=normalizedNext?normalized[1]*b.sigma:nyw;wz=normalizedNext?normalized[2]*b.sigma:nzw;x=xn;y=yn;z=zn;t=tn;a=b;coag=cb;hazard=nextHazard;screened=nextScreened;steps++;
   if(populationPath){
    const last=populationPath.at(-1),keys=['powderConcentration','G','eps','sigma','TL','p','T'],changed=(u,v)=>keys.some(k=>Math.abs((u[k]||0)-(v[k]||0))>(p.populationTolerance??.02)*Math.max(Math.abs(u[k]||0),Math.abs(v[k]||0),1e-20));
    if(hit||t-last.t>=(p.populationStep??20)*1e-6||changed(a,last.a)){
     // Preserve both sides of rapid field transitions; coarse uniform sampling
     // biases nonlinear kernels even when parcel integration conserves mass.
     if(changed(a,populationPrevious.a)&&populationPrevious.t>last.t)populationPath.push(populationPrevious);
     populationPath.push({t,a});
    }
   }
   if(hit){const q=geo.dimension===3?geo.closest(x,y,z):null,s=q?q.branch.length-q.s:-geo.cs*x+geo.sn*y;if(x>=geo.end-d/2-2e-6&&Math.hypot(y,z)<geo.radius(geo.end))fate='outlet';else if(q?q.s<=d/2+2e-6:s>=geo.cap-d/2-2e-6)fate='inlet';else fate='wall';break;}
  }
  if(steps>=(p.maxSteps??2000000)&&t<horizon)fate='numerical';
  if(keep&&points.at(-1)[0]!==t)points.push(geo.dimension===3?[t,x,y,a.G,z]:[t,x,y,a.G]);
  // A finite-size particle centre stops just inside the outlet plane.
  if(fate==='outlet')stationSamples[sections.length-1]={t,x,y,z,G:a.G,hazard,screened,probability:coagProbability(hazard,p.alpha??1),...components};
  hazardSeries.fill(hazard,nextBin);screenedSeries.fill(screened,nextBin);for(const key of doseKeys)componentSeries[key].fill(components[key],nextBin);
  if(populationPath&&populationPath.at(-1).t<t)populationPath.push({t,a});
  const growth=populationPath?evolvePopulation(populationPath,p,d,stationSamples):null;
  return {growth,maxMach,t,half,quarter,peak,maxRe,steps,maxDt,minDt:Number.isFinite(minDt)?minDt:0,meanDt:steps?t/steps:0,fate,numericalReason,contact:fate==='wall',deposited:null,branchIndex,dimension:geo.dimension||2,series,points,stationSamples,hazard,screened,probability:coagProbability(hazard,p.alpha??1),shearDose,epsilonDose,maxSt,maxSize,maxMt,maxConcentration,maxVolume,maxCollisionRe,...components,componentSeries,hazardSeries,screenedSeries};
 }
 function collisionOptions(p,slip=0){return {legacy:p.collisionModel==='legacy',tracerOnly:p.collisionModel==='tracer-brownian',slip,convention:p.kernelConvention||'zhou2001',extended:p.extendedKernels,clustering:p.clustering};}
 function evolvePopulation(path,p,d,stations=[]){
  const mixture=p.background==='mixture',ds=mixture?[p.small*1e-6,p.large*1e-6]:[d],shares=mixture?[p.smallShare??.5,1-(p.smallShare??.5)]:[1],pop=population.create(ds,shares,p.density,{bins:36});
  const stationGrowth=stations.map(()=>null);if(stations[0])stationGrowth[0]=population.summary(pop);const diameters=pop.m.map(m=>Math.cbrt(6*m/(p.density*Math.PI)));
  for(let j=1;j<path.length;j++){const start=path[j-1].t,end=path[j].t,breaks=[...stations.map((s,i)=>s&&!s.diagnosticOnly&&s.t>start&&s.t<=end?{t:s.t,index:i}:null).filter(Boolean),{t:end,index:null}].sort((a,b)=>a.t-b.t);let previous=start;for(const cut of breaks){
   for(let i=0;i<stations.length;i++){const station=stations[i];if(!station?.diagnosticOnly||station.t<=previous||station.t>cut.t)continue;const probe=structuredClone(pop),a=path[j-1].a,b=path[j].a,q={...a},fraction=(.5*(previous+station.t)-start)/Math.max(end-start,1e-30);for(const key of ['p','T','G','eps','sigma','TL','nu','rho','powderConcentration'])q[key]=a[key]+fraction*(b[key]-a[key]);population.advance(probe,Math.max(0,q.powderConcentration),station.t-previous,(i,k)=>collisions.pair(q,diameters[i],diameters[k],p.density,collisionOptions(p)).rate,p.alpha??1);stationGrowth[i]=population.summary(probe);}
   const a=path[j-1].a,b=path[j].a,q={...a},fraction=(.5*(previous+cut.t)-start)/Math.max(end-start,1e-30);for(const key of ['p','T','G','eps','sigma','TL','nu','rho','powderConcentration'])q[key]=a[key]+fraction*(b[key]-a[key]);population.advance(pop,Math.max(0,q.powderConcentration),cut.t-previous,(i,k)=>collisions.pair(q,diameters[i],diameters[k],p.density,collisionOptions(p)).rate,p.alpha??1);if(cut.index!==null)stationGrowth[cut.index]=population.summary(pop);previous=cut.t;}}
  return {...population.summary(pop),stations:stationGrowth,scope:'Sectional growth along a fixed transport parcel; growing sizes do not change the carrier field or parcel path.'};
 }
 function quantile(sorted,q){if(!sorted.length)return null;const x=(sorted.length-1)*q,i=Math.floor(x);return sorted[i]+(sorted[Math.min(i+1,sorted.length-1)]-sorted[i])*(x-i);}
 function stats(a){if(!a.length)return {n:0,mean:null,sd:null,se:null,ci:null,p05:null,p25:null,median:null,p75:null,p95:null,p99:null,p10:null,p90:null,cv:null,min:null,max:null};const n=a.length,mean=a.reduce((s,x)=>s+x,0)/n,sd=n>1?Math.sqrt(a.reduce((s,x)=>s+(x-mean)**2,0)/(n-1)):0,sorted=[...a].sort((a,b)=>a-b),se=sd/Math.sqrt(n),df=n-1,small=[0,12.7062,4.3027,3.1824,2.7764,2.5706,2.4469,2.3646,2.306,2.2622,2.2281],t=df<=10?small[df]:1.959964+2.372272/df+2.8225/(df*df)+2.556/(df**3);return {n,mean,sd,se,ci:n>1?[mean-t*se,mean+t*se]:null,p05:quantile(sorted,.05),p25:quantile(sorted,.25),median:quantile(sorted,.5),p75:quantile(sorted,.75),p95:quantile(sorted,.95),p99:quantile(sorted,.99),p10:quantile(sorted,.1),p90:quantile(sorted,.9),cv:mean!==0?sd/Math.abs(mean):null,min:sorted[0],max:sorted.at(-1)};}
 function weightedStats(values,weights){
  if(!values.length)return stats([]);
  const total=weights.reduce((a,b)=>a+b,0),w=weights.map(v=>v/total),mean=values.reduce((a,v,i)=>a+v*w[i],0),w2=w.reduce((a,b)=>a+b*b,0),effectiveN=1/w2,sd=Math.sqrt(values.reduce((a,v,i)=>a+w[i]*(v-mean)**2,0)/Math.max(1e-12,1-w2)),se=sd/Math.sqrt(effectiveN),sorted=values.map((v,i)=>({v,w:w[i]})).sort((a,b)=>a.v-b.v);
  const q=x=>{let c=0;for(const a of sorted){c+=a.w;if(c>=x)return a.v;}return sorted.at(-1).v;};
  const df=effectiveN-1,t=df>=1?1.959964+2.372272/df+2.8225/df**2+2.556/df**3:0;
  return {n:values.length,effectiveN,mean,sd,se,ci:effectiveN>1?[mean-t*se,mean+t*se]:null,p05:q(.05),p25:q(.25),median:q(.5),p75:q(.75),p95:q(.95),p99:q(.99),p10:q(.1),p90:q(.9),cv:mean?sd/Math.abs(mean):null,min:sorted[0].v,max:sorted.at(-1).v};
 }
 // A frozen-background encounter model. Full dose retains explicitly flagged
 // extrapolation; screened dose is a subtotal, never a corrected total.
 const collisionKernel=collisions.saffmanTurner;
 function collisionState(a,p,d,slip=0){
  const concentration=Math.max(0,a.powderConcentration||0),share=clamp(p.smallShare??.5,0,1),background=p.background==='mixture'?[[p.small*1e-6,share],[p.large*1e-6,1-share]]:[[d,1]],legacy=p.collisionModel==='legacy';
  const out={rate:0,screened:0,...Object.fromEntries(doseKeys.map(k=>[k,0])),maxSt:0,maxSize:0,maxCollisionRe:0,Mt:0,concentration,volumeFraction:concentration/p.density};
  for(const [dj,weight] of background){if(!weight)continue;const n=concentration*weight/(p.density*Math.PI*dj**3/6),q=collisions.pair(a,d,dj,p.density,collisionOptions(p,slip));
   for(const key of ['rate','screened',...doseKeys])out[key]+=q[key]*n;
   out.maxSt=Math.max(out.maxSt,q.maxSt);out.maxSize=Math.max(out.maxSize,q.maxSize);out.maxCollisionRe=Math.max(out.maxCollisionRe,q.maxRe);out.Mt=Math.max(out.Mt,q.Mt);
  }
  return out;
 }
 const coagProbability=(hazard,alpha=1)=>-Math.expm1(-clamp(alpha,0,1)*Math.max(0,hazard));
 const integrateRate=(a,b,dt,f=1)=>dt*(a*f+(b-a)*f*f/2);

 function wilson(k,n){const z=1.95996398454,d=1+z*z/n,c=(k/n+z*z/(2*n))/d,h=z*Math.sqrt((k/n)*(1-k/n)/n+z*z/(4*n*n))/d;return [Math.max(0,c-h),Math.min(1,c+h)];}
 function* simulation(g,p,f,options={}){
  const field=options.field||makeField(g,p,f),Gref=p.reference==='slide'?3829961.607121514:field.peak;
  if(!(Gref>0&&Number.isFinite(Gref)))return {ok:false,reason:'No positive turbulence reference: increase the turbulence intensity or provide an energetically feasible flow.'};
  const rows=[],offsets=[0],sizes=[p.small*1e-6,p.large*1e-6];
  for(let species=0;species<2;species++)for(let o=0;o<offsets.length;o++){
   if(options.groups&&!options.groups.includes(species*offsets.length+o))continue;
   const d=sizes[species],offset=offsets[o],adaptive=p.noise&&p.sampling==='precision',limit=p.noise?(adaptive?(p.maxHistories??4096):p.n):(field.geo.branches?.length||1),runs=[],mean=new Float64Array(241),sq=new Float64Array(241);let first;
   const fractions=field.branchFractions||[1],counts=fractions.map(()=>0),paired=p.noise&&p.antithetic===true;
   const seedFor=j=>(p.seed+species*103309+o*1009+(paired?Math.floor(j/2):j)*7919)>>>0;
   for(let j=0;j<limit;j++){
    const position=(((paired?Math.floor(j/2):j)+.5)*.618033988749895)%1;let c=0;const branch=!p.noise?j:(j<(paired?2:1)*fractions.length?Math.floor(j/(paired?2:1)):fractions.findIndex((v,i)=>(c+=v)>=position||i===fractions.length-1));
    const r=track(field,{...p,antitheticSign:paired&&j%2?-1:1},d,offset,Gref,seedFor(j),p.noise,counts[branch]===0,branch);counts[branch]++;runs.push(r);if(j===0)first=r;
    if((j+1)%(p.populationBalance||p.wallModel==='cad'?1:16)===0)yield {completed:(species*offsets.length+o)*limit+j+1,total:2*offsets.length*limit};
    if(adaptive&&(j+1)%128===0&&counts.every((n,i)=>n>0||fractions[i]===0)){
     runs.forEach(r=>r.weight=fractions[r.branchIndex]/counts[r.branchIndex]);
     const check=numerics.precision(runs,p.targetHalfWidth??.02);if(paired){check.halfWidth=(numerics.wilson(runs.reduce((s,r)=>s+(r.contact?r.weight:0),0),check.effectiveN/2)[1]-numerics.wilson(runs.reduce((s,r)=>s+(r.contact?r.weight:0),0),check.effectiveN/2)[0])/2;}
     if(check.halfWidth<=(p.targetHalfWidth??.02))break;
    }
   }
   const n=runs.length;for(const r of runs){r.weight=fractions[r.branchIndex]/counts[r.branchIndex];for(let b=0;b<241;b++){mean[b]+=r.series[b]*r.weight;sq[b]+=r.series[b]**2*r.weight;}}
   const summary=(items,fn)=>p.modelVersion===2?weightedStats(items.map(fn),items.map(r=>r.weight)):stats(items.map(fn));
   const sd=Array.from(mean,(v,i)=>n>1?Math.sqrt(Math.max(0,(sq[i]-v*v)/(p.modelVersion===2?1-runs.reduce((s,r)=>s+r.weight*r.weight,0):(n-1)/n))):0),wall=runs.filter(r=>r.fate==='wall').length;
   const sampleKeys=['t','half','quarter','peak','maxRe','maxDt','meanDt','maxMach','fate','hazard','screened','probability','shearDose','epsilonDose','maxSt','maxSize','maxMt','maxConcentration','maxVolume','maxCollisionRe',...doseKeys];
   const samples=runs.map((r,j)=>Object.fromEntries([['index',j+1],['seed',seedFor(j)],['antitheticSign',paired&&j%2?-1:1],['stationSamples',r.stationSamples],['branch',r.branchIndex],['powderFluxWeight',r.weight],['contact',r.contact],['deposited',r.deposited],['growth',r.growth],['numericalReason',r.numericalReason],...sampleKeys.map(k=>[k,r[k]])]));
   const stationShear=stationSections(field).map((section,i)=>({...section,total:n,...summary(runs.filter(r=>r.stationSamples[i]),r=>r.stationSamples[i].G)}));
   const stationCoagulation=Object.fromEntries(cumulativeKeys.map(key=>[key,stationSections(field).map((section,i)=>({...section,total:n,...summary(runs.filter(r=>r.stationSamples[i]),r=>r.stationSamples[i][key])}))]));
   const distribution=Object.fromEntries(sampleKeys.filter(k=>k!=='fate').map(k=>[k,summary(runs,r=>r[k])]));distribution.outletTime=summary(runs.filter(r=>r.fate==='outlet'),r=>r.t);
   const residence=numerics.residence(samples),precision=numerics.precision(runs,p.targetHalfWidth??.02);if(paired){const ci=numerics.wilson(runs.reduce((s,r)=>s+(r.contact?r.weight:0),0),precision.effectiveN/2);Object.assign(precision,{effectiveN:precision.effectiveN/2,ci,halfWidth:(ci[1]-ci[0])/2,met:(ci[1]-ci[0])/2<=precision.target,method:'Conservative pair-count Wilson interval'});}
   const cumulative={},traceValues=(r,key)=>key==='probability'?Array.from(r.hazardSeries,v=>coagProbability(v,p.alpha??1)):Array.from(key==='screened'?r.screenedSeries:key==='hazard'?r.hazardSeries:r.componentSeries[key]);
   for(const key of cumulativeKeys){
    const c=Object.fromEntries(['mean','median','sd','se','p05','p25','p75','p95','ciLow','ciHigh'].map(k=>[k,[]])),histories=runs.map(r=>traceValues(r,key));
    for(let b=0;b<241;b++){const q=p.modelVersion===2?weightedStats(histories.map(v=>v[b]),runs.map(r=>r.weight)):stats(histories.map(v=>v[b]));for(const k of ['mean','median','sd','se','p05','p25','p75','p95'])c[k].push(q[k]);c.ciLow.push(q.ci?.[0]??null);c.ciHigh.push(q.ci?.[1]??null);if(b%60===0)yield {completed:(species*offsets.length+o+1)*n,total:2*offsets.length*n,phase:'statistics'};}
    c.traces=histories.slice(0,12);c.convergence=[...new Set([Math.max(1,Math.floor(n/4)),Math.max(1,Math.floor(n/2)),n])].map(count=>({count,...summary(runs.slice(0,count),r=>r[key])}));cumulative[key]=c;
   }
   const fates=Object.fromEntries(['wall','outlet','inlet','censored','numerical'].map(k=>{const subset=runs.filter(r=>r.fate===k),count=subset.length,fraction=p.modelVersion===2?subset.reduce((s,r)=>s+r.weight,0):count/n,neff=(p.modelVersion===2?1/runs.reduce((s,r)=>s+r.weight*r.weight,0):n)/(paired?2:1);return [k,{count,fraction,ci:wilson(fraction*neff,neff)}];}));
   rows.push({species,d,offset,n,residence,precision,branchFractions:fractions,branchCounts:counts,paths:runs.filter(r=>r.points).map(r=>({points:r.points,branch:r.branchIndex,fate:r.fate})),samples,stationShear,stationCoagulation,distribution,cumulative,fates,time:summary(runs,r=>r.t),half:summary(runs,r=>r.half),quarter:summary(runs,r=>r.quarter),nh:summary(runs,r=>turnover(r.half,.5,Gref,field.peakNu,d)),nq:summary(runs,r=>turnover(r.quarter,.25,Gref,field.peakNu,d)),growth:runs.some(r=>r.growth)?{numberLoss:summary(runs,r=>r.growth.numberLoss),meanDiameter:summary(runs,r=>r.growth.meanDiameter),maxMassError:Math.max(...runs.map(r=>r.growth.massError)),maxOverflow:Math.max(...runs.map(r=>r.growth.overflowMassFraction))}:null,wall:runs.reduce((s,r)=>s+(r.contact?r.weight:0),0),wallCI:precision.ci,outlet:fates.outlet.fraction,inlet:fates.inlet.fraction,censored:fates.censored.fraction,maxRe:Math.max(...runs.map(r=>r.maxRe)),steps:Math.max(...runs.map(r=>r.steps)),curve:Array.from(mean),sd,pathFate:first.fate,path:first.points,meanPath:p.noise?track(field,p,d,offset,Gref,1,false,true).points:first.points});
  }
  return {ok:true,field,Gref,rows,etaHalf:Math.sqrt(field.peakNu/(Gref/2)),etaQuarter:Math.sqrt(field.peakNu/(Gref/4)),pressure:f,settings:{...p}};
 }
 function simulate(g,p,f,options={}){const it=simulation(g,p,f,options);let step;do{step=it.next();}while(!step.done);return step.value;}
 async function simulateAsync(g,p,f,{cancelled=()=>false,onProgress=()=>{},field,groups}={}){const it=simulation(g,p,f,{field,groups});while(true){if(cancelled()){it.return();return {ok:false,cancelled:true};}const step=it.next();if(step.done)return step.value;onProgress(step.value);await new Promise(resolve=>{if(typeof MessageChannel==='undefined'){setTimeout(resolve,0);return;}const c=new MessageChannel();c.port1.onmessage=()=>{c.port1.close();c.port2.close();resolve();};c.port2.postMessage(0);});}}
 return {numerics,dispersion,population,evolvePopulation,particle,nozzle,flow,geometry,makeField,stationSections,sampleStationCrossings,track,crossing,turnover,wilson,simulate,simulateAsync,stats,weightedStats,quantile,collisionKernel,collisionState,collisions,coagProbability,integrateRate};
})();
if(typeof module!=='undefined')module.exports=Exposure;

// One independent reference, including frozen views. Never enters the solver.
