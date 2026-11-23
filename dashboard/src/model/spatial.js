/* CAD-centreline 3D tube union + conservative dilute scalar reconstruction.
 * The internal junctions and turbulent field remain reduced approximations.
 */
const EductorSpatial=(()=>{
 const WALL=typeof EductorWalls!=='undefined'?EductorWalls:require('./walls.js');
 const JET=typeof EductorJet!=='undefined'?EductorJet:require('./jet.js');
 const NUM=typeof ExposureNumerics!=='undefined'?ExposureNumerics:require('../exposure/numerics.js');
 const U=typeof EductorUpgrades!=='undefined'?EductorUpgrades:require('./physics.js');
 const {clamp,area}=U,R=287.05,cp=1004.675;
 const smooth=x=>(x=clamp(x,0,1),x*x*(3-2*x));
 const mu=T=>1.716e-5*(T/273.15)**1.5*(273.15+110.4)/(T+110.4);
 function geometry(g,p,legacy){
  const source=U.topology(g,p.release),cad=!!U.CAD[g.key],convert=q=>[q[0]/1000,q[2]/1000,q[1]/1000];
  const branches=source.branches.map((b,i)=>{
   let pts=b.points.filter((v,j)=>j===0||j===b.points.length-1||j%3===0).map(convert);
   if(source.common.lengthMM>0)pts=[convert(source.common.points[0]),...pts];
   const segments=[];let length=0;
   for(let j=0;j<pts.length-1;j++){const a=pts[j],b=pts[j+1],delta=b.map((v,k)=>v-a[k]),L=Math.hypot(...delta);if(L<1e-10)continue;segments.push({a,b,L,s:length,tangent:delta.map(v=>v/L),radius:(j===0&&source.common.lengthMM>0?source.common.radius:source.branches[i].radius)/1000});length+=L;}
   const last=segments.at(-1),entry=last.b;
   return {id:i,segments,length,entry,angle:Math.atan2(last.a[2],last.a[1]),radius:b.radius/1000};
  });
  function suctionRadius(y){const pts=g.suctionProfile;if(!pts||p.flowPhysics!=='phase2')return null;const mm=clamp(y*1000,pts[0][0],pts.at(-1)[0]);let i=1;while(i<pts.length-1&&pts[i][0]<mm)i++;const a=pts[i-1],b=pts[i],slope=(b[1]-a[1])/(b[0]-a[0]);return {radius:(a[1]+slope*(mm-a[0]))/1000,slope:y*1000<=pts[0][0]||y*1000>=pts.at(-1)[0]?0:slope};}
  function closest(x,y,z=0,only=null){
   let best=null,bestClear=-Infinity;
   for(let bi=0;bi<branches.length;bi++){if(only!==null&&bi!==only)continue;const branch=branches[bi];for(const segment of branch.segments){
    const a=segment.a,v=segment.tangent,t=clamp((x-a[0])*v[0]+(y-a[1])*v[1]+(z-a[2])*v[2],0,segment.L),cx=a[0]+t*v[0],cy=a[1]+t*v[1],cz=a[2]+t*v[2],distance=Math.hypot(x-cx,y-cy,z-cz),shape=suctionRadius(cy),radius=shape?.radius??segment.radius,clearance=radius-distance;
    if(clearance>bestClear){bestClear=clearance;best={branch,segment:shape?{...segment,radius}:segment,radiusSlope:shape?shape.slope*segment.tangent[1]:0,t,cx,cy,cz,distance,clearance,s:segment.s+t};}
   }}
   best.center=[best.cx,best.cy,best.cz];best.offset=[x-best.cx,y-best.cy,z-best.cz];return best;
  }
  function clearance(x,y,dp=0,z=0){
   const rr=legacy.radius(x),main=Math.min(x,legacy.end-x,rr-Math.hypot(y,z));
   if(main>dp/2+1e-8)return main-dp/2;
   let tube=-Infinity;
   for(const b of branches){const q=closest(x,y,z,b.id),first=b.segments[0],inlet=[x,y,z].reduce((s,v,k)=>s+(v-first.a[k])*first.tangent[k],0);tube=Math.max(tube,Math.min(q.clearance,inlet));}
   return Math.min(Math.max(main,tube),legacy.end-x)-dp/2;
  }
  function release(off,branch=0){
   const b=branches[branch],distance=cad&&p.releaseMode!=='distance'?b.length-.0005:Math.min(p.release/1000,b.length-.0005),s=b.length-distance,seg=b.segments.find(v=>s<=v.s+v.L)||b.segments.at(-1),t=s-seg.s;
   // A stable transverse frame for each curved passage. Offset never becomes
   // an axial displacement on a non-planar branch.
   let e=[-seg.tangent[1],seg.tangent[0],0],norm=Math.hypot(...e);if(norm<1e-5){e=[0,-seg.tangent[2],seg.tangent[1]];norm=Math.hypot(...e);}
   return seg.a.map((v,k)=>v+t*seg.tangent[k]+off*e[k]/norm);
  }
  return {...legacy,branches,closest,clearance,release,cad,source,dimension:3};
 }
 function makeField(g,p,f,legacyGeometry){
  const candidate=p.flowPhysics==='phase2',notional=JET.notional(f,g.dn/1000),jetSpeed=candidate&&notional.valid?notional.u:f.uJet,baseGeo=geometry(g,p,legacyGeometry),geo=p.wallModel==='cad'?WALL.attach(baseGeo,g.key):baseGeo,dn=geo.dn,b0=(candidate&&notional.valid?notional.diameter:dn)/Math.sqrt(2),bleedX=g.xh!=null?(g.xh-(g.xn||0))/1000:geo.end*.5;
  const smoothJunction=candidate&&geo.branches.every(b=>Math.abs(b.segments.at(-1).tangent[0])<.1);
  const fractions=geo.branches.map((b,i)=>f.circuit?.flowUsed?f.circuit.branches[i]?.fraction||0:1/geo.branches.length);
  const entries=geo.branches.map(b=>b.entry[0]-b.radius);
  // Linked diffusivity closure: changing turbulence intensity or mixing length
  // also changes spreading. The Mach reduction concerns the shear layer only.
  const Mc=Math.abs(jetSpeed-f.ms/(f.inletP/(R*f.Tin)*area(geo.db)))/(2*Math.sqrt(1.4*R*f.Tin));
  const spread=clamp(p.spread*(g.spreadScale??1)*(p.intensity/.18)*(p.lengthScale/.12)*(p.spreadingLaw==='compressible-fit'?.23+.77*Math.exp(-3.4*Mc*Mc):1/Math.sqrt(1+.25*Mc*Mc)),.02,.35),h=dn*.035;
  const referenceRadii=Array.from({length:101},(_,i)=>geo.radius(geo.end*i/100)),confinement=Object.freeze({...JET.confinement(f,g.d/2000),constantArea:Math.max(...referenceRadii)-Math.min(...referenceRadii)<1e-6});
  let noiseScale=1;
  const ingress=x=>fractions.reduce((s,w,i)=>s+w*smooth((x-entries[i]+geo.branches[i].radius)/(2*geo.branches[i].radius)),0);
  const axialCache=new Map();
  function axial(x){
   x=clamp(x,0,geo.end);if(axialCache.has(x))return axialCache.get(x);const r=geo.radius(x),b=b0+spread*x,F=ingress(x),H=smooth((x-bleedX+dn)/(2*dn)),mass=f.mp+f.ms*F+f.mh*H,powder=f.md*F;
   const totalT=(f.mp*f.T0+f.ms*F*f.Tin+f.mh*H*293.15)/mass,split=clamp(bleedX,dn,geo.end*.8);
   const middle=f.network?.bleedPressure||f.throatP,pressure=x<split?f.portP+(middle-f.portP)*x/split:middle+(f.po-middle)*(x-split)/(geo.end-split);
   const ec=(1+powder/mass)/(2*cp)*(mass*R/(pressure*Math.PI*r*r))**2,temp=2*totalT/(1+Math.sqrt(1+4*ec*totalT));
   const rho=pressure/(R*temp),mean=mass/(rho*Math.PI*r*r),fb=b*b/(r*r)*(-Math.expm1(-r*r/(b*b))),amp0=Math.min(jetSpeed*b0/b,mean/Math.max(fb,1e-9)*.98);
   const recirc=candidate?JET.recirculation({x,r,b,rho,mean,fb,amp:amp0,F,H},f,notional.valid?notional.diameter:dn,spread,confinement):{amp:amp0,reverse:0,Ct:null,demand:0},amp=recirc.amp;
   // Analytic Gaussian moments; normalization uses gas flux, not area alone.
   const moment=k=>b*b/(k*r*r)*(-Math.expm1(-k*r*r/(b*b))),e2=moment(2),e3=moment(3),variance=e2-fb*fb;
   const beta=1+(amp/mean)**2*variance,alpha=1+3*(amp/mean)**2*variance+(amp/mean)**3*(e3-3*fb*e2+2*fb**3);
   const result=Object.freeze({recirculation:Object.freeze(recirc),x,r,b,mass,powder,temp,pressure,rho,mean,fb,amp,turbulenceAmp:amp0,beta,alpha,fluxGaussian:fb+amp/mean*variance,F,H});
   if(axialCache.size>=2048)axialCache.delete(axialCache.keys().next().value);axialCache.set(x,result);return result;
  }
  const junctionLimit=Math.max(0,...geo.branches.map(b=>b.entry[0]+b.radius)),junctionScale=[],junctionPowderScale=[];
  function sectionScale(values,x){if(!smoothJunction||!values.length||x>=junctionLimit)return 1;const u=clamp(x/junctionLimit,0,1)*(values.length-1),i=Math.min(values.length-2,Math.floor(u));return values[i]+(values[i+1]-values[i])*(u-i);}
  const mainScale=x=>sectionScale(junctionScale,x);
  const integral=(a,y)=>a.rho*(a.mean*y*y/2+a.amp*(a.b*a.b/2*(-Math.expm1(-y*y/(a.b*a.b)))-a.fb*y*y/2));
  function sample(x,y,z=0,branchIndex=null){
   const a=axial(x),rad=Math.hypot(y,z),eta=rad/a.b,gauss=Math.exp(-eta*eta),ux0=a.mean+a.amp*(gauss-a.fb),am=axial(x-h),ap=axial(x+h),dx=Math.max(ap.x-am.x,1e-12);
   const ur=rad>1e-8?-(integral(ap,rad)-integral(am,rad))/(dx*a.rho*rad):0;
   const insideMain=x>=0&&x<=geo.end&&rad<a.r;
   const q=smoothJunction?geo.closest(x,y,z,branchIndex):insideMain?{branch:geo.branches[0],segment:geo.branches[0].segments[0],s:0,clearance:0}:geo.closest(x,y,z,branchIndex);
   // Keep suction-directed transport through the main-bore overlap. The old
   // hard insideMain switch turned particles sideways at the top of the bore,
   // before they could reach the jet. Join smoothly at the branch centreline end.
   const w=smoothJunction?smooth(q.clearance/(q.segment.radius*.12))*smooth((q.branch.length-q.s)/q.segment.radius):insideMain?0:clamp(q.clearance/(q.segment.radius*.12),0,1);
   const branchFlow=f.ms*fractions[q.branch.id],common=q.segment===q.branch.segments[0]&&geo.source.common.lengthMM>0,branchMass=common?f.ms:branchFlow;
   const inletPressure=f.inletP+(f.portP-f.inletP)*q.s/q.branch.length,branchRho=inletPressure/(R*f.Tin),us=branchMass/(branchRho*area(q.segment.radius*2));
   const contraction=candidate?us*(q.radiusSlope||0)/q.segment.radius:0,branchU=q.segment.tangent.map((v,i)=>us*v+contraction*(q.offset?.[i]||0));
   const ux=candidate?(1-w)*ux0*mainScale(x)+w*branchU[0]:(1-w)*ux0+w*us*q.segment.tangent[0],uy=candidate?(1-w)*(rad?ur*y/rad:0)+w*branchU[1]:(1-w)*(rad?ur*y/rad:0)+w*us*q.segment.tangent[1],uz=candidate?(1-w)*(rad?ur*z/rad:0)+w*branchU[2]:(1-w)*(rad?ur*z/rad:0)+w*us*q.segment.tangent[2];
   const development=Math.sqrt(1-Math.exp(-((Math.max(0,x)/((g.developmentD??.7)*dn))**2))),sigmaMain=(g.turbulenceScale??1)*p.intensity*a.turbulenceAmp*development*Math.exp(-(((eta-.7)/.65)**2));
   const sigma=noiseScale*(candidate?Math.hypot(Math.sqrt(1-w)*sigmaMain,Math.sqrt(w)*.15*p.intensity*us):(1-w)*sigmaMain+w*.15*p.intensity*us),ell=Math.max(dn*.015,(1-w)*p.lengthScale*Math.min(a.b,a.r)+w*.1*q.segment.radius*2),k=1.5*sigma*sigma,eps=.9*k**1.5/ell;
   const pressure=(1-w)*a.pressure+w*inletPressure,temp=(1-w)*a.temp+w*f.Tin,nu=mu(temp)*R*temp/pressure;
   const meanEntry=entries.reduce((s,v,i)=>s+v*fractions[i],0),contrast=.9*Math.exp(-Math.max(0,x-meanEntry)/(Math.max(1e-5,geo.db*2)));
   const theta=Math.atan2(z,y),azimuth=fractions.reduce((s,v,i)=>s+v*Math.cos(theta-geo.branches[i].angle),0);
   const scalarShape=a.recirculation.reverse>0?1:(1+contrast*(a.fluxGaussian-gauss))*(1+.7*contrast*azimuth);
   const mainConcentration=a.rho*a.powder/a.mass*scalarShape,branchConcentration=branchRho*f.md/f.ms,baseConcentration=(1-w)*mainConcentration+w*branchConcentration,lambda=sectionScale(junctionPowderScale,x)-1,concentration=candidate?baseConcentration*Math.exp(lambda*(1-w)):(1-w)*a.rho*a.powder/a.mass*scalarShape+w*branchRho*f.md/f.ms;
   return {ux,uy,uz,baseConcentration,mainConcentration,ux0,mainPowderDelta:mainConcentration*((1-w)*ux-mainScale(x)*ux0),branchPowderFlux:w*branchConcentration*ux,inletWeight:w,mainMassFlux:pressure/(R*temp)*(1-w)*ux0,branchMassFlux:pressure/(R*temp)*w*branchU[0],k,eps,G:Math.sqrt(eps/nu),meanGradient:Math.abs(2*a.amp*rad/a.b**2*gauss),sigma,ell,nu,p:pressure,T:temp,rho:pressure/(R*temp),powderConcentration:concentration,bulkPowderConcentration:a.rho*a.powder/a.mass,scalarShape,b:a.b,r:a.r,TL:NUM.lagrangianTime(sigma,eps)};
  }
  // Normalize the blended axial mass flux section by section. The branch
  // direction is retained; only the main-flow component is scaled. Uniform
  // area quadrature and linear interpolation are checked independently below.
  if(smoothJunction&&junctionLimit>0){
   const ns=Math.max(24,Math.ceil(junctionLimit/(dn*.04))),scales=[];
   for(let ix=0;ix<=ns;ix++){const x=junctionLimit*ix/ns,a=axial(x),nr=96,nt=128;let main=0,branch=0;
    for(let ir=0;ir<nr;ir++)for(let it=0;it<nt;it++){const r=a.r*Math.sqrt((ir+.5)/nr),theta=2*Math.PI*(it+.5)/nt,v=sample(x,r*Math.cos(theta),r*Math.sin(theta));main+=v.mainMassFlux-a.rho*(a.mean+a.amp*(Math.exp(-((r/a.b)**2))-a.fb));branch+=v.branchMassFlux;}
    const cell=Math.PI*a.r*a.r/(nr*nt),scale=(a.mass-branch*cell)/(a.mass+main*cell);if(!(scale>0&&Number.isFinite(scale)))throw Error('Junction flux normalization has no forward-flow solution.');scales.push(scale);
   }
   junctionScale.push(...scales);
   if(f.md>0){const cs=[];for(let ix=0;ix<=ns;ix++){
    const x=junctionLimit*ix/ns,a=axial(x),nr=96,nt=128,cell=Math.PI*a.r*a.r/(nr*nt),terms=[];if(!(a.powder>0)){cs.push(1);continue;}
    let baseline=0;for(let ir=0;ir<nr;ir++)for(let it=0;it<nt;it++){const r=a.r*Math.sqrt((ir+.5)/nr),theta=2*Math.PI*(it+.5)/nt,v=sample(x,r*Math.cos(theta),r*Math.sin(theta));terms.push([v.baseConcentration*v.ux*cell,1-v.inletWeight]);baseline+=v.mainConcentration*v.ux0*mainScale(x)*cell;}
    // Positive exponential scalar normalization preserves the undiluted branch
    // (w=1). This avoids a negative main-stream concentration at high loading.
    let lambda=0;const correction=a.powder*mainScale(x)-baseline;
    for(let j=0;j<16;j++){let flux=correction*Math.exp(lambda),derivative=flux;for(const [v,q] of terms){const f=v*Math.exp(lambda*q);flux+=f;derivative+=q*f;}if(Math.abs(flux/a.powder-1)<1e-10)break;lambda-=clamp((flux-a.powder)/derivative,-2,2);}
    if(!Number.isFinite(lambda))throw Error('Junction scalar normalization is non-finite.');cs.push(1+lambda);
   }junctionPowderScale.push(...cs);}


  }
  let lossPower=0;const nx=60,nr=24,dx=geo.end/nx;
  for(let i=0;i<nx;i++){const x=(i+.5)*dx,r=geo.radius(x),dr=r/nr;for(let j=0;j<nr;j++){const y=(j+.5)*dr,v=sample(x,-y);lossPower+=v.rho*v.eps*2*Math.PI*y*dr*dx;}}
  const out=axial(geo.end),kinetic=.5*(f.mg+f.md)*out.mean**2*out.alpha,available=Math.max(0,(candidate&&notional.valid?f.power-.5*f.mp*f.uJet**2+notional.power:f.power)-kinetic);
  if(lossPower>available&&lossPower>0)noiseScale=(available/lossPower)**(1/3);
  let peak=0,peakX=0,peakY=0,peakNu=mu(f.T)*R*f.T/f.po;
  for(let i=0;i<=180;i++){const x=dn*1.001+(geo.end-dn*1.001)*i/180,r=geo.radius(x);for(let j=0;j<=48;j++){const y=-r+2*r*j/48,v=sample(x,y);if(v.G>peak){peak=v.G;peakX=x;peakY=y;peakNu=v.nu;}}}
  return {smoothJunction,junctionPowderScale,junctionScale,junctionLimit,candidate,notional,confinement,geo,sample,axial,peak,peakX,peakY,peakNu,noiseScale,lossPower:lossPower*noiseScale**3,available,kinetic,spread,convectiveMach:Mc,branchFractions:fractions,modelVersion:2,scalarClosure:'Flux-normalized radial/azimuthal mixing reconstruction; passive powder, no inferred deposition.'};
 }
 return {geometry,makeField};
})();
if(typeof module!=='undefined')module.exports=EductorSpatial;

'use strict';
// Geometric encounter kernels, SI units. No sticking, growth or depletion here.
// Fuchs transition kernel: WRF/MOSAIC module_mosaic_coag.F, lines 894–930.
// WZW relative velocity: Ayala, Rosa & Wang (2008), equations 22–24;
// theta = tau / (0.4 u'^2 / epsilon), following Wang et al. (2000).
