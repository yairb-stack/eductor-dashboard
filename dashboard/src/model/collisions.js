const ExposureCollisions=(()=>{
 const NUM=typeof ExposureNumerics!=='undefined'?ExposureNumerics:require('./numerics.js');
 const kB=1.380649e-23, Rgas=287.05;
 const SCREEN={radiusEta:.1,turbulentMach:.3,particleRe:.1,densityRatio:100,tracerSt:.1};
 const KEYS=['brownian','turbulent','inertial','unsupported','tracer'];
 function properties(d,rhop,p,T){
  if(![d,rhop,p,T].every(v=>Number.isFinite(v)&&v>0))throw Error('Collision kernels require positive finite diameter, density, absolute pressure and temperature.');
  const {mu,rho,lambda}=NUM.gas(p,T),Kn=2*lambda/d,Cc=1+Kn*(1.257+.4*Math.exp(-1.1/Kn)),mass=rhop*Math.PI*d**3/6;
  const D=kB*T*Cc/(3*Math.PI*mu*d),c=Math.sqrt(8*kB*T/(Math.PI*mass)),ell=8*D/(Math.PI*c),q=ell/d;
  // Algebraically rationalized Fuchs length. Avoid subtracting nearly equal
  // cubes for q << 1 or q >> 1; the small-q series avoids final cancellation.
  const h=Math.hypot(1,q),delta=d*(q<1e-4?q/2+q*q/3-q**3/8:2*((1+q)**2+(1+q)*h+h*h)/(3*(1+q+h))-1);
  return {d,mass,D,c,delta,Kn,Cc,lambda,mu,rho,nu:mu/rho,tau:rhop*d*d*Cc/(18*mu)};
 }
 function fuchs(a,b){const R=(a.d+b.d)/2,D=a.D+b.D,c=Math.hypot(a.c,b.c);return 4*Math.PI*R*D/(R/(R+Math.hypot(a.delta,b.delta))+4*D/(R*c));}
 function brownian(d1,d2,rhop,p,T){return fuchs(properties(d1,rhop,p,T),properties(d2,rhop,p,T));}
 function saffmanTurner(d1,d2,G){return Math.sqrt(8*Math.PI/15)*((d1+d2)/2)**3*Math.max(0,G);}
 function accelerationVariance(tau1,tau2,sigma,eps,nu,{TL:givenTL,convention='zhou2001'}={}){
  if(!(sigma>0&&eps>0&&nu>0&&tau1>0&&tau2>0))return 0;
  const TL=givenTL??NUM.lagrangianTime(sigma,eps),a=tau1/TL,b=tau2/TL,s=a+b,ab=a*b,phi=Math.max(a/b,b/a),gamma=(convention==='restatement2025'?1:phi)*.183*sigma*sigma/Math.sqrt(eps*nu),A=(1+a)*(1+b),B=(1+gamma*a)*(1+gamma*b);
  // WZW correction to the Kruis–Kusters acceleration variance. Rationalizing
  // gamma/(gamma-1)*(1/A-1/B) removes the removable singularity at gamma=1.
  const decorrelation=((a-b)**2+4*ab*(ab/(A*(1+Math.sqrt((1+s)/A)))))/s;
  return (convention==='restatement2025'?1/3:1)*(1+.6*Math.exp(-((phi-1)**1.5)))*sigma*sigma*decorrelation*gamma*(s+(gamma+1)*ab)/(A*B);
 }
 // Zhou, Wexler & Wang (2001), p.102. DNS fit at contact R/eta=1;
 // at other contact radii this is a sensitivity upper bound, not validation.
 // Below St=.1 the empirical tanh tails are tapered to recover g→1.
 // This asymptotic guard is our extension, not a fitted Zhou coefficient.
 function rdfMono(alpha,Rlambda){const z0=.5*(1+Math.tanh((alpha-.5)/.25)),z1=.5*(1+Math.tanh((alpha-1.25)/.1)),z2=.5*(1+Math.tanh((alpha-6.5)/2.5)),y0=18*alpha*alpha,y1=.36*alpha**2.5*Math.exp(-(alpha**2.5)),y2=.24*Math.exp(-.5*alpha),y3=.013*Math.exp(-.07*alpha);const excess=y0*(1-z0*z0)+Rlambda*z0*z0*(y1*(1-z1)+y2*z1+y3*z2);return 1+excess*Math.min(1,(alpha/.1)**2);}
 function rdf(tau1,tau2,G,sigma,nu){if(!(G>0&&sigma>0))return 1;const psi=Math.max(tau1/tau2,tau2/tau1),correlation=Math.min(1,2.6*Math.exp(-psi)+.205*Math.exp(-.0206*psi)*.5*(1+Math.tanh(psi-3))),Rlambda=Math.sqrt(15)*sigma*sigma/(nu*G);return 1+correlation*Math.sqrt(Math.max(0,rdfMono(tau1*G,Rlambda)-1)*Math.max(0,rdfMono(tau2*G,Rlambda)-1));}
 function abrahamson(R,tau1,tau2,sigma,TL){return Math.sqrt(8*Math.PI)*R*R*sigma*Math.sqrt(1/(1+tau1/TL)+1/(1+tau2/TL));}
 function pair(a,d1,d2,rhop,{legacy=false,tracerOnly=false,slip=0,convention='zhou2001',extended=false,clustering=false}={}){
  const i=properties(d1,rhop,a.p,a.T),j=d1===d2?i:properties(d2,rhop,a.p,a.T),G=Math.max(0,a.G||0),nu=i.nu,eps=G*G*nu,sigma=Math.max(0,a.sigma||0),R=(d1+d2)/2,eta=G>0?Math.sqrt(nu/G):Infinity;
  const st=Math.max(i.tau,j.tau)*G,size=R/eta,Mt=Math.sqrt(3)*sigma/Math.sqrt(1.4*Rgas*a.T),TL=a.TL??NUM.lagrangianTime(sigma,eps);
  const relativeSlip=q=>sigma>0?Math.sqrt(3)*sigma*Math.sqrt(q.tau/(TL+q.tau)):0;
  // Estimated turbulent slip for both species, plus measured tracked-particle
  // slip when available. The partner slip remains a local statistical estimate.
  const maxRe=Math.max(i.d*Math.max(Math.abs(slip),relativeSlip(i)),j.d*relativeSlip(j))/nu;
  const tracer=saffmanTurner(d1,d2,G),variance=legacy||tracerOnly?0:accelerationVariance(i.tau,j.tau,sigma,eps,nu,{TL,convention}),inertialKernel=Math.sqrt(8*Math.PI)*R*R*Math.sqrt(variance);
  // Combine turbulent shear and acceleration variances, not their kernels.
  const kinetic=!legacy&&!tracerOnly&&extended&&(st>100||size>1),g12=!legacy&&!tracerOnly&&clustering&&!kinetic?rdf(i.tau,j.tau,G,sigma,nu):1,unclustered=kinetic?abrahamson(R,i.tau,j.tau,sigma,TL):Math.hypot(tracer,inertialKernel),turbulent=unclustered*g12,brownianKernel=legacy?0:fuchs(i,j),rate=brownianKernel+turbulent;
  const flags=[];if(!kinetic&&size>SCREEN.radiusEta)flags.push('finite-size');if(Mt>SCREEN.turbulentMach)flags.push('compressible-turbulence');if(!kinetic&&maxRe>SCREEN.particleRe)flags.push('non-Stokes-drag');if(rhop/i.rho<SCREEN.densityRatio)flags.push('density-ratio');if(G>0&&!(sigma>0)&&!legacy&&!tracerOnly)flags.push('missing-turbulence-variance');if((legacy||tracerOnly)&&st>SCREEN.tracerSt)flags.push('tracer-inertia');
  if(kinetic&&Math.min(i.tau,j.tau)/TL<10)flags.push('kinetic-decorrelation-unverified');if(!legacy&&!tracerOnly&&clustering&&!kinetic&&Math.abs(size-1)>.1)flags.push('RDF-contact-radius-extrapolation');const unsupported=flags.length?turbulent:0;
  return {TL,convention,g12,kinetic,unclustered,rate,brownian:brownianKernel,turbulent,inertial:Math.max(0,turbulent-(kinetic?0:tracer)),tracer:kinetic?0:tracer,screened:brownianKernel+(flags.length?0:turbulent),unsupported,maxSt:st,maxSize:size,maxRe,Mt,flags,regime:st<=SCREEN.tracerSt?'tracer':'inertial'};
 }
 return {rdfMono,rdf,abrahamson,properties,fuchs,brownian,saffmanTurner,accelerationVariance,pair,SCREEN,KEYS,label:model=>model==='legacy'?'Legacy Saffman–Turner (no Brownian)':model==='tracer-brownian'?'Fuchs Brownian + Saffman–Turner':'Fuchs Brownian + WZW turbulence',version:'fuchs-zhou2001-st-reference-v2',conventions:{zhou2001:'Zhou et al. 2001 pp.101–102: gamma includes phi; radial acceleration variance has no extra 1/3.',restatement2025:'Sensitivity only: 2025 restatement without phi and with 1/3. Not the original Zhou convention.'}};
})();
if(typeof module!=='undefined')module.exports=ExposureCollisions;

'use strict';
