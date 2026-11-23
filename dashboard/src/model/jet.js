const EductorJet=(()=>{
 const R=287.05,cp=1004.675;
 function notional(flow,dn){
  const thrust=flow.network?.motiveThrust??flow.mp*flow.uJet,u=flow.mp?thrust/flow.mp:0,T=flow.T0-u*u/(2*cp),rho=flow.portP/(R*T),A=u>0?flow.mp/(rho*u):0;
  return {u,T,rho,A,diameter:Math.sqrt(4*A/Math.PI),thrust,power:.5*flow.mp*u*u,machDisk:dn*.67*Math.sqrt(flow.p0/flow.portP),valid:T>0&&Number.isFinite(A)&&A>0,source:'Mass + momentum + adiabatic energy equivalent source; Birch-type momentum formulation (1987). No shock-cell resolution.'};
 }
 function annulusMass(amp,a){
  const base=a.mean-amp*a.fb;if(base>=0)return 0;
  const y0=Math.min(1,-a.b*a.b*Math.log(-base/amp)/(a.r*a.r));
  return -a.rho*Math.PI*a.r*a.r*(base*(1-y0)+amp*a.b*a.b/(a.r*a.r)*(Math.exp(-y0*a.r*a.r/(a.b*a.b))-Math.exp(-a.r*a.r/(a.b*a.b))));
 }
 // One gas-only inlet diagnostic per passage: never recompute onset from
 // the expanding local area. Powder affects the solved flow, not this gas Ct.
 function confinement(f,r){
  const rho=f.portP/(R*f.Tin),rhoA=rho*Math.PI*r*r,uk=f.mg/rhoA;
  const momentum=(f.network?.motiveThrust??f.mp*f.uJet)+f.ms*(f.network?.secondaryVelocity||0)+f.mh*(f.network?.bleedVelocity||0),ud2=momentum/rhoA;
  return {Ct:ud2>.5*uk*uk?uk/Math.sqrt(ud2-.5*uk*uk):null,rho,uk,ud2,radius:r,convention:'Gas-only, fixed reference bore; all supplied gas mass and axial momentum. Candidate threshold, not universal.'};
 }
 function recirculation(a,f,dEq,spread,reference){
  const fixed=reference||confinement(f,a.r),Ct=fixed.Ct;
  const rhoA=a.rho*Math.PI*a.r*a.r,ud2=((f.network?.motiveThrust??f.mp*f.uJet)+f.ms*(f.network?.secondaryVelocity||0)+f.mh*(f.network?.bleedVelocity||0))/rhoA;
  // The cylindrical confined-jet closure has no validated diffuser/junction
  // extension. Keep its diagnostic, but do not impose reverse flow there.
  const supported=fixed.constantArea!==false;
  const entrained=.32*f.mp*a.x/Math.max(dEq,1e-9),supply=f.ms*a.F+f.mh*a.H,closure=Math.max(0,Math.min(1,(1-a.b/a.r)/.2)),demand=supported&&Ct!==null&&Ct<.98?Math.max(0,entrained-supply)*closure:0;
  if(!demand)return {amp:a.amp,Ct,supported,reason:supported?null:'Variable-area passage: diffuser recirculation not modeled',demand:0,reverse:annulusMass(a.amp,a),momentumLimited:false};
  let lo=a.mean/Math.max(a.fb,1e-12),hi=Math.max(lo*2,a.amp);while(annulusMass(hi,a)<demand&&hi<1e6)hi*=2;
  for(let k=0;k<45;k++){const mid=(lo+hi)/2;if(annulusMass(mid,a)<demand)lo=mid;else hi=mid;}
  const variance=a.b*a.b/(2*a.r*a.r)*(-Math.expm1(-2*a.r*a.r/(a.b*a.b)))-a.fb*a.fb;
  const budget=Math.max(0,ud2-a.mean*a.mean),cap=Math.sqrt(budget/Math.max(variance,1e-16)),amp=Math.max(a.amp,Math.min((lo+hi)/2,cap));
  return {amp,Ct,supported,demand,reverse:annulusMass(amp,a),momentumLimited:cap<(lo+hi)/2,closureX:(a.r-dEq/Math.sqrt(2))/Math.max(spread,1e-9)};
 }
 return {notional,annulusMass,confinement,recirculation};
})();
if(typeof module!=='undefined')module.exports=EductorJet;
