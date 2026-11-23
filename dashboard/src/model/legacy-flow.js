const Legacy=(()=>{
"use strict";
/* ---------------- flow model (mirrors eductor_model.py) ---------------- */
const GAM=1.4, RG=287.05, T0=293.15, ATM=101325, MU=1.81e-5;
let PA=ATM, PF=ATM, RHO=PA/(RG*T0), RHOF=RHO, NU=MU/RHO;
function setBoundary(suctionAtm=1,outletAtm=1){
  PF=suctionAtm*ATM; PA=outletAtm*ATM; RHO=PA/(RG*T0); RHOF=PF/(RG*T0); NU=MU/RHO;
}
const CD_SHARP=0.82, CD=CD_SHARP, CD_CONE=0.94, CS=0.12, FD=0.02;   // fitted to the PSD database; textbook values 0.84 and 0.16
const DN=3.5e-3, D1=22.1e-3, D15=34.8e-3, LM1=61.7e-3, LM15=110.3e-3, LLEG=87.7e-3, LDROP=71.3e-3, DIN=60e-3;
const area=d=>Math.PI*d*d/4;
function nozzle(p0,CD=CD_SHARP,Tm=T0,dn=DN){      // Tm: stagnation temperature of the motive air, K
  const g=GAM, a=area(dn), crit=Math.pow(2/(g+1),g/(g-1)), n={p0, Tm, aEff:CD*a};
  if(p0<=PA){ return Object.assign(n,{mdot:0,S:PA*CD*a,fNet:0,uEq:0,dEq:0,mj:0,uj:0,power:0,choked:false,vExit:0}); }
  n.choked = PA/p0<=crit;
  if(n.choked){
    n.mdot=CD*a*p0*Math.sqrt(g/(RG*Tm))*Math.pow(2/(g+1),(g+1)/(2*(g-1)));
    n.S=CD*a*crit*p0*(1+g); n.vExit=Math.sqrt(2*g/(g+1)*RG*Tm);
  }else{
    const m2=2/(g-1)*(Math.pow(PA/p0,-(g-1)/g)-1), t=Tm/(1+(g-1)/2*m2), v=Math.sqrt(m2*g*RG*t);
    n.mdot=CD*a*PA/(RG*t)*v; n.S=n.mdot*v+CD*a*PA; n.vExit=v;
  }
  n.fNet=n.S-PA*n.aEff; n.uEq=n.fNet/n.mdot; n.dEq=2*n.mdot/Math.sqrt(Math.PI*RHO*n.fNet);
  const npr=p0/PA; n.mj=Math.sqrt(2/(g-1)*(Math.pow(npr,(g-1)/g)-1));
  n.uj=n.mj*Math.sqrt(g*RG*Tm/(1+(g-1)/2*n.mj*n.mj));
  n.power=n.mdot*g/(g-1)*RG*Tm*(1-Math.pow(npr,-(g-1)/g));
  return n;
}
const betaExit=(d,L)=>Math.max(1,Math.pow(d/2/(CS*L),2));
const tMix=(st,ms)=>{ const n=st.n, mg=n.mdot+ms; return mg>1e-9? (n.mdot*n.Tm+ms*(st.Ts||T0))/mg : T0; };   // mixed stream temperature, equal cp
function exitP(st,ms,md,p1){
  const a=area(st.d), ab=area(st.db), n=st.n, vb=ms/((st.rhoS||RHO)*ab), mg=n.mdot+ms, rhoM=RHO*T0/tMix(st,ms);
  const k=betaExit(st.d,st.L)+FD*st.L/(2*st.d);
  return (n.S+p1*(a-n.aEff)+(ms+md)*vb*Math.cos(st.th*Math.PI/180)-(k*mg+md)*mg/(rhoM*a))/a;
}
function suction(pt,ms,md,db,k,vup=0,rhoS=RHO){
  const ab=area(db), vb=ms/(rhoS*ab), q=.5*rhoS*vb*vb;
  return pt-Math.sign(ms)*(1+k)*q-md*(vb-vup)/ab;
}
function bisect(f,lo,hi){ let flo=f(lo), fhi=f(hi);
  if(!Number.isFinite(flo)||!Number.isFinite(fhi)) return NaN;
  if(Math.abs(flo)<1e-9) return lo; if(Math.abs(fhi)<1e-9) return hi;
  if((flo>0)===(fhi>0)) return NaN; for(let i=0;i<70;i++){ const mid=.5*(lo+hi), fm=f(mid); if((fm>0)===(flo>0)){lo=mid;flo=fm}else hi=mid; } return .5*(lo+hi); }
let KS_A=Math.pow(D1/DIN,4)+0.05*Math.pow(D1/D15,4)+0.20+FD*LLEG/D1;
let KS_1=Math.pow(D15/DIN,4)+0.05+FD*LDROP/D15;
let K_RED=0.20+FD*LLEG/D1;
function report(st,ms,md,p1,pe,ks){
  const n=st.n, a=area(st.d), mg=n.mdot+ms, Tmix=tMix(st,ms), rhoM=RHO*T0/Tmix, vExit=mg/(rhoM*a), bExit=Math.min(st.d/2,CS*st.L);
  const vJet=mg/(rhoM*Math.PI*bExit*bExit);
  const lEnt=n.mdot>0? mg/(CS*Math.sqrt(Math.PI*RHO*n.fNet)) : Infinity;
  const eps=Math.max(n.power-.5*mg*vJet*vJet,0)/(rhoM*a*st.L), shear=Math.sqrt(eps/NU), tRes=st.L/vExit;
  return {md,n,st,mp:n.mdot,ms,mg,omega:ms/n.mdot,p1:p1-PA,pe:pe-PA,vb:ms/((st.rhoS||RHO)*area(st.db)),vExit,vJet,beta:betaExit(st.d,st.L),
    lAttach:st.d/2/CS,lEnt,loading:md/mg,conc:md/(mg/rhoM),eps,shear,tRes,ks,Tmix};
}
function valveState(dpv,v,shut){          // equivalent loss coefficient and an indicative ball angle for the chosen pressure drop
  if(!(dpv>0)) return {K:0,ang:0,shut:false};
  if(shut||!(v>0.05)) return {K:Infinity,ang:90,shut:true};
  const K=dpv/(.5*RHO*v*v), T=[[0.05,0],[0.3,10],[1.6,20],[5.5,30],[17,40],[53,50],[206,60],[750,70],[1e4,80],[1e7,88]];
  let ang=88; if(K<=T[0][0]) ang=0; else for(let i=1;i<T.length;i++) if(K<=T[i][0]){ const t=Math.log(K/T[i-1][0])/Math.log(T[i][0]/T[i-1][0]); ang=T[i-1][1]+t*(T[i][1]-T[i-1][1]); break; }
  return {K,ang,shut:false};
}
function solveA(pg,feed,ext,dpv=0,Tm=T0){
  const md=feed/3600, n=nozzle(pg+PA,CD_SHARP,Tm), st={n,d:D1,L:LM1+ext,th:45,db:D1,rhoS:RHOF};
  const p1of=ms=>suction(PF-dpv,ms,md,D1,KS_A,0,RHOF);
  let ms=bisect(m=>exitP(st,m,md,p1of(m))-PA,-0.02,0.2), p1=p1of(ms), shut=false;
  if(dpv>0&&ms<=1e-6){ shut=true; ms=0; p1=bisect(p=>exitP(st,0,0,p)-PA,PA-5e4,PA+5e4); }   // the valve holds more than the jet can pull: no flow
  const r=report(st,ms,shut?0:md,p1,PA,KS_A); r.valve=valveState(dpv,r.vb,shut); return r;
}
function solveB(pg1,pg2,feed,ext,dpv=0,Tm=T0,dpf=0){      // dpv: valve between the stages; dpf: valve in the feeder tube above the tee
  const n1=nozzle(pg1+PA,CD_SHARP,Tm), n2=nozzle(pg2+PA,CD_SHARP,Tm);
  const s1={n:n1,d:D15,L:LM15,th:90,db:D15,rhoS:RHOF}, s2={n:n2,d:D1,L:LM1+ext,th:45,db:D1};
  const chain=(ms1,md,p11o)=>{ const p11=p11o!=null? p11o : suction(PF-Math.sign(ms1)*dpf,ms1,md,D15,KS_1,0,RHOF), p2=exitP(s1,ms1,md,p11), m1=n1.mdot+ms1, T1=tMix(s1,ms1), rho1=RHO*T0/T1, v2=m1/(rho1*area(D15));
    const pt2=p2+betaExit(D15,LM15)*.5*rho1*v2*Math.abs(v2), p12=suction(pt2-Math.sign(m1)*dpv,m1,md,D1,K_RED,v2,rho1); return {p11,p2,pt2,p12,m1,s2:Object.assign({},s2,{Ts:T1,rhoS:rho1})}; };
  const res=(c,md)=>exitP(c.s2,c.m1,md,c.p12)-PA;
  let md=feed/3600, ms1=bisect(m=>res(chain(m,md),md),-0.02,0.2), shut=false, c;
  if(dpf>0&&Math.abs(ms1)<1e-6){ shut=true; ms1=0; md=0; c=chain(0,0,bisect(p=>res(chain(0,0,p),0),Math.max(1,PF-dpf),Math.max(PA,PF)+dpf+5e4)); }   // the feeder valve holds: no flow, the pressure under it floats
  else c=chain(ms1,md);
  const r1=report(s1,ms1,md,c.p11,c.p2,KS_1), r2=report(c.s2,c.m1,md,c.p12,PA,K_RED);
  r1.pt=c.pt2-PA; r1.valve=valveState(dpf,Math.abs(r1.vb),shut); r2.valve=valveState(dpv,r2.vb,false); return [r1,r2];
}
/* British eductor: 16 mm port, 8.2 mm converge-diverge insert, 8 x 8 mm bleed holes, 13.7 mm throat, 10 degree diffuser */
const BR={dPort:16e-3,rT:4.09e-3,xT:14.3e-3,rC:5.2e-3,xC:21.9e-3,aH:8*Math.PI*16e-6,d2:13.7e-3,x2in:32.3e-3,x2out:56.3e-3,dE:19.34e-3,KH:1.5,KP:1.3,K1:0.59,KE2:0.15,eta:0.5};
function solveC(pg,feed,dpv=0,Tm=T0){
  const md=feed/3600, n=nozzle(pg+PA,CD_CONE,Tm), sq=v=>v*Math.abs(v);
  const jetR=x=>n.mdot>0? Math.max(n.dEq/2,CS*x):0, ent=x=>n.mdot>0? n.mdot*Math.max(0,2*CS*x/n.dEq-1):0;
  const aP=area(BR.dPort), aC=Math.PI*BR.rC*BR.rC-Math.PI*Math.pow(Math.min(jetR(BR.xC),BR.rC),2), a2=area(BR.d2), aA2=a2-Math.PI*Math.pow(Math.min(jetR(BR.x2in),BR.d2/2),2), aE=area(BR.dE);
  const beta=betaExit(BR.d2,BR.x2out), k=beta+FD*(BR.x2out-BR.x2in)/(2*BR.d2);
  const pGap=mh=>PA-BR.KH*.5*RHO*sq(mh/(RHO*BR.aH));
  const feeder=pg_=>{ const e1=ent(BR.xC); return bisect(ms=>{ const vp=ms/(RHOF*aP), vc=(ms-e1)/(RHO*aC); return PF-dpv-BR.KP*.5*RHOF*sq(vp)-md*vp/aP-pg_-(1+BR.K1)*.5*RHO*sq(vc); },-0.02,0.1); };
  const chain=mh=>{ const pg_=pGap(mh); let ms=feeder(pg_), shut=false; if(dpv>0&&ms<=1e-6){ ms=0; shut=true; }
    const co=ms+mh-ent(BR.x2in), va=co/(RHO*aA2), p2s=pg_-(1+BR.KE2)*.5*RHO*sq(va), mg=n.mdot+ms+mh, mdd=shut?0:md;
    const Tmix=(n.mdot*n.Tm+(ms+mh)*T0)/mg, rhoM=RHO*T0/Tmix;
    const p2m=(n.S+p2s*(a2-n.aEff)+co*va-(k*mg+mdd)*mg/(rhoM*a2))/a2, v2=mg/(rhoM*a2), pe=p2m+BR.eta*(1-Math.pow(a2/aE,2))*.5*rhoM*v2*v2;
    return {pg:pg_,ms,shut,p2s,p2m,pe,mg,v2,mdd,Tmix,rhoM}; };
  const mh=bisect(m=>chain(m).pe-PA,-0.02,0.2), c=chain(mh), vp=c.ms/(RHOF*aP);
  const r={md:c.mdd,n,mp:n.mdot,ms:c.ms,mh,mg:c.mg,pCh:PF-PA-dpv-(1+BR.KP)*.5*RHOF*sq(vp),pGap:c.pg-PA,pThroat:c.p2m-PA,vPort:vp,vThroat:c.v2,vExit:c.mg/(c.rhoM*aE),Tmix:c.Tmix,beta,
    lAttach:BR.d2/2/CS,loading:c.mdd/c.mg,loadFeed:c.ms>0? c.mdd/c.ms:Infinity,conc:c.mdd/(c.mg/c.rhoM),vAnn1:(c.ms-ent(BR.xC))/(RHO*aC),vHole:mh/(RHO*BR.aH)};
  r.valve=valveState(dpv,vp*aP/area(D1),c.shut); return r;
}
// sensitivity: scale every assumed loss coefficient by m, run fn, restore
const K0={a:KS_A,b:KS_1,r:K_RED,KH:BR.KH,KP:BR.KP,K1:BR.K1,KE2:BR.KE2};
function withLoss(m,fn){
  KS_A=K0.a*m; KS_1=K0.b*m; K_RED=K0.r*m; BR.KH=K0.KH*m; BR.KP=K0.KP*m; BR.K1=K0.K1*m; BR.KE2=K0.KE2*m;
  try{ return fn(); } finally { KS_A=K0.a; KS_1=K0.b; K_RED=K0.r; BR.KH=K0.KH; BR.KP=K0.KP; BR.K1=K0.K1; BR.KE2=K0.KE2; }
}
function exposure(q,r){ if(!(r.mp>0)) return 0; const xs=r.n.uEq*r.n.dEq/(2*CS)*Math.sqrt(RHO/(2*q)); return Math.min(1,xs/r.lEnt); }
const collisions=(r,rhoC=1000)=>7.8/Math.PI*(r.conc/rhoC)*r.shear*r.tRes;


function assess(r,kind){
  const pressures=kind==='C'?[r.pCh,r.pGap,r.pThroat]:[r.p1,r.pe];
  let reason='';
  if(![r.ms,r.mg,...pressures].every(Number.isFinite)) reason='No bracketed flow solution';
  else if(r.ms<=1e-7 || !(r.mg>0)) reason='No forward powder transport';
  else if(pressures.some(v=>v+PA<=0)) reason='Nonphysical pressure: compressible solution needed';
  else if(kind==='C' && r.mh<0) reason='Bleed outflow: powder split is unresolved';
  else if(kind==='C' && r.n.dEq/2>=BR.rT) reason='Expanded jet fills the insert: compressible solution needed';
  const bulkMach=Math.max(Math.abs(r.vb||r.vPort||0)/Math.sqrt(GAM*RG*T0),Math.abs(r.vExit||0)/Math.sqrt(GAM*RG*(r.Tmix||T0)),kind==='C'?Math.abs(r.vThroat)/Math.sqrt(GAM*RG*r.Tmix):0);
  if(!reason&&bulkMach>=1) reason='Sonic bulk flow: compressible solution needed';
  const chokeCap=area(kind==='C'?BR.dPort:r.st.db)*PF*Math.sqrt(GAM/(RG*T0))*Math.pow(2/(GAM+1),(GAM+1)/(2*(GAM-1)));
  if(!reason&&kind!=='B2'&&r.ms>chokeCap) reason='Feeder choking: compressible solution needed';
  r.supported=!reason; r.reason=reason; r.bulkMach=bulkMach;
  r.caution=bulkMach>0.3||pressures.some(v=>Math.abs(v)/PA>0.1)||Math.abs(PF/PA-1)>0.1;
  return r;
}
function prepareResults(A,B,C){
  assess(A,'A'); assess(B[0],'B1'); assess(B[1],'B2'); assess(C,'C');
  if(!B[0].supported||!B[1].supported){const why=B.find(r=>!r.supported).reason;B.forEach(r=>{r.supported=false;r.reason=why;});}
  return [A,B[0],B[1],C];
}
function safeDrawing(r){
  if(r.supported)return r;
  const copy={...r,n:{...r.n},valve:{...r.valve,shut:true}};
  for(const k of Object.keys(copy))if(typeof copy[k]==='number')copy[k]=0;
  for(const k of Object.keys(copy.n))if(typeof copy.n[k]==='number')copy.n[k]=0;
  copy.Tmix=T0; copy.beta=1;copy.n.Tm=T0;copy.n.p0=PA;
  return copy;
}
function stationRows(A,B1,B2,C){
  const row=(name,air,powder,pAbs,T,detail)=>({name,air,powder,pAbs,T,ratio:air>0?powder/air:null,fraction:air+powder>0?powder/(air+powder):null,detail});
  const station=(r,name,air,pAbs,T,detail)=>row(name,air,r.md,pAbs,T,detail);
  return [
    {name:'1-stage',r:A,rows:[row('Motive nozzle',A.mp,0,A.n.p0,A.n.Tm,'Pure motive air; supply pressure'),station(A,'Feeder / powder port',A.ms,PA+A.p1,T0,'Before motive air joins'),station(A,'After the powder junction',A.mg,PA+A.p1,A.Tmix,'All incoming streams; section-average ratio'),station(A,'Mixing tube outlet',A.mg,PA,A.Tmix,'Same ratio: no further air inlet')]},
    {name:'2-stage',r:B1,rows:[station(B1,'Feeder / first tee',B1.ms,PA+B1.p1,T0,'Before the first motive jet'),station(B1,'1.5″ mixing tube outlet',B1.mg,PA+B1.pe,B1.Tmix,'Feeder air + first motive jet'),station(B1,'Interstage / second powder port',B2.ms,PA+B2.p1,B1.Tmix,'The valve changes flow; it adds no air'),station(B2,'After second motive jet',B2.mg,PA+B2.p1,B2.Tmix,'Both motive jets + feeder air'),station(B2,'Final outlet',B2.mg,PA,B2.Tmix,'Same ratio as after the second jet')]},
    {name:'British',r:C,rows:[row('Motive nozzle',C.mp,0,C.n.p0,C.n.Tm,'Pure motive air; supply pressure'),station(C,'Powder port / suction chamber',C.ms,PA+C.pCh,T0,'Powder-bearing feed stream only'),station(C,'Before bleed-air entry',C.mp+C.ms,PA+C.pGap,(C.mp*C.n.Tm+C.ms*T0)/(C.mp+C.ms),'Section includes motive air; mixing may be incomplete'),station(C,'After the bleed holes',C.mg,PA+C.pGap,C.Tmix,'Includes all inward bleed air'),station(C,'13.7 mm mixing throat',C.mg,PA+C.pThroat,C.Tmix,'No additional air source'),station(C,'Diffuser outlet',C.mg,PA,C.Tmix,'Same mass ratio despite area expansion')]}
  ];
}


function solveCase(g,p){
 setBoundary(p.suction,p.outlet);
 const mm=x=>x/1000, md=p.feed/3600, tm=p.temp+273.15, pg=p.pg*1e5;
 let r;
 if(g.family==='british'){
  const backup={...BR};
  Object.assign(BR,{dPort:mm(g.db),rT:mm(g.rt),rC:mm(g.rc),xC:mm(g.xc-g.xn),aH:g.nh*area(mm(g.dh)),d2:mm(g.d),x2in:mm(g.xin-g.xn),x2out:mm(g.xe-g.xn),dE:mm(g.de),KH:1.5*p.loss,KP:1.3*p.loss,K1:(g.insert?.59:.05)*p.loss,KE2:.15*p.loss,eta:g.insert?.5:0});
  // solveC uses this per-case nozzle through a temporary adapter; no case shares state after return.
  r=solveBritish(pg,p.feed,tm,g.dn/1000);
  Object.assign(r,{st:{d:mm(g.d),db:mm(g.db),L:mm(g.L)},family:g.family});
  assess(r,'C');
  Object.assign(BR,backup);
 } else {
  const n=nozzle(pg+PA,CD_SHARP,tm,mm(g.dn)),a=area(mm(g.d)),af=area(mm(g.df)),ap=g.nports*area(mm(g.db)),ao=area(mm(g.de));
  const st={n,d:mm(g.d),db:mm(g.db)*Math.sqrt(g.nports),L:mm(g.L),th:45,rhoS:RHOF};
  // A common feed supplies the manifold: port flows are equal only by the symmetry assumption.
  const ps=s=>PF-.5*RHOF*(p.loss*2*(s/(RHOF*af))**2+(1+.3*p.loss)*(s/(RHOF*ap))**2)-md*s/(RHOF*ap*ap);
  const outletDrop=s=>{const mg=n.mdot+s,rho=RHO*T0/tMix(st,s),vm=mg/(rho*a),vo=mg/(rho*ao);return .5*rho*(vo*vo-vm*vm)+(g.de<g.d?.2*p.loss*.5*rho*vo*vo:0)+md*(vo-vm)/ao;};
  const balance=s=>exitP(st,s,md,ps(s))-outletDrop(s)-PA;
  const s=bisect(balance,0,.2);
  r=report(st,s,md,ps(s),PA+outletDrop(s),2*p.loss);r.vMix=r.vExit;r.vExit=r.mg/(RHO*T0/r.Tmix*ao);r.pMix=r.pe;r.pe=0;r.family=g.family;
  assess(r,'A');
  // Diagnose the existing balance without inventing a carrier-flow fallback.
  r.balance={zeroSecondaryResidual:balance(0),residual:Number.isFinite(s)?balance(s):null,
    zeroSecondaryVesselPressure:PF-balance(0)*a/(a-n.aEff)};
  if(!Number.isFinite(s)&&r.balance.zeroSecondaryResidual<0){r.reason='No forward suction: outlet pressure exceeds the reduced-model pressure recovery';}
  const cap=af*PF*Math.sqrt(GAM/(RG*T0))*Math.pow(2/(GAM+1),(GAM+1)/(2*(GAM-1)));
  if(r.ms>cap){r.supported=false;r.reason='Shared feeder reaches choking; compressible model needed';}
  if(r.supported&&r.n.dEq>mm(g.de)){r.supported=false;r.reason='Expanded jet exceeds the exit bore; compressible model needed';}
 }
 r.geometry=g;r.md=md;r.loading=md/r.mg;r.feedLoading=md/r.ms;r.portAir=r.ms/g.nports;r.portPowder=md/g.nports;
 r.rows=[{name:'Motive nozzle',air:r.mp,powder:0,pressure:r.n.p0,detail:'Pure air; pressure is stagnation supply pressure'},
 {name:g.nports>1?'Common feeder before split':'Powder feeder',air:r.ms,powder:md,pressure:PF,detail:'Upstream suction vessel'},
 {name:g.nports>1?'Each suction port ('+g.nports+' equal branches)':'Suction port',air:r.portAir,powder:r.portPowder,pressure:PA+(r.pCh??r.p1),detail:g.nports>1?'Equal gas and powder split assumed':'Powder-bearing stream only'},
 {name:'After motive air joins',air:r.mp+r.ms,powder:md,pressure:PA+(r.pCh??r.p1),detail:'Section-average loading; local mixing may be incomplete'}];
 if(g.family==='british')r.rows.push({name:'After 8 bleed holes',air:r.mg,powder:md,pressure:PA+r.pGap,detail:'Bleed air comes from the selected outlet surroundings'}, {name:'Mixing throat / bore',air:r.mg,powder:md,pressure:PA+r.pThroat,detail:'Mass ratio remains constant downstream of the last inlet'});
 else r.rows.push({name:'End of mixing bore',air:r.mg,powder:md,pressure:PA+r.pMix,detail:g.de<g.d?'Before the outlet contraction':'No further air inlet'});
 r.rows.push({name:'Outlet',air:r.mg,powder:md,pressure:PA,detail:'Conserved powder and air mass flow'});
 for(const row of r.rows){row.ratio=row.air>0?row.powder/row.air:null;row.fraction=row.air+row.powder>0?row.powder/(row.air+row.powder):null;}
 return r;
}

function solveBritish(pg,feed,Tm,dn){
  const md=feed/3600, dpv=0, n=nozzle(pg+PA,CD_CONE,Tm,dn), sq=v=>v*Math.abs(v);
  const jetR=x=>n.mdot>0? Math.max(n.dEq/2,CS*x):0, ent=x=>n.mdot>0? n.mdot*Math.max(0,2*CS*x/n.dEq-1):0;
  const aP=area(BR.dPort), aC=Math.PI*BR.rC*BR.rC-Math.PI*Math.pow(Math.min(jetR(BR.xC),BR.rC),2), a2=area(BR.d2), aA2=a2-Math.PI*Math.pow(Math.min(jetR(BR.x2in),BR.d2/2),2), aE=area(BR.dE);
  const beta=betaExit(BR.d2,BR.x2out), k=beta+FD*(BR.x2out-BR.x2in)/(2*BR.d2);
  const pGap=mh=>PA-BR.KH*.5*RHO*sq(mh/(RHO*BR.aH));
  const feeder=pg_=>{ const e1=ent(BR.xC); return bisect(ms=>{ const vp=ms/(RHOF*aP), vc=(ms-e1)/(RHO*aC); return PF-dpv-BR.KP*.5*RHOF*sq(vp)-md*vp/aP-pg_-(1+BR.K1)*.5*RHO*sq(vc); },-0.02,0.1); };
  const chain=mh=>{ const pg_=pGap(mh); let ms=feeder(pg_), shut=false; if(dpv>0&&ms<=1e-6){ ms=0; shut=true; }
    const co=ms+mh-ent(BR.x2in), va=co/(RHO*aA2), p2s=pg_-(1+BR.KE2)*.5*RHO*sq(va), mg=n.mdot+ms+mh, mdd=shut?0:md;
    const Tmix=(n.mdot*n.Tm+(ms+mh)*T0)/mg, rhoM=RHO*T0/Tmix;
    const p2m=(n.S+p2s*(a2-n.aEff)+co*va-(k*mg+mdd)*mg/(rhoM*a2))/a2, v2=mg/(rhoM*a2), pe=p2m+BR.eta*(1-Math.pow(a2/aE,2))*.5*rhoM*v2*v2;
    return {pg:pg_,ms,shut,p2s,p2m,pe,mg,v2,mdd,Tmix,rhoM}; };
  const mh=bisect(m=>chain(m).pe-PA,-0.02,0.2), c=chain(mh), vp=c.ms/(RHOF*aP);
  const r={md:c.mdd,n,mp:n.mdot,ms:c.ms,mh,mg:c.mg,pCh:PF-PA-dpv-(1+BR.KP)*.5*RHOF*sq(vp),pGap:c.pg-PA,pThroat:c.p2m-PA,vPort:vp,vThroat:c.v2,vExit:c.mg/(c.rhoM*aE),Tmix:c.Tmix,beta,
    lAttach:BR.d2/2/CS,loading:c.mdd/c.mg,loadFeed:c.ms>0? c.mdd/c.ms:Infinity,conc:c.mdd/(c.mg/c.rhoM),vAnn1:(c.ms-ent(BR.xC))/(RHO*aC),vHole:mh/(RHO*BR.aH)};
  r.valve=valveState(dpv,vp*aP/area(D1),c.shut); return r;
}

if(typeof module!=="undefined")module.exports={solveCase};

// UI entry point: shared boundaries plus an automatically selected compressible estimate.
function solveCaseEstimate(g,p){
 const net=typeof EductorNetwork!=='undefined'?EductorNetwork:require('../exposure/network.js');
 p={...p,outlet:p.suction};
 const base=solveCase(g,p);if(p.modelVersion===2){const u=typeof EductorUpgrades!=='undefined'?EductorUpgrades:require('../upgrades/physics.js'),geom={...g,key:'nov'+g.id,angle:g.family==='british'?90:45},po=p.suction*101325,f=u.solve(geom,{ps:po/1e5,p0:po/1e5+p.pg,temp:p.temp,feed:p.feed,loss:p.loss});return net.legacyResult(f,geom,base,po);}
 if(!net.assess(base,g,{ps:p.suction*1.01325}).needed)return base;
 const geom={...g,key:'nov'+g.id,angle:g.family==='british'?90:45},po=p.suction*101325;
 const flow=net.solve(geom,{ps:po/1e5,p0:po/1e5+p.pg,temp:p.temp,feed:p.feed,loss:p.loss});
 return net.legacyResult(flow,geom,base,po);
}

if(typeof module!=='undefined')module.exports.solveCaseEstimate=solveCaseEstimate;
;return {setBoundary,solveCase,solveA,solveB,solveC,prepareResults,stationRows};})();const GEOMETRIES=[{"key":"revA","name":"Rev 2 \u00b7 1-stage","dn":3.5,"db":22.1,"d":22.1,"de":22.1,"L":61.7,"end":61.7,"angle":45,"nports":1,"chamber":22.1,"family":"branched","insert":false,"profile":[[0,11.05],[61.7,11.05]],"source":"Existing Rev 2 dashboard nominal dimensions."},{"key":"revB1","name":"Rev 2 \u00b7 2-stage: first stage","dn":3.5,"db":34.8,"d":34.8,"de":34.8,"L":110.3,"end":110.3,"angle":90,"nports":1,"chamber":34.8,"family":"branched","insert":false,"profile":[[0,17.4],[110.3,17.4]],"source":"Existing Rev 2 dashboard nominal dimensions. Stage-local release; interstage particle-history transfer is not modeled."},{"key":"revB2","name":"Rev 2 \u00b7 2-stage: second stage","dn":3.5,"db":22.1,"d":22.1,"de":22.1,"L":61.7,"end":61.7,"angle":45,"nports":1,"chamber":22.1,"family":"branched","insert":false,"profile":[[0,11.05],[61.7,11.05]],"source":"Existing Rev 2 dashboard nominal dimensions. Stage-local release; interstage particle-history transfer is not modeled."},{"key":"revC","name":"British \u00b7 slide geometry","family":"british","dn":3.5,"db":16,"d":13.7,"de":19.3424,"L":56.3,"end":72.2871884,"angle":90,"nports":1,"chamber":22.5,"rt":4.09,"rc":5.2,"xc":21.9,"xn":0,"xh":28,"nh":8,"dh":8,"insert":true,"profile":[[0.0,11.0],[0.3012,11.0],[0.60239,11.0],[0.90359,11.0],[1.20479,11.0],[1.50598,11.0],[1.80718,11.0],[2.10838,11.0],[2.40957,11.0],[2.71077,11.0],[3.01197,11.0],[3.31316,11.0],[3.61436,11.0],[3.91556,11.0],[4.21675,11.0],[4.51795,11.0],[4.81915,11.0],[5.12034,11.0],[5.42154,11.0],[5.72274,11.0],[6.02393,11.0],[6.32513,11.0],[6.62633,11.0],[6.92752,11.0],[7.22872,11.0],[7.52992,11.0],[7.83111,11.0],[8.13231,9.04256],[8.43351,8.0895],[8.7347,7.44895],[9.0359,6.94733],[9.3371,6.53105],[9.63829,6.17577],[9.93949,5.86708],[10.24069,5.59632],[10.54188,5.35733],[10.84308,5.14612],[11.14427,4.95915],[11.44547,4.79437],[11.74667,4.64965],[12.04786,4.52362],[12.34906,4.41512],[12.65026,4.32324],[12.95145,4.24722],[13.25265,4.18648],[13.55385,4.14056],[13.85504,4.10914],[14.15624,4.0922],[14.45744,4.08942],[14.75863,4.10077],[15.05983,4.12635],[15.36103,4.16633],[15.66222,4.21404],[15.96342,4.26177],[16.26462,4.30948],[16.56581,4.35718],[16.86701,4.40489],[17.16821,4.45259],[17.4694,4.5003],[17.7706,4.548],[18.0718,4.59571],[18.37299,4.64341],[18.67419,4.69112],[18.97539,4.73882],[19.27658,4.78653],[19.57778,4.83423],[19.87898,4.88194],[20.18017,4.92964],[20.48137,4.97735],[20.78257,5.02505],[21.08376,5.07275],[21.38496,5.12046],[21.68616,5.16811],[21.98735,11.25],[22.28855,11.25],[22.58975,11.25],[22.89094,11.25],[23.19214,11.25],[23.49334,11.25],[23.79453,11.25],[24.09573,11.25],[24.39693,11.25],[24.69812,11.25],[24.99932,11.25],[25.30052,11.25],[25.60171,11.25],[25.90291,11.25],[26.20411,11.25],[26.5053,11.25],[26.8065,11.25],[27.1077,11.25],[27.40889,11.25],[27.71009,9.5401],[28.01129,9.36317],[28.31248,9.18614],[28.61368,9.00911],[28.91488,8.83208],[29.21607,8.65505],[29.51727,8.47802],[29.81847,8.30099],[30.11966,8.12396],[30.42086,7.94693],[30.72206,7.7699],[31.02325,7.59287],[31.32445,7.41584],[31.62564,7.23881],[31.92684,7.06175],[32.22804,6.88463],[32.52923,6.84987],[32.83043,6.84987],[33.13163,6.84987],[33.43282,6.84987],[33.73402,6.84987],[34.03522,6.84987],[34.33641,6.84987],[34.63761,6.84987],[34.93881,6.84987],[35.24,6.84987],[35.5412,6.84987],[35.8424,6.84987],[36.14359,6.84987],[36.44479,6.84987],[36.74599,6.84987],[37.04718,6.84987],[37.34838,6.84987],[37.64958,6.84987],[37.95077,6.84987],[38.25197,6.84987],[38.55317,6.84987],[38.85436,6.84987],[39.15556,6.84987],[39.45676,6.84987],[39.75795,6.84987],[40.05915,6.84987],[40.36035,6.84987],[40.66154,6.84987],[40.96274,6.84987],[41.26394,6.84987],[41.56513,6.84987],[41.86633,6.84987],[42.16753,6.84987],[42.46872,6.84987],[42.76992,6.84987],[43.07112,6.84987],[43.37231,6.84987],[43.67351,6.84987],[43.97471,6.84987],[44.2759,6.84987],[44.5771,6.84987],[44.8783,6.84987],[45.17949,6.84987],[45.48069,6.84987],[45.78189,6.84987],[46.08308,6.84987],[46.38428,6.84987],[46.68548,6.84987],[46.98667,6.84987],[47.28787,6.84987],[47.58907,6.84987],[47.89026,6.84987],[48.19146,6.84987],[48.49266,6.84987],[48.79385,6.84987],[49.09505,6.84987],[49.39625,6.84987],[49.69744,6.84987],[49.99864,6.84987],[50.29984,6.84987],[50.60103,6.84987],[50.90223,6.84987],[51.20343,6.84987],[51.50462,6.84987],[51.80582,6.84987],[52.10701,6.84987],[52.40821,6.84987],[52.70941,6.84987],[53.0106,6.84987],[53.3118,6.84987],[53.613,6.84987],[53.91419,6.84987],[54.21539,6.84987],[54.51659,6.84987],[54.81778,6.84987],[55.11898,6.84987],[55.42018,6.84987],[55.72137,6.84987],[56.02257,6.84987],[56.32377,6.85633],[56.62496,6.90948],[56.92616,6.96264],[57.22736,7.01578],[57.52855,7.06889],[57.82975,7.122],[58.13095,7.17511],[58.43214,7.22822],[58.73334,7.28133],[59.03454,7.33444],[59.33573,7.38755],[59.63693,7.44066],[59.93813,7.49376],[60.23932,7.54687],[60.54052,7.59998],[60.84172,7.65309],[61.14291,7.7062],[61.44411,7.75931],[61.74531,7.81242],[62.0465,7.86553],[62.3477,7.91864],[62.6489,7.97175],[62.95009,8.02486],[63.25129,8.07796],[63.55249,8.13107],[63.85368,8.18418],[64.15488,8.23729],[64.45608,8.2904],[64.75727,8.34351],[65.05847,8.39662],[65.35967,8.44973],[65.66086,8.50284],[65.96206,8.55595],[66.26326,8.60906],[66.56445,8.66216],[66.86565,8.71527],[67.16685,8.76838],[67.46804,8.82149],[67.76924,8.8746],[68.07044,8.92771],[68.37163,8.98082],[68.67283,9.03393],[68.97403,9.08704],[69.27522,9.14015],[69.57642,9.19326],[69.87762,9.24636],[70.17881,9.29947],[70.48001,9.35258],[70.78121,9.40569],[71.0824,9.4588],[71.3836,9.51191],[71.6848,9.56497],[71.98599,9.61801],[72.28719,9.61801]],"source":"Slide geometry centre-plane section, simplified to a closed radial envelope.","suctionProfile":[[12.0,7.999835649441105],[13.000000000000002,7.999837419901852],[14.0,7.999819081093637],[15.0,7.999870621159801],[16.0,7.999876937558089],[17.0,7.999888765502138],[18.000000000000004,7.999914724648536],[19.0,7.999952100919642],[20.0,8.0],[21.0,8.011936994380152],[22.0,8.047082723602852],[23.0,8.105495293375293],[24.0,8.187273812890652],[25.0,8.292559074745249],[26.000000000000004,8.421534529758784],[27.0,8.574427573333363],[28.0,8.751672598241823],[29.0,8.942492348217858],[30.000000000000004,9.133903453202446],[31.0,9.32539879462792],[32.0,9.516982216799251],[33.0,9.708657569309377],[34.0,9.900428705827615],[35.0,10.091912846983096],[36.00000000000001,10.283330994637657],[37.00000000000001,10.474844935296025],[38.00000000000001,10.666458515890715],[39.0,10.85817558662884],[40.0,11.05]],"suctionProfileSource":"OpenFOAM reference closed CAD envelope: circular sections normal to y, 1 mm spacing, 12\u201340 mm upstream; no fit to particle outcomes."},{"name":"British 1\u2033 \u00b7 straight bore","family":"british","dn":3.5,"db":16,"d":22,"de":22,"nports":1,"xn":30.3983,"xe":114,"xh":66.4,"dh":8,"nh":8,"rt":11,"rc":11,"xc":62.4,"xin":70.4,"insert":false,"id":0,"file":"British Eductor 1'' - w-o converge-diverge.step","hash":"b9c8c52bc3705774c7eacd6cc4ef87c254e0533a7784cb32dbcc6c931af8f95a","L":83.6017,"key":"nov0","chamber":22,"angle":90,"end":83.6017,"profile":[[0.0,11.0],[0.34834,11.0],[0.69668,11.0],[1.04502,11.0],[1.39336,11.0],[1.7417,11.0],[2.09004,11.0],[2.43838,11.0],[2.78672,11.0],[3.13506,11.0],[3.4834,11.0],[3.83174,11.0],[4.18009,11.0],[4.52843,11.0],[4.87677,11.0],[5.22511,11.0],[5.57345,11.0],[5.92179,11.0],[6.27013,11.0],[6.61847,11.0],[6.96681,11.0],[7.31515,11.0],[7.66349,11.0],[8.01183,11.0],[8.36017,11.0],[8.70851,11.0],[9.05685,11.0],[9.40519,11.0],[9.75353,11.0],[10.10187,11.0],[10.45021,11.0],[10.79855,11.0],[11.14689,11.0],[11.49523,11.0],[11.84357,11.0],[12.19191,11.0],[12.54025,11.0],[12.8886,11.0],[13.23694,11.0],[13.58528,11.0],[13.93362,11.0],[14.28196,11.0],[14.6303,11.0],[14.97864,11.0],[15.32698,11.0],[15.67532,11.0],[16.02366,11.0],[16.372,11.0],[16.72034,11.0],[17.06868,11.0],[17.41702,11.0],[17.76536,11.0],[18.1137,11.0],[18.46204,11.0],[18.81038,11.0],[19.15872,11.0],[19.50706,11.0],[19.8554,11.0],[20.20374,11.0],[20.55208,11.0],[20.90042,11.0],[21.24877,11.0],[21.59711,11.0],[21.94545,11.0],[22.29379,11.0],[22.64213,11.0],[22.99047,11.0],[23.33881,11.0],[23.68715,11.0],[24.03549,11.0],[24.38383,11.0],[24.73217,11.0],[25.08051,11.0],[25.42885,11.0],[25.77719,11.0],[26.12553,11.0],[26.47387,11.0],[26.82221,11.0],[27.17055,11.0],[27.51889,11.0],[27.86723,11.0],[28.21557,11.0],[28.56391,11.0],[28.91225,11.0],[29.26059,11.0],[29.60894,11.0],[29.95728,11.0],[30.30562,11.0],[30.65396,11.0],[31.0023,11.0],[31.35064,11.0],[31.69898,11.0],[32.04732,11.0],[32.39566,11.0],[32.744,11.0],[33.09234,11.0],[33.44068,11.0],[33.78902,11.0],[34.13736,11.0],[34.4857,11.0],[34.83404,11.0],[35.18238,11.0],[35.53072,11.0],[35.87906,11.0],[36.2274,11.0],[36.57574,11.0],[36.92408,11.0],[37.27242,11.0],[37.62076,11.0],[37.96911,11.0],[38.31745,11.0],[38.66579,11.0],[39.01413,11.0],[39.36247,11.0],[39.71081,11.0],[40.05915,11.0],[40.40749,11.0],[40.75583,11.0],[41.10417,11.0],[41.45251,11.0],[41.80085,11.0],[42.14919,11.0],[42.49753,11.0],[42.84587,11.0],[43.19421,11.0],[43.54255,11.0],[43.89089,11.0],[44.23923,11.0],[44.58757,11.0],[44.93591,11.0],[45.28425,11.0],[45.63259,11.0],[45.98094,11.0],[46.32928,11.0],[46.67762,11.0],[47.02596,11.0],[47.3743,11.0],[47.72264,11.0],[48.07098,11.0],[48.41932,11.0],[48.76766,11.0],[49.116,11.0],[49.46434,11.0],[49.81268,11.0],[50.16102,11.0],[50.50936,11.0],[50.8577,11.0],[51.20604,11.0],[51.55438,11.0],[51.90272,11.0],[52.25106,11.0],[52.5994,11.0],[52.94774,11.0],[53.29608,11.0],[53.64442,11.0],[53.99276,11.0],[54.3411,11.0],[54.68945,11.0],[55.03779,11.0],[55.38613,11.0],[55.73447,11.0],[56.08281,11.0],[56.43115,11.0],[56.77949,11.0],[57.12783,11.0],[57.47617,11.0],[57.82451,11.0],[58.17285,11.0],[58.52119,11.0],[58.86953,11.0],[59.21787,11.0],[59.56621,11.0],[59.91455,11.0],[60.26289,11.0],[60.61123,11.0],[60.95957,11.0],[61.30791,11.0],[61.65625,11.0],[62.00459,11.0],[62.35293,11.0],[62.70127,11.0],[63.04962,11.0],[63.39796,11.0],[63.7463,11.0],[64.09464,11.0],[64.44298,11.0],[64.79132,11.0],[65.13966,11.0],[65.488,11.0],[65.83634,11.0],[66.18468,11.0],[66.53302,11.0],[66.88136,11.0],[67.2297,11.0],[67.57804,11.0],[67.92638,11.0],[68.27472,11.0],[68.62306,11.0],[68.9714,11.0],[69.31974,11.0],[69.66808,11.0],[70.01642,11.0],[70.36476,11.0],[70.7131,11.0],[71.06144,11.0],[71.40979,11.0],[71.75813,11.0],[72.10647,11.0],[72.45481,11.0],[72.80315,11.0],[73.15149,11.0],[73.49983,11.0],[73.84817,11.0],[74.19651,11.0],[74.54485,11.0],[74.89319,11.0],[75.24153,11.0],[75.58987,11.0],[75.93821,11.0],[76.28655,11.0],[76.63489,11.0],[76.98323,11.0],[77.33157,11.0],[77.67991,11.0],[78.02825,11.0],[78.37659,11.0],[78.72493,11.0],[79.07327,11.0],[79.42161,11.0],[79.76996,11.0],[80.1183,11.0],[80.46664,11.0],[80.81498,11.0],[81.16332,11.0],[81.51166,11.0],[81.86,11.0],[82.20834,11.0],[82.55668,11.0],[82.90502,11.0],[83.25336,11.0],[83.6017,11.0]],"source":"20Nov STEP dimensions; axisymmetric envelope closes side holes. One representative suction branch."},{"name":"British 1\u2033 \u00b7 C\u2013D","family":"british","dn":3.5,"db":16,"d":13.7,"de":19.3424,"nports":1,"xn":30.3983,"xe":98,"xh":66.4,"dh":8,"nh":8,"rt":4.089,"rc":5.2,"xc":63.6,"xin":74,"insert":true,"id":1,"file":"British Eductor 1'' - with converge-diverge.step","hash":"37591a07f7a1e3917abab0396dc8b8592e9b8c2540d445da05dd5a66dd0b497c","L":67.6017,"key":"nov1","chamber":22,"angle":90,"end":83.6017,"profile":[[0.0,11.0],[0.34834,11.0],[0.69668,11.0],[1.04502,11.0],[1.39336,11.0],[1.7417,11.0],[2.09004,11.0],[2.43838,11.0],[2.78672,11.0],[3.13506,11.0],[3.4834,11.0],[3.83174,11.0],[4.18009,11.0],[4.52843,11.0],[4.87677,11.0],[5.22511,11.0],[5.57345,11.0],[5.92179,11.0],[6.27013,11.0],[6.61847,11.0],[6.96681,11.0],[7.31515,11.0],[7.66349,11.0],[8.01183,11.0],[8.36017,11.0],[8.70851,11.0],[9.05685,11.0],[9.40519,11.0],[9.75353,11.0],[10.10187,11.0],[10.45021,11.0],[10.79855,11.0],[11.14689,11.0],[11.49523,11.0],[11.84357,11.0],[12.19191,11.0],[12.54025,11.0],[12.8886,11.0],[13.23694,11.0],[13.58528,11.0],[13.93362,11.0],[14.28196,11.0],[14.6303,11.0],[14.97864,11.0],[15.32698,11.0],[15.67532,11.0],[16.02366,11.0],[16.372,11.0],[16.72034,11.0],[17.06868,11.0],[17.41702,11.0],[17.76536,11.0],[18.1137,11.0],[18.46204,11.0],[18.81038,11.0],[19.15872,11.0],[19.50706,8.82764],[19.8554,7.85375],[20.20374,7.20803],[20.55208,6.69416],[20.90042,6.24522],[21.24877,5.90005],[21.59711,5.57746],[21.94545,5.31142],[22.29379,5.08184],[22.64213,4.85902],[22.99047,4.70424],[23.33881,4.54947],[23.68715,4.41614],[24.03549,4.32605],[24.38383,4.23597],[24.73217,4.1701],[25.08051,4.1388],[25.42885,4.1075],[25.77719,4.09961],[26.12553,4.12553],[26.47387,4.15146],[26.82221,4.18912],[27.17055,4.24432],[27.51889,4.29953],[27.86723,4.35474],[28.21557,4.40995],[28.56391,4.46516],[28.91225,4.52036],[29.26059,4.5755],[29.60894,4.63078],[29.95728,4.68599],[30.30562,4.7412],[30.65396,4.7964],[31.0023,4.85161],[31.35064,4.90672],[31.69898,4.96185],[32.04732,11.0],[32.39566,11.0],[32.744,11.0],[33.09234,11.0],[33.44068,11.0],[33.78902,11.0],[34.13736,11.0],[34.4857,11.0],[34.83404,11.0],[35.18238,11.0],[35.53072,11.0],[35.87906,11.0],[36.2274,11.0],[36.57574,11.0],[36.92408,11.0],[37.27242,11.0],[37.62076,11.0],[37.96911,11.0],[38.31745,11.0],[38.66579,11.0],[39.01413,11.0],[39.36247,11.0],[39.71081,11.0],[40.05915,8.93221],[40.40749,8.72747],[40.75583,8.52272],[41.10417,8.31797],[41.45251,8.11323],[41.80085,7.90848],[42.14919,7.70373],[42.49753,7.49899],[42.84587,7.29424],[43.19421,7.08949],[43.54255,6.88475],[43.89089,6.85],[44.23923,6.85],[44.58757,6.85],[44.93591,6.85],[45.28425,6.85],[45.63259,6.85],[45.98094,6.85],[46.32928,6.85],[46.67762,6.85],[47.02596,6.85],[47.3743,6.85],[47.72264,6.85],[48.07098,6.85],[48.41932,6.85],[48.76766,6.85],[49.116,6.85],[49.46434,6.85],[49.81268,6.85],[50.16102,6.85],[50.50936,6.85],[50.8577,6.85],[51.20604,6.85],[51.55438,6.85],[51.90272,6.85],[52.25106,6.85],[52.5994,6.85],[52.94774,6.85],[53.29608,6.85],[53.64442,6.85],[53.99276,6.85],[54.3411,6.85],[54.68945,6.85],[55.03779,6.85],[55.38613,6.85],[55.73447,6.85],[56.08281,6.85],[56.43115,6.85],[56.77949,6.85],[57.12783,6.85],[57.47617,6.85],[57.82451,6.85],[58.17285,6.85],[58.52119,6.85],[58.86953,6.85],[59.21787,6.85],[59.56621,6.85],[59.91455,6.85],[60.26289,6.85],[60.61123,6.85],[60.95957,6.85],[61.30791,6.85],[61.65625,6.85],[62.00459,6.85],[62.35293,6.85],[62.70127,6.85],[63.04962,6.85],[63.39796,6.85],[63.7463,6.85],[64.09464,6.85],[64.44298,6.85],[64.79132,6.85],[65.13966,6.85],[65.488,6.85],[65.83634,6.85],[66.18468,6.85],[66.53302,6.85],[66.88136,6.85],[67.2297,6.85],[67.57804,6.85],[67.92638,6.90713],[68.27472,6.96842],[68.62306,7.02971],[68.9714,7.091],[69.31974,7.1525],[69.66808,7.21407],[70.01642,7.27564],[70.36476,7.33717],[70.7131,7.39859],[71.06144,7.45991],[71.40979,7.5212],[71.75813,7.58248],[72.10647,7.64389],[72.45481,7.70546],[72.80315,7.76703],[73.15149,7.8285],[73.49983,7.88992],[73.84817,7.95134],[74.19651,8.01268],[74.54485,8.07397],[74.89319,8.13527],[75.24153,8.19684],[75.58987,8.25842],[75.93821,8.31984],[76.28655,8.38125],[76.63489,8.44267],[76.98323,8.50409],[77.33157,8.56546],[77.67991,8.62675],[78.02825,8.68823],[78.37659,8.74975],[78.72493,8.81117],[79.07327,8.87259],[79.42161,8.934],[79.76996,8.99542],[80.1183,9.05684],[80.46664,9.11823],[80.81498,9.17962],[81.16332,9.24109],[81.51166,9.30251],[81.86,9.36392],[82.20834,9.42534],[82.55668,9.48676],[82.90502,9.54817],[83.25336,9.60959],[83.6017,9.60959]],"source":"20Nov STEP dimensions; axisymmetric envelope closes side holes. One representative suction branch."},{"name":"British original \u00b7 C\u2013D","family":"british","dn":2,"db":10,"d":8.56,"de":12.0866,"nports":1,"xn":24.144,"xe":61.2,"xh":41.45,"dh":5,"nh":8,"rt":2.552,"rc":3.25,"xc":39.7,"xin":46.2,"insert":true,"id":2,"file":"British Eductor - original size - with converge-diverge.step","hash":"f4c5f9e89805e8536577b542324c9ef3095e631ca88cc4650f3c3275c808db77","L":37.056,"key":"nov2","chamber":14,"angle":90,"end":47.056000000000004,"profile":[[0.0,7.0],[0.19607,7.0],[0.39213,7.0],[0.5882,7.0],[0.78427,7.0],[0.98033,7.0],[1.1764,7.0],[1.37247,7.0],[1.56853,7.0],[1.7646,7.0],[1.96067,7.0],[2.15673,7.0],[2.3528,7.0],[2.54887,7.0],[2.74493,7.0],[2.941,7.0],[3.13707,7.0],[3.33313,7.0],[3.5292,7.0],[3.72527,7.0],[3.92133,7.0],[4.1174,7.0],[4.31347,7.0],[4.50953,7.0],[4.7056,7.0],[4.90167,7.0],[5.09773,7.0],[5.2938,7.0],[5.48987,7.0],[5.68593,7.0],[5.882,7.0],[6.07807,7.0],[6.27413,7.0],[6.4702,7.0],[6.66627,7.0],[6.86233,6.14957],[7.0584,5.25985],[7.25447,4.79015],[7.45053,4.41711],[7.6466,4.15372],[7.84267,3.89369],[8.03873,3.70615],[8.2348,3.51861],[8.43087,3.36291],[8.62693,3.23019],[8.823,3.09748],[9.01907,2.99729],[9.21513,2.90764],[9.4112,2.81799],[9.60727,2.75259],[9.80333,2.7],[9.9994,2.64741],[10.19547,2.60735],[10.39153,2.58909],[10.5876,2.57084],[10.78367,2.55259],[10.97973,2.56594],[11.1758,2.58034],[11.37187,2.59475],[11.56793,2.61852],[11.764,2.64948],[11.96007,2.68033],[12.15613,2.71123],[12.3522,2.74219],[12.54827,2.77364],[12.74433,2.80473],[12.9404,2.83576],[13.13647,2.86678],[13.33253,2.89767],[13.5286,2.92862],[13.72467,2.95988],[13.92073,2.9909],[14.1168,3.02193],[14.31287,3.05296],[14.50893,3.08399],[14.705,3.11502],[14.90107,7.0],[15.09713,7.0],[15.2932,7.0],[15.48927,7.0],[15.68533,7.0],[15.8814,7.0],[16.07747,7.0],[16.27353,7.0],[16.4696,7.0],[16.66567,7.0],[16.86173,7.0],[17.0578,7.0],[17.25387,7.0],[17.44993,7.0],[17.646,7.0],[17.84207,7.0],[18.03813,7.0],[18.2342,7.0],[18.43027,7.0],[18.62633,7.0],[18.8224,7.0],[19.01847,7.0],[19.21453,7.0],[19.4106,7.0],[19.60667,7.0],[19.80273,7.0],[19.9988,5.48871],[20.19487,5.37351],[20.39093,5.25831],[20.587,5.14311],[20.78307,5.02791],[20.97913,4.91271],[21.1752,4.79751],[21.37127,4.68231],[21.56733,4.56711],[21.7634,4.45191],[21.95947,4.33671],[22.15553,4.28],[22.3516,4.28],[22.54767,4.28],[22.74373,4.28],[22.9398,4.28],[23.13587,4.28],[23.33193,4.28],[23.528,4.28],[23.72407,4.28],[23.92013,4.28],[24.1162,4.28],[24.31227,4.28],[24.50833,4.28],[24.7044,4.28],[24.90047,4.28],[25.09653,4.28],[25.2926,4.28],[25.48867,4.28],[25.68473,4.28],[25.8808,4.28],[26.07687,4.28],[26.27293,4.28],[26.469,4.28],[26.66507,4.28],[26.86113,4.28],[27.0572,4.28],[27.25327,4.28],[27.44933,4.28],[27.6454,4.28],[27.84147,4.28],[28.03753,4.28],[28.2336,4.28],[28.42967,4.28],[28.62573,4.28],[28.8218,4.28],[29.01787,4.28],[29.21393,4.28],[29.41,4.28],[29.60607,4.28],[29.80213,4.28],[29.9982,4.28],[30.19427,4.28],[30.39033,4.28],[30.5864,4.28],[30.78247,4.28],[30.97853,4.28],[31.1746,4.28],[31.37067,4.28],[31.56673,4.28],[31.7628,4.28],[31.95887,4.28],[32.15493,4.28],[32.351,4.28],[32.54707,4.28],[32.74313,4.28],[32.9392,4.28],[33.13527,4.28],[33.33133,4.28],[33.5274,4.28],[33.72347,4.28],[33.91953,4.28],[34.1156,4.28],[34.31167,4.28],[34.50773,4.28],[34.7038,4.28],[34.89987,4.28],[35.09593,4.28],[35.292,4.28],[35.48807,4.28],[35.68413,4.28],[35.8802,4.28],[36.07627,4.28],[36.27233,4.28],[36.4684,4.28],[36.66447,4.28],[36.86053,4.28],[37.0566,4.28011],[37.25267,4.31462],[37.44873,4.34913],[37.6448,4.38364],[37.84087,4.41815],[38.03693,4.45274],[38.233,4.48747],[38.42907,4.52208],[38.62513,4.55664],[38.8212,4.59121],[39.01727,4.62577],[39.21333,4.66034],[39.4094,4.69491],[39.60547,4.72947],[39.80153,4.76404],[39.9976,4.79861],[40.19367,4.83317],[40.38973,4.86774],[40.5858,4.90231],[40.78187,4.93682],[40.97793,4.97133],[41.174,5.00584],[41.37007,5.04036],[41.56613,5.07487],[41.7622,5.10952],[41.95827,5.14421],[42.15433,5.17884],[42.3504,5.21341],[42.54647,5.24797],[42.74253,5.28254],[42.9386,5.31711],[43.13467,5.35167],[43.33073,5.38624],[43.5268,5.4208],[43.72287,5.45531],[43.91893,5.48982],[44.115,5.52434],[44.31107,5.55885],[44.50713,5.59355],[44.7032,5.6282],[44.89927,5.66277],[45.09533,5.69734],[45.2914,5.7319],[45.48747,5.76647],[45.68353,5.80104],[45.8796,5.8356],[46.07567,5.87017],[46.27173,5.90474],[46.4678,5.9393],[46.66387,5.97387],[46.85993,6.00844],[47.056,6.00844]],"source":"20Nov STEP dimensions; axisymmetric envelope closes side holes. One representative suction branch."},{"name":"British original \u00b7 straight bore","family":"british","dn":2,"db":10,"d":14,"de":14,"nports":1,"xn":24.144,"xe":71.2,"xh":41.45,"dh":5,"nh":8,"rt":7,"rc":7,"xc":38.95,"xin":43.95,"insert":false,"id":3,"file":"British Eductor - original size w-o converge-diverge.step","hash":"4a72e84eec0c32029692e3d7335a547c0d571fd3f80a0fcca6915538b0fa36fe","L":47.056,"key":"nov3","chamber":14,"angle":90,"end":47.056000000000004,"profile":[[0.0,7.0],[0.19607,7.0],[0.39213,7.0],[0.5882,7.0],[0.78427,7.0],[0.98033,7.0],[1.1764,7.0],[1.37247,7.0],[1.56853,7.0],[1.7646,7.0],[1.96067,7.0],[2.15673,7.0],[2.3528,7.0],[2.54887,7.0],[2.74493,7.0],[2.941,7.0],[3.13707,7.0],[3.33313,7.0],[3.5292,7.0],[3.72527,7.0],[3.92133,7.0],[4.1174,7.0],[4.31347,7.0],[4.50953,7.0],[4.7056,7.0],[4.90167,7.0],[5.09773,7.0],[5.2938,7.0],[5.48987,7.0],[5.68593,7.0],[5.882,7.0],[6.07807,7.0],[6.27413,7.0],[6.4702,7.0],[6.66627,7.0],[6.86233,7.0],[7.0584,7.0],[7.25447,7.0],[7.45053,7.0],[7.6466,7.0],[7.84267,7.0],[8.03873,7.0],[8.2348,7.0],[8.43087,7.0],[8.62693,7.0],[8.823,7.0],[9.01907,7.0],[9.21513,7.0],[9.4112,7.0],[9.60727,7.0],[9.80333,7.0],[9.9994,7.0],[10.19547,7.0],[10.39153,7.0],[10.5876,7.0],[10.78367,7.0],[10.97973,7.0],[11.1758,7.0],[11.37187,7.0],[11.56793,7.0],[11.764,7.0],[11.96007,7.0],[12.15613,7.0],[12.3522,7.0],[12.54827,7.0],[12.74433,7.0],[12.9404,7.0],[13.13647,7.0],[13.33253,7.0],[13.5286,7.0],[13.72467,7.0],[13.92073,7.0],[14.1168,7.0],[14.31287,7.0],[14.50893,7.0],[14.705,7.0],[14.90107,7.0],[15.09713,7.0],[15.2932,7.0],[15.48927,7.0],[15.68533,7.0],[15.8814,7.0],[16.07747,7.0],[16.27353,7.0],[16.4696,7.0],[16.66567,7.0],[16.86173,7.0],[17.0578,7.0],[17.25387,7.0],[17.44993,7.0],[17.646,7.0],[17.84207,7.0],[18.03813,7.0],[18.2342,7.0],[18.43027,7.0],[18.62633,7.0],[18.8224,7.0],[19.01847,7.0],[19.21453,7.0],[19.4106,7.0],[19.60667,7.0],[19.80273,7.0],[19.9988,7.0],[20.19487,7.0],[20.39093,7.0],[20.587,7.0],[20.78307,7.0],[20.97913,7.0],[21.1752,7.0],[21.37127,7.0],[21.56733,7.0],[21.7634,7.0],[21.95947,7.0],[22.15553,7.0],[22.3516,7.0],[22.54767,7.0],[22.74373,7.0],[22.9398,7.0],[23.13587,7.0],[23.33193,7.0],[23.528,7.0],[23.72407,7.0],[23.92013,7.0],[24.1162,7.0],[24.31227,7.0],[24.50833,7.0],[24.7044,7.0],[24.90047,7.0],[25.09653,7.0],[25.2926,7.0],[25.48867,7.0],[25.68473,7.0],[25.8808,7.0],[26.07687,7.0],[26.27293,7.0],[26.469,7.0],[26.66507,7.0],[26.86113,7.0],[27.0572,7.0],[27.25327,7.0],[27.44933,7.0],[27.6454,7.0],[27.84147,7.0],[28.03753,7.0],[28.2336,7.0],[28.42967,7.0],[28.62573,7.0],[28.8218,7.0],[29.01787,7.0],[29.21393,7.0],[29.41,7.0],[29.60607,7.0],[29.80213,7.0],[29.9982,7.0],[30.19427,7.0],[30.39033,7.0],[30.5864,7.0],[30.78247,7.0],[30.97853,7.0],[31.1746,7.0],[31.37067,7.0],[31.56673,7.0],[31.7628,7.0],[31.95887,7.0],[32.15493,7.0],[32.351,7.0],[32.54707,7.0],[32.74313,7.0],[32.9392,7.0],[33.13527,7.0],[33.33133,7.0],[33.5274,7.0],[33.72347,7.0],[33.91953,7.0],[34.1156,7.0],[34.31167,7.0],[34.50773,7.0],[34.7038,7.0],[34.89987,7.0],[35.09593,7.0],[35.292,7.0],[35.48807,7.0],[35.68413,7.0],[35.8802,7.0],[36.07627,7.0],[36.27233,7.0],[36.4684,7.0],[36.66447,7.0],[36.86053,7.0],[37.0566,7.0],[37.25267,7.0],[37.44873,7.0],[37.6448,7.0],[37.84087,7.0],[38.03693,7.0],[38.233,7.0],[38.42907,7.0],[38.62513,7.0],[38.8212,7.0],[39.01727,7.0],[39.21333,7.0],[39.4094,7.0],[39.60547,7.0],[39.80153,7.0],[39.9976,7.0],[40.19367,7.0],[40.38973,7.0],[40.5858,7.0],[40.78187,7.0],[40.97793,7.0],[41.174,7.0],[41.37007,7.0],[41.56613,7.0],[41.7622,7.0],[41.95827,7.0],[42.15433,7.0],[42.3504,7.0],[42.54647,7.0],[42.74253,7.0],[42.9386,7.0],[43.13467,7.0],[43.33073,7.0],[43.5268,7.0],[43.72287,7.0],[43.91893,7.0],[44.115,7.0],[44.31107,7.0],[44.50713,7.0],[44.7032,7.0],[44.89927,7.0],[45.09533,7.0],[45.2914,7.0],[45.48747,7.0],[45.68353,7.0],[45.8796,7.0],[46.07567,7.0],[46.27173,7.0],[46.4678,7.0],[46.66387,7.0],[46.85993,7.0],[47.056,7.0]],"source":"20Nov STEP dimensions; axisymmetric envelope closes side holes. One representative suction branch."},{"name":"0.75\u2033 \u00b7 short \u00b7 printed","family":"branched","dn":3.5,"db":16,"d":22.5,"de":16,"nports":1,"df":16,"xn":-62.16,"xe":-39.5326,"xout":-32.4542,"insert":false,"id":4,"file":"0.75'' - SHORT - 3D printed.step","hash":"831ec6fd67f0e776d9952f5b0efdd555a9312f04215daa8fb6acbb2ce2046b55","L":22.6274,"key":"nov4","chamber":22.5,"angle":45,"end":29.705799999999996,"profile":[[0,11.25],[22.6274,11.25],[29.705799999999996,8.0]],"source":"20Nov STEP dimensions; axisymmetric envelope closes side holes. One representative suction branch."},{"name":"1\u2033 \u00b7 dual suction \u00b7 short","family":"branched","dn":3.5,"db":15.9,"d":22.5,"de":22.5,"nports":2,"df":15.9,"xn":-57.493,"xe":-29.9437,"xout":-29.9437,"insert":false,"id":5,"file":"1'' Eductor dual suction SHORT.step","hash":"2f5a509a8bb37b38b9f4f648afff39246f9455002406ad875de0ac547ca52d07","L":27.5493,"key":"nov5","chamber":22.5,"angle":45,"end":27.549300000000002,"profile":[[0,11.25],[27.5493,11.25],[27.549300000000002,11.25]],"source":"20Nov STEP dimensions; axisymmetric envelope closes side holes. One representative suction branch."},{"name":"1\u2033 \u00b7 four suction","family":"branched","dn":3.5,"db":14,"d":22.5,"de":22.5,"nports":4,"df":14,"xn":-57.493,"xe":0,"xout":0,"insert":false,"id":6,"file":"1'' Eductor 4 suction.step","hash":"cc059c87f555a10f6d866e7317e8a388f83d9d7866b57793e737025241618f60","L":57.493,"key":"nov6","chamber":22.5,"angle":45,"end":57.493,"profile":[[0,11.25],[57.493,11.25],[57.493,11.25]],"source":"20Nov STEP dimensions; axisymmetric envelope closes side holes. One representative suction branch."},{"name":"1\u2033 \u00b7 four suction \u00b7 short","family":"branched","dn":3.5,"db":14,"d":22.5,"de":22.5,"nports":4,"df":14,"xn":-57.493,"xe":-34.7721,"xout":-34.7721,"insert":false,"id":7,"file":"1'' Eductor 4 suction SHORT.step","hash":"49263cce89b82f2d275148346730af0a94a94e372471c3a5de3c63b693002bc0","L":22.7209,"key":"nov7","chamber":22.5,"angle":45,"end":22.7209,"profile":[[0,11.25],[22.7209,11.25],[22.7209,11.25]],"source":"20Nov STEP dimensions; axisymmetric envelope closes side holes. One representative suction branch."}]
;const AUDIT=[
  {
    "key": "revA",
    "name": "Rev 2 · 1-stage",
    "states": [
      {
        "atm": 1,
        "status": "screened",
        "reason": "Existing reduced flow estimate passes its screening checks.",
        "air": 0.023752889009181863
      },
      {
        "atm": 0.5,
        "status": "screened",
        "reason": "Existing reduced flow estimate passes its screening checks.",
        "air": 0.016322601458980332
      },
      {
        "atm": 0.2,
        "status": "compressible",
        "reason": "Automatic compressible mixing estimate; uncalibrated complete-mixing closure.",
        "air": 0.01557414125018068
      }
    ]
  },
  {
    "key": "revB1",
    "name": "Rev 2 · 2-stage: first stage",
    "states": [
      {
        "atm": 1,
        "status": "screened",
        "reason": "Existing reduced flow estimate passes its screening checks.",
        "air": 0.020861061907178645
      },
      {
        "atm": 0.5,
        "status": "screened",
        "reason": "Existing reduced flow estimate passes its screening checks.",
        "air": 0.013299184283503209
      },
      {
        "atm": 0.2,
        "status": "compressible",
        "reason": "Automatic compressible mixing estimate; uncalibrated complete-mixing closure.",
        "air": 0.012873871223080497
      }
    ]
  },
  {
    "key": "revB2",
    "name": "Rev 2 · 2-stage: second stage",
    "states": [
      {
        "atm": 1,
        "status": "screened",
        "reason": "Existing reduced flow estimate passes its screening checks.",
        "air": 0.02758981184228046
      },
      {
        "atm": 0.5,
        "status": "screened",
        "reason": "Existing reduced flow estimate passes its screening checks.",
        "air": 0.019084475300492256
      },
      {
        "atm": 0.2,
        "status": "compressible",
        "reason": "Automatic compressible mixing estimate; uncalibrated complete-mixing closure.",
        "air": 0.018093086889201894
      }
    ]
  },
  {
    "key": "revC",
    "name": "British · slide geometry",
    "states": [
      {
        "atm": 1,
        "status": "compressible",
        "reason": "Automatic compressible mixing estimate; uncalibrated complete-mixing closure.",
        "air": 0.026866443242541168
      },
      {
        "atm": 0.5,
        "status": "compressible",
        "reason": "Automatic choked-mixing / internal-shock estimate; uncalibrated homogeneous flow.",
        "air": 0.01803073610832394
      },
      {
        "atm": 0.2,
        "status": "compressible",
        "reason": "Automatic choked-mixing / internal-shock estimate; uncalibrated homogeneous flow.",
        "air": 0.01001426898096781
      }
    ]
  },
  {
    "key": "nov0",
    "name": "British 1″ · straight bore",
    "states": [
      {
        "atm": 1,
        "status": "screened",
        "reason": "Existing reduced flow estimate passes its screening checks.",
        "air": 0.032947123897068255
      },
      {
        "atm": 0.5,
        "status": "screened",
        "reason": "Existing reduced flow estimate passes its screening checks.",
        "air": 0.022804502981520675
      },
      {
        "atm": 0.2,
        "status": "compressible",
        "reason": "Automatic compressible mixing estimate; uncalibrated complete-mixing closure.",
        "air": 0.014538490107987877
      }
    ]
  },
  {
    "key": "nov1",
    "name": "British 1″ · C–D",
    "states": [
      {
        "atm": 1,
        "status": "compressible",
        "reason": "Automatic compressible mixing estimate; uncalibrated complete-mixing closure.",
        "air": 0.026710950758591864
      },
      {
        "atm": 0.5,
        "status": "compressible",
        "reason": "Automatic compressible mixing estimate; uncalibrated complete-mixing closure.",
        "air": 0.017952022240875734
      },
      {
        "atm": 0.2,
        "status": "compressible",
        "reason": "Automatic choked-mixing / internal-shock estimate; uncalibrated homogeneous flow.",
        "air": 0.009966734680507546
      }
    ]
  },
  {
    "key": "nov2",
    "name": "British original · C–D",
    "states": [
      {
        "atm": 1,
        "status": "compressible",
        "reason": "Automatic compressible mixing estimate; uncalibrated complete-mixing closure.",
        "air": 0.009097377228985213
      },
      {
        "atm": 0.5,
        "status": "compressible",
        "reason": "Automatic compressible mixing estimate; uncalibrated complete-mixing closure.",
        "air": 0.006046924130753884
      },
      {
        "atm": 0.2,
        "status": "compressible",
        "reason": "Automatic choked-mixing / internal-shock estimate; uncalibrated homogeneous flow.",
        "air": 0.0031404998065799102
      }
    ]
  },
  {
    "key": "nov3",
    "name": "British original · straight bore",
    "states": [
      {
        "atm": 1,
        "status": "screened",
        "reason": "Existing reduced flow estimate passes its screening checks.",
        "air": 0.010365870259041458
      },
      {
        "atm": 0.5,
        "status": "screened",
        "reason": "Existing reduced flow estimate passes its screening checks.",
        "air": 0.0070693774234393296
      },
      {
        "atm": 0.2,
        "status": "compressible",
        "reason": "Automatic compressible mixing estimate; uncalibrated complete-mixing closure.",
        "air": 0.004762751871662772
      }
    ]
  },
  {
    "key": "nov4",
    "name": "0.75″ · short · printed",
    "states": [
      {
        "atm": 1,
        "status": "screened",
        "reason": "Existing reduced flow estimate passes its screening checks.",
        "air": 0.008235240154784246
      },
      {
        "atm": 0.5,
        "status": "compressible",
        "reason": "Automatic compressible mixing estimate; uncalibrated complete-mixing closure.",
        "air": 0.012208614585769897
      },
      {
        "atm": 0.2,
        "status": "compressible",
        "reason": "Automatic compressible mixing estimate; uncalibrated complete-mixing closure.",
        "air": 0.008055199467608794
      }
    ]
  },
  {
    "key": "nov5",
    "name": "1″ · dual suction · short",
    "states": [
      {
        "atm": 1,
        "status": "screened",
        "reason": "Existing reduced flow estimate passes its screening checks.",
        "air": 0.010483336835682151
      },
      {
        "atm": 0.5,
        "status": "screened",
        "reason": "Existing reduced flow estimate passes its screening checks.",
        "air": 0.007347370531561289
      },
      {
        "atm": 0.2,
        "status": "compressible",
        "reason": "Automatic compressible mixing estimate; uncalibrated complete-mixing closure.",
        "air": 0.010245130189473466
      }
    ]
  },
  {
    "key": "nov6",
    "name": "1″ · four suction",
    "states": [
      {
        "atm": 1,
        "status": "screened",
        "reason": "Existing reduced flow estimate passes its screening checks.",
        "air": 0.016424944145573595
      },
      {
        "atm": 0.5,
        "status": "screened",
        "reason": "Existing reduced flow estimate passes its screening checks.",
        "air": 0.0119788598510981
      },
      {
        "atm": 0.2,
        "status": "screened",
        "reason": "Existing reduced flow estimate passes its screening checks.",
        "air": 0.008233282230204819
      }
    ]
  },
  {
    "key": "nov7",
    "name": "1″ · four suction · short",
    "states": [
      {
        "atm": 1,
        "status": "screened",
        "reason": "Existing reduced flow estimate passes its screening checks.",
        "air": 0.008759307129248294
      },
      {
        "atm": 0.5,
        "status": "screened",
        "reason": "Existing reduced flow estimate passes its screening checks.",
        "air": 0.006116566456976122
      },
      {
        "atm": 0.2,
        "status": "compressible",
        "reason": "Automatic compressible mixing estimate; uncalibrated complete-mixing closure.",
        "air": 0.009240718975670643
      }
    ]
  }
]
;const SLIDE={"Gref":3829961.607121514,"etaHalfUm":2.69,"p0bar":3,"poutbar":1,"rows":[["300 nm",-7,1.47,0.22,18,36,35,69,64,118,61,113,79],["300 nm",0,1.65,0.17,91,35,174,67,266,68,255,66,12],["300 nm",7,1.47,0.18,61,30,117,58,219,65,210,62,10],["10 \u00b5m",-7,1.13,0.01,0,0,0,0,0,0,0,0,100],["10 \u00b5m",0,1.62,0.011,0.0085,0.075,0.0068,0.06,27,18,14,9.4,0],["10 \u00b5m",7,1.55,0.0084,43,5.4,35,4.3,99,24,50,12,0]],"note":"Rounded reference values from the attached slide; CFD-field stochastic tracking, not experimental validation. Analytical calculations are independent.","definitions":"Accumulated time above q Gref; N_l=q Gref t_q min[1,(eta_q/d_p)^(2/3)]; eta_q=sqrt(nu_q/(q Gref)).","provenance":["eductor-quarter-eddy-counts/final-reviewed/build-manifest.json","eductor-half-maximum-exposure/deck/british/data/case-normalized-stochastic.json"]};// UI powder feed is g/s; existing solver inputs and archived URLs use kg/h.
// Keep exact model values through rounded preset displays until the user edits them.
