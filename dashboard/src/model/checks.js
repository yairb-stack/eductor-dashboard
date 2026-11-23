const ExposureChecks=(()=>{
 function classify(f,p){
  const reason=f.reason||'The calculation did not return a result.';
  if(f.stage==='particles')return {code:'particles',title:'Particle model unavailable',explanation:reason+' The gas calculation and the particle calculation have separate validity limits.'};
  if(/positive|required|above outlet|above the shared|Prescribed mode/i.test(reason))return {code:'inputs',title:'Check the operating inputs',explanation:reason};
  if(/chok|sonic|compressible|expanded jet/i.test(reason))return {code:'compressible',title:'Compressible estimate unresolved',explanation:reason+' The automatic estimate cannot resolve every shock, restriction or mixing regime; this is not proof that the hardware cannot operate.'};
  if(/Bleed outflow/i.test(reason))return {code:'reverse-bleed',title:'Bleed flow reverses',explanation:'The model predicts outward flow through the bleed holes. Powder splitting through those openings is not modeled, so particle and dilution predictions are withheld.'};
  if(/No forward/i.test(reason))return {code:'no-forward',title:'No forward suction at these settings',explanation:p.po>p.ps?'The selected outlet is above the vessel pressure. This reduced model cannot drive a forward powder-carrying stream against that pressure rise.':'The balance does not support forward suction at these settings. Short mixing lengths or outlet restrictions can exhaust the modeled pressure recovery.'};
  return {code:'unresolved',title:'Pressure balance is unresolved',explanation:'The reduced solver found no admissible root. This is a model limit, not proof that the physical eductor cannot operate.'};
 }
 function diagnose(g,p,f,engine,legacy){
  const result={...classify(f,p),reason:f.reason,settings:{p0:p.p0,ps:p.ps,po:p.po},alternatives:[],threshold:null};
  const threshold=f.raw?.balance?.zeroSecondaryVesselPressure;
  if(Number.isFinite(threshold)&&threshold>0)result.threshold=threshold/1e5;
  if(p.mode!=='auto'||f.stage==='particles'||f.stage==='inputs')return result;
  // Prefer checked motive-only alternatives; do not silently reset the vacuum.
  const candidates=[.5,.75,1,1.25,1.5,2,2.5,3,4,5,6].filter(v=>v>p.ps&&v!==p.p0).sort((a,b)=>Math.abs(Math.log(a/p.p0))-Math.abs(Math.log(b/p.p0)));
  const options=candidates.map(p0=>({label:'Keep vessel = outlet '+p.ps+' bar · motive '+p0+' bar',changes:{p0,ps:p.ps,po:p.ps},sameBoundary:true}));
  options.push({label:'Use motive 3 bar · vessel = outlet 1 bar',changes:{p0:3,ps:1,po:1}});
  const seen=new Set();
  for(const option of options){
   if(option.sameBoundary&&result.alternatives.filter(v=>v.sameBoundary).length>=2)continue;
   const q={...p,...option.changes},key=[q.p0,q.ps,q.po].join('/');if(seen.has(key)||key===[p.p0,p.ps,p.po].join('/'))continue;seen.add(key);
   const flow=engine.flow(g,q,legacy);if(!flow.ok)continue;
   const field=engine.makeField(g,q,flow);if(!(field.peak>0&&Number.isFinite(field.peak)))continue;
   const releaseOK=[q.small,q.large].every(d=>field.geo.clearance(...field.geo.release(0),d*1e-6)>0);
   if(!releaseOK)continue;
   result.alternatives.push({...option,air:flow.ms*1000,chamber:flow.portP/1e5});
  }
  return result;
 }
 return {classify,diagnose};
})();
if(typeof module!=='undefined')module.exports=ExposureChecks;

const $=id=>document.getElementById(id),COLORS=['#008678','#08678b','#bf5015'];
PowderRate.init($('feed'));
