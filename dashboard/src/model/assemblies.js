const ExposureAssemblies=(()=>{
 function geometry(g,p={}){
  const mm=1000*(p.ext||0);
  if(!(mm>0)||!['revA','revB2'].includes(g.key))return g;
  return {...g,L:g.L+mm,end:g.end+mm,profile:[...g.profile.map(a=>[...a]),[g.end+mm,g.profile.at(-1)[1]]],extensionMM:mm,source:g.source+' Added straight outlet: '+mm+' mm (analytical geometry).'};
 }
 function parameters(key,inputs){
  if(!key.startsWith('rev'))return {};
  const two=key.startsWith('revB');
  return {ext:key==='revC'?0:inputs.extension/1000,feedDrop:1000*(two?inputs.feederDrop:inputs.valveDrop),interDrop:two?1000*inputs.valveDrop:0,...(two?{p01:inputs.p01}:{})};
 }
 return {geometry,parameters};
})();
if(typeof module!=='undefined')module.exports=ExposureAssemblies;

'use strict';
// Execution-only acceleration. Cache keys include the build, complete geometry,
// every setting, and the solved flow. Never round coordinates or physics inputs.
