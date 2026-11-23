const ExposureOverlay=(()=>{
 let saved=[],enabled=false,seq=1;
 const clone=value=>structuredClone(value);
 function snapshot(run,points){if(!run?.r?.ok)return null;return {geometry:run.g.key,name:run.g.name,settings:run.p,Gref:run.r.Gref,codeSHA256:run.r.manifest?.codeSHA256,capturedAt:new Date().toISOString(),points,summary:{Gref:run.r.Gref,etaHalf:run.r.etaHalf,outletAir:run.r.pressure?.mg,turbulentPower:run.r.field?.lossPower},rows:run.r.rows.map(({species,d,offset,n,stationShear,stationCoagulation})=>({species,d,offset,n,stationShear,stationCoagulation}))};}
 function store(value){if(!value)return null;const item=clone(value);item.id=seq++;item.colorIndex=saved.length;saved.push(item);enabled=true;return item;}
 function keep(run,points,visuals){const q=snapshot(run,points);if(!q)return null;if(visuals)q.visuals=visuals;return store(q);}
 function list(){return saved;}
 function get(){return enabled?saved.slice():[];}
 function toggle(value){enabled=!!value&&saved.length>0;return enabled;}
 function remove(id){saved=saved.filter(item=>item.id!==id);if(!saved.length)enabled=false;saved.forEach((item,i)=>item.colorIndex=i);}
 return {snapshot,store,keep,get,peek:()=>saved.at(-1)||null,list,toggle,has:()=>saved.length>0,clear:()=>{saved=[];enabled=false;},remove};
})();
if(typeof module!=='undefined')module.exports=ExposureOverlay;

// UI configuration for the nominal earlier assemblies. Pressure drops are Pa;
// the outlet extension is metres internally. STEP geometries are not modified.
