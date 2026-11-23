const ExposureOverlay=(()=>{
 let saved=null,enabled=false;
 function snapshot(run,points){if(!run?.r?.ok)return null;return structuredClone({geometry:run.g.key,name:run.g.name,settings:run.p,Gref:run.r.Gref,codeSHA256:run.r.manifest?.codeSHA256,capturedAt:new Date().toISOString(),points,summary:{Gref:run.r.Gref,etaHalf:run.r.etaHalf,outletAir:run.r.pressure?.mg,turbulentPower:run.r.field?.lossPower},rows:run.r.rows.map(({species,d,offset,n,stationShear,stationCoagulation})=>({species,d,offset,n,stationShear,stationCoagulation}))});}
 function store(value){if(!value)return null;saved=structuredClone(value);enabled=true;return saved;}
 function keep(run,points,visuals){const q=snapshot(run,points);if(q&&visuals)q.visuals=visuals;return store(q);}
 function get(){return enabled?saved:null;}
 function toggle(value){enabled=!!value&&!!saved;return enabled;}
 return {snapshot,store,keep,get,peek:()=>saved,toggle,has:()=>!!saved,clear:()=>{saved=null;enabled=false;}};
})();
if(typeof module!=='undefined')module.exports=ExposureOverlay;

// UI configuration for the nominal earlier assemblies. Pressure drops are Pa;
// the outlet extension is metres internally. STEP geometries are not modified.
