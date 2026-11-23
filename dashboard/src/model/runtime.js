const ExposureRuntime=(()=>{
 function key(value){if(Array.isArray(value))return '['+value.map(key).join(',')+']';if(value&&typeof value==='object')return '{'+Object.keys(value).sort().map(k=>JSON.stringify(k)+':'+key(value[k])).join(',')+'}';return value===undefined?'undefined':Object.is(value,-0)?'-0':JSON.stringify(value);}
 class LRU{
  constructor(limit=3,bytes=48*1024*1024){this.limit=limit;this.budget=bytes;this.bytes=0;this.map=new Map();}
  get(k){const entry=this.map.get(k);if(!entry)return;this.map.delete(k);this.map.set(k,entry);return entry.value;}
  set(k,value,bytes=1){if(this.map.has(k)){this.bytes-=this.map.get(k).bytes;this.map.delete(k);}if(bytes>this.budget)return false;while(this.map.size&&(this.map.size>=this.limit||this.bytes+bytes>this.budget)){const k=this.map.keys().next().value;this.bytes-=this.map.get(k).bytes;this.map.delete(k);}this.map.set(k,{value,bytes});this.bytes+=bytes;return true;}
  clear(){this.map.clear();this.bytes=0;}
 }
 function estimate(result){let bytes=0;for(const row of result.rows){const {samples,...rest}=row;bytes+=JSON.stringify(rest).length*4;const sizes=[samples[0],samples.at(-1),samples[Math.floor(samples.length/2)]].filter(Boolean).map(s=>JSON.stringify(s).length);bytes+=Math.max(0,...sizes)*samples.length*4;}return bytes;}
 function create(options={}){
  const E=options.model||(typeof Exposure!=='undefined'?Exposure:null),W=options.WorkerClass??globalThis.Worker,build=options.buildHash||'unversioned',url=options.workerURL||'/model-assets/exposure-worker.js?v='+build;
  const cores=Math.max(1,Math.floor(options.hardwareConcurrency??globalThis.navigator?.hardwareConcurrency??2)),memory=options.deviceMemory??globalThis.navigator?.deviceMemory??8,maxWorkers=Math.max(1,Math.min(4,cores-1||1,memory<=4?2:4,options.maxWorkers??4));
  const cache=new LRU(options.cacheEntries??3,options.cacheBytes??(memory<=4?48:96)*1024*1024),prepared=new LRU(2),slots=[];let active=null,sequence=0,workersEnabled=!!W,closed=false;
  const counters={runs:0,cacheHits:0,workersStarted:0,cancelled:0,preparedHits:0};
  function remove(slot){clearTimeout(slot.timer);slot.worker.terminate();const i=slots.indexOf(slot);if(i>=0)slots.splice(i,1);}
  function cancel(){sequence++;if(!active)return;const a=active;active=null;counters.cancelled++;for(const slot of [...slots])if(slot.job?.run===a)remove(slot);a.resolve({ok:false,cancelled:true});}
  function fail(run,error,infrastructure=false){if(active!==run)return;active=null;for(const slot of [...slots])if(slot.job?.run===run)remove(slot);error.infrastructure=infrastructure;run.reject(error);}
  function dispatch(slot){const job=slot.job;if(!job||!slot.ready||job.sent)return;job.sent=true;const {run,group}=job,message={id:run.id,group,g:run.g,p:run.p,f:run.f,key:run.key,build};
   if(run.wallData&&!slot.cad.has(run.g.key)){message.wallData=run.wallData;slot.cad.add(run.g.key);}try{slot.worker.postMessage(message);}catch(error){fail(run,error,true);}
  }
  function assign(slot,run){if(active!==run||!run.queue.length)return;slot.job={run,group:run.queue.shift(),sent:false};dispatch(slot);}
  function spawn(){const slot={worker:new W(url),ready:false,job:null,cad:new Set()};slots.push(slot);counters.workersStarted++;
   slot.timer=setTimeout(()=>{if(slot.job)fail(slot.job.run,Error('Parallel engine did not start'),true);else remove(slot);},15000);
   slot.worker.onerror=e=>{e.preventDefault?.();if(slot.job)fail(slot.job.run,Error(e.message||'Parallel engine unavailable'),true);else remove(slot);};
   slot.worker.onmessageerror=()=>{if(slot.job)fail(slot.job.run,Error('Parallel result could not be transferred'),true);};
   slot.worker.onmessage=event=>{const m=event.data;if(m.ready){clearTimeout(slot.timer);if(m.build!==build){if(slot.job)fail(slot.job.run,Error('Model version changed; reload for parallel processing'),true);else remove(slot);return;}slot.ready=true;dispatch(slot);return;}
    const job=slot.job,run=job?.run;if(!run||active!==run||m.id!==run.id||m.group!==job.group)return;
    if(m.error){fail(run,Error(m.error));return;}
    if(m.progress){run.counts[m.group]=m.progress.completed;run.onProgress({completed:run.counts.reduce((a,b)=>a+b,0),total:6*run.limit,phase:m.progress.phase,workers:run.workers});return;}
    if(!m.result?.ok){fail(run,Error(m.result?.reason||'Particle model did not return a valid result'));return;}
    run.rows[m.group]=m.result.rows[0];run.meta=m.result;run.counts[m.group]=m.result.rows[0].n;slot.job=null;run.done++;
    if(run.done===6){active=null;const {rows,meta}=run;run.resolve({ok:true,rows,Gref:meta.Gref,etaHalf:meta.etaHalf,etaQuarter:meta.etaQuarter,workers:run.workers});}
    else assign(slot,run);
   };return slot;
  }
  async function simulate(g,p,f,{onProgress=()=>{},cancelled=()=>false}={}){
   if(closed)throw Error('Calculation engine was disposed');cancel();const started=performance.now(),id=++sequence,k=key({build,g,p,f});counters.runs++;
   let field=prepared.get(k);if(field)counters.preparedHits++;else{field=E.makeField(g,p,f);prepared.set(k,field);}
   const hit=cache.get(k);if(hit){counters.cacheHits++;return {...structuredClone(hit),field,pressure:f,settings:{...p},execution:{mode:'cache',workers:0,cacheHit:true,milliseconds:performance.now()-started}};}
   let result,fallbackReason=null;
   if(workersEnabled){try{const count=Math.min(maxWorkers,p.wallModel==='cad'?2:4);while(slots.length>count)remove(slots.at(-1));result=await new Promise((resolve,reject)=>{const run={id,key:k,g,p,f,wallData:options.wallData?.(g.key),queue:[0,1,2,3,4,5],rows:Array(6),counts:Array(6).fill(0),limit:p.noise?(p.sampling==='precision'?(p.maxHistories??4096):p.n):(field.geo.branches?.length||1),done:0,workers:count,resolve,reject,onProgress};active=run;try{while(slots.length<count)spawn();for(const slot of slots)assign(slot,run);}catch(error){fail(run,error,true);}});}catch(error){if(!error.infrastructure)throw error;workersEnabled=false;fallbackReason=error.message;}}
   if(result?.cancelled||cancelled()||id!==sequence)return {ok:false,cancelled:true};
   if(!result){result=await E.simulateAsync(g,p,f,{field,cancelled:()=>cancelled()||id!==sequence,onProgress});if(result.cancelled||cancelled()||id!==sequence)return {ok:false,cancelled:true};}
   if(!result.ok)return result;const data={ok:true,rows:result.rows,Gref:result.Gref,etaHalf:result.etaHalf,etaQuarter:result.etaQuarter};const bytes=estimate(data);if(bytes<=cache.budget)cache.set(k,structuredClone(data),bytes);
   return {...data,field,pressure:f,settings:{...p},execution:{mode:result.workers?'parallel':'single',workers:result.workers||1,cacheHit:false,milliseconds:performance.now()-started,...(fallbackReason?{fallbackReason}:{})}};
  }
  function dispose(){cancel();closed=true;for(const slot of [...slots])remove(slot);cache.clear();prepared.clear();}
  return {simulate,cancel,dispose,clearCache:()=>cache.clear(),stats:()=>({...counters,workers:slots.length,maxWorkers,cachedRuns:cache.map.size,cacheBytes:cache.bytes})};
 }
 return {key,LRU,estimate,create};
})();
if(typeof module!=='undefined')module.exports=ExposureRuntime;
