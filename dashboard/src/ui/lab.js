let current=null,playing=false,animation=0,started=0,updateTimer=0,runVersion=0;
const fmt=(v,n=2)=>Number.isFinite(v)?v.toLocaleString('en-US',{minimumFractionDigits:n,maximumFractionDigits:n}):'—';
const sizeLabel=d=>d<1e-6?fmt(d*1e9,0)+' nm':fmt(d*1e6,d*1e6%1?1:0)+' µm';
const stat=(v,scale=1,digits=1)=>fmt(v.mean*scale,digits)+' ± '+fmt(v.sd*scale,digits);
function syncCollisionControls(){const hybrid=$('collisionModel').value==='hybrid';$('kernelRange').disabled=!hybrid;$('kernelConvention').disabled=!hybrid;}
$('collisionModel').addEventListener('input',syncCollisionControls);syncCollisionControls();
function settings(){$('po').value=$('ps').value;const p={modelVersion:2,releaseMode:'inlet'};for(const id of ['p0','ps','po','feed','temp','small','large','density','air','bleed','chamber','release','offset','intensity','lengthScale','spread','n','dt','horizon','seed','loss'])p[id]=id==='feed'?PowderRate.read($(id)):+$(id).value;p.alpha=+$('alpha').value;p.smallShare=+$('smallShare').value/100;p.background=$('background').value;p.collisionModel=$('collisionModel').value;p.mode=$('flow-mode').value;p.reference=$('reference').value;p.noise=$('noise').checked;p.sampling=$('sampling').value;p.targetHalfWidth=.02;p.maxHistories=4096;p.adaptiveDt=$('adaptiveDt').checked;p.antithetic=$('antithetic').checked;p.commonNoise=$('commonNoise').checked;p.kernelConvention=$('kernelConvention').value;p.extendedKernels=p.collisionModel==='hybrid'&&$('kernelRange').value!=='baseline';p.clustering=p.collisionModel==='hybrid'&&$('kernelRange').value==='clustered';p.populationBalance=$('populationBalance').checked;p.wallModel=$('wallModel').value;p.phaseModel=$('phase-model').value;if(p.phaseModel!=='stable'){p.flowPhysics='phase2';p.dispersion='well-mixed';}if(p.phaseModel==='phase3'){p.trainClosure='staged';p.calibrationCandidate=true;p.spreadingLaw='compressible-fit';}Object.assign(p,ExposureAssemblies.parameters($('geometry').value,Object.fromEntries(['p01','valveDrop','feederDrop','extension'].map(id=>[id,+$(id).value]))));return p;}
function syncFlowControls(){const key=$('geometry').value,two=key.startsWith('revB');$('assembly-config').hidden=!key.startsWith('rev');document.querySelector('label[for=p0]').textContent=two?'Second motive · bar abs':'Motive · bar abs';$('first-motive-field').hidden=$('feeder-drop-field').hidden=!two;$('extension-field').hidden=key==='revC';$('valve-scope').textContent=two?'Between the two stages.':'Powder supply line.';$('assembly-scope').textContent=two?'Coupled two-stage gas balance; particle histories are local to the selected stage.':'Nominal earlier-assembly geometry.';const cad=$('geometry').value.startsWith('nov');$('release').disabled=cad;$('release').title=cad?'CAD cases release at the common suction inlet; passage length comes from STEP.':'Distance upstream of the nominal junction';const manual=$('flow-mode').value==='prescribed';$('prescribed-fields').hidden=!manual;$('loss').closest('.field').hidden=manual||!$('geometry').value.startsWith('nov');const out=$('solved-flow');if(out)out.hidden=manual;}
function markDirty(){ExposureCompute.cancel();document.querySelector('label[for=dt]').textContent=$('adaptiveDt').checked?'Base timestep · µs (up to 20×)':'Maximum timestep · µs';runVersion++;clearTimeout(updateTimer);document.body.classList.add('model-pending');$('busy-note').textContent='Updating automatically…';$('csv').disabled=$('json').disabled=true;const raw=$('samples-csv');if(raw)raw.disabled=true;syncFlowControls();$('smallShare').closest('.field').hidden=$('background').value!=='mixture';stop();updateTimer=setTimeout(run,350);}

function stop(){playing=false;cancelAnimationFrame(animation);$('play').textContent='Play histories';}
function getGeometry(){return ExposureAssemblies.geometry(GEOMETRIES.find(g=>g.key===$('geometry').value),settings());}
function unavailable(g,p,f){
 drawnField=null;current=null;$('csv').disabled=$('json').disabled=true;
 const issue=ExposureChecks.diagnose(g,p,f,Exposure,Legacy);
 $('flow-status').className='notice error';$('flow-status').textContent=issue.title+' · details and checked alternatives are shown in the result panel.';
 $('flow-status').dataset.details=issue.explanation;$('busy-note').textContent=issue.title;
 window.exposureUnavailable?.(g,issue);
}
function run(){
 ExposureCompute.cancel();const version=++runVersion;clearTimeout(updateTimer);stop();syncFlowControls();$('smallShare').closest('.field').hidden=$('background').value!=='mixture';
 const p=settings(),g=getGeometry();
 for(const el of document.querySelectorAll('.experiment input[type=number],main>section input[type=number]'))if(!($('flow-mode').value==='auto'&&$('prescribed-fields').contains(el))&&(!el.checkValidity()||el.value==='')){
  unavailable(g,p,{stage:'inputs',reason:'A valid '+(el.labels?.[0]?.textContent||'value')+' is required.'});document.body.classList.remove('model-pending');$('run').disabled=false;return;
 }
 $('run').disabled=true;$('busy-note').textContent='Computing particle histories…';
 setTimeout(async()=>{if(version!==runVersion)return;let gasSolved=false;try{
  if(p.wallModel==='cad'&&g.key.startsWith('nov')&&g.key!=='nov4'&&!EductorWalls.has(g.key)){const response=await fetch('../model-assets/walls/'+g.key+'.json');if(!response.ok)throw Error('Native CAD wall data could not be loaded. Select Reduced passages or retry.');EductorWalls.register(g.key,await response.json());if(version!==runVersion)return;}

  const f=Exposure.flow(g,p,Legacy);if(!f.ok){unavailable(g,p,f);return;}gasSolved=true;
  const r=await ExposureCompute.simulate(g,p,f,{cancelled:()=>version!==runVersion,onProgress:q=>{if(version===runVersion)$('busy-note').textContent=(q.phase==='statistics'?'Summarizing':'Computing')+' · '+q.completed.toLocaleString()+' / '+q.total.toLocaleString()+' histories…';}});if(version!==runVersion||r.cancelled)return;
  if(!r.ok){unavailable(g,p,{...r,stage:'particles'});return;}r.manifest={...EDUCTOR_BUILD,created:new Date().toISOString(),runtime:navigator.userAgent,execution:r.execution,geometry:g.key,settings:{...p},histories:r.rows.map(x=>x.n),seed:p.seed,powderInputUnit:'kg/h internal; g/s display'};current={g,p,r};render();$('busy-note').textContent='Updated · '+r.rows.reduce((s,row)=>s+row.n,0)+' particle histories · '+(r.execution.cacheHit?'cached':r.execution.workers+' CPU worker'+(r.execution.workers===1?'':'s'))+' · '+fmt(r.execution.milliseconds/1000,2)+' s.';$('csv').disabled=$('json').disabled=false;
 }catch(e){if(version!==runVersion)return;unavailable(g,p,{reason:e.message,stage:gasSolved?'particles':'calculation'});}finally{if(version===runVersion){$('run').disabled=false;document.body.classList.remove('model-pending');}}},20);
}

function render(){
 const {g,p,r}=current,f=r.pressure;$('results').hidden=false;$('flow-status').className='notice';
 const text=[f.reason,f.manual?'Chamber pressure and carrier flows are prescribed inputs.':'Feeder flow '+fmt(f.ms*1000,3)+' g/s; bleed flow '+fmt(f.mh*1000,3)+' g/s; chamber pressure '+fmt(f.portP/1e5,5)+' bar absolute. Gas mass balance: outlet = motive + incoming + bleed. The reduced model includes bulk powder momentum; the spatial field is not a fully coupled compressible solution.',...(f.compressible?[f.gate.reason,f.network.assumptions,'Outlet pressure residual: '+fmt(f.network.pressureResidual,6)+' Pa.']:[]),'Particle histories are 3D, one-way coupled analytical estimates. First wall contact is recorded; rebound/deposition is unknown, not inferred. Local concentration is flux-normalized; bulk powder is conserved.'];
 if(f.md/f.ms>.1)text.push('Powder / incoming air = '+fmt(f.md/f.ms,2)+': particle feedback on the spatial flow can matter.');
 if(g.key.startsWith('revB'))text.push('Stage-local histories: particles are released anew here, not propagated through the whole two-stage train.');
 if(g.nports>1)text.push(g.nports+' passages are sampled by predicted powder flux. Curved CAD centrelines set transport and passage lengths.');
 if(f.circuit?.flowUsed)text.push('Suction branches: '+f.circuit.branches.map((b,i)=>`${i+1}: ${(b.m*1000).toFixed(3)} g/s, ${(100*b.fraction).toFixed(1)}% of powder, ${b.lengthMM.toFixed(1)} mm path`).join('; ')+'. Common feeder losses are included.');
 $('flow-status').dataset.details=[...text,...(p.phaseModel!=='stable'?['Candidate field: phase 2/3 gates have not all passed. Consult Model review before interpreting differences.']:[])].join('\n\n');$('flow-status').innerHTML=`<span class="flow-mode-label">${f.manual?'Prescribed-flow scenario':f.compressible?(f.network.shock?'Choked flow · diffuser shock':'Compressible estimate active'):'Low-Mach balance solved'}</span><span class="flow-loading">Powder / inlet air <b>${fmt(f.md/f.ms,2)}</b> kg/kg</span><span class="flow-scope">${f.compressible?'Mass + momentum + energy balanced':'One-way particle model'}</span><button type="button" class="flow-details" aria-haspopup="dialog">Assumptions</button>`;$('result-title').textContent=g.name;
 $('run-caption').textContent=`${g.key.startsWith('revB')?'Motives: first '+fmt(p.p01??p.p0,3)+' / second '+fmt(p.p0,3):'Motive '+fmt(p.p0,3)} bar abs · vessel = outlet ${fmt(p.ps,3)} bar abs · powder ${fmt(PowderRate.fromModel(p.feed),3)} g/s · ${p.noise?'seed '+p.seed:'deterministic mean field'} · ${p.mode==='auto'?'predicted gas flow':'prescribed gas-flow scenario'}`;
 $('p0').parentElement.querySelector('small').hidden=!g.key.startsWith('revB');$('geometry-source').textContent=g.source+' '+(r.field.geo.wallStatus||'Reduced wall model.')+' Revision 2: CAD-centreline tube unions for 20Nov; nominal passages for historical geometries. Junction shape and wall contact remain approximations. Bleed-hole escape is unresolved.';
 $('kpis').innerHTML=[['Reference shear · s⁻¹',fmt(r.Gref/1e6,2)+' ×10⁶','s⁻¹ · '+(p.reference==='slide'?'fixed slide reference':'analytic maximum at x/D = '+fmt(r.field.peakX/(g.dn/1000),2))],['Half scale · µm',fmt(r.etaHalf*1e6,2),'µm · ν at analytical jet maximum'],['Outlet air · g/s',fmt(f.mg*1000,2),'g/s · '+(p.feed?fmt(f.mg/f.md,2)+' kg air/kg powder':'no powder mass feed')],['Turbulent power · W',fmt(r.field.lossPower,1),'W · '+fmt(r.field.available,1)+' W available ceiling']].map(x=>`<div class="kpi"><span>${x[0]}</span><strong>${x[1]}</strong><small>${x[2]}</small></div>`).join('');
 $('path-legend').innerHTML=r.rows.slice(0,3).map((row,i)=>`<span><i class="key" style="color:${COLORS[i]}"></i>${fmt(row.offset*1000,2)} mm</span>`).join('');
 MAP.xMax=Math.max(120,g.end+12);MAP.scale=Math.min(790/(MAP.xMax+60),374/85);if(!$('field').closest('.plot-wrap').hidden)drawField();drawPaths();drawHistories();
 $('metrics').innerHTML=r.rows.map((row,i)=>`<tr><td>${sizeLabel(row.d)}</td><td style="color:${COLORS[i%3]}">${fmt(row.offset*1000,2)}</td><td>${row.n}</td><td>${stat(row.time,1000,2)}</td><td>${stat(row.half,1e6,1)}</td><td>${stat(row.nh,1,1)}</td><td>${stat(row.quarter,1e6,1)}</td><td>${stat(row.nq,1,1)}</td><td>${fmt(row.wall*100,1)}${p.noise?`<br><small>95% MC: ${fmt(row.wallCI[0]*100,0)}–${fmt(row.wallCI[1]*100,0)}</small>`:''}</td><td>${fmt(row.outlet*100,1)}</td><td>${fmt(row.censored*100,1)}</td></tr>`).join('');
 const censor=r.rows.reduce((s,row)=>s+row.censored*row.n,0),escape=r.rows.reduce((s,row)=>s+row.inlet*row.n,0);
 $('diagnostics').textContent=`Maximum particle Re: ${fmt(Math.max(...r.rows.map(row=>row.maxRe)),1)}. Inlet escapes: ${Math.round(escape)}. Censored: ${Math.round(censor)}. Dissipation amplitude factor: ${fmt(r.field.noiseScale,3)}. ${p.noise?'Wall intervals are Wilson Monte Carlo intervals for these chosen assumptions; they exclude model uncertainty.':'Fluctuations are off: one deterministic path per size/offset and no sampling interval.'} ${censor?'Censored transit times are truncated at the tracking horizon; increase it before comparing completed residence times.':''}`;
 drawBenchmark();$('time-control').max=p.horizon;$('time-control').value=p.horizon;drawTime();window.exposureRendered?.();
}
const MAP={scale:790/180,xMax:120,x:x=>55+(x*1000+60)*MAP.scale,y:y=>30+(60-y*1000)*MAP.scale};
let drawnField=null;
function drawField(){
 if(!current||drawnField===current.r.field)return;drawnField=current.r.field;
 const {r}=current,{geo}=r.field,c=$('field'),ctx=c.getContext('2d');ctx.clearRect(0,0,c.width,c.height);ctx.fillStyle='#fff';ctx.fillRect(0,0,900,510);
 const scale=MAP.scale,step=Math.max(1.15,(MAP.xMax+60)/180);
 for(let x=-60;x<MAP.xMax;x+=step)for(let y=-25;y<60;y+=step){if(geo.clearance(x/1000,y/1000)>0){const v=r.field.sample(x/1000,y/1000).ux,t=Math.min(1,Math.max(0,v/500));ctx.fillStyle=`rgb(${Math.round(224-180*t)},${Math.round(237-131*t)},${Math.round(246-104*t)})`;ctx.fillRect(MAP.x(x/1000),MAP.y((y+step)/1000),step*scale+1,step*scale+1);}}
 ctx.strokeStyle='#d5dee4';ctx.lineWidth=1;ctx.font='14px system-ui';ctx.fillStyle='#566773';
 for(let x=-60;x<=MAP.xMax;x+=Math.ceil((MAP.xMax+60)/180)*30){ctx.beginPath();ctx.moveTo(MAP.x(x/1000),30);ctx.lineTo(MAP.x(x/1000),404);ctx.stroke();ctx.fillText(x,MAP.x(x/1000)-10,424);}
 for(let y=-20;y<=60;y+=20){ctx.beginPath();ctx.moveTo(55,MAP.y(y/1000));ctx.lineTo(845,MAP.y(y/1000));ctx.stroke();ctx.fillText(y,17,MAP.y(y/1000)+5);}
 // Draw the boundary of the same union envelope used for first-contact checks.
 ctx.strokeStyle='#516774';ctx.lineWidth=1.3;ctx.beginPath();
 const ds=Math.max(.0007,(MAP.xMax+60)/1000/260);
 for(let x=-.060;x<MAP.xMax/1000;x+=ds)for(let y=-.025;y<.060;y+=ds){const v=[[x,y],[x+ds,y],[x+ds,y+ds],[x,y+ds]],d=v.map(q=>geo.clearance(...q)),hits=[];for(let j=0;j<4;j++){const k=(j+1)%4;if((d[j]>0)!==(d[k]>0)){const t=d[j]/(d[j]-d[k]);hits.push([v[j][0]+t*(v[k][0]-v[j][0]),v[j][1]+t*(v[k][1]-v[j][1])]);}}for(let j=0;j+1<hits.length;j+=2){ctx.moveTo(MAP.x(hits[j][0]),MAP.y(hits[j][1]));ctx.lineTo(MAP.x(hits[j+1][0]),MAP.y(hits[j+1][1]));}}ctx.stroke();
 ctx.fillText('x (mm)',420,451);ctx.fillText('y (mm)',8,18);ctx.fillText('Axial gas velocity Ux (m/s)',55,482);
 const grad=ctx.createLinearGradient(310,0,650,0);grad.addColorStop(0,'rgb(241,246,249)');grad.addColorStop(1,'rgb(44,106,142)');ctx.fillStyle=grad;ctx.fillRect(310,467,340,13);ctx.fillStyle='#566773';ctx.fillText('0',307,500);ctx.fillText('250',460,500);ctx.fillText('500',630,500);
}
function pathD(points){return points.map((p,i)=>(i?'L':'M')+MAP.x(p[1]).toFixed(2)+','+MAP.y(p[2]).toFixed(2)).join('');}
function drawPaths(){
 const {r}=current,geo=r.field.geo;let s='<defs><clipPath id="plotclip"><rect x="55" y="30" width="790" height="374"/></clipPath></defs><g clip-path="url(#plotclip)">';
 for(const row of r.rows)s+=`<path d="${pathD(row.meanPath)}" fill="none" stroke="#637581" stroke-width="1.3" stroke-dasharray="6 5" opacity=".6"/>`;
 for(let i=0;i<r.rows.length;i++){const row=r.rows[i];s+=`<path d="${pathD(row.path)}" fill="none" stroke="${COLORS[i%3]}" stroke-width="2.2" ${row.species?'stroke-dasharray="2 4"':''} opacity=".85"/>`;const p=row.path[0];s+=`<circle cx="${MAP.x(p[1])}" cy="${MAP.y(p[2])}" r="4" fill="${COLORS[i%3]}"/>`;}
 s+='<g id="markers"></g></g>';$('paths').innerHTML=s;
}
function drawHistories(){
 if(!current)return;const {r}=current,combined=matchMedia('(max-height:780px) and (min-width:651px)').matches;
 if($('lab-size-legend'))$('lab-size-legend').textContent='Particle charts: solid '+sizeLabel(r.rows[0].d)+' · dashed '+sizeLabel(r.rows[3].d);
 for(let species=0;species<2;species++){
  const target=$(species?'history-large':'history-small'),rows=combined?r.rows:r.rows.filter(row=>row.species===species),width=Math.max(150,target.clientWidth||300),height=Math.max(130,target.clientHeight||220);
  target.setAttribute('viewBox',`0 0 ${width} ${height}`);
  target.setAttribute('aria-label',(combined?'Both particle sizes':sizeLabel(rows[0].d))+' shear at twelve numbered locations, G / Gref from 0 to 1');
  target.innerHTML=ExposureStationChart.svg({rows,Gref:r.Gref,width,height,active:ExposureDilution.selectedIndex(),title:combined?'Both particle sizes':sizeLabel(rows[0].d),combined,id:target.id,overlay:ExposureOverlay.get()});
 }
 drawStationCoagulation();
}
function drawStationCoagulation(){
 if(!current||!$('lab-coag-quantity'))return;const {r}=current,quantity=$('lab-coag-quantity').value,combined=matchMedia('(max-height:780px) and (min-width:651px)').matches;
 for(let species=0;species<2;species++){
  const target=$(species?'lab-coag-large':'lab-coag-small'),rows=combined?r.rows:r.rows.filter(row=>row.species===species),width=Math.max(150,target.clientWidth||300),height=Math.max(130,target.clientHeight||220);
  target.setAttribute('viewBox',`0 0 ${width} ${height}`);
  target.setAttribute('aria-label',(combined?'Both particle sizes':sizeLabel(rows[0].d))+' accumulated '+(quantity==='probability'?'first-sticking probability':quantity==='screened'?'screened collision dose':'collision dose')+' at twelve numbered locations; mean and P05–P95 among reaching histories');
  target.innerHTML=ExposureStationChart.svg({rows,quantity,width,height,active:ExposureDilution.selectedIndex(),title:combined?'Both particle sizes':sizeLabel(rows[0].d),combined,id:target.id,overlay:ExposureOverlay.get()});
 }
}

function drawTime(){if(!current)return;const t=+$('time-control').value/1000;$('time-label').textContent=fmt(t*1000,2)+' ms';let s='';for(let i=0;i<current.r.rows.length;i++){const row=current.r.rows[i],ps=row.path;if(t>ps.at(-1)[0])continue;let k=ps.findIndex(v=>v[0]>=t);if(k<0)continue;const a=ps[Math.max(0,k-1)],b=ps[k],f=(t-a[0])/Math.max(b[0]-a[0],1e-12),x=a[1]+f*(b[1]-a[1]),y=a[2]+f*(b[2]-a[2]);s+=`<circle cx="${MAP.x(x)}" cy="${MAP.y(y)}" r="${row.species?5:3.5}" fill="${row.species?'white':COLORS[i%3]}" stroke="${COLORS[i%3]}" stroke-width="2"/>`;}$('markers').innerHTML=s;window.exposureTime?.(t);}
function play(){if(!current)return;if(playing){stop();return;}playing=true;started=performance.now();$('play').textContent='Pause histories';function frame(now){if(!playing)return;const f=Math.min(1,(now-started)/6000);$('time-control').value=f*current.p.horizon;drawTime();if(f<1)animation=requestAnimationFrame(frame);else stop();}animation=requestAnimationFrame(frame);}
function drawBenchmark(){
 const {g,p,r}=current,match=g.key==='revC'&&Math.abs(p.p0-3)<1e-8&&Math.abs(p.ps-1)<1e-8&&Math.abs(p.po-1)<1e-8&&Math.abs(p.small-.3)<1e-8&&Math.abs(p.large-10)<1e-8&&Math.abs(p.offset-.875)<1e-8&&Math.abs(p.release-30)<1e-8&&Math.abs(p.density-2200)<1e-8&&p.mode==='auto'&&p.temp===20;
 $('benchmark-status').textContent=match?('Geometry, pressures, sizes and release positions match the slide preset. The spatial field, gas closure and turbulence model differ; only gas-field coefficients were calibrated; particle residence and contact were not fitting targets. '+(p.reference==='model'?'Half-exposure thresholds also differ because each model uses its own Gref. Select the fixed slide reference for a common threshold.':'The half-exposure comparison uses the same fixed Gref as the slide.')):'Different configuration: the slide values are contextual only. Use the slide preset for matching geometry, pressures, sizes and release positions.';
 $('benchmark').innerHTML=SLIDE.rows.map((s,i)=>{const row=r.rows[i];return `<tr><td>${s[0]} / ${s[1]>0?'+':''}${s[1]} mm</td><td>${match?fmt(row.time.mean*1000,2):'—'} / ${fmt(s[2],2)}</td><td>${match?fmt(row.half.mean*1e6,2):'—'} / ${fmt(s[4],s[4]<.1?4:0)}</td><td>${match?fmt(row.wall*100,1):'—'} / ${s[12]}</td></tr>`;}).join('');
}
function download(kind){if(!current)return;const {g,p,r}=current;let data,type;
 if(kind==='json'){data=JSON.stringify({version:'analytical-exposure-phases1-3-coagulation',manifest:r.manifest,collisionKernel:{version:Exposure.collisions.version,model:p.collisionModel,label:Exposure.collisions.label(p.collisionModel),screen:Exposure.collisions.SCREEN,clusteringRDF:p.clustering?'Zhou contact RDF sensitivity':1,extendedKernels:p.extendedKernels,populationBalance:p.populationBalance,interpretation:'Encounter plots use a frozen background; optional fixed-pivot population balance reports separate growth/depletion on fixed transport paths. Screened dose is a subtotal.'},geometry:g,settings:p,flow:{motive:r.pressure.mp,feeder:r.pressure.ms,bleed:r.pressure.mh,outlet:r.pressure.mg,source:p.mode,circuit:r.pressure.circuit,network:r.pressure.network,calibration:r.pressure.calibration},Gref:r.Gref,nuRef:r.field.peakNu,etaHalf:r.etaHalf,limitations:'Unvalidated reduced field; 3D reduced tube-union walls; hydraulic branch split; contact is not deposition; stage-local release; one-way particle coupling; dispersion and wall model identified by settings; candidates retain open validation gates; coagulation uses the selected collision kernel in a frozen background, with explicit finite-size/drag/Mach screening and optional sectional growth/depletion; no growth feedback on transport or breakup; sticking probability is a first-event estimate.',rows:r.rows,overlay:ExposureOverlay.get()},null,2);type='application/json';}
 else{const head=['geometry','flow_source','motive_bar_abs','vessel_bar_abs','outlet_bar_abs','seed','diameter_um','offset_mm','n','tracked_time_mean_ms','tracked_time_sd_ms','half_mean_us','half_sd_us','N_half_mean','N_half_sd','quarter_mean_us','quarter_sd_us','N_quarter_mean','N_quarter_sd','wall_pct','outlet_pct','inlet_pct','censored_pct','Gref_s_inv'];const records=r.rows.map(x=>[g.name,p.mode,p.p0,p.ps,p.po,p.seed,x.d*1e6,x.offset*1000,x.n,x.time.mean*1000,x.time.sd*1000,x.half.mean*1e6,x.half.sd*1e6,x.nh.mean,x.nh.sd,x.quarter.mean*1e6,x.quarter.sd*1e6,x.nq.mean,x.nq.sd,x.wall*100,x.outlet*100,x.inlet*100,x.censored*100,r.Gref]);data=[head,...records].map(row=>row.map(x=>'"'+String(x).replaceAll('"','""')+'"').join(',')).join('\r\n');type='text/csv';}
 const a=document.createElement('a'),url=URL.createObjectURL(new Blob([data],{type}));a.href=url;a.download='eductor-exposure-'+g.key+'.'+kind;a.click();setTimeout(()=>URL.revokeObjectURL(url),1000);
}
$('geometry').innerHTML=[['Earlier assemblies',g=>g.key.startsWith('rev')],['20Nov STEP cases',g=>g.key.startsWith('nov')]].map(([label,filter])=>`<optgroup label="${label}">`+GEOMETRIES.filter(filter).map(g=>`<option value="${g.key}">${g.name}</option>`).join('')+'</optgroup>').join('');$('geometry').value='revC';
if($('coverage'))$('coverage').innerHTML=AUDIT.map(row=>`<tr><td>${row.name}</td>${row.states.map(s=>`<td class="${s.status}" title="${s.reason}">${s.status==='screened'?'Screened':s.status==='sensitive'?'Sensitive':s.status==='compressible'?'Compressible estimate':'Unavailable'}${s.status==='unavailable'?'<br><small>'+s.reason+'</small>':''}</td>`).join('')}</tr>`).join('');
for(const id of ['p01','valveDrop','feederDrop','extension','geometry','flow-mode','p0','ps','po','feed','temp','small','large','density','air','bleed','chamber','release','offset','intensity','lengthScale','spread','n','dt','horizon','seed','loss','alpha','smallShare','background','collisionModel','reference','noise','sampling','phase-model','wallModel','kernelConvention','adaptiveDt','antithetic','commonNoise','kernelRange','populationBalance'])$(id).addEventListener('input',markDirty);
$('run').onclick=run;$('csv').onclick=()=>download('csv');$('json').onclick=()=>download('json');$('time-control').oninput=()=>{stop();drawTime();};$('play').onclick=play;
$('vacuum').onclick=()=>{$('ps').value=$('po').value=.20265;markDirty();};
$('slide-preset').onclick=()=>{const defaults={geometry:'revC','flow-mode':'auto',p01:3,valveDrop:0,feederDrop:0,extension:0,p0:3,ps:1,po:1,feed:6,temp:20,small:.3,large:10,density:2200,release:30,offset:.875,intensity:.18,lengthScale:.12,spread:.12,reference:'slide',n:128,alpha:1,smallShare:50,background:'mono',collisionModel:'hybrid',dt:.25,horizon:6,seed:42,loss:1};for(const [k,v] of Object.entries(defaults))k==='feed'?PowderRate.write($(k),v):$(k).value=v;$('noise').checked=true;$('phase-model').value='stable';$('sampling').value='fixed';$('adaptiveDt').checked=$('antithetic').checked=false;$('wallModel').value='reduced';$('kernelConvention').value='zhou2001';syncFlowControls();syncCollisionControls();run();};
// Deep links from a STEP case preserve physical absolute pressures.
const linkedSettings=new URLSearchParams(location.search),linkedGeometry=linkedSettings.get('geometry');
if(GEOMETRIES.some(g=>g.key===linkedGeometry)){
 $('geometry').value=linkedGeometry;
 for(const id of ['p0','ps','temp','loss','p01','valveDrop','feederDrop','extension','small','large','density','alpha','offset','seed','smallShare'])if(linkedSettings.has(id)){const v=Number(linkedSettings.get(id));if(Number.isFinite(v)&&v>=Number($(id).min)&&(!$(id).max||v<=Number($(id).max)))$(id).value=v;}
 for(const id of ['dt','background','collisionModel','phase-model','sampling','reference','wallModel','flow-mode'])if(linkedSettings.has(id)){const value=linkedSettings.get(id);if([...$(id).options].some(o=>o.value===value))$(id).value=value;}
 for(const [id,min,max] of [['horizon',6,100],['n',4,4096]])if(linkedSettings.has(id)){const v=Number(linkedSettings.get(id));if(Number.isFinite(v)&&v>=min&&v<=max&&(id!=='n'||Number.isInteger(v))){if(![...$(id).options].some(o=>o.value===String(v)))$(id).add(new Option(String(v),String(v)));$(id).value=String(v);}}
 for(const id of ['commonNoise','antithetic','adaptiveDt'])if(['true','false'].includes(linkedSettings.get(id)))$(id).checked=linkedSettings.get(id)==='true';
 const feedG=linkedSettings.has('feed_g_s')?Number(linkedSettings.get('feed_g_s')):linkedSettings.has('feed')?PowderRate.fromModel(Number(linkedSettings.get('feed'))):null;
 if(feedG!==null&&Number.isFinite(feedG)&&feedG>=0&&feedG<=Number($('feed').max))PowderRate.write($('feed'),PowderRate.toModel(feedG));
 $('po').value=$('ps').value;
 syncCollisionControls();
}
run();

window.addEventListener('pagehide',()=>ExposureCompute.cancel());

// Display the boundary of the same main-bore / representative-branch union as
// Exposure.geometry. Units here are mm. This is a visual mesh, not a flow solve.
const ExposureSurface=(()=>{
 function domain(g,release){
  const a=g.angle*Math.PI/180,sn=Math.sin(a),cs=Math.cos(a),cap=release+8,rb=g.db/2;
  function radius(x){x=Math.max(0,Math.min(g.end,x));const p=g.profile;let lo=0,hi=p.length-1;while(hi-lo>1){const m=(lo+hi)>>1;if(p[m][0]<x)lo=m;else hi=m;}const f=Math.max(0,Math.min(1,(x-p[lo][0])/(p[hi][0]-p[lo][0]||1)));return p[lo][1]+f*(p[hi][1]-p[lo][1]);}
  function field(x,y,z){const s=-cs*x+sn*z,c=sn*x+cs*z;return Math.max(Math.min(x,g.end-x,radius(x)-Math.hypot(y,z)),Math.min(s,cap-s,rb-Math.hypot(c,y)));}
  const r=Math.max(...g.profile.map(p=>p[1]));
  return {field,radius,sn,cs,cap,rb,bounds:[[-cs*cap-sn*rb-1,g.end+1],[-Math.max(r,rb)-1,Math.max(r,rb)+1],[-Math.max(r,cs*rb)-1,Math.max(r,sn*cap+cs*rb)+1]]};
 }
 function mesh(g,release,cut=true){
  const d=domain(g,release),step=Math.max(.38,Math.min(.72,g.db/28)),bounds=d.bounds.map(b=>[...b]);
  // The cut plane is exactly y=0, so trajectories remain in the open section.
  if(cut)bounds[1][0]=0;
  const n=bounds.map(b=>Math.ceil((b[1]-b[0])/step)),h=bounds.map((b,i)=>(b[1]-b[0])/n[i]),stride=n[0]+1,plane=stride*(n[1]+1),values=new Float32Array(plane*(n[2]+1));
  for(let k=0;k<=n[2];k++)for(let j=0;j<=n[1];j++)for(let i=0;i<=n[0];i++)values[k*plane+j*stride+i]=d.field(bounds[0][0]+i*h[0],bounds[1][0]+j*h[1],bounds[2][0]+k*h[2]);
  const offsets=[[0,0,0],[1,0,0],[1,1,0],[0,1,0],[0,0,1],[1,0,1],[1,1,1],[0,1,1]],tets=[[0,5,1,6],[0,1,2,6],[0,2,3,6],[0,3,7,6],[0,7,4,6],[0,4,5,6]],edges=[[0,1],[0,2],[0,3],[1,2],[1,3],[2,3]],vertices=[],normals=[],section=[];
  const gradient=p=>{const e=.025,v=[(d.field(p[0]+e,p[1],p[2])-d.field(p[0]-e,p[1],p[2])),(d.field(p[0],p[1]+e,p[2])-d.field(p[0],p[1]-e,p[2])),(d.field(p[0],p[1],p[2]+e)-d.field(p[0],p[1],p[2]-e))],l=Math.hypot(...v)||1;return v.map(x=>-x/l);};
  function triangle(a,b,c){const v=[a,b,c],center=a.map((x,i)=>(x+b[i]+c[i])/3),normal=gradient(center),u=b.map((x,i)=>x-a[i]),w=c.map((x,i)=>x-a[i]),cross=[u[1]*w[2]-u[2]*w[1],u[2]*w[0]-u[0]*w[2],u[0]*w[1]-u[1]*w[0]];if(Math.hypot(...cross)<1e-10)return;if(cross.reduce((s,x,i)=>s+x*normal[i],0)<0)v.reverse();for(const p of v){vertices.push(...p);normals.push(...gradient(p));}if(cut){const on=v.filter(p=>Math.abs(p[1])<1e-7);if(on.length===2)section.push(...on[0],...on[1]);}}
  for(let k=0;k<n[2];k++)for(let j=0;j<n[1];j++)for(let i=0;i<n[0];i++){
   const vs=offsets.map(o=>values[(k+o[2])*plane+(j+o[1])*stride+i+o[0]]);if(vs.every(v=>v>=0)||vs.every(v=>v<0))continue;
   const ps=offsets.map(o=>[bounds[0][0]+(i+o[0])*h[0],bounds[1][0]+(j+o[1])*h[1],bounds[2][0]+(k+o[2])*h[2]]);
   for(const tet of tets){const pts=[];for(const [a,b] of edges){const ia=tet[a],ib=tet[b];if((vs[ia]>=0)===(vs[ib]>=0))continue;let lo=0,hi=1;for(let q=0;q<11;q++){const f=(lo+hi)/2,p=ps[ia].map((x,k)=>x+f*(ps[ib][k]-x));if((d.field(...p)>=0)===(vs[ia]>=0))lo=f;else hi=f;}const f=(lo+hi)/2;pts.push(ps[ia].map((x,q)=>x+f*(ps[ib][q]-x)));}if(pts.length===3)triangle(...pts);else if(pts.length===4){
     // Sort the convex tetrahedron section around its centroid before triangulation.
     const c=[0,1,2].map(q=>pts.reduce((s,p)=>s+p[q],0)/4),normal=gradient(c),u=pts[0].map((x,q)=>x-c[q]),v=[normal[1]*u[2]-normal[2]*u[1],normal[2]*u[0]-normal[0]*u[2],normal[0]*u[1]-normal[1]*u[0]],angle=p=>Math.atan2(p.reduce((s,x,q)=>s+(x-c[q])*v[q],0),p.reduce((s,x,q)=>s+(x-c[q])*u[q],0));pts.sort((a,b)=>angle(a)-angle(b));triangle(pts[0],pts[1],pts[2]);triangle(pts[0],pts[2],pts[3]);
   }}
  }
  return {positions:new Float32Array(vertices),normals:new Float32Array(normals),section:new Float32Array(section),step};
 }
 return {domain,mesh};
})();
if(typeof module!=='undefined')module.exports=ExposureSurface;

const ExposureViewer=(()=>{
 const cad=CAD_MESHES;let renderer,scene,camera,orbit,root,paths,markers,geom,currentRun,mode='envelope',view='iso',displayRows=[];const surfaceCache=new Map();
 const host=()=>document.getElementById('three-view'),by=id=>document.getElementById(id);
 function paint(){if(renderer&&host().clientWidth&&host().clientHeight)renderer.render(scene,camera);}
 function dispose(group){if(!group)return;group.traverse(o=>{o.geometry?.dispose();if(o.material)(Array.isArray(o.material)?o.material:[o.material]).forEach(m=>{m.map?.dispose();m.dispose();});});scene.remove(group);}
 function shell(){return new THREE.MeshStandardMaterial({color:0x537d9d,metalness:.05,roughness:.72,transparent:by('view-ghost').checked,opacity:by('view-ghost').checked?.28:1,side:THREE.DoubleSide,depthWrite:!by('view-ghost').checked,clippingPlanes:by('view-cut').checked?[new THREE.Plane(new THREE.Vector3(...(geom?.key==='nov5'?[0,1,1]:['nov6','nov7'].includes(geom?.key)?[0,1,-1]:[0,1,0])).normalize(),0)]:[]});}
 function label(text,position){const canvas=document.createElement('canvas');canvas.width=384;canvas.height=64;const ctx=canvas.getContext('2d');ctx.font='600 46px system-ui';ctx.textAlign='center';ctx.fillStyle='#324f63';ctx.fillText(text,192,48);const map=new THREE.CanvasTexture(canvas),sprite=new THREE.Sprite(new THREE.SpriteMaterial({map,transparent:true,depthTest:false}));sprite.position.set(...position);sprite.scale.set(30,5,1);root.add(sprite);}

 function frame(cameraMode=view){if(!root)return;view=cameraMode;const box=new THREE.Box3();root.traverse(o=>{if(o.isMesh&&o.geometry){o.updateWorldMatrix(true,false);box.union(new THREE.Box3().setFromObject(o));}});box.expandByScalar(5);const center=box.getCenter(new THREE.Vector3()),size=box.getSize(new THREE.Vector3()),span=Math.max(size.x,size.y,size.z,40);orbit.target.copy(center);camera.up.set(0,0,1);const dir=(cameraMode==='side'?new THREE.Vector3(0,-1,0):cameraMode==='end'?new THREE.Vector3(1,0,0):new THREE.Vector3(...(['nov5','nov6','nov7'].includes(geom.key)?[-.53,-1,.80]:geom.key.startsWith('nov')?[-.52,-1,.55]:[.12,-1.8,.42]))).normalize(),right=new THREE.Vector3().crossVectors(new THREE.Vector3(0,0,1),dir).normalize(),up=new THREE.Vector3().crossVectors(dir,right).normalize(),tan=Math.tan(camera.fov*Math.PI/360);let dist=0;for(const x of [-1,1])for(const y of [-1,1])for(const z of [-1,1]){const corner=new THREE.Vector3(x*size.x/2,y*size.y/2,z*size.z/2);dist=Math.max(dist,corner.dot(dir)+Math.abs(corner.dot(right))/(tan*camera.aspect),corner.dot(dir)+Math.abs(corner.dot(up))/tan);}dist*=1.15;camera.position.copy(center).addScaledVector(dir,dist);camera.near=.1;camera.far=span*30;camera.updateProjectionMatrix();orbit.update();paint();}
 function build(){if(!renderer||!geom)return;dispose(root);dispose(paths);dispose(markers);root=new THREE.Group();paths=new THREE.Group();markers=new THREE.Group();const asset=cad.find(g=>g.key===geom.key);mode=asset?by('view-surface').value:'envelope';by('view-surface').disabled=!asset;by('view-surface').value=mode;
  if(mode==='cad'){
   const decode=s=>Uint8Array.from(atob(s),c=>c.charCodeAt(0)).buffer,q=new Uint16Array(decode(asset.mesh.v)),pos=new Float32Array(q.length);for(let i=0;i<q.length;i++)pos[i]=asset.mesh.lo[i%3]+q[i]/65535*asset.mesh.span-(i%3===0?(geom.xn||0):0);const geo=new THREE.BufferGeometry();geo.setAttribute('position',new THREE.BufferAttribute(pos,3));geo.setIndex(new THREE.BufferAttribute(new Uint32Array(decode(asset.mesh.t)),1));geo.computeVertexNormals();const mat=shell();mat.opacity=by('view-ghost').checked?.28:1;mat.transparent=by('view-ghost').checked;mat.depthWrite=!mat.transparent;root.add(new THREE.Mesh(geo,mat));by('viewer-provenance').textContent='STEP mesh · '+geom.dn+' mm jet · '+geom.d+' mm bore · '+geom.nports+' suction port'+(geom.nports>1?'s':'');by('viewer-note').textContent='Complete STEP assembly. 3D paths follow extracted passage centrelines; junctions are approximate. Use Ghost to inspect all branches. × marks first contact, not deposition.';
  }else{
   if(currentRun?.r.field.geo.dimension===3){
    const domain=currentRun.r.field.geo,mat=shell();mat.color.setHex(0x9bbbd0);mat.transparent=true;mat.opacity=.25;mat.depthWrite=false;
    const profile=geom.profile.map(q=>new THREE.Vector2(q[1],q[0])),main=new THREE.Mesh(new THREE.LatheGeometry(profile,48),mat);main.rotation.z=-Math.PI/2;root.add(main);
    for(const branch of domain.branches){const segments=branch.segments,points=[segments[0].a,...segments.map(s=>s.b)].map(p=>new THREE.Vector3(p[0]*1000,p[2]*1000,p[1]*1000)),lengths=[0];for(let i=1;i<points.length;i++)lengths.push(lengths.at(-1)+points[i].distanceTo(points[i-1]));const curve=new THREE.Curve();curve.getPoint=t=>{const d=t*lengths.at(-1);let j=lengths.findIndex(v=>v>=d);j=Math.max(1,j);return points[j-1].clone().lerp(points[j],(d-lengths[j-1])/(lengths[j]-lengths[j-1]));};root.add(new THREE.Mesh(new THREE.TubeGeometry(curve,160,branch.radius*1000,12,false),mat.clone()));}
   }else{
   const release=currentRun?.p.release||30,key=[geom.key,geom.end,release,by('view-cut').checked].join(':');
   if(!surfaceCache.has(key)){if(surfaceCache.size>6)surfaceCache.delete(surfaceCache.keys().next().value);surfaceCache.set(key,ExposureSurface.mesh(geom,release,by('view-cut').checked));}
   const surface=surfaceCache.get(key),geo=new THREE.BufferGeometry();geo.setAttribute('position',new THREE.BufferAttribute(surface.positions,3));geo.setAttribute('normal',new THREE.BufferAttribute(surface.normals,3));
   const mat=shell();mat.clippingPlanes=[];mat.color.setHex(0x9bbbd0);mat.opacity=by('view-ghost').checked?.18:1;mat.transparent=by('view-ghost').checked;mat.depthWrite=!mat.transparent;root.add(new THREE.Mesh(geo,mat));
   // One union boundary, with no cylinder walls inside the junction.
   if(surface.section.length){const edge=new THREE.BufferGeometry();edge.setAttribute('position',new THREE.BufferAttribute(surface.section,3));root.add(new THREE.LineSegments(edge,new THREE.LineBasicMaterial({color:0x365c76})));}
   const d=ExposureSurface.domain(geom,release),cap=d.cap;
   root.add(new THREE.ArrowHelper(new THREE.Vector3(1,0,0),new THREE.Vector3(-18,0,0),15,0x64758b,2.4,1.3));
   root.add(new THREE.ArrowHelper(new THREE.Vector3(1,0,0),new THREE.Vector3(geom.end+3,0,0),9,0x08678b,2.4,1.3));
   label('Motive air',[-18,0,-5.5]);label('Outlet',[geom.end+5,0,-geom.profile.at(-1)[1]-4]);label('Powder + air',[-d.cs*cap,0,d.sn*cap+d.cs*d.rb+6]);

   }
   by('viewer-provenance').textContent='Particle domain · '+(geom.nports>1?geom.nports+' branches':'connected passage');by('viewer-note').textContent='Solid: small particles · dotted: large · ×: first wall contact. Paths are 3D histories through the reduced passage domain.'+(geom.nports>1?' Full port geometry is in STEP view.':'');

  }
   displayRows=currentRun?currentRun.r.rows.flatMap((row,i)=>(row.paths?.length?row.paths:[{points:row.path,fate:row.pathFate}]).map(path=>({...row,path:path.points,pathFate:path.fate,colorIndex:i%3}))):[];displayRows.forEach((row,i)=>{const pts=row.path.map(p=>new THREE.Vector3(p[1]*1000,(p[4]||0)*1000,p[2]*1000)),geo=new THREE.BufferGeometry().setFromPoints(pts),mat=row.species?new THREE.LineDashedMaterial({color:COLORS[row.colorIndex],dashSize:.55,gapSize:.55}):new THREE.LineBasicMaterial({color:COLORS[row.colorIndex]});const line=new THREE.Line(geo,mat);line.computeLineDistances();paths.add(line);if(!row.species){const curve=new THREE.Curve();curve.getPoint=t=>{const q=t*(pts.length-1),j=Math.min(pts.length-2,Math.floor(q));return pts[j].clone().lerp(pts[j+1],q-j);};paths.add(new THREE.Mesh(new THREE.TubeGeometry(curve,Math.max(2,pts.length-1),.12,5,false),new THREE.MeshBasicMaterial({color:COLORS[row.colorIndex]})));}
    const sphere=new THREE.Mesh(new THREE.SphereGeometry(row.species?.48:.36,12,8),new THREE.MeshBasicMaterial({color:row.species?0xffffff:COLORS[row.colorIndex]}));markers.add(sphere);
    const start=new THREE.Mesh(new THREE.SphereGeometry(.48,12,8),new THREE.MeshBasicMaterial({color:COLORS[row.colorIndex]}));start.position.copy(pts[0]);paths.add(start);
    // Show first-contact endpoints rather than leaving unexplained broken paths.
    if(row.pathFate==='wall'){const p=pts.at(-1),end=new THREE.Group(),material=new THREE.LineBasicMaterial({color:COLORS[row.colorIndex]});for(const sign of [-1,1])end.add(new THREE.Line(new THREE.BufferGeometry().setFromPoints([new THREE.Vector3(-.65,0,sign*-.65),new THREE.Vector3(.65,0,sign*.65)]),material));end.position.copy(p);paths.add(end);}
   });
  scene.add(root,paths,markers);frame();time(+(by('time-control')?.value||0)/1000);
 }
 function time(t){if(!currentRun||!markers)return;displayRows.forEach((row,i)=>{const marker=markers.children[i];if(!marker)return;const ps=row.path;marker.visible=t<=ps.at(-1)[0];if(!marker.visible)return;let k=ps.findIndex(v=>v[0]>=t);if(k<0){marker.visible=false;return;}const a=ps[Math.max(0,k-1)],b=ps[k],f=(t-a[0])/Math.max(b[0]-a[0],1e-12);marker.position.set((a[1]+f*(b[1]-a[1]))*1000,((a[4]||0)+f*((b[4]||0)-(a[4]||0)))*1000,(a[2]+f*(b[2]-a[2]))*1000);});paint();}
 function update(run,g){const changed=geom?.key!==(g||run.g).key;currentRun=run;geom=g||run.g;if(changed){by('view-surface').value=cad.some(a=>a.key===geom.key)?'cad':'envelope';by('view-ghost').checked=false;by('view-cut').checked=!cad.some(a=>a.key===geom.key);}if(!host().hidden)build();}
 function resize(){if(!renderer)return;const w=host().clientWidth,h=host().clientHeight;if(!w||!h)return;renderer.setSize(w,h,false);camera.aspect=w/h;camera.updateProjectionMatrix();if(root)frame(view);else paint();}
 function init(){try{renderer=new THREE.WebGLRenderer({antialias:true,alpha:true});renderer.localClippingEnabled=true;renderer.setPixelRatio(Math.min(devicePixelRatio,2));scene=new THREE.Scene();camera=new THREE.PerspectiveCamera(35,1,.1,10000);scene.add(new THREE.HemisphereLight(0xffffff,0x66869b,.85));const light=new THREE.DirectionalLight(0xffffff,.8);light.position.set(-60,-120,150);scene.add(light);host().appendChild(renderer.domElement);renderer.domElement.setAttribute('aria-label','Interactive 3D eductor geometry');renderer.domElement.setAttribute('role','img');orbit=new THREE.OrbitControls(camera,renderer.domElement);orbit.enableDamping=false;orbit.addEventListener('change',paint);new ResizeObserver(resize).observe(host());for(const id of ['view-cut','view-ghost'])by(id).addEventListener('change',build);by('view-surface').addEventListener('change',()=>{by('view-ghost').checked=false;build();});const showPaths=()=>{if(mode==='cad'&&currentRun){by('view-ghost').checked=true;by('view-cut').checked=false;build();}};by('play').addEventListener('click',showPaths);by('time-control').addEventListener('input',showPaths);for(const mode of ['iso','side','end'])by('camera-'+mode).onclick=()=>frame(mode);resize();}catch(e){by('viewer-provenance').textContent='3D unavailable on this device. Use the 2D section.';by('viewer-note').textContent=e.message;renderer=null;}}
 function capture(){if(!renderer)throw Error('3D image is unavailable. Switch to Blender or 2D and keep again.');paint();return renderer.domElement.toDataURL('image/png');}
 return {init,update,time,resize,capture};
})();

const BLENDER_PROJECTIONS=[
  {
    "key": "revA",
    "width": 960,
    "height": 640,
    "matrix": [
      [
        8.215341567993164,
        0.8215427398681641,
        0.0,
        364.09071922302246
      ],
      [
        0.27400970458984375,
        -2.7400779724121094,
        -7.783527374267578,
        422.9713439941406
      ]
    ],
    "suctionAnchor": [
      -21.213203435596427,
      0,
      21.213203435596423
    ],
    "kind": "Analytical particle domain",
    "blender": "4.5.0"
  },
  {
    "key": "revB1",
    "width": 960,
    "height": 640,
    "matrix": [
      [
        6.181068420410156,
        0.6180953979492188,
        0.0,
        187.51097202301025
      ],
      [
        0.206146240234375,
        -2.061595916748047,
        -5.856227874755859,
        388.67847442626953
      ]
    ],
    "suctionAnchor": [
      -1.83697019872103e-15,
      0,
      30
    ],
    "kind": "Analytical particle domain",
    "blender": "4.5.0"
  },
  {
    "key": "revB2",
    "width": 960,
    "height": 640,
    "matrix": [
      [
        8.215341567993164,
        0.8215427398681641,
        0.0,
        364.09071922302246
      ],
      [
        0.27400970458984375,
        -2.7400779724121094,
        -7.783527374267578,
        422.9713439941406
      ]
    ],
    "suctionAnchor": [
      -21.213203435596427,
      0,
      21.213203435596423
    ],
    "kind": "Analytical particle domain",
    "blender": "4.5.0"
  },
  {
    "key": "revC",
    "width": 960,
    "height": 640,
    "matrix": [
      [
        9.846482276916504,
        0.984649658203125,
        0.0,
        157.96043872833252
      ],
      [
        0.32840728759765625,
        -3.284130096435547,
        -9.328956604003906,
        452.68930435180664
      ]
    ],
    "suctionAnchor": [
      -1.83697019872103e-15,
      0,
      30
    ],
    "kind": "Analytical particle domain",
    "blender": "4.5.0"
  },
  {
    "key": "nov0",
    "width": 960,
    "height": 640,
    "matrix": [
      [
        4.55841064453125,
        -2.3703575134277344,
        0.0,
        456.4072895050049
      ],
      [
        -1.0395050048828125,
        -1.999053955078125,
        -4.617443084716797,
        359.0880584716797
      ]
    ],
    "suctionAnchor": [
      11.301700000000004,
      -8.938606635269468e-16,
      38.5
    ],
    "kind": "Complete STEP assembly \u00b7 exterior",
    "blender": "4.5.0",
    "cameraDirection": [
      -0.41462242603302,
      -0.7973508238792419,
      0.4385429620742798
    ],
    "cutaway": false
  },
  {
    "key": "nov1",
    "width": 960,
    "height": 640,
    "matrix": [
      [
        4.55841064453125,
        -2.3703575134277344,
        0.0,
        456.4072895050049
      ],
      [
        -1.0395050048828125,
        -1.999053955078125,
        -4.617443084716797,
        359.0880584716797
      ]
    ],
    "suctionAnchor": [
      11.301700000000004,
      -8.938606635269468e-16,
      38.5
    ],
    "kind": "Complete STEP assembly \u00b7 exterior",
    "blender": "4.5.0",
    "cameraDirection": [
      -0.41462242603302,
      -0.7973508238792419,
      0.4385429620742798
    ],
    "cutaway": false
  },
  {
    "key": "nov2",
    "width": 960,
    "height": 640,
    "matrix": [
      [
        7.280187606811523,
        -3.78570556640625,
        0.0,
        476.15572929382324
      ],
      [
        -1.6601943969726562,
        -3.1926727294921875,
        -7.374515533447266,
        374.70951080322266
      ]
    ],
    "suctionAnchor": [
      1.8560000000000016,
      -8.938606635269468e-16,
      38.5
    ],
    "kind": "Complete STEP assembly \u00b7 exterior",
    "blender": "4.5.0",
    "cameraDirection": [
      -0.41462242603302,
      -0.7973508238792419,
      0.4385429620742798
    ],
    "cutaway": false
  },
  {
    "key": "nov3",
    "width": 960,
    "height": 640,
    "matrix": [
      [
        7.280187606811523,
        -3.78570556640625,
        0.0,
        476.15572929382324
      ],
      [
        -1.6601943969726562,
        -3.1926727294921875,
        -7.374515533447266,
        374.70951080322266
      ]
    ],
    "suctionAnchor": [
      1.8560000000000016,
      -8.938606635269468e-16,
      38.5
    ],
    "kind": "Complete STEP assembly \u00b7 exterior",
    "blender": "4.5.0",
    "cameraDirection": [
      -0.41462242603302,
      -0.7973508238792419,
      0.4385429620742798
    ],
    "cutaway": false
  },
  {
    "key": "nov4",
    "width": 960,
    "height": 640,
    "matrix": [
      [
        4.730072021484375,
        -2.4596786499023438,
        0.0,
        660.7865524291992
      ],
      [
        -1.0786628723144531,
        -2.0743560791015625,
        -4.791374206542969,
        384.10400390625
      ]
    ],
    "suctionAnchor": [
      -27.782294321497375,
      -8.93860663526947e-16,
      50.34600282048214
    ],
    "kind": "Complete STEP assembly \u00b7 exterior",
    "blender": "4.5.0",
    "cameraDirection": [
      -0.41462242603302,
      -0.7973508238792419,
      0.4385429620742798
    ],
    "cutaway": false
  },
  {
    "key": "nov5",
    "width": 960,
    "height": 640,
    "matrix": [
      [
        4.198894500732422,
        -2.2254180908203125,
        0.0,
        562.5022315979004
      ],
      [
        -1.2845230102539062,
        -2.4236679077148438,
        -3.880615234375,
        188.797607421875
      ]
    ],
    "suctionAnchor": [
      -91.08927522873023,
      -65.56014756989518,
      2.60140027202554e-15
    ],
    "kind": "Complete STEP assembly \u00b7 exterior",
    "blender": "4.5.0",
    "cameraDirection": [
      -0.3824049234390259,
      -0.7215187549591064,
      0.5772150158882141
    ],
    "cutaway": false
  },
  {
    "key": "nov6",
    "width": 960,
    "height": 640,
    "matrix": [
      [
        3.451995849609375,
        -1.8296241760253906,
        0.0,
        545.0947952270508
      ],
      [
        -1.056060791015625,
        -1.9925689697265625,
        -3.1903457641601562,
        207.3047637939453
      ]
    ],
    "suctionAnchor": [
      -114.70979579551081,
      -85.30330085889939,
      -2.4477358023231794e-16
    ],
    "kind": "Complete STEP assembly \u00b7 exterior",
    "blender": "4.5.0",
    "cameraDirection": [
      -0.3824049234390259,
      -0.7215187549591064,
      0.5772150158882141
    ],
    "cutaway": false
  },
  {
    "key": "nov7",
    "width": 960,
    "height": 640,
    "matrix": [
      [
        3.562488555908203,
        -1.8881034851074219,
        0.0,
        596.645393371582
      ],
      [
        -1.0898590087890625,
        -2.0563125610351562,
        -3.2924270629882812,
        188.56277465820312
      ]
    ],
    "suctionAnchor": [
      -114.70979579551081,
      -85.30330085889939,
      -2.4477358023231794e-16
    ],
    "kind": "Complete STEP assembly \u00b7 exterior",
    "blender": "4.5.0",
    "cameraDirection": [
      -0.3824049234390259,
      -0.7215187549591064,
      0.5772150158882141
    ],
    "cutaway": false
  }
];
// Shear at first forward passage through the same twelve sections as dilution.
// Statistics condition on particles reaching a section; missing is not zero.
const ExposureStationChart=(()=>{
 const colors=['#008678','#08678b','#bf5015'],ink='#365568';
 const escape=s=>String(s).replaceAll('&','&amp;').replaceAll('<','&lt;').replaceAll('"','&quot;');
 const num=x=>Number.isFinite(x)?x.toLocaleString('en-US',{maximumFractionDigits:3}):'—';
 function svg({rows,Gref=1,width,height,active=2,title='Shear',id='shear',quantity=null,combined=false,overlay=null}){
  const coag=quantity!==null,values=row=>coag?row.stationCoagulation[quantity]:row.stationShear,scale=coag?1:Gref;
  const metric=coag?({hazard:'Cumulative collision dose Λ',probability:'First-sticking probability P',screened:'Within-screen subtotal Λ',brownian:'Brownian collision dose Λ',turbulent:'Turbulent collision dose Λ',unsupported:'Outside-screen dose Λ'}[quantity]):'G / Gref';
  const bounds=s=>coag?[s.p05,s.p95]:[Math.max(0,s.mean-s.sd),s.mean+s.sd];
  const W=Math.max(150,width),H=Math.max(130,height),left=W<240?38:50,right=W-20,
   dense=(right-left-20)/(values(rows[0]).length-1)<18,short=H<240,top=short?28:44,base=H-(W>=680?(short?65:84):dense?66:46),
   xx=i=>left+10+i*(right-left-20)/(values(rows[0]).length-1),yy=v=>base-v*(base-top),clip=id+'-clip',
   sections=values(rows[0]),labels=sections.map(s=>s.label.split(/ |-/));
  let out=`<defs><clipPath id="${clip}"><rect x="${left}" y="${top}" width="${right-left}" height="${base-top}"/></clipPath></defs><text x="${left}" y="17" font-size="12" font-weight="650" fill="#173c51">${escape(title)}</text>`;
  if(!short)out+=`<text x="${left}" y="32" font-size="${W<280?8:10}" fill="${ink}">${coag?'Mean · P05–P95':'First passage · mean ± SD'}</text>`;
  if(rows.some(row=>values(row).some(s=>s.n&&Math.max(s.mean,bounds(s)[1])/scale>1)))out+=`<text x="${right}" y="17" text-anchor="end" font-size="9" fill="#9a5f1a"><title>Mean or spread exceeds the fixed 0–1 axis; full values remain in JSON exports and tooltips.</title>↑ above 1</text>`;
  out+=`<rect x="${xx(active)-10}" y="${top}" width="20" height="${base-top}" fill="#b64f12" opacity=".07"/>`;
  for(let i=0;i<=4;i++){const v=i/4;out+=`<line x1="${left}" x2="${right}" y1="${yy(v)}" y2="${yy(v)}" stroke="#dce5eb"/><text x="${left-7}" y="${yy(v)+3}" text-anchor="end" font-size="9" fill="${ink}">${v.toFixed(2)}</text>`;}
  if(!coag)for(const q of [.25,.5])out+=`<line x1="${left}" x2="${right}" y1="${yy(q)}" y2="${yy(q)}" stroke="#8299a6" stroke-dasharray="${q===.5?'5 3':'2 3'}"/>`;
  out+=`<path d="M${left},${top}V${base}H${right}" fill="none" stroke="#8fa3af"/><g clip-path="url(#${clip})">`;
  if(overlay){const previous=overlay.rows.filter(row=>combined||row.species===rows[0].species),denominator=coag?1:overlay.Gref;previous.forEach((row,k)=>{let path='',last=false;values(row).forEach((s,i)=>{if(s.n){path+=(last?'L':'M')+xx(i)+','+yy(s.mean/denominator);last=true;}else last=false;});out+=`<path class="saved-experiment-curve" d="${path}" fill="none" stroke="${colors[k%3]}" stroke-width="2" stroke-dasharray="${row.species&&combined?'2 3 7 3':'7 4'}" opacity=".5"><title>Saved experiment · ${escape(overlay.name)} · ${num(row.d*1e6)} µm · release ${num(row.offset*1000)} mm</title></path>`;});}
  rows.forEach((row,k)=>{
   let path='',previous=false;const c=colors[k%3];
   values(row).forEach((s,i)=>{if(s.n){path+=(previous?'L':'M')+xx(i)+','+yy(s.mean/scale);previous=true;}else previous=false;});
   out+=`<path class="station-${coag?'coagulation':'shear'}-curve" d="${path}" fill="none" stroke="${c}" stroke-width="2" ${combined&&row.species?'stroke-dasharray="5 3"':''}/>`;
   values(row).forEach((s,i)=>{
    if(!s.n)return;const x=xx(i),[low,high]=bounds(s),lo=yy(low/scale),hi=yy(high/scale);
    if(s.n>1)out+=`<path d="M${x},${lo}V${hi}M${x-3},${lo}H${x+3}M${x-3},${hi}H${x+3}" fill="none" stroke="${c}" stroke-opacity=".5" stroke-width="1.4"/>`;
   });
  });out+='</g>';
  rows.forEach((row,k)=>values(row).forEach((s,i)=>{
   const c=colors[k%3],value=s.mean/scale,spread=coag?' · P05–P95 '+num(s.p05)+'–'+num(s.p95)+'; median '+num(s.median)+'; SE '+num(s.se):' ± '+num(s.sd/scale),description=`${s.id}. ${s.label}; ${num(row.d*1e6)} µm; release ${num(row.offset*1000)} mm: ${s.n?metric+' '+num(value)+spread+'; '+s.n+'/'+s.total+' histories reached':'not crossed by any of '+s.total+' histories'}`;
   out+=`<g role="button" tabindex="0" data-station="${i}" aria-label="${escape(description)}" aria-pressed="${i===active}"><title>${escape(description)}. Conditional on first forward crossing; wall stops and unfinished histories are excluded from later sections.</title>`;
   if(s.n)out+=`<circle class="station-${coag?'coagulation':'shear'}-point" data-value="${value}" data-count="${s.n}" data-total="${s.total}" cx="${xx(i)}" cy="${yy(Math.min(1,value))}" r="${i===active?4:3}" fill="${c}" stroke="white" stroke-width="1"/>`;
   else out+=`<text class="station-shear-missing" x="${xx(i)}" y="${base-5-k*9}" text-anchor="middle" font-size="10" fill="${c}">—</text>`;
   out+='</g>';
  }));
  sections.forEach((s,i)=>{const x=xx(i),selected=i===active,c=selected?'#b64f12':'#08678b',r=dense?7:9,badgeY=base+17+(dense&&i%2?17:0);
   out+=`<g role="button" tabindex="0" data-station="${i}" aria-label="${s.id}. ${escape(s.label)}" aria-pressed="${selected}"><title>${s.id}. ${escape(s.label)} · same section as dilution and Blender</title><circle class="station-badge" cx="${x}" cy="${badgeY}" r="${r}" fill="${selected?c:'white'}" stroke="${c}"/><text x="${x}" y="${badgeY+3}" text-anchor="middle" font-size="9" font-weight="700" fill="${selected?'white':c}">${s.id}</text>`;
   if(W>=680)labels[i].forEach((line,j)=>out+=`<text x="${x}" y="${base+(short?32:40)+j*(short?10:12)}" text-anchor="middle" font-size="${W<440?9:11}" fill="${ink}">${line}</text>`);
   out+='</g>';
  });
  out+=`<text x="${(left+right)/2}" y="${H-7}" text-anchor="middle" font-size="${W<280?8:10}" fill="${ink}">${W<680?'Locations 1–12 · select for details':'Location along flow · matches dilution'}</text>`;
  return out;
 }
 return {svg};
})();
if(typeof module!=='undefined')module.exports=ExposureStationChart;

// Bulk dilution at sections of the selected stage. Render annotations use the
// camera projection exported by Blender, in nozzle-relative millimetres.
const ExposureDilution=(()=>{
 const model=typeof Exposure!=='undefined'?Exposure:require('./model.js');
 let latest=null,active=2;
 const names=model.stationSections({peakX:0,geo:{end:1}}).map(s=>s.label);
 function stations(run){
  if(!run?.r?.ok)return [];
  const f=run.r.pressure,field=run.r.field,sections=model.stationSections(field),px=sections[2].x,end=sections.at(-1).x;
  const point=(label,air,pressure,temp,detail,x=null,powder=f.md)=>({label,air,powder,pressure,temp,detail,x,
   powderAir:air>0?powder/air:null,concentration:air>0?1000*pressure/(287.05*temp)*powder/air:null});
  const section=(label,x,detail)=>{const a=field.axial(x);return point(label,a.mass,a.pressure,a.temp,detail,x,a.powder??f.md);};
  const inlet=point(names[0],f.ms,f.inletP,f.Tin,'Powder-carrying suction stream, before this stage’s motive jet.');
  const points=sections.map((s,i)=>i===0?inlet:i===sections.length-1?point(s.label,f.mg,f.po,f.T,`Selected-stage outlet; x = ${(end*1000).toFixed(1)} mm.`,end):
   i===1?section(s.label,s.x,'Nozzle-tip section, x = 0 mm; CAD suction may enter further downstream.'):
   field.peak>0?section(s.label,s.x,`${i===2?'Modeled shear maximum':Number((100*(s.x-px)/(end-px)).toFixed(1))+'% of the distance from peak shear to outlet'}; x = ${(s.x*1000).toFixed(1)} mm.`):null);
  return points.map((p,i)=>p?{...p,id:i+1,ratioRelative:inlet.powderAir>0?p.powderAir/inlet.powderAir:null,
   concentrationRelative:inlet.concentration>0?p.concentration/inlet.concentration:null}:null);
 }
 function project(mapping,xyz){return mapping.matrix.map(row=>row[3]+xyz.reduce((s,v,i)=>s+v*row[i],0));}
 const by=id=>document.getElementById(id);
 const number=v=>!Number.isFinite(v)?'—':v!==0&&(Math.abs(v)<.01||Math.abs(v)>=10000)?v.toExponential(1):v.toLocaleString('en-US',{maximumFractionDigits:2});
 const compact=v=>Number.isFinite(v)?number(Number(v.toPrecision(3))):'—';
 const escape=s=>String(s).replaceAll('&','&amp;').replaceAll('<','&lt;').replaceAll('"','&quot;');
 function clear(){latest=null;by('lab-dilution-chart').innerHTML='<text x="50%" y="50%" text-anchor="middle" font-size="12" fill="#647785">No valid flow solution</text>';by('dilution-overlay').innerHTML='';by('dilution-station-key').innerHTML='';by('dilution-station-detail').textContent='No station values without a supported flow solution.';by('lab-dilution-selection').textContent='';}
 function render(run){latest=run;draw();}
 function select(index){if(!Number.isInteger(index)||index<0||index>=names.length)return;active=index;draw();window.exposureStationSelected?.();}
 function init(){
  for(const id of ['dilution-overlay','dilution-station-key','lab-dilution-plots','shear-history-panel','coagulation-history-panel']){
   const el=by(id),choose=e=>{const t=e.target.closest('[data-station]');if(t)select(Number(t.dataset.station));};el.addEventListener('click',choose);el.addEventListener('keydown',e=>{if(e.key==='Enter'||e.key===' '){e.preventDefault();choose(e);}});
  }
  by('dilution-quantity').addEventListener('change',draw);
  by('blender-view').addEventListener('load',drawOverlay);
  new ResizeObserver(drawOverlay).observe(by('dilution-overlay').parentElement);
 }
 function drawOverlay(){
  const target=by('dilution-overlay');if(!latest||!target)return;
  const mapping=BLENDER_PROJECTIONS.find(x=>x.key===latest.g.key);if(!mapping){target.innerHTML='';return;}
  const points=stations(latest),W=by('blender-view').clientWidth-10,H=by('blender-view').clientHeight-25,scale=Math.min(W/960,H/640),ui=Math.min(2.3,Math.max(1,1/Math.max(.2,scale)));
  const anchors=points.map((p,i)=>p?project(mapping,i===0?mapping.suctionAnchor:[p.x*1000,0,0]):null);
  const downstream=anchors.slice(1).filter(Boolean),laneY=Math.min(610-12*ui,Math.max(...downstream.map(p=>p[1]))+48*ui),a=downstream[0],z=downstream.at(-1),width=Math.min(900,Math.max(z[0]-a[0],23*ui*(points.length-2))),center=Math.min(940-width/2,Math.max(20+width/2,(a[0]+z[0])/2));
  let svg='';anchors.forEach((q,i)=>{if(!q)return;const p=points[i],selected=i===active,color=selected?'#b64f12':'#075d7a',bx=i===0?Math.max(20*ui,q[0]-26*ui):center-width/2+(i-1)*width/(points.length-2),byy=i===0?Math.max(25*ui,q[1]-30*ui):laneY,r=8*ui;
   svg+=`<g role="button" tabindex="0" data-station="${i}" aria-label="${i+1}. ${escape(p.label)}: ${number(p.powderAir)} kg powder per kg air" aria-pressed="${selected}"><title>${i+1}. ${escape(p.label)}. ${escape(p.detail)}</title><path d="M${q[0]},${q[1]}L${bx},${byy}" stroke="white" stroke-width="${3.5*ui}" fill="none"/><path d="M${q[0]},${q[1]}L${bx},${byy}" stroke="${color}" stroke-width="${1.25*ui}" fill="none"/><circle cx="${q[0]}" cy="${q[1]}" r="${3.4*ui}" fill="${color}" stroke="white" stroke-width="${1.5*ui}"/><circle class="station-badge" cx="${bx}" cy="${byy}" r="${r}" fill="${selected?color:'white'}" stroke="${color}" stroke-width="${1.5*ui}"/><text x="${bx}" y="${byy+3.8*ui}" text-anchor="middle" font-size="${11*ui}" font-weight="700" fill="${selected?'white':color}">${i+1}</text></g>`;
  });target.innerHTML=latest.g.extensionMM?'':svg;
  by('dilution-station-key').innerHTML=names.map((name,i)=>`<button type="button" data-station="${i}" aria-pressed="${i===active}" title="${i+1}. ${name}"><b>${i+1}</b>${name}</button>`).join('');
  const p=points[active];by('dilution-station-detail').innerHTML=p?`<strong>${p.id} · ${escape(p.label)}</strong> · ${number(p.powderAir)} kg/kg · ${number(p.concentration)} g/m³<span>${escape(p.detail)}${active===0&&latest.g.nports>1?' Common suction inlet; branch flows follow hydraulic resistance.':''}</span>`:'No modeled shear maximum for this station.';by('lab-dilution-selection').innerHTML=p?`<b>${p.id} · ${escape(p.label)}</b> · ${number(p.powderAir)} kg/kg · ${number(p.concentration)} g/m³`:'';
 }
 // A zero-based linear axis in the SAME units as the plotted labels.
 function axisScale(values){
  const peak=Math.max(0,...values.filter(Number.isFinite));
  if(!peak)return {max:1,ticks:[0,.25,.5,.75,1]};
  const raw=peak/4,power=10**Math.floor(Math.log10(raw)),step=[1,2,2.5,5,10].find(v=>v*power>=raw)*power;
  const count=Math.ceil(peak*1.08/step);
  return {max:count*step,ticks:Array.from({length:count+1},(_,i)=>i*step)};
 }
 function draw(){
  if(!latest)return;
  const previous=typeof ExposureOverlay!=='undefined'?ExposureOverlay.get():null;const points=stations(latest),quantity=by('dilution-quantity').value,
   concentration=quantity==='concentration',key=concentration?'concentration':'powderAir',
   title=concentration?'Powder concentration':'Powder / air',unit=concentration?'g/m³':'kg/kg',color=concentration?'#008678':'#08678b';
  const target=by('lab-dilution-chart'),W=Math.max(150,target.clientWidth||300),H=Math.max(150,target.clientHeight||300),
   left=W<240?38:50,right=W-20,dense=(right-left-20)/(points.length-1)<18,short=H<240,top=short?28:52,base=H-(W>=680?(short?65:84):dense?66:52),xx=i=>left+10+i*(right-left-20)/(points.length-1),
   scale=axisScale([...points,...(previous?.points||[])].map(p=>p?.[key])),yy=v=>base-v/scale.max*(base-top);
  target.setAttribute('viewBox',`0 0 ${W} ${H}`);
  target.setAttribute('aria-label',title+' in '+unit+' at twelve numbered flow locations; zero-based linear axis');
  let svg=`<text x="${left}" y="17" font-size="${W<240?10:12}" font-weight="650" fill="#173c51">${title} (${unit})</text>`;
  if(!short)svg+=`<text x="${left}" y="33" font-size="${W<240?9:10}" fill="#536b78">Lower = more dilute</text>`;
  for(const v of scale.ticks)svg+=`<line x1="${left}" x2="${right}" y1="${yy(v)}" y2="${yy(v)}" stroke="#dce5eb"/><text class="dilution-tick" x="${left-7}" y="${yy(v)+4}" text-anchor="end" font-size="10" fill="#536b78">${number(v)}</text>`;
  svg+=`<path d="M${left},${top}V${base}H${right}" fill="none" stroke="#8fa3af"/>`;
  let path='';points.forEach((p,i)=>{if(Number.isFinite(p?.[key]))path+=(i&&Number.isFinite(points[i-1]?.[key])?'L':'M')+xx(i)+','+yy(p[key]);});
  if(previous){let oldPath='';previous.points.forEach((p,i)=>{if(Number.isFinite(p?.[key]))oldPath+=(i&&Number.isFinite(previous.points[i-1]?.[key])?'L':'M')+xx(i)+','+yy(p[key]);});svg+=`<path class="saved-experiment-curve" d="${oldPath}" fill="none" stroke="${color}" stroke-width="2" stroke-dasharray="7 4" opacity=".5"><title>Saved experiment · ${escape(previous.name)}</title></path>`;previous.points.forEach((p,i)=>{if(Number.isFinite(p?.[key]))svg+=`<circle class="saved-experiment-point" cx="${xx(i)}" cy="${yy(p[key])}" r="3" fill="white" stroke="${color}" opacity=".65"><title>Saved ${p.id}. ${escape(p.label)}: ${number(p[key])} ${unit}</title></circle>`;});}
  svg+=`<path class="dilution-curve" d="${path}" fill="none" stroke="${color}" stroke-width="2.5"/>`;
  // Prioritize the selected and end values when all twelve labels cannot fit.
  // Every point still exposes its exact value through selection and its tooltip.
  const labelFont=W<240?9:11,labelBoxes=[],labels=new Map();
  const candidates=[active,0,points.length-1,...(W>=600?points.map((_,i)=>i):[])];
  for(const i of new Set(candidates)){const v=points[i]?.[key];if(!Number.isFinite(v))continue;const text=number(v),width=text.length*labelFont*.62+6,x=Math.max(left+width/2,Math.min(right-width/2,xx(i))),y=yy(v)-10,box={l:x-width/2,r:x+width/2,t:y-labelFont,b:y+3};if(labelBoxes.some(b=>box.l<b.r&&box.r>b.l&&box.t<b.b&&box.b>b.t))continue;labelBoxes.push(box);labels.set(i,{x,y,text});}
  points.forEach((p,i)=>{
   const x=xx(i),v=p?.[key],valid=Number.isFinite(v),y=valid?yy(v):base,selected=i===active,c=selected?'#b64f12':color,badgeY=base+19+(dense&&i%2?17:0);
   svg+=`<g role="button" tabindex="0" data-station="${i}" aria-label="${i+1}. ${escape(names[i])}: ${number(v)} ${unit}" aria-pressed="${selected}"><title>${i+1}. ${escape(names[i])}: ${number(v)} ${unit}. ${escape(p?.detail||'No current value')}</title><rect x="${x-10}" y="${top-8}" width="20" height="${base-top+37}" fill="${c}" opacity="${selected?.08:0}"/>`;
   if(valid)svg+=`<circle class="dilution-point" data-value="${v}" cx="${x}" cy="${y}" r="${selected?5:3.5}" fill="${c}" stroke="white" stroke-width="1.2"/>`;
   const label=labels.get(i);if(label)svg+=`<text class="dilution-value-label" x="${label.x}" y="${label.y}" text-anchor="middle" font-size="${labelFont}" font-weight="600" fill="${c}">${label.text}</text>`;
   svg+=`<circle cx="${x}" cy="${badgeY}" r="${dense?7:9}" fill="${selected?c:'white'}" stroke="${c}"/><text x="${x}" y="${badgeY+3}" text-anchor="middle" font-size="10" font-weight="700" fill="${selected?'white':c}">${i+1}</text>`;
   if(W>=680){const lines=names[i].replace('Suction inlet','Suction inlet').split(/ |-/);lines.forEach((line,j)=>svg+=`<text x="${x}" y="${base+(short?34:42)+j*(short?10:12)}" text-anchor="middle" font-size="10" fill="#536b78">${line}</text>`);}
   svg+='</g>';
  });
  svg+=`<text x="${(left+right)/2}" y="${H-8}" text-anchor="middle" font-size="${W<330?8:10}" fill="#536b78">${W<680?'Locations 1–12 · select for details':'Location along flow · spacing is schematic'}</text>`;
  if(points.every(p=>!p||p[key]===0))svg+=`<text x="${(left+right)/2}" y="${(top+base)/2}" text-anchor="middle" font-size="11" fill="#647785">No powder feed · all values are zero</text>`;
  target.innerHTML=svg;drawOverlay();
 }

 return {stations,project,axisScale,render,draw,clear,init,select,selectedIndex:()=>active};
})();
if(typeof module!=='undefined')module.exports=ExposureDilution;

const ExposureAnalytics=(()=>{
 let latest=null;
 const titles={hazard:'Total encounter dose Λ',probability:'First-sticking probability',screened:'Within-screen subtotal Λ',brownian:'Brownian encounter dose',turbulent:'Turbulent encounter dose',unsupported:'Outside-screen dose'};
 const metricInfo={t:['Tracked time · ms',1000],outletTime:['Outlet transit · ms',1000],half:['Time above ½ reference · µs',1e6],quarter:['Time above ¼ reference · µs',1e6],shearDose:['Accumulated shear ∫G dt',1],epsilonDose:['Accumulated dissipation · m²/s²',1],peak:['Peak shear · 10⁶ s⁻¹',1e-6],screened:['Within-screen subtotal',1],brownian:['Brownian dose',1],turbulent:['Turbulent dose',1],inertial:['Inertial increment over tracer dose',1],unsupported:['Outside-screen dose',1],maxCollisionRe:['Maximum collision-model slip Reynolds number',1],maxSt:['Maximum Stokes number Stη',1],maxRe:['Maximum particle Reynolds number',1],maxConcentration:['Maximum powder concentration · kg/m³',1],hazard:['Total encounter dose Λ',1],probability:['At least one sticking event · %',100]};
 const f=(v,n=3)=>v===null||!Number.isFinite(v)?'—':Math.abs(v)>0&&(Math.abs(v)<.001||Math.abs(v)>=1e5)?v.toExponential(2):fmt(v,n);
 function init(pages){
 pages.coagulation.innerHTML=`<div class="page-head"><div><h2>Brownian & turbulent collision exposure</h2><p>Integrate collision exposure along every history, until wall contact, exit or the tracking horizon.</p></div><select id="coag-quantity" aria-label="Coagulation quantity"><option value="hazard">Total encounter dose Λ</option><option value="probability">Probability of ≥1 sticking event</option><option value="screened">Within-screen subtotal Λ</option><option value="brownian">Brownian collision dose</option><option value="turbulent">Turbulent collision dose</option><option value="unsupported">Outside-screen dose · diagnostic</option></select></div><div class="coag-band-tools"><label>Band <select id="coag-band"><option value="spread">P05–P95 spread</option><option value="ci">95% CI of mean</option></select></label><label><input type="checkbox" id="coag-traces"> Show 12 example tracks per group</label><span>Axes fixed at 0–1; ↑ indicates values above view. Probability is a fraction.</span></div><div id="coag-summary" class="analytics-summary"></div><div class="analytics-charts"><section class="card"><svg id="coag-small" viewBox="0 0 600 270" role="img" aria-label="Small particle cumulative coagulation"></svg></section><section class="card"><svg id="coag-large" viewBox="0 0 600 270" role="img" aria-label="Large particle cumulative coagulation"></svg></section></div><div class="analytics-legend" id="coag-legend"></div><p id="coag-screen" class="notice"></p><section class="card analytics-table"><div class="table-wrap"><table><thead><tr><th>Size / offset</th><th>n</th><th>Mean Λ ± SD</th><th>Median Λ</th><th>P05–P95 Λ</th><th>95% CI · mean Λ</th><th>Relative SE</th><th>Max Λ</th><th>Mean sticking P</th><th>Within-screen share</th></tr></thead><tbody id="coag-table"></tbody></table></div></section><section class="card analytics-table"><h3>Mechanisms and applicability</h3><p class="note">Brownian + turbulent = total dose. Inertial increment is included in turbulence, never added a second time. Outside-screen dose is an extrapolation; the within-screen subtotal does not estimate the missing physics.</p><div class="table-wrap"><table><thead><tr><th>Size / offset</th><th>Brownian Λ</th><th>Turbulent Λ</th><th>Inertial increment</th><th>Outside screen</th><th>Max Stη</th><th>Max R/η</th><th>Max slip Re</th></tr></thead><tbody id="coag-mechanisms"></tbody></table></div></section><section class="card analytics-table"><h3>Sample-size stability · collision dose Λ</h3><p class="note">Same seeded histories at N/4, N/2 and N. A stable mean alone does not prove convergence; use the confidence interval and compare seeds.</p><div class="table-wrap"><table><thead><tr><th>Size / offset</th><th>N/4 mean</th><th>N/2 mean</th><th>N mean</th><th>Change N/2 → N</th></tr></thead><tbody id="coag-convergence"></tbody></table></div></section><details class="coag-methods"><summary>Equations, concentration and applicability</summary><p><b>Brownian:</b> Kᴮ = 4πR(Dᵢ+Dⱼ) / [R/(R+√(δᵢ²+δⱼ²)) + 4(Dᵢ+Dⱼ)/(R√(cᵢ²+cⱼ²))]. Fuchs transition interpolation includes Cunningham slip, local absolute pressure and temperature. D = kᴮTCc/(3πμd); c = √(8kᴮT/πm). It approaches continuum diffusion and free-molecular encounter limits.</p><p><b>Turbulent:</b> Kᵀ = √(8π) R² √(R²G²/15 + w²accel). WZW (Wang–Zhou–Wexler) supplies the inertial acceleration variance using local fluctuation variance, dissipation and both particle response times. The zero-inertia limit recovers Saffman–Turner. Shear and acceleration are combined as velocity variances, not added collision kernels. The baseline assumes g(R) = 1. Optional Zhou RDF is an upper-bound sensitivity away from its fitted contact radius. The extended option uses Abrahamson kinetic theory at large inertia or contact radius, with a separate decorrelation screen.</p><p><b>Along each history:</b> Λ = ∫Σⱼ(Kᴮ+Kᵀ)nⱼ dt and P = 1 − exp(−αΛ). Adding Brownian and turbulent encounter rates is a reduced independent-mechanism approximation. P is a first-event estimate in a frozen background, not coagulated mass or final size. Optional sectional aggregation separately calculates growth and depletion. Breakup, hydrodynamic deflection, electrostatics and differential settling collisions are omitted. α is the selected capture × sticking efficiency; 1 is the ideal limit.</p><p><b>Applicability:</b> conservative turbulence screening requires R/η ≤ 0.1, turbulent Mach ≤ 0.3, estimated slip Re ≤ 0.1 and particle/gas density ≥ 100. The partner slip is a local statistical estimate. Legacy mode additionally requires Stη ≤ 0.1. WZW includes inertia but remains a locally isotropic, homogeneous, heavy point-particle approximation. Passing these checks does not validate its use in a confined jet. Finite-size, nonlinear drag and compressible-turbulence intervals remain explicit extrapolations. Brownian dose remains included in the within-screen subtotal. That subtotal is neither a corrected total nor a rigorous lower bound.</p><p>Local powder concentration is the flux-conservative analytical scalar field. nⱼ = Cwⱼ/(ρparticle πdⱼ³/6); branch powder shares follow solved gas shares. Same-size mode means independent monodisperse scenarios; mixture mode uses fixed mass fractions. Histories end at first contact, exit or the tracking limit. Station averages condition on reaching each section; later sections can contain different histories. P05–P95 shows history spread; 95% mean intervals include sampling error only.</p><p><a href="https://www2.mmm.ucar.edu/people/wrfhelp/wrfvar/code/trunk/module_mosaic_coag.html" target="_blank" rel="noopener">Fuchs implementation in WRF/MOSAIC</a> · <a href="https://research.me.udel.edu/lwang/reprints/njp8_7_075016.pdf" target="_blank" rel="noopener">WZW formulation, equations 22–24</a> · <a href="https://authors.library.caltech.edu/records/mvczm-fey53" target="_blank" rel="noopener">Saffman & Turner</a></p></details>`;
 pages.statistics.innerHTML=`<div class="page-head"><div><h2>History statistics</h2><p id="stats-caption"></p></div><button id="samples-csv">Export every history · CSV</button></div><div class="statistics-tools"><label for="stats-metric">Distribution</label><select id="stats-metric">${Object.entries(metricInfo).map(([k,v])=>`<option value="${k}">${v[0]}</option>`).join('')}</select><span>More samples: change Histories in Experiment → Model.</span></div><section class="card distribution-card"><svg id="distribution-chart" viewBox="0 0 920 310" role="img" aria-label="Distribution by particle size and release offset"></svg><p class="note">Whiskers: P05–P95 · box: P25–P75 · line: median · dot: mean. Statistics are weighted by predicted branch powder flux. These show spread between histories, not uncertainty in the model.</p></section><section class="card analytics-table"><div class="table-wrap"><table><thead><tr><th>Size / offset</th><th>n</th><th>Mean</th><th>SD</th><th>SE</th><th>Approx. 95% CI of mean</th><th>P05</th><th>Median</th><th>P95</th><th>P99</th><th>CV · SD/mean</th></tr></thead><tbody id="stats-table"></tbody></table></div></section><section class="card fate-card"><h3>Where histories end</h3><div class="fate-legend"><span style="--fate:#bf5015">Wall contact</span><span style="--fate:#008678">Outlet</span><span style="--fate:#8061ad">Inlet escape</span><span style="--fate:#a2adb5">Censored</span><span style="--fate:#d32f2f">Numerical stop</span></div><div id="fate-bars"></div><p id="stats-caution" class="note"></p></section>`;
 for(const id of ['coag-band','coag-traces','coag-quantity'])$(id).onchange=()=>latest&&drawCoag();$('stats-metric').onchange=()=>latest&&drawStatistics();$('samples-csv').onclick=exportSamples;
 }
 function renderQuality({p,r}){
 let panel=$('numerical-quality');if(!panel){panel=document.createElement('section');panel.id='numerical-quality';panel.className='card analytics-table';$('page-statistics').append(panel);}
 panel.innerHTML='<h3>Numerical precision & outlet residence</h3><p class="note">Outlet statistics include successful exits only. A missing value means there were no outlet histories. Wilson intervals describe Monte Carlo precision, not physical-model accuracy. Export JSON includes code hashes, settings, seed and runtime.</p><div class="table-wrap"><table><thead><tr><th>Size / offset</th><th>Histories / effective n</th><th>Wall half-width</th><th>±2 pp target</th><th>Max step · µs</th><th>Outlet n</th><th>Outlet mean · ms</th><th>P10 / P50 / P90 · ms</th></tr></thead><tbody>'+r.rows.map(x=>{const t=x.distribution.outletTime;return `<tr><td>${sizeLabel(x.d)} / ${fmt(x.offset*1000,1)}</td><td>${x.n} / ${fmt(x.precision.effectiveN,0)}</td><td>±${fmt(100*x.precision.halfWidth,2)} pp</td><td>${p.noise?(x.precision.met?'Reached':'Not reached'):'Deterministic'}</td><td>${fmt(x.distribution.maxDt.max*1e6,4)}</td><td>${t.n}</td><td>${t.mean===null?'—':fmt(t.mean*1000,3)}</td><td>${[t.p10,t.median,t.p90].map(v=>v===null?'—':fmt(v*1000,3)).join(' / ')}</td></tr>`;}).join('')+'</tbody></table></div><p class="note">Code: '+r.manifest.codeSHA256.slice(0,16)+' · '+p.kernelConvention+' · '+(p.phaseModel||'stable')+'</p>';
 }
 function renderPopulation({p,r}){
  let panel=$('population-results');if(!panel){panel=document.createElement('section');panel.id='population-results';panel.className='card analytics-table';$('page-coagulation').append(panel);}
  if(!p.populationBalance){panel.innerHTML='<h3>Coagulation growth and depletion</h3><p>Enable “Evolve size distribution” in Experiment → Particles to calculate number depletion and the size distribution along each history. The plots above remain encounter-dose diagnostics.</p>';return;}
  const chart=species=>{const rows=r.rows.filter(row=>row.species===species),lines=[];let svg='<svg viewBox="0 0 560 270" role="img" aria-label="Particle number remaining by location"><text x="48" y="22" font-size="16" fill="#173748">'+sizeLabel(rows[0].d)+' · particle number remaining</text>';for(let i=0;i<5;i++){const y=43+i*42;svg+=`<path d="M48 ${y}H540" stroke="#dbe5eb"/><text x="40" y="${y+4}" text-anchor="end" font-size="11">${(1-i/4).toFixed(2)}</text>`;}
  rows.forEach((row,k)=>{let points=[];for(let i=0;i<r.rows[0].stationShear.length;i++){const samples=row.samples.filter(s=>s.growth?.stations[i]),q=Exposure.weightedStats(samples.map(s=>s.growth.stations[i].numberFraction),samples.map(s=>s.powderFluxWeight)),x=48+i*492/(r.rows[0].stationShear.length-1);if(!q.n){if(points.length)lines.push(`<polyline points="${points.join(' ')}" fill="none" stroke="${COLORS[k]}" stroke-width="2"/>`);points=[];continue;}const y=211-168*q.mean;points.push(x+','+y);svg+=`<circle cx="${x}" cy="${y}" r="3" fill="${COLORS[k]}"><title>${i+1}: ${fmt(q.mean,4)}; ${q.n} histories</title></circle>`;}if(points.length)lines.push(`<polyline points="${points.join(' ')}" fill="none" stroke="${COLORS[k]}" stroke-width="2"/>`);});for(let i=0;i<r.rows[0].stationShear.length;i++)svg+=`<text x="${48+i*492/(r.rows[0].stationShear.length-1)}" y="234" text-anchor="middle" font-size="12">${i+1}</text>`;return svg+lines.join('')+'<text x="290" y="258" text-anchor="middle" font-size="12">Location · same numbered sections as shear and dilution</text></svg>';};
  panel.innerHTML='<h3>Coagulation growth and depletion</h3><p>Fixed-pivot aggregation along each parcel. Dilution sets concentration; aggregation conserves powder mass. Capture × sticking efficiency: '+fmt(p.alpha,2)+'. Histories end at wall contact, outlet or the horizon. Transport follows the release-size path; aggregate shape, breakup and growth feedback on the path are excluded.</p><div class="analytics-charts"><div>'+chart(0)+'</div><div>'+chart(1)+'</div></div><p class="note">Fraction of initial particle count, 0–1. Means condition on reaching each location. Hover over points for counts. Colors match release offsets; this is distinct from first-event probability.</p><div class="table-wrap"><table><thead><tr><th>Release size / offset</th><th>Number depletion · mean ± SD</th><th>Mean volume diameter · µm</th><th>Max mass error</th><th>Overflow mass</th></tr></thead><tbody>'+r.rows.map(row=>`<tr><td>${sizeLabel(row.d)} / ${fmt(row.offset*1000,1)} mm</td><td>${stat(row.growth.numberLoss,100,2)}%</td><td>${stat(row.growth.meanDiameter,1e6,3)}</td><td>${row.growth.maxMassError.toExponential(2)}</td><td>${fmt(row.growth.maxOverflow*100,4)}%</td></tr>`).join('')+'</tbody></table></div><details><summary>Final size distributions · separate diameter grids</summary>'+r.rows.map(row=>'<details><summary>'+sizeLabel(row.d)+' / '+fmt(row.offset*1000,1)+' mm</summary><table><thead><tr><th>Diameter · µm</th><th>Mass fraction</th><th>Number per kg powder</th></tr></thead><tbody>'+row.samples[0].growth.distribution.map((bin,i)=>'<tr><td>'+fmt(bin.diameter*1e6,4)+'</td><td>'+fmt(row.samples.reduce((v,s)=>v+s.powderFluxWeight*s.growth.distribution[i].massFraction,0),6)+'</td><td>'+row.samples.reduce((v,s)=>v+s.powderFluxWeight*s.growth.distribution[i].numberPerKg,0).toExponential(3)+'</td></tr>').join('')+'</tbody></table></details>').join('')+'<p>Overflow mass is conserved but no longer collides. Material overflow requires a wider size grid. JSON export includes all distributions and station values.</p></details>';
 }
 function clear(){latest=null;if($('population-results'))$('population-results').remove();if($('numerical-quality'))$('numerical-quality').remove();for(const id of ['coag-mechanisms','coag-convergence','coag-summary','coag-table','coag-legend','stats-table','fate-bars','coag-small','coag-large','distribution-chart','lab-shear-legend'])$(id).innerHTML='';$('coag-screen').textContent='No valid particle histories for the current settings.';$('stats-caption').textContent='No current result.';$('samples-csv').disabled=true;$('lab-sampling-summary').textContent='No current particle statistics';if($('lab-kernel-status'))$('lab-kernel-status').textContent='';}
 function render(run){latest=run;renderQuality(run);renderPopulation(run);$('samples-csv').disabled=false;drawCoag();drawInline();drawStatistics();}
 function drawCoag(){const {p,r}=latest,rows=r.rows,key=$('coag-quantity').value,alpha=p.alpha??1,mult=key==='probability'?100:1,title=titles[key];
 $('coag-summary').innerHTML=[['Background',p.background==='mixture'?fmt(p.smallShare*100,0)+'% small / '+fmt((1-p.smallShare)*100,0)+'% large by mass':'Independent same-size scenarios'],['Capture × sticking',fmt(alpha,2)+(alpha===1?' · ideal limit':'')],['Histories',rows.reduce((s,x)=>s+x.n,0)+' · '+rows[0].n+' per group']].map(([a,b])=>`<div><span>${a}</span><strong>${b}</strong></div>`).join('');
 drawCumulative(['coag-small','coag-large'],key);
 $('coag-legend').innerHTML=rows.slice(0,3).map((x,i)=>`<span><i class="key" style="color:${COLORS[i]}"></i>${f(x.offset*1000,2)} mm</span>`).join('')+`<span>Solid: mean · dashed: median · ${$('coag-band').value==='ci'?'pointwise 95% CI of mean':'P05–P95 spread'} · accumulation ends with each track</span>`;
 const total=rows.reduce((a,x)=>a+x.distribution.hazard.mean*x.n,0),screened=rows.reduce((a,x)=>a+x.distribution.screened.mean*x.n,0),coverage=total>0?screened/total:null,volume=Math.max(...rows.map(x=>x.distribution.maxVolume.max)),maxDose=Math.max(...rows.map(x=>x.distribution.hazard.max));
 const label=Exposure.collisions.label(p.collisionModel)+(!['legacy','tracer-brownian'].includes(p.collisionModel)?(p.extendedKernels?' / kinetic extension':'')+(p.clustering?' + RDF sensitivity':''):'');
 const bySize=[0,1].map(species=>{const group=rows.filter(x=>x.species===species),total=group.reduce((s,x)=>s+x.distribution.hazard.mean,0),outside=group.reduce((s,x)=>s+x.distribution.unsupported.mean,0);return sizeLabel(group[0].d)+': '+(total?f(outside/total*100,1)+'% outside screen':'zero dose');}).join(' · ');
 const warning=total===0?'Zero encounter dose at these settings.':bySize+'. Dose shares pool the three sampled releases within each size; they are applicability diagnostics, not physical mixture averages. Outside-screen contributions are extrapolated; the screened subtotal is not a corrected total.';
 const numericalStops=rows.reduce((n,row)=>n+row.samples.filter(s=>s.fate==='numerical').length,0);
 $('coag-screen').textContent=label+'. '+(numericalStops?numericalStops+' histories stopped numerically. Their partial dose is incomplete, not a completed residence or encounter estimate. ':'')+warning+(maxDose>.1?(p.populationBalance?' Encounter plots retain a frozen background; the separate population balance includes number depletion.':' Some histories have Λ > 0.1; enable population balance to examine growth and depletion.'):'' )+(volume>.001?' Local powder volume fraction exceeds 0.1%; dilute one-way assumptions are questionable.':'');
 if($('lab-kernel-status')){$('lab-kernel-status').textContent=(numericalStops?numericalStops+' numerical stops · incomplete exposure; do not use terminal statistics. ':'')+label+' · '+bySize+' · first-event estimate';$('lab-kernel-status').title=warning+' See Coagulation for mechanism and per-size statistics.';}
 $('coag-mechanisms').innerHTML=rows.map((x,i)=>{const d=x.distribution;return `<tr><td style="color:${COLORS[i%3]}">${sizeLabel(x.d)} / ${f(x.offset*1000,2)} mm</td>${['brownian','turbulent','inertial'].map(k=>`<td>${f(d[k].mean)}</td>`).join('')}<td>${d.hazard.mean?f(100*d.unsupported.mean/d.hazard.mean,1)+'%':'—'}</td>${['maxSt','maxSize','maxCollisionRe'].map(k=>`<td>${f(d[k].max)}</td>`).join('')}</tr>`;}).join('');
 $('coag-table').innerHTML=rows.map((x,i)=>{const h=x.distribution.hazard;return `<tr><td style="color:${COLORS[i%3]}">${sizeLabel(x.d)} / ${f(x.offset*1000,2)} mm</td><td>${x.n}</td><td>${f(h.mean)} ± ${f(h.sd)}</td><td>${f(h.median)}</td><td>${f(h.p05)}–${f(h.p95)}</td><td>${h.ci?h.ci.map(v=>f(v)).join('–'):'—'}</td><td>${h.mean&&h.n>1?f(h.se/h.mean*100,2)+'%':'—'}</td><td>${f(h.max)}</td><td>${f(x.distribution.probability.mean*100,2)}%</td><td>${h.mean?f(x.distribution.screened.mean/h.mean*100,2)+'%':'—'}</td></tr>`;}).join('');
 $('coag-convergence').innerHTML=rows.map((x,i)=>{const c=x.cumulative.hazard.convergence,a=c.at(-2),z=c.at(-1);return `<tr><td style="color:${COLORS[i%3]}">${sizeLabel(x.d)} / ${f(x.offset*1000,2)} mm</td>${[c[0],c[Math.min(1,c.length-1)],z].map(v=>`<td>${f(v.mean)} <small>(n=${v.count})</small></td>`).join('')}<td>${a&&z.mean?f((z.mean-a.mean)/z.mean*100,2)+'%':'—'}</td></tr>`;}).join('');

 }
 function drawCumulative(ids,key,responsive=false){
  if(!latest)return;const {p,r}=latest,title=titles[key],ci=$('coag-band').value==='ci',traces=$('coag-traces').checked;
  for(let species=0;species<2;species++){
   const target=$(ids[species]),group=r.rows.filter(x=>x.species===species),W=responsive?Math.max(150,target.clientWidth||300):600,H=responsive?Math.max(115,target.clientHeight||220):270,left=responsive?(W<220?32:40):65,right=W-(responsive?10:25),bottom=responsive?34:52,base=H-bottom,top=responsive?28:45,ticks=W<320?2:4,xx=t=>left+(right-left)*t/p.horizon,yy=v=>base-v*(base-top),clip=ids[species]+'-clip';
   target.setAttribute('viewBox',`0 0 ${W} ${H}`);let svg=`<defs><clipPath id="${clip}"><rect x="${left}" y="${top}" width="${right-left}" height="${base-top}"/></clipPath></defs><text x="${left}" y="${responsive?17:22}" font-size="${responsive?12:17}" font-weight="650" fill="#173c51">${sizeLabel(group[0].d)}</text>`;
   const over=group.some(row=>{const c=row.cumulative[key];return c.mean.some(v=>v>1)||c.median.some(v=>v>1)||(ci?c.ciHigh:c.p95).some(v=>v>1)||(traces&&c.traces.some(t=>t.some(v=>v>1)));});
   if(!responsive||over)svg+=`<text x="${right}" y="${responsive?17:22}" text-anchor="end" font-size="${responsive?9:12}" fill="${over?'#9a5f1a':'#536b78'}"><title>Fixed display range only; full values are retained in tables and exports.</title>${responsive?'':title}${over?' ↑ above 1':''}</text>`;
   for(let j=0;j<=4;j++){const v=j/4;svg+=`<line x1="${left}" x2="${right}" y1="${yy(v)}" y2="${yy(v)}" stroke="#dce5eb"/><text x="${left-5}" y="${yy(v)+3}" text-anchor="end" font-size="${responsive?9:11}" fill="#536b78">${f(v,2)}</text>`;}
   for(let j=0;j<=ticks;j++)svg+=`<text x="${xx(p.horizon*j/ticks)}" y="${H-(responsive?21:31)}" text-anchor="middle" font-size="${responsive?9:11}" fill="#536b78">${f(p.horizon*j/ticks,1)}</text>`;
   svg+=`<g clip-path="url(#${clip})">`;group.forEach((row,i)=>{const c=row.cumulative[key],point=(v,j)=>xx(j*p.horizon/240)+','+yy(v),low=ci?c.ciLow:c.p05,high=ci?c.ciHigh:c.p95;
    if(!ci||row.n>1)svg+=`<path d="M${high.map(point).join('L')}L${low.map(point).reverse().join('L')}Z" fill="${COLORS[i]}" opacity=".14"/>`;
    if(traces)c.traces.forEach(t=>{svg+=`<path d="M${t.map(point).join('L')}" fill="none" stroke="${COLORS[i]}" stroke-width=".7" opacity=".14"/>`;});
    svg+=`<path d="M${c.median.map(point).join('L')}" fill="none" stroke="${COLORS[i]}" stroke-width="1" stroke-dasharray="4 3"/><path d="M${c.mean.map(point).join('L')}" fill="none" stroke="${COLORS[i]}" stroke-width="${responsive?1.7:2.4}"/>`;});svg+=`</g><text x="${(left+right)/2}" y="${H-5}" text-anchor="middle" font-size="${responsive?9:12}" fill="#536b78">Time from release (ms)</text>`;target.innerHTML=svg;
  }
 }
 function drawInline(){if(!latest)return;const rows=latest.r.rows,errors=rows.map(x=>x.distribution.shearDose).filter(s=>s.n>1&&s.mean>0).map(s=>s.se/s.mean*100);$('lab-shear-legend').innerHTML=rows.slice(0,3).map((x,i)=>`<span><i class="key" style="color:${COLORS[i]}"></i>${f(x.offset*1000,2)} mm</span>`).join('');$('lab-sampling-summary').textContent=`${rows.reduce((n,x)=>n+x.n,0).toLocaleString()} histories · ${rows[0].n}/group · ${errors.length?'shear dose max relative SE '+f(Math.max(...errors),1)+'%':'relative SE unavailable'} · statistics ↗`;}
 function drawStatistics(){const {p,r}=latest,key=$('stats-metric').value,[title,scale]=metricInfo[key],values=r.rows.map(x=>x.distribution[key]),max=Math.max(1e-12,...values.filter(x=>x.n).map(x=>Math.max(x.p95,x.mean)*scale))*1.1,xx=v=>175+700*v/max;
 $('stats-caption').textContent=`${r.rows.reduce((s,x)=>s+x.n,0)} histories · seed ${p.seed} · ${p.noise?'seeded turbulent fluctuations':'deterministic field'} · ${key==='outletTime'?'only histories reaching the outlet':'all histories, including wall contacts and censored tracks'}`;
 let svg=`<text x="175" y="20" font-size="14" fill="#173c51">${title}</text>`;for(let j=0;j<=4;j++){const v=max*j/4;svg+=`<line x1="${xx(v)}" x2="${xx(v)}" y1="31" y2="278" stroke="#e1e8ed"/><text x="${xx(v)}" y="300" text-anchor="middle" font-size="12" fill="#536b78">${f(v,2)}</text>`;}
 r.rows.forEach((row,i)=>{const y=51+i*40,s=values[i],c=COLORS[i%3];svg+=`<text x="163" y="${y+4}" text-anchor="end" font-size="12" fill="${c}">${sizeLabel(row.d)} / ${f(row.offset*1000,1)} mm</text>`;if(!s.n){svg+=`<text x="190" y="${y+4}" font-size="12" fill="#647785">No completed outlet histories</text>`;return;}const x=k=>xx(s[k]*scale);svg+=`<line x1="${x('p05')}" x2="${x('p95')}" y1="${y}" y2="${y}" stroke="${c}" stroke-width="2"/><path d="M${x('p05')} ${y-6}V${y+6}M${x('p95')} ${y-6}V${y+6}" stroke="${c}"/><rect x="${x('p25')}" y="${y-10}" width="${Math.max(1,x('p75')-x('p25'))}" height="20" fill="${c}" fill-opacity=".17" stroke="${c}"/><line x1="${x('median')}" x2="${x('median')}" y1="${y-10}" y2="${y+10}" stroke="${c}" stroke-width="2"/><circle cx="${x('mean')}" cy="${y}" r="3.5" fill="${c}"/>`;});$('distribution-chart').innerHTML=svg;
 $('stats-table').innerHTML=r.rows.map((x,i)=>{const s=values[i];return `<tr><td style="color:${COLORS[i%3]}">${sizeLabel(x.d)} / ${f(x.offset*1000,1)} mm</td><td>${s.n}</td>${['mean','sd','se'].map(k=>`<td>${f(s[k]===null?null:s[k]*scale)}</td>`).join('')}<td>${s.ci?s.ci.map(v=>f(v*scale)).join('–'):'—'}</td>${['p05','median','p95','p99'].map(k=>`<td>${f(s[k]===null?null:s[k]*scale)}</td>`).join('')}<td>${f(s.cv,2)}</td></tr>`;}).join('');
 const colors={wall:'#bf5015',outlet:'#008678',inlet:'#8061ad',censored:'#a2adb5',numerical:'#d32f2f'};$('fate-bars').innerHTML=r.rows.map((row,i)=>`<div class="fate-row"><span style="color:${COLORS[i%3]}">${sizeLabel(row.d)} / ${f(row.offset*1000,1)} mm</span><div class="fate-stack">${Object.entries(row.fates).filter(([k,v])=>v.count).map(([k,v])=>`<span style="width:${v.fraction*100}%;background:${colors[k]}" title="${k}: ${v.count}/${row.n}; 95% Wilson interval ${v.ci.map(x=>f(x*100,1)).join('–')}%">${v.fraction>=.1?f(v.fraction*100,0)+'%':''}</span>`).join('')}</div><span>${row.n} tracks</span></div>`).join('');
 $('stats-caution').textContent='Hover over fate segments for counts and 95% Wilson intervals. Mean intervals use a Student-t approximation; they describe sampling uncertainty only, not uncertainty in the gas field or collision model. Tracked time includes wall stops and the finite horizon; outlet transit conditions on successful exit. '+(p.noise?'Increase histories and compare seeds to assess sampling stability.':'With fluctuations off, each group has one history; a sampling interval cannot be estimated.');
 }
 function exportSamples(){if(!latest)return;const {g,p,r}=latest,columns=['group','history','seed','branch','powder_flux_weight','contact','deposited','diameter_um','offset_mm','fate','time_s','half_s','quarter_s','peak_G_s_inv','shear_integral','epsilon_integral_m2_s2','collision_dose','screened_dose','sticking_probability','max_St_eta','max_collision_radius_over_eta','max_turbulent_Mach','max_powder_kg_m3','max_solids_volume_fraction','brownian_dose','turbulent_dose','inertial_increment','outside_screen_dose','max_collision_slip_Re'],records=r.rows.flatMap((row,i)=>row.samples.map(s=>[i+1,s.index,s.seed,s.branch,s.powderFluxWeight,s.contact,s.deposited===null?'unknown':s.deposited,row.d*1e6,row.offset*1000,s.fate,s.t,s.half,s.quarter,s.peak,s.shearDose,s.epsilonDose,s.hazard,s.screened,s.probability,s.maxSt,s.maxSize,s.maxMt,s.maxConcentration,s.maxVolume,s.brownian,s.turbulent,s.inertial,s.unsupported,s.maxCollisionRe]));const data=[['geometry',g.name],['collision_model',p.collisionModel||'hybrid'],['kernel_version',Exposure.collisions.version],['background',p.background],['small_mass_fraction',p.smallShare],['capture_sticking_efficiency',p.alpha],['motive_bar_abs',p.p0],['vessel_bar_abs',p.ps],['outlet_bar_abs',p.po],['powder_kg_h',p.feed],['caution','Frozen-background encounter estimate with explicit applicability screening; no growth/depletion; see JSON export'],[],columns,...records].map(a=>a.map(x=>'"'+String(x).replaceAll('"','""')+'"').join(',')).join('\r\n'),url=URL.createObjectURL(new Blob([data],{type:'text/csv'})),a=document.createElement('a');a.href=url;a.download='eductor-histories-'+g.key+'.csv';a.click();setTimeout(()=>URL.revokeObjectURL(url),1000);}
 return {init,render,clear,drawInline,drawStatistics};
})();

// Browser-only picture capture. No external service, re-simulation or wall changes.
const ExposureSnapshots=(()=>{
 const by=id=>document.getElementById(id),NS='http://www.w3.org/2000/svg';
 const escape=s=>String(s).replaceAll('&','&amp;').replaceAll('<','&lt;').replaceAll('"','&quot;');
 const uri=svg=>'data:image/svg+xml;charset=utf-8,'+encodeURIComponent(svg);
 const shell=(w,h,content,bg='white')=>`<svg xmlns="${NS}" width="${w}" height="${h}" viewBox="0 0 ${w} ${h}"><rect width="100%" height="100%" fill="${bg}"/><g font-family="system-ui, sans-serif">${content}</g></svg>`;
 function svgImage(el){const copy=el.cloneNode(true);copy.querySelectorAll('.saved-experiment-curve,.saved-experiment-point').forEach(n=>n.remove());copy.setAttribute('xmlns',NS);const box=el.viewBox.baseVal;copy.setAttribute('width',box.width||600);copy.setAttribute('height',box.height||400);copy.style.cssText='font-family:system-ui,sans-serif;background:white';return uri(new XMLSerializer().serializeToString(copy));}
 async function raster(source,maxWidth=1400){const img=new Image();img.src=source;await img.decode();const scale=Math.min(1,maxWidth/img.naturalWidth),c=document.createElement('canvas');c.width=Math.max(1,Math.round(img.naturalWidth*scale));c.height=Math.max(1,Math.round(img.naturalHeight*scale));const ctx=c.getContext('2d');ctx.fillStyle='#fff';ctx.fillRect(0,0,c.width,c.height);ctx.drawImage(img,0,0,c.width,c.height);return c.toDataURL('image/png');}
 function geometrySource(){
  if(!by('three-view').hidden)return {mode:'3D · camera and playback frozen',source:ExposureViewer.capture()};
  if(!by('blender-view').hidden){const img=by('blender-view');if(!img.complete||!img.naturalWidth)throw Error('Wait for the geometry image to finish loading, then keep again.');const w=Math.max(200,img.clientWidth),h=Math.max(200,img.clientHeight),overlay=by('dilution-overlay').innerHTML;return {mode:'Blender · numbered locations frozen',source:uri(shell(w,h,`<image x="5" y="20" width="${w-10}" height="${h-25}" href="${escape(img.src)}" preserveAspectRatio="xMidYMid meet"/><svg x="5" y="20" width="${w-10}" height="${h-25}" viewBox="0 0 960 640">${overlay}</svg>`,'#edf4f8'))};}
  const c=by('field');return {mode:'2D · field and paths frozen',source:uri(shell(c.width,c.height,`<image width="${c.width}" height="${c.height}" href="${c.toDataURL('image/png')}"/><svg width="${c.width}" height="${c.height}" viewBox="0 0 900 510">${by('paths').innerHTML}</svg>`))};
 }
 function pairedChart(run,quantity){const w=620,h=300;let body='';for(let species=0;species<2;species++){const rows=run.r.rows.filter(r=>r.species===species);body+=`<svg y="${species*h}" width="${w}" height="${h}" viewBox="0 0 ${w} ${h}">${ExposureStationChart.svg({rows,Gref:run.r.Gref,width:w,height:h,active:ExposureDilution.selectedIndex(),title:sizeLabel(rows[0].d),id:'frozen-'+(quantity||'shear')+'-'+species,quantity})}</svg>`;}return uri(shell(w,h*2,body));}
 async function capture(run){
  // Everything is collected synchronously before decoding images. Controls may
  // change while decoding; the captured geometry, data and labels stay coherent.
  const geometry=geometrySource(),quantity=by('lab-coag-quantity').value,
   summary=[...by('kpis').querySelectorAll('.kpi')].map(el=>({label:el.querySelector('span').textContent,value:el.querySelector('strong').textContent,detail:el.querySelector('small').textContent})),
   sources={geometry:geometry.source,shear:pairedChart(run,null),dilution:svgImage(by('lab-dilution-chart')),coagulation:pairedChart(run,quantity)},
   labels={geometry:geometry.mode,shear:'Shear · G / Gref',dilution:by('dilution-quantity').selectedOptions[0].textContent,coagulation:by('lab-coag-quantity').selectedOptions[0].textContent},
   releaseLegend=run.r.rows.slice(0,3).map(r=>fmt(r.offset*1000,2)+' mm').join(' / '),timeMS=+by('time-control').value;
  const images=Object.fromEntries(await Promise.all(Object.entries(sources).map(async([k,v])=>[k,await raster(v)])));
  return {images,labels,summary,releaseLegend,timeMS};
 }
 let dialog,previews={},initialized=false;
 const titles={geometry:'Geometry & paths',shear:'Shear',dilution:'Dilution',coagulation:'Coagulation'};
 function init(){if(initialized)return;initialized=true;
  for(const [key,host] of [['geometry',by('dilution-station-key').closest('.path-card')],['shear',by('shear-history-panel')],['dilution',by('dilution-history-panel')],['coagulation',by('coagulation-history-panel')]]){
   const button=document.createElement('button');button.type='button';button.className='saved-picture';button.dataset.savedPanel=key;button.hidden=true;button.innerHTML='<span>Saved · '+titles[key]+' <i>Expand ↗</i></span><img alt="Saved experiment '+titles[key]+'">';button.onclick=()=>open(key);if(key==='geometry')host.querySelector('.visual-stage').after(button);else host.append(button);previews[key]=button;
  }
  dialog=document.createElement('dialog');dialog.className='saved-view-dialog';dialog.id='saved-view-dialog';dialog.innerHTML='<header><h2 id="saved-view-title">Saved experiment</h2><button type="button" id="close-saved-view">Close</button></header><p id="saved-view-caption"></p><div id="saved-view-summary"></div><div id="saved-view-images"></div><p class="saved-view-note">Frozen pictures of the kept experiment. These pictures do not update when settings change.</p>';dialog.setAttribute('aria-labelledby','saved-view-title');document.body.append(dialog);by('close-saved-view').onclick=()=>dialog.close();dialog.addEventListener('click',e=>{if(e.target===dialog)dialog.close();});
 }
 function refresh(){if(!initialized)return;const q=ExposureOverlay.get(),v=q?.visuals;document.body.classList.toggle('has-saved-pictures',!!v);
  for(const [key,button] of Object.entries(previews)){button.hidden=!v;if(v){const img=button.querySelector('img');if(img.getAttribute('src')!==v.images[key])img.src=v.images[key];button.title='Saved '+q.name+' · '+v.labels[key]+' · click to expand';button.setAttribute('aria-label','Expand saved '+titles[key]+' from '+q.name);}}
  by('kpis').querySelectorAll('.saved-kpi').forEach(n=>n.remove());if(v)[...by('kpis').children].forEach((card,i)=>{const s=v.summary[i];if(!s)return;const b=document.createElement('button');b.type='button';b.className='saved-kpi';b.textContent='Saved: '+s.value;b.title=q.name+' · '+s.label+' · '+s.detail;b.onclick=()=>open('summary');card.append(b);});
  if(by('saved-panels'))by('saved-panels').hidden=!ExposureOverlay.peek()?.visuals;if(!ExposureOverlay.has()&&dialog.open)dialog.close();
 }
 function open(key='all'){const q=ExposureOverlay.peek(),v=q?.visuals;if(!v)return;by('saved-view-title').textContent='Saved · '+q.name;by('saved-view-caption').textContent=`Motive ${fmt(q.settings.p0,3)} bar abs · vessel = outlet ${fmt(q.settings.ps,3)} bar abs · powder ${fmt(PowderRate.fromModel(q.settings.feed),3)} g/s · release positions ${v.releaseLegend} · playback ${fmt(v.timeMS,2)} ms`;
  const summary=by('saved-view-summary');summary.replaceChildren();for(const s of v.summary){const card=document.createElement('div'),label=document.createElement('span'),value=document.createElement('strong'),detail=document.createElement('small');label.textContent=s.label;value.textContent=s.value;detail.textContent=s.detail;card.append(label,value,detail);summary.append(card);}
  const images=by('saved-view-images');images.replaceChildren();dialog.dataset.focus=key;for(const k of key==='all'?Object.keys(titles):Object.hasOwn(titles,key)?[key]:[]){const figure=document.createElement('figure'),caption=document.createElement('figcaption'),img=document.createElement('img');caption.textContent=titles[k]+' · '+v.labels[k];img.src=v.images[k];img.alt='Saved '+titles[k]+' for '+q.name;figure.append(caption,img);images.append(figure);}if(!dialog.open)dialog.showModal();
 }
 return {capture,init,refresh,open};
})();

// Reuse the scientific controls and tables in a fixed-height, tabbed workspace.
(()=>{
 const main=document.querySelector('main'),nav=main.querySelector('nav'),header=main.querySelector('header'),controls=$('controls-title').closest('section'),plots=main.querySelector('.plots'),results=$('results'),metrics=$('metrics').closest('section'),comparison=$('slide-comparison'),audit=$('model-audit'),methods=$('methods');
 const make=(tag,cls,html='')=>{const el=document.createElement(tag);el.className=cls;el.innerHTML=html;return el;};
 const tabs=make('div','workspace-tabs');tabs.setAttribute('role','tablist');tabs.setAttribute('aria-label','Exposure workspace');const pages=[['lab','Exposure lab'],['coagulation','Coagulation'],['statistics','Statistics'],['metrics','Metrics'],['slide','Slide comparison'],['compressible','Compressible flow'],['review','Model review']];tabs.innerHTML=pages.map(([id,name])=>`<button role="tab" id="tab-${id}" aria-controls="page-${id}" aria-selected="${id==='lab'}" data-page="${id}">${name}</button>`).join('');
 const grid=make('div','app-grid'),stage=make('div','stage');controls.className='card experiment';const pageEls={};for(const [id] of pages){const page=make('section','workspace-page');page.dataset.page=id;page.id='page-'+id;page.setAttribute('role','tabpanel');page.setAttribute('aria-labelledby','tab-'+id);page.hidden=id!=='lab';pageEls[id]=page;stage.append(page);}
 const heading=make('div','result-heading');heading.append($('result-title'),$('run-caption'));const kpis=$('kpis'),status=$('flow-status');status.classList.add('flow-bar');const mobileSwitch=make('div','mobile-chart-switch','<button id="show-geometry" aria-selected="true">Geometry</button><button id="show-histories" aria-selected="false">Histories</button>');
 results.replaceChildren(heading,kpis,mobileSwitch,plots);pageEls.lab.append(results,status);
 const exportRow=make('div','export-row');$('csv').textContent='Export CSV';$('json').textContent='Export JSON';exportRow.append($('csv'),$('json'));const metricsHead=make('div','page-head','<div><h2>Exposure metrics</h2><p>Six release groups · mean ± sample SD</p></div>');metricsHead.append(exportRow);pageEls.metrics.append(metricsHead,metrics);pageEls.slide.append(comparison);pageEls.review.append(audit,methods);
 const controlHead=make('div','section-head','<h2 id="controls-title">Experiment</h2><button class="mobile-close" id="close-controls">Done</button>'),controlTabs=make('div','control-tabs');controlTabs.setAttribute('role','tablist');controlTabs.setAttribute('aria-label','Experiment settings');const sets=[['operating','Operating',['geometry','flow-mode','p0','ps','po','feed','temp']],['particles','Particles',['collisionModel','kernelRange','small','large','density','release','offset','noise','background','smallShare','alpha','populationBalance']],['numerics','Model',['intensity','lengthScale','spread','reference','n','dt','horizon','seed','loss','sampling','phase-model','wallModel','kernelConvention','adaptiveDt','antithetic','commonNoise']]];
 const controlPages=[];for(const [id,title,fields] of sets){const btn=make('button','',title);btn.type='button';btn.setAttribute('role','tab');btn.setAttribute('aria-selected',id==='operating');btn.id='control-tab-'+id;btn.setAttribute('aria-controls','control-page-'+id);btn.dataset.control=id;controlTabs.append(btn);const page=make('div','control-page');page.id='control-page-'+id;page.hidden=id!=='operating';page.setAttribute('role','tabpanel');page.setAttribute('aria-labelledby',btn.id);const fieldGrid=make('div','fields');for(const field of fields){const el=$(field).closest('.field')||$(field).closest('label');if(['geometry','flow-mode','density','release','offset','noise','reference','loss','background','smallShare','alpha','collisionModel','kernelRange','populationBalance'].includes(field))el.classList.add('span2');fieldGrid.append(el);}page.append(fieldGrid);if(id==='operating'){page.append($('assembly-config'));page.append($('prescribed-fields'));page.append(make('div','solved-flow','<h3>Calculated airflow</h3><div id="solved-flow-values">Solving pressure balance…</div><small>Low-Mach or compressible balance, selected automatically.</small>'));page.lastElementChild.id='solved-flow';page.append(make('p','note','Vessel and outlet always share this absolute pressure. Chamber and interstage pressures are calculated.'));}if(id==='particles')page.append(make('p','note','Sizes are compact-sphere diameters. CAD cases start at the common suction inlet; the release-distance input applies to nominal geometries. Offsets scale with port radius.'));controlPages.push(page);}
 const actions=make('div','experiment-actions'),runrow=make('div','toolbar');$('run').hidden=true;runrow.append($('run'));const preset=make('div','preset-row');$('slide-preset').textContent='Slide · 3 / 1 bar';$('vacuum').textContent='Vacuum · 0.2 atm';preset.append($('slide-preset'),$('vacuum'));actions.append(make('div','auto-update-label','● Automatic updates'),preset,$('busy-note'));controls.replaceChildren(controlHead,controlTabs,...controlPages,actions,runrow);runrow.hidden=true;
 for(const [id,text] of Object.entries({p0:'Motive · bar abs',ps:'Vessel = outlet · bar abs',po:'Outlet · bar abs',feed:'Powder · g/s',temp:'Temperature · °C'}))controls.querySelector(`label[for="${id}"]`).textContent=text;$('po').closest('.field').hidden=true;$('ps').closest('.field').classList.add('span2');const open=make('button','mobile-actions','Experiment');open.id='open-controls';header.append(open);grid.append(controls,stage);main.replaceChildren(nav,header,tabs,grid);const lock=document.querySelector('body>form[action="/logout"]');if(lock){lock.removeAttribute('style');lock.className='lock-form';const button=lock.querySelector('button');button.removeAttribute('style');button.textContent='Lock';button.setAttribute('aria-label','Lock dashboard');header.append(lock);main.classList.add('has-lock');}
 // Larger 3D/2D scene and synchronized shear plots.
 const pathCard=plots.children[0],historyCard=plots.children[1];pathCard.classList.add('path-card');historyCard.classList.add('history-card');const plotWrap=pathCard.querySelector('.plot-wrap'),playbar=pathCard.querySelector('.toolbar'),pathNote=pathCard.querySelector('p.note'),legend=$('path-legend');pathNote.remove();pathCard.querySelector('h3').remove();const top=make('div','plot-topline','<h3>Geometry & particle paths</h3><div class="visual-tabs" role="tablist" aria-label="Geometry view"><button id="view-render" role="tab" aria-selected="true">Blender</button><button id="view-3d" role="tab" aria-selected="false">3D</button><button id="view-2d" role="tab" aria-selected="false">2D</button></div>');
 const visual=make('div','visual-stage','<div id="three-view" hidden></div><img id="blender-view" alt="Blender-rendered eductor geometry"><svg id="dilution-overlay" viewBox="0 0 960 640" role="group" aria-label="Dilution reference stations on the Blender geometry"></svg><div class="viewer-badge" id="viewer-provenance"></div>');plotWrap.hidden=true;visual.append(plotWrap);const viewtools=make('div','viewer-toolbar','<select id="view-surface" aria-label="3D surface"><option value="cad">STEP geometry</option><option value="envelope">Particle domain</option></select><button id="camera-iso">Iso</button><button id="camera-side">Side</button><button id="camera-end">End</button><label><input type="checkbox" id="view-cut">Cutaway</label><label><input type="checkbox" id="view-ghost">Ghost</label>');playbar.className='playback';const viewerNote=make('p','plot-note');viewerNote.id='viewer-note';const stationControls=make('div','dilution-stations','<div id="dilution-station-key" aria-label="Dilution stations"></div><p id="dilution-station-detail" aria-live="polite"></p>');pathCard.replaceChildren(top,legend,visual,stationControls,viewtools,playbar,viewerNote);
 const overlayTools=make('div','overlay-tools','<button id="keep-experiment" type="button" disabled>Keep experiment</button><label><input id="show-overlay" type="checkbox" disabled> Overlay saved</label><button id="saved-panels" type="button" hidden>Saved panels</button><button id="clear-overlay" type="button" hidden>Clear</button><span id="overlay-description">Keep a run to compare settings on these charts.</span>');
 function overlayDescription(){const q=ExposureOverlay.get();const label=q?'Dashed: '+q.name+' · '+fmt(q.settings.p0,2)+' / '+fmt(q.settings.ps,2)+' bar abs · '+fmt(PowderRate.fromModel(q.settings.feed),3)+' g/s · '+q.settings.n+'/group. '+(q.geometry!==current?.g.key?'Different geometry: locations match by role. ':'')+'Shear uses each run’s own Gref. Saved in this tab.':ExposureOverlay.has()?'Saved experiment hidden.':'Keep a run to compare settings on these charts.';$('overlay-description').textContent=label;$('overlay-description').title=label;}
 let keepCapture=0;
 overlayTools.querySelector('#keep-experiment').onclick=async()=>{
  if(!current||document.body.classList.contains('model-pending'))return;const token=++keepCapture,button=$('keep-experiment'),run=current,snapshot=ExposureOverlay.snapshot(run,ExposureDilution.stations(run)),wasEnabled=!!ExposureOverlay.get();stop();button.disabled=true;button.textContent='Saving pictures…';
  try{ExposureOverlay.toggle(false);drawHistories();ExposureDilution.draw();const pending=ExposureSnapshots.capture(run);ExposureOverlay.toggle(wasEnabled);drawHistories();ExposureDilution.draw();snapshot.visuals=await pending;if(token!==keepCapture)return;ExposureOverlay.store(snapshot);$('show-overlay').disabled=false;$('show-overlay').checked=true;$('clear-overlay').hidden=false;overlayDescription();if(current)drawHistories();ExposureDilution.draw();ExposureSnapshots.refresh();}
  catch(error){ExposureOverlay.toggle(wasEnabled);$('overlay-description').textContent='Could not keep pictures: '+error.message;if(current)drawHistories();ExposureDilution.draw();}
  finally{button.textContent='Keep experiment';button.disabled=!current||document.body.classList.contains('model-pending');}
 };
 overlayTools.querySelector('#show-overlay').onchange=()=>{ExposureOverlay.toggle($('show-overlay').checked);overlayDescription();if(current)drawHistories();ExposureDilution.draw();ExposureSnapshots.refresh();};
 overlayTools.querySelector('#clear-overlay').onclick=()=>{keepCapture++;ExposureOverlay.clear();$('show-overlay').checked=false;$('show-overlay').disabled=true;$('clear-overlay').hidden=true;overlayDescription();if(current)drawHistories();ExposureDilution.draw();ExposureSnapshots.refresh();};
 overlayTools.querySelector('#saved-panels').onclick=()=>ExposureSnapshots.open();
 const histSmall=$('history-small'),histLarge=$('history-large'),histBox=make('div','histories');histBox.append(histSmall,histLarge);
 const historyHead=make('div','history-heading','<h3>Shear, dilution & coagulation</h3>');
 const sampleTools=make('div','history-sampling','<label>Histories/group <select id="lab-n" aria-label="Histories per group">'+$('n').innerHTML+'</select></label><span>Bulk dilution applies to both particle sizes.</span>');$('n').value='128';
 const historyGrid=make('div','history-comparison'),shearPanel=make('div','history-panel','<h4>Shear by location · G / Gref</h4>'),dilutionPanel=make('div','history-panel','<div class="dilution-heading"><h4>Dilution by location</h4><label for="dilution-quantity">Quantity</label><select id="dilution-quantity"><option value="powderAir">Powder / air · kg/kg</option><option value="concentration">Concentration · g/m³</option></select></div><div class="histories" id="lab-dilution-plots"><svg id="lab-dilution-chart" class="hist" viewBox="0 0 300 440" role="img" aria-label="Powder to air mass ratio in kg/kg at twelve numbered locations"></svg></div>');shearPanel.id='shear-history-panel';dilutionPanel.id='dilution-history-panel';shearPanel.append(histBox);const coagPanel=make('div','history-panel','<div class="dilution-heading"><h4>Coagulation by location</h4><label for="lab-coag-quantity">Quantity</label><select id="lab-coag-quantity" aria-label="Coagulation by location quantity"><option value="probability">First-sticking probability P</option><option value="hazard">Cumulative collision dose Λ</option><option value="screened">Within-screen subtotal Λ</option><option value="brownian">Brownian collision dose</option><option value="turbulent">Turbulent collision dose</option><option value="unsupported">Outside-screen dose · diagnostic</option></select></div><div class="histories" id="lab-coag-plots"><svg id="lab-coag-small" class="hist" role="img" aria-label="Small particle coagulation by location"></svg><svg id="lab-coag-large" class="hist" role="img" aria-label="Large particle coagulation by location"></svg></div>');coagPanel.id='coagulation-history-panel';historyGrid.append(shearPanel,dilutionPanel,coagPanel);sampleTools.querySelector('#lab-n').value=$('n').value;
 const mobileMetrics=make('div','mobile-metric-switch','<button data-metric="shear" aria-pressed="true">Shear</button><button data-metric="dilution" aria-pressed="false">Dilution</button><button data-metric="coagulation" aria-pressed="false">Coagulation</button>');mobileMetrics.setAttribute('aria-label','Location chart');historyGrid.dataset.mobileMetric='shear';mobileMetrics.onclick=e=>{const button=e.target.closest('button[data-metric]');if(!button)return;historyGrid.dataset.mobileMetric=button.dataset.metric;for(const b of mobileMetrics.children)b.setAttribute('aria-pressed',b===button);if(current)drawHistories();ExposureDilution.draw();};
 const historyLegend=make('div','legend');historyLegend.id='lab-shear-legend';const sizeLegend=make('div','compact-species-note');sizeLegend.id='lab-size-legend';const selectedStation=make('div','dilution-current');selectedStation.id='lab-dilution-selection';const historyNote=make('p','plot-note','First crossing: shear mean ± SD; cumulative coagulation mean / P05–P95. Gaps = not reached; hover for counts. Particle axes 0–1; ↑ overflow. Selected collision kernel; first-event estimate. Growth, if enabled, is reported separately.');const samplingSummary=make('button','sampling-summary','Particle statistics');samplingSummary.id='lab-sampling-summary';samplingSummary.onclick=()=>{$('stats-metric').value='shearDose';ExposureAnalytics.drawStatistics();setView('statistics');};const kernelStatus=make('p','plot-note');kernelStatus.id='lab-kernel-status';historyCard.replaceChildren(historyHead,kernelStatus,sampleTools,overlayTools,mobileMetrics,historyLegend,sizeLegend,historyGrid,selectedStation,historyNote,samplingSummary);$('lab-n').onchange=()=>{$('n').value=$('lab-n').value;markDirty();};

 const trackingNotice=make('div','tracking-notice');trackingNotice.id='tracking-notice';trackingNotice.hidden=true;historyCard.append(trackingNotice);
 const recovery=make('div','flow-recovery');recovery.id='flow-recovery';recovery.setAttribute('role','status');recovery.hidden=true;historyCard.append(recovery);
 const regularHistory=[kernelStatus,sampleTools,mobileMetrics,historyLegend,sizeLegend,historyGrid,selectedStation,historyNote,samplingSummary,trackingNotice];
 function showRecovery(issue){
  for(const el of regularHistory)el.hidden=!!issue;recovery.hidden=!issue;historyCard.classList.toggle('has-recovery',!!issue);if(!issue)return;plots.classList.add('show-history');$('show-geometry').setAttribute('aria-selected','false');$('show-histories').setAttribute('aria-selected','true');
  const escape=s=>String(s).replaceAll('&','&amp;').replaceAll('<','&lt;').replaceAll('"','&quot;');
  recovery.innerHTML=`<span class="recovery-tag">${issue.code==='particles'?'Gas solved · particle limit':'Operating-point check'}</span><h3>${escape(issue.title)}</h3><p>${escape(issue.explanation)}</p><div class="pressure-triplet">${[['Motive',issue.settings.p0],['Vessel',issue.settings.ps],['Outlet',issue.settings.po]].map(([k,v])=>`<div><span>${k}</span><strong>${fmt(v,3)}</strong><small>bar absolute</small></div>`).join('')}</div>${issue.threshold?`<p class="threshold-note">Approximate zero-suction vessel threshold: <b>${fmt(issue.threshold,3)} bar abs</b> with the selected motive, outlet and powder feed. This is a diagnostic of the reduced balance, not a measured operating limit.</p>`:''}<div class="recovery-options"><h4>${issue.alternatives.length?'Checked alternatives':'Next step'}</h4>${issue.alternatives.length?'<p>Each option passes the gas and field checks. Choose one to change only the stated pressures.</p>':'<p>'+escape(issue.code==='particles'?'Review the particle release, turbulence settings and model assumptions. A gas solution alone does not guarantee valid particle histories.':issue.code==='compressible'?'The automatic compressible estimate also found no admissible branch. Resolved shock / mixing losses or measured flow inputs are needed for this condition.':'Adjust the pressure boundaries or review the model assumptions. No manual airflow has been inserted.')+'</p>'}</div>`;
  for(const option of issue.alternatives){const button=make('button','recovery-option');button.type='button';button.innerHTML='<b>'+escape(option.label)+'</b><span>Incoming air '+fmt(option.air,2)+' g/s · chamber '+fmt(option.chamber,3)+' bar</span>';button.onclick=()=>{for(const [key,value] of Object.entries(option.changes))$(key).value=value;markDirty();};recovery.querySelector('.recovery-options').append(button);}
 }
 new ResizeObserver(()=>ExposureDilution.draw()).observe($('lab-dilution-plots'));
 $('lab-coag-quantity').onchange=drawStationCoagulation;new ResizeObserver(drawStationCoagulation).observe($('lab-coag-plots'));
 function showTrackingNotice(){
  const {g,p,r}=current,total=r.rows.reduce((s,row)=>s+row.n,0),unfinished=r.rows.reduce((s,row)=>s+row.censored*row.n,0)/total;
  trackingNotice.hidden=unfinished<=.05;if(trackingNotice.hidden)return;
  const [x,y,z=0]=r.field.geo.release(0),entry=r.field.sample(x,y,z),speed=Math.max(.001,Math.hypot(entry.ux,entry.uy,entry.uz||0)),out=r.field.axial(r.field.geo.end);
  const travelMs=(r.field.geo.branches?Math.max(...r.field.geo.branches.map(b=>b.length))*1000:p.release)/speed+1000*r.field.geo.end/Math.max(.001,out.mean),target=Math.max(p.horizon*1.5,travelMs*1.5),options=Array.from($('horizon').options,o=>+o.value),next=options.find(v=>v>=target)||options.findLast(v=>v>p.horizon);
  trackingNotice.replaceChildren(make('span','',fmt(unfinished*100,0)+'% of histories are still active at '+p.horizon+' ms.'));
  if(next){const extend=make('button','','Extend tracking to '+next+' ms');extend.type='button';extend.onclick=()=>{$('horizon').value=next;markDirty();};trackingNotice.append(extend);}
 }
 function setView(id){for(const [key] of pages){pageEls[key].hidden=key!==id;const b=$('tab-'+key);b.setAttribute('aria-selected',key===id);b.tabIndex=key===id?0:-1;}if(id==='compressible')drawCompressible();ExposureViewer.resize();}
 tabs.addEventListener('click',e=>{const b=e.target.closest('button[data-page]');if(b)setView(b.dataset.page);});
 function keyboardTabs(container,fn){container.addEventListener('keydown',e=>{if(!['ArrowLeft','ArrowRight','Home','End'].includes(e.key))return;const buttons=[...container.querySelectorAll('[role=tab]')],i=buttons.indexOf(document.activeElement);if(i<0)return;e.preventDefault();const j=e.key==='Home'?0:e.key==='End'?buttons.length-1:(i+(e.key==='ArrowRight'?1:-1)+buttons.length)%buttons.length;buttons[j].focus();fn(buttons[j]);});}
 keyboardTabs(tabs,b=>setView(b.dataset.page));const setControl=b=>{for(const button of controlTabs.children){const active=button===b;button.setAttribute('aria-selected',active);button.tabIndex=active?0:-1;$(button.getAttribute('aria-controls')).hidden=!active;}};controlTabs.addEventListener('click',e=>{if(e.target.dataset.control)setControl(e.target);});keyboardTabs(controlTabs,setControl);
 for(const id of ['open-controls','close-controls'])$(id).onclick=()=>grid.classList.toggle('controls-open',id==='open-controls');$('run').addEventListener('click',()=>grid.classList.remove('controls-open'));
 $('show-geometry').onclick=()=>{plots.classList.remove('show-history');$('show-geometry').setAttribute('aria-selected','true');$('show-histories').setAttribute('aria-selected','false');ExposureViewer.resize();};$('show-histories').onclick=()=>{plots.classList.add('show-history');$('show-histories').setAttribute('aria-selected','true');$('show-geometry').setAttribute('aria-selected','false');};
 function sceneMode(dimension){const entering3D=dimension===3&&$('three-view').hidden;$('blender-view').hidden=dimension!=='render';$('dilution-overlay').toggleAttribute('hidden',dimension!=='render');stationControls.hidden=dimension!=='render';$('view-render').setAttribute('aria-selected',dimension==='render');plotWrap.hidden=dimension!==2;$('three-view').hidden=dimension!==3;$('viewer-provenance').hidden=dimension===2;viewtools.hidden=dimension!==3;$('view-3d').setAttribute('aria-selected',dimension===3);$('view-2d').setAttribute('aria-selected',dimension===2);viewerNote.hidden=dimension===2;if(entering3D)ExposureViewer.update(current,current?.g||getGeometry());ExposureViewer.resize();} $('view-render').onclick=()=>{sceneMode('render');renderCaption();};$('view-3d').onclick=()=>sceneMode(3);$('view-2d').onclick=()=>{sceneMode(2);drawField();};
 function renderCaption(){const g=getGeometry();$('blender-view').src=BLENDER_RENDERS[g.key]||'';$('blender-view').alt='Blender render · '+g.name;$('viewer-provenance').textContent='Blender / Cycles · '+(g.key.startsWith('nov')?'complete STEP exterior':'analytical domain');viewerNote.textContent=g.key.startsWith('nov')?'Complete CAD assembly. Numbered markers project internal model sections onto the exterior. Use 3D → Cutaway to see inside; particle histories use 3D passage paths with approximate junctions.':g.extensionMM?'Nominal Blender image; the added '+g.extensionMM+' mm tube is shown in 3D. Select a numbered location below for its calculated values.':'Numbered dots mark shear, dilution and coagulation sections; click a number to inspect its values. This is the analytical particle domain.';}
 $('play').addEventListener('click',()=>sceneMode(3));$('time-control').addEventListener('input',()=>sceneMode(3));sceneMode('render');renderCaption();
 pageEls.compressible.innerHTML=`<div class="page-head"><div><h2>Automatic compressible estimation</h2><p>The active flow balance is selected automatically. The ideal nozzle example remains separate.</p></div></div><section class="card" id="active-flow-network" aria-live="polite"></section><div class="compress-layout"><section class="card"><h3>Backpressure selects the regime</h3><p class="note">Uses the selected motive diameter, reservoir pressure and temperature. This illustrative C–D nozzle is not the eductor insert or a solved eductor flow.</p><div class="compress-controls"><div class="field"><label for="comp-back">Local nozzle backpressure · bar abs</label><input id="comp-back" type="number" min="0.001" step="any" value="1"></div><div class="field"><label for="comp-area">Demonstration exit / throat area</label><input id="comp-area" type="number" min="1" max="20" step="0.1" value="2"></div><div class="field"><label for="comp-sweep">Backpressure / total pressure</label><input id="comp-sweep" type="range" min="0.02" max="0.999" step="0.001" value="0.333"></div></div><div id="comp-gate" class="notice" role="status"></div><div id="comp-regime" class="regime" role="status"></div><svg id="comp-chart" class="compress-chart" viewBox="0 0 720 340" role="img" aria-label="Quasi-one-dimensional nozzle Mach and static pressure profiles"></svg><div id="comp-kpis" class="compress-kpis"></div><p id="comp-description" class="note"></p><p class="note">Ideal gas γ = 1.4, adiabatic, frictionless except for one normal shock; discharge coefficient = 1. Pressure thrust and external adjustment are not modeled here. The displayed nozzle shape is illustrative.</p><div class="reference-strip"><a href="https://www.grc.nasa.gov/www/k-12/airplane/isentrop.html" target="_blank" rel="noopener">NASA · isentropic flow</a><a href="https://www.grc.nasa.gov/www/k-12/airplane/normal.html" target="_blank" rel="noopener">NASA · normal shocks</a></div></section><section class="card"><h3>What the automatic estimate solves</h3><div class="network"><span>Motive reservoir → nozzle</span><span>Shared vessel → suction</span><span class="wide">Mass + axial momentum + total energy</span><span class="wide">Outlet → shared boundary pressure</span></div><ol class="compact-steps"><li><strong>Only when needed.</strong> The reduced balance is tried first. Bulk Mach ≥ 0.3, significant pressure variation or a failed validity check activates compressible mixing. Motive-nozzle choking alone does not activate the extra model.</li><li><strong>Local inlet restrictions.</strong> Ideal gas inlet laws cap motive, secondary and bleed flows at their sonic limits. Effective areas account for assumed branch and feeder losses; the insert restricts the secondary passage.</li><li><strong>Coupled balance.</strong> Solve chamber pressure, gas flows and mixed total temperature together. Axial momentum includes pressure thrust, powder inertia and wall loss. Reject nonphysical temperatures, entropy decrease and unbracketed roots.</li><li><strong>Match the outlet.</strong> The subsonic area-change branch must match the shared external pressure. If it cannot and the bore has an expanding diffuser, try critical mixed flow (M = 1) and solve an internal normal-shock area against that same outlet pressure. For two stages, transfer mass and total enthalpy through a calculated interstage pressure.</li><li><strong>Update the plots.</strong> Calculated air flows and local density feed dilution, drag and exposure automatically. No manual air-flow substitution is made.</li></ol><p class="note">This is an uncalibrated complete-mixing estimate, particularly uncertain in short tubes. Bleed is lumped in the balance. An ideal internal diffuser shock is now estimated for the critical branch. External shocks, compound choking, separation, jet blockage and reverse-flow powder routing remain unresolved. Particle histories use a smooth analytical pressure field; they do not resolve shock-particle dynamics.</p><div class="reference-strip"><a href="https://www.grc.nasa.gov/www/k-12/airplane/mflchk.html" target="_blank" rel="noopener">NASA · inlet mass flow</a><a href="https://doi.org/10.1016/j.energy.2022.126551" target="_blank" rel="noopener">1D ejector network research</a></div></section></div>`;
 function drawCompressible(){try{const f=Exposure.flow(getGeometry(),settings(),Legacy);$('active-flow-network').innerHTML='<h3>'+(f.ok?(f.compressible?'Compressible mixing estimate active':f.manual?'Prescribed airflow scenario':'Low-Mach balance active'):'No admissible flow estimate')+'</h3><p class="note">'+(f.gate?.reason||f.reason)+'</p>'+(f.ok?'<div class="compress-kpis">'+[['Incoming air',fmt(f.ms*1000,2)+' g/s'],['Chamber pressure',fmt(f.portP/1e5,4)+' bar abs'],['Outlet air',fmt(f.mg*1000,2)+' g/s']].map(([a,b])=>'<div>'+a+'<strong>'+b+'</strong></div>').join('')+'</div>':'')+(f.network?'<p class="note"><b>'+f.network.regime+'</b>'+(f.network.shock?' · shock Mach '+fmt(f.network.shock.before.M,2)+' → '+fmt(f.network.shock.after.M,2)+' · equivalent diameter '+fmt(f.network.shock.diameter*1000,2)+' mm':'')+'</p><p class="note">Outlet residual '+Math.abs(f.network.pressureResidual).toExponential(2)+' Pa · mass residual '+f.network.residuals.mass.toExponential(2)+' kg/s · momentum residual '+f.network.residuals.momentum.toExponential(2)+' N · energy residual '+f.network.residuals.energy.toExponential(2)+' W. These check the equations, not experimental accuracy.</p>':'');const g=getGeometry(),p=settings(),input={p0:p.p0*1e5,pb:+$('comp-back').value*1e5,T0:p.temp+273.15,diameter:g.dn/1000,ratio:+$('comp-area').value},gate=Compressible.assess(input);$('comp-sweep').value=+$('comp-back').value/p.p0;$('comp-gate').textContent=gate.reason+(gate.valid&&!gate.choked?' Estimated peak nozzle Mach: '+fmt(gate.peakMach,3)+'.':'');if(!gate.valid||!gate.needed){$('comp-regime').textContent=gate.valid?'Low-Mach regime · compressible module inactive':'No admissible forward flow';$('comp-chart').innerHTML='';$('comp-kpis').innerHTML='';$('comp-description').textContent='The extended nozzle/shock solver was not called. This demonstration never replaces the current exposure flow model.';return;}const r=Compressible.solve(input);$('comp-sweep').value=+$('comp-back').value/p.p0;if(!r.ok)throw Error(r.reason);$('comp-regime').textContent=r.regime;$('comp-kpis').innerHTML=[['Ideal motive flow',fmt(r.mdot*1000,2)+' g/s'],['Exit Mach',fmt(r.exit.M,2)],['Total pressure retained',fmt(r.p02/(p.p0*1e5)*100,1)+'%']].map(([a,b])=>`<div>${a}<strong>${b}</strong></div>`).join('');
 const xx=x=>60+600*x,maxM=Math.max(1.2,...r.samples.map(s=>s.M))*1.1,ym=m=>150-105*m/maxM,yp=v=>285-90*v;let s='<text x="60" y="22" font-size="15" fill="#18394b">Mach number</text><text x="60" y="179" font-size="15" fill="#18394b">Static pressure / reservoir total pressure</text>';
 for(let i=0;i<=4;i++){let y=ym(i*maxM/4);s+=`<line x1="60" x2="660" y1="${y}" y2="${y}" stroke="#e0e8ed"/><text x="49" y="${y+4}" text-anchor="end" font-size="12" fill="#566773">${fmt(i*maxM/4,1)}</text>`;y=yp(i/4);s+=`<line x1="60" x2="660" y1="${y}" y2="${y}" stroke="#e0e8ed"/><text x="49" y="${y+4}" text-anchor="end" font-size="12" fill="#566773">${fmt(i/4,2)}</text>`;}
 const line=(key,y,color)=>`<path d="${r.samples.map((v,i)=>(i?'L':'M')+xx(v.x).toFixed(2)+','+y(v[key]).toFixed(2)).join('')}" fill="none" stroke="${color}" stroke-width="2.7"/>`;s+=line('M',ym,'#08678b')+line('p',v=>yp(v/(p.p0*1e5)),'#008678');s+=`<line x1="60" x2="660" y1="${ym(1)}" y2="${ym(1)}" stroke="#a76b1f" stroke-dasharray="4 4"/><line x1="${xx(.4)}" x2="${xx(.4)}" y1="35" y2="285" stroke="#9badb8" stroke-dasharray="3 4"/><text x="${xx(.4)}" y="310" text-anchor="middle" font-size="13">Throat</text><text x="60" y="310" font-size="13">Inlet</text><text x="660" y="310" text-anchor="end" font-size="13">Exit</text><text x="360" y="334" text-anchor="middle" font-size="12" fill="#566773">Illustrative axial coordinate, not CAD distance</text>`;
 if(r.shockX!==null)s+=`<line x1="${xx(r.shockX)}" x2="${xx(r.shockX)}" y1="35" y2="285" stroke="#bf5015" stroke-width="1.5"/><text x="${xx(r.shockX)}" y="31" text-anchor="middle" font-size="13" fill="#bf5015">Normal shock</text>`;$('comp-chart').innerHTML=s;
 $('comp-description').textContent=`${g.dn} mm throat; ${fmt(p.p0,3)} bar total; ${fmt(p.temp,0)} °C. Choking starts at pb/p0 = ${fmt(r.onset,4)} for this nozzle area ratio; shock-at-exit boundary = ${fmt(r.exitShock,4)}. `+(r.shockX!==null?'Internal shock position is solved to match the specified backpressure.':r.choked?'Exit static pressure is '+fmt(r.exit.p/1e5,3)+' bar; any remaining pressure mismatch adjusts outside this ideal nozzle.':'The subsonic exit pressure matches the backpressure.');
 }catch(e){$('comp-regime').textContent=e.message;$('comp-chart').innerHTML='';$('comp-kpis').innerHTML='';$('comp-description').textContent='';}}
 for(const id of ['p0','ps','feed','loss','temp','geometry','flow-mode'])$(id).addEventListener('input',()=>{if(!pageEls.compressible.hidden)drawCompressible();});for(const id of ['comp-back','comp-area'])$(id).addEventListener('input',drawCompressible);$('comp-sweep').oninput=()=>{$('comp-back').value=(+$('comp-sweep').value*+$('p0').value).toFixed(5);drawCompressible();};
 window.exposureRendered=()=>{showRecovery(null);$('keep-experiment').disabled=false;overlayDescription();if(!pageEls.compressible.hidden)drawCompressible();const q=current.r.pressure;$('solved-flow').hidden=q.manual;$('solved-flow-values').innerHTML=[['Incoming air',fmt(q.ms*1000,2)+' g/s'],['Bleed air',fmt(q.mh*1000,2)+' g/s'],['Chamber',fmt(q.portP/1e5,4)+' bar abs']].map(([a,b])=>'<span>'+a+'</span><b>'+b+'</b>').join('');if(q.circuit?.flowUsed&&q.circuit.branches.length>1)$('solved-flow-values').insertAdjacentHTML('beforeend','<small style="grid-column:1/-1">Branch shares: '+q.circuit.branches.map(b=>fmt(b.fraction*100,1)+'%').join(' / ')+'</small>');$('lab-n').value=$('n').value;$('play').disabled=false;status.classList.add('flow-bar');if(current.g.extensionMM)sceneMode(3);ExposureViewer.update(current);if(!$('blender-view').hidden)renderCaption();ExposureAnalytics.render(current);ExposureDilution.render(current);showTrackingNotice();ExposureSnapshots.refresh();};window.exposureUnavailable=(g,issue)=>{showRecovery(issue);$('keep-experiment').disabled=true;$('p0').parentElement.querySelector('small').hidden=!g.key.startsWith('revB');$('solved-flow-values').textContent=issue.code==='particles'?'Gas solved · particle model unavailable':'No supported forward flow';ExposureAnalytics.clear();ExposureDilution.clear();results.hidden=false;$('play').disabled=true;$('diagnostics').textContent='';$('path-legend').innerHTML='';$('geometry-source').textContent=g.source;$('result-title').textContent=g.name;$('run-caption').textContent=issue.title+' · geometry remains available';kpis.innerHTML='';for(const id of ['history-small','history-large','lab-coag-small','lab-coag-large'])$(id).innerHTML='<text x="50%" y="50%" text-anchor="middle" font-size="12" fill="#647785">No valid flow solution</text>';$('paths').innerHTML='';const ctx=$('field').getContext('2d');ctx.clearRect(0,0,900,510);ctx.font='20px system-ui';ctx.fillStyle='#647785';ctx.fillText('No flow solution. Inspect the 3D geometry.',180,250);$('metrics').innerHTML='<tr><td colspan="11">No valid flow solution for the selected inputs.</td></tr>';$('benchmark').innerHTML='';$('benchmark-status').textContent='No current particle result.';status.classList.add('flow-bar');ExposureViewer.update(null,g);if(!$('blender-view').hidden)renderCaption();viewerNote.textContent='Geometry reference only. '+issue.title+'.';ExposureSnapshots.refresh();};
 ExposureSnapshots.init();window.exposureStationSelected=()=>{if(current)drawHistories();};ExposureAnalytics.init(pageEls);ExposureDilution.init();$('smallShare').closest('.field').hidden=$('background').value!=='mixture';window.exposureTime=t=>ExposureViewer.time(t);ExposureViewer.init();ExposureViewer.update(null,getGeometry());renderCaption();new ResizeObserver(()=>{if(current)drawHistories();}).observe(histBox);$('geometry').addEventListener('change',run);
 const requestedView=new URLSearchParams(location.search).get('view');if(pages.some(([id])=>id===requestedView))setView(requestedView);
 const info=make('dialog','info-dialog','<button id="close-info">Close</button><h2>Model assumptions</h2><p id="info-text"></p>');document.body.append(info);info.setAttribute('aria-label','Model assumptions');status.onclick=e=>{if(!e.target.closest('.flow-details'))return;$('info-text').textContent=status.dataset.details||status.textContent;info.showModal();};$('close-info').onclick=()=>info.close();info.addEventListener('click',e=>{if(e.target===info)info.close();});
})();
