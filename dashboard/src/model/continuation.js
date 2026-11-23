const EductorContinuation=(()=>{
 function trace(residual,a,b,{step=.005,count=80,tolerance=1e-9}={}){
  const points=[a.slice(),b.slice()];let ds=step;
  for(let n=0;n<count;n++){
   const u=points.at(-2),v=points.at(-1),delta=v.map((x,i)=>x-u[i]),norm=Math.hypot(...delta);if(!norm)break;
   const t=delta.map(x=>x/norm),pred=v.map((x,i)=>x+ds*t[i]);let x=pred.slice(),good=false;
   for(let k=0;k<16;k++){
    const f=residual(...x),arc=(x[0]-pred[0])*t[0]+(x[1]-pred[1])*t[1];if(!Number.isFinite(f))break;
    if(Math.abs(f)<tolerance&&Math.abs(arc)<tolerance){good=true;break;}
    const h=1e-6,fx=(residual(x[0]+h,x[1])-residual(x[0]-h,x[1]))/(2*h),fy=(residual(x[0],x[1]+h)-residual(x[0],x[1]-h))/(2*h),det=fx*t[1]-fy*t[0];if(!Number.isFinite(det)||Math.abs(det)<1e-12)break;
    x=[x[0]-(f*t[1]-fy*arc)/det,x[1]-(fx*arc-f*t[0])/det];
   }
   if(good){points.push(x);ds=Math.min(step,ds*1.4);}else{ds*=.5;if(ds<step/128)break;n--;}
  }
  return points;
 }
 return {trace};
})();
if(typeof module!=='undefined')module.exports=EductorContinuation;

'use strict';
