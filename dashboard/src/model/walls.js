/* STEP triangle signed distance, with a lazy 0.1 mm trilinear grid (0.05 mm for the short dual-suction case). The native
 * surface supplies distance/sign; the existing inlet/outlet planes remain open.
 * No invalid STEP is auto-repaired. Fluid-volume certification is separate. */
const EductorWalls=(()=>{
 const DATA=typeof EDUCTOR_WALLS!=='undefined'?EDUCTOR_WALLS:typeof require==='function'?require('../phases/walls.json'):{};
 const cache=new Map(),sub=(a,b)=>a.map((v,i)=>v-b[i]),dot=(a,b)=>a.reduce((s,v,i)=>s+v*b[i],0),cross=(a,b)=>[a[1]*b[2]-a[2]*b[1],a[2]*b[0]-a[0]*b[2],a[0]*b[1]-a[1]*b[0]];
 const bytes=s=>typeof Buffer!=='undefined'?new Uint8Array(Buffer.from(s,'base64')):Uint8Array.from(atob(s),c=>c.charCodeAt(0));
 function nearest(p,a,b,c){
  const ab=sub(b,a),ac=sub(c,a),ap=sub(p,a),d1=dot(ab,ap),d2=dot(ac,ap);if(d1<=0&&d2<=0)return a;
  const bp=sub(p,b),d3=dot(ab,bp),d4=dot(ac,bp);if(d3>=0&&d4<=d3)return b;
  const vc=d1*d4-d3*d2;if(vc<=0&&d1>=0&&d3<=0)return a.map((v,i)=>v+d1/(d1-d3)*ab[i]);
  const cp=sub(p,c),d5=dot(ab,cp),d6=dot(ac,cp);if(d6>=0&&d5<=d6)return c;
  const vb=d5*d2-d1*d6;if(vb<=0&&d2>=0&&d6<=0)return a.map((v,i)=>v+d2/(d2-d6)*ac[i]);
  const va=d3*d6-d5*d4;if(va<=0&&d4-d3>=0&&d5-d6>=0){const w=(d4-d3)/(d4-d3+d5-d6);return b.map((v,i)=>v+w*(c[i]-v));}
  const v=vb/(va+vb+vc),w=vc/(va+vb+vc);return a.map((q,i)=>q+ab[i]*v+ac[i]*w);
 }
 function build(key){
  if(cache.has(key))return cache.get(key);const data=DATA[key];if(!data)return null;
  const vb=bytes(data.v),tb=bytes(data.t),vs=new Uint16Array(vb.buffer,vb.byteOffset,vb.byteLength/2),ts=new Uint32Array(tb.buffer,tb.byteOffset,tb.byteLength/4),v=Array.from({length:vs.length/3},(_,i)=>[0,1,2].map(k=>(data.lo[k]+data.span*vs[i*3+k]/65535)/1000));
  const tris=Array.from({length:ts.length/3},(_,i)=>{const a=v[ts[i*3]],b=v[ts[i*3+1]],c=v[ts[i*3+2]],n=cross(sub(b,a),sub(c,a)),length=Math.hypot(...n);return {a,b,c,n:n.map(v=>v/(length||1)),lo:a.map((q,k)=>Math.min(q,b[k],c[k])),hi:a.map((q,k)=>Math.max(q,b[k],c[k]))};});
  function node(list){const lo=[Infinity,Infinity,Infinity],hi=[-Infinity,-Infinity,-Infinity];for(const t of list)for(let k=0;k<3;k++){lo[k]=Math.min(lo[k],t.lo[k]);hi[k]=Math.max(hi[k],t.hi[k]);}if(list.length<=10)return {lo,hi,list};const span=hi.map((q,k)=>q-lo[k]),axis=span.indexOf(Math.max(...span));list.sort((a,b)=>a.lo[axis]+a.hi[axis]-b.lo[axis]-b.hi[axis]);const mid=list.length>>1;return {lo,hi,left:node(list.slice(0,mid)),right:node(list.slice(mid))};}
  const root=node(tris),grid=new Map(),h=key==='nov5'?.00005:.0001;
  const box=(p,n)=>p.reduce((s,q,i)=>s+(q<n.lo[i]?n.lo[i]-q:q>n.hi[i]?q-n.hi[i]:0)**2,0);
  function exact(p){let best=Infinity,sign=1,normal=[0,0,0];function search(n){if(box(p,n)>best)return;if(n.list){for(const t of n.list){const q=nearest(p,t.a,t.b,t.c),d=sub(p,q),r=dot(d,d);if(r<best){best=r;sign=dot(d,t.n)>=0?1:-1;normal=t.n;}}}else{const l=box(p,n.left),r=box(p,n.right);if(l<r){search(n.left);search(n.right);}else{search(n.right);search(n.left);}}}search(root);
   // Oriented ray winding gives the sign of the UNION of closed hardware
   // solids. A nearest-face normal alone is discontinuous at concave edges
   // and at assembly-part overlaps; it is not a valid signed-distance sign.
   const direction=[.9238795325,.331413574,.191341716];let winding=0;
   function raybox(n){let lo=0,hi=Infinity;for(let i=0;i<3;i++){const a=(n.lo[i]-p[i])/direction[i],b=(n.hi[i]-p[i])/direction[i];lo=Math.max(lo,Math.min(a,b));hi=Math.min(hi,Math.max(a,b));}return hi>=lo;}
   function ray(n){if(!raybox(n))return;if(n.list){for(const t of n.list){const e1=sub(t.b,t.a),e2=sub(t.c,t.a),h=cross(direction,e2),det=dot(e1,h);if(Math.abs(det)<1e-18)continue;const q=sub(p,t.a),u=dot(q,h)/det;if(u<0||u>1)continue;const v=cross(q,e1),w=dot(direction,v)/det;if(w<0||u+w>1)continue;const distance=dot(e2,v)/det;if(distance>1e-10)winding+=dot(direction,t.n)>0?1:-1;}}else{ray(n.left);ray(n.right);}}
   ray(root);sign=winding!==0?-1:1;return {distance:sign*Math.sqrt(best),normal};}
  function distance(x,y,z){const q=[x/h,y/h,z/h],i=q.map(Math.floor),f=q.map((v,k)=>v-i[k]);let value=0,reference=null,sharp=false;for(let mask=0;mask<8;mask++){const j=i.map((v,k)=>v+((mask>>k)&1)),key=j.join(','),weight=f.reduce((s,v,k)=>s*(((mask>>k)&1)?v:1-v),1);let d=grid.get(key);if(d===undefined){d=exact(j.map(v=>v*h));if(grid.size>=100000)grid.delete(grid.keys().next().value);grid.set(key,d);}if(reference&&dot(reference,d.normal)<.99)sharp=true;reference??=d.normal;value+=weight*d.distance;}return sharp?exact([x,y,z]).distance:value;}
  function normal(x,y,z){const q=[x,y,z],d=q.map((_,i)=>{const a=q.slice(),b=q.slice();a[i]-=h/2;b[i]+=h/2;return distance(...b)-distance(...a);}),n=Math.hypot(...d);return d.map(v=>v/(n||1));}
  const out={distance,exact,normal,gridSpacing:h,source:data,triangles:tris.length};if(cache.size>=2)cache.delete(cache.keys().next().value);cache.set(key,out);return out;
 }
 function attach(geo,key){const sdf=build(key);if(!sdf)return {...geo,wallStatus:'No valid native STEP SDF; reduced walls retained.'};
  const clearance=(x,y,dp=0,z=0)=>{const inlet=Math.max(...geo.branches.map(b=>{const a=b.segments[0];return dot(sub([x,y,z],a.a),a.tangent);}));return Math.min(sdf.distance(x,y,z),geo.end-x,inlet)-dp/2;};
  return {...geo,clearance,sdf,wallStatus:'Native STEP surface distance · '+(sdf.gridSpacing*1000)+' mm grid · fluid-volume gate pending',wallNormal:sdf.normal};
 }
 return {build,attach,nearest,register:(key,data)=>{DATA[key]=data;cache.delete(key);},has:key=>!!DATA[key]};
})();
if(typeof module!=='undefined')module.exports=EductorWalls;

'use strict';
// Predictor/corrector pseudo-arclength continuation in scaled (P, boundary P).
