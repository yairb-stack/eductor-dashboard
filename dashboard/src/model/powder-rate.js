var PowderRate=(()=>{
 const exact=new WeakMap(),fromModel=value=>value/3.6,toModel=value=>value*3.6;
 function write(input,kgPerHour){input.value=String(+fromModel(kgPerHour).toFixed(6));exact.set(input,{value:input.value,kgPerHour});}
 function read(input){const saved=exact.get(input);return saved?.value===input.value?saved.kgPerHour:input.value===''?NaN:toModel(Number(input.value));}
 function init(input,kgPerHour=6){write(input,kgPerHour);input.addEventListener('input',()=>exact.delete(input));}
 return {fromModel,toModel,read,write,init};
})();
if(typeof module!=='undefined')module.exports=PowderRate;
