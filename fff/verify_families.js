"use strict";
// Verify every affine component of every retained support representative.
const fs=require("node:fs"),path=require("node:path");
const {matrices,snf,mm,seed}=require("./fff.js");
const ROOT=process.argv[2]||__dirname;
const supports=JSON.parse(fs.readFileSync(path.join(ROOT,"fff_supports.json"))),classes=JSON.parse(fs.readFileSync(path.join(ROOT,"fff_classes.json"))).classes;
function assert(x,m){if(!x)throw Error(m);}
function branches(d){const out=[];function rec(a){if(a.length===d.length){out.push(a.slice());return;}for(let k=0;k<d[a.length];k++){a.push(k);rec(a);a.pop();}}rec([]);return out;}
const report={verified:true,support_classes:classes.length,orientations:0,components:0,by_dimension:{},seed_orders:{},examples:{}},certificates=[];
for(let z=0;z<classes.length;z++){
 const c=classes[z],q=supports.orbits[c.vi].representative,r=supports.partitions[c.wi],mat=matrices(q,r),s=snf(mat.A),uf=mm(s.U,mat.f.map(x=>[x])).map(x=>x[0]),ub=mm(s.U,mat.B),bb=branches(s.diag),allowed=[];
 assert(8-s.rank===c.dimension,"dimension");assert(s.diag.reduce((a,b)=>a*b,1)===c.torsion,"torsion");assert(uf.slice(s.rank).every(x=>x%2===0),"binary");
 for(let mask=0;mask<8;mask++){
  const sig=Array.from({length:3},(_,j)=>mask>>j&1?1:-1);
  if(!ub.slice(s.rank).every(row=>row.reduce((t,x,j)=>t+x*sig[j],0)%3===0))continue;
  allowed.push(mask);report.orientations++;
  for(const branch of bb){
   const ex=seed(q,r,s,mat,mask,branch);
   assert(24%ex.cross_ratio_order===0,"seed order divides 24");
   report.components++;report.by_dimension[c.dimension]=(report.by_dimension[c.dimension]||0)+1;
   report.seed_orders[ex.cross_ratio_order]=(report.seed_orders[ex.cross_ratio_order]||0)+1;
   if(!report.examples["d"+c.dimension])report.examples["d"+c.dimension]={class_index:z,vi:c.vi,wi:c.wi,...ex};
   if(!report.examples["r"+ex.cross_ratio_order])report.examples["r"+ex.cross_ratio_order]={class_index:z,vi:c.vi,wi:c.wi,...ex};
  }
 }
 assert(allowed.length===c.orientations,"orientations");
 certificates.push({class_index:z,vi:c.vi,wi:c.wi,rank:s.rank,diagonal:s.diag,U:s.U,V:s.V,orientations:allowed});
 if((z+1)%2000===0)console.error("Verified class families",z+1);
}
fs.writeFileSync(path.join(ROOT,"family_verification.json"),JSON.stringify(report,null,2)+"\n");
fs.writeFileSync(path.join(ROOT,"affine_certificates.json"),JSON.stringify(certificates)+"\n");
console.log(JSON.stringify({...report,examples:Object.keys(report.examples)},null,2));

