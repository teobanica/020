"use strict";
// Exact chosen-support reduction for F=3+2+2+2+2, width 11.
// Node.js, no external packages. Run: node fff.js [output-directory].
const fs=require("node:fs"),path=require("node:path");
const OUT=process.argv[2]||__dirname;
fs.mkdirSync(OUT,{recursive:true});
function assert(x,msg="assertion"){if(!x)throw Error(msg);}
const SAFE={maxTerm:0,maxU:0,maxV:0,matrices:0};
function checked(x){assert(Number.isSafeInteger(x),"unsafe integer");SAFE.maxTerm=Math.max(SAFE.maxTerm,Math.abs(x));return x;}
function combinations(a,k){const out=[];function rec(s,p){if(p.length===k){out.push(p.slice());return;}for(let i=s;i<=a.length-k+p.length;i++){p.push(a[i]);rec(i+1,p);p.pop();}}rec(0,[]);return out;}
function matchings(a){if(!a.length)return [[]];const out=[];for(let i=1;i<a.length;i++){const b=[a[0],a[i]];for(const r of matchings(a.slice(1).filter(x=>x!==a[i])))out.push([b,...r]);}return out;}
function cmp(a,b){for(let i=0;i<a.length;i++){if(Array.isArray(a[i])){const c=cmp(a[i],b[i]);if(c)return c;}else if(a[i]!==b[i])return a[i]-b[i];}return 0;}
function key(p){return p.map(b=>b.join(",")).join("|");}
const BASE=[[0,1,2],[3,4],[5,6],[7,8],[9,10]];
function partitions(){const all=Array.from({length:11},(_,i)=>i),out=[];for(const tri of combinations(all,3)){const s=new Set(tri);for(const m of matchings(all.filter(x=>!s.has(x))))out.push([tri,...m]);}return out.sort(cmp);}
function action(p,g){return [p[0].map(x=>g[x]).sort((a,b)=>a-b),...p.slice(1).map(b=>b.map(x=>g[x]).sort((a,b)=>a-b)).sort(cmp)];}
function orbits(parts){
 const index=new Map(parts.map((p,i)=>[key(p),i])),gens=[];
 for(const [a,b] of [[0,1],[1,2],[3,4],[5,6],[7,8],[9,10]]){const g=Array.from({length:11},(_,i)=>i);[g[a],g[b]]=[g[b],g[a]];gens.push(g);}
 for(let k=0;k<3;k++){const g=Array.from({length:11},(_,i)=>i),a=3+2*k;[g[a],g[a+2]]=[g[a+2],g[a]];[g[a+1],g[a+3]]=[g[a+3],g[a+1]];gens.push(g);}
 const seen=new Uint8Array(parts.length),out=[];
 for(let i=0;i<parts.length;i++)if(!seen[i]){
  seen[i]=1;const todo=[i];
  for(let h=0;h<todo.length;h++)for(const g of gens){const j=index.get(key(action(parts[todo[h]],g)));assert(j!==undefined,"action");if(!seen[j]){seen[j]=1;todo.push(j);}}
  out.push({index:i,representative:parts[i],size:todo.length});
 }return out;
}
function maps(p){const block=Array(11),pos=Array(11),eps=Array(11).fill(0);p.forEach((b,k)=>b.forEach((j,h)=>{block[j]=k;pos[j]=h;if(k)eps[j]=h;}));return {block,pos,eps};}
function matrices(q,r){
 const p=maps(BASE),qq=maps(q),rr=maps(r),w=[],off=[];
 for(let j=0;j<11;j++){
  const a=Array(8).fill(0),b=Array(3).fill(0);
  if(p.block[j])a[p.block[j]-1]=1;else b[0]=-p.pos[j];
  if(qq.block[j])a[qq.block[j]+3]=-1;else b[1]=qq.pos[j];
  if(!rr.block[j])b[2]=rr.pos[j];
  w.push(a);off.push(b);
 }
 const h=Array.from({length:11},(_,j)=>rr.eps[j]-p.eps[j]+qq.eps[j]),A=[],f=[],B=[];
 for(const block of r){const ref=block[0];for(const j of block.slice(1)){A.push(w[j].map((x,i)=>x-w[ref][i]));f.push(h[j]-h[ref]);B.push(off[j].map((x,i)=>x-off[ref][i]));}}
 return {A,f,B};
}
function ident(n){return Array.from({length:n},(_,i)=>Array.from({length:n},(_,j)=>+(i===j)));}
function rowadd(a,i,j,q){for(let k=0;k<a[0].length;k++)a[i][k]=checked(a[i][k]+checked(q*a[j][k]));}
function coladd(a,i,j,q){for(const r of a)r[i]=checked(r[i]+checked(q*r[j]));}
function colswap(a,i,j){for(const r of a)[r[i],r[j]]=[r[j],r[i]];}
function snf(input){
 const m=input.length,n=input[0].length,a=input.map(r=>r.slice()),u=ident(m),v=ident(n);let k=0;
 for(;k<Math.min(m,n);k++){
  let pi=-1,pj=-1;
  for(let i=k;i<m;i++)for(let j=k;j<n;j++)if(a[i][j]&&(pi<0||Math.abs(a[i][j])<Math.abs(a[pi][pj]))){pi=i;pj=j;}
  if(pi<0)break;[a[k],a[pi]]=[a[pi],a[k]];[u[k],u[pi]]=[u[pi],u[k]];colswap(a,k,pj);colswap(v,k,pj);
  for(let loop=0;;loop++){
   assert(loop<500,"Smith loop");let changed=false;
   for(let i=k+1;i<m;i++)if(a[i][k]){
    const q=Math.trunc(a[i][k]/a[k][k]);rowadd(a,i,k,-q);rowadd(u,i,k,-q);
    if(a[i][k]){[a[i],a[k]]=[a[k],a[i]];[u[i],u[k]]=[u[k],u[i]];}changed=true;break;
   }if(changed)continue;
   for(let j=k+1;j<n;j++)if(a[k][j]){
    const q=Math.trunc(a[k][j]/a[k][k]);coladd(a,j,k,-q);coladd(v,j,k,-q);
    if(a[k][j]){colswap(a,j,k);colswap(v,j,k);}changed=true;break;
   }if(changed)continue;
   let bi=-1;for(let i=k+1;i<m;i++)for(let j=k+1;j<n;j++)if(a[i][j]%a[k][k])bi=i;
   if(bi<0)break;rowadd(a,k,bi,1);rowadd(u,k,bi,1);
  }
  if(a[k][k]<0){a[k]=a[k].map(x=>-x);u[k]=u[k].map(x=>-x);}
 }
 return {D:a,U:u,V:v,rank:k,diag:Array.from({length:k},(_,i)=>a[i][i])};
}
function det(input){
 const a=input.map(r=>r.slice()),n=a.length;let old=1,sign=1;
 for(let k=0;k<n-1;k++){
  if(!a[k][k]){let j=k+1;while(j<n&&!a[j][k])j++;if(j===n)return 0;[a[k],a[j]]=[a[j],a[k]];sign=-sign;}
  const z=a[k][k];for(let i=k+1;i<n;i++)for(let j=k+1;j<n;j++){
   const num=checked(checked(z*a[i][j])-checked(a[i][k]*a[k][j]));assert(num%old===0,"Bareiss divisibility");a[i][j]=checked(num/old);
  }for(let i=k+1;i<n;i++)a[i][k]=0;old=z;
 }return checked(sign*a[n-1][n-1]);
}
function mm(a,b){const out=Array.from({length:a.length},()=>Array(b[0].length).fill(0));for(let i=0;i<a.length;i++)for(let k=0;k<b.length;k++)if(a[i][k])for(let j=0;j<b[0].length;j++)out[i][j]=checked(out[i][j]+checked(a[i][k]*b[k][j]));return out;}
function mod(x,n){return (x%n+n)%n;}
function gcd(a,b){return b?gcd(b,a%b):a;}
function lcm(a,b){return a/gcd(a,b)*b;}
function popcount(x){let n=0;while(x){n+=x&1;x>>>=1;}return n;}
function add(obj,k,n=1){obj[k]=(obj[k]||0)+n;}
function branches(diag){const out=[];function rec(a){if(a.length===diag.length){out.push(a.slice());return;}for(let k=0;k<diag[a.length];k++){a.push(k);rec(a);a.pop();}}rec([]);return out;}
function seed(q,r,s,mat,mask,branch){
 const M=6*(s.diag.at(-1)||1),sig=Array.from({length:3},(_,i)=>mask>>i&1?1:-1);
 const rhs=mat.f.map((f,i)=>3*f+2*mat.B[i].reduce((t,x,k)=>t+x*sig[k],0)),ur=mm(s.U,rhs.map(x=>[x])).map(x=>x[0]);
 const theta=Array.from({length:8},(_,i)=>mod(s.diag.reduce((t,d,j)=>t+s.V[i][j]*(M/(6*d))*(ur[j]+6*branch[j]),0),M));
 const rows=[],basis=[];
 for(const [part,start,si] of [[BASE,0,0],[q,4,1]]){
  const m=maps(part),e=[],b=[];
  for(let j=0;j<11;j++){
   e.push(mod(m.block[j]?theta[start+m.block[j]-1]+M/2*m.pos[j]:M/3*m.pos[j]*sig[si],M));
   b.push(m.block[j]?s.V[start+m.block[j]-1].slice(s.rank):Array(8-s.rank).fill(0));
  }rows.push(e);basis.push(b);
 }
 const w=rows[0].map((x,j)=>mod(x-rows[1][j],M)),wb=basis[0].map((x,j)=>x.map((y,k)=>y-basis[1][j][k]));
 for(const [e,b,p] of [[rows[0],basis[0],BASE],[rows[1],basis[1],q],[w,wb,r]])for(const block of p){
  const first=block[0],delta=block.map(j=>mod(e[j]-e[first],M)).sort((a,b)=>a-b);
  assert(delta.every((v,k)=>v===k*M/block.length),"seed cycle");
  assert(block.every(j=>b[j].every((v,k)=>v===b[first][k])),"affine cycle");
 }
 let order=1;for(const row of rows)for(const x of row)order=lcm(order,M/gcd(M,mod(x-row[0],M)));
 return {modulus:M,exponents:rows,basis,supports:[BASE,q,r],cross_ratio_order:order,orientation:mask,branch};
}
function tensor(q,r){const p=maps(BASE).block,qm=maps(q).block,rm=maps(r).block,t=Array(125).fill(0);for(let j=0;j<11;j++)t[25*p[j]+5*qm[j]+rm[j]]++;return t.join("");}
function run(){
 const parts=partitions(),reps=orbits(parts);assert(parts.length===17325);assert(reps.reduce((t,r)=>t+r.size,0)===parts.length);
 const state={partitions:parts.length,relative_orbits:reps.length,presentations:parts.length*reps.length,binary:0,feasible:0,orientations:0,components:0,smith:{},byDim:{},supportDim:{},maxd:{},orders:{},orderByMax:{},examples:{}},raw=new Map(),records=[];
 console.error("Partitions",parts.length,"relative orbits",reps.length);
 for(let vi=0;vi<reps.length;vi++){
  const q=reps[vi].representative;
  for(let wi=0;wi<parts.length;wi++){
   const r=parts[wi],mat=matrices(q,r),s=snf(mat.A);
   assert(Math.abs(det(s.U))===1&&Math.abs(det(s.V))===1,"unimodular");
   const check=mm(mm(s.U,mat.A),s.V);
   for(let i=0;i<6;i++)for(let j=0;j<8;j++)assert(check[i][j]===(i===j&&i<s.rank?s.diag[i]:0),"Smith identity");
   assert(s.diag.every((x,i)=>x>0&&(!i||x%s.diag[i-1]===0)),"Smith divisibility");
   SAFE.matrices++;SAFE.maxU=Math.max(SAFE.maxU,...s.U.flat().map(Math.abs));SAFE.maxV=Math.max(SAFE.maxV,...s.V.flat().map(Math.abs));
   add(state.smith,s.diag.join(","));
   const uf=mm(s.U,mat.f.map(x=>[x])).map(x=>x[0]);
   if(!uf.slice(s.rank).every(x=>x%2===0))continue;state.binary++;
   const ub=mm(s.U,mat.B);let allowed=0;
   for(let mask=0;mask<8;mask++){const sig=Array.from({length:3},(_,j)=>mask>>j&1?1:-1);if(ub.slice(s.rank).every(row=>row.reduce((t,x,j)=>t+x*sig[j],0)%3===0))allowed|=1<<mask;}
   if(!allowed)continue;
   const dim=8-s.rank,tor=s.diag.reduce((a,b)=>a*b,1),maxd=s.diag.at(-1)||1,n=popcount(allowed);
   state.feasible++;state.orientations+=n;state.components+=n*tor;
   add(state.byDim,dim,n*tor);add(state.supportDim,dim);add(state.maxd,maxd);
   records.push([vi,wi,dim,maxd,tor,allowed]);
   const t=tensor(q,r),old=raw.get(t),inv=[dim,maxd,tor,n].join(",");
   if(old){assert(old.invariant===inv,"raw invariant");old.presentations++;}
   else raw.set(t,{vi,wi,dimension:dim,largest_smith:maxd,torsion:tor,orientations:n,presentations:1,invariant:inv});
   const bb=branches(s.diag);
   for(let mask=0;mask<8;mask++)if(allowed>>mask&1){
    const sig=Array.from({length:3},(_,j)=>mask>>j&1?1:-1),ug=ub.map(row=>row.reduce((t,x,j)=>t+x*sig[j],0));
    for(const b of bb){
     let order=6;
     for(let j=0;j<s.rank;j++){const den=6*s.diag[j],num=3*uf[j]+2*ug[j]+6*b[j];order=lcm(order,den/gcd(den,mod(num,den)));}
     add(state.orders,order);state.orderByMax[maxd]??={};add(state.orderByMax[maxd],order);
     if(!state.examples["d"+dim]||!state.examples["r"+order]){
      const ex={vi,wi,dimension:dim,...seed(q,r,s,mat,mask,b)};assert(ex.cross_ratio_order===order,"root formula");
      state.examples["d"+dim]??=ex;state.examples["r"+order]??=ex;
     }
    }
   }
  }
  console.error("Orbit",vi+1,"/",reps.length,"feasible",state.feasible);
 }
 state.integer_audit=SAFE;state.raw_tensors=raw.size;
 fs.writeFileSync(path.join(OUT,"fff_summary.json"),JSON.stringify(state,null,2)+"\n");
 fs.writeFileSync(path.join(OUT,"fff_supports.json"),JSON.stringify({base:BASE,partitions:parts,orbits:reps})+"\n");
 fs.writeFileSync(path.join(OUT,"fff_feasible.json"),JSON.stringify(records)+"\n");
 fs.writeFileSync(path.join(OUT,"fff_raw.json"),JSON.stringify(Array.from(raw))+"\n");
 fs.writeFileSync(path.join(OUT,"class_input.txt"),Array.from(raw,([t,r])=>[t,r.vi,r.wi,r.dimension,r.largest_smith,r.torsion,r.orientations,r.presentations].join(" ")).join("\n")+"\n");
 console.log(JSON.stringify({...state,examples:Object.keys(state.examples)},null,2));
}
if(require.main===module)run();
module.exports={partitions,orbits,maps,matrices,snf,det,mm,seed,tensor,BASE,SAFE};

