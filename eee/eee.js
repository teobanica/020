"use strict";
// Exact EEE census at width 11. No dependencies. Run: node eee.js --write

function combinations(a,k){
 const out=[]; function rec(start,p){if(p.length===k){out.push(p.slice());return;}for(let i=start;i<=a.length-(k-p.length);i++){p.push(a[i]);rec(i+1,p);p.pop();}}rec(0,[]);return out;
}
function triples(rem){if(!rem.length)return [[]];const out=[];for(const more of combinations(rem.slice(1),2)){const b=[rem[0],...more],s=new Set(b);for(const r of triples(rem.filter(x=>!s.has(x))))out.push([b,...r]);}return out;}
function cmp(a,b){for(let i=0;i<a.length;i++){if(Array.isArray(a[i])){const c=cmp(a[i],b[i]);if(c)return c;}else if(a[i]!==b[i])return a[i]-b[i];}return 0;}
function key(p){return p.map(b=>b.join(",")).join("|");}
function partitions(){const all=Array.from({length:11},(_,i)=>i),out=[];for(const pair of combinations(all,2)){const s=new Set(pair);for(const t of triples(all.filter(x=>!s.has(x))))out.push([pair,...t]);}return out.sort(cmp);}
const BASE=[[0,1],[2,3,4],[5,6,7],[8,9,10]];
function action(p,g){return [p[0].map(x=>g[x]).sort((a,b)=>a-b),...p.slice(1).map(b=>b.map(x=>g[x]).sort((a,b)=>a-b)).sort(cmp)];}
function orbits(parts){
 const index=new Map(parts.map((p,i)=>[key(p),i])),gens=[];
 for(const [a,b] of [[0,1],[2,3],[3,4],[5,6],[6,7],[8,9],[9,10]]){const g=Array.from({length:11},(_,i)=>i);[g[a],g[b]]=[g[b],g[a]];gens.push(g);}
 gens.push([0,1,5,6,7,2,3,4,8,9,10]);gens.push([0,1,2,3,4,8,9,10,5,6,7]);
 const seen=new Uint8Array(parts.length),out=[];
 for(let i=0;i<parts.length;i++)if(!seen[i]){
  seen[i]=1;const todo=[i];for(let h=0;h<todo.length;h++)for(const g of gens){const j=index.get(key(action(parts[todo[h]],g)));if(j===undefined)throw Error("bad action");if(!seen[j]){seen[j]=1;todo.push(j);}}
  out.push({index:i,representative:parts[i],size:todo.length});
 }return out;
}
function maps(p){
 const block=Array(11),pos=Array(11),eps=Array(11).fill(0);
 p.forEach((b,k)=>b.forEach((j,h)=>{block[j]=k;pos[j]=h;if(k===0)eps[j]=h;}));
 return {block,pos,eps};
}
function matrices(q,r){
 const p=maps(BASE),qq=maps(q),rr=maps(r),w=[],off=[];
 for(let j=0;j<11;j++){
  const a=Array(6).fill(0),b=Array(9).fill(0);
  if(p.block[j]){a[p.block[j]-1]=1;b[p.block[j]-1]=-p.pos[j];}
  if(qq.block[j]){a[qq.block[j]+2]=-1;b[qq.block[j]+2]=qq.pos[j];}
  if(rr.block[j])b[rr.block[j]+5]=rr.pos[j];
  w.push(a);off.push(b);
 }
 const h=Array.from({length:11},(_,j)=>rr.eps[j]-p.eps[j]+qq.eps[j]),A=[],f=[],B=[];
 for(const block of r){const ref=block[0];for(const j of block.slice(1)){A.push(w[j].map((x,i)=>x-w[ref][i]));f.push(h[j]-h[ref]);B.push(off[j].map((x,i)=>x-off[ref][i]));}}
 return {A,f,B};
}

function assert(x,msg="assertion"){if(!x)throw Error(msg);}
function ident(n){return Array.from({length:n},(_,i)=>Array.from({length:n},(_,j)=>+(i===j)));}
function rowadd(a,i,j,q){for(let k=0;k<a[0].length;k++)a[i][k]+=q*a[j][k];}
function coladd(a,i,j,q){for(const r of a)r[i]+=q*r[j];}
function colswap(a,i,j){for(const r of a)[r[i],r[j]]=[r[j],r[i]];}
function snf(input){
 const a=input.map(r=>r.slice()),u=ident(7),v=ident(6);let k=0;
 for(;k<6;k++){
  let pi=-1,pj=-1;
  for(let i=k;i<7;i++)for(let j=k;j<6;j++)if(a[i][j]&&(pi<0||Math.abs(a[i][j])<Math.abs(a[pi][pj]))){pi=i;pj=j;}
  if(pi<0)break;[a[k],a[pi]]=[a[pi],a[k]];[u[k],u[pi]]=[u[pi],u[k]];colswap(a,k,pj);colswap(v,k,pj);
  for(let loop=0;;loop++){
   assert(loop<500,"Smith loop");let changed=false;
   for(let i=k+1;i<7;i++)if(a[i][k]){
    const q=Math.trunc(a[i][k]/a[k][k]);rowadd(a,i,k,-q);rowadd(u,i,k,-q);
    if(a[i][k]){[a[i],a[k]]=[a[k],a[i]];[u[i],u[k]]=[u[k],u[i]];}changed=true;break;
   }if(changed)continue;
   for(let j=k+1;j<6;j++)if(a[k][j]){
    const q=Math.trunc(a[k][j]/a[k][k]);coladd(a,j,k,-q);coladd(v,j,k,-q);
    if(a[k][j]){colswap(a,j,k);colswap(v,j,k);}changed=true;break;
   }if(changed)continue;
   let bi=-1;for(let i=k+1;i<7;i++)for(let j=k+1;j<6;j++)if(a[i][j]%a[k][k])bi=i;
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
  const z=a[k][k];
  for(let i=k+1;i<n;i++)for(let j=k+1;j<n;j++){
   const num=z*a[i][j]-a[i][k]*a[k][j];assert(Number.isSafeInteger(num)&&num%old===0,"det integer");a[i][j]=num/old;
  }for(let i=k+1;i<n;i++)a[i][k]=0;old=z;
 }return sign*a[n-1][n-1];
}
function mm(a,b){
 const out=Array.from({length:a.length},()=>Array(b[0].length).fill(0));
 for(let i=0;i<a.length;i++)for(let k=0;k<b.length;k++)if(a[i][k])for(let j=0;j<b[0].length;j++)out[i][j]+=a[i][k]*b[k][j];
 return out;
}
function mod(x,n){return (x%n+n)%n;}
function popcount(x){x=x-((x>>>1)&0x55555555);x=(x&0x33333333)+((x>>>2)&0x33333333);return (((x+(x>>>4))&0x0F0F0F0F)*0x01010101)>>>24;}
function zeroMask(coeff,cache){
 let c=coeff.map(x=>mod(x,3));const first=c.find(x=>x);if(first===2)c=c.map(x=>x*2%3);
 const id=c.join("");if(cache[id])return cache[id];
 const out=Array(16).fill(0);let sum=-c.reduce((a,b)=>a+b,0),prev=0;
 for(let k=0;k<512;k++){
  const gray=k^(k>>>1);
  if(k){const bit=31-Math.clz32(gray^prev);sum+=(gray>>bit&1?2:-2)*c[bit];}
  if(mod(sum,3)===0)out[gray>>>5]=(out[gray>>>5]|(1<<(gray&31)))>>>0;prev=gray;
 }
 // Independent direct dot products verify every bit of the Gray-code mask.
 for(let mask=0;mask<512;mask++){let s=0;for(let j=0;j<9;j++)s+=c[j]*(mask>>j&1?1:-1);assert(((out[mask>>>5]>>>(mask&31))&1)===+(mod(s,3)===0),"orientation mask");}
 cache[id]=out;return out;
}
function seed(q,r,s,mat,mask,branch){
 const dmax=s.diag.length?s.diag[s.diag.length-1]:1,M=6*dmax;
 const sig=Array.from({length:9},(_,i)=>mask>>i&1?1:-1);
 const rhs=mat.f.map((f,i)=>3*f+2*mat.B[i].reduce((t,x,k)=>t+x*sig[k],0));
 const ur=mm(s.U,rhs.map(x=>[x])).map(x=>x[0]);
 const theta=Array.from({length:6},(_,i)=>mod(s.diag.reduce((t,d,j)=>t+s.V[i][j]*(M/(6*d))*(ur[j]+6*branch[j]),0),M));
 const rows=[],basis=[];
 for(const [part,start] of [[BASE,0],[q,3]]){
  const m=maps(part),e=[],b=[];
  for(let j=0;j<11;j++){
   e.push(mod(m.block[j]?theta[start+m.block[j]-1]+M/3*m.pos[j]*sig[start+m.block[j]-1]:M/2*m.pos[j],M));
   b.push(m.block[j]?s.V[start+m.block[j]-1].slice(s.rank):Array(6-s.rank).fill(0));
  }rows.push(e);basis.push(b);
 }
 const w=rows[0].map((x,j)=>mod(x-rows[1][j],M));
 const wb=basis[0].map((x,j)=>x.map((y,k)=>y-basis[1][j][k]));
 for(const [e,b,p] of [[rows[0],basis[0],BASE],[rows[1],basis[1],q],[w,wb,r]])for(const block of p){
  const first=block[0],delta=block.map(j=>mod(e[j]-e[first],M)).sort((a,b)=>a-b);
  assert(delta.every((v,k)=>v===k*M/block.length),"seed cycle");
  assert(block.every(j=>b[j].every((v,k)=>v===b[first][k])),"basis cycle");
 }
 const gcd=(a,b)=>b?gcd(b,a%b):a;
 let order=1;for(const row of rows)for(const x of row){const o=M/gcd(M,mod(x-row[0],M));order=order/gcd(order,o)*o;}
 return {modulus:M,exponents:rows,basis,supports:[BASE,q,r],cross_ratio_order:order,orientation:mask,branch};
}
function runOrbit(vi,parts,reps,state){
 const q=reps[vi].representative,codes=[];
 for(let wi=0;wi<parts.length;wi++){
  const r=parts[wi],mat=matrices(q,r),s=snf(mat.A);
  assert(s.U.flat().every(x=>Math.abs(x)<=1024)&&s.V.flat().every(x=>Math.abs(x)<=1024),"matrix arithmetic bound");
  SAFE.matrices++;SAFE.maxU=Math.max(SAFE.maxU,...s.U.flat().map(Math.abs));SAFE.maxV=Math.max(SAFE.maxV,...s.V.flat().map(Math.abs));
  assert(Math.abs(det(s.U))===1&&Math.abs(det(s.V))===1,"unimodular");
  const check=mm(mm(s.U,mat.A),s.V);
  for(let i=0;i<7;i++)for(let j=0;j<6;j++)assert(check[i][j]===(i===j&&i<s.rank?s.diag[i]:0),"Smith identity");
  assert(s.diag.every((x,i)=>Number.isSafeInteger(x)&&x>0&&(!i||x%s.diag[i-1]===0)),"Smith divisibility");
  const dk=s.diag.join(",");state.smith[dk]=(state.smith[dk]||0)+1;
  const uf=mm(s.U,mat.f.map(x=>[x])).map(x=>x[0]);
  const binary=uf.slice(s.rank).every(x=>x%2===0);
  let count=0,allowed=null;
  if(binary){
   state.binary++;
   const UB=mm(s.U,mat.B);allowed=Array(16).fill(0xffffffff);
   for(let i=s.rank;i<7;i++){const m=zeroMask(UB[i],state.cache);for(let j=0;j<16;j++)allowed[j]=(allowed[j]&m[j])>>>0;}
   count=allowed.reduce((t,x)=>t+popcount(x),0);
  }
  const tor=s.diag.reduce((a,b)=>a*b,1),maxd=s.diag.length?s.diag[s.rank-1]:1,dim=6-s.rank;
  assert(tor<256&&maxd<256);
  // All entries remain exact integers; this code is arithmetic, not bitwise packing.
  const code=count?1+dim+8*maxd+2048*tor+524288*count:0;
  codes.push(code);
  if(count){
   state.feasible++;state.orientations+=count;state.components+=count*tor;
   state.byDim[dim]=(state.byDim[dim]||0)+count*tor;
   state.supportDim[dim]=(state.supportDim[dim]||0)+1;
   state.maxd[maxd]=(state.maxd[maxd]||0)+1;
   let mask=0;while(!((allowed[mask>>>5]>>>(mask&31))&1))mask++;
   if(!state.examples["d"+dim]||!state.examples["m"+maxd]){
    const branches=[];function rec(a){if(a.length===s.rank){branches.push(a.slice());return;}for(let k=0;k<s.diag[a.length];k++){a.push(k);rec(a);a.pop();}}rec([]);
    for(const b of branches){
     const ex={vi,wi,dimension:dim,...seed(q,r,s,mat,mask,b)};
     state.examples["r"+ex.cross_ratio_order]??=ex;
     state.examples["d"+dim]??=ex;
     state.examples["m"+maxd]??=ex;
    }
   }
  }
 }
 state.done=vi+1;return codes;
}

const SAFE={matrices:0,maxTerm:0,maxU:0,maxV:0};
function checked(x){assert(Number.isSafeInteger(x),"unsafe integer");SAFE.maxTerm=Math.max(SAFE.maxTerm,Math.abs(x));return x;}
function rowadd(a,i,j,q){for(let k=0;k<a[0].length;k++)a[i][k]=checked(a[i][k]+checked(q*a[j][k]));}
function coladd(a,i,j,q){for(const r of a)r[i]=checked(r[i]+checked(q*r[j]));}
function det(input){
 const a=input.map(r=>r.slice()),n=a.length;let old=1,sign=1;
 for(let k=0;k<n-1;k++){
  if(!a[k][k]){let j=k+1;while(j<n&&!a[j][k])j++;if(j===n)return 0;[a[k],a[j]]=[a[j],a[k]];sign=-sign;}
  const z=a[k][k];
  for(let i=k+1;i<n;i++)for(let j=k+1;j<n;j++){
   const num=checked(checked(z*a[i][j])-checked(a[i][k]*a[k][j]));assert(num%old===0);a[i][j]=checked(num/old);
  }for(let i=k+1;i<n;i++)a[i][k]=0;old=z;
 }return checked(sign*a[n-1][n-1]);
}
function auditSafe(parts,reps,start,end){
 for(let vi=start;vi<end;vi++)for(const r of parts){
  const s=snf(matrices(reps[vi].representative,r).A);
  assert(Math.abs(det(s.U))===1&&Math.abs(det(s.V))===1);
  SAFE.maxU=Math.max(SAFE.maxU,...s.U.flat().map(Math.abs));SAFE.maxV=Math.max(SAFE.maxV,...s.V.flat().map(Math.abs));SAFE.matrices++;
 }
 assert(7*6*2*SAFE.maxU*SAFE.maxV<Number.MAX_SAFE_INTEGER,"matrix product bound");
 return SAFE;
}

function rootAudit(parts,reps,codesByOrbit){
 const gcd=(a,b)=>b?gcd(b,a%b):a,orders={},byMax={},examples={};let cases=0,components=0;
 for(let vi=0;vi<reps.length;vi++)for(let wi=0;wi<parts.length;wi++){
  const code=codesByOrbit[vi][wi];if(!code)continue;
  const maxd=Math.floor((code-1)/8)%256;if(maxd<3)continue;
  cases++;
  const q=reps[vi].representative,r=parts[wi],mat=matrices(q,r),s=snf(mat.A),uf=mm(s.U,mat.f.map(x=>[x])).map(x=>x[0]),UB=mm(s.U,mat.B);
  const branches=[];function rec(a){if(a.length===s.rank){branches.push(a.slice());return;}for(let k=0;k<s.diag[a.length];k++){a.push(k);rec(a);a.pop();}}rec([]);
  let counted=0;
  for(let mask=0;mask<512;mask++){
   const sig=Array.from({length:9},(_,j)=>mask>>j&1?1:-1),ug=UB.map(row=>row.reduce((t,x,j)=>t+x*sig[j],0));
   if(!ug.slice(s.rank).every(x=>x%3===0))continue;
   for(const b of branches){
    let order=6;
    for(let j=0;j<s.rank;j++){
     const den=6*s.diag[j],num=3*uf[j]+2*ug[j]+6*b[j],v=den/gcd(den,mod(num,den));
     order=order/gcd(order,v)*v;
    }
    orders[order]=(orders[order]||0)+1;byMax[maxd]??={};byMax[maxd][order]=(byMax[maxd][order]||0)+1;counted++;components++;
    if(!examples[order])examples[order]={vi,wi,dimension:6-s.rank,...seed(q,r,s,mat,mask,b)};
   }
  }
  assert(counted===Math.floor((code-1)/524288)*(Math.floor((code-1)/2048)%256),"root count");
 }
 return {cases,components,orders,byMax,examples};
}

function perms3(){
 const out=[];for(let a=0;a<3;a++)for(let b=0;b<3;b++)for(let c=0;c<3;c++)if(a!==b&&a!==c&&b!==c)out.push([a,b,c]);return out;
}
function partialTransforms(){
 const pp=perms3(),out=[];
 for(const axis of pp)for(const pa of pp)for(const pb of pp){
  const a=[0,...pa.map(x=>x+1)],b=[0,...pb.map(x=>x+1)],m=[];
  for(let i=0;i<4;i++)for(let j=0;j<4;j++)for(let k=0;k<4;k++){
   const q=[];q[axis[0]]=a[i];q[axis[1]]=b[j];q[axis[2]]=k;m.push(16*q[0]+4*q[1]+q[2]);
  }out.push(m);
 }return out;
}
function canonical(t,transforms){
 let best=t;
 for(const m of transforms){
  const order=[1,2,3];
  order.sort((a,b)=>{for(let c=0;c<16;c++){const d=t.charCodeAt(m[4*c+a])-t.charCodeAt(m[4*c+b]);if(d)return d;}return 0;});
  order.unshift(0);
  let cmp=0;
  for(let j=0;j<64;j++){
   const c=t.charCodeAt(m[(j&~3)+order[j&3]])-best.charCodeAt(j);
   if(c){cmp=c;break;}
  }
  if(cmp<0){const v=Array(64);for(let j=0;j<64;j++)v[j]=t.charCodeAt(m[(j&~3)+order[j&3]]);best=String.fromCharCode(...v);}
 }
 return best;
}
function classifyBatch(raw,start,end,entries){
 const tr=partialTransforms(),classes=new Map(entries);
 for(let i=start;i<end;i++){
  const [key,r]=raw[i],c=canonical(key,tr),old=classes.get(c);
  if(old){assert(old.code===r.code,"class invariant");old.presentations+=r.presentations;old.raw_tensors++;}
  else classes.set(c,{...r,raw_tensors:1});
 }
 return Array.from(classes);
}

function fullTransforms(){
 const pp=perms3(),out=[];
 for(const axis of pp)for(const pa of pp)for(const pb of pp)for(const pc of pp){
  const maps=[[0,...pa.map(x=>x+1)],[0,...pb.map(x=>x+1)],[0,...pc.map(x=>x+1)]],m=[];
  for(let i=0;i<4;i++)for(let j=0;j<4;j++)for(let k=0;k<4;k++){
   const q=[];q[axis[0]]=maps[0][i];q[axis[1]]=maps[1][j];q[axis[2]]=maps[2][k];m.push(16*q[0]+4*q[1]+q[2]);
  }assert(new Set(m).size===64);out.push(m);
 }assert(out.length===1296);return out;
}
function verifyClasses(parts,reps,entries){
 const tr=fullTransforms(),pm=maps(BASE).block,rm=parts.map(x=>maps(x).block),qm=reps.map(x=>maps(x.representative).block);
 let checked=0,covered=0,raw=0;const seen=new Set();
 for(const [key,r] of entries){
  assert(!seen.has(key),"duplicate class");seen.add(key);
  const t=Array(64).fill(0);for(let j=0;j<11;j++)t[16*pm[j]+4*qm[r.vi][j]+rm[r.wi][j]]++;
  let equal=false;
  for(const m of tr){
   let cmp=0;for(let j=0;j<64;j++){cmp=t[m[j]]-key.charCodeAt(j);if(cmp)break;}
   assert(cmp>=0,"noncanonical support");if(cmp===0)equal=true;
  }assert(equal,"unreachable representative");checked++;covered+=r.presentations;raw+=r.raw_tensors;
 }
 return {verified:true,classes:checked,presentations:covered,raw_tensors:raw,group_size:tr.length};
}

function buildRaw(parts,reps,codesByOrbit){
 const p=maps(BASE).block,rs=parts.map(x=>maps(x).block),qs=reps.map(x=>maps(x.representative).block),raw=new Map();
 for(let vi=0;vi<reps.length;vi++)for(let wi=0;wi<parts.length;wi++){
  const code=codesByOrbit[vi][wi];if(!code)continue;
  const a=Array(64).fill(0);for(let j=0;j<11;j++)a[16*p[j]+4*qs[vi][j]+rs[wi][j]]++;
  const k=String.fromCharCode(...a),r=raw.get(k);
  if(r){assert(r.code===code,"raw invariant");r.presentations++;}
  else raw.set(k,{code,vi,wi,presentations:1});
 }return Array.from(raw);
}
function reproduceEEE(){
 const parts=partitions(),reps=orbits(parts);
 assert(parts.length===15400&&reps.length===42);
 assert(reps.reduce((s,r)=>s+r.size,0)===15400);
 const state={done:0,smith:{},binary:0,feasible:0,orientations:0,components:0,byDim:{},supportDim:{},maxd:{},examples:{},cache:{}},codes=[];
 for(let vi=0;vi<reps.length;vi++){codes.push(runOrbit(vi,parts,reps,state));console.error("Support orbit",vi+1,"of",reps.length);}
 assert(state.binary===152132&&state.feasible===147642);
 assert(state.orientations===19460244&&state.components===27688632);
 const roots=rootAudit(parts,reps,codes);
 assert(Object.keys(roots.orders).map(Number).sort((a,b)=>a-b).join(",")==="6,12,18,24,30,36,42");
 const raw=buildRaw(parts,reps,codes),classes=classifyBatch(raw,0,raw.length,[]),classCheck=verifyClasses(parts,reps,classes);
 assert(classes.length===6904&&classCheck.presentations===147642&&raw.length===64833);
 const dimensions={},maxSmith={};
 for(const [key,r] of classes){const d=(r.code-1)%8,m=Math.floor((r.code-1)/8)%256;dimensions[d]=(dimensions[d]||0)+1;maxSmith[m]=(maxSmith[m]||0)+1;}
 assert(JSON.stringify(dimensions)==='{"0":5230,"1":1544,"2":121,"3":8,"4":1}');
 const summary={verified:true,partitions:parts.length,relative_orbits:reps.length,support_presentations:parts.length*reps.length,binary_feasible:state.binary,feasible_presentations:state.feasible,feasible_orientation_choices:state.orientations,phase_components:state.components,phase_components_by_dimension:state.byDim,feasible_support_classes:classes.length,classes_by_dimension:dimensions,max_smith_by_class:maxSmith,seed_orders:[6,12,18,24,30,36,42],class_check:classCheck,integer_audit:SAFE};
 console.log(JSON.stringify(summary,null,2));
 if(typeof process!=="undefined"&&process.argv.includes("--write")){
  const fs=require("node:fs");
  fs.writeFileSync("eee_verification.json",JSON.stringify(summary,null,2)+"\n");
  fs.writeFileSync("eee_support_classes.json",JSON.stringify(classes.map(([key,r])=>({tensor:Array.from(key,x=>x.charCodeAt(0)),...r})),null,2)+"\n");
  fs.writeFileSync("eee_examples.json",JSON.stringify({four_parameters:state.examples.d4,root_examples:roots.examples},null,2)+"\n");
 }
 return summary;
}
reproduceEEE();

