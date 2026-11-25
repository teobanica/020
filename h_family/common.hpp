#pragma once
#include <algorithm>
#include <array>
#include <cassert>
#include <chrono>
#include <fstream>
#include <iostream>
#include <numeric>
#include <set>
#include <vector>
using namespace std;
constexpr int M=420,L=352800;
#ifdef AUDIT
constexpr int P=2116801;
#else
constexpr int P=1764001;
#endif

const int H[5][11]={
{0,0,0,0,0,0,0,0,0,0,0},
{0,4,8,12,16,0,10,15,5,5,15},
{0,8,16,4,12,0,15,10,15,5,5},
{0,12,4,16,8,0,5,15,10,15,5},
{0,16,12,8,4,0,5,5,15,10,15}};
int mod(int a,int m){a%=m;return a<0?a+m:a;}
struct Field{
 vector<int>pw,iv,log;int z;
 int mul(int a,int b)const{return int(1LL*a*b%P);}
 int power(int a,int n)const{int b=1;for(;n;n>>=1,a=mul(a,a))if(n&1)b=mul(a,b);return b;}
 Field():pw(L),iv(P),log(P,-1){
  for(int d=2;d*d<=P;d++)assert(P%d);
  vector<int>fac;int n=P-1;for(int d=2;d*d<=n;d++)if(n%d==0){fac.push_back(d);while(n%d==0)n/=d;}if(n>1)fac.push_back(n);
  int g=2;for(;;g++){bool ok=true;for(int d:fac)if(power(g,(P-1)/d)==1)ok=false;if(ok)break;}
  z=power(g,(P-1)/L);pw[0]=1;for(int j=1;j<L;j++)pw[j]=mul(pw[j-1],z);assert(mul(pw.back(),z)==1);
  for(int j=0;j<L;j++){assert(log[pw[j]]==-1);log[pw[j]]=j;}
  iv[1]=1;for(int j=2;j<P;j++)iv[j]=P-int(1LL*(P/j)*iv[P%j]%P);
 }
}f;
using Poly=vector<int>;
Poly trim(Poly a){while(!a.empty()&&!a.back())a.pop_back();return a;}
Poly add(Poly a,const Poly&b,int s=1){a.resize(max(a.size(),b.size()));for(int j=0;j<(int)b.size();j++)a[j]=mod(a[j]+s*b[j],P);return trim(a);}
Poly mul(const Poly&a,const Poly&b){if(a.empty()||b.empty())return {};Poly c(a.size()+b.size()-1);for(int i=0;i<(int)a.size();i++)for(int j=0;j<(int)b.size();j++)c[i+j]=(c[i+j]+f.mul(a[i],b[j]))%P;return trim(c);}
Poly rem(Poly a,const Poly&b){assert(!b.empty());
#ifdef AUDIT
 Poly original=a,q(max(0,int(a.size())-int(b.size())+1));
#endif
 int iv=f.iv[b.back()];while(a.size()>=b.size()){int d=a.size()-b.size(),v=f.mul(a.back(),iv);
#ifdef AUDIT
 q[d]=v;
#endif
 for(int j=0;j<(int)b.size();j++)a[d+j]=mod(a[d+j]-f.mul(v,b[j]),P);a=trim(a);}
#ifdef AUDIT
 assert(trim(original)==add(mul(trim(q),b),a));
#endif
 return a;}

Poly gcdp(Poly a,Poly b){while(!b.empty()){auto c=rem(a,b);a=b;b=c;}if(!a.empty()){int v=f.iv[a.back()];for(int &x:a)x=f.mul(x,v);}return a;}
Poly powerp(Poly a,int n,const Poly&g){Poly out=rem({1},g);a=rem(a,g);for(;n;n>>=1,a=rem(mul(a,a),g))if(n&1)out=rem(mul(out,a),g);return out;}
vector<int> roots(const Poly&a){
 if(a.size()==1)return {};if(a.size()==2)return {mod(-f.mul(a[0],f.iv[a[1]]),P)};
 for(int v=0;v<P;v++){
  auto d=gcdp(a,add(powerp({v,1},(P-1)/2,a),{1},-1));
  if(d.size()>1&&d.size()<a.size()){
   Poly b=a,q(a.size()-d.size()+1);while(b.size()>=d.size()){int e=b.size()-d.size(),c=f.mul(b.back(),f.iv[d.back()]);q[e]=c;for(int j=0;j<(int)d.size();j++)b[e+j]=mod(b[e+j]-f.mul(c,d[j]),P);b=trim(b);}assert(b.empty());
   auto r=roots(d),s=roots(trim(q));r.insert(r.end(),s.begin(),s.end());return r;
  }
 }assert(false);return {};
}
struct State{
 array<int,5>parent{},a{},e{};vector<int>qs;bool fixed=false;
 State(){iota(parent.begin(),parent.end(),0);}
 array<int,3>find(int b)const{
#ifdef AUDIT
 return {parent[b],a[b],e[b]};
#endif
int x=0,y=0;while(parent[b]!=b){x+=a[b];y+=e[b];b=parent[b];}return {b,x,mod(y,M)};}
 // u_b/u_c = q^d zeta_420^v.
 bool link(int b,int c,int d,int v){
  auto B=find(b),C=find(c);d-=B[1]-C[1];v=mod(v-B[2]+C[2],M);
  if(B[0]!=C[0]){
#ifdef AUDIT
   if(B[0]<C[0]){for(int i=0;i<5;i++)if(parent[i]==C[0]){parent[i]=B[0];a[i]-=d;e[i]=mod(e[i]-v,M);}}
   else {for(int i=0;i<5;i++)if(parent[i]==B[0]){parent[i]=C[0];a[i]+=d;e[i]=mod(e[i]+v,M);}}
#else
   parent[B[0]]=C[0];a[B[0]]=d;e[B[0]]=v;
#endif
   return true;
  }
  // 1=q^d*zeta^v.
  if(!d)return !v;
  assert(abs(d)<=5);
  if(!fixed){fixed=true;int dd=abs(d),vv=mod(-v*(d>0?1:-1),M);assert((L/M)%dd==0);
   for(int k=0;k<dd;k++)qs.push_back((vv+M*k)*(L/M/dd));return true;
  }
  vector<int>out;for(int q:qs)if(mod(d*q+v*(L/M),L)==0)out.push_back(q);qs=out;return !qs.empty();
 }
 int components(int k)const{int n=0;for(int j=0;j<k;j++)n+=parent[j]==j;return n;}
};
vector<int>sizes,block(11),exponent(11);vector<vector<int>>blocks;
long long tests=0,limit=100,leaves=0,nodes=0,unconnected=0,freezero=0,rootcandidates=0,rootpassed=0;
array<long long,7>bydim{};
ofstream families,points,inputs;set<vector<int>> seen;long long duplicates=0,actual=0;char profile='F';
bool normcheck(const vector<int>&x){
 for(int r=1;r<5;r++){
  int a=0,b=0,ac=0,bc=0;for(int j=0;j<11;j++){
   int e=mod(x[j]-H[r][j]*(L/20),L);if(j<5){a=(a+f.pw[e])%P;ac=(ac+f.pw[mod(-e,L)])%P;}else{b=(b+f.pw[e])%P;bc=(bc+f.pw[mod(-e,L)])%P;}
  }if(f.mul(a,ac)!=f.mul(b,bc)){
#ifdef AUDIT
   int check=0;for(int j=0;j<11;j++)for(int k=0;k<11;k++)if((j<5)==(k<5)){
    int e=mod(x[j]-x[k]-(H[r][j]-H[r][k])*(L/20),L);check=mod(check+(j<5?1:-1)*f.pw[e],P);
   }assert(check==mod(f.mul(a,ac)-f.mul(b,bc),P)&&check);
#endif
   return false;
  }
 }return true;
}

bool covered(const vector<int>&x,int mask,vector<int>&cache){
 if(!mask)return true;if(cache[mask]>=0)return cache[mask];
 int j=__builtin_ctz((unsigned)mask),rest=mask^(1<<j),n=__builtin_popcount((unsigned)mask);
 for(int prime:{2,3,5,7})if(prime<=n){
  vector<int>choices={1<<j};
  for(int h=1;h<prime;h++){
   int target=mod(x[j]+L/prime*h,L);vector<int>next;
   for(int k=0;k<11;k++)if(((rest>>k)&1)&&x[k]==target)for(int used:choices)next.push_back(used|(1<<k));
   choices=next;if(choices.empty())break;
  }
  for(int used:choices)if(covered(x,mask^used,cache))return cache[mask]=1;
 }return cache[mask]=0;
}

#ifdef AUDIT
bool point_cover(int mask,bool fixed,int phase,const vector<int>&y,int &answer);
bool point_cycle(int mask,int anchor,int prime,int label,bool fixed,int phase,const vector<int>&y,int &answer){
 if(label==prime)return point_cover(mask,fixed,phase,y,answer);
 for(int j=0;j<11;j++)if((mask>>j)&1){
  int d=int(j>=5)-int(anchor>=5),rhs=mod(L/prime*label-y[j]+y[anchor],L);
  if(!d){if(rhs)continue;if(point_cycle(mask^(1<<j),anchor,prime,label+1,fixed,phase,y,answer))return true;}
  else{int want=mod(d*rhs,L);if(fixed&&phase!=want)continue;if(point_cycle(mask^(1<<j),anchor,prime,label+1,true,want,y,answer))return true;}
 }return false;
}
bool point_cover(int mask,bool fixed,int phase,const vector<int>&y,int &answer){
 if(!mask){answer=fixed?phase:0;return true;}int n=__builtin_popcount((unsigned)mask),anchor=__builtin_ctz((unsigned)mask);
 for(int prime:{2,3,5,7})if(prime<=n&&(n==prime||n-prime>=2))if(point_cycle(mask^(1<<anchor),anchor,prime,1,fixed,phase,y,answer))return true;
 return false;
}
#endif
bool regular_at_row(const vector<int>&x,int r,int &phase){
 vector<int>y(11);for(int j=0;j<11;j++)y[j]=mod(x[j]-H[r][j]*(L/20),L);
#ifdef AUDIT
 return point_cover(2047,false,0,y,phase);
#endif

 set<int>possible={0};
 for(int j=0;j<5;j++)for(int k=5;k<11;k++)for(int prime:{2,3,5,7})for(int h=1;h<prime;h++)possible.insert(mod(y[j]-y[k]-L/prime*h,L));
 // Same-root labels cannot occur twice in a prime cycle, so h=0 is not needed for cross pairs.
 for(int q:possible){vector<int>z=y;for(int j=5;j<11;j++)z[j]=mod(z[j]+q,L);vector<int>cache(2048,-1);if(covered(z,2047,cache)){phase=q;return true;}}
 return false;
}
void checkpoint_x(const vector<int>&x,int q){
 rootcandidates++;if(!normcheck(x))return;rootpassed++;points<<tests<<' '<<q;for(int e:x)points<<' '<<e;
 array<int,4>phases{};bool yes=true;for(int r=1;r<5;r++)if(!regular_at_row(x,r,phases[r-1])){yes=false;break;}
 points<<' '<<yes;for(int e:phases)points<<' '<<e;points<<'\n';if(yes){actual++;cerr<<"ACTUAL "<<tests<<endl;}
}
void univariate(vector<int>aa,vector<int>ee,int coeff_conductor,int q_saved){
 int range=*max_element(aa.begin(),aa.end())-*min_element(aa.begin(),aa.end());assert(range<=4);Poly g;
 for(int r=1;r<5;r++){
  Poly pol(2*range+1);
#ifdef AUDIT
  array<Poly,4>factor;for(auto &v:factor)v.resize(range+1);int low=*min_element(aa.begin(),aa.end());
  for(int j=0;j<11;j++){
   int index=j<5?0:2,exp=mod(ee[j]-(L/20)*H[r][j],L),power=aa[j]-low;
   factor[index][power]=(factor[index][power]+f.pw[exp])%P;
   factor[index+1][range-power]=(factor[index+1][range-power]+f.pw[mod(-exp,L)])%P;
  }
  pol=add(mul(trim(factor[0]),trim(factor[1])),mul(trim(factor[2]),trim(factor[3])),-1);
#else
  for(int j=0;j<11;j++)for(int k=0;k<11;k++)if((j<5)==(k<5)){
   int d=aa[j]-aa[k]+range,e=mod(ee[j]-ee[k]-(L/20)*(H[r][j]-H[r][k]),L);
   pol[d]=mod(pol[d]+(j<5?1:-1)*f.pw[e],P);
  }
#endif
  g=gcdp(g,trim(pol));if(g.size()==1)return;
 }
 if(g.empty()){freezero++;families<<tests<<" 1 FREE_NORM_ZERO "<<coeff_conductor<<' '<<q_saved;for(int j=0;j<11;j++)families<<' '<<aa[j]<<' '<<ee[j];families<<'\n';return;}
 assert((coeff_conductor==420&&g.size()<=9)||(coeff_conductor==25200&&g.size()<=3));
 auto d=gcdp(g,add(powerp({0,1},L,g),{1},-1));
 for(int rt:roots(d)){
  assert(f.log[rt]>=0);int exponent=f.log[rt];vector<int>x(11);
  for(int j=0;j<11;j++)x[j]=mod(ee[j]+aa[j]*exponent,L);
  checkpoint_x(x,q_saved<0?exponent:q_saved);
 }
}
