#include <algorithm>
#include <array>
#include <cassert>
#include <fstream>
#include <iostream>
#include <numeric>
#include <set>
#include <vector>
using namespace std;constexpr int M=120;using V=array<int,4>;using Key=array<long long,67>;
long long mod(long long x,long long m=M){x%=m;return x<0?x+m:x;}
// Independent solvability test in two old parameters. The primitive integer
// kernel of at most three character equations gives every consistency condition.
bool independent(const array<V,3>& base,V a){
 array<V,3> nonzero;int n=0;
 for(const auto&r:base)if(r[0]||r[1]){assert(r[2]==0);nonzero[n++]=r;}
 assert(a[2]==0);a[3]=mod(a[3]);
 if(a[0]==0&&a[1]==0)return a[3]==0;
 nonzero[n++]=a;if(n<=1)return true;assert(n<=3);
 auto det=[](V u,V v)->long long{return (long long)u[0]*v[1]-(long long)u[1]*v[0];};
 if(n==2){auto u=nonzero[0],v=nonzero[1];if(det(u,v))return true;
  int k=u[0]?0:1;long long d=gcd(llabs(u[k]),llabs(v[k]));return mod((long long)v[k]/d*u[3]-(long long)u[k]/d*v[3])==0;}
 auto u=nonzero[0],v=nonzero[1],w=nonzero[2];
 long long A=det(v,w),B=-det(u,w),C=det(u,v);long long g=gcd(gcd(llabs(A),llabs(B)),llabs(C));
 assert(g>0);return mod((A/g)*u[3]+(B/g)*v[3]+(C/g)*w[3])==0;
}

struct Lattice{array<V,3>r{};
 bool add_impl(V a){a[3]=mod(a[3]);for(int i=0;i<3;i++){
  while(a[i]){if(!r[i][i]){if(a[i]<0)for(auto&x:a)x=-x;a[3]=mod(a[3]);r[i]=a;return true;}
   long long q=a[i]/r[i][i];for(int k=i;k<4;k++)a[k]-=q*r[i][k];a[3]=mod(a[3]);
   if(a[i]){swap(a,r[i]);if(r[i][i]<0)for(auto&x:r[i])x=-x;r[i][3]=mod(r[i][3]);}
  }
 }return mod(a[3])==0;}
 bool add(V a){
#ifdef AUDIT
 bool expected=independent(r,a);
#endif
 bool answer=add_impl(a);
#ifdef AUDIT
 assert(answer==expected);for(auto &v:r)for(auto x:v)assert(llabs(x)<1000000000);
#endif
 return answer;}
 void canonical(){for(int j=0;j<3;j++)if(r[j][j])for(int i=0;i<j;i++)if(r[i][i]){
  long long q=(r[i][j]-mod(r[i][j],r[j][j]))/r[j][j];for(int k=j;k<4;k++)r[i][k]-=q*r[j][k];r[i][3]=mod(r[i][3]);}}
 V reduce(V a)const{for(int i=0;i<3;i++)if(r[i][i]){long long q=(a[i]-mod(a[i],r[i][i]))/r[i][i];for(int k=i;k<3;k++)a[k]-=q*r[i][k];a[3]+=q*r[i][3];}a[3]=mod(a[3]);return a;}
};
struct State{array<int,5>p;array<V,5>e{};Lattice l;State(){iota(p.begin(),p.end(),0);}
 pair<int,V>find(int a)const{
#ifdef AUDIT
 return {p[a],e[a]};
#else
 V x{};while(p[a]!=a){for(int k=0;k<4;k++)x[k]+=e[a][k];a=p[a];}x[3]=mod(x[3]);return {a,x};
#endif
 }
 bool link(int a,int b,V x){auto [ra,ea]=find(a);auto [rb,eb]=find(b);for(int k=0;k<4;k++)x[k]+=-ea[k]+eb[k];x[3]=mod(x[3]);if(ra==rb){x[3]=mod(-x[3]);return l.add(x);}
#ifdef AUDIT
 if(ra<rb){for(int j=0;j<5;j++)if(p[j]==rb){p[j]=ra;for(int k=0;k<4;k++)e[j][k]-=x[k];e[j][3]=mod(e[j][3]);}}
 else{for(int j=0;j<5;j++)if(p[j]==ra){p[j]=rb;for(int k=0;k<4;k++)e[j][k]+=x[k];e[j][3]=mod(e[j][3]);}}
#else
 p[ra]=rb;e[ra]=x;
#endif
 return true;}

};
struct Family{array<int,11>b;array<V,11>x;Lattice l;};
array<array<V,11>,2>old;
using PairKey=array<int,35>;
set<PairKey>pairs;set<Key>finals;long long tests=0,nodes[2]{},leaves[2]{},limit=1000000000;vector<int>sizes;vector<vector<int>>blocks;Family initial;int chosen_profile;ofstream output;
Family normalized(const Family&f,const State&s){Family o;o.l=s.l;o.l.canonical();array<int,5>grp;grp.fill(-1);array<V,5>anchor;int nc=0;
 for(int j=0;j<11;j++){auto [g,v]=s.find(f.b[j]);for(int k=0;k<4;k++)v[k]+=f.x[j][k];if(grp[g]<0){grp[g]=nc++;anchor[g]=v;}o.b[j]=grp[g];for(int k=0;k<4;k++)v[k]-=anchor[g][k];o.x[j]=o.l.reduce(v);}return o;}
Key key(const Family&f){Key a;int p=0;for(auto&r:f.l.r)for(auto x:r)a[p++]=x;for(int j=0;j<11;j++){a[p++]=f.b[j];for(auto x:f.x[j])a[p++]=x;}assert(p==67);return a;}
PairKey pair_key(const Family&f){
 PairKey k;assert(f.l.r[0][1]==0&&f.l.r[0][2]==0);for(int i=1;i<3;i++)for(auto x:f.l.r[i])assert(x==0);
 k[0]=f.l.r[0][0];k[1]=f.l.r[0][3];int p=2;
 for(int j=0;j<11;j++){assert(f.x[j][1]==0&&f.x[j][2]==0);k[p++]=f.b[j];k[p++]=f.x[j][0];k[p++]=f.x[j][3];}assert(p==35);return k;
}
void cover(int mask,const State&s,const Family&f,int stage);
void leaf(const State&s,const Family&f,int stage){leaves[stage]++;auto g=normalized(f,s);auto k=key(g);
 if(stage==0){if(pairs.insert(pair_key(g)).second){State t;t.l=g.l;assert(t.l.add({0,120,0,0}));cover(2047,t,g,1);}}
 else if(finals.insert(k).second){for(auto x:k)output<<x<<' ';output<<'\n';}
}
void cycle(int mask,int anchor,int prime,int label,const State&s,const Family&f,int stage){nodes[stage]++;if(label==prime){cover(mask,s,f,stage);return;}
 for(int j=0;j<11;j++)if(mask>>j&1){V v;for(int k=0;k<4;k++)v[k]=f.x[anchor][k]-f.x[j][k]+old[stage][j][k]-old[stage][anchor][k];v[3]+=M/prime*label;State t=s;if(t.link(f.b[j],f.b[anchor],v))cycle(mask^(1<<j),anchor,prime,label+1,t,f,stage);}}
void cover(int mask,const State&s,const Family&f,int stage){if(!mask){leaf(s,f,stage);return;}int n=__builtin_popcount((unsigned)mask);
#ifdef AUDIT
 int j=__builtin_ctz((unsigned)mask);vector<int>primes={5,3,2};
#else
 int j=__builtin_ctz((unsigned)mask);vector<int>primes={2,3,5};
#endif
 for(int p:primes)if(p<=n&&(p==n||n-p>=2))cycle(mask^(1<<j),j,p,1,s,f,stage);}
void test(){if(tests>=limit)return;tests++;State s;assert(s.l.add({120,0,0,0}));cover(2047,s,initial,0);if(tests%10000==0)cerr<<tests<<" pair "<<pairs.size()<<" final "<<finals.size()<<" nodes "<<nodes[0]<<' '<<nodes[1]<<endl;}
void labels(int k){if(tests>=limit)return;if(k==(int)blocks.size()){test();return;}auto&B=blocks[k];vector<int>p(B.size()-1);iota(p.begin(),p.end(),1);do{
#ifdef REDUCE_SEVENTH
if(B.size()==7 && p[0]!=1)continue;
#endif
initial.x[B[0]][3]=0;for(int j=1;j<(int)B.size();j++)initial.x[B[j]][3]=p[j-1]*M/B.size();labels(k+1);}while(next_permutation(p.begin(),p.end()));}
void partition(int k,vector<int>rem){if(tests>=limit)return;if(k==(int)sizes.size()){labels(0);return;}for(int mask=0;mask<(1<<(int)rem.size());mask++)if(__builtin_popcount((unsigned)mask)==sizes[k]){vector<int>b,r;for(int j=0;j<(int)rem.size();j++)((mask>>j&1)?b:r).push_back(rem[j]);if(k&&sizes[k]==sizes[k-1]&&b[0]<blocks.back()[0])continue;for(int j:b)initial.b[j]=k;blocks.push_back(b);partition(k+1,r);blocks.pop_back();}}
#ifndef NO_MAIN
int main(int argc,char**argv){assert(argc>=4);ifstream in(argv[1]);for(auto&r:old)for(auto&x:r){in>>x[3];x[3]%=M;for(int k=0;k<3;k++)in>>x[k];}assert(in);chosen_profile=argv[2][0]-'B';assert(1<=chosen_profile&&chosen_profile<5);if(argc>4)limit=stoll(argv[4]);string final_path=argv[3];output.open(final_path+".writing");assert(output);sizes=vector<vector<int>>{{7,2,2},{5,3,3},{5,2,2,2},{3,3,3,2},{3,2,2,2,2}}[chosen_profile];vector<int>rem(11);iota(rem.begin(),rem.end(),0);partition(0,rem);output.flush();assert(output);output.close(); assert(!rename((final_path+".writing").c_str(),final_path.c_str()));cout<<"tests "<<tests<<" pairs "<<pairs.size()<<" final "<<finals.size()<<" nodes "<<nodes[0]<<' '<<nodes[1]<<" leaves "<<leaves[0]<<' '<<leaves[1]<<endl;}
#endif
