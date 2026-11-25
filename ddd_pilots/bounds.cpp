// Orthogonality upper bounds on affine parameter strata.
#include <algorithm>
#include <array>
#include <cassert>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <map>
#include <numeric>
#include <set>
#include <vector>
using namespace std;
constexpr int L=3360,N=11;
#ifdef AUDIT
constexpr int P=20161;
#else
constexpr int P=3361;
#endif
using V=array<int,4>; // row phase group, two powers, fixed root
using Row=array<V,N>;
int mod(long long x){x%=L;return x<0?x+L:x;}
array<int,L> roots;
void table(){
 for(int k=2;k*k<=P;k++)assert(P%k);
 for(int g=2;g<P;g++){
  roots[0]=1;for(int j=1;j<L;j++)roots[j]=(long long)roots[j-1]*g%P;
  if((long long)roots.back()*g%P==1&&roots[L/2]!=1&&roots[L/3]!=1&&roots[L/5]!=1&&roots[L/7]!=1)return;
 }
 assert(false);
}
struct Raw{int h;Row v;};
Row norm(Row row,int mask){
 array<V,2>anchor{};array<bool,2>seen{};
 for(int j=0;j<N;j++){
  int g=mask>>j&1;row[j][0]=g;
  if(!seen[g]){seen[g]=true;anchor[g]=row[j];}
  for(int k=1;k<4;k++)row[j][k]-=anchor[g][k];
  row[j][3]=mod(row[j][3]);
 }
 return row;
}
vector<Row> reduce(const set<Row>&all){
 map<int,set<Row>> circles;vector<Row>keep;
 for(const Row&row:all){
  int mask=0;for(int j=0;j<N;j++)mask|=row[j][0]<<j;
  if(mask){circles[mask].insert(norm(row,mask));keep.push_back(row);}
 }
 for(const Row&row:all){
  if(any_of(row.begin(),row.end(),[](auto v){return v[0]!=0;}))continue;
  bool covered=false;for(auto&[mask,cs]:circles)if(cs.count(norm(row,mask))){covered=true;break;}
  if(!covered)keep.push_back(row);
 }
 return keep;
}
bool edge(const Row&x,const Row&y){
 array<array<int,4>,2>cells{};int n=0;array<long long,2>pos{},neg{};
 for(int j=0;j<N;j++){
  array<int,4>c{x[j][1]-y[j][1],x[j][2]-y[j][2],x[j][0],y[j][0]};
  int k=0;while(k<n&&cells[k]!=c)k++;
  if(k==n){if(n==2)return true;cells[n++]=c;}
  int e=mod(x[j][3]-y[j][3]);pos[k]+=roots[e];neg[k]+=roots[mod(-e)];
 }
 if(n==1)return pos[0]%P==0&&neg[0]%P==0;
 return (pos[0]*neg[0]-pos[1]*neg[1])%P==0;
}
struct Graph{
 vector<vector<char>>adj;int best=0;
 // Independent sets in each color class give a rigorous bound in every prefix.
 void expand(const vector<int>&v,int k){
  vector<int>rest=v,order,color;int n=0;
  while(!rest.empty()){n++;vector<int>used,next;
   for(int i:rest){
    if(any_of(used.begin(),used.end(),[&](int j){return adj[i][j];}))next.push_back(i);
    else{used.push_back(i);order.push_back(i);color.push_back(n);}
   }rest=move(next);
  }
  for(int i=int(order.size())-1;i>=0;i--){
   if(k+color[i]<=best||best>=8)return;
   int x=order[i];best=max(best,k+1);vector<int>w;
   for(int j=0;j<i;j++)if(adj[x][order[j]])w.push_back(order[j]);
   if(!w.empty())expand(w,k+1);
  }
 }
#ifdef AUDIT
 // Ordered common-neighbour enumeration, without coloring.
 void enumerate(const vector<int>&v,int k){
  best=max(best,k);
  for(int i=0;i<(int)v.size()&&best<8;i++){
   if(k+(int)v.size()-i<=best)return;
   vector<int>next;for(int j=i+1;j<(int)v.size();j++)if(adj[v[i]][v[j]])next.push_back(v[j]);
   enumerate(next,k+1);
  }
 }
#endif
};
int main(int argc,char**argv){
 assert(argc==4);table();ifstream in(argv[1]),domains(argv[2]);ofstream out(argv[3]);int h;
 vector<Raw> data;
 while(in>>h){
  int dummy;for(int j=0;j<12;j++)in>>dummy;Raw r;r.h=h;
  for(V&v:r.v){int third;in>>v[0]>>v[1]>>v[2]>>third>>v[3];assert(v[0]<=1&&!third);}
  assert(in);data.push_back(r);
 }
 int d,n;array<int,2>b,c,f;
 int processed=0;
 while(domains>>h>>d>>b[0]>>b[1]>>c[0]>>c[1]>>f[0]>>f[1]>>n){
  ++processed;set<Row>unique;for(int t=0;t<n;t++){
   int i;domains>>i;assert(data.at(i).h==h);Row r=data[i].v;
   for(V&v:r){int a=v[1],s=v[2];v[1]=a*c[0]+s*c[1];v[2]=a*f[0]+s*f[1];v[3]=mod(4*v[3]+a*b[0]+s*b[1]);}
   unique.insert(r);
  }
  auto v=reduce(unique);Graph g;int nv=v.size(),edges=0;g.adj.assign(nv,vector<char>(nv));
  for(int i=0;i<nv;i++)for(int j=0;j<i;j++)if(edge(v[i],v[j])){g.adj[i][j]=g.adj[j][i]=1;edges++;}
  vector<int>order(nv);iota(order.begin(),order.end(),0);
#ifdef AUDIT
  g.enumerate(order,0);
#else
  g.expand(order,0);
#endif
  out<<h<<' '<<d<<' '<<b[0]<<' '<<b[1]<<' '<<c[0]<<' '<<c[1]<<' '<<f[0]<<' '<<f[1]<<' '<<n<<' '<<nv<<' '<<edges<<' '<<g.best<<'\n';
 }
 out.flush();assert(out);out.close();assert(out);cerr<<"graphs "<<processed<<"\n";
}
