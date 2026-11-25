#define NO_MAIN
#include "../full.cpp"
#include <map>
#include <sstream>
namespace triple_check {
using Relation=array<int,12>;
struct RowFamily{array<int,11>g,e;int c;};
array<V,11>quotient;
set<Relation>relations;
long long visited=0;
Relation relkey(Lattice h){h.canonical();Relation k;int p=0;for(auto&r:h.r)for(int x:r)k[p++]=x;return k;}
void cover(int mask,const Lattice&h);
void cycle(int mask,int anchor,int p,int label,const Lattice&h){
 visited++;if(label==p){cover(mask,h);return;}
 for(int j=0;j<11;j++)if(mask>>j&1){
  V equation;for(int k=0;k<3;k++)equation[k]=quotient[j][k]-quotient[anchor][k];
  equation[3]=mod(M/p*label-quotient[j][3]+quotient[anchor][3]);
  Lattice next=h;if(next.add(equation))cycle(mask^(1<<j),anchor,p,label+1,next);
 }
}
void cover(int mask,const Lattice&h){
 if(!mask){relations.insert(relkey(h));return;}
 int count=__builtin_popcount((unsigned)mask);
#ifdef AUDIT
 int anchor=31-__builtin_clz((unsigned)mask);vector<int>primes={7,5,3,2};
#else
 int anchor=__builtin_ctz((unsigned)mask);vector<int>primes={2,3,5,7};
#endif
 for(int p:primes)if(p<=count&&(count==p||count-p>=2))cycle(mask^(1<<anchor),anchor,p,1,h);
}
set<Relation>pair_loci(const RowFamily&a,const RowFamily&b){
 for(int j=0;j<11;j++)quotient[j]={a.g[j],-b.g[j],0,(int)mod(a.e[j]-b.e[j])};
 relations.clear();Lattice h;cover(2047,h);return relations;
}
bool insert_relation(Lattice&h,const Relation&r,int a,int b){
 for(int i=0;i<3;i++){assert(r[4*i+2]==0);V v{};v[a]=r[4*i];v[b]=r[4*i+1];v[3]=r[4*i+3];if(!h.add(v))return false;}return true;
}
}
int main(int argc,char**argv){
 using namespace triple_check;assert(argc==4);ifstream in(argv[1]);map<int,vector<RowFamily>>data;
 int id;RowFamily row;while(in>>id>>row.c){for(int&v:row.g)in>>v;for(int&v:row.e)in>>v;assert(in&&row.c<=2);data[id].push_back(row);}
 ifstream graph(argv[2]);string final=argv[3];ofstream out(final+".writing");assert(in.eof()&&graph&&out);string l;
 long long total_triangles=0,total_pairs=0,total_combinations=0,possible=0,point_count=0;
 while(getline(graph,l)){
  istringstream s(l);int id,n,m,k;s>>id;string type;s>>type>>n>>m>>k;assert(s&&(type=="FINITE"||type=="PARAMETRIC"));int ignored;for(int j=0;j<k;j++)s>>ignored;
  vector<vector<bool>>adj(n,vector<bool>(n));for(int j=0;j<m;j++){int a,b;s>>a>>b;assert(s);adj[a][b]=adj[b][a]=true;}
  vector<array<int,3>>triangles;set<pair<int,int>>edges;
  for(int a=0;a<n;a++)for(int b=a+1;b<n;b++)if(adj[a][b])for(int c=b+1;c<n;c++)if(adj[a][c]&&adj[b][c]){
   triangles.push_back({a,b,c});edges.insert({a,b});edges.insert({a,c});edges.insert({b,c});
  }
  if(triangles.empty())continue;point_count++;assert((int)data[id].size()==n);map<pair<int,int>,set<Relation>>loci;
  for(auto [a,b]:edges){auto found=pair_loci(data[id][a],data[id][b]);assert(!found.empty());loci[{a,b}]=found;total_pairs++;
   for(auto&r:found){out<<"PAIR "<<id<<' '<<a<<' '<<b;for(int x:r)out<<' '<<x;out<<'\n';}
  }
  for(auto [a,b,c]:triangles){long long combinations=0;set<Relation>consistent;
   for(auto&ab:loci[{a,b}])for(auto&ac:loci[{a,c}])for(auto&bc:loci[{b,c}]){
    combinations++;Lattice h;if(insert_relation(h,ab,0,1)&&insert_relation(h,ac,0,2)&&insert_relation(h,bc,1,2))consistent.insert(relkey(h));
   }
   out<<"TRIANGLE "<<id<<' '<<a<<' '<<b<<' '<<c<<' '<<loci[{a,b}].size()<<' '<<loci[{a,c}].size()<<' '<<loci[{b,c}].size()<<' '<<combinations<<' '<<consistent.size()<<'\n';
   for(auto&r:consistent){out<<"SURVIVOR "<<id<<' '<<a<<' '<<b<<' '<<c;for(int x:r)out<<' '<<x;out<<'\n';}
   total_triangles++;total_combinations+=combinations;if(!consistent.empty())possible++;
  }
 }
 out.flush();assert(out);out.close();assert(!rename((final+".writing").c_str(),final.c_str()));
 cout<<"points "<<point_count<<" triangles "<<total_triangles<<" pairs "<<total_pairs<<" combinations "<<total_combinations<<" feasible_triangles "<<possible<<" nodes "<<visited<<endl;
}
