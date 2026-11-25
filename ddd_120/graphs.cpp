#include <algorithm>
#include <array>
#include <cassert>
#include <fstream>
#include <iostream>
#include <map>
#include <numeric>
#include <set>
#include <vector>
using namespace std;constexpr int M=120;
#ifdef AUDIT
constexpr int P=601;
#else
constexpr int P=241;
#endif
using Row=array<int,11>;struct Family{Row g,e;int c;};
int mod(int v){v%=M;return v<0?v+M:v;}
int power(int a,int n){int b=1;for(;n;n>>=1,a=1LL*a*a%P)if(n&1)b=1LL*b*a%P;return b;}
array<int,M>pw;array<long long,2048>cache{};array<bool,2048>value{};long long stamp=0,pair_tests=0;
void tables(){for(int d=2;d*d<=P;d++)assert(P%d);int g=2;while(power(g,(P-1)/2)==1||power(g,(P-1)/3)==1||power(g,(P-1)/5)==1)g++;int z=power(g,(P-1)/M);pw[0]=1;for(int i=1;i<M;i++)pw[i]=pw[i-1]*z%P;assert(pw[M/2]!=1&&pw[M/3]!=1&&pw[M/5]!=1&&pw.back()*z%P==1);}
bool cover(const Row&x,int mask){
 if(!mask)return true;if(cache[mask]==stamp)return value[mask];cache[mask]=stamp;
#ifdef AUDIT
 int a=31-__builtin_clz((unsigned)mask);array<int,3>primes={5,3,2};
#else
 int a=__builtin_ctz((unsigned)mask);array<int,3>primes={2,3,5};
#endif
 int n=__builtin_popcount((unsigned)mask);
 for(int p:primes)if(p<=n&&(p==n||n-p>=2)){
  vector<int>choices={1<<a};
  for(int k=1;k<p;k++){int target=mod(x[a]+M/p*k);vector<int>next;for(int j=0;j<11;j++)if((mask>>j&1)&&x[j]==target)for(int m:choices)next.push_back(m|1<<j);choices=move(next);if(choices.empty())break;}
  for(int used:choices)if(cover(x,mask^used))return value[mask]=true;
 }
 return value[mask]=false;
}
bool compatible(const Row&a,const Row&b){pair_tests++;Row x;int sum=0;for(int j=0;j<11;j++){x[j]=mod(a[j]-b[j]);sum+=pw[x[j]];}if(sum%P)return false;stamp++;return cover(x,2047);}
vector<vector<bool>>adj;vector<int>best,current;
void expand(vector<int>candidates){
 if(current.size()+candidates.size()<=best.size())return;
 vector<int>order,colors,remaining=candidates;int color=0;
 while(!remaining.empty()){color++;vector<int>independent,next;for(int v:remaining){bool ok=true;for(int w:independent)if(adj[v][w]){ok=false;break;}if(ok){independent.push_back(v);order.push_back(v);colors.push_back(color);}else next.push_back(v);}remaining=move(next);}
 for(int i=int(order.size())-1;i>=0;i--){if(current.size()+colors[i]<=best.size())return;int v=order[i];current.push_back(v);vector<int>next;for(int j=0;j<i;j++)if(adj[v][order[j]])next.push_back(order[j]);if(next.empty()){if(current.size()>best.size())best=current;}else expand(next);current.pop_back();}
}
void phases(const Family&f,int g,array<int,5>&t,set<Row>&rows){
 if(g==f.c){Row x;for(int j=0;j<11;j++)x[j]=mod(f.e[j]+t[f.g[j]]);assert(x[0]==0);rows.insert(x);return;}
 for(int x=0;x<M;x++){t[g]=x;phases(f,g+1,t,rows);}
}
int main(int argc,char**argv){
 assert(argc==4);tables();ifstream in(argv[1]);map<int,vector<Family>>data;int id;Family f;
 while(in>>id>>f.c){assert(1<=f.c&&f.c<=5);for(int&v:f.g)in>>v;for(int&v:f.e)in>>v;assert(in);data[id].push_back(f);}
 string dest=argv[2],rdest=argv[3];ofstream out(dest+".writing"),rout(rdest+".writing");assert(in.eof()&&out&&rout);
 map<int,int>hist,nhist;long long total_rows=0,total_edges=0;int maxn=0;
 for(auto&[id,fs]:data){set<Row>unique;array<int,5>t{};for(auto&f:fs)phases(f,1,t,unique);vector<Row>rows(unique.begin(),unique.end());int n=rows.size();maxn=max(maxn,n);nhist[n]++;total_rows+=n;
  for(auto&r:rows){rout<<id;for(int e:r)rout<<' '<<e;rout<<'\n';}
  adj.assign(n,vector<bool>(n));vector<pair<int,int>>edges;
  for(int i=0;i<n;i++)for(int j=0;j<i;j++)if(compatible(rows[i],rows[j])){adj[i][j]=adj[j][i]=true;edges.push_back({j,i});}
  best.clear();current.clear();vector<int>candidates(n);iota(candidates.begin(),candidates.end(),0);expand(candidates);hist[best.size()]++;total_edges+=edges.size();
  out<<id<<' '<<n<<' '<<edges.size()<<' '<<best.size();for(int x:best)out<<' '<<x;for(auto[a,b]:edges)out<<' '<<a<<' '<<b;out<<'\n';
  if(best.size()>=3){cerr<<"SURVIVOR point "<<id<<" added_rows "<<best.size()<<" vertices "<<n<<endl;}
 }
 out.flush();rout.flush();assert(out&&rout);out.close();rout.close();assert(!rename((dest+".writing").c_str(),dest.c_str()));assert(!rename((rdest+".writing").c_str(),rdest.c_str()));
 cout<<"points "<<data.size()<<" rows "<<total_rows<<" max_vertices "<<maxn<<" edges "<<total_edges<<" pair_tests "<<pair_tests<<" maximum_added_rows_histogram";for(auto[k,n]:hist)cout<<' '<<k<<':'<<n;cout<<" row_count_histogram";for(auto[k,n]:nhist)cout<<' '<<k<<':'<<n;cout<<endl;
}
