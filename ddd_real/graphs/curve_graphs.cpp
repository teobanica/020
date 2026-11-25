#include <algorithm>
#include <array>
#include <cassert>
#include <fstream>
#include <iostream>
#include <vector>
#include <map>
#include <set>
using namespace std;constexpr int M=840;
#ifdef AUDIT
constexpr int P=3361;
#else
constexpr int P=2521;
#endif
using Row=array<int,11>;struct Family{Row k,e;};
int mod(int v){v%=M;return v<0?v+M:v;}
array<int,M>pw;array<int,2048>cache{};array<bool,2048>value{};int stamp=0;long long tests=0;
int power(int a,int n){int b=1;for(;n;n>>=1,a=1LL*a*a%P)if(n&1)b=1LL*b*a%P;return b;}
void tables(){for(int d=2;d*d<=P;d++)assert(P%d);int g=2;while(power(g,(P-1)/2)==1||power(g,(P-1)/3)==1||power(g,(P-1)/5)==1||power(g,(P-1)/7)==1)g++;int z=power(g,(P-1)/M);pw[0]=1;for(int i=1;i<M;i++)pw[i]=1LL*pw[i-1]*z%P;assert(1LL*pw.back()*z%P==1);}
bool cover(const Row&x,int mask){if(!mask)return true;if(cache[mask]==stamp)return value[mask];cache[mask]=stamp;
#ifdef AUDIT
 int a=31-__builtin_clz((unsigned)mask);array<int,4>primes={7,5,3,2};
#else
 int a=__builtin_ctz((unsigned)mask);array<int,4>primes={2,3,5,7};
#endif
 int n=__builtin_popcount((unsigned)mask);for(int p:primes)if(p<=n&&(n==p||n-p>=2)){
  vector<int>choices={1<<a};for(int k=1;k<p;k++){int target=mod(x[a]+M/p*k);vector<int>next;for(int j=0;j<11;j++)if((mask>>j&1)&&x[j]==target)for(int m:choices)next.push_back(m|1<<j);choices=move(next);if(choices.empty())break;}
  for(int used:choices)if(cover(x,mask^used))return value[mask]=true;
 }return value[mask]=false;}
bool possible(const Family&a,const Family&b,int phase){tests++;Row x,k;for(int j=0;j<11;j++){k[j]=a.k[j]-b.k[j];x[j]=mod(a.e[j]-b.e[j]+(phase<0?0:phase*k[j]));}
 if(phase>=0){int s=0;for(int z:x)s+=pw[z];if(s%P)return false;stamp++;return cover(x,2047);}
 map<int,int>groups;for(int j=0;j<11;j++)groups[k[j]]|=1<<j;stamp++;
 for(auto &[power,mask]:groups){int s=0;for(int j=0;j<11;j++)if(mask>>j&1)s+=pw[x[j]];if(s%P||!cover(x,mask))return false;}return true;}
int main(int argc,char**argv){assert(argc==3);tables();ifstream in(argv[1]);string dest=argv[2];ofstream out(dest+".writing");int count;in>>count;long long graphs=0,triangles=0;map<int,int>edge_hist;int maxn=0;
 for(int cc=0;cc<count;cc++){int id,A,B,q,n,D;in>>id>>A>>B>>q>>n>>D;assert(in&&D<=2);maxn=max(maxn,n);vector<Family>rows(n);for(auto&r:rows){for(int&v:r.k)in>>v;for(int&v:r.e){in>>v;assert(v%14==0);}}assert(in);
  for(int t=-1;t<M;t++){if(t>=0&&t%14==0)continue;vector<vector<bool>>adj(n,vector<bool>(n));vector<pair<int,int>>edges;array<int,3>tri={-1,-1,-1};
   for(int i=0;i<n;i++)for(int j=0;j<i;j++)if(possible(rows[i],rows[j],t)){adj[i][j]=adj[j][i]=true;edges.push_back({j,i});if(tri[0]<0)for(int h=0;h<j;h++)if(adj[i][h]&&adj[j][h]){tri={h,j,i};break;}}
   out<<id<<' '<<t<<' '<<n<<' '<<edges.size()<<' '<<tri[0]<<' '<<tri[1]<<' '<<tri[2];for(auto [i,j]:edges)out<<' '<<i<<' '<<j;out<<'\n';graphs++;edge_hist[edges.size()]++;if(tri[0]>=0){triangles++;cerr<<"TRIANGLE "<<id<<' '<<t<<' '<<tri[0]<<' '<<tri[1]<<' '<<tri[2]<<endl;}
  }
  if(cc%10==0)cerr<<"curves "<<cc+1<<" graphs "<<graphs<<" triangle_graphs "<<triangles<<endl;
 }
 out.flush();assert(out);out.close();assert(!rename((dest+".writing").c_str(),dest.c_str()));cout<<"curves "<<count<<" graphs "<<graphs<<" pair_tests "<<tests<<" max_rows "<<maxn<<" triangle_graphs "<<triangles<<" edge_hist";for(auto[k,v]:edge_hist)cout<<' '<<k<<':'<<v;cout<<endl;return 0;
}
