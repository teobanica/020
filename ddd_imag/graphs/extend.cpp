#define main old_pair_main
#include "pair.cpp"
#undef main
#include <cstdint>
using Row=array<int,11>;
int stamp=0;array<int,2048>cache_stamp{};array<bool,2048>cache_value{};
bool point_cover(const Row&y,int mask){
 if(!mask)return true;if(cache_stamp[mask]==stamp)return cache_value[mask];cache_stamp[mask]=stamp;
 int anchor=__builtin_ctz((unsigned)mask);
#ifdef AUDIT
 anchor=31-__builtin_clz((unsigned)mask);
#endif
 int n=__builtin_popcount((unsigned)mask),rest=mask^(1<<anchor);
 for(int p:{2,3,5,7})if(p<=n&&(n==p||n-p>=2)){
  vector<int>choices={1<<anchor};for(int h=1;h<p;h++){int target=mod(y[anchor]+M/p*h);vector<int>next;for(int j=0;j<11;j++)if(((rest>>j)&1)&&y[j]==target)for(int c:choices)next.push_back(c|(1<<j));choices=move(next);if(choices.empty())break;}
  for(int used:choices)if(point_cover(y,mask^used))return cache_value[mask]=true;
 }return cache_value[mask]=false;
}
namespace stage2 {
set<Key>result;long long branches=0,complete=0;
void leaf(const State&s){complete++;Key key;array<int,5>groups,anchor;groups.fill(-1);int c=0;
 for(int j=0;j<11;j++){auto v=s.find(block[j]);if(groups[v.first]<0){groups[v.first]=c++;anchor[v.first]=mod(v.second+exponent[j]);}key[j]=groups[v.first];key[11+j]=mod(v.second+exponent[j]-anchor[v.first]);}result.insert(key);
}
void cover(int mask,const State&s);
void cycle(int mask,int anchor,int prime,int label,const State&s){branches++;if(label==prime){stage2::cover(mask,s);return;}
 for(int j=0;j<11;j++)if((mask>>j)&1){State t=s;int v=mod(label*(M/prime)+exponent[anchor]-exponent[j]+u[j]-u[anchor]);if(t.link(block[j],block[anchor],v))stage2::cycle(mask^(1<<j),anchor,prime,label+1,t);}
}
void cover(int mask,const State&s){if(!mask){stage2::leaf(s);return;}int n=__builtin_popcount((unsigned)mask),j=__builtin_ctz((unsigned)mask);mask^=1<<j;for(int p:{2,3,5,7})if(p<=n&&(n==p||n-p>=2))stage2::cycle(mask,j,p,1,s);}
}
#ifdef AUDIT
constexpr int P=3361;
#else
constexpr int P=2521;
#endif
array<array<int,40>,840>positive{},negative{};
void tables(){
 for(int d=2;d*d<=P;d++)assert(P%d);
 auto pow=[](int a,int n){int b=1;for(;n;n>>=1,a=1LL*a*a%P)if(n&1)b=1LL*b*a%P;return b;};
 int g=2;while(pow(g,(P-1)/2)==1||pow(g,(P-1)/3)==1||pow(g,(P-1)/5)==1||pow(g,(P-1)/7)==1)g++;
 int z=pow(g,(P-1)/840);array<int,840>pw;pw[0]=1;for(int j=1;j<840;j++)pw[j]=1LL*pw[j-1]*z%P;assert(1LL*pw.back()*z%P==1);
 for(int x=0;x<840;x++)for(int h=0;h<40;h++){positive[x][h]=pw[mod(x-21*h)];negative[x][h]=pw[mod(-x+21*h)];}
}
bool necessary(int c,const Key&key,const Row&v){
 if(c==1){int a=0;for(int j=0;j<11;j++)a+=positive[key[11+j]][v[j]];return a%P==0;}
 if(c==2){array<int,2>a{},b{};for(int j=0;j<11;j++){a[key[j]]+=positive[key[11+j]][v[j]];b[key[j]]+=negative[key[11+j]][v[j]];}return (1LL*a[0]*b[0]-1LL*a[1]*b[1])%P==0;}
 return true;
}
#ifndef EXTEND_NO_MAIN
int main(int argc,char**argv){tables();assert(argc==5);int group=stoi(argv[1]);ifstream in(argv[2]);int n;in>>n;vector<pair<int,Row>>oldrows;for(int r=0;r<n;r++){int id,nv;in>>id>>nv;Row a;for(int&x:a)in>>x;for(int k=0;k<nv;k++){int idx;in>>idx;for(int&x:a)in>>x;if(id==group)oldrows.push_back({idx,a});}}assert(!oldrows.empty());
 ifstream ff(argv[3]);vector<pair<int,Key>>families;char profile;long long index;int c;while(ff>>profile>>index>>c){Key key;for(int&x:key)ff>>x;assert(ff);families.push_back({c,key});}
 string temp=string(argv[4])+".writing";ofstream out(temp);assert(out);long long rows_with_extensions=0,extensions=0,pointtests=0,familytests=0;array<long long,6>dims{};
 for(auto&[id,old]:oldrows){stage2::result.clear();for(int j=0;j<11;j++)u[j]=21*old[j];
  for(auto&[c,key]:families){
   if(!necessary(c,key,old))continue;
   if(c==1){pointtests++;Row y;for(int j=0;j<11;j++)y[j]=mod(key[11+j]-u[j]);stamp++;if(point_cover(y,2047))stage2::result.insert(key);}
   else{familytests++;for(int j=0;j<11;j++){block[j]=key[j];exponent[j]=key[11+j];}State s;stage2::cover(2047,s);}
  }
  if(!stage2::result.empty())rows_with_extensions++;for(auto&key:stage2::result){int c=*max_element(key.begin(),key.begin()+11)+1;dims[c-1]++;out<<id<<' '<<c;for(int x:key)out<<' '<<x;out<<'\n';}extensions+=stage2::result.size();cout<<"seed "<<id<<" families "<<stage2::result.size()<<endl;
 }
 out.flush();assert(out);out.close();assert(!out.fail());assert(rename(temp.c_str(),argv[4])==0);
 cout<<"summary group "<<group<<" seeds "<<oldrows.size()<<" extendible "<<rows_with_extensions<<" extension_families "<<extensions<<" pointtests "<<pointtests<<" familytests "<<familytests<<" branches "<<stage2::branches<<" covers "<<stage2::complete<<" dimensions";for(auto x:dims)cout<<' '<<x;cout<<endl;
 return 0;
}
#endif
