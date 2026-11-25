#include <algorithm>
#include <array>
#include <cassert>
#include <fstream>
#include <iostream>
#include <map>
#include <set>
#include <vector>
using namespace std;constexpr int M=840;using Input=array<int,67>;using Key=array<int,22>;
int mod(int x){x%=M;return x<0?x+M:x;}
Key point_row(const Input&r,int a,int b){Key out;array<int,5>labels,anchors;labels.fill(-1);int count=0;
 for(int j=0;j<11;j++){int g=r[12+5*j],e=mod(r[16+5*j]+a*r[13+5*j]+b*r[14+5*j]);assert(r[15+5*j]==0);if(labels[g]<0){labels[g]=count++;anchors[g]=e;}out[j]=labels[g];out[11+j]=mod(e-anchors[g]);}return out;}
int main(int argc,char**argv){assert(argc==4);ifstream in(argv[1]);vector<Input>records;Input r;
 while(in>>r[0]){for(int i=1;i<67;i++)in>>r[i];assert(in);records.push_back(r);}
 ifstream pf(argv[2]);map<int,set<Key>>given;int id,c;Key k;while(pf>>id>>c){for(int&x:k)pf>>x;assert(pf);assert(c==1+*max_element(k.begin(),k.begin()+11));given[id].insert(k);}assert(given.size()==3600);
 long long total=0;
 for(int a=0;a<60;a++)for(int b=0;b<60;b++){set<Key>expected;
  for(const auto&r:records){bool ok=true;for(int i=0;i<12;i+=4)if(mod(14*a*r[i]+14*b*r[i+1]-r[i+3]))ok=false;if(ok)expected.insert(point_row(r,14*a,14*b));}
  assert(expected==given[60*a+b]);total+=expected.size();
 }
 ifstream cf(argv[3]);int nc;cf>>nc;assert(nc==180);long long curve_rows=0;
 for(int cc=0;cc<nc;cc++){int id,A,B,q,n,D;cf>>id>>A>>B>>q>>n>>D;set<Key>stored;for(int i=0;i<n;i++){Key x;for(int&v:x)cf>>v;stored.insert(x);}assert(cf&&stored.size()==(unsigned)n);
  int a0=0,b0=0,da=0,db=0;if(A==1&&B==0){a0=14*q;db=1;}else if(A==0&&B==1){b0=14*q;da=1;}else{assert(A==1&&B==-1);a0=14*q;da=db=1;}
  set<Key>expected;
  for(const auto&r:records){bool fits=true;for(int i=0;i<12;i+=4)if(r[i]*da+r[i+1]*db||mod(r[i]*a0+r[i+1]*b0-r[i+3]))fits=false;if(!fits)continue;Key x;
   for(int j=0;j<11;j++){assert(r[12+5*j]==0&&r[15+5*j]==0);x[j]=r[13+5*j]*da+r[14+5*j]*db;x[11+j]=mod(r[13+5*j]*a0+r[14+5*j]*b0+r[16+5*j]);assert(x[11+j]%14==0);}expected.insert(x);
  }
  assert(expected==stored);
  for(auto&x:stored)for(auto&y:stored){int lo=100,hi=-100;for(int j=0;j<11;j++){lo=min(lo,x[j]-y[j]);hi=max(hi,x[j]-y[j]);}assert(hi-lo<=2);}
  curve_rows+=stored.size();
 }
 cout<<"records "<<records.size()<<" point_graphs 3600 point_families "<<total<<" curves "<<nc<<" curve_families "<<curve_rows<<" exact_expansion_verified"<<endl;
}
