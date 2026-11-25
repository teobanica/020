#include <algorithm>
#include <array>
#include <cassert>
#include <fstream>
#include <iostream>
#include <map>
#include <set>
#include <vector>
using namespace std;using Record=array<int,67>;using Key=array<int,22>;using Row=array<int,11>;constexpr int M=120;
int mod(int x){x%=M;return x<0?x+M:x;}
Key substitute(const Record&r,int a,int b){Key out;array<int,5>labels,anchors;labels.fill(-1);int count=0;
 for(int j=0;j<11;j++){int g=r[12+5*j],e=mod(r[16+5*j]+a*r[13+5*j]+b*r[14+5*j]);assert(r[15+5*j]==0);if(labels[g]<0){labels[g]=count++;anchors[g]=e;}out[j]=labels[g];out[11+j]=mod(e-anchors[g]);}return out;}
int main(int argc,char**argv){assert(argc==4);ifstream in(argv[1]);vector<Record>records;Record r;
 while(in>>r[0]){for(int i=1;i<67;i++)in>>r[i];assert(in);records.push_back(r);}assert(in.eof());
 ifstream pf(argv[2]);map<int,set<Key>>given;int id,c;Key k;while(pf>>id>>c){for(int&x:k)pf>>x;assert(pf);assert(c==1+*max_element(k.begin(),k.begin()+11));given[id].insert(k);}assert(given.size()==M*M);
 ifstream rf(argv[3]);map<int,set<Row>>actual;Row row;while(rf>>id){for(int&x:row)rf>>x;assert(rf);actual[id].insert(row);}assert(actual.size()==M*M);
 long long total_families=0,total_rows=0;
 for(int a=0;a<M;a++)for(int b=0;b<M;b++){set<Key>expected;
  for(const auto&r:records){bool ok=true;for(int i=0;i<12;i+=4)if(mod(a*r[i]+b*r[i+1]-r[i+3]))ok=false;if(ok)expected.insert(substitute(r,a,b));}
  assert(expected==given[M*a+b]);total_families+=expected.size();set<Row>expanded;
  for(const auto&f:expected){int free=*max_element(f.begin(),f.begin()+11);assert(free<=1);
   for(int phase=0;phase<(free?M:1);phase++){Row x;for(int j=0;j<11;j++)x[j]=mod(f[11+j]+phase*f[j]);expanded.insert(x);}}
  assert(expanded==actual[M*a+b]);total_rows+=expanded.size();
 }
 cout<<"records "<<records.size()<<" points "<<M*M<<" point_families "<<total_families<<" actual_rows "<<total_rows<<" exact_expansion_verified"<<endl;
}
