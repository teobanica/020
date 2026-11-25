// Exact canonicalization of 5x5x5 intersection tensors.
// First minimize a two-axis marginal, then minimize the tensor among
// transformations realizing that marginal. Final-axis pair slices are sorted.
#include <algorithm>
#include <array>
#include <fstream>
#include <iostream>
#include <map>
#include <unordered_map>
#include <vector>
#include <string>
#include <cassert>
using namespace std;
vector<array<int,3>> axes;
vector<array<int,5>> perms;
struct PairCanon { string key; vector<array<int,2>> choices; };
unordered_map<string,PairCanon> cache;
PairCanon& pc(const string&m){
 auto it=cache.find(m);if(it!=cache.end())return it->second;
 PairCanon r;r.key=string(25,'9');
 for(int a=0;a<24;a++)for(int b=0;b<24;b++){
  string k(25,'0');for(int i=0;i<5;i++)for(int j=0;j<5;j++)k[5*i+j]=m[5*perms[a][i]+perms[b][j]];
  if(k<r.key){r.key=k;r.choices.clear();r.choices.push_back({a,b});}
  else if(k==r.key)r.choices.push_back({a,b});
 }
 return cache.emplace(m,move(r)).first->second;
}
int ix(const array<int,3>&ax,int i,int j,int k){int q[3];q[ax[0]]=i;q[ax[1]]=j;q[ax[2]]=k;return 25*q[0]+5*q[1]+q[2];}
pair<string,int> canon(const string&t){
 string bestM(25,'9'),best(125,'9');int witness=-1;
 for(int ai=0;ai<6;ai++){auto ax=axes[ai];
  string m(25,'0');for(int i=0;i<5;i++)for(int j=0;j<5;j++)for(int k=0;k<5;k++)m[5*i+j]+=t[ix(ax,i,j,k)]-'0';
  PairCanon&r=pc(m);
  if(r.key>bestM)continue;
  if(r.key<bestM){bestM=r.key;best=string(125,'9');}
  for(auto ab:r.choices){
   auto a=perms[ab[0]],b=perms[ab[1]];array<string,5>cols;for(auto&s:cols)s.reserve(25);
   for(int i=0;i<5;i++)for(int j=0;j<5;j++)for(int k=0;k<5;k++)cols[k].push_back(t[ix(ax,a[i],b[j],k)]);
   array<int,5> order{0,1,2,3,4};sort(order.begin()+1,order.end(),[&](int x,int y){return cols[x]<cols[y];});string v;v.reserve(125);
   for(int i=0;i<25;i++)for(int k=0;k<5;k++)v.push_back(cols[order[k]][i]);
   if(v<best){best=v;int ci=find(perms.begin(),perms.end(),order)-perms.begin();assert(ci<24);witness=ai*13824+ab[0]*576+ab[1]*24+ci;}
  }
 }
 assert(witness>=0);return {best,witness};
}
struct Record{int vi,wi,dim,maxd,tor,ori,count,raw;};
int main(int argc,char**argv){
 if(argc!=4)return 2;
 array<int,3>a{0,1,2};do{axes.push_back(a);}while(next_permutation(a.begin(),a.end()));
 array<int,5>b{0,1,2,3,4};do{perms.push_back(b);}while(next_permutation(b.begin()+1,b.end()));assert(perms.size()==24);
 ifstream in(argv[1]);ofstream out(argv[2]),witnesses(argv[3]);map<string,Record>classes;string t;Record r;int z=0;
 while(in>>t>>r.vi>>r.wi>>r.dim>>r.maxd>>r.tor>>r.ori>>r.count){
  assert(t.size()==125);r.raw=1;auto [c,w]=canon(t);witnesses<<w<<'\n';auto q=classes.find(c);
  if(q==classes.end())classes[c]=r;
  else{auto&s=q->second;assert(s.dim==r.dim&&s.maxd==r.maxd&&s.tor==r.tor&&s.ori==r.ori);s.count+=r.count;s.raw++;}
  if(++z%20000==0)cerr<<"Raw tensors "<<z<<" classes "<<classes.size()<<" pair cache "<<cache.size()<<"\n";
 }
 map<int,int>dims,tors,maxds,oris;int count=0;
 out<<"{\"classes\":[\n";int n=0;
 for(auto&[t,r]:classes){
  out<<"{\"tensor\":\""<<t<<"\",\"vi\":"<<r.vi<<",\"wi\":"<<r.wi<<",\"dimension\":"<<r.dim<<",\"largest_smith\":"<<r.maxd<<",\"torsion\":"<<r.tor<<",\"orientations\":"<<r.ori<<",\"presentations\":"<<r.count<<",\"raw_tensors\":"<<r.raw<<"}"<<(++n<(int)classes.size()?",":"")<<"\n";
  dims[r.dim]++;tors[r.tor]++;maxds[r.maxd]++;oris[r.ori]++;count+=r.count;
 }
 out<<"],\"feasible_support_classes\":"<<classes.size()<<",\"raw_tensors\":"<<z<<",\"presentations\":"<<count<<"}\n";
 cerr<<"Classes "<<classes.size()<<" presentations "<<count<<"\n";
 for(auto[k,v]:dims)cerr<<"dimension "<<k<<": "<<v<<"\n";
 for(auto[k,v]:maxds)cerr<<"largest Smith "<<k<<": "<<v<<"\n";
}

