// Independent full-group minimum and transformation-witness verification.
// Group: S3 on quotient axes and S4 on the pair labels of each axis.
// A marginal comparison prunes transformations that cannot be minimal;
// all 24 final-axis permutations are explicitly tested when needed.
#include <algorithm>
#include <array>
#include <fstream>
#include <iostream>
#include <map>
#include <string>
#include <vector>
#include <cassert>
using namespace std;
vector<array<int,3>> axes;
vector<array<int,5>> perms;
int loc(array<int,3> ax,int i,int j,int k){array<int,3>q;q[ax[0]]=i;q[ax[1]]=j;q[ax[2]]=k;return 25*q[0]+5*q[1]+q[2];}
string transform(const string&t,int w){
 int ai=w/13824;w%=13824;int a=w/576;w%=576;int b=w/24,c=w%24;
 string s(125,'0');for(int i=0;i<5;i++)for(int j=0;j<5;j++)for(int k=0;k<5;k++)s[25*i+5*j+k]=t[loc(axes[ai],perms[a][i],perms[b][j],perms[c][k])];return s;
}
string independent(const string&t){
 string bestM(25,'9'),bestT(125,'9');
 for(auto axis:axes){
  int marginal[5][5]{};for(int i=0;i<5;i++)for(int j=0;j<5;j++)for(int k=0;k<5;k++)marginal[i][j]+=t[loc(axis,i,j,k)]-'0';
  for(auto pa:perms)for(auto pb:perms){
   string m(25,'0');for(int i=0;i<5;i++)for(int j=0;j<5;j++)m[5*i+j]+=marginal[pa[i]][pb[j]];
   if(m>bestM)continue;
   if(m<bestM){bestM=m;bestT=string(125,'9');}
   for(auto pc:perms){
    string u(125,'0');for(int i=0;i<5;i++)for(int j=0;j<5;j++)for(int k=0;k<5;k++)u[25*i+5*j+k]=t[loc(axis,pa[i],pb[j],pc[k])];
    if(u<bestT)bestT=u;
   }
  }
 }
 return bestT;
}
struct Rec{int dim,maxd,tor,ori,n,raw;long long seen=0,seenRaw=0;};
int main(int argc,char**argv){
 if(argc!=5)return 2;
 array<int,3>a{0,1,2};do{axes.push_back(a);}while(next_permutation(a.begin(),a.end()));
 array<int,5>b{0,1,2,3,4};do{perms.push_back(b);}while(next_permutation(b.begin()+1,b.end()));
 ifstream classes(argv[1]),raw(argv[2]),witness(argv[3]);string t;Rec r;map<string,Rec>dict;map<int,int>dims;int z=0;
 while(classes>>t>>r.dim>>r.maxd>>r.tor>>r.ori>>r.n>>r.raw){
  assert(t.size()==125);assert(independent(t)==t);assert(!dict.count(t));dict[t]=r;dims[r.dim]++;
  if(++z%2000==0)cerr<<"Full-group minima checked "<<z<<"\n";
 }
 long long count=0,rawCount=0;
 int vi,wi,dim,maxd,tor,ori,n,w;
 while(raw>>t>>vi>>wi>>dim>>maxd>>tor>>ori>>n){
  assert(witness>>w);assert(w>=0&&w<82944);string c=transform(t,w);auto it=dict.find(c);assert(it!=dict.end());
  Rec&v=it->second;assert(v.dim==dim&&v.maxd==maxd&&v.tor==tor&&v.ori==ori);v.seen+=n;v.seenRaw++;count+=n;rawCount++;
 }
 assert(!(witness>>w));for(auto&[s,v]:dict)assert(v.seen==v.n&&v.seenRaw==v.raw);
 ofstream out(argv[4]);out<<"{\"verified\":true,\"group_order\":82944,\"classes\":"<<dict.size()<<",\"raw_tensors\":"<<rawCount<<",\"presentations\":"<<count<<",\"dimensions\":{";int k=0;for(auto[d,n]:dims)out<<(k++?",":"")<<"\""<<d<<"\":"<<n;out<<"}}\n";
 cerr<<"Verified "<<dict.size()<<" classes, "<<rawCount<<" transformation witnesses, "<<count<<" presentations\n";
}

