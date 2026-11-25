// Canonical 4x4x4 support tensors under quotient permutations and pair relabeling.
// Optimize the final-axis permutation by lexicographically sorting its 3 slices.
#include <algorithm>
#include <array>
#include <fstream>
#include <iostream>
#include <map>
#include <unordered_map>
#include <vector>
#include <string>
using namespace std;
vector<array<int,3>> perms;
string canon(const string&t){
 string best(64,char(127));
 for(auto axis:perms)for(auto pa:perms)for(auto pb:perms){
  int a[4]={0,pa[0]+1,pa[1]+1,pa[2]+1},b[4]={0,pb[0]+1,pb[1]+1,pb[2]+1};
  array<string,4> cols;for(auto&s:cols)s.reserve(16);
  for(int i=0;i<4;i++)for(int j=0;j<4;j++)for(int k=0;k<4;k++){
   int q[3];q[axis[0]]=a[i];q[axis[1]]=b[j];q[axis[2]]=k;
   cols[k].push_back(t[16*q[0]+4*q[1]+q[2]]);
  }
  sort(cols.begin()+1,cols.end());string v;v.reserve(64);
  for(int i=0;i<16;i++)for(int k=0;k<4;k++)v.push_back(cols[k][i]);
  if(v<best)best=v;
 }return best;
}
struct Record{int vi,wi,dim,torsion,ok,count;};
int main(int argc,char**argv){
 if(argc!=3)return 2;array<int,3>a{0,1,2};do{perms.push_back(a);}while(next_permutation(a.begin(),a.end()));
 ifstream in(argv[1]);ofstream out(argv[2]);int n;if(!(in>>n))return 3;
 unordered_map<string,string>cache;map<string,Record>classes;
 for(int z=0;z<n;z++){
  Record r;if(!(in>>r.vi>>r.wi>>r.dim>>r.torsion>>r.ok))return 4;r.count=1;
  string t(64,0);for(int j=0;j<11;j++){int c;if(!(in>>c)||c<0||c>=64)return 5;t[c]++;}
  auto p=cache.find(t);if(p==cache.end())p=cache.emplace(t,canon(t)).first;
  auto q=classes.find(p->second);
  if(q==classes.end())classes[p->second]=r;
  else{auto&s=q->second;if(s.dim!=r.dim||s.torsion!=r.torsion||s.ok!=r.ok)return 6;s.count++;}
 }
 map<int,int>dims,torsions;int feasible=0;for(auto&[t,r]:classes)if(r.ok){feasible++;dims[r.dim]++;torsions[r.torsion]++;}
 out<<"{\"total_support_classes\":"<<classes.size()<<",\"feasible_support_classes\":"<<feasible<<",\"classes\":[\n";int z=0;
 for(auto&[t,r]:classes){
  out<<"{\"tensor\":[";for(int j=0;j<64;j++)out<<(j?",":"")<<(int)(unsigned char)t[j];
  out<<"],\"representative\":["<<r.vi<<","<<r.wi<<"],\"dimension\":"<<r.dim<<",\"torsion_order\":"<<r.torsion<<",\"feasible\":"<<(r.ok?"true":"false")<<",\"presentations\":"<<r.count<<"}"<<(++z<(int)classes.size()?",":"")<<"\n";
 }out<<"]}\n";
 cerr<<"Classes="<<classes.size()<<" feasible="<<feasible<<" raw="<<cache.size()<<"\n";
 for(auto[d,c]:dims)cerr<<"dimension "<<d<<": "<<c<<"\n";
 for(auto[d,c]:torsions)cerr<<"torsion "<<d<<": "<<c<<"\n";
}
