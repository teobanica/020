// Integer Smith reduction with explicit unimodular row and column transforms.
#include <algorithm>
#include <array>
#include <fstream>
#include <iostream>
#include <vector>
#include <cstdint>
#include <cstdlib>
using namespace std;using Mat=vector<vector<long long>>;
Mat ident(int n){Mat a(n,vector<long long>(n));for(int i=0;i<n;i++)a[i][i]=1;return a;}
void rowadd(Mat&a,int i,int j,long long q){for(int k=0;k<(int)a[0].size();k++)a[i][k]+=q*a[j][k];}
void coladd(Mat&a,int i,int j,long long q){for(auto&r:a)r[i]+=q*r[j];}
void colswap(Mat&a,int i,int j){for(auto&r:a)swap(r[i],r[j]);}
struct SNF{Mat d,u,v;int rank;};
SNF snf(Mat a){
 Mat u=ident(8),v=ident(4);int k=0;
 for(;k<4;k++){
  int pi=-1,pj=-1;for(int i=k;i<8;i++)for(int j=k;j<4;j++)if(a[i][j]&&(pi<0||llabs(a[i][j])<llabs(a[pi][pj]))){pi=i;pj=j;}
  if(pi<0)break;swap(a[k],a[pi]);swap(u[k],u[pi]);colswap(a,k,pj);colswap(v,k,pj);
  while(true){
   bool changed=false;
   for(int i=k+1;i<8;i++)if(a[i][k]){
    long long q=a[i][k]/a[k][k];rowadd(a,i,k,-q);rowadd(u,i,k,-q);
    if(a[i][k]){swap(a[i],a[k]);swap(u[i],u[k]);}changed=true;break;
   }if(changed)continue;
   for(int j=k+1;j<4;j++)if(a[k][j]){
    long long q=a[k][j]/a[k][k];coladd(a,j,k,-q);coladd(v,j,k,-q);
    if(a[k][j]){colswap(a,j,k);colswap(v,j,k);}changed=true;break;
   }if(changed)continue;
   int bi=-1;for(int i=k+1;i<8;i++)for(int j=k+1;j<4;j++)if(a[i][j]%a[k][k])bi=i;
   if(bi<0)break;rowadd(a,k,bi,1);rowadd(u,k,bi,1);
  }
  if(a[k][k]<0){for(auto&x:a[k])x=-x;for(auto&x:u[k])x=-x;}
 }
 return {a,u,v,k};
}
void printmat(ostream&out,const Mat&a){out<<"[";for(int i=0;i<(int)a.size();i++){out<<(i?",":"")<<"[";for(int j=0;j<(int)a[i].size();j++)out<<(j?",":"")<<a[i][j];out<<"]";}out<<"]";}
int main(int argc,char**argv){
 if(argc!=3)return 2;ifstream in(argv[1]);ofstream out(argv[2]);int n;if(!(in>>n))return 3;out<<"[\n";long long supports=0,orients=0,components=0;
 for(int z=0;z<n;z++){
  int vi,wi;if(!(in>>vi>>wi))return 4;Mat a(8,vector<long long>(4)),b(8,vector<long long>(6));
  for(int i=0;i<8;i++){for(auto&x:a[i])if(!(in>>x))return 5;for(auto&x:b[i])if(!(in>>x))return 5;}
  auto s=snf(a);Mat ub(8,vector<long long>(6));
  for(int i=0;i<8;i++)for(int j=0;j<6;j++)for(int k=0;k<8;k++)ub[i][j]+=s.u[i][k]*b[k][j];
  // Exact transformation identity, including all off-diagonal entries.
  for(int i=0;i<8;i++)for(int j=0;j<4;j++){long long q=0;for(int k=0;k<8;k++)for(int l=0;l<4;l++)q+=s.u[i][k]*a[k][l]*s.v[l][j];if(q!=s.d[i][j])return 6;}
  vector<int> masks;long long torsion=1;for(int i=0;i<s.rank;i++)torsion*=s.d[i][i];
  for(int mask=0;mask<64;mask++){
   bool ok=true;for(int i=s.rank;i<8;i++){long long q=0;for(int j=0;j<6;j++)q+=ub[i][j]*(mask>>j&1?1:-1);if(q%3)ok=false;}
   if(ok)masks.push_back(mask);
  }
  supports+=!masks.empty();orients+=masks.size();components+=masks.size()*torsion;
  out<<"{\"vi\":"<<vi<<",\"wi\":"<<wi<<",\"rank\":"<<s.rank<<",\"diagonal\":[";
  for(int i=0;i<s.rank;i++)out<<(i?",":"")<<s.d[i][i];out<<"],\"U\":";printmat(out,s.u);out<<",\"V\":";printmat(out,s.v);
  out<<",\"orientations\":[";for(int i=0;i<(int)masks.size();i++)out<<(i?",":"")<<masks[i];out<<"]}"<<(z+1<n?",":"")<<"\n";
 }out<<"]\n";cerr<<"Feasible supports="<<supports<<" orientations="<<orients<<" components="<<components<<"\n";
}
