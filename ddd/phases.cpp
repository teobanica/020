// Exact Smith reduction, with an independent calculation from integer minors.
#include <algorithm>
#include <array>
#include <fstream>
#include <iostream>
#include <numeric>
#include <vector>
#include <cstdint>
#include <cstdlib>
using namespace std;using Mat=vector<vector<long long>>;
constexpr int M=7,N=6;
Mat ident(int n){Mat a(n,vector<long long>(n));for(int i=0;i<n;i++)a[i][i]=1;return a;}
void rowadd(Mat&a,int i,int j,long long q){for(int k=0;k<(int)a[0].size();k++)a[i][k]+=q*a[j][k];}
void coladd(Mat&a,int i,int j,long long q){for(auto&r:a)r[i]+=q*r[j];}
void colswap(Mat&a,int i,int j){for(auto&r:a)swap(r[i],r[j]);}
struct SNF{Mat d,u,v;int rank;};
SNF snf(Mat a){
 Mat u=ident(M),v=ident(N);int k=0;
 for(;k<N;k++){
  int pi=-1,pj=-1;for(int i=k;i<M;i++)for(int j=k;j<N;j++)if(a[i][j]&&(pi<0||llabs(a[i][j])<llabs(a[pi][pj]))){pi=i;pj=j;}
  if(pi<0)break;swap(a[k],a[pi]);swap(u[k],u[pi]);colswap(a,k,pj);colswap(v,k,pj);
  while(true){
   bool changed=false;
   for(int i=k+1;i<M;i++)if(a[i][k]){
    long long q=a[i][k]/a[k][k];rowadd(a,i,k,-q);rowadd(u,i,k,-q);
    if(a[i][k]){swap(a[i],a[k]);swap(u[i],u[k]);}changed=true;break;
   }if(changed)continue;
   for(int j=k+1;j<N;j++)if(a[k][j]){
    long long q=a[k][j]/a[k][k];coladd(a,j,k,-q);coladd(v,j,k,-q);
    if(a[k][j]){colswap(a,j,k);colswap(v,j,k);}changed=true;break;
   }if(changed)continue;
   int bi=-1;for(int i=k+1;i<M;i++)for(int j=k+1;j<N;j++)if(a[i][j]%a[k][k])bi=i;
   if(bi<0)break;rowadd(a,k,bi,1);rowadd(u,k,bi,1);
  }
  if(a[k][k]<0){for(auto&x:a[k])x=-x;for(auto&x:u[k])x=-x;}
 }return {a,u,v,k};
}
long long det(Mat a){
 int n=a.size();if(!n)return 1;long long prev=1,sign=1;
 for(int k=0;k<n-1;k++){
  if(!a[k][k]){int j=k+1;while(j<n&&!a[j][k])j++;if(j==n)return 0;swap(a[k],a[j]);sign=-sign;}
  long long p=a[k][k];
  for(int i=k+1;i<n;i++)for(int j=k+1;j<n;j++){
   long long num=p*a[i][j]-a[i][k]*a[k][j];if(num%prev){cerr<<"Nonexact Bareiss division\n";exit(10);}a[i][j]=num/prev;
  }
  for(int i=k+1;i<n;i++)a[i][k]=0;prev=p;
 }return sign*a[n-1][n-1];
}
vector<vector<int>> combos(int n,int k){vector<vector<int>>out;for(int mask=0;mask<(1<<n);mask++)if(__builtin_popcount((unsigned)mask)==k){vector<int>v;for(int j=0;j<n;j++)if(mask>>j&1)v.push_back(j);out.push_back(v);}return out;}
long long minor_gcd(const Mat&a,int k){
 if(k==0)return 1;long long g=0;
 for(auto r:combos(M,k))for(auto c:combos(N,k)){
  Mat b(k,vector<long long>(k));for(int i=0;i<k;i++)for(int j=0;j<k;j++)b[i][j]=a[r[i]][c[j]];
  g=gcd(g,llabs(det(b)));if(g==1)return 1;
 }return g;
}
void printmat(ostream&out,const Mat&a){out<<"[";for(int i=0;i<(int)a.size();i++){out<<(i?",":"")<<"[";for(int j=0;j<(int)a[i].size();j++)out<<(j?",":"")<<a[i][j];out<<"]";}out<<"]";}
int main(int argc,char**argv){
 if(argc!=3)return 2;ifstream in(argv[1]);ofstream out(argv[2]);int count;if(!(in>>count))return 3;out<<"[\n";long long feasible=0,components=0;
 for(int z=0;z<count;z++){
  int vi,wi;if(!(in>>vi>>wi))return 4;Mat a(M,vector<long long>(N));vector<long long>g(M);
  for(int i=0;i<M;i++){for(auto&x:a[i])if(!(in>>x))return 5;if(!(in>>g[i]))return 5;}
  auto s=snf(a);vector<long long>ug(M);for(int i=0;i<M;i++)for(int k=0;k<M;k++)ug[i]+=s.u[i][k]*g[k];
  for(int i=0;i<M;i++)for(int j=0;j<N;j++){long long q=0;for(int k=0;k<M;k++)for(int l=0;l<N;l++)q+=s.u[i][k]*a[k][l]*s.v[l][j];if(q!=s.d[i][j])return 6;}
  if(llabs(det(s.u))!=1||llabs(det(s.v))!=1)return 7;
  int rank=N;long long minor=minor_gcd(a,rank);while(rank&&!minor){rank--;minor=minor_gcd(a,rank);}
  long long maxd=rank?minor/minor_gcd(a,rank-1):1;long long torsion=1;
  for(int i=0;i<s.rank;i++)torsion*=s.d[i][i];
  if(rank!=s.rank||minor!=torsion||maxd!=(rank?s.d[rank-1][rank-1]:1))return 8;
  bool ok=true;for(int i=rank;i<M;i++)if(ug[i]%2)ok=false;
  feasible+=ok;components+=ok*torsion;
  out<<"{\"vi\":"<<vi<<",\"wi\":"<<wi<<",\"rank\":"<<rank<<",\"binary_feasible\":"<<(ok?"true":"false")<<",\"diagonal\":[";
  for(int i=0;i<rank;i++)out<<(i?",":"")<<s.d[i][i];out<<"],\"U\":";printmat(out,s.u);out<<",\"V\":";printmat(out,s.v);out<<"}"<<(z+1<count?",":"")<<"\n";
  if(z%50000==0)cerr<<z<<"/"<<count<<" binary feasible="<<feasible<<"\n";
 }out<<"]\n";cerr<<"Binary feasible supports="<<feasible<<" components="<<components<<"\n";
}
