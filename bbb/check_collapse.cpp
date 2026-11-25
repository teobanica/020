// Independent finite-grid check, with Smith invariant bounds from exact minors.
#include <algorithm>
#include <fstream>
#include <iostream>
#include <numeric>
#include <vector>
#include <cstdlib>
using namespace std;
using Mat=vector<vector<long long>>;
long long det(Mat a){
 int n=a.size();if(!n)return 1;long long prev=1,sgn=1;
 for(int k=0;k<n-1;k++){
  if(!a[k][k]){int j=k+1;while(j<n&&!a[j][k])j++;if(j==n)return 0;swap(a[k],a[j]);sgn=-sgn;}
  long long p=a[k][k];
  for(int i=k+1;i<n;i++)for(int j=k+1;j<n;j++)a[i][j]=(p*a[i][j]-a[i][k]*a[k][j])/prev;
  for(int i=k+1;i<n;i++)a[i][k]=0;
  prev=p;
 }
 return sgn*a[n-1][n-1];
}
vector<vector<int>> combos(int n,int k){
 vector<vector<int>> out;for(int m=0;m<(1<<n);m++)if(__builtin_popcount((unsigned)m)==k){
  vector<int> v;for(int j=0;j<n;j++)if(m>>j&1)v.push_back(j);out.push_back(v);
 }return out;
}
long long minor_gcd(const Mat&a,int k){
 if(k==0)return 1;long long g=0;
 for(auto r:combos(a.size(),k))for(auto c:combos(4,k)){
  Mat b(k,vector<long long>(k));for(int i=0;i<k;i++)for(int j=0;j<k;j++)b[i][j]=a[r[i]][c[j]];
  g=gcd(g,llabs(det(b)));if(g==1)return 1;
 }return g;
}
int main(int argc,char**argv){
 if(argc!=3)return 2;ifstream in(argv[1]);ofstream out(argv[2]);int n;in>>n;out<<"[\n";
 for(int z=0;z<n;z++){
  int vi,wi,m;in>>vi>>wi>>m;Mat a(m,vector<long long>(4));vector<int>b(m);
  for(int i=0;i<m;i++){for(auto&v:a[i])in>>v;in>>b[i];}
  int rank=4;long long g=minor_gcd(a,rank);while(rank&&!g){rank--;g=minor_gcd(a,rank);}
  long long maxd=rank?g/minor_gcd(a,rank-1):1;
  if(4%maxd){cerr<<"Grid bound not established "<<maxd<<"\n";return 3;}
  int count=0;
  for(int k=0;k<4096;k++){
   int x[4]={k%8,k/8%8,k/64%8,k/512};bool ok=true;
   for(int i=0;i<m;i++){long long s=-4*b[i];for(int j=0;j<4;j++)s+=a[i][j]*x[j];if(s%8){ok=false;break;}}
   count+=ok;
  }
  out<<"{\"vi\":"<<vi<<",\"wi\":"<<wi<<",\"rank\":"<<rank<<",\"max_smith_invariant\":"<<maxd<<",\"grid_count\":"<<count<<"}"<<(z+1<n?",":"")<<"\n";
 }out<<"]\n";cerr<<"Checked "<<n<<" systems on eighth-root grid with exact Smith bounds.\n";
}
