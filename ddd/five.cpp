#include <array>
#include <algorithm>
#include <fstream>
#include <iostream>
#include <vector>
#include <cstdint>
using namespace std;
constexpr int NVAR=24,NAUG=25;
using D=array<unsigned,NVAR>;
struct Problem {int vi,wi; array<array<int,3>,11> eq;};
uint64_t nodes;
int pc(unsigned x){return __builtin_popcount(x);}
vector<int> linear_collision(const Problem&p){
 vector<array<int,NAUG>> a;
 for(auto e:p.eq){array<int,NAUG> r{};r[e[0]]=1;r[e[1]]=4;r[e[2]]=4;a.push_back(r);}
 for(int b:{0,8,16}){array<int,NAUG>r{};for(int j=0;j<5;j++)r[b+j]=1;a.push_back(r);}
 for(auto f:vector<pair<int,int>>{{0,0},{8,0},{1,1}}){array<int,NAUG>r{};r[f.first]=1;r[NVAR]=f.second;a.push_back(r);}
 vector<int> piv;int r=0;int inv[5]={0,1,3,2,4};
 for(int c=0;c<NVAR;c++){
  int k=r;while(k<(int)a.size()&&!a[k][c])k++;if(k==(int)a.size())continue;swap(a[r],a[k]);
  int z=inv[a[r][c]];for(int j=0;j<NAUG;j++)a[r][j]=a[r][j]*z%5;
  for(int i=0;i<(int)a.size();i++)if(i!=r&&a[i][c]){z=a[i][c];for(int j=0;j<NAUG;j++)a[i][j]=(a[i][j]-z*a[r][j]+25)%5;}
  piv.push_back(c);r++;
 }
 for(auto row:a){bool zero=true;for(int j=0;j<NVAR;j++)zero &= row[j]==0;if(zero&&row[NVAR])return {-1,-1};}
 array<array<int,NAUG>,NVAR> form{};
 for(int c=0;c<NVAR;c++)form[c][c]=1;
 for(int i=0;i<r;i++){
  int c=piv[i];form[c].fill(0);form[c][NVAR]=a[i][NVAR];
  for(int j=0;j<NVAR;j++)if(j!=c)form[c][j]=(5-a[i][j])%5;
 }
 for(int b:{0,8,16})for(int i=0;i<5;i++)for(int j=i+1;j<5;j++)if(form[b+i]==form[b+j])return {b+i,b+j};
 return {};
}
bool propagate(const Problem&p,D&d){
 bool change=true;
 while(change){change=false;
  for(auto e:p.eq){unsigned a=0,b=0,c=0;
   for(int y=0;y<5;y++)if(d[e[1]]>>y&1)for(int z=0;z<5;z++)if(d[e[2]]>>z&1){int x=(y+z)%5;if(d[e[0]]>>x&1){a|=1u<<x;b|=1u<<y;c|=1u<<z;}}
   if(!a||!b||!c)return false;
   for(int k=0;k<3;k++){unsigned m=k==0?a:k==1?b:c;if(m!=d[e[k]]){d[e[k]]=m;change=true;}}
  }
  for(int base:{0,8,16}){array<unsigned,32>un{};
   for(int mask=1;mask<32;mask++){int k=__builtin_ctz((unsigned)mask);un[mask]=un[mask^(1<<k)]|d[base+k];int n=pc(mask),m=pc(un[mask]);
    if(m<n)return false;
    if(m==n&&mask!=31)for(int j=0;j<5;j++)if(!(mask>>j&1)){unsigned q=d[base+j]&~un[mask];if(!q)return false;if(q!=d[base+j]){d[base+j]=q;change=true;}}
   }
  }
 }return true;
}
bool solve(const Problem&p,D&d){
 nodes++;if(!propagate(p,d))return false;int v=-1,best=6;
 for(int i=0;i<NVAR;i++)if(pc(d[i])>1&&pc(d[i])<best){best=pc(d[i]);v=i;}
 if(v<0)return true;
 for(int x=0;x<5;x++)if(d[v]>>x&1){D ch=d;ch[v]=1u<<x;if(solve(p,ch)){d=ch;return true;}}
 return false;
}
int main(int argc,char**argv){
 if(argc!=3)return 2;ifstream in(argv[1]);ofstream out(argv[2]);int n;if(!(in>>n))return 3;
 int sat=0,linear=0,unsat=0;out<<"[\n";
 for(int z=0;z<n;z++){
  Problem p;if(!(in>>p.vi>>p.wi))return 4;for(auto&e:p.eq)for(auto&j:e)if(!(in>>j)||j<0||j>=NVAR)return 5;
  auto cert=linear_collision(p);D d;d.fill(31);d[0]=1;d[8]=1;d[1]=2;nodes=0;
  bool ok=cert.empty()&&solve(p,d);if(ok)sat++;else if(!cert.empty())linear++;else unsat++;
  out<<"{\"vi\":"<<p.vi<<",\"wi\":"<<p.wi<<",\"satisfiable\":"<<(ok?"true":"false");
  if(!cert.empty())out<<",\"collision\":["<<cert[0]<<","<<cert[1]<<"]";
  if(ok){out<<",\"labels\":[";for(int i=0;i<NVAR;i++)out<<(i?",":"")<<__builtin_ctz(d[i]);out<<"]";}
  out<<",\"nodes\":"<<nodes<<"}"<<(z+1<n?",":"")<<"\n";
  if(z%10000==0)cerr<<z<<"/"<<n<<" sat="<<sat<<" linear="<<linear<<" search_unsat="<<unsat<<"\n";
 }out<<"]\n";cerr<<"DONE sat="<<sat<<" linear="<<linear<<" search_unsat="<<unsat<<"\n";
}
