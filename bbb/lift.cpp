#include <array>
#include <algorithm>
#include <fstream>
#include <iostream>
#include <vector>
#include <cstdint>
using namespace std;
using D=array<unsigned,27>;
struct Problem {int vi,wi; array<array<int,3>,11> eq; int u0,v0,u1;};
uint64_t nodes;
int pc(unsigned x){return __builtin_popcount(x);}
bool propagate(const Problem& p,D& d){
 bool change=true;
 while(change){
  change=false;
  for(auto e:p.eq){
   unsigned a=0,b=0,c=0;
   for(int y=0;y<7;y++) if(d[e[1]]>>y&1)
    for(int z=0;z<7;z++) if(d[e[2]]>>z&1){
     int x=(y+z)%7;
     if(d[e[0]]>>x&1){a|=1u<<x;b|=1u<<y;c|=1u<<z;}
    }
   if(!a||!b||!c)return false;
   for(int k=0;k<3;k++){
    unsigned m=k==0?a:k==1?b:c;
    if(m!=d[e[k]]){d[e[k]]=m;change=true;}
   }
  }
  for(int base:{0,9,18}){
   // Hall sets: every subset of these seven domains needs at least its size
   // in possible residues. Equality reserves those residues for that subset.
   array<unsigned,128> un{};
   for(int mask=1;mask<128;mask++){
    int k=__builtin_ctz((unsigned)mask);
    un[mask]=un[mask^(1<<k)]|d[base+k];
    int n=pc(mask),m=pc(un[mask]);
    if(m<n)return false;
    if(m==n&&mask!=127){
     for(int j=0;j<7;j++)if(!(mask>>j&1)){
      unsigned q=d[base+j]&~un[mask];
      if(!q)return false;
      if(q!=d[base+j]){d[base+j]=q;change=true;}
     }
    }
   }
  }
 }
 return true;
}
bool solve(const Problem& p,D& d){
 nodes++;
 if(!propagate(p,d))return false;
 int v=-1,best=8;
 for(int i=0;i<27;i++)if(pc(d[i])>1&&pc(d[i])<best){best=pc(d[i]);v=i;}
 if(v<0)return true;
 unsigned choices=d[v];
 for(int x=0;x<7;x++)if(choices>>x&1){D child=d;child[v]=1u<<x;
  if(solve(p,child)){d=child;return true;}
 }
 return false;
}
int main(int argc,char**argv){
 if(argc!=3)return 2;
 ifstream in(argv[1]);ofstream out(argv[2]);int count;if(!(in>>count)||count<0)return 3;
 int sat=0,unsat=0;uint64_t total=0;
 out<<"[\n";
 for(int k=0;k<count;k++){
  Problem p;if(!(in>>p.vi>>p.wi>>p.u0>>p.v0>>p.u1))return 3;
  for(auto& e:p.eq){if(!(in>>e[0]>>e[1]>>e[2]))return 3;for(int i:e)if(i<0||i>=27)return 3;}
  for(int i:{p.u0,p.v0,p.u1})if(i<0||i>=27)return 3;
  D d;d.fill(127);d[p.u0]=1;d[p.v0]=1;d[p.u1]=2;
  nodes=0;bool ok=solve(p,d);total+=nodes;(ok?sat:unsat)++;
  out<<"{\"vi\":"<<p.vi<<",\"wi\":"<<p.wi<<",\"satisfiable\":"<<(ok?"true":"false")<<",\"nodes\":"<<nodes;
  if(ok){out<<",\"labels\":[";for(int i=0;i<27;i++)out<<(i?",":"")<<__builtin_ctz(d[i]);out<<"]";}
  out<<"}"<<(k+1<count?",":"")<<"\n";out.flush();
  if(k%100==0)cerr<<k<<"/"<<count<<" sat="<<sat<<" unsat="<<unsat<<" nodes="<<total<<"\n";
 }
 out<<"]\n";cerr<<"DONE sat="<<sat<<" unsat="<<unsat<<" nodes="<<total<<"\n";
}
