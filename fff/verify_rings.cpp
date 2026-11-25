// Independent finite-ring verification of every support presentation.
// Solve A*z=12*(3*f+2*B*sigma) over Z/72 using Z/8 and Z/9.
// All Smith divisors from the separately certified calculation divide 12.
#include <array>
#include <vector>
#include <fstream>
#include <iostream>
#include <cassert>
#include <algorithm>
#include <map>
using namespace std;
using Part=array<int,11>;
struct Info{array<int,11>block,pos;};
Info info(const Part&p){
 Info r;for(int k=0;k<11;k++){int b=k<3?0:1+(k-3)/2;int h=k<3?k:(k-3)%2;r.block[p[k]]=b;r.pos[p[k]]=h;}return r;
}
int mod(int x,int m){x%=m;return x<0?x+m:x;}
struct Ring{int mask,rank,tor,maxd;};
Ring solve(array<array<int,8>,6>A,array<array<int,8>,6>B,int p,int m){
 for(auto&r:A)for(int&x:r)x=mod(x,m);for(auto&r:B)for(int&x:r)x=mod(x,m);
 int rank=0,tor=1,maxd=1;
 for(int k=0;k<6;k++){
  int pi=-1,pj=-1,best=m;
  for(int i=k;i<6;i++)for(int j=k;j<8;j++)if(A[i][j]){
   int d=1,x=A[i][j];while(x%p==0){d*=p;x/=p;}
   if(d<best){best=d;pi=i;pj=j;}
  }
  if(pi<0)break;
  swap(A[k],A[pi]);swap(B[k],B[pi]);for(auto&r:A)swap(r[k],r[pj]);
  int unit=A[k][k]/best,inv=1;while(mod(unit*inv,m)!=1)inv++;assert(inv<m);
  for(int&x:A[k])x=mod(x*inv,m);for(int&x:B[k])x=mod(x*inv,m);assert(A[k][k]==best);
  for(int i=k+1;i<6;i++){
   assert(A[i][k]%best==0);int q=A[i][k]/best;
   for(int j=0;j<8;j++){A[i][j]=mod(A[i][j]-q*A[k][j],m);B[i][j]=mod(B[i][j]-q*B[k][j],m);}
  }
  for(int j=k+1;j<8;j++){
   assert(A[k][j]%best==0);int q=A[k][j]/best;
   for(int i=0;i<6;i++)A[i][j]=mod(A[i][j]-q*A[i][k],m);
  }
  rank++;tor*=best;maxd=best;
 }
 int mask=0;
 for(int c=0;c<8;c++){
  bool ok=true;for(int i=0;i<6;i++){if(i<rank){if(B[i][c]%A[i][i])ok=false;}else if(B[i][c])ok=false;}
  if(ok)mask|=1<<c;
 }
 return {mask,rank,tor,maxd};
}
struct Expect{int dim=0,maxd=0,tor=0,mask=0;};
int main(int argc,char**argv){
 if(argc!=3)return 2;ifstream in(argv[1]);int nq,nr,nf;in>>nq>>nr>>nf;
 vector<Part>qs(nq),rs(nr);for(auto&p:qs)for(int&x:p)in>>x;for(auto&p:rs)for(int&x:p)in>>x;
 vector<Expect>ex(nq*nr);for(int z=0;z<nf;z++){int vi,wi;Expect e;in>>vi>>wi>>e.dim>>e.maxd>>e.tor>>e.mask;assert(!ex[vi*nr+wi].mask);ex[vi*nr+wi]=e;}
 Part base{0,1,2,3,4,5,6,7,8,9,10};Info ip=info(base);long long feasible=0,orientations=0,components=0;map<int,int>dims;
 for(int vi=0;vi<nq;vi++){
  Info iq=info(qs[vi]);
  for(int wi=0;wi<nr;wi++){
   Info ir=info(rs[wi]);array<array<int,8>,11>rows{},offset{};
   // Build the phase of u/v directly, for each triangle orientation.
   for(int j=0;j<11;j++){
    if(ip.block[j])rows[j][ip.block[j]-1]++;
    if(iq.block[j])rows[j][iq.block[j]+3]--;
    for(int c=0;c<8;c++){
     int sp=(c&1)?1:-1,sq=(c&2)?1:-1,sr=(c&4)?1:-1;
     int up=ip.block[j]?3*ip.pos[j]:2*sp*ip.pos[j];
     int vq=iq.block[j]?3*iq.pos[j]:2*sq*iq.pos[j];
     int wr=ir.block[j]?3*ir.pos[j]:2*sr*ir.pos[j];
     offset[j][c]=wr-up+vq;
    }
   }
   array<array<int,8>,6>A{},B{};int row=0;
   for(int start=0;start<11;){
    int size=start==0?3:2,ref=rs[wi][start];
    for(int h=1;h<size;h++){
     int j=rs[wi][start+h];for(int k=0;k<8;k++){A[row][k]=rows[j][k]-rows[ref][k];B[row][k]=12*(offset[j][k]-offset[ref][k]);}row++;
    }start+=size;
   }
   assert(row==6);Ring a=solve(A,B,2,8),b=solve(A,B,3,9);assert(a.rank==b.rank);
   int mask=a.mask&b.mask;Expect e=ex[vi*nr+wi];assert(mask==e.mask);
   if(mask){
    assert(8-a.rank==e.dim&&a.tor*b.tor==e.tor&&a.maxd*b.maxd==e.maxd);
    feasible++;int n=__builtin_popcount((unsigned)mask);orientations+=n;components+=(long long)n*e.tor;dims[e.dim]++;
   }
  }
 }
 assert(feasible==nf);ofstream out(argv[2]);out<<"{\"verified\":true,\"method\":\"independent elimination over Z/8 and Z/9\",\"presentations\":"<<nq*nr<<",\"feasible_presentations\":"<<feasible<<",\"orientations\":"<<orientations<<",\"phase_components\":"<<components<<"}\n";
 cerr<<"Verified "<<nq*nr<<" presentations; feasible "<<feasible<<" orientations "<<orientations<<" components "<<components<<"\n";
}

