#include <algorithm>
#include <array>
#include <cassert>
#include <fstream>
#include <iostream>
#include <numeric>
#include <vector>
using namespace std;
using Mat = array<array<int,5>,4>;
const int M=420;
int H[5][11]={
 {0,0,0,0,0,0,0,0,0,0,0},
 {0,4,8,12,16,0,10,15,5,5,15},
 {0,8,16,4,12,0,15,10,15,5,5},
 {0,12,4,16,8,0,5,15,10,15,5},
 {0,16,12,8,4,0,5,5,15,10,15}};
struct Field {
 int p,g,z; vector<int> pw,inv;
 int mul(int a,int b) const {return int(1LL*a*b%p);}
 int power(int a,int n) const {int b=1;for(;n;n>>=1,a=mul(a,a))if(n&1)b=mul(a,b);return b;}
 Field(int P):p(P),pw(M),inv(P){
  for(int d=2;d*d<=p;d++) assert(p%d);
  assert((p-1)%M==0); vector<int> fac;int n=p-1;
  for(int d=2;d*d<=n;d++)if(n%d==0){fac.push_back(d);while(n%d==0)n/=d;}
  if(n>1)fac.push_back(n);
  for(g=2;;g++){bool ok=true;for(int d:fac)if(power(g,(p-1)/d)==1)ok=false;if(ok)break;}
  z=power(g,(p-1)/M);assert(power(z,M)==1);
  for(int d:{2,3,5,7})assert(power(z,M/d)!=1);
  pw[0]=1;for(int j=1;j<M;j++)pw[j]=mul(pw[j-1],z);
  inv[1]=1;for(int j=2;j<p;j++)inv[j]=p-int(1LL*(p/j)*inv[p%j]%p);
 }
 int rank(Mat a,int k) const {
  int r=0;for(int c=0;c<k&&r<4;c++){
   int s=r;while(s<4&&!a[s][c])s++;if(s==4)continue;swap(a[s],a[r]);
   int v=inv[a[r][c]];for(int j=c;j<k;j++)a[r][j]=mul(a[r][j],v);
   for(int i=r+1;i<4;i++)if(a[i][c]){v=a[i][c];for(int j=c;j<k;j++)a[i][j]=(a[i][j]-mul(v,a[r][j])+p)%p;}r++;
  }return r;
 }
 int det4(Mat a,int skip) const {
  int b[4][4];for(int i=0;i<4;i++){int t=0;for(int j=0;j<5;j++)if(j!=skip)b[i][t++]=a[i][j];}
  int ans=1;for(int c=0;c<4;c++){
   int s=c;while(s<4&&!b[s][c])s++;if(s==4)return 0;
   if(s!=c){for(int j=c;j<4;j++)swap(b[s][j],b[c][j]);ans=(p-ans)%p;}
   int pivot=b[c][c];ans=mul(ans,pivot);int v=inv[pivot];
   for(int i=c+1;i<4;i++){int w=mul(b[i][c],v);for(int j=c+1;j<4;j++)b[i][j]=(b[i][j]-mul(w,b[c][j])+p)%p;}
  }return ans;
 }
 int minor_det(const Mat &a,int k,int rowmask,int skip) const {
  int b[4][4],r=0;for(int i=0;i<4;i++)if((rowmask>>i)&1){int t=0;for(int j=0;j<k;j++)if(j!=skip)b[r][t++]=a[i][j];r++;}
  assert(r==k-1);int ans=1;for(int c=0;c<r;c++){
   int s=c;while(s<r&&!b[s][c])s++;if(s==r)return 0;
   if(s!=c){for(int j=c;j<r;j++)swap(b[s][j],b[c][j]);ans=(p-ans)%p;}
   int pivot=b[c][c];ans=mul(ans,pivot);int v=inv[pivot];
   for(int i=c+1;i<r;i++){int w=mul(b[i][c],v);for(int j=c+1;j<r;j++)b[i][j]=(b[i][j]-mul(w,b[c][j])+p)%p;}
  }return ans;
 }
 // Separate determinant check: explicit Leibniz expansion, no elimination.
 int expanded_minor(const Mat &a,int rowmask,int colmask) const {
  vector<int> rows,cols;for(int i=0;i<4;i++)if((rowmask>>i)&1)rows.push_back(i);
  for(int j=0;j<5;j++)if((colmask>>j)&1)cols.push_back(j);
  assert(rows.size()==cols.size());
#ifdef AUDIT
  int n=rows.size(),b[4][4]{};for(int i=0;i<n;i++)for(int j=0;j<n;j++)b[i][j]=a[rows[i]][cols[j]];
  int det=1;for(int c=0;c<n;c++){
   int r=c;while(r<n&&!b[r][c])r++;if(r==n)return 0;
   if(r!=c){for(int j=0;j<n;j++)swap(b[r][j],b[c][j]);det=(p-det)%p;}
   det=mul(det,b[c][c]);int iv=inv[b[c][c]];
   for(int i=c+1;i<n;i++){int factor=mul(b[i][c],iv);for(int j=c;j<n;j++)b[i][j]=(b[i][j]-mul(factor,b[c][j])+p)%p;}
  }return det;
#endif
  int ans=0;vector<int> s(cols.size());iota(s.begin(),s.end(),0);
  do{int v=1,parity=0;for(int i=0;i<(int)s.size();i++){
    v=mul(v,a[rows[i]][cols[s[i]]]);for(int j=i+1;j<(int)s.size();j++)parity^=(s[i]>s[j]);
   }ans=(ans+(parity?p-v:v))%p;
  }while(next_permutation(s.begin(),s.end()));return ans;
 }
};
vector<Field> fields;
vector<int> sizes,block(11),exponent(11);vector<vector<int>> blocks;
long long partitions=0,tested=0,rank_excluded=0,norm_excluded=0,unresolved=0;
array<long long,3> byfield{};ofstream pending;
char profile;
// Bihomogeneous necessary conditions on P^1 x P^1.
struct Span16 {
 const Field &f;array<array<int,16>,16> basis{};array<bool,16> used{};int rank=0;
 Span16(const Field &F):f(F){}
 void add(array<int,16> row){
#ifdef AUDIT
  reverse(row.begin(),row.end());
#endif
  for(int c=0;c<16;c++)if(row[c]){
   if(!used[c]){int inv=f.inv[row[c]];for(int j=c;j<16;j++)row[j]=f.mul(row[j],inv);basis[c]=row;used[c]=true;rank++;return;}
   int v=row[c];for(int j=c;j<16;j++)row[j]=(row[j]-f.mul(v,basis[c][j])+f.p)%f.p;
  }
 }
 void polynomial(const array<int,9>&p,int dq,int ds){
  for(int a=0;a<=3-dq;a++)for(int b=0;b<=3-ds;b++){
   array<int,16> row{};
   for(int i=0;i<=dq;i++)for(int j=0;j<=ds;j++)row[4*(a+i)+b+j]=p[3*i+j];
   add(row);if(rank==16)return;
  }
 }
};
using UP=vector<int>;
UP trim(UP a){while(!a.empty()&&!a.back())a.pop_back();return a;}
UP addp(const Field&f,UP a,const UP&b,int sign=1){a.resize(max(a.size(),b.size()));for(int i=0;i<(int)b.size();i++)a[i]=(a[i]+(sign==1?b[i]:f.p-b[i]))%f.p;return trim(a);}
UP mulp(const Field&f,const UP&a,const UP&b){if(a.empty()||b.empty())return {};UP c(a.size()+b.size()-1);for(int i=0;i<(int)a.size();i++)for(int j=0;j<(int)b.size();j++)c[i+j]=(c[i+j]+f.mul(a[i],b[j]))%f.p;return trim(c);}
UP remp(const Field&f,UP a,const UP&b){assert(!b.empty());int iv=f.inv[b.back()];while(a.size()>=b.size()){
 int d=a.size()-b.size(),v=f.mul(a.back(),iv);for(int j=0;j<(int)b.size();j++)a[d+j]=(a[d+j]-f.mul(v,b[j])+f.p)%f.p;a=trim(a);
 }return a;}
UP gcdp(const Field&f,UP a,UP b){while(!b.empty()){UP c=remp(f,a,b);a=b;b=c;}if(!a.empty()){int v=f.inv[a.back()];for(int &x:a)x=f.mul(x,v);}return a;}
UP powx(const Field&f,int n,const UP&g){UP a=remp(f,{0,1},g),out=remp(f,{1},g);for(;n;n>>=1,a=remp(f,mulp(f,a,a),g))if(n&1)out=remp(f,mulp(f,out,a),g);return out;}
struct BP{array<int,9> p{};int dq=0,ds=0;};
UP resultant(const Field&f,const BP&a,const BP&b,bool transpose){
 int da=transpose?a.dq:a.ds,db=transpose?b.dq:b.ds;
 auto coef=[&](const BP &x,int j){UP c(3);for(int i=0;i<3;i++)c[i]=x.p[transpose?3*j+i:3*i+j];return trim(c);};
 if(!da)return coef(a,0);if(!db)return coef(b,0);
#ifdef AUDIT
 int n=da+db;vector<vector<UP>>S(n,vector<UP>(n));
 for(int i=0;i<db;i++)for(int j=0;j<=da;j++)S[i][i+j]=coef(a,da-j);
 for(int i=0;i<da;i++)for(int j=0;j<=db;j++)S[db+i][i+j]=coef(b,db-j);
 vector<int>perm(n);iota(perm.begin(),perm.end(),0);UP out;
 do{UP term={1};int sign=1;for(int i=0;i<n;i++){term=mulp(f,term,S[i][perm[i]]);for(int j=i+1;j<n;j++)if(perm[i]>perm[j])sign=-sign;}out=addp(f,out,term,sign);}while(next_permutation(perm.begin(),perm.end()));return out;
#endif
 UP a0=coef(a,0),a1=coef(a,1),b0=coef(b,0),b1=coef(b,1);
 if(da==1&&db==1)return addp(f,mulp(f,a1,b0),mulp(f,a0,b1),-1);
 if(da==2&&db==1)return resultant(f,b,a,transpose);
 UP b2=coef(b,2);
 if(da==1)return addp(f,addp(f,mulp(f,b0,mulp(f,a1,a1)),mulp(f,b1,mulp(f,a0,a1)),-1),mulp(f,b2,mulp(f,a0,a0)));
 UP a2=coef(a,2),d=addp(f,mulp(f,a2,b0),mulp(f,a0,b2),-1);
 UP e=addp(f,mulp(f,a2,b1),mulp(f,a1,b2),-1),h=addp(f,mulp(f,a1,b0),mulp(f,a0,b1),-1);
 return addp(f,mulp(f,d,d),mulp(f,e,h),-1);
}
bool torsion_excludes(const Field&f,const vector<BP>&polys,bool transpose){
 UP g;
 for(int i=0;i<(int)polys.size();i++){
  const BP&a=polys[i];int ds=transpose?a.dq:a.ds;
  if(!ds){UP u(3);for(int n=0;n<3;n++)u[n]=a.p[transpose?n:3*n];u=trim(u);if(!u.empty())g=gcdp(f,g,u);}
  for(int j=0;j<i;j++){UP r=resultant(f,polys[j],a,transpose);if(!r.empty())g=gcdp(f,g,r);}
  if(g.size()==1)return true;
 }
 if(g.empty())return false;
 assert(g.size()<=9);
 UP tors=addp(f,powx(f,352800,g),{1},-1);UP d=gcdp(f,g,tors);
 return d.size()==1 || remp(f,{f.p-1,0,0,0,1},d).empty();
}
