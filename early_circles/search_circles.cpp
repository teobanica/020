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
  for(int j=1;j<p;j++)inv[j]=power(j,p-2);
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
  assert(rows.size()==cols.size());int ans=0;vector<int> s(cols.size());iota(s.begin(),s.end(),0);
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
// Necessary polynomials in q=conjugate(t_4), with the other t_r fixed.
// Rank four of their degree-three homogeneous multiples excludes every q in P^1.
struct Span {
 const Field &f; array<array<int,4>,4> basis{}; array<bool,4> used{};Mat witnesses{};int rank=0;
 Span(const Field &field):f(field){}
 void add(array<int,4> row){
  const auto original=row;
  for(int c=0;c<4;c++)if(row[c]){
   if(!used[c]){int inv=f.inv[row[c]];for(int j=c;j<4;j++)row[j]=f.mul(row[j],inv);basis[c]=row;used[c]=true;for(int j=0;j<4;j++)witnesses[rank][j]=original[j];rank++;return;}
   int x=row[c];for(int j=c;j<4;j++)row[j]=(row[j]-f.mul(x,basis[c][j])+f.p)%f.p;
  }
 }
 void linear(int a,int b){add({a,b,0,0});add({0,a,b,0});add({0,0,a,b});}
 void quadratic(int a,int b,int c){add({a,b,c,0});add({0,a,b,c});}
};
void test(){
 tested++;int k=sizes.size();
 for(int fi=0;fi<(int)fields.size();fi++){
  const auto &f=fields[fi];Mat a0{},a1{},b0{},b1{};
  for(int r=0;r<4;r++)for(int j=0;j<11;j++){
   int e=(exponent[j]-21*H[r+1][j]+M)%M;
   bool variable=(r==3&&j>=5);
   auto &a=variable?a1:a0;auto &b=variable?b1:b0;
   a[r][block[j]]=(a[r][block[j]]+f.pw[e])%f.p;
   b[r][block[j]]=(b[r][block[j]]+f.pw[(M-e)%M])%f.p;
  }
  for(int r=0;r<3;r++)for(int c=0;c<k;c++){a1[r][c]=a0[r][c];b1[r][c]=b0[r][c];}
  Span span(f);bool exclude=false;
  if(k<=4)for(int mask=0;mask<16;mask++)if(__builtin_popcount((unsigned)mask)==k){
   int d0=f.expanded_minor(a0,mask,(1<<k)-1);
   if(!(mask&8)){if(d0){exclude=true;break;}}
   else {int d1=f.expanded_minor(a1,mask,(1<<k)-1);span.linear(d0,d1);}
  }
  if(!exclude)for(int mask=0;mask<16;mask++)if(__builtin_popcount((unsigned)mask)==k-1){
   array<array<int,3>,5> pol{};
   for(int j=0;j<k;j++){
    int cm=((1<<k)-1)^(1<<j);
    int u0=f.expanded_minor(a0,mask,cm),v0=f.expanded_minor(b0,mask,cm);
    if(!(mask&8)){pol[j]={0,f.mul(u0,v0),0};}
    else {
     int u1=f.expanded_minor(a1,mask,cm),v1=f.expanded_minor(b1,mask,cm);
     pol[j]={f.mul(u0,v1),(f.mul(u0,v0)+f.mul(u1,v1))%f.p,f.mul(u1,v0)};
    }
   }
   for(int j=1;j<k;j++){
    array<int,3> d{};for(int s=0;s<3;s++)d[s]=(pol[j][s]-pol[0][s]+f.p)%f.p;
    if(!(mask&8)){if(d[1]){exclude=true;break;}}
    else span.quadratic(d[0],d[1],d[2]);
   }
   if(exclude||span.rank==4)break;
  }
  if(exclude||span.rank==4){if(span.rank==4)assert(f.expanded_minor(span.witnesses,15,15));rank_excluded++;byfield[fi]++;return;}
 }
 unresolved++;pending<<profile;
 for(int x:block)pending<<' '<<x;for(int x:exponent)pending<<' '<<x;pending<<'\n';
}

void labels(int k){
 if(k==(int)blocks.size()){test();return;}
 const auto &s=blocks[k];vector<int> perm(s.size()-1);iota(perm.begin(),perm.end(),1);
 do{exponent[s[0]]=0;for(int j=1;j<(int)s.size();j++)exponent[s[j]]=perm[j-1]*M/s.size();labels(k+1);}while(next_permutation(perm.begin(),perm.end()));
}
void partition(int k,vector<int> rem){
 if(k==(int)sizes.size()){assert(rem.empty());partitions++;labels(0);return;}
 int n=rem.size(),s=sizes[k];
 for(int mask=0;mask<(1<<n);mask++)if(__builtin_popcount((unsigned)mask)==s){
  vector<int> chosen,rest;for(int j=0;j<n;j++)((mask>>j)&1?chosen:rest).push_back(rem[j]);
  if(k&&sizes[k]==sizes[k-1]&&chosen[0]<blocks.back()[0])continue;
  for(int j:chosen)block[j]=k;blocks.push_back(chosen);partition(k+1,rest);blocks.pop_back();
 }
}
int main(int argc,char **argv){
 assert(argc==6);int variable=stoi(argv[1]);assert(variable>=0&&variable<4);
 int temp[5][11];for(int r=0;r<5;r++)for(int j=0;j<11;j++)temp[r][j]=H[r][j];
 int row=1;for(int r=1;r<5;r++)if(r!=variable+1){
  int shift=stoi(argv[r+1]);assert(shift>=0&&shift<20);
  for(int j=0;j<11;j++)H[row][j]=(temp[r][j]+(j>=5?shift:0))%20;row++;
 }
 for(int j=0;j<11;j++)H[4][j]=temp[variable+1][j];

 fields.emplace_back(421);fields.emplace_back(3361);fields.emplace_back(4621);
 for(auto &f:fields)cerr<<"field "<<f.p<<" generator "<<f.g<<" root420 "<<f.z<<'\n';
 pending.open("unresolved.txt");
 vector<vector<int>> profiles={{7,2,2},{5,3,3},{5,2,2,2},{3,3,3,2},{3,2,2,2,2}};
 const long long expected_partitions[]={990,4620,6930,15400,17325};
 const long long expected_labels[]={712800,443520,166320,123200,34650};
 cout<<"profile partitions labelled_presentations rank_excluded norm_excluded unresolved field1 field2 field3\n";
 for(int i=0;i<5;i++){
  profile='B'+i;sizes=profiles[i];partitions=tested=rank_excluded=norm_excluded=unresolved=0;byfield.fill(0);
  vector<int> rem(11);iota(rem.begin(),rem.end(),0);partition(0,rem);
  assert(partitions==expected_partitions[i]&&tested==expected_labels[i]);
  cout<<profile<<' '<<partitions<<' '<<tested<<' '<<rank_excluded<<' '<<norm_excluded<<' '<<unresolved;
  for(auto n:byfield)cout<<' '<<n;cout<<endl;
 }
}
