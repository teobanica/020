#include "common.hpp"
long long macaulay=0,torsion_q=0,torsion_s=0;
void test(){
 tested++;int k=sizes.size();
 for(int fi=0;fi<(int)fields.size();fi++){
  const Field &f=fields[fi];array<Mat,4>A{},B{};
  for(int m=0;m<4;m++)for(int r=0;r<4;r++)for(int j=0;j<11;j++){
   if(r==2 && ((m&1)!=0)!=(j>=5))continue;
   if(r==3 && ((m&2)!=0)!=(j>=5))continue;
   int e=(exponent[j]-21*H[r+1][j]+M)%M;
   A[m][r][block[j]]=(A[m][r][block[j]]+f.pw[e])%f.p;
   B[m][r][block[j]]=(B[m][r][block[j]]+f.pw[(M-e)%M])%f.p;
  }
  Span16 span(f);vector<BP> polys;
  if(k<=4)for(int mask=0;mask<16;mask++)if(__builtin_popcount((unsigned)mask)==k){
   BP pol;pol.dq=(mask>>2)&1;pol.ds=(mask>>3)&1;
   for(int i=0;i<=pol.dq;i++)for(int j=0;j<=pol.ds;j++)pol.p[3*i+j]=f.expanded_minor(A[i+2*j],mask,(1<<k)-1);
   polys.push_back(pol);span.polynomial(pol.p,pol.dq,pol.ds);if(span.rank==16)break;
  }
  if(span.rank!=16)for(int mask=0;mask<16;mask++)if(__builtin_popcount((unsigned)mask)==k-1){
   int dq=(mask>>2)&1,ds=(mask>>3)&1;array<array<int,9>,5>P{};
   for(int omit=0;omit<k;omit++){
    array<int,4>a{},b{};int cm=((1<<k)-1)^(1<<omit);
    for(int i=0;i<=dq;i++)for(int j=0;j<=ds;j++){
     a[i+2*j]=f.expanded_minor(A[i+2*j],mask,cm);
     b[i+2*j]=f.expanded_minor(B[i+2*j],mask,cm);
    }
    for(int i=0;i<=dq;i++)for(int j=0;j<=ds;j++)for(int u=0;u<=dq;u++)for(int v=0;v<=ds;v++){
     int pos=3*(i-u+dq)+(j-v+ds);
     P[omit][pos]=(P[omit][pos]+f.mul(a[i+2*j],b[u+2*v]))%f.p;
    }
   }
   for(int omit=1;omit<k;omit++){
    BP pol;pol.dq=2*dq;pol.ds=2*ds;for(int n=0;n<9;n++)pol.p[n]=(P[omit][n]-P[0][n]+f.p)%f.p;
    polys.push_back(pol);if(k<5)span.polynomial(pol.p,pol.dq,pol.ds);if(span.rank==16)break;
   }
   if(span.rank==16)break;
  }
  if(span.rank==16){rank_excluded++;macaulay++;byfield[fi]++;return;}
  if(torsion_excludes(f,polys,false)){rank_excluded++;torsion_q++;byfield[fi]++;return;}
  if(torsion_excludes(f,polys,true)){rank_excluded++;torsion_s++;byfield[fi]++;return;}
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
 assert(argc==7);int v0=stoi(argv[1]),v1=stoi(argv[2]);assert(0<=v0&&v0<v1&&v1<4);
 int temp[5][11];for(int r=0;r<5;r++)for(int j=0;j<11;j++)temp[r][j]=H[r][j];
 int row=1;for(int r=1;r<5;r++)if(r!=v0+1&&r!=v1+1){
  int shift=stoi(argv[r+2]);assert(shift>=0&&shift<20);
  for(int j=0;j<11;j++)H[row][j]=(temp[r][j]+(j>=5?shift:0))%20;row++;
 }
 for(int j=0;j<11;j++){H[3][j]=temp[v0+1][j];H[4][j]=temp[v1+1][j];}

 fields.emplace_back(43261);fields.emplace_back(55021);fields.emplace_back(62581);
 for(auto &f:fields)cerr<<"field "<<f.p<<" generator "<<f.g<<" root420 "<<f.z<<'\n';
 pending.open("unresolved.txt");
 vector<vector<int>> profiles={{7,2,2},{5,3,3},{5,2,2,2},{3,3,3,2},{3,2,2,2,2}};
 const long long expected_partitions[]={990,4620,6930,15400,17325};
 const long long expected_labels[]={712800,443520,166320,123200,34650};
 cout<<"profile partitions labelled_presentations excluded reserved_zero unresolved field1 field2 field3 macaulay torsion_q torsion_s\n";
 for(int i=0;i<5;i++){
  profile='B'+i;sizes=profiles[i];partitions=tested=rank_excluded=norm_excluded=unresolved=0;byfield.fill(0);macaulay=torsion_q=torsion_s=0;
  vector<int> rem(11);iota(rem.begin(),rem.end(),0);partition(0,rem);
  assert(partitions==expected_partitions[i]&&tested==expected_labels[i]);
  cout<<profile<<' '<<partitions<<' '<<tested<<' '<<rank_excluded<<' '<<norm_excluded<<' '<<unresolved;
  for(auto n:byfield)cout<<' '<<n;cout<<' '<<macaulay<<' '<<torsion_q<<' '<<torsion_s<<endl;
 }
}
