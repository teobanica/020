#include "common.hpp"
long long macaulay=0,torsion_q=0,torsion_s=0;
bool test(){
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
  if(span.rank==16){rank_excluded++;macaulay++;byfield[fi]++;return true;}
  if(torsion_excludes(f,polys,false)){rank_excluded++;torsion_q++;byfield[fi]++;return true;}
  if(torsion_excludes(f,polys,true)){rank_excluded++;torsion_s++;byfield[fi]++;return true;}
 }
 return false;
}

vector<int> oldblock,oldexp;int target;long long branches=0,failed=0;vector<int>potential;vector<bool>assigned;vector<bool>used;
bool cycle(int prime,int label){
 if(label==prime){
  vector<int>mapped(5,-1);int k=1;for(int b=0;b<5;b++)if(assigned[b])mapped[b]=0;else mapped[b]=k++;
  assert(k<=4);block=oldblock;exponent=oldexp;
  for(int j=0;j<11;j++){block[j]=mapped[oldblock[j]];exponent[j]=(exponent[j]+(assigned[oldblock[j]]?potential[oldblock[j]]:0))%M;}
  sizes.assign(k,1);branches++;bool ok=test();if(!ok)failed++;
  return ok;
 }
 int et=(oldexp[target]-21*H[1][target]+M)%M;
 for(int col=0;col<11;col++)if(!used[col]){
  int b=oldblock[col],ec=(oldexp[col]-21*H[1][col]+M)%M;
  int d=(label*(M/prime)+et-ec+M)%M;
  bool was=assigned[b];if(was&&potential[b]!=d)continue;
  assigned[b]=true;potential[b]=d;used[col]=true;
  bool ok=cycle(prime,label+1);used[col]=false;assigned[b]=was;
  if(!ok)return false;
 }
 return true;
}
int main(int argc,char**argv){
 assert(argc==7);int v0=stoi(argv[1]),v1=stoi(argv[2]);
 int temp[5][11];for(int r=0;r<5;r++)for(int j=0;j<11;j++)temp[r][j]=H[r][j];
 int row=1;for(int r=1;r<5;r++)if(r!=v0+1&&r!=v1+1){int shift=stoi(argv[r+2]);for(int j=0;j<11;j++)H[row][j]=(temp[r][j]+(j>=5?shift:0))%20;row++;}
 for(int j=0;j<11;j++){H[3][j]=temp[v0+1][j];H[4][j]=temp[v1+1][j];}
 fields.emplace_back(43261);fields.emplace_back(55021);fields.emplace_back(62581);
 ifstream input("unresolved.txt"),retained("curve_unresolved.txt");vector<int>ids;int id;while(retained>>id)ids.push_back(id);
 int idx=0;while(input>>profile){for(int&v:block)input>>v;for(int&v:exponent)input>>v;
  if(find(ids.begin(),ids.end(),idx)!=ids.end()){
   oldblock=block;oldexp=exponent;target=find(block.begin(),block.end(),0)-block.begin();branches=failed=0;
   vector<long long>counts;bool ok=true;for(int prime:{2,3,5,7}){
    long long before=branches;
    assigned.assign(5,false);potential.assign(5,0);used.assign(11,false);assigned[0]=true;used[target]=true;
    if(!cycle(prime,1)){ok=false;break;}counts.push_back(branches-before);
   }
   cout<<idx<<' '<<ok<<' '<<branches<<' '<<failed;for(auto n:counts)cout<<' '<<n;cout<<endl;
  }idx++;
 }
}
