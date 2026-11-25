#include "common.hpp"
struct Context {vector<BP> base;array<BP,5> cof;};
vector<Context> contexts;
Context build_context(const Field&f){
 Context ctx;array<Mat,4>A{},B{};
 for(int m=0;m<4;m++)for(int r=0;r<4;r++)for(int j=0;j<11;j++){
  if(r==2&&((m&1)!=0)!=(j>=5))continue;
  if(r==3&&((m&2)!=0)!=(j>=5))continue;
  int e=(exponent[j]-21*H[r+1][j]+M)%M;
  A[m][r][block[j]]=(A[m][r][block[j]]+f.pw[e])%f.p;
  B[m][r][block[j]]=(B[m][r][block[j]]+f.pw[(M-e)%M])%f.p;
 }
 array<array<int,9>,5>N{};
 for(int c=0;c<5;c++){
  int mask=31^(1<<c);array<int,4>a{},b{};ctx.cof[c].dq=ctx.cof[c].ds=1;
  for(int i=0;i<2;i++)for(int j=0;j<2;j++){
   a[i+2*j]=f.expanded_minor(A[i+2*j],15,mask);b[i+2*j]=f.expanded_minor(B[i+2*j],15,mask);
   ctx.cof[c].p[3*i+j]=(c%2?f.p-a[i+2*j]:a[i+2*j])%f.p;
  }
  for(int i=0;i<2;i++)for(int j=0;j<2;j++)for(int u=0;u<2;u++)for(int v=0;v<2;v++){
   int e=3*(i-u+1)+(j-v+1);N[c][e]=(N[c][e]+f.mul(a[i+2*j],b[u+2*v]))%f.p;
  }
 }
 // Retain a subset of original necessary polynomials spanning their reduced coefficient space.
 array<array<int,9>,9>basis{};array<bool,9>used{};
 for(int c=1;c<5;c++){
  BP pol;pol.dq=pol.ds=2;for(int e=0;e<9;e++)pol.p[e]=(N[c][e]-N[0][e]+f.p)%f.p;
  auto row=pol.p;bool take=false;
  for(int e=0;e<9;e++)if(row[e]){
   if(!used[e]){int iv=f.inv[row[e]];for(int n=e;n<9;n++)row[n]=f.mul(row[n],iv);basis[e]=row;used[e]=true;take=true;break;}
   int z=row[e];for(int n=e;n<9;n++)row[n]=(row[n]-f.mul(z,basis[e][n])+f.p)%f.p;
  }
  if(take)ctx.base.push_back(pol);
 }
 return ctx;
}

vector<int> isolate(const Field&f,const vector<BP>&polys,bool t){
 UP g;for(int i=0;i<(int)polys.size();i++)for(int j=0;j<i;j++){auto r=resultant(f,polys[i],polys[j],t);if(!r.empty())g=gcdp(f,g,r);}
 if(g.empty())return {-1};assert(g.size()<=9);
 UP d=gcdp(f,g,addp(f,powx(f,352800,g),{1},-1));vector<int>roots;UP product={1};
 for(int e=0;e<M;e++)if(remp(f,d,{(f.p-f.pw[e])%f.p,1}).empty()){roots.push_back(e);product=mulp(f,product,{(f.p-f.pw[e])%f.p,1});}
 assert(remp(f,product,d).empty());return roots;
}
int main(int argc,char**argv){
 assert(argc==7);int v0=stoi(argv[1]),v1=stoi(argv[2]);
 int temp[5][11];for(int r=0;r<5;r++)for(int j=0;j<11;j++)temp[r][j]=H[r][j];
 int row=1;for(int r=1;r<5;r++)if(r!=v0+1&&r!=v1+1){int shift=stoi(argv[r+2]);for(int j=0;j<11;j++)H[row][j]=(temp[r][j]+(j>=5?shift:0))%20;row++;}
 for(int j=0;j<11;j++){H[3][j]=temp[v0+1][j];H[4][j]=temp[v1+1][j];}
 Field f(43261);ifstream input("unresolved.txt"),retained("rank2_unresolved.txt");vector<int>ids;int id;while(retained>>id)ids.push_back(id);
 ofstream curves("curve_unresolved.txt");int idx=0;while(input>>profile){for(int&v:block)input>>v;for(int&v:exponent)input>>v;
  if(find(ids.begin(),ids.end(),idx)!=ids.end()){
   auto c=build_context(f);auto q=isolate(f,c.base,false),s=isolate(f,c.base,true);assert(q.size()==1&&s.size()==1);
   if(q[0]<0||s[0]<0)curves<<idx<<'\n';else cout<<idx<<' '<<q[0]<<' '<<s[0]<<'\n';
  }idx++;
 }
}
