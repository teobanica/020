#include "common.hpp"
vector<pair<int,int>> dependencies;
int method=0;
bool run_case(const Field&f){
 dependencies.clear();
 array<Mat,4>A{},B{};
 for(int m=0;m<4;m++)for(int r=0;r<4;r++)for(int j=0;j<11;j++){
  if(r==2&&((m&1)!=0)!=(j>=5))continue;
  if(r==3&&((m&2)!=0)!=(j>=5))continue;
  int e=(exponent[j]-21*H[r+1][j]+M)%M;
  A[m][r][block[j]]=(A[m][r][block[j]]+f.pw[e])%f.p;
  B[m][r][block[j]]=(B[m][r][block[j]]+f.pw[(M-e)%M])%f.p;
 }
 vector<BP>polys;
 for(int a=0;a<5;a++)for(int b=a+1;b<5;b++){
  int cm=31^(1<<a)^(1<<b);bool zero=true;
  for(int rows=0;rows<16;rows++)if(__builtin_popcount((unsigned)rows)==3)
   for(int i=0;i<4;i++)if(f.expanded_minor(A[i],rows,cm))zero=false;
  if(!zero)continue;
  dependencies.push_back({a,b});
  // Exact characteristic-zero certification is supplied separately.
  for(int omit=0;omit<5;omit++)if((cm>>omit)&1){
   int cc=cm^(1<<omit);
   for(int rows=0;rows<16;rows++)if(__builtin_popcount((unsigned)rows)==3){
    int dq=(rows>>2)&1,ds=(rows>>3)&1;
    array<array<int,9>,2>N{};
    for(int v=0;v<2;v++){
     int cols=cc|(1<<(v?b:a));array<int,4>u{},w{};
     for(int i=0;i<=dq;i++)for(int j=0;j<=ds;j++){
      u[i+2*j]=f.expanded_minor(A[i+2*j],rows,cols);
      w[i+2*j]=f.expanded_minor(B[i+2*j],rows,cols);
     }
     for(int i=0;i<=dq;i++)for(int j=0;j<=ds;j++)for(int h=0;h<=dq;h++)for(int l=0;l<=ds;l++){
      int e=3*(i-h+dq)+(j-l+ds);N[v][e]=(N[v][e]+f.mul(u[i+2*j],w[h+2*l]))%f.p;
     }
    }
    BP pol;pol.dq=2*dq;pol.ds=2*ds;for(int e=0;e<9;e++)pol.p[e]=(N[0][e]-N[1][e]+f.p)%f.p;
    polys.push_back(pol);
   }
  }
 }
 Span16 span(f);for(auto&pol:polys)span.polynomial(pol.p,pol.dq,pol.ds);
 if(span.rank==16){method=1;return true;}
 if(torsion_excludes(f,polys,false)){method=2;return true;}
 if(torsion_excludes(f,polys,true)){method=3;return true;}
 return false;
}
int main(int argc,char**argv){
 assert(argc==7);int v0=stoi(argv[1]),v1=stoi(argv[2]);assert(0<=v0&&v0<v1&&v1<4);
 int temp[5][11];for(int r=0;r<5;r++)for(int j=0;j<11;j++)temp[r][j]=H[r][j];
 int row=1;for(int r=1;r<5;r++)if(r!=v0+1&&r!=v1+1){int shift=stoi(argv[r+2]);for(int j=0;j<11;j++)H[row][j]=(temp[r][j]+(j>=5?shift:0))%20;row++;}
 for(int j=0;j<11;j++){H[3][j]=temp[v0+1][j];H[4][j]=temp[v1+1][j];}
 fields.emplace_back(43261);fields.emplace_back(55021);fields.emplace_back(62581);
 ifstream input("unresolved.txt");ofstream output("rank2_unresolved.txt"),depout("dependencies.txt");int index=0,passed=0;
 while(input>>profile){assert(profile=='F');for(int&v:block)input>>v;for(int&v:exponent)input>>v;
  bool ok=false;for(const auto&f:fields)if(run_case(f)){ok=true;break;}
  if(ok)passed++;else output<<index<<'\n';
  if(ok){depout<<index;for(auto [a,b]:dependencies)depout<<' '<<a<<' '<<b;depout<<'\n';}
  cout<<index++<<' '<<ok<<' '<<(ok?method:0)<<endl;
 }
 cout<<"TOTAL "<<index<<" excluded "<<passed<<endl;
}
