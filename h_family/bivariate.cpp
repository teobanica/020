#include "common.hpp"
#include <sstream>
struct BP{array<int,9>p{};int dq=2,ds=2;};
Poly resultant(const BP&a,const BP&b,bool flip){
 auto coefficient=[&](const BP&x,int j){Poly z(3);for(int i=0;i<3;i++)z[i]=x.p[flip?3*j+i:3*i+j];return trim(z);};
#ifdef AUDIT
 vector<vector<Poly>>S(4,vector<Poly>(4));for(int i=0;i<2;i++)for(int j=0;j<=2;j++){S[i][i+j]=coefficient(a,2-j);S[i+2][i+j]=coefficient(b,2-j);}
 vector<int>permutation={0,1,2,3};Poly out;
 do{Poly term={1};int sign=1;for(int i=0;i<4;i++){term=mul(term,S[i][permutation[i]]);for(int j=i+1;j<4;j++)if(permutation[i]>permutation[j])sign=-sign;}out=add(out,term,sign);}while(next_permutation(permutation.begin(),permutation.end()));return out;
#endif
 auto a0=coefficient(a,0),a1=coefficient(a,1),a2=coefficient(a,2),b0=coefficient(b,0),b1=coefficient(b,1),b2=coefficient(b,2);
 auto d=add(mul(a2,b0),mul(a0,b2),-1),e=add(mul(a2,b1),mul(a1,b2),-1),h=add(mul(a1,b0),mul(a0,b1),-1);
 return add(mul(d,d),mul(e,h),-1);
}
struct Span{array<array<int,16>,16>b{};array<bool,16>used{};int rank=0;
 void addrow(array<int,16>v){
#ifdef AUDIT
 reverse(v.begin(),v.end());
#endif
for(int c=0;c<16;c++)if(v[c]){if(!used[c]){int iv=f.iv[v[c]];for(int j=c;j<16;j++)v[j]=f.mul(v[j],iv);b[c]=v;used[c]=true;rank++;return;}int w=v[c];for(int j=c;j<16;j++)v[j]=mod(v[j]-f.mul(w,b[c][j]),P);}}
 void polynomial(const BP&p){for(int a=0;a<2;a++)for(int b=0;b<2;b++){array<int,16>v{};for(int i=0;i<3;i++)for(int j=0;j<3;j++)v[4*(a+i)+b+j]=p.p[3*i+j];addrow(v);}}
};
Poly projection(const vector<BP>&ps,bool flip){Poly g;for(int i=0;i<(int)ps.size();i++)for(int j=0;j<i;j++){auto r=resultant(ps[i],ps[j],flip);if(!r.empty())g=gcdp(g,r);}return g;}
int main(int argc,char**argv){
 assert(argc==2);string path=argv[1];ifstream in(path+"/refine_unresolved.txt");ofstream left(path+"/bivariate_unresolved.txt");points.open(path+"/bivariate_points.txt");string line;int count=0,spanexcluded=0,torsionexcluded=0,isolated=0,unresolved=0;
 while(getline(in,line)){
  istringstream row(line);string kind;int active,c;row>>tests>>kind>>active>>c;assert(kind=="GENUINE_Q_FREE"&&c==2);
  vector<int>a(11),b(11),e(11);for(int j=0;j<11;j++)row>>b[j]>>a[j]>>e[j];assert(row);
  assert(*max_element(a.begin(),a.end())-*min_element(a.begin(),a.end())==1);
  vector<BP>ps(4);for(int r=1;r<5;r++)for(int j=0;j<11;j++)for(int k=0;k<11;k++)if((j<5)==(k<5)){
   int pos=3*(a[j]-a[k]+1)+b[j]-b[k]+1;
   int exp=mod(e[j]-e[k]-(L/20)*(H[r][j]-H[r][k]),L);
   ps[r-1].p[pos]=mod(ps[r-1].p[pos]+(j<5?1:-1)*f.pw[exp],P);
  }
  count++;Span span;for(auto p:ps)span.polynomial(p);if(span.rank==16){spanexcluded++;continue;}
  auto g=projection(ps,false),h=projection(ps,true);
  if(!g.empty())g=gcdp(g,add(powerp({0,1},L,g),{1},-1));
  if(!h.empty())h=gcdp(h,add(powerp({0,1},L,h),{1},-1));
  if(g.size()==1||h.size()==1){torsionexcluded++;continue;}
  if(g.empty()||h.empty()){unresolved++;left<<line<<'\n';continue;}
  isolated++;
  for(int qr:roots(g))for(int sr:roots(h)){
   int qe=f.log[qr],se=f.log[sr];assert(qe>=0&&se>=0);vector<int>x(11);
   for(int j=0;j<11;j++)x[j]=mod(a[j]*qe+b[j]*se+e[j],L);
   checkpoint_x(x,qe);
  }
 }
 cout<<"count "<<count<<" span "<<spanexcluded<<" torsion "<<torsionexcluded<<" isolated "<<isolated<<" unresolved "<<unresolved<<" rootpairs "<<rootcandidates<<" normpassed "<<rootpassed<<" actual "<<actual<<endl;
}
