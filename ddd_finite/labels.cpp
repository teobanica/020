#include <algorithm>
#include <array>
#include <fstream>
#include <iostream>
#include <map>
#include <numeric>
#include <set>
#include <vector>
#include <cassert>
using namespace std;
int mod(int a,int n){a%=n;return a<0?a+n:a;}
using Mat=array<array<int,11>,3>;using Key=array<int,22>;
Key canon(const Mat&e){
#ifdef AUDIT
 Key best;best.fill(100);
 array<array<int,11>,6>left{},right{};
 for(int j=0;j<11;j++){
  int u=e[1][j]-e[0][j],v=e[2][j]-e[0][j];
  left[0][j]=u;right[0][j]=v;left[1][j]=v;right[1][j]=u;
  left[2][j]=-u;right[2][j]=v-u;left[3][j]=v-u;right[3][j]=-u;
  left[4][j]=-v;right[4][j]=u-v;left[5][j]=u-v;right[5][j]=-v;
 }
 for(int t=0;t<6;t++)for(int b=0;b<11;b++){
  array<pair<int,int>,11>cols;for(int j=0;j<11;j++){cols[j]={mod(left[t][j]-left[t][b],40),mod(right[t][j]-right[t][b],40)};for(int k=j;k>0&&cols[k]<cols[k-1];k--)swap(cols[k],cols[k-1]);}
  Key key;for(int j=0;j<11;j++){key[2*j]=cols[j].first;key[2*j+1]=cols[j].second;}best=min(best,key);
 }return best;
#else
 Key best;best.fill(100);
 for(int a=0;a<3;a++)for(int b=0;b<3;b++)if(a!=b){int c=3-a-b;
  for(int base=0;base<11;base++){
   array<pair<int,int>,11>cols;
   for(int j=0;j<11;j++)cols[j]={mod(e[b][j]-e[a][j]-e[b][base]+e[a][base],40),mod(e[c][j]-e[a][j]-e[c][base]+e[a][base],40)};
   sort(cols.begin(),cols.end());Key key;for(int j=0;j<11;j++){key[2*j]=cols[j].first;key[2*j+1]=cols[j].second;}best=min(best,key);
  }
 }return best;
#endif
}

int main(int argc,char**argv){
 assert(argc==2);string root=argv[1];ifstream in(root+"/labels_input.txt");ofstream stats(root+"/labels_counts.txt"),families(root+"/families.txt"),zero(root+"/zero_seeds.txt"),special(root+"/exceptional_seeds.txt");int n;in>>n;assert(n==452);
 vector<array<int,5>>perms;array<int,4>p={1,2,3,4};do{perms.push_back({0,p[0],p[1],p[2],p[3]});}while(next_permutation(p.begin(),p.end()));
 map<Key,set<int>>allzero;long long total=0;int inv[5]={0,1,3,2,4};
 for(int ci=0;ci<n;ci++){
  int id,dim,rank;in>>id>>dim>>rank;vector<int>diag(rank);int components=1;for(int&d:diag){in>>d;components*=d;}
  array<array<int,11>,3>support;for(auto&r:support)for(int&x:r)in>>x;
  int U[7][7],V[6][6],g[7],Ug[7]={};for(auto&r:U)for(int&x:r)in>>x;for(auto&r:V)for(int&x:r)in>>x;for(int&x:g)in>>x;
  for(int i=0;i<7;i++)for(int j=0;j<7;j++)Ug[i]+=U[i][j]*g[j];
  int which[2][11],signs[2][11]={};for(auto&r:which)fill(r,r+11,-1);
  for(int r=0;r<2;r++)for(int b=0;b<3;b++)for(int k=0;k<2;k++){int j=support[r][5+2*b+k];which[r][j]=3*r+b;signs[r][j]=k;}
  int RR[7][6]={},T[7][7]={};for(int i=0;i<7;i++)T[i][i]=1;int rr=0;
  for(int k=1;k<5;k++){int j=support[2][k],a=support[2][0];for(int r=0;r<2;r++){int sign=r?-1:1;if(which[r][j]>=0)RR[rr][which[r][j]]+=sign;if(which[r][a]>=0)RR[rr][which[r][a]]-=sign;}rr++;}
  for(int k=0;k<3;k++){int j=support[2][6+2*k],a=support[2][5+2*k];for(int r=0;r<2;r++){int sign=r?-1:1;if(which[r][j]>=0)RR[rr][which[r][j]]+=sign;if(which[r][a]>=0)RR[rr][which[r][a]]-=sign;}rr++;}
  for(auto&r:RR)for(int&x:r)x=mod(x,5);int rrank=0;vector<int>piv;
  for(int col=0;col<6;col++){int k=rrank;while(k<7&&!RR[k][col])k++;if(k==7)continue;for(int j=0;j<6;j++)swap(RR[k][j],RR[rrank][j]);for(int j=0;j<7;j++)swap(T[k][j],T[rrank][j]);int z=inv[RR[rrank][col]];for(int&x:RR[rrank])x=x*z%5;for(int&x:T[rrank])x=x*z%5;for(int i=0;i<7;i++)if(i!=rrank){int a=RR[i][col];for(int j=0;j<6;j++)RR[i][j]=mod(RR[i][j]-a*RR[rrank][j],5);for(int j=0;j<7;j++)T[i][j]=mod(T[i][j]-a*T[rrank][j],5);}piv.push_back(col);rrank++;}assert(rrank==rank);
  long long count=0;set<Key>distinct;int maxlevel=20;
  for(int pa=0;pa<24;pa++)for(int pb=0;pb<24;pb++)for(int pc=0;pc<24;pc++){
   int labs[3][11]={};for(int r=0;r<3;r++)for(int k=0;k<5;k++)labs[r][support[r][k]]=perms[r==0?pa:r==1?pb:pc][k];
   int f[7],Uf[7]={},line=0;for(int k=1;k<5;k++){int j=support[2][k],a=support[2][0];f[line++]=labs[2][j]-labs[0][j]+labs[1][j]-labs[2][a]+labs[0][a]-labs[1][a];}
   for(int b=0;b<3;b++){int a=support[2][5+2*b],j=support[2][6+2*b];f[line++]=-labs[0][j]+labs[1][j]+labs[0][a]-labs[1][a];}
   assert(line==7);bool ok=true;for(int i=0;i<7;i++){for(int j=0;j<7;j++)Uf[i]+=U[i][j]*f[j];Uf[i]=mod(Uf[i],5);if(i>=rank&&Uf[i])ok=false;}
#ifdef AUDIT
   int Tf[7]={};bool gauss=true;for(int i=0;i<7;i++){for(int j=0;j<7;j++)Tf[i]+=T[i][j]*f[j];Tf[i]=mod(Tf[i],5);if(i>=rank&&Tf[i])gauss=false;}assert(gauss==ok);ok=gauss;
#endif
   if(!ok)continue;
   count++;
   for(int k=0;k<components;k++){
    int code=k;array<int,6>ys{},xs{};for(int i=0;i<rank;i++){int ki=code%diag[i];code/=diag[i];ys[i]=8*mod(inv[diag[i]%5]*Uf[i],5)+(20/diag[i])*(Ug[i]+2*ki);}
    for(int i=0;i<6;i++)for(int j=0;j<rank;j++)xs[i]+=V[i][j]*ys[j];
#ifdef AUDIT
    array<int,6>original5{},gaussian5{};for(int i=0;i<6;i++)for(int j=0;j<rank;j++)original5[i]+=V[i][j]*mod(inv[diag[j]%5]*Uf[j],5);for(int i=0;i<rank;i++)gaussian5[piv[i]]=Tf[i];
    if(dim==0){for(int i=0;i<6;i++)assert(mod(original5[i]-gaussian5[i],5)==0);}
    // Free-coordinate choices give the same affine coset. Keep the original
    // representative for stable data bytes; the independent solve certifies it.
#endif

    Mat e{};for(int r=0;r<2;r++)for(int j=0;j<11;j++)e[r+1][j]=mod(8*labs[r][j]+20*signs[r][j]+(which[r][j]>=0?xs[which[r][j]]:0),40);
    if(dim==0){auto key=canon(e);distinct.insert(key);allzero[key].insert(id);if(id==826){special<<pa<<' '<<pb<<' '<<pc<<' '<<k;for(auto x:key)special<<' '<<x;special<<'\n';}}
    else{families<<id<<' '<<dim<<' '<<pa<<' '<<pb<<' '<<pc<<' '<<k;for(int r=1;r<3;r++)for(int j=0;j<11;j++)families<<' '<<e[r][j];for(int r=0;r<2;r++)for(int j=0;j<11;j++)for(int d=rank;d<6;d++)families<<' '<<(which[r][j]>=0?V[which[r][j]][d]:0);families<<'\n';}
   }
  }
  stats<<id<<' '<<dim<<' '<<components<<' '<<count<<' '<<count*components<<' '<<distinct.size()<<'\n';total+=count*components;
  if(ci%20==0)cerr<<ci<<" classes; "<<total<<" labelled components; "<<allzero.size()<<" fixed triples"<<endl;
 }
 for(auto&[key,classes]:allzero){for(int x:key)zero<<x<<' ';zero<<classes.size();for(int id:classes)zero<<' '<<id;zero<<'\n';}
 cout<<"classes "<<n<<" labelled_components "<<total<<" fixed_triples "<<allzero.size()<<endl;
}
