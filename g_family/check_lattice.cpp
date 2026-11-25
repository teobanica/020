#define NO_MAIN
#include "search.cpp"
#include <random>
int main(){mt19937 rng(112233);long long checked=0,good=0;for(int z=0;z<30000;z++){vector<V>a;int n=1+rng()%10;array<int,3>w;for(auto&x:w)x=rng()%M;for(int i=0;i<n;i++){V v;for(int k=0;k<3;k++)v[k]=int(rng()%9)-4;v[3]=z%2?rng()%M:mod(v[0]*w[0]+v[1]*w[1]+v[2]*w[2]);a.push_back(v);}Lattice l;bool feasible=true;for(auto v:a)if(!l.add(v)){feasible=false;break;}assert(feasible==independent(a));if(feasible){good++;l.canonical();auto key=l.r;for(int j=0;j<3;j++){shuffle(a.begin(),a.end(),rng);Lattice x;for(auto v:a)assert(x.add(v));x.canonical();assert(x.r==key);}}checked++;}cout<<"independent diagonal reduction comparisons "<<checked<<" feasible and order-invariant "<<good<<endl;}
