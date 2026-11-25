// Independent upper bound, using a second split prime and direct K4 detection.
#include <algorithm>
#include <array>
#include <bitset>
#include <cassert>
#include <fstream>
#include <iostream>
#include <map>
#include <set>
#include <vector>
using namespace std;
constexpr int P=3361,M=840,L=5040,MAXN=300;
using Row=array<int,11>;
using Family=array<Row,3>; // free group, exponent of z, fixed root
int mod(int x,int n){return (x%n+n)%n;}
array<int,M> powers;
void setup(){
    for(int g=2;g<P;g++){
        int t=1;set<int>seen;
        for(int i=0;i<M;i++){powers[i]=t;seen.insert(t);t=t*g%P;}
        if(t==1&&seen.size()==M)return;
    }
    assert(false);
}
bool allowed(const Family&x,const Family&y){
    int first=-1,other=-1;array<long long,2>sum{},conjugate{};
    for(int j=0;j<11;j++){
        // Equal triples have the same unknown phase in the ratio x/y.
        int cell=8*(x[1][j]-y[1][j]+4)+2*x[0][j]+y[0][j];
        if(first<0)first=cell;
        if(cell!=first){if(other<0)other=cell;else if(other!=cell)return true;}
        int k=mod(x[2][j]-y[2][j],M),slot=int(cell!=first);
        sum[slot]+=powers[k];conjugate[slot]+=powers[mod(-k,M)];
    }
    if(other<0)return sum[0]%P==0 && conjugate[0]%P==0;
    return (sum[0]*conjugate[0]-sum[1]*conjugate[1])%P==0;
}
int main(int argc,char**argv){
    assert(argc==3);setup();ifstream in(argv[1]);ofstream out(argv[2]);assert(in&&out);
    int id,q;Family f;map<int,map<int,vector<Family>>> cases;
    while(in>>id>>q){for(auto&r:f)for(int&x:r)in>>x;assert(in);cases[id][q].push_back(f);}
    for(auto&[id,part]:cases){
        part[-1];for(auto&[q,rows]:part){
            set<Family> vertices;
            for(auto f:part[-1]){
                for(int j=0;j<11;j++){
                    if(q>=0){f[2][j]=mod(f[2][j]+q*f[1][j],L);f[1][j]=0;}
                    assert(f[2][j]%6==0);f[2][j]/=6;
                }
                vertices.insert(f);
            }
            if(q>=0)for(auto f:rows){for(int&e:f[2]){assert(e%6==0);e/=6;}vertices.insert(f);}
            vector<Family> v(vertices.begin(),vertices.end());int n=v.size();assert(n<=MAXN);
            vector<bitset<MAXN>> higher(n);int edges=0,maximum=(n>0);
            for(int i=0;i<n;i++)for(int j=i+1;j<n;j++)if(allowed(v[i],v[j])){higher[i].set(j);edges++;maximum=2;}
            for(int i=0;i<n;i++)for(int j=i+1;j<n;j++)if(higher[i][j]){
                auto common=higher[i]&higher[j];
                if(common.any())maximum=3;
                for(int k=j+1;k<n;k++)if(common[k])assert((common&higher[k]).none());
            }
            out<<id<<' '<<q<<' '<<n<<' '<<edges<<' '<<maximum<<'\n';
        }
    }
}
