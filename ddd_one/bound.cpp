// Necessary orthogonality tests. A coloring or a clique search bounds the
// number of simultaneously addable rows even when individual phases differ.
#include <algorithm>
#include <array>
#include <cassert>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <map>
#include <numeric>
#include <set>
#include <tuple>
#include <vector>
using namespace std;
using Row=array<int,11>;
constexpr int M=840,L=5040,P=2521;
int mod(int x,int n){return (x%n+n)%n;}
struct Family{
    Row g,a,e;
    bool operator<(const Family &r)const{return tie(g,a,e)<tie(r.g,r.a,r.e);}
};
array<int,M>pw;
void table(){
    for(int g=2;;g++){pw[0]=1;for(int j=1;j<M;j++)pw[j]=int64_t(pw[j-1])*g%P;
        if(int64_t(pw.back())*g%P==1&&pw[M/2]!=1&&pw[M/3]!=1&&pw[M/5]!=1&&pw[M/7]!=1)break;}
}
bool edge(const Family&x,const Family&y){
    array<int,40> plus{},minus{};int first=-1,second=-1;
    for(int j=0;j<11;j++){
        int cell=4*(4+x.a[j]-y.a[j])+2*x.g[j]+y.g[j];assert(cell>=0&&cell<40);
        if(first==-1)first=cell;
        else if(cell!=first&&second==-1)second=cell;
        else if(cell!=first&&cell!=second)return true;
        int k=mod(x.e[j]-y.e[j],M);plus[cell]+=pw[k];minus[cell]+=pw[mod(-k,M)];
    }
    if(second==-1)return plus[first]%P==0;
    return (int64_t(plus[first])*minus[first]-int64_t(plus[second])*minus[second])%P==0;
}
struct Graph{
    vector<vector<char>> adj;int best=0;
    int greedy(const vector<int>&vv,vector<int>&order,vector<int>&col){
        vector<int>rest=vv;int n=0;
        while(!rest.empty()){
            vector<int>used,next;n++;
            for(int x:rest){
                if(any_of(used.begin(),used.end(),[&](int y){return adj[x][y];}))next.push_back(x);
                else{used.push_back(x);order.push_back(x);col.push_back(n);}
            }
            rest=move(next);
        }
        return n;
    }
    void expand(const vector<int>&vv,int size){
        vector<int>order,col;greedy(vv,order,col);
        for(int i=int(order.size())-1;i>=0;i--){
            if(size+col[i]<=best||best>=8)return;
            int v=order[i];vector<int>next;for(int j=0;j<i;j++)if(adj[v][order[j]])next.push_back(order[j]);
            best=max(best,size+1);if(!next.empty())expand(next,size+1);
        }
    }
};
int main(int argc,char**argv){
    assert(argc==4);table();ifstream in(argv[1]);ofstream out(argv[2]);assert(in&&out);int only=stoi(argv[3]);
    map<int,map<int,vector<Family>>> all;int id,z;Family f;
    while(in>>id>>z){for(auto p:{&f.g,&f.a,&f.e})for(int&v:*p)in>>v;assert(in&&*max_element(f.g.begin(),f.g.end())<=1);all[id][z].push_back(f);}
    for(auto &[id,by]:all){
        if(only>=0&&id!=only)continue;
        vector<int>points{-1};for(auto &[k,fs]:by)if(k>=0)points.push_back(k);
        for(int z:points){
            set<Family>rows;
            for(auto f:by[-1]){for(int j=0;j<11;j++){if(z>=0){f.e[j]=mod(f.e[j]+f.a[j]*z,L);f.a[j]=0;}assert(f.e[j]%6==0);f.e[j]/=6;}rows.insert(f);}
            if(z>=0)for(auto f:by[z]){for(int&v:f.e){assert(v%6==0);v/=6;}rows.insert(f);}
            vector<Family>v(rows.begin(),rows.end());int n=v.size();Graph graph;graph.adj.assign(n,vector<char>(n));int edges=0;
            for(int i=0;i<n;i++)for(int j=0;j<i;j++)if(edge(v[i],v[j])){graph.adj[i][j]=graph.adj[j][i]=1;edges++;}
            vector<int>vv(n),order,col;iota(vv.begin(),vv.end(),0);int colors=graph.greedy(vv,order,col);
            graph.expand(vv,0);out<<id<<' '<<z<<' '<<n<<' '<<edges<<' '<<colors<<' '<<graph.best<<'\n';
        }
        cout<<id<<" graphs "<<points.size()<<endl;
    }
}
