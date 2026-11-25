// Exact generic fourth-row scan for the common-block CCC triples.
// V=(s,r,t,b,root), with root denominator 210. All phase variables are free.
#include <algorithm>
#include <array>
#include <cassert>
#include <fstream>
#include <iostream>
#include <numeric>
#include <set>
#include <vector>
using namespace std;
constexpr int M=210, N=11, D=5;
using V=array<int,D>;
using Key=array<int,N*(D+1)>;
int mod(int x){x%=M;return x<0?x+M:x;}
array<set<V>,2> exceptional;
array<long long,2> nodes{},leaves{};
bool identity(V v,int stage){
    if(any_of(v.begin(),v.end()-1,[](int x){return x!=0;})){
        // The closing equation is monomial(v)=1; record character=root.
        v[D-1]=-v[D-1];
        int j=0;while(!v[j])++j;
        if(v[j]<0)for(auto &x:v)x=-x;
        v[D-1]=mod(v[D-1]);exceptional[stage].insert(v);return false;
    }
    return mod(v[D-1])==0;
}
// Audit implementation: maintain every component label and potential explicitly.
struct State{
    array<int,5> component;
    array<V,5> potential{};
    State(){iota(component.begin(),component.end(),0);}
    pair<int,V> find(int j)const{return {component[j],potential[j]};}
    bool link(int a,int b,V v,int stage){
        int ca=component[a],cb=component[b];
        for(int k=0;k<D;++k)v[k]+=-potential[a][k]+potential[b][k];
        v[D-1]=mod(v[D-1]);
        if(ca==cb)return identity(v,stage);
        // Attach the larger component label to the smaller one.
        if(ca<cb){
            for(int j=0;j<5;++j)if(component[j]==cb){
                component[j]=ca;
                for(int k=0;k<D;++k)potential[j][k]-=v[k];
                potential[j][D-1]=mod(potential[j][D-1]);
            }
        }else{
            for(int j=0;j<5;++j)if(component[j]==ca){
                component[j]=cb;
                for(int k=0;k<D;++k)potential[j][k]+=v[k];
                potential[j][D-1]=mod(potential[j][D-1]);
            }
        }
        return true;
    }
};
struct Family{array<int,N> group{};array<V,N> value{};};
array<array<V,N>,2> old;
vector<array<int,5>> heads;
set<Key> pairs;
vector<set<Key>> finals;
vector<int> sizes;
vector<vector<int>> blocks;
Family initial;
int current=0;
long long tests=0;
ofstream output;

Family normalized(const Family& f,const State& s){
    Family g;array<int,5> number;number.fill(-1);array<V,5> anchor{};int count=0;
    for(int j=0;j<N;++j){
        auto [r,v]=s.find(f.group[j]);
        for(int k=0;k<D;++k)v[k]+=f.value[j][k];
        if(number[r]<0){number[r]=count++;anchor[r]=v;}
        g.group[j]=number[r];
        for(int k=0;k<D;++k)g.value[j][k]=v[k]-anchor[r][k];
        g.value[j][D-1]=mod(g.value[j][D-1]);
    }
    return g;
}
Key key(const Family& f){
    Key result;int k=0;
    for(int j=0;j<N;++j){result[k++]=f.group[j];for(int x:f.value[j])result[k++]=x;}
    assert(k==(int)result.size());return result;
}
void cover(int mask,const State& s,const Family& f,int stage);
void leaf(const State& s,const Family& f,int stage){
    ++leaves[stage];auto g=normalized(f,s);auto k=key(g);
    if(stage==0){
        if(!pairs.insert(k).second)return;
        for(current=0;current<(int)heads.size();++current){
            for(int j=0;j<N;++j){
                if(j<5)old[1][j]={0,0,0,0,42*heads[current][j]};
                else if(j<8)old[1][j]={0,0,1,0,mod(140*(j-5))};
                else old[1][j]={0,0,0,1,mod(140*(j-8))};
            }
            State next;cover((1<<N)-1,next,g,1);
        }
    }else if(finals[current].insert(k).second){
        output<<current<<' ';for(int x:k)output<<x<<' ';output<<'\n';
    }
}
void cycle(int mask,int anchor,int prime,int label,const State& s,const Family& f,int stage){
    ++nodes[stage];if(label==prime){cover(mask,s,f,stage);return;}
    for(int j=0;j<N;++j)if(mask>>j&1){
        V v;
        for(int k=0;k<D;++k)v[k]=f.value[anchor][k]-f.value[j][k]+old[stage][j][k]-old[stage][anchor][k];
        v[D-1]+=M/prime*label;State next=s;
        if(next.link(f.group[j],f.group[anchor],v,stage))
            cycle(mask^(1<<j),anchor,prime,label+1,next,f,stage);
    }
}
void cover(int mask,const State& s,const Family& f,int stage){
    if(!mask){leaf(s,f,stage);return;}
    int n=__builtin_popcount((unsigned)mask),j=31-__builtin_clz((unsigned)mask);
    for(int p:{7,5,3,2})if(p<=n&&(p==n||n-p>=2))
        cycle(mask^(1<<j),j,p,1,s,f,stage);
}
void test(){
    ++tests;State s;cover((1<<N)-1,s,initial,0);
    if(tests%10000==0)cerr<<tests<<" pairs "<<pairs.size()<<" nodes "<<nodes[0]<<' '<<nodes[1]<<'\n';
}
void labels(int k){
    if(k==(int)blocks.size()){test();return;}
    const auto& block=blocks[k];vector<int> p(block.size()-1);iota(p.begin(),p.end(),1);
    do{
        initial.value[block[0]][D-1]=0;
        for(int j=1;j<(int)block.size();++j)initial.value[block[j]][D-1]=p[j-1]*M/block.size();
        labels(k+1);
    }while(next_permutation(p.begin(),p.end()));
}
void partition(int k,const vector<int>& rem){
    if(k==(int)sizes.size()){labels(0);return;}
    for(int mask=0;mask<(1<<(int)rem.size());++mask)if(__builtin_popcount((unsigned)mask)==sizes[k]){
        vector<int> block,rest;
        for(int j=0;j<(int)rem.size();++j)((mask>>j&1)?block:rest).push_back(rem[j]);
        if(k&&sizes[k]==sizes[k-1]&&block[0]<blocks.back()[0])continue;
        for(int j:block)initial.group[j]=k;
        blocks.push_back(block);partition(k+1,rest);blocks.pop_back();
    }
}
int main(int argc,char** argv){
    assert(argc==5);ifstream in(argv[1]);for(auto& v:old[0])for(auto& x:v)in>>x;assert(in);
    ifstream hf(argv[2]);array<int,5> h;
    while(hf>>h[0]){for(int j=1;j<5;++j)hf>>h[j];assert(hf);heads.push_back(h);}
    assert(heads.size()==3);finals.resize(heads.size());
    int profile=argv[3][0]-'B';assert(0<=profile&&profile<5);
    string path=argv[4];output.open(path+".writing");assert(output);
    sizes=vector<vector<int>>{{7,2,2},{5,3,3},{5,2,2,2},{3,3,3,2},{3,2,2,2,2}}[profile];
    vector<int> rem(N);iota(rem.begin(),rem.end(),0);partition(0,rem);
    output.flush();assert(output);output.close();assert(!rename((path+".writing").c_str(),path.c_str()));
    ofstream pairfile(path+".pairs");for(const auto& k:pairs){for(int x:k)pairfile<<x<<' ';pairfile<<'\n';}
    for(int stage=0;stage<2;++stage){
        ofstream ex(path+".exceptional"+to_string(stage));
        for(auto v:exceptional[stage]){for(int x:v)ex<<x<<' ';ex<<'\n';}
    }
    size_t total=0;for(const auto& f:finals)total+=f.size();
    cout<<"tests "<<tests<<" pairs "<<pairs.size()<<" final "<<total<<" nodes "<<nodes[0]<<' '<<nodes[1]
        <<" leaves "<<leaves[0]<<' '<<leaves[1]<<" exceptional "<<exceptional[0].size()<<' '<<exceptional[1].size()<<'\n';
    for(int i=0;i<(int)heads.size();++i)cout<<"head "<<i<<" fourth_rows "<<finals[i].size()<<'\n';
}
