// All prime-cycle covers of x/v(z) for a certified row family x of the fixed pair.
#include <algorithm>
#include <array>
#include <cassert>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <numeric>
#include <set>
#include <tuple>
#include <vector>
using namespace std;
constexpr int N=11,M=840,L=5040;
#ifdef AUDIT
constexpr int P=3361;
#else
constexpr int P=2521;
#endif
using Row=array<int,N>;
int mod(int x,int n=L){x%=n;return x<0?x+n:x;}
struct Record{int z;Row g,a,e;bool operator<(const Record&r)const{return std::tie(z,g,a,e)<std::tie(r.z,r.g,r.a,r.e);}};
struct Pair{int n;Row g,e;};
struct Target{int id;Row e,b;};
struct Potential{int group,a,e;};
#ifdef AUDIT
// A forest on the at most three phase groups. Each non-tree edge supplies
// one linear congruence on the old phase. Keep its exact solution coset.
struct State{
    array<int,3> parent{0,1,2},slope{},offset{},roots{};
    int nz=-1;
    Potential get(int j)const{
        int a=0,e=0;
        while(parent[j]!=j){a+=slope[j];e+=offset[j];j=parent[j];}
        return {j,a,mod(e)};
    }
    bool link(int j,int h,int b,int v){
        auto left=get(j),right=get(h);
        int a=left.a-right.a-b;
        int e=mod(v+right.e-left.e);
        if(left.group!=right.group){
            parent[right.group]=left.group;
            slope[right.group]=a;offset[right.group]=mod(-e);
            return true;
        }
        if(!a)return e==0;
        assert(abs(a)<=3 && e%6==0);
        if(nz<0){
            // Solve a*t=e in turns with denominator L. The input constants
            // are in mu_840 and L is divisible by 840*abs(a).
            int n=abs(a),step=L/n;
            int base=mod((a>0?e:-e)/n,step);nz=n;
            for(int i=0;i<n;i++)roots[i]=base+i*step;
        }else{
            int n=0;for(int i=0;i<nz;i++)if(mod(a*roots[i]-e)==0)roots[n++]=roots[i];
            nz=n;
        }
        return nz>0;
    }
};
#else
struct State{
    array<int,3> group{0,1,2},a{},e{};
    int nz=-1;array<int,3> roots{};
    Potential get(int j)const{return {group[j],a[j],e[j]};}
    bool restrict_z(int k,int v){
        if(!k)return mod(v)==0;
        assert(abs(k)<=3);if(nz<0){
            if(k<0){k=-k;v=-v;}
            v=mod(v);assert(v%k==0);
            nz=k;for(int i=0;i<k;i++)roots[i]=(v+i*L)/k;
        }else{
            int n=0;for(int i=0;i<nz;i++)if(mod(k*roots[i]-v)==0)roots[n++]=roots[i];
            nz=n;
        }
        return nz>0;
    }
    bool link(int j,int h,int b,int v){
        // phase[j]/phase[h] = z^b times the root with exponent v.
        b+=a[h]-a[j];v=mod(v+e[h]-e[j]);
        int cj=group[j],ch=group[h];
        if(cj==ch)return restrict_z(b,-v);
        for(int i=0;i<3;i++)if(group[i]==cj){group[i]=ch;a[i]+=b;e[i]=mod(e[i]+v);}
        return true;
    }
};
#endif
array<int,M> pw;
void make_table(){
    for(int g=2;;g++){
        pw[0]=1;for(int j=1;j<M;j++)pw[j]=int64_t(pw[j-1])*g%P;
        if(int64_t(pw.back())*g%P==1&&pw[M/2]!=1&&pw[M/3]!=1&&pw[M/5]!=1&&pw[M/7]!=1)break;
    }
}
bool possible_point(const Pair&x,const Target&t){
    if(x.n>1)return true;
    array<int,2>plus{},minus{};
    for(int j=0;j<N;j++){
        int k=mod(x.e[j]-21*t.e[j],M);plus[t.b[j]]+=pw[k];minus[t.b[j]]+=pw[mod(-k,M)];
    }
    return (int64_t(plus[0])*minus[0]-int64_t(plus[1])*minus[1])%P==0;
}
struct Search{
    const Pair&x;const Target&t;set<Record>&output;
    void leaf(const State&s){
        for(int i=0;i<(s.nz<0?1:s.nz);i++){
            Record r{};r.z=s.nz<0?-1:s.roots[i];array<int,3> num{-1,-1,-1},aa{},ee{};int count=0;
            for(int j=0;j<N;j++){
                auto potential=s.get(x.g[j]);
                int g=potential.group,a=potential.a,e=mod(6*x.e[j]+potential.e);
                if(r.z>=0){e=mod(e+a*r.z);a=0;}
                if(num[g]<0){num[g]=count++;aa[g]=a;ee[g]=e;}
                r.g[j]=num[g];r.a[j]=a-aa[g];r.e[j]=mod(e-ee[g]);
            }
            output.insert(r);
        }
    }
    void cover(int mask,const State&s){
        if(!mask){leaf(s);return;}
        int n=__builtin_popcount((unsigned)mask);
#ifdef AUDIT
        int h=31-__builtin_clz((unsigned)mask);
        constexpr array<int,4> primes{7,5,3,2};
#else
        int h=__builtin_ctz((unsigned)mask);
        constexpr array<int,4> primes{2,3,5,7};
#endif
        for(int p:primes)if(p<=n&&(p==n||n-p>=2))cycle(mask^(1<<h),h,p,1,s);
    }
    void cycle(int mask,int h,int p,int k,const State&s){
        if(k==p){cover(mask,s);return;}
        for(int j=0;j<N;j++)if(mask>>j&1){
            int v=L/p*k+6*(x.e[h]-x.e[j])+126*(t.e[j]-t.e[h]);
            State next=s;
            if(next.link(x.g[j],x.g[h],t.b[j]-t.b[h],v))cycle(mask^(1<<j),h,p,k+1,next);
        }
    }
};
int main(int argc,char**argv){
    assert(argc==4);make_table();ifstream xf(argv[1]),tf(argv[2]);ofstream out(argv[3]);assert(xf&&tf&&out);
    vector<Pair> xs;char profile;long long index;Pair x;
    while(xf>>profile>>index>>x.n){for(int&g:x.g)xf>>g;for(int&e:x.e)xf>>e;assert(xf&&x.n<=3);xs.push_back(x);}
    assert(xs.size()==339008);Target t;while(tf>>t.id){
        for(int&e:t.e)tf>>e;for(int&b:t.b)tf>>b;assert(tf);
        set<Record> records;for(auto&x:xs)if(possible_point(x,t))Search{x,t,records}.cover(2047,State{});
        set<int> points;int free=0;
        for(const Record&r:records){
            if(r.z>=0)points.insert(r.z);else free++;
            out<<t.id<<' '<<r.z;for(int g:r.g)out<<' '<<g;for(int a:r.a)out<<' '<<a;for(int e:r.e)out<<' '<<e;out<<'\n';
        }
        cout<<t.id<<" records "<<records.size()<<" generic "<<free<<" points "<<points.size()<<endl;
    }
}
