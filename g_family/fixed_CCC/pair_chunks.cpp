#include <algorithm>
#include <array>
#include <cassert>
#include <fstream>
#include <iostream>
#include <numeric>
#include <set>
#include <vector>
#include <cstdio>
#include <cerrno>
#include <cstring>
using namespace std;
constexpr int M=840;
int mod(int a){a%=M;return a<0?a+M:a;}
struct State{
 array<int,5>parent{},e{};State(){iota(parent.begin(),parent.end(),0);}
 pair<int,int>find(int b)const{
#ifdef AUDIT
  return {parent[b],e[b]};
#else
  int v=0;while(parent[b]!=b){v+=e[b];b=parent[b];}return{b,mod(v)};
#endif
 }
 bool link(int b,int c,int v){auto B=find(b),C=find(c);v=mod(v-B.second+C.second);if(B.first==C.first)return v==0;
#ifdef AUDIT
  if(B.first<C.first){for(int j=0;j<5;j++)if(parent[j]==C.first){parent[j]=B.first;e[j]=mod(e[j]-v);}}
  else{for(int j=0;j<5;j++)if(parent[j]==B.first){parent[j]=C.first;e[j]=mod(e[j]+v);}}
#else
  parent[B.first]=C.first;e[B.first]=v;
#endif
  return true;
 }
};
vector<int>sizes,block(11),exponent(11),u(11);vector<vector<int>>blocks;
using Key=array<int,22>;set<Key>seen;ofstream out;long long start_index=1,tested=0,limit=1000000000,nodes=0,leaves=0;array<long long,6>dimensions{};char profile;
void leaf(const State&s){
 leaves++;Key key;array<int,5>groups,anchor;groups.fill(-1);int c=0;
 for(int j=0;j<11;j++){auto v=s.find(block[j]);if(groups[v.first]<0){groups[v.first]=c++;anchor[v.first]=mod(v.second+exponent[j]);}key[j]=groups[v.first];key[11+j]=mod(v.second+exponent[j]-anchor[v.first]);}
 if(!seen.insert(key).second)return;dimensions[c-1]++;out<<profile<<' '<<tested<<' '<<c;for(int x:key)out<<' '<<x;out<<'\n';
}
void cover(int mask,const State&s);
void cycle(int mask,int anchor,int prime,int label,const State&s){nodes++;if(label==prime){cover(mask,s);return;}
 for(int j=0;j<11;j++)if((mask>>j)&1){State t=s;int v=mod(label*(M/prime)+exponent[anchor]-exponent[j]+u[j]-u[anchor]);if(t.link(block[j],block[anchor],v))cycle(mask^(1<<j),anchor,prime,label+1,t);}
}
void cover(int mask,const State&s){if(!mask){leaf(s);return;}int n=__builtin_popcount((unsigned)mask),j=__builtin_ctz((unsigned)mask);mask^=1<<j;for(int p:{2,3,5,7})if(p<=n&&(n==p||n-p>=2))cycle(mask,j,p,1,s);}
void test(){if(tested>=limit)return;tested++;if(tested<start_index)return;State s;cover(2047,s);if(tested%100000==0)cerr<<profile<<' '<<tested<<" nodes "<<nodes<<" unique "<<seen.size()<<endl;}
void labels(int k){if(tested>=limit)return;if(k==(int)blocks.size()){test();return;}auto&s=blocks[k];vector<int>p(s.size()-1);iota(p.begin(),p.end(),1);do{exponent[s[0]]=0;for(int j=1;j<(int)s.size();j++)exponent[s[j]]=p[j-1]*M/s.size();labels(k+1);}while(next_permutation(p.begin(),p.end()));}
void partition(int k,vector<int>rem){if(tested>=limit)return;if(k==(int)sizes.size()){labels(0);return;}int n=rem.size();for(int mask=0;mask<(1<<n);mask++)if(__builtin_popcount((unsigned)mask)==sizes[k]){vector<int>chosen,rest;for(int j=0;j<n;j++)((mask>>j)&1?chosen:rest).push_back(rem[j]);if(k&&sizes[k]==sizes[k-1]&&chosen[0]<blocks.back()[0])continue;for(int j:chosen)block[j]=k;blocks.push_back(chosen);partition(k+1,rest);blocks.pop_back();}}
int main(int argc,char**argv){assert(argc>=4);int group=stoi(argv[1]);ifstream in(argv[2]);int n;in>>n;bool found=false;for(int r=0;r<n;r++){int id,nv;in>>id>>nv;vector<int>v(11);for(int&x:v)in>>x;if(id==group){for(int j=0;j<11;j++)u[j]=7*v[j];found=true;}for(int k=0;k<nv;k++){int idx;in>>idx;for(int&x:v)in>>x;}}assert(found);string final_path=argv[3],temporary_path=final_path+".writing";out.open(temporary_path);assert(out);if(argc>4)limit=stoll(argv[4]);if(argc>5)start_index=stoll(argv[5]);int selected=argc>6?argv[6][0]-'B':-1;vector<vector<int>>sz={{7,2,2},{5,3,3},{5,2,2,2},{3,3,3,2},{3,2,2,2,2}};
 for(int p=0;p<5;p++){if(selected>=0&&p!=selected)continue;profile='B'+p;sizes=sz[p];tested=nodes=leaves=0;dimensions.fill(0);auto before=seen.size();vector<int>rem(11);iota(rem.begin(),rem.end(),0);partition(0,rem);cout<<"profile "<<profile<<" tests "<<tested<<" nodes "<<nodes<<" leaves "<<leaves<<" unique "<<seen.size()-before<<" dimensions";for(auto x:dimensions)cout<<' '<<x;cout<<endl;}
 out.flush();if(!out){cerr<<"WRITE FAILURE "<<strerror(errno)<<endl;return 7;}out.close();if(!out)return 8;if(rename(temporary_path.c_str(),final_path.c_str()))return 9;return 0;
}
