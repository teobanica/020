#include "common.hpp"
void fixedstate(const State&s,int q){
 int c=s.components(sizes.size());auto base=s.find(0);vector<int>groups(5,-1),ee(11),aa(11);groups[base[0]]=0;int next=1;
 for(int j=0;j<11;j++){
  auto v=s.find(block[j]);if(groups[v[0]]<0)groups[v[0]]=next++;
  aa[j]=groups[v[0]];ee[j]=mod((v[1]-base[1])*q+(v[2]-base[2]+exponent[j])*(L/M),L);
 }
 vector<int>key={c};key.insert(key.end(),aa.begin(),aa.end());key.insert(key.end(),ee.begin(),ee.end());if(!seen.insert(key).second){duplicates++;return;}
 if(c==1){checkpoint_x(ee,q);return;}
 if(c==2){assert((1LL*q*25200)%L==0);univariate(aa,ee,25200,q);return;}
 unconnected++;families<<tests<<' '<<c<<" FIXED "<<q;for(int j=0;j<11;j++)families<<' '<<aa[j]<<' '<<ee[j];families<<'\n';
}

void leaf(const State&s){
 leaves++;int c=s.components(sizes.size());bydim[c-1+!s.fixed]++;
 if(s.fixed){for(int q:s.qs)fixedstate(s,q);return;}
 vector<int>aa(11),ee(11);auto base=s.find(0);
 for(int j=0;j<11;j++){auto v=s.find(block[j]);aa[j]=v[1]-base[1];ee[j]=mod(v[2]-base[2]+exponent[j],M)*(L/M);}
 if(c==1){vector<int>key={0};key.insert(key.end(),aa.begin(),aa.end());key.insert(key.end(),ee.begin(),ee.end());if(!seen.insert(key).second){duplicates++;return;}univariate(aa,ee,420,-1);return;}
 unconnected++;families<<tests<<' '<<c<<" FREE";for(int b=0;b<(int)sizes.size();b++){auto v=s.find(b);families<<' '<<v[0]<<' '<<v[1]<<' '<<v[2];}families<<'\n';
}
bool sizes_possible(int rem){return rem==0||rem>=2;}
void cover(int mask,const State&s);
void cycle(int mask,int anchor,int prime,int label,const State&s){
 nodes++;
 if(label==prime){cover(mask,s);return;}
 for(int j=0;j<11;j++)if((mask>>j)&1){
  State t=s;int d=int(anchor>=5)-int(j>=5);
  int v=mod(label*(M/prime)+exponent[anchor]-exponent[j]+21*(H[1][j]-H[1][anchor]),M);
  if(t.link(block[j],block[anchor],d,v))cycle(mask^(1<<j),anchor,prime,label+1,t);
 }
}
void cover(int mask,const State&s){
 if(!mask){leaf(s);return;}
 int n=__builtin_popcount((unsigned)mask),j=__builtin_ctz((unsigned)mask);mask^=1<<j;
 for(int p:{2,3,5,7})if(p<=n&&sizes_possible(n-p))cycle(mask,j,p,1,s);
}
void test(){
 if(tests>=limit)return;tests++;seen.clear();long long before=unconnected+freezero;auto start=chrono::steady_clock::now();State s;cover(2047,s);if(unconnected+freezero>before){inputs<<tests<<' '<<profile;for(int b:block)inputs<<' '<<b;for(int e:exponent)inputs<<' '<<e;inputs<<'\n';}
 if(tests%1000==0)cerr<<tests<<" leaves "<<leaves<<" nodes "<<nodes<<" family "<<unconnected<<" free "<<freezero<<" roots "<<rootcandidates<<" passed "<<rootpassed<<endl;
}
void labels(int k){if(tests>=limit)return;if(k==(int)blocks.size()){test();return;}auto &s=blocks[k];vector<int>p(s.size()-1);iota(p.begin(),p.end(),1);do{exponent[s[0]]=0;for(int j=1;j<(int)s.size();j++)exponent[s[j]]=p[j-1]*M/s.size();labels(k+1);}while(next_permutation(p.begin(),p.end()));}
void partition(int k,vector<int>rem){if(tests>=limit)return;if(k==(int)sizes.size()){labels(0);return;}int n=rem.size();for(int mask=0;mask<(1<<n);mask++)if(__builtin_popcount((unsigned)mask)==sizes[k]){vector<int>chosen,rest;for(int j=0;j<n;j++)((mask>>j)&1?chosen:rest).push_back(rem[j]);if(k&&sizes[k]==sizes[k-1]&&chosen[0]<blocks.back()[0])continue;for(int j:chosen)block[j]=k;blocks.push_back(chosen);partition(k+1,rest);blocks.pop_back();}}
int main(int argc,char**argv){if(argc>1)limit=stoll(argv[1]);if(argc>2)profile=argv[2][0];vector<vector<int>>sz={{7,2,2},{5,3,3},{5,2,2,2},{3,3,3,2},{3,2,2,2,2}};sizes=sz[profile-'B'];families.open("families.txt");points.open("points.txt");inputs.open("inputs.txt");vector<int>rem(11);iota(rem.begin(),rem.end(),0);partition(0,rem);cout<<"tests "<<tests<<" leaves "<<leaves<<" nodes "<<nodes<<" families "<<unconnected<<" free "<<freezero<<" roots "<<rootcandidates<<" passed "<<rootpassed<<" actual "<<actual<<" duplicates "<<duplicates<<" dimensions";for(auto n:bydim)cout<<' '<<n;cout<<endl;return 0;}
