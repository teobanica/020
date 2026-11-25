#include "common.hpp"
#include <map>
#include <sstream>
int active_row=1;bool block_gauge_allowed=false;long long block_gauges=0;
void refine_family(const vector<int>&groups,const vector<int>&coeff,int count);
void fixedstate(const State&s,int q){
 int c=s.components(sizes.size());auto base=s.find(0);vector<int>groups(5,-1),ee(11),aa(11);groups[base[0]]=0;int next=1;
 for(int j=0;j<11;j++){
  auto v=s.find(block[j]);if(groups[v[0]]<0)groups[v[0]]=next++;
  aa[j]=groups[v[0]];ee[j]=mod((v[1]-base[1])*q+(v[2]-base[2]+exponent[j])*(L/M),L);
 }
 vector<int>key={c};key.insert(key.end(),aa.begin(),aa.end());key.insert(key.end(),ee.begin(),ee.end());if(!seen.insert(key).second){duplicates++;return;}
 if(c==1){checkpoint_x(ee,q);return;}
 if(c==2){assert((1LL*q*25200)%L==0);univariate(aa,ee,25200,q);return;}
 refine_family(aa,ee,c);
}

void leaf(const State&s){
 leaves++;int c=s.components(sizes.size());bydim[c-1+!s.fixed]++;
 if(s.fixed){for(int q:s.qs)fixedstate(s,q);return;}
 vector<int>aa(11),ee(11);auto base=s.find(0);
 for(int j=0;j<11;j++){auto v=s.find(block[j]);aa[j]=v[1]-base[1];ee[j]=mod(v[2]-base[2]+exponent[j],M)*(L/M);}
 if(c==1){vector<int>key={0};key.insert(key.end(),aa.begin(),aa.end());key.insert(key.end(),ee.begin(),ee.end());if(!seen.insert(key).second){duplicates++;return;}univariate(aa,ee,420,-1);return;}

 vector<int>gg(11);array<int,5>group,anchor_a,anchor_e;group.fill(-1);int ng=0;
 for(int j=0;j<11;j++){
  auto v=s.find(block[j]);if(group[v[0]]<0){group[v[0]]=ng++;anchor_a[v[0]]=v[1];anchor_e[v[0]]=v[2]+exponent[j];}
  gg[j]=group[v[0]];aa[j]=v[1]-anchor_a[v[0]];ee[j]=mod(v[2]+exponent[j]-anchor_e[v[0]],M)*(L/M);
 }
 bool depends=any_of(aa.begin(),aa.end(),[](int a){return a!=0;});
 if(depends&&block_gauge_allowed){
  for(int shift=-4;shift<=4;shift++){
   array<int,5>values;values.fill(999);bool ok=true;
   for(int j=0;j<11;j++){int value=aa[j]-shift*int(j>=5);if(values[gg[j]]==999)values[gg[j]]=value;else if(values[gg[j]]!=value)ok=false;}
   if(ok){depends=false;block_gauges++;break;}
  }
 }
 if(depends){
  unconnected++;families<<tests<<" GENUINE_Q_FREE "<<active_row<<' '<<c;for(int j=0;j<11;j++)families<<' '<<gg[j]<<' '<<aa[j]<<' '<<ee[j];families<<'\n';return;
 }
 if(c==2){univariate(gg,ee,420,0);return;}
 refine_family(gg,ee,c);

}
bool sizes_possible(int rem){return rem==0||rem>=2;}
void cover(int mask,const State&s);
void cycle(int mask,int anchor,int prime,int label,const State&s){
 nodes++;
 if(label==prime){cover(mask,s);return;}
 for(int j=0;j<11;j++)if((mask>>j)&1){
  State t=s;int d=int(anchor>=5)-int(j>=5);
  int v=mod(label*(M/prime)+exponent[anchor]-exponent[j]+21*(H[active_row][j]-H[active_row][anchor]),M);
  if(t.link(block[j],block[anchor],d,v))cycle(mask^(1<<j),anchor,prime,label+1,t);
 }
}
void cover(int mask,const State&s){
 if(!mask){leaf(s);return;}
 int n=__builtin_popcount((unsigned)mask),j=__builtin_ctz((unsigned)mask);mask^=1<<j;
 for(int p:{2,3,5,7})if(p<=n&&sizes_possible(n-p))cycle(mask,j,p,1,s);
}

long long refined=0,refine_skipped=0;set<vector<int>> refined_keys;
void refine_family(const vector<int>&groups,const vector<int>&coeff,int count){
 vector<int>key={active_row,count};key.insert(key.end(),groups.begin(),groups.end());key.insert(key.end(),coeff.begin(),coeff.end());
 if(!refined_keys.insert(key).second){refine_skipped++;return;}
 if(active_row==4||any_of(coeff.begin(),coeff.end(),[](int x){return x%(L/M)!=0;})){
  unconnected++;families<<tests<<" RETAINED "<<active_row<<' '<<count;for(int j=0;j<11;j++)families<<' '<<groups[j]<<' '<<coeff[j];families<<'\n';return;
 }
 auto oldblock=block,oldexponent=exponent,oldsizes=sizes;auto oldseen=move(seen);seen.clear();int oldrow=active_row;
 block=groups;for(int j=0;j<11;j++)exponent[j]=coeff[j]/(L/M);sizes.assign(count,1);active_row++;refined++;State t;cover(2047,t);
 active_row=oldrow;block=oldblock;exponent=oldexponent;sizes=oldsizes;seen=move(oldseen);
}
int main(int argc,char**argv){
 assert(argc==2);string path=argv[1];
 ifstream in(path+"/inputs.txt"),ff(path+"/families.txt");string line;map<long long,pair<vector<int>,vector<int>>> input;
 while(getline(in,line)){istringstream row(line);long long id;char pr;row>>id>>pr;vector<int>b(11),e(11);for(int&x:b)row>>x;for(int&x:e)row>>x;input[id]={b,e};}
 families.open(path+"/refine_unresolved.txt");points.open(path+"/refine_points.txt");long long records=0;
 while(getline(ff,line)){
  istringstream row(line);int c;string type;row>>tests>>c>>type;active_row=1;seen.clear();refined_keys.clear();
  assert(input.count(tests));auto original=input[tests].first;block_gauge_allowed=true;
  for(int j=0;j<5;j++)for(int k=5;k<11;k++)if(original[j]==original[k])block_gauge_allowed=false;
  if(type=="FIXED"){
   int q;row>>q;vector<int>groups(11),coeff(11);for(int j=0;j<11;j++)row>>groups[j]>>coeff[j];
   if(!row)continue;records++;refine_family(groups,coeff,c);
  }else if(type=="FREE"){
   assert(input.count(tests));auto data=input[tests];block=data.first;exponent=data.second;int k=*max_element(block.begin(),block.end())+1;sizes.assign(k,1);State st;
   for(int b=0;b<k;b++){int root,a,e;row>>root>>a>>e;st.parent[b]=root;st.a[b]=a;st.e[b]=e;}
   if(!row)continue;records++;leaf(st);
  }else {cerr<<"unknown "<<line<<endl;return 2;}
 }
 cout<<"records "<<records<<" refined "<<refined<<" skipped "<<refine_skipped<<" leaves "<<leaves<<" nodes "<<nodes<<" block_gauges "<<block_gauges<<" unresolved "<<unconnected<<" free "<<freezero<<" roots "<<rootcandidates<<" passed "<<rootpassed<<" actual "<<actual<<endl;
}
