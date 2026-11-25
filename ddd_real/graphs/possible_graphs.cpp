#define EXTEND_NO_MAIN
#include "extend.cpp"
#include <map>
struct Family{int c;Row groups,e;};
namespace compatible{
Row cells,coeff;array<int,4>cell_group;int count;long long nodes=0;
bool relation(const State&s){
 if(count<4)return true;array<int,5>charge{};int phase=0;
 for(int cell=0;cell<4;cell++){auto v=s.find(cell_group[cell]);int sign=(cell==0||cell==3)?1:-1;charge[v.first]+=sign;phase+=sign*v.second;}
 return any_of(charge.begin(),charge.end(),[](int x){return x!=0;})||mod(phase)==0;
}
bool cover(int mask,const State&s);
bool cycle(int mask,int anchor,int prime,int label,const State&s){nodes++;if(label==prime)return compatible::cover(mask,s);
 for(int j=0;j<11;j++)if((mask>>j)&1){State t=s;int v=mod(M/prime*label+coeff[anchor]-coeff[j]);if(t.link(cells[j],cells[anchor],v)&&compatible::cycle(mask^(1<<j),anchor,prime,label+1,t))return true;}
 return false;
}
bool cover(int mask,const State&s){if(!mask)return relation(s);int anchor=__builtin_ctz((unsigned)mask),n=__builtin_popcount((unsigned)mask);for(int p:{2,3,5,7})if(p<=n&&(n==p||n-p>=2))if(compatible::cycle(mask^(1<<anchor),anchor,p,1,s))return true;return false;}
#ifdef AUDIT
Row xgroups,ygroups;int witness_q=0,witness_s=0;
int mod2(int x){x%=2*M;return x<0?x+2*M:x;}
struct Congruences{
 int rank=0,a=0,b=0,v=0;vector<pair<int,int>>solutions;
 bool add(int aa,int bb,int vv){
  vv=mod(vv);if(!aa&&!bb)return vv==0;
  if(!rank){rank=1;a=aa;b=bb;v=vv;return true;}
  if(rank==1){int det=a*bb-b*aa;if(!det){int sign=(aa==a&&bb==b)?1:-1;assert(aa==sign*a&&bb==sign*b);return mod(vv-sign*v)==0;}
   rank=2;int d=abs(det),sgn=det>0?1:-1;assert(d<=2);
   set<pair<int,int>>list;for(int k=0;k<d;k++)for(int l=0;l<d;l++){
    int q=sgn*(bb*(v+M*k)-b*(vv+M*l))*(2/d);
    int z=sgn*(-aa*(v+M*k)+a*(vv+M*l))*(2/d);list.insert({mod2(q),mod2(z)});
   }solutions.assign(list.begin(),list.end());assert(int(solutions.size())==d);return true;
  }
  vector<pair<int,int>>keep;for(auto [q,z]:solutions)if(mod2(aa*q+bb*z-2*vv)==0)keep.push_back({q,z});solutions=move(keep);return !solutions.empty();
 }
};
bool direct_cover(int mask,const Congruences&s);
bool direct_cycle(int mask,int anchor,int prime,int label,const Congruences&s){nodes++;if(label==prime)return direct_cover(mask,s);
 for(int j=0;j<11;j++)if((mask>>j)&1){Congruences t=s;int a=xgroups[j]-xgroups[anchor],b=ygroups[anchor]-ygroups[j],v=mod(M/prime*label-coeff[j]+coeff[anchor]);if(t.add(a,b,v)&&direct_cycle(mask^(1<<j),anchor,prime,label+1,t))return true;}return false;
}
bool direct_cover(int mask,const Congruences&s){if(!mask){
 witness_q=witness_s=0;if(s.rank==1){if(s.a)witness_q=mod2(2*s.v*s.a);else witness_s=mod2(2*s.v*s.b);}else if(s.rank==2){witness_q=s.solutions[0].first;witness_s=s.solutions[0].second;}return true;}int anchor=31-__builtin_clz((unsigned)mask),n=__builtin_popcount((unsigned)mask);for(int p:{7,5,3,2})if(p<=n&&(n==p||n-p>=2))if(direct_cycle(mask^(1<<anchor),anchor,p,1,s))return true;return false;}
#endif
bool possible(const Family&a,const Family&b){
#ifdef AUDIT
 witness_q=witness_s=0;
#endif
 cell_group.fill(-1);count=0;for(int j=0;j<11;j++){int cell=2*a.groups[j]+b.groups[j];if(cell_group[cell]<0)cell_group[cell]=count++;cells[j]=cell_group[cell];coeff[j]=mod(a.e[j]-b.e[j]);}
 if(count==1){int sum=0;for(int e:coeff)sum+=positive[e][0];if(sum%P)return false;stamp++;return point_cover(coeff,2047);}
 if(count==2){array<int,2>v{},vc{};for(int j=0;j<11;j++){v[cells[j]]+=positive[coeff[j]][0];vc[cells[j]]+=negative[coeff[j]][0];}if((1LL*v[0]*vc[0]-1LL*v[1]*vc[1])%P)return false;}
#ifdef AUDIT
 xgroups=a.groups;ygroups=b.groups;Congruences constraints;return direct_cover(2047,constraints);
#else
 State s;return compatible::cover(2047,s);
#endif
}
}
vector<vector<bool>>adj;vector<int>best,current;
void expand(vector<int>candidates){
 if(current.size()+candidates.size()<=best.size())return;
 vector<int>order,colors,remaining=candidates;int color=0;
 while(!remaining.empty()){color++;vector<int>independent,next;for(int v:remaining){bool ok=true;for(int w:independent)if(adj[v][w]){ok=false;break;}if(ok){independent.push_back(v);order.push_back(v);colors.push_back(color);}else next.push_back(v);}remaining=move(next);}
 for(int i=int(order.size())-1;i>=0;i--){if(current.size()+colors[i]<=best.size())return;int v=order[i];current.push_back(v);vector<int>next;for(int j=0;j<i;j++)if(adj[v][order[j]])next.push_back(order[j]);if(next.empty()){if(current.size()>best.size())best=current;}else expand(next);current.pop_back();}
}
int main(int argc,char**argv){assert(argc==3);tables();ifstream in(argv[1]);string temporary=string(argv[2])+".writing";ofstream out(temporary);assert(in&&out);
#ifdef AUDIT
 string wfinal=string(argv[2])+".witnesses",wtemporary=wfinal+".writing";ofstream witnesses(wtemporary);assert(witnesses);
#endif
 map<int,vector<Family>>data;int id;Family family;while(in>>id>>family.c){for(int&x:family.groups)in>>x;for(int&x:family.e)in>>x;assert(in);data[id].push_back(family);}
 int finite=0,parametric=0,skipped=0;map<int,int>hist;
 for(auto&[id,rows]:data){if(any_of(rows.begin(),rows.end(),[](auto&r){return r.c>2;})){skipped++;out<<id<<" SKIPPED "<<rows.size()<<'\n';continue;}bool param=any_of(rows.begin(),rows.end(),[](auto&r){return r.c==2;});if(param)parametric++;else finite++;
  int n=rows.size();adj.assign(n,vector<bool>(n));vector<pair<int,int>>edges;
  for(int i=0;i<n;i++)for(int j=0;j<i;j++)if(compatible::possible(rows[i],rows[j])){adj[i][j]=adj[j][i]=true;
#ifdef AUDIT
 if(edges.empty())witnesses<<id<<' '<<j<<' '<<i<<' '<<compatible::witness_s<<' '<<compatible::witness_q<<'\n';
#endif
 edges.push_back({j,i});}
  best.clear();current.clear();vector<int>candidates(n);iota(candidates.begin(),candidates.end(),0);expand(candidates);hist[best.size()]++;
  out<<id<<(param?" PARAMETRIC ":" FINITE ")<<n<<' '<<edges.size()<<' '<<best.size();for(int v:best)out<<' '<<v;for(auto [a,b]:edges)out<<' '<<a<<' '<<b;out<<'\n';
 }
 out.flush();assert(out);out.close();assert(!out.fail());assert(rename(temporary.c_str(),argv[2])==0);
#ifdef AUDIT
 witnesses.flush();assert(witnesses);witnesses.close();assert(!witnesses.fail());assert(rename(wtemporary.c_str(),wfinal.c_str())==0);
#endif
 cout<<"finite "<<finite<<" parametric "<<parametric<<" skipped "<<skipped<<" nodes "<<compatible::nodes<<" max_clique_counts";for(auto [a,b]:hist)cout<<' '<<a<<':'<<b;cout<<endl;
}
