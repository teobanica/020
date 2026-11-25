#define main curve_main
#include "curve_graphs.cpp"
#undef main
int main(int argc,char**argv){assert(argc==3);tables();ifstream in(argv[1]);ofstream out(argv[2]);vector<Family>rows;vector<pair<int,int>>chars;int x;
 while(in>>x){array<int,67>z;z[0]=x;for(int k=1;k<67;k++)in>>z[k];assert(in);for(int k=0;k<12;k++)assert(z[k]==0);Family r;r.k.fill(0);pair<int,int>ch={z[13+35],z[14+35]};assert(ch==make_pair(0,0)||ch==make_pair(1,0)||ch==make_pair(0,1));
  for(int j=0;j<11;j++){assert(z[12+5*j]==0&&z[15+5*j]==0);r.e[j]=z[16+5*j];assert(r.e[j]%14==0);assert(z[13+5*j]==(j<7?0:ch.first)&&z[14+5*j]==(j<7?0:ch.second));}rows.push_back(r);chars.push_back(ch);}
 int n=rows.size(),outside=0;vector<vector<bool>>adj(n,vector<bool>(n));vector<pair<int,int>>edges;
 for(int i=0;i<n;i++)for(int j=0;j<i;j++){auto a=rows[i],b=rows[j];bool variable=chars[i]!=chars[j];if(variable)for(int k=7;k<11;k++)a.k[k]=1;bool generic=possible(a,b,-1);if(generic){adj[i][j]=adj[j][i]=true;edges.push_back({j,i});}
  vector<int>roots;if(!generic&&variable)for(int t=0;t<840;t++)if(possible(a,b,t)){roots.push_back(t);if(t%14)outside++;}
  out<<j<<' '<<i<<' '<<generic<<' '<<chars[i].first-chars[j].first<<' '<<chars[i].second-chars[j].second<<' '<<roots.size();for(int t:roots)out<<' '<<t;out<<'\n';
 }
 int triangles=0;for(int i=0;i<n;i++)for(int j=0;j<i;j++)for(int k=0;k<j;k++)if(adj[i][j]&&adj[i][k]&&adj[j][k])triangles++;
 out.flush();assert(out);cout<<"rows "<<n<<" generic_edges "<<edges.size()<<" generic_triangles "<<triangles<<" exceptional_roots_outside_mu60 "<<outside<<endl;
}
