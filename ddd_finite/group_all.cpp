#include <array>
#include <cassert>
#include <string>
#include <algorithm>
#include <fstream>
#include <iostream>
#include <map>
#include <vector>
using namespace std;int mod(int a){a%=40;return a<0?a+40:a;}using Row=array<int,11>;using Key=array<int,22>;
int main(int argc,char**argv){assert(argc==2);string root=argv[1];ifstream in(root+"/zero_seeds.txt");ofstream out(root+"/zero_groups.txt"),mapout(root+"/zero_index.txt");assert(in&&out&&mapout);map<Row,vector<pair<int,Row>>>groups;int id=0;Key raw;
 while(in>>raw[0]){for(int j=1;j<22;j++)in>>raw[j];int n;in>>n;vector<int>classes(n);for(int&x:classes)in>>x;array<Row,3>e{};for(int j=0;j<11;j++){e[1][j]=raw[2*j];e[2][j]=raw[2*j+1];}
  Key best;best.fill(100);for(int a=0;a<3;a++)for(int b=0;b<3;b++)if(a!=b){int c=3-a-b;for(int base=0;base<11;base++){array<pair<int,int>,11>cols;for(int j=0;j<11;j++)cols[j]={mod(e[b][j]-e[a][j]-e[b][base]+e[a][base]),mod(e[c][j]-e[a][j]-e[c][base]+e[a][base])};sort(cols.begin(),cols.end());Key key;for(int j=0;j<11;j++){key[j]=cols[j].first;key[11+j]=cols[j].second;}best=min(best,key);}}
  Row u,v;copy(best.begin(),best.begin()+11,u.begin());copy(best.begin()+11,best.end(),v.begin());groups[u].push_back({id,v});mapout<<id<<' '<<n;for(int c:classes)mapout<<' '<<c;for(int x:best)mapout<<' '<<x;mapout<<'\n';id++;
 }
 out<<groups.size()<<'\n';int group=0;for(auto&[u,vs]:groups){out<<group++<<' '<<vs.size();for(int x:u)out<<' '<<x;out<<'\n';for(auto&[id,v]:vs){out<<id;for(int x:v)out<<' '<<x;out<<'\n';}}
 cout<<"groups "<<groups.size()<<" triples "<<id<<endl;
}
