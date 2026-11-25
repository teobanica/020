// Refine the verified group-21 pair records against three two-parameter third rows.
#define NO_MAIN
#include "lattice.cpp"
int main(int argc,char**argv){
    assert(argc==4);
    ifstream xf(argv[1]),tf(argv[2]);ofstream out(argv[3]);
    assert(xf&&tf&&out);
    vector<Family> rows;char profile;int index,count;Family f{};
    while(xf>>profile>>index>>count){
        for(int&g:f.b)xf>>g;
        for(auto&v:f.x){v={};xf>>v[3];}
        assert(xf&&count<=3);rows.push_back(f);
    }
    assert(rows.size()==339008);
    int id;
    while(tf>>id){
        for(auto&v:old[1]){v={};tf>>v[3];v[3]*=21;}
        for(auto&v:old[1]){int g;tf>>g;if(g)v[g-1]=1;}
        assert(tf);finals.clear();
        for(const Family& row:rows) {State s;cover(2047,s,row,1);}
        for(const auto& rec:finals){out<<id;for(auto x:rec)out<<' '<<x;out<<'\n';}
        cout<<id<<" records "<<finals.size()<<endl;
    }
    out.flush();assert(out);out.close();assert(out);
}
