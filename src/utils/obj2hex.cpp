#include <sstream>
#include<iostream>
#include <cassert>
#include<fstream>
#include <iomanip>
#include<string>
#include<vector>
#include<map>
#include<set>
#include <utility>
#include<algorithm>
#include <boost/property_tree/ptree.hpp>
#include <jtflib/mesh/io.h>

using namespace std;
using namespace zjucad::matrix;

int n=0;
vector<pair<int,int> > preedge;
map<pair<int,int>,vector<pair<int,int> > >::iterator coniter;
vector<int> ver_store;

void build_edge(zjucad::matrix::matrix<size_t> &fc,
                map<pair<int,int>,vector<pair<int,int> > > & o_edge)
{
  pair<int,int> one_edge, oppo_edge;
  for(size_t fi = 0; fi < fc.size(2); ++fi){
      for(size_t pi = 0; pi < 2; ++pi){
          one_edge.first = fc(pi,fi);
          one_edge.second = fc((pi+1)%4,fi);
          if(one_edge.first > one_edge.second) swap(one_edge.first, one_edge.second);

          oppo_edge.first = fc((pi+2)%4,fi);
          oppo_edge.second = fc((pi+3)%4,fi);
          if(oppo_edge.first > oppo_edge.second) swap(oppo_edge.first, oppo_edge.second);

          o_edge[one_edge].push_back(oppo_edge);
          o_edge[oppo_edge].push_back(one_edge);
        }
    }
}

void find_cube(map<pair<int,int>,vector<pair<int,int> > >::iterator iter,
               map<pair<int,int>,vector<pair<int,int> > > & t_edge,
               vector<vector<int> > &vcube,
               zjucad::matrix::matrix<double> &vtable,
               set<set<int> > &cmp_ver){

  map<pair<int,int>,vector<pair<int,int> > >::iterator niter;
  vector<pair<int,int> >::iterator viter;
  set<int> ver_set;
  preedge.push_back(iter->first);
  matrix<double> p(3,1),p0(3,1);
  if(n==3){
      int flagh=0;
      vector<pair<int,int> >::iterator iterh;
      for(iterh=iter->second.begin();!flagh&&iterh!=iter->second.end();iterh++)
        if((*iterh)==coniter->first)
          flagh=1;
      if(flagh){
          vector<pair<int,int> >::iterator iter7;
          ver_set.clear();
          for(iter7=preedge.begin();iter7!=preedge.end();iter7++){
              ver_set.insert(iter7->first);
              ver_set.insert(iter7->second);
            }
          int presum=cmp_ver.size();
          cmp_ver.insert(ver_set);
          if(cmp_ver.size()>(size_t)presum){
              vector<pair<int,int> >::iterator iterpre=preedge.begin();
              vector<pair<int,int> >::iterator iternext=preedge.begin();
              vector<int> s_store;
              ++iternext;
              int a=iterpre->first;
              int b=iterpre->second;

              p0  =vtable(colon(),a) - vtable(colon(),b);

              s_store.clear();
              s_store.push_back(a);
              s_store.push_back(b);
              for(;iternext!=preedge.end();iternext++){
                  int c=iternext->first;
                  int d=iternext->second;

                  p = vtable(colon(),c) -vtable(colon(),d);
                  if(dot(p,p0)<0){
                      s_store.push_back(d);
                      s_store.push_back(c);
                    }
                  else{
                      s_store.push_back(c);
                      s_store.push_back(d);
                    }
                }
              a=s_store[6];
              b=s_store[7];
              s_store[6]=s_store[4];
              s_store[7]=s_store[5];
              s_store[4]=a;
              s_store[5]=b;

              vcube.push_back(s_store);//put the vertex into vector
            }
        }
      n--;
      preedge.pop_back();
    }
  else{
      for(viter=iter->second.begin();viter!=iter->second.end();viter++){
          int flag1=1;
          vector<pair<int,int> >::iterator iter6;
          for(iter6=preedge.begin();flag1&&iter6!=preedge.end();iter6++)
            if((*iter6)==(*viter))
              flag1=0;
          if(flag1){
              niter=t_edge.find(*viter);
              if(niter!=t_edge.end()){
                  n++;
                  find_cube(niter,t_edge,vcube,vtable,cmp_ver);
                }
            }
        }
      n--;
      preedge.pop_back();
    }
}
void record_cube(vector<vector<int> > &cube_v, zjucad::matrix::matrix<double> &v,
                 const std::string& out_hex_file){
  ofstream out_file(out_hex_file.c_str());
  if (out_file.fail()) {
      std::cerr << "# can not open file: " << out_file << std::endl;
      return;
    }
  vector<vector<int> >::iterator iter;
  vector<int>::iterator iterv;
  cout <<"Size:"<<cube_v.size()<<endl;
  out_file << "vertex " << v.size(2) << endl;
  out_file << "hex " << cube_v.size() << endl;
  for(size_t i = 0; i < v.size(2); ++i){
      for(size_t di = 0; di < v.size(1); ++di)
        out_file << v(di,i) << " ";
      out_file << endl;
    }
  for(iter=cube_v.begin();iter!=cube_v.end();iter++){
      for(iterv=iter->begin();iterv!=iter->end();iterv++)
        out_file <<*iterv<<" ";
      out_file <<endl;
    }
  out_file.close();
}

int obj2hex(int argc, char * argv[])
{
  if(argc != 3){
      cerr << "# [usage] obj2hex in_obj out_hex" << endl;
      return __LINE__;
    }

  map<pair<int,int>,vector<pair<int,int> > >  opp_edge;
  map<pair<int,int>,vector<pair<int,int> > >::iterator itered;

  vector<vector<int> > cube_ver;
  set<set<int> > cmp_set;

  std::cout << "# [ Read obj file ]" << std::endl;
  matrix<size_t> face;
  matrix<double> vertex;

  if(jtf::mesh::load_obj(argv[1], face, vertex)){
      cerr << "# [error] can not open obj." << endl;
      return __LINE__;
    }
  std::cout << "# [ Build Edge ] Start..." << std::endl;
  build_edge(face,opp_edge);
  std::cout << "# [ Build Edge ] END " << std::endl;

  cout <<"vunm:"<<vertex.size(2)<<endl;
  cout <<"The size of oppsite-edge:"<<opp_edge.size()<<endl;
  for(itered=opp_edge.begin();itered!=opp_edge.end();++itered){
      coniter=itered;
      ver_store.clear();
      preedge.clear();
      n=0;
      find_cube(itered,opp_edge,cube_ver,vertex,cmp_set);
    }
  record_cube(cube_ver,vertex, argv[2]);
  return 0;
}
