#include <fstream>
#include <iostream>
#include <vector>
#include <list>
#include "../hexmesh/io.h"

using namespace std;
using namespace zjucad::matrix;
using namespace jtf::hexmesh;

int merge_vertex_pair(const char* hex_name,
                      const char * merge_vertex_file,
                      const char * dump_out_hex_file,
                      size_t hex_format)
{
  matrixd node;
  matrixst hex;
  if(hex_mesh_read_from_wyz(hex_name, hex, node, hex_format))
    return __LINE__;

  list<pair<size_t,size_t> > merge_vertex;
  ifstream ifs(merge_vertex_file);
  if(ifs.fail()){
    cerr << "# [error] can not open merge_vertex file." << endl;
    return __LINE__;
  }
  size_t vertex_pair_num = 0;
  ifs >> vertex_pair_num;
  assert(vertex_pair_num > 0);
  size_t vertex_from = 0, vertex_to = 0;

  for(size_t t = 0; t < vertex_pair_num; ++t){
    ifs >> vertex_from >> vertex_to;
    merge_vertex.push_back(make_pair(vertex_from,vertex_to));
  }

  // remove link in merge_vertex
  vector<pair<size_t,size_t> > merge_vertex_;
  while(!merge_vertex.empty()){
    vector<pair<size_t,size_t> > chain;
    chain.push_back(merge_vertex.front());
    merge_vertex.pop_front();

    for(list<pair<size_t,size_t> >::iterator lit = merge_vertex.begin();
        lit != merge_vertex.end();){
      if(lit->first == chain.back().second){
        chain.push_back(*lit);
        merge_vertex.erase(lit++);
      }else
        ++lit;
    }
    if(chain.front().first == chain.back().second) // it's a loop
      chain.pop_back();

    for(size_t t = 0; t < chain.size(); ++t){
      chain[t].second = chain.back().second;
    }
    merge_vertex_.insert(merge_vertex_.end(),chain.begin(),chain.end());
  }

  for(size_t t = 0; t < hex.size(2); ++t){

    for(size_t i = 0; i < merge_vertex_.size(); ++i){
      const size_t & from = merge_vertex_[i].first;
      const size_t & to = merge_vertex_[i].second;
      // if this hex contain "from" and no "to"
      if(find(hex(colon(),t).begin(),hex(colon(),t).end(),from)
         != hex(colon(),t).end()){

        if(find(hex(colon(),t).begin(),hex(colon(),t).end(),to)
           == hex(colon(),t).end()){
          //          if(new_hex_idx.back() != t)
          //            new_hex_idx.push_back(t);
         // cerr << "# [info] merge a vertex. " << endl;
          for(size_t j = 0; j < hex.size(1); ++j){
            if(hex(j,t)  == from) hex(j,t) = to;

          }
        }
      }
    }
  }

  if(dump_out_hex_file){
    ofstream ofs(dump_out_hex_file);
    if(ofs.fail()){
      cerr << "# [error] can not open nex hex file." << endl;
      return __LINE__;
    }
    ofs << "vertex_num " << node.size(2) << endl;
    ofs << "hex_num " << hex.size(2) << endl;
    for(size_t t = 0; t < node.size(2); ++t){
      ofs << node(0,t) << " " << node(1,t) << " " << node(2,t) << endl;
    }
    for(size_t t = 0; t < hex.size(2); ++t){
      for(size_t i = 0; i < hex.size(1); ++i){
        ofs << hex(i,t) << " ";
      }
      ofs << endl;
    }
  }
  return 0;
}


int merge_hex_vertex(int argc, char *argv[])
{
  if(argc != 5){
    cerr << argc << endl;
    cerr << "Useage: merge_hex_vertex hex format[default 1] input_hex input_merge_vertex_file output_hex" << endl;
    return __LINE__;
  }

  const size_t hex_format = atoi(argv[1]);

  if(merge_vertex_pair(argv[2],argv[3],argv[4],hex_format))
    return __LINE__;
  return 0;
}
