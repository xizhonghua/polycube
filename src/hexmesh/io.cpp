#include "io.h"
#include <string>
#include <iostream>
#include <fstream>
#include <vector>
#include <set>
#include <sstream>
#include <zjucad/matrix/io.h>

#include <boost/unordered_set.hpp>
using namespace std;
using namespace zjucad::matrix;

namespace jtf{
namespace hexmesh{

int hex_mesh_read_from_wyz(const char *path,
                           matrixst &hex,
                           matrixd &node,
                           const size_t hex_format,
                           bool is_remove_invalid_hex)
{
  if(hex_format > 3 || hex_format == 0) {
    cerr << "# [error] error hex_format." << endl;
    return __LINE__;
  }
  string call_function_name = "hex_mesh_read_from_wyz_";
  ostringstream os ;
  os << hex_format;
  call_function_name += os.str();
  if(is_remove_invalid_hex)
    call_function_name += "_check_valid";

#define CALL_SUB_PROG(prog)						\
  int prog(const char *path,                    \
  matrixst &hex,          \
  matrixd &node);		\
  if(call_function_name == #prog)				\
  return prog(path, hex, node);

  CALL_SUB_PROG(hex_mesh_read_from_wyz_1);
  CALL_SUB_PROG(hex_mesh_read_from_wyz_2);
  CALL_SUB_PROG(hex_mesh_read_from_wyz_3);
  CALL_SUB_PROG(hex_mesh_read_from_wyz_1_check_valid);
  CALL_SUB_PROG(hex_mesh_read_from_wyz_2_check_valid);
  CALL_SUB_PROG(hex_mesh_read_from_wyz_3_check_valid);

  return __LINE__;
}

int hex_mesh_read_from_wyz_1(const char *path,
                             matrixst &hex,
                             matrixd &node)
{
  ifstream ifs(path);
  if(ifs.fail()){
    cerr << "# [error] can not open file." << endl;
    return __LINE__;
  }
  string temp;
  size_t vertex_num;
  size_t hex_num;
  ifs >> temp >> vertex_num;
  ifs >> temp >> hex_num;
  node.resize(3,vertex_num);
  hex.resize(8,hex_num);

  for(size_t t = 0; t < vertex_num * 3; ++t)
    ifs >> node[t];

  for(size_t t = 0; t < hex_num * 8; ++t)
    ifs >> hex[t];

  return 0;
}

int hex_mesh_read_from_wyz_1_check_valid(const char *path,
                                         matrixst &hex,
                                         matrixd &node)
{
  ifstream ifs(path);
  if(ifs.fail()){
    cerr << "# [error] can not open file." << endl;
    return __LINE__;
  }
  string temp;
  size_t vertex_num;
  size_t hex_num;
  ifs >> temp >> vertex_num;
  ifs >> temp >> hex_num;
  node.resize(3,vertex_num);
  //hex.resize(8,hex_num);
  vector<size_t> hex_vec;
  hex_vec.reserve(8 * hex_num);

  for(size_t t = 0; t < vertex_num ; ++t){
    for(size_t i = 0; i < 3; ++i)
      ifs >> node[t*3+i];
  }

  vector<size_t> one_hex(8);
  std::set<size_t> one_hex_set;
  boost::unordered_set<set<size_t> > hex_record;
  for(size_t t = 0; t < hex_num; ++t){
    for(size_t i = 0; i < 8; ++i){
      ifs >> one_hex[i];
    }
    {
      one_hex_set.clear();
      one_hex_set.insert(one_hex.begin(), one_hex.end());
      if(one_hex_set.size() == one_hex.size()){
        if(hex_record.find(one_hex_set) == hex_record.end()){
          hex_vec.insert(hex_vec.end(),one_hex.begin(),one_hex.end());
          hex_record.insert(one_hex_set);
        }
      }
    }
  }
  hex.resize(8,hex_vec.size()/8);
  copy(hex_vec.begin(),hex_vec.end(),hex.begin());
  return 0;
}

int hex_mesh_read_from_wyz_2_check_valid(const char *path,
                                         matrixst &hex,
                                         matrixd &node)
{
  ifstream ifs(path);
  if(ifs.fail()){
    cerr << "# [error] can not open file." << endl;
    return __LINE__;
  }
  string temp;
  size_t vertex_num;
  size_t hex_num;
  ifs >> temp >> vertex_num;
  ifs >> temp >> hex_num;
  node.resize(3,vertex_num);
  //hex.resize(8,hex_num);
  vector<size_t> hex_vec;
  hex_vec.reserve(8 * hex_num);
  double trash;

  for(size_t t = 0; t < vertex_num ; ++t){
    for(size_t i = 0; i < 3; ++i)
      ifs >> node[t*3+i];
    ifs >> trash >> trash >> trash;
  }

  vector<size_t> one_hex(8);
  std::set<size_t> one_hex_set;
  boost::unordered_set<set<size_t> > hex_record;
  for(size_t t = 0; t < hex_num; ++t){
    for(size_t i = 0; i < 8; ++i){
      ifs >> one_hex[i];
    }
    {
      one_hex_set.clear();
      one_hex_set.insert(one_hex.begin(), one_hex.end());
      if(one_hex_set.size() == one_hex.size()){
        if(hex_record.find(one_hex_set) == hex_record.end()){
          hex_vec.insert(hex_vec.end(),one_hex.begin(),one_hex.end());
          hex_record.insert(one_hex_set);
        }
      }
    }
  }
  hex.resize(8,hex_vec.size()/8);
  copy(hex_vec.begin(),hex_vec.end(),hex.begin());
  return 0;
}

int hex_mesh_read_from_wyz_2(const char *path,
                             matrixst &hex,
                             matrixd &node)
{
  ifstream ifs(path);
  if(ifs.fail()){
    cerr << "# [error] can not open file." << endl;
    return __LINE__;
  }
  string temp;
  size_t vertex_num;
  size_t hex_num;
  ifs >> temp >> vertex_num;
  ifs >> temp >> hex_num;
  node.resize(3,vertex_num);
  hex.resize(8,hex_num);
  double trash;

  for(size_t t = 0; t < vertex_num ; ++t)
  {
    for(size_t i = 0; i < 3; ++i)
      ifs >> node[t*3+i];
    ifs >> trash >> trash >> trash;
  }

  for(size_t t = 0; t < hex_num * 8; ++t)
    ifs >> hex[t];

  return 0;
}


int hex_mesh_read_from_wyz_3(const char *path,
                             matrixst &hex,
                             matrixd &node)
{
  ifstream ifs(path);
  if(ifs.fail()){
    cerr << "# [error] can not open file." << endl;
    return __LINE__;
  }
  string temp;
  size_t vertex_num;
  size_t hex_num;
  ifs >> temp >> vertex_num;
  ifs >> temp >> hex_num;
  node.resize(3,vertex_num);
  //hex.resize(8,hex_num);
  vector<size_t> hex_vec;
  hex_vec.reserve(8 * hex_num);
  double trash;

  for(size_t t = 0; t < vertex_num ; ++t)
  {
    for(size_t i = 0; i < 3; ++i)
      ifs >> node[t*3+i];
    ifs >> trash >> trash >> trash;
  }

  vector<size_t> one_hex(8);
  for(size_t t = 0; t < hex_num; ++t){
    for(size_t i = 0; i < 8; ++i){
      ifs >> one_hex[i];
    }
    if(find(one_hex.begin(),one_hex.end(),-1) == one_hex.end()){
      hex_vec.insert(hex_vec.end(),one_hex.begin(),one_hex.end());
    }
  }
  hex.resize(8,hex_vec.size()/8);
  copy(hex_vec.begin(),hex_vec.end(),hex.begin());
  return 0;
}

int hex_mesh_read_from_wyz_3_check_valid(const char *path,
                                         matrixst &hex,
                                         matrixd &node)
{
  ifstream ifs(path);
  if(ifs.fail()){
    cerr << "# [error] can not open file." << endl;
    return __LINE__;
  }
  string temp;
  size_t vertex_num;
  size_t hex_num;
  ifs >> temp >> vertex_num;
  ifs >> temp >> hex_num;
  node.resize(3,vertex_num);
  //hex.resize(8,hex_num);
  vector<size_t> hex_vec;
  hex_vec.reserve(8 * hex_num);
  double trash;

  for(size_t t = 0; t < vertex_num ; ++t){
    for(size_t i = 0; i < 3; ++i)
      ifs >> node[t*3+i];
    ifs >> trash >> trash >> trash;
  }

  vector<size_t> one_hex(8);
  set<size_t> one_hex_set;
  boost::unordered_set<set<size_t> > hex_record;
  for(size_t t = 0; t < hex_num; ++t){
    for(size_t i = 0; i < 8; ++i){
      ifs >> one_hex[i];
    }
    if(find(one_hex.begin(),one_hex.end(),-1) == one_hex.end()){
      one_hex_set.clear();
      one_hex_set.insert(one_hex.begin(), one_hex.end());
      if(one_hex_set.size() == one_hex.size()){
        // for
        if(find(one_hex_set.begin(), one_hex_set.end(),0) != one_hex_set.end())
          continue;

        if(hex_record.find(one_hex_set) == hex_record.end()){
          hex_record.insert(one_hex_set);
          hex_vec.insert(hex_vec.end(),one_hex.begin(),one_hex.end());
        }
      }
    }
  }
  hex.resize(8,hex_vec.size()/8);
  copy(hex_vec.begin(),hex_vec.end(),hex.begin());
  return 0;
}

int hex_mesh_write_to_wyz(const char *path,
                          matrixst &hex,
                          matrixd &node)
{
  ofstream ofs(path);
  if(ofs.fail()){
    cerr << "# [error] can not open file." << endl;
    return __LINE__;
  }

  ofs << "vertex_num " << node.size(2) << endl;
  ofs << "hex_num " << hex.size(2) << endl;

  for(size_t t = 0; t < node.size(2); ++t)  ofs << node(0,t) << " " << node(1,t) << " " << node(2,t) << endl;

  for(size_t t = 0; t < hex.size(2); ++t){
    for(size_t i = 0; i < hex.size(1); ++i){
      ofs << hex(i,t) << " ";
    }
    ofs << endl;
  }
  return 0;
}

int hex_mesh_read_from_vtk(
    const char *path,
    matrixd *node,
    matrixst *hex)
{
  ifstream ifs(path);
  if(ifs.fail()) {
    cerr << "[info] " << "can not open file" << path << endl;
    return __LINE__;
  }

  matrixd node0;
  matrix<int> hex1;

  string str;
  int point_num = 0,cell_num = 0;

  while(!ifs.eof()){
    ifs >> str;
    if(str == "POINTS"){
      ifs >> point_num >> str;
      node0.resize(3, point_num);
      for(size_t i = 0;i < point_num; ++i){
        for(size_t j = 0;j < 3; ++j)
          ifs >> node0(j, i);
      }
      continue;
    }
    if(str == "CELLS"){
      ifs >> cell_num >> str;
      int point_number_of_cell = 0;
      vector<size_t> hex_temp;
      for(size_t ci = 0; ci < cell_num; ++ci){
        ifs >> point_number_of_cell;
        if(point_number_of_cell != 8){
          for(size_t i = 0; i < point_number_of_cell; ++i)
            ifs >> str;
        }else{
          int p;
          for(size_t i = 0; i < point_number_of_cell; ++i){
            ifs >> p;
            hex_temp.push_back(p);
          }
        }
      }
      hex1.resize(8, hex_temp.size()/8);
      copy(hex_temp.begin(), hex_temp.end(), hex1.begin());
    }
  }

  vector<size_t> one_hex(hex1.size(1));
  for(size_t hi = 0; hi < hex1.size(2); ++hi){
      copy(hex1(colon(),hi).begin(), hex1(colon(), hi).end(), one_hex.begin());
      hex1(0,hi) = one_hex[6];
      hex1(1,hi) = one_hex[5];
      hex1(2,hi) = one_hex[7];
      hex1(3,hi) = one_hex[4];
      hex1(4,hi) = one_hex[2];
      hex1(5,hi) = one_hex[1];
      hex1(6,hi) = one_hex[3];
      hex1(7,hi) = one_hex[0];
    }

  *node = node0;
  *hex = hex1;
  return 0;
}
}
}
