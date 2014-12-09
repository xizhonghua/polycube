#include "util.h"
#include <fstream>
#include <iostream>

using namespace std;

int load_box_constraint(
    const char * box_file,
    std::vector<std::tuple<matrixd, boost::unordered_set<size_t> > > & box_vec)
{
  ifstream ifs(box_file);
  if(ifs.fail()){
    cerr << "# [error] can not open box_file" << endl;
    return __LINE__;
  }

  size_t box_num = 0;
  ifs >> box_num;
  box_vec.resize(box_num);
  size_t tet_num = 0;
  size_t tet_idx = -1;
  for(size_t bi = 0; bi < box_num; ++bi){
    get<0>(box_vec[bi]).resize(3,3);
    for(size_t ri = 0; ri < 9; ++ri)
      ifs >> get<0>(box_vec[bi])[ri];
    ifs >> tet_num;
    for(size_t ti = 0; ti < tet_num; ++ti){
      ifs >> tet_idx;
      get<1>(box_vec[bi]).insert(tet_idx);
    }
  }

  return 0;
}

int load_plane_constraint(
    const char * plane_file,
    std::vector<std::tuple<matrixd, boost::unordered_set<size_t> > > & plane_vec)
{
  ifstream ifs(plane_file);
  if(ifs.fail()){
    cerr << "# [error] can not open plane_file" << endl;
    return __LINE__;
  }

  size_t plane_num = 0;
  ifs >> plane_num;
  plane_vec.resize(plane_num);
  size_t tet_num = 0;
  size_t tet_idx = -1;
  for(size_t bi = 0; bi < plane_num; ++bi){
    get<0>(plane_vec[bi]).resize(3,1);
    for(size_t ri = 0; ri < 3; ++ri)
      ifs >> get<0>(plane_vec[bi])[ri];
    ifs >> tet_num;
    for(size_t ti = 0; ti < tet_num; ++ti){
      ifs >> tet_idx;
      get<1>(plane_vec[bi]).insert(tet_idx);
    }
  }

  return 0;
}

int load_line_constraint(
    const char * line_file,
    std::vector<std::tuple<matrixd, boost::unordered_set<size_t> > > & line_vec)
{
  ifstream ifs(line_file);
  if(ifs.fail()){
    cerr << "# [error] can not open line_file" << endl;
    return __LINE__;
  }

  size_t line_num = 0;
  ifs >> line_num;
  line_vec.resize(line_num);
  size_t tet_num = 0;
  size_t tet_idx = -1;
  for(size_t bi = 0; bi < line_num; ++bi){
    get<0>(line_vec[bi]).resize(3,1);
    for(size_t ri = 0; ri < 3; ++ri)
      ifs >> get<0>(line_vec[bi])[ri];
    ifs >> tet_num;
    for(size_t ti = 0; ti < tet_num; ++ti){
      ifs >> tet_idx;
      get<1>(line_vec[bi]).insert(tet_idx);
    }
  }

  return 0;
}

