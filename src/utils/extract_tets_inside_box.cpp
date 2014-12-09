#include <fstream>
#include "../common/def.h"
#include "../tetmesh/hex_io.h"
#include "../common/IO.h"

using namespace std;
using namespace zjucad::matrix;
int read_box(
    const char * box_file,
    matrixd & box)
{
  ifstream ifs(box_file);
  if(ifs.fail())
    return __LINE__;

  size_t col,row;
  ifs >> row >> col;
  box.resize(row, col);

  for(size_t t = 0; t < row * col; ++t){
    ifs >> box[t];
  }
  return 0;
}

int extract_tets_inside_box(int argc, char * argv[])
{
  if(argc < 4){
    cerr << "# [usage] extract_tets_inside_box tet frame box" << endl;
    return __LINE__;
  }

  matrixst tet;
  matrixd node;
  matrixd zyz;

  if(jtf::mesh::tet_mesh_read_from_zjumat(argv[1], &node, &tet)){
    cerr << "# [error] can not read tet mesh." << endl;
    return __LINE__;
  }

  jtf::mesh::read_matrix(argv[2], zyz);
  if(zyz.size(2) != tet.size(2)){
    cerr << "# [error] wrong zyz files." << endl;
    return __LINE__;
  }

  matrixd box;
  if(read_box(argv[3], box)){
    cerr << "# [error] can not read box file" << endl;
    return __LINE__;
  }

  matrixd max = zeros<double>(3,1),min = zeros<double>(3,1);
  for(size_t pi = 0; pi < box.size(1); ++pi){
    max[pi] = *max_element(box(pi,colon()).begin(),box(pi, colon()).end());
    min[pi] = *min_element(box(pi,colon()).begin(),box(pi, colon()).end());
  }

  set<size_t> inside_node;
  for(size_t pi = 0; pi < node.size(2); ++pi){
    size_t ai = 0;
    for(; ai < node.size(1); ++ai){
      if(node(ai, pi) > min[ai] && node(ai, pi) < max[ai]) continue;
      else
        break;
    }
    if(ai == node.size(1))
      inside_node.insert(pi);
  }

  vector<size_t> new_tet_idx;
  for(size_t ti = 0; ti < tet.size(2); ++ti){
    size_t inside_num = 0;
    for(size_t pi = 0; pi < tet.size(1); ++pi){
      if(find(inside_node.begin(), inside_node.end(),tet(pi,ti)) == inside_node.end())
        continue;
      ++inside_num;
    }
    if(inside_num == 4)
      new_tet_idx.push_back(ti);
  }

  // output new tet and new frame
  matrixst new_tet(4, new_tet_idx.size());
  matrixd new_zyz(3, new_tet_idx.size());
  for(size_t ti = 0; ti < new_tet_idx.size(); ++ti){
    new_tet(colon(),ti) = tet(colon(), new_tet_idx[ti]);
    new_zyz(colon(),ti) = zyz(colon(), new_tet_idx[ti]);
  }

  jtf::mesh::tet_mesh_write_to_zjumat("new_tet.tet", &node, &new_tet);
  jtf::mesh::write_matrix("new_zyz.zyz", new_zyz);
  return 0;
}
