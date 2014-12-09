#include "tetmesh.h"

#include <fstream>
#include <iostream>
#include <memory>
#include <stack>

#include <zjucad/matrix/io.h>
#include <hjlib/sparse/sparse.h>
#include <jtflib/mesh/mesh.h>
#include "../common/IO.h"

using namespace std;
using namespace zjucad::matrix;
using namespace hj::sparse;

int tet2tet_adj_matrix(
    const matrixst &tet,
    matrixst &ptr,
    matrixst &idx,
    const jtf::mesh::face2tet_adjacent *fa,
    int type
    )
{
  if(!fa) {
    cerr << "needjtf::mesh::face2tet_adjacent in the current implementation." << endl;
    return __LINE__;
  }
  const size_t tet_num = tet.size(2);
  csc_by_vm<bool, size_t, map_by_sorted_vector> csc(tet_num, tet_num);
  matrixst tet_i(4);
  for(size_t ti = 0; ti < tet_num; ++ti) { // for each tet
    tet_i = tet(colon(), ti);
    for(size_t fi = 0; fi < 4; ++fi) { // for each face
      rotate(tet_i.begin(), tet_i.begin()+1, tet_i.end());
      const pair<size_t, size_t> nbt = fa->query(&tet_i[0]);
      const size_t other = (nbt.first == ti)?nbt.second:nbt.first;
      if(other == -1) continue; // surface triangle face
      if(type == 1 && other < ti) continue;
      if(type ==-1 && other > ti) continue;
      csc[ti][other] = true;
    }
  }
  hj::sparse::csc<bool> rtn;
  convert(csc, rtn);
  ptr = rtn.ptr();
  idx = rtn.idx();
  return 0;
}

int face2face_adj_matrix(
    const matrixst &tet,
    matrixst &ptr,
    matrixst &idx,
    const jtf::mesh::face2tet_adjacent *fa,
    int type
    )
{
  if(!fa) {
    cerr << "needjtf::mesh::face2tet_adjacent in the current implementation." << endl;
    return __LINE__;
  }
  const size_t tet_num = tet.size(2), face_num = fa->faces_.size();
  csc_by_vm<bool, size_t, map_by_sorted_vector> csc(face_num, face_num);
  matrixst tet_i(4);
  size_t face_idx[4];
  for(size_t ti = 0; ti < tet_num; ++ti) { // for each tet
    tet_i = tet(colon(), ti);
    for(int fi = 0; fi < 4; ++fi) {// for each face
      face_idx[fi] = fa->get_face_idx(tet_i[fi], tet_i[(fi+1)%4], tet_i[(fi+2)%4]);
      if(face_idx[fi] >= fa->faces_.size())
        cerr << "strange and wrong face." << endl;
    }
    for(int fi = 0; fi < 4; ++fi) {
      for(int fj = 0; fj < 4; ++fj) {
        if(fi == fj) continue;
        if(type == 1 && face_idx[fi] > face_idx[fj]) continue;
        if(type == -1 && face_idx[fi] < face_idx[fj]) continue;
        csc[face_idx[fi]][face_idx[fj]] = true;
      }
    }
  }
  hj::sparse::csc<bool> rtn;
  convert(csc, rtn);
  ptr = rtn.ptr();
  idx = rtn.idx();
  return 0;
}
