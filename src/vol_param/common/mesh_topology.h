#ifndef MESH_TOPOLOGY_H
#define MESH_TOPOLOGY_H

#include <zjucad/ptree/ptree.h>
#include <zjucad/matrix/matrix.h>
#include <jtflib/mesh/mesh.h>
#include <vector>

template <typename T1, typename T2>
void construct_original_face_in_cut_mesh(
    const zjucad::matrix::matrix<size_t> & uncut_tet,
    const zjucad::matrix::matrix<double> & uncut_node,
    const jtf::mesh::face2tet_adjacent & uncut_fa,
    const zjucad::matrix::matrix_expression<T1> & cut_tet,
    const jtf::mesh::face2tet_adjacent & cut_fa,
    zjucad::matrix::matrix_expression<T2> &orig_face_in_cut,
    const zjucad::matrix::matrix<size_t> & cut_tet2tet)
{
  using namespace std;
  using namespace zjucad::matrix;

  matrix<size_t> outside_face_cut;
  jtf::mesh::get_outside_face(cut_fa, outside_face_cut, true);


  vector<size_t> orig_face_idx_in_cut; // it store original surface in outside_face_cut

  for(size_t fi = 0; fi < outside_face_cut.size(2); ++fi){
      const size_t face_idx  = uncut_fa.get_face_idx(
            cut_tet2tet[outside_face_cut(0,fi)],
          cut_tet2tet[outside_face_cut(1,fi)],
          cut_tet2tet[outside_face_cut(2,fi)]);
      assert(face_idx != -1);
      if(uncut_fa.is_outside_face(uncut_fa.face2tet_[face_idx]))
        orig_face_idx_in_cut.push_back(fi);
    }

  orig_face_in_cut().resize(3, orig_face_idx_in_cut.size());
  for(size_t i = 0; i < orig_face_idx_in_cut.size(); ++i){
      orig_face_in_cut()(colon(),i) =
          outside_face_cut(colon(),orig_face_idx_in_cut[i]);
    }
}

#endif // MESH_TOPOLOGY_H
