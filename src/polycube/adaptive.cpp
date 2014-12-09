#include "adaptive.h"

#include <map>
#include <set>
#include <vector>

#include <zjucad/matrix/io.h>

#include "../tetmesh/tetmesh.h"
#include "../tetmesh/subdivide_tet.h"

using namespace std;
using namespace zjucad::matrix;

int boundary_uniform_subdivide(const zjucad::matrix::matrix<double> &parent_node,
                               const zjucad::matrix::matrix<size_t> &parent_tet,
                               zjucad::matrix::matrix<double> &child_node,
                               zjucad::matrix::matrix<size_t> &child_tet,
                               zjucad::matrix::matrix<size_t> &new_node_parent_id //nx2 to ori node
                               )
{
  unique_ptr<jtf::mesh::face2tet_adjacent> fa(jtf::mesh::face2tet_adjacent::create(parent_tet));
  if(!fa.get()){
    cerr  << "# [error] can not buildjtf::mesh::face2tet_adjacent." << endl;
    return __LINE__;
  }

  matrix<size_t> faces;
  get_outside_face(*fa, faces);

  const static size_t face2edge[] = {
    0, 1, 0, 2, 1, 2
  };
  matrix<size_t> &split_edge = new_node_parent_id;
  create_ordered_table(faces, itr_matrix<const size_t *>(2, 3, face2edge), new_node_parent_id);

  time_t beg = clock();
  matrix<size_t> before_validate = new_node_parent_id;
  validate_edge_nodes(parent_tet, before_validate, new_node_parent_id);  
  cout << "validate: " << (clock()-beg)/double(CLOCKS_PER_SEC) << endl;
  cerr << "# validate subdivide: " << before_validate.size(1)
       << " " << new_node_parent_id.size(1) << endl;

  beg = clock();
  subdivide_tetmesh(parent_tet, parent_node.size(2), new_node_parent_id, child_tet);
  cout << "subdivide: " << (clock()-beg)/double(CLOCKS_PER_SEC) << endl;

  subdivide_top2geo(parent_node, new_node_parent_id, child_node);

  cerr << "sub tet size: " << child_node.size(2) << " " << child_tet.size(2) << endl;
  return 0;
}
