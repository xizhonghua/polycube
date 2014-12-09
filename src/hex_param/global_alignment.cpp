#include "global_alignment.h"
#include "hex_param.h"
#include "../common/graph.h"
#include "../common/zyz.h"
#include "../spherical_harmonics/rot_cubic_f_SH.h"
#include "../common/transition.h"

#include <iostream>

using namespace zjucad::matrix;
using namespace std;

void hj_frame_alignemt(const matrixst & tet,
                       const jtf::mesh::face2tet_adjacent &fa,
                       const matrixd &tet_zyz,
                       matrixd &new_frame_inner)
{
  matrix<matrixd > frame_inner(tet_zyz.size(2),1);

  fix_graph fg;
  tet2dual_graph(tet,fg,fa);

  matrixd tet_cf_sh(9, tet_zyz.size(2));
  for(size_t fi = 0; fi < tet_cf_sh.size(2); ++fi)
    calc_rot_cubic_f_sh_(&tet_cf_sh(0, fi), &tet_zyz(0, fi));

  std::vector<double> err(fg.edge_num());
  for(size_t ni = 0; ni < fg.node_num(); ++ni) {
    for(size_t ei = fg.ptr_[ni]; ei < fg.ptr_[ni+1]; ++ei)
      err[ei] = norm(tet_cf_sh(colon(), ni) - tet_cf_sh(colon(), fg.idx_[ei]));
  }

  std::vector<size_t> spanning_tree;
  fix_graph spanning_tree_graph;
  build_minimum_spanning_tree(fg, err, spanning_tree);
  subgraph(fg, spanning_tree, spanning_tree_graph);

  matrixd orient_err(spanning_tree_graph.edge_num());
  orient_frame(tet_zyz,frame_inner,spanning_tree_graph,orient_err);

  new_frame_inner.resize(3,tet.size(2));
  for(size_t t = 0; t < tet.size(2); ++t)
    rotation_matrix_2_zyz_angle(&frame_inner[t][0], &new_frame_inner(0,t),0);
}
