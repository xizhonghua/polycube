#include <fstream>
#include <iostream>
#include <vector>
#include <algorithm>
#include <memory>

#include "../common/zyz.h"
#include "../common/IO.h"
#include "../numeric/util.h"
#include "../tetmesh/tetmesh.h"

#include <zjucad/matrix/io.h>

using namespace std;
using namespace zjucad::matrix;

int init_prism(int argc, char *argv[])
{
  if(argc < 3) {
      cerr << "tet2vtk tet output" << endl;
      return __LINE__;
    }

  jtf::mesh::meshes tm;
  if(jtf::mesh::tet_mesh_read_from_zjumat(
       argv[1],
       &tm.node_, &tm.mesh_, 0))
    return __LINE__;

  // triangle bottom is at x-y plain, z is from 0 to 0.5
  matrixd xy_corners = zeros<double>(2, 3);
  const double tan_30 = tan(My_PI()/6);
  cerr << "# tan_30: " << tan_30 << endl;
  xy_corners(0, 0) = -tan_30+0.2;
  xy_corners(1, 0) = 1;
  xy_corners(0, 1) = tan_30;
  xy_corners(1, 1) = 1+0.2;

  unique_ptr<jtf::mesh::face2tet_adjacent> fa(jtf::mesh::face2tet_adjacent::create(tm.mesh_));
  matrixst node_idx(3);
  matrixd dist(3), zyz(3, fa->faces_.size()), zyz_frame = zeros<double>(3, 3);
  zyz_frame(2, 0) = (90+15)/180.0*My_PI();
  zyz_frame(2, 1) = (-30+15)/180.0*My_PI();
  zyz_frame(2, 2) = (-150+15)/180.0*My_PI();
  cerr << xy_corners << zyz_frame << endl;
  for(size_t i = 0; i < fa->faces_.size(); ++i) {
      copy(fa->faces_[i].begin(), fa->faces_[i].end(), node_idx.begin());
      const matrixd pos = tm.node_(colon(), node_idx)*ones<double>(3, 1)/3.0;
      for(size_t ci = 0; ci < 3; ++ci)
        dist[ci] = norm(pos(colon(0, 1))-xy_corners(colon(), ci));
      const size_t voronoi_id = min_element(dist.begin(), dist.end())-dist.begin();
      zyz(colon(), i) = zyz_frame(colon(), voronoi_id)+rand<double>(3, 1)*0.2;
      if(i < 100) {
          copy(fa->faces_[i].begin(), fa->faces_[i].end(), ostream_iterator<size_t>(cout, " "));
          cout << endl;
          cout << "node_idx: " << trans(node_idx) << trans(pos) << trans(dist) << voronoi_id << endl;
        }
    }

  ofstream ofs(argv[2], ofstream::binary);
  jtf::mesh::write_matrix(ofs, zyz);

  return 0;
}
