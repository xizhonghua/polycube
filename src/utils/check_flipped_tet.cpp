#include <fstream>
#include <jtflib/mesh/io.h>
#include <jtflib/mesh/util.h>

#include "../tetmesh/tetmesh.h"
#include "../tetmesh/hex_io.h"

#include "../common/vtk.h"

using namespace std;
using namespace zjucad::matrix;

int check_flipped_tet(int argc, char *argv[])
{
  if(argc != 2){
    cerr << "# [info] check_flipped_tet tet" << endl;
    return __LINE__;
  }

  jtf::mesh::meshes tm;
  if(jtf::mesh::tet_mesh_read_from_zjumat(argv[1], &tm.node_, &tm.mesh_))
    return __LINE__;
  vector<size_t> flipped_tet_idx;
  for(size_t ti = 0; ti < tm.mesh_.size(2); ++ti){
    const double v = jtf::mesh::cal_tet_vol(tm.node_(colon(), tm.mesh_(colon(),ti)));
    if(v < 0)
      flipped_tet_idx.push_back(ti);
//      cerr << "# [info] tet " << tm.mesh_(0,ti) << " "
//           << tm.mesh_(1,ti) << " " << tm.mesh_(2,ti) << " "  << tm.mesh_(3,ti) << " "
//           << " is flipped." << endl;
  }

  {
    cerr << "# [info] flipped tet number: " << flipped_tet_idx.size() << endl;
    ofstream ofs("flipped_tet.vtk");
    matrixst selected_tet(4, flipped_tet_idx.size());
    for(size_t ti = 0; ti <flipped_tet_idx.size(); ++ti){
      selected_tet(colon(), ti) = tm.mesh_(colon(), flipped_tet_idx[ti]);
    }
    tet2vtk(ofs, &tm.node_[0], tm.node_.size(2), &selected_tet[0], selected_tet.size(2));
  }

  return 0;
}
