#include <iostream>
#include <jtflib/mesh/io.h>

#include "../tetmesh/util.h"

using namespace std;
using namespace jtf::tetmesh;

int improve_tet(int argc, char *argv[])
{
  if(argc != 3){
    cerr << "# [usage] improve_tet input_tet output_tet." << endl;
    return __LINE__;
  }

  jtf::mesh::meshes tm;
  if(jtf::mesh::tet_mesh_read_from_zjumat(argv[1], &tm.node_, &tm.mesh_))
    return __LINE__;

  tetmesh_quality_improver tqi(tm.mesh_,tm.node_);
  tqi.improve("fix_surface");

  if(jtf::mesh::tet_mesh_write_to_zjumat(argv[2], &tm.node_, &tm.mesh_))
    return __LINE__;

  return 0;
}
