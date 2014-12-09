#include "../hexmesh/io.h"
#include <jtflib/mesh/mesh.h>
#include "../hexmesh/util.h"
#include <iostream>

using namespace std;
using namespace jtf::hexmesh;

int vtk2hex(int argc, char * argv[])
{
  if(argc != 3) {
    cerr << "# [usage] vtk2hex vtk hex."<< endl;
    return __LINE__;
  }

  jtf::mesh::meshes hm;
  if(hex_mesh_read_from_vtk(argv[1], &hm.node_, &hm.mesh_))
    return __LINE__;
  cerr << "# [info] node " << hm.node_.size(2) << " hex " << hm.mesh_.size(2) << endl;
  {
//    orient_hex(hm.mesh_, hm.node_);
//    if(jtf::hexmesh::hex_mesh_write_to_wyz("fixed.hex", hm.mesh_, hm.node_))
//      return __LINE__;
  }
  if(hex_mesh_write_to_wyz(argv[2], hm.mesh_, hm.node_))
    return __LINE__;
  return 0;
}
