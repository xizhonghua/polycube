#include <iostream>
#include <fstream>
#include <jtflib/mesh/mesh.h>
#include "../hexmesh/io.h"
#include "../hexmesh/util.h"
#include "../common/vtk.h"

using namespace std;
using namespace jtf::hexmesh;

int subdivide_hexmesh(int argc, char * argv[])
{
  if(argc != 4){
    cerr << "# [usage] subdivide_hexmesh input_hex hex_format sub_num" << endl;
    return __LINE__;
  }

  jtf::mesh::meshes hm;
  if(hex_mesh_read_from_wyz(argv[1], hm.mesh_, hm.node_, atoi(argv[2])))
    return __LINE__;

  const size_t iter = atoi(argv[3]);
  for(size_t i = 0; i < iter; ++i)
    jtf::hexmesh::subdivide_hexmesh(hm.mesh_, hm.node_);

  ofstream ofs("sub_hex.vtk");
  hex2vtk(ofs, &hm.node_[0], hm.node_.size(2), &hm.mesh_[0], hm.mesh_.size(2));

  if(hex_mesh_write_to_wyz("sub_hex.hex", hm.mesh_, hm.node_))
    return __LINE__;
  return 0;
}
