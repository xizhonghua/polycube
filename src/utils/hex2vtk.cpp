#include <fstream>
#include <iostream>
#include <vector>
#include <algorithm>

#include <zjucad/matrix/matrix.h>
#include "../common/vtk.h"

#include "../hexmesh/util.h"
#include <jtflib/mesh/mesh.h>
#include "../hexmesh/io.h"
using namespace std;
using namespace jtf::hexmesh;
using namespace zjucad::matrix;

int hex2vtk(int argc, char *argv[])
{
  if(argc < 3) {
      cerr << "useage: hex2vtk hex_format[1/2/3,defaul 1] hex" << endl;
      return __LINE__;
    }

  jtf::mesh::meshes hm;
  const size_t hex_format = atoi(argv[1]);

  if(hex_mesh_read_from_wyz(argv[2], hm.mesh_, hm.node_, hex_format))
    return __LINE__;

  jtf::hexmesh::orient_hex(hm.mesh_, hm.node_);

  hex2vtk(cout, &hm.node_[0], hm.node_.size(2), &hm.mesh_[0], hm.mesh_.size(2));

  return 0;
}

