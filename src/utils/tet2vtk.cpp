#include <fstream>
#include <iostream>
#include <vector>
#include <algorithm>

#include <jtflib/mesh/io.h>

#include "../common/vtk.h"
#include "../tetmesh/tetmesh.h"

using namespace std;



int tet2vtk(int argc, char *argv[])
{
  if(argc < 2) {
      cerr << "tet2vtk tet" << endl;
      return __LINE__;
    }

  jtf::mesh::meshes tm;
  if(jtf::mesh::tet_mesh_read_from_zjumat(
       argv[1],
       &tm.node_, &tm.mesh_, 0))
    return __LINE__;

  tet2vtk(cout,
          &tm.node_[0], tm.node_.size(2),
      &tm.mesh_[0], tm.mesh_.size(2));

  return 0;
}
