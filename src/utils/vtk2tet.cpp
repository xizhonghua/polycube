#include <iostream>
#include <jtflib/mesh/io.h>
#include <jtflib/mesh/util.h>
#include "../tetmesh/hex_io.h"

using namespace std;
int vtk2tet(int argc, char * argv[])
{
  if(argc != 3 && argc != 4) {
      cerr << "# [usage] vtk2tet vtk tet [orient_tet?1:Yes,0:No]"<< endl;
      return __LINE__;
    }

  jtf::mesh::meshes tm;
  if(jtf::mesh::tet_mesh_read_from_vtk(argv[1], &tm.node_, &tm.mesh_))
    return __LINE__;
  cerr << "# [info] node " << tm.node_.size(2) << " tet " << tm.mesh_.size(2) << endl;
  if(argc == 4){
      int need_orient = atoi(argv[3]);
      if(need_orient)
        orient_tet(tm.node_, tm.mesh_);
    }
  if(jtf::mesh::tet_mesh_write_to_zjumat(argv[2], &tm.node_, &tm.mesh_))
    return __LINE__;
  return 0;
}
