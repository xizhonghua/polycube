#include "../tetmesh/tetmesh.h"
#include "../tetmesh/hex_io.h"
#include <iostream>

#include <jtflib/mesh/mesh.h>
#include <jtflib/mesh/io.h>


using namespace std;

int validate_tetmesh(int argc, char * argv[])
{
  if(argc != 2){
    cerr << "# [usage] validate_tetmesh input_tet." << endl;
    return __LINE__;
  }

  jtf::mesh::meshes tm;
  if(jtf::mesh::tet_mesh_read_from_zjumat(argv[1], &tm.node_, &tm.mesh_))
    return __LINE__;

  unique_ptr<jtf::mesh::face2tet_adjacent> fa(jtf::mesh::face2tet_adjacent::create(tm.mesh_));
  if(!fa.get()){
    cerr <<  "# [error] can not buildjtf::mesh::face2tet_adjacent." << endl;
    return __LINE__;
  }

  matrixst outside_face;
  get_outside_face(*fa, outside_face);

  unique_ptr<jtf::mesh::edge2cell_adjacent> ea(
        jtf::mesh::edge2cell_adjacent::create(outside_face));
  if(!ea.get()){
    cerr << "# [error] can not build edge2cell_adjacent." << endl;
    return __LINE__;
  }

  cerr << "# [info] valid." << endl;
  return 0;
}
