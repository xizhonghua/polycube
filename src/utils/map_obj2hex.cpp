#include <iostream>
#include <jtflib/mesh/io.h>
#include <jtflib/mesh/mesh.h>
#include "../hexmesh/hexmesh.h"

using namespace std;
using namespace zjucad::matrix;

int map_obj2hex(int argc, char *argv[])
{
  if(argc != 4 && argc != 5){
      cerr << "# [usage] map_obj2hex obj hex new_hex [s2v]" << endl;
      return __LINE__;
    }

  jtf::mesh::meshes obj;
  if(jtf::mesh::load_obj(argv[1], obj.mesh_, obj.node_))
    return __LINE__;

  jtf::hex_mesh hm(argv[2]);
  if(argc == 4)
    hm.hexmesh_.node_ = obj.node_;
  else if(argc == 5){
      matrix<int> s2v;
      jtf::mesh::read_matrix(argv[4], s2v);
      for(size_t i = 0; i < s2v.size(); ++i){
          hm.hexmesh_.node_(colon(), s2v[i]) = obj.node_(colon(),i);
        }
    }

  jtf::hexmesh::hex_mesh_write_to_wyz(argv[3], hm.hexmesh_.mesh_, hm.hexmesh_.node_);

  return 0;
}
