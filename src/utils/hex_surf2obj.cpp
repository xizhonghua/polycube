#include <fstream>
#include <iostream>
#include <vector>
#include <algorithm>
using namespace std;

#include <jtflib/mesh/mesh.h>
#include "../common/vtk.h"
#include "../common/util.h"
#include "../common/IO.h"
#include <jtflib/mesh/mesh.h>
#include "../hexmesh/io.h"

using namespace zjucad::matrix;
//using namespace jtf::hexmesh;

int hex_surf2obj(int argc, char *argv[])
{
  if(argc != 3 && argc != 4) {
      cerr << "# [usage] hex_surf2obj hex obj [s2v]" << endl;
      return __LINE__;
    }

  jtf::mesh::meshes hm;
  matrixst quad;

  if(jtf::hexmesh::hex_mesh_read_from_wyz(argv[1], hm.mesh_, hm.node_,1))
    return __LINE__;

  cerr << "# hex number " << hm.mesh_.size(2) << endl;
  cerr << "# hex point " << hm.node_.size(2) << endl;

  unique_ptr<jtf::mesh::face2hex_adjacent> fa(jtf::mesh::face2hex_adjacent::create(hm.mesh_));
  if(!fa.get()){
      cerr << "# [error] can not build face2hex_adjacent." << endl;
      return __LINE__;
    }

  jtf::mesh::get_outside_face(*fa, quad);

  if(argc == 3)
    remove_extra_node(quad, hm.node_);
  if(argc == 4){
      matrix<size_t> orig2new_mapping;
      remove_extra_node(quad, hm.node_, &orig2new_mapping);
      matrix<int> s2v(hm.node_.size(2),1);
      for(size_t i = 0; i < orig2new_mapping.size(); ++i){
          if(orig2new_mapping[i] != -1)
            s2v[orig2new_mapping[i]] = i;
        }
      jtf::mesh::write_matrix(argv[3], s2v);
    }

  jtf::mesh::save_obj(argv[2], quad, hm.node_);

  return 0;
}
