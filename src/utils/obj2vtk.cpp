#include <zjucad/matrix/matrix.h>
#include <zjucad/matrix/itr_matrix.h>

#include <jtflib/mesh/mesh.h>
#include <jtflib/mesh/util.h>
#include <jtflib/mesh/io.h>
#include "../common/vtk.h"

#include <fstream>
#include <iostream>

using namespace std;
using namespace zjucad::matrix;

int obj2vtk(int argc, char * argv[])
{
  if(argc != 2){
      cerr << "# [usage] obj2vtk obj." << endl;
      return __LINE__;
    }

  jtf::mesh::meshes trm;
  if(jtf::mesh::load_obj(argv[1], trm.mesh_, trm.node_))
    return __LINE__;

  string vtk_name = argv[1];
  vtk_name += ".vtk";
  ofstream ofs(vtk_name.c_str());
  tri2vtk(ofs, &trm.node_[0], trm.node_.size(2), &trm.mesh_[0], trm.mesh_.size(2));
  return 0;
}

int obj2obj(int argc, char * argv[])
{
  if(argc != 2){
      cerr << "# [usage] obj2obj obj." << endl;
      return __LINE__;
    }

  jtf::mesh::meshes trm;
  if(jtf::mesh::load_obj(argv[1], trm.mesh_, trm.node_))
    return __LINE__;
  if(jtf::mesh::save_obj("test.obj", trm.mesh_, trm.node_))
    return __LINE__;
  return 0;
}
