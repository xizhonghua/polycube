#include <iostream>
#include <fstream>
#include <jtflib/mesh/io.h>
#include <jtflib/mesh/mesh.h>

#include "../tetmesh/tetmesh.h"
using namespace std;
using namespace jtf::mesh;
using namespace zjucad::matrix;

int load_off(const char * file,
             matrix<size_t> & mesh,
             matrix<double> & node)
{
  ifstream ifs(file);
  if(ifs.fail()){
      cerr << "# [error] can not open off file." << endl;
      return __LINE__;
    }

  string trash;
  ifs >> trash; // "OFF"
  size_t point_num, cell_num;
  ifs >> point_num >> cell_num >> trash;
  node.resize(3, point_num);
  for(size_t i = 0; i < point_num; ++i)
    ifs >> node(0,i) >> node(1,i) >> node(2,i);
  size_t cell_point_num = 0;
  ifs >> cell_point_num; // to see how many points a cell should have

  mesh.resize(cell_point_num, cell_num);
  for(size_t i = 0; i < cell_num; ++i){
      for(size_t j = 0; j < cell_point_num; ++j)
        ifs >> mesh(j,i);
      ifs >> trash;
    }

  return 0;
}

int off2cell(int argc, char * argv[])
{
  if(argc != 4){
      cerr << "# [usage] off2cell tet/obj input_cell output_cell" << endl;
      return __LINE__;
    }

  jtf::mesh::meshes mesh;
  const string type = argv[1];
  const string input_cell_name = argv[2];
  const string output_cell_name = argv[3];

  if(load_off(input_cell_name.c_str(), mesh.mesh_, mesh.node_))
    return __LINE__;

  if(type == "obj"){
      if(jtf::mesh::save_obj(output_cell_name.c_str(), mesh.mesh_, mesh.node_))
        return __LINE__;
    }else if(type == "tet") {
      if(jtf::mesh::tet_mesh_write_to_zjumat(output_cell_name.c_str(), &mesh.node_, &mesh.mesh_))
        return __LINE__;
    }else {
      cerr << "# [error] do not support this cell type " << type << endl;
      return __LINE__;
    }

  return 0;
}
