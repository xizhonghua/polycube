#include "../hexmesh/io.h"
#include "../hexmesh/util.h"
#include <jtflib/mesh/io.h>
#include <iostream>
#include "../common/def.h"

using namespace std;
using namespace jtf::hexmesh;

int load_quad(const char * file, matrixst & quad, matrixd & node)
{
  if(jtf::mesh::load_obj(file, quad, node)){
    matrixst hex;
    if(hex_mesh_read_from_wyz(file, hex, node, 1)){
      cerr << "# [error] can not load quad mesh." << endl;
      return __LINE__;
    }
    std::unique_ptr<jtf::mesh::face2hex_adjacent> fa(jtf::mesh::face2hex_adjacent::create(hex));
    if(!fa.get()){
      cerr << "# [error] invalid quad mesh." << endl;
      return __LINE__;
    }
    jtf::mesh::get_outside_face(*fa, quad);
  }
  return 0;
}


int extract_quad_polyedge(int argc, char *argv[])
{
  if(argc != 2){
    cerr << "# [usage] extract_quad_polyedge hex[quad]"  << endl;
    return __LINE__;
  }

  matrixst quad;
  matrixd node;

  if(load_quad(argv[1], quad, node))
    return __LINE__;

  cerr << "# [error] not finished. " << endl;

  return 0;
}
