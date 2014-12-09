#include <iostream>
#include "../hexmesh/util.h"
#include "../hexmesh/io.h"

using namespace std;
using namespace jtf::hexmesh;

int improve_hex(int argc, char *argv[])
{
  if(argc != 3 && argc != 4 ){
    cerr << "# [usage] improve_hex input_hex output_hex strategy[shape/lap/size]." << endl;
    return __LINE__;
  }

  jtf::mesh::meshes hm;
  if(hex_mesh_read_from_wyz(argv[1], hm.mesh_, hm.node_, 1))
    return __LINE__;

  hexmesh_quality_improver hqi(hm.mesh_,hm.node_);

  string strategy;
  if(argc == 3) strategy = "shape";
  else strategy = argv[3];
  hqi.improve(strategy);

  if(hex_mesh_write_to_wyz(argv[2], hm.mesh_, hm.node_))
    return __LINE__;

  return 0;
}
