#include <fstream>
#include <iostream>

#include <jtflib/mesh/io.h>
#include <jtflib/mesh/mesh.h>

using namespace std;
using namespace jtf::mesh;

int tet2ascii(int argc, char * argv[])
{
  if(argc != 3){
      cerr << "# [error] tet2ascii tet_file output_ascii." << endl;
      return __LINE__;
    }

  jtf::mesh::meshes tm;

  if(tet_mesh_read_from_zjumat(argv[1], &tm.node_, &tm.mesh_))
    return __LINE__;

  ofstream ofs(argv[2]);
  ofs << tm.node_.size(1) << " " << tm.node_.size(2) << endl;
  for(size_t ni = 0; ni < tm.node_.size(2); ++ni)
    ofs << tm.node_(0,ni) << " " << tm.node_(1,ni) << " " << tm.node_(2,ni) << endl;
  ofs << tm.mesh_.size(1) << " " << tm.mesh_.size(2) << endl;
  for(size_t ni = 0; ni < tm.mesh_.size(2); ++ni)
    ofs << tm.mesh_(0,ni) << " " << tm.mesh_(1,ni) << " " << tm.mesh_(2,ni) << " " << tm.mesh_(3,ni) << endl;

  return 0;
}
