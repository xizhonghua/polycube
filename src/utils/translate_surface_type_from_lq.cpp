#include <fstream>
#include <jtflib/mesh/io.h>

#include "../tetmesh/hex_io.h"

using namespace std;

int translate_surface_type_from_lq(int argc, char * argv[])
{
  if(argc != 4){
    cerr << "# [usage] translate_surface_type_from_lq tet surface_type_lq surface_type." << endl;
    return __LINE__;
  }

  jtf::mesh::meshes tm;
  if(jtf::mesh::tet_mesh_read_from_zjumat(argv[1], &tm.node_, &tm.mesh_))
    return __LINE__;

  unique_ptr<jtf::mesh::face2tet_adjacent> fa(jtf::mesh::face2tet_adjacent::create(tm.mesh_));

  if(!fa.get()){
    cerr << "# [error] can  not buildjtf::mesh::face2tet_adjacent." << endl;
    return __LINE__;
  }
\
  ifstream ifs(argv[2]);
  if(ifs.fail()){
    cerr << "# [error] can not open surface_type_lq file." << endl;
    return __LINE__;
  }

  ofstream ofs(argv[3]);
  if(ofs.fail()){
    cerr << "# [error] can not open output surface type." << endl;
    return __LINE__;
  }

  size_t face[3], type;
  while(!ifs.eof()){
    ifs >> face[0] >> face[1] >> face[2] >> type;
    if(ifs.eof())
      break;
    const size_t face_idx = fa->get_face_idx(&face[0]);
    if(face_idx == -1){
      cerr << "# [error] can not find face in tet." << endl;
      return __LINE__;
    }
    ofs << face_idx << " " << type << endl;
  }

  return 0;
}
