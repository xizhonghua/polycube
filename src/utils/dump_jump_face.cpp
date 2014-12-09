#include <map>
#include <set>
#include <fstream>

#include <jtflib/mesh/io.h>

#include "../tetmesh/tetmesh.h"
#include "../hex_param/hex_param.h"
#include "../hex_param/io.h"
#include "../common/transition_type.h"
#include "../common/vtk.h"

using namespace std;

int dump_jump_face_to_vtk(int argc, char * argv[])
{
  if(argc != 3){
      cerr << "# [error] wrong num arguments: dump_jump_face_to_vtk tet inner_face_jump_type." << endl;
      return __LINE__;
    }

  jtf::mesh::meshes tm;
  if(jtf::mesh::tet_mesh_read_from_zjumat(argv[1],&tm.node_,&tm.mesh_)){
      cerr << "# [error] read tet mesh error" << endl;
      return __LINE__;
    }

  boost::unordered_map<pair<size_t,size_t>,size_t> inner_face_jump_type;
  if(load_inner_face_jump_type(argv[2], inner_face_jump_type))
    return __LINE__;

  unique_ptr<jtf::mesh::face2tet_adjacent> fa(jtf::mesh::face2tet_adjacent::create(tm.mesh_));
  if(!fa.get()){
      cerr << "# [error] can not buildjtf::mesh::face2tet_adjacent." << endl;
      return __LINE__;
    }

  vector<size_t> jump_face;
  vector<size_t> jump_type;

  vector<size_t> one_face_vec(3);
  for(const auto & one_face : inner_face_jump_type){
      if(one_face.second == TRIVIAL_TYPE) continue;
      const pair<size_t,size_t> & two_tets = one_face.first;
      jtf::mesh::find_common_face(tm.mesh_(colon(),two_tets.first), tm.mesh_(colon(),two_tets.second),
                                  &one_face_vec[0]);
      jump_face.insert(jump_face.end(), one_face_vec.begin(), one_face_vec.end());
      jump_type.push_back(one_face.second);
    }
  {// visualization
    tri2vtk(cout,&tm.node_[0],tm.node_.size(2),&jump_face[0],jump_face.size()/3);
    cell_data(cout,&jump_type[0],jump_type.size(),"uvw");
  }
  return 0;
}
