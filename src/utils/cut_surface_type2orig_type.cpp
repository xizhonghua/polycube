#include <jtflib/mesh/io.h>

#include "../tetmesh/tetmesh.h"
#include "../tetmesh/hex_io.h"
#include "../hex_param/io.h"

using namespace std;
using namespace zjucad::matrix;

int cut_surface_type2orig_type(int argc, char * argv[])
{
  if(argc != 4){
    cerr << "# [usage] cut_surface_type2orig_type tet cut_tet cut_surface_type." << endl;
    return __LINE__;
  }

  jtf::mesh::meshes tm, cut_tm;
  if(jtf::mesh::tet_mesh_read_from_zjumat(argv[1], &tm.node_, &tm.mesh_))
    return __LINE__;
  if(jtf::mesh::tet_mesh_read_from_zjumat(argv[2], &cut_tm.node_, &cut_tm.mesh_))
    return __LINE__;

  if(tm.mesh_.size(2) != cut_tm.mesh_.size(2)){
    cerr << "# [error] orig tet is not compatiable with cut_tet." << endl;
    return __LINE__;
  }

  boost::unordered_map<size_t,size_t> surface_type;
  if(load_surface_type(argv[3], surface_type))
    return __LINE__;
  unique_ptr<jtf::mesh::face2tet_adjacent> fa_cut(jtf::mesh::face2tet_adjacent::create(cut_tm.mesh_));
  unique_ptr<jtf::mesh::face2tet_adjacent> fa(jtf::mesh::face2tet_adjacent::create(tm.mesh_));
  if(!fa_cut.get() || !fa.get()){
    cerr << "# [error] can not buildjtf::mesh::face2tet_adjacent." << endl;
    return __LINE__;
  }

  matrixst cut_tet2tet(max(cut_tm.mesh_)+1);
  cut_tet2tet(cut_tm.mesh_) = tm.mesh_(colon());

  for(boost::unordered_map<size_t,size_t>::const_iterator cit
      = surface_type.begin(); cit != surface_type.end(); ++cit){
    const size_t & face_idx = cit->first;
    const vector<size_t> & face_vec = fa_cut->faces_[face_idx];
    const size_t & face_idx_orig =
        fa->get_face_idx(cut_tet2tet[face_vec[0]], cut_tet2tet[face_vec[1]],
                         cut_tet2tet[face_vec[2]]);

    if(face_idx_orig == -1){
      cerr << "# [error] strange face " << cut_tet2tet[face_vec[0]] << " "
           << cut_tet2tet[face_vec[1]]  << " " << cut_tet2tet[face_vec[2]]
           << " is not inside the orig tetmesh" << endl;
      return __LINE__;
    }
    cout << face_idx_orig << " " << cit->second << endl;
  }

  return 0;
}
