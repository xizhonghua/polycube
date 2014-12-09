#include <iostream>
#include <jtflib/mesh/io.h>
#include <fstream>
#include <sstream>

#include "../common/util.h"
#include "../tetmesh/tetmesh.h"
#include "../tetmesh/hex_io.h"
#include "../tetmesh/util.h"

using namespace std;
using namespace zjucad::matrix;

int load_surface_restricted_type_static(
    const char * filename,
    boost::unordered_map<size_t,size_t> & result_surface_type)
{
  ifstream ifs(filename);
  if(ifs.fail()){
    cerr << "# [error] can not open surface restricted type." << endl;
    return __LINE__;
  }
  result_surface_type.clear();
  size_t face_idx, face_type;
  while(!ifs.eof()){
    ifs >> face_idx >> face_type;
    if(face_type > 2) {
      cerr << "# [error] restrcited type should only be 0/1/2." << endl;
      return __LINE__;
    }
    result_surface_type[face_idx] = face_type;
  }
  return 0;
}

int sft2obj(int argc, char * argv[])
{
  if(argc != 4 ){
    cerr << "# [usage] sft2obj tet surface_type output_obj." << endl;
    return __LINE__;
  }

  jtf::mesh::meshes tm;
  if(jtf::mesh::tet_mesh_read_from_zjumat(argv[1], &tm.node_, &tm.mesh_))
    return __LINE__;

  boost::unordered_map<size_t,size_t> surface_type;
  if(load_surface_restricted_type_static(argv[2], surface_type))
    return __LINE__;

  unique_ptr<jtf::mesh::face2tet_adjacent> fa(jtf::mesh::face2tet_adjacent::create(tm.mesh_));
  if(!fa.get()){
    cerr << "# [error] can not buildjtf::mesh::face2tet_adjacent." << endl;
    return __LINE__;
  }

  matrix<size_t> outside_face, outside_face_idx;
  get_outside_face(*fa, outside_face);
  get_outside_face_idx(*fa, outside_face_idx);

  vector<vector<size_t> > face_patch(3); // u,v,w
  for(size_t fi = 0; fi < outside_face_idx.size(); ++fi){
    boost::unordered_map<size_t,size_t>::const_iterator cit =
        surface_type.find(outside_face_idx[fi]);
    if(cit == surface_type.end() || cit->second > 2){
      cerr << "# [error] invalid surface type." << endl;
      return __LINE__;
    }
    face_patch[cit->second].insert(
          face_patch[cit->second].end(), outside_face(colon(),fi).begin(),
          outside_face(colon(),fi).end());
  }

  matrix<double> node;
  matrix<size_t> face;
  string output_name = argv[3];
  string obj_name;
  for(size_t di = 0; di < face_patch.size(); ++di){
    const vector<size_t> & one_group = face_patch[di];
    face.resize(3, one_group.size()/3);
    copy(one_group.begin(), one_group.end(), face.begin());
    node = tm.node_;
    remove_extra_node(face, node);
    ostringstream os;
    os << "_" << di << ".obj";
    obj_name = output_name + os.str();
    jtf::mesh::save_obj(obj_name.c_str(), face, node);
  }
  return 0;
}
