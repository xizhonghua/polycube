#include <fstream>
#include "io.h"

#include "../common/vtk.h"
#include "../common/IO.h"
#include <jtflib/mesh/mesh.h>
#include <jtflib/mesh/io.h>

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

int dump_surface_restricted_type_to_vtk_static(
    const char * filename,
    const std::string & type_name,
    const matrixd & node,
    const jtf::mesh::face2tet_adjacent & fa,
    const boost::unordered_map<size_t,size_t> & surface_type)
{
  ofstream ofs(filename);
  if(ofs.fail()){
    cerr << "# [error] can not open file." << endl;
    return __LINE__;
  }

  vector<size_t> face_vec;
  face_vec.reserve(3 * surface_type.size());
  vector<size_t> face_type;
  face_type.reserve(surface_type.size());
  for(boost::unordered_map<size_t,size_t>::const_iterator cit = surface_type.begin();
      cit != surface_type.end(); ++cit){
    const vector<size_t> & one_face = fa.faces_.at(cit->first);
    face_vec.insert(face_vec.end(), one_face.begin(), one_face.end());
    face_type.push_back(cit->second);
  }

  tri2vtk(ofs, &node[0], node.size(2), &face_vec[0], face_vec.size()/3);
  cell_data(ofs, &face_type[0], face_type.size(), const_cast<char*>(type_name.c_str()));
  return 0;
}

int load_fnode_group_static(
    const char * filename,
    const matrixst & cut_tet,
    std::vector<boost::unordered_set<size_t> > & fnode_group)
{
  ifstream ifs(filename);
  if(ifs.fail()){
    cerr << "# [error] can not open fnode_group file." << endl;
    return __LINE__;
  }

  fnode_group.clear();
  boost::unordered_set<size_t> one_group;
  size_t group_idx, group_size, type_trash;
  string trash;
  size_t f_node = -1;
  while(!ifs.eof()){
    ifs >> trash >> group_idx >> group_size >> type_trash;
    if(trash.empty()) break;
    one_group.clear();
    for(size_t i = 0; i < group_size; ++i){
      ifs >> f_node;
      one_group.insert(f_node);
    }
    fnode_group.push_back(one_group);
    trash.clear();
  }

  return 0;
}



int load_point_diffused_value(
    const char * filename,
    zjucad::matrix::matrix<double> & point_diffused)
{
  ifstream ifs(filename);
  if(ifs.fail()){
    cerr << "# [error] can not open point_diffused_value file." << endl;
    return __LINE__;
  }
  
  return  jtf::mesh::read_matrix(ifs, point_diffused);
}

int load_surface_patch(
    const char * surface_patch_file,
    const char * obj_file,
    const char * s2v_file,
    std::vector<matrix<size_t> > &patch_faces)
{
  ifstream ifs_patch(surface_patch_file);
  if(ifs_patch.fail()){
    cerr << "# [error] can not open patch surface file." << endl;
    return __LINE__;
  }

  size_t patch_num = 0;
  ifs_patch >> patch_num;

  cerr << "# [info] patch number " << patch_num << endl;

  if(patch_num == 0) return 0;

  jtf::mesh::meshes trm;
  if(jtf::mesh::load_obj(obj_file, trm.mesh_, trm.node_))
    return __LINE__;

  matrix<int32_t> surf_node2vol_nod;
  if(jtf::mesh::read_matrix(s2v_file, surf_node2vol_nod))
    return __LINE__;

  trm.mesh_(colon()) = surf_node2vol_nod(trm.mesh_);

  patch_faces.clear();

  patch_faces.resize(patch_num);
  size_t face_num = 0, face_idx, trash_1 ,trash_2;

  for(size_t pi = 0; pi < patch_num; ++pi){
    ifs_patch >> face_num >> trash_1;
    matrix<size_t> one_patch(3, face_num);
    for(size_t fi = 0; fi < face_num; ++fi){
      ifs_patch >> face_idx;
      one_patch(colon(),fi) = trm.mesh_(colon(), face_idx);
    }
    patch_faces[pi] = one_patch;
    for(size_t li = 0; li < trash_1; ++li)
      ifs_patch >> trash_2;
  }
  return 0;
}

int save_to_uv(const char * filename,
               const vector<size_t> &uv_basis_point,
               const vector<size_t> &boundary_nodes,
               const matrix<double> &polycube_node)
{
  if(uv_basis_point.size() != 3){
    cerr << "# [error] wrong uv_basis_point." << endl;
    return __LINE__;
  }

  ofstream ofs(filename);
  if(ofs.fail()){
    cerr << "# [error] can not open uv file." << endl;
    return __LINE__;
  }

  for(size_t pi = 0; pi < uv_basis_point.size(); ++pi){
    ofs << "# " << uv_basis_point[pi] << " ";
    for(size_t di = 0; di < polycube_node.size(1); ++di)
      ofs << polycube_node(di, uv_basis_point[pi]) << " ";
    ofs << endl;
  }

  matrix<double> e1 = zeros<double>(3,1),e2 = zeros<double>(3,1);
  e1 = polycube_node(colon(), uv_basis_point[1]) - polycube_node(colon(), uv_basis_point[0]);
  e2 = polycube_node(colon(), uv_basis_point[2]) - polycube_node(colon(), uv_basis_point[0]);

  const double len1 = norm(e1), len2 = norm(e2);
  if(len1 < 1e-6 || len2 < 1e-6){
    cerr << "# [error] degenerated basis." << endl;
    return __LINE__;
  }

  e1 /= len1;
  e2 /= len2;

  ofs << "# boundary" << endl;
  for(size_t pi = 0; pi < boundary_nodes.size(); ++pi){
    ofs << boundary_nodes[pi] + 1 << " "; // index in obj
    ofs << dot(e1, polycube_node(colon(), boundary_nodes[pi])
               - polycube_node(colon(), uv_basis_point[0])) << " ";
    ofs << dot(e2, polycube_node(colon(), boundary_nodes[pi])
               - polycube_node(colon(), uv_basis_point[0])) << endl;
  }

  ofs << "# init" << endl;
  for(size_t pi = 0; pi < polycube_node.size(2); ++pi){
    ofs << pi + 1 << " "; // index in obj
    ofs << dot(e1, polycube_node(colon(), pi)
               - polycube_node(colon(), uv_basis_point[0])) << " ";
    ofs << dot(e2, polycube_node(colon(), pi)
               - polycube_node(colon(), uv_basis_point[0])) << endl;
  }
  return 0;
}


int load_from_uv(const char * filename,
                 const jtf::mesh::meshes & trm,
                 zjucad::matrix::matrix<double> & uv,
                 zjucad::matrix::matrix<double> * orig_node_ptr,
                 zjucad::matrix::matrix<double> * uv_basis_ptr)
{
  ifstream ifs(filename);
  if(ifs.fail()){
    cerr << "# [error] can not open uv file." << endl;
    return __LINE__;
  }

  matrix<double> uv_basis(3,2);
  string string; size_t idx;
  matrix<double> point(3,3);
  for(size_t p = 0; p < 3; ++p){
    ifs >> string >> idx >> point(0,p) >> point(1,p) >> point(2,p);
    if(string != "#") {
      cerr << "# [error] wrong uv file." << endl;
      return __LINE__;
    }
  }

  uv_basis(colon(),0) = point(colon(),1) - point(colon(),0);
  uv_basis(colon(),1) = point(colon(),2) - point(colon(),0);

  double len_0 = norm(uv_basis(colon(),0));
  if(len_0 < 1e-6) len_0 = 1.0;
  uv_basis(colon(),0) /= len_0;
  double len_1 = norm(uv_basis(colon(),1));
  if(len_1 < 1e-6) len_1 = 1.0;
  uv_basis(colon(),1) /= len_1;

  uv = zeros<double>(2, trm.node_.size(2));
  for(size_t p = 0; p < trm.node_.size(2); ++p){
    uv(0,p) = dot(trm.node_(colon(),p) - point(colon(),0),
                  uv_basis(colon(),0));
    uv(1,p) = dot(trm.node_(colon(),p) - point(colon(),0),
                  uv_basis(colon(),1));
  }

  if(orig_node_ptr){
    *orig_node_ptr = point(colon(),0);
  }

  if(uv_basis_ptr){
    *uv_basis_ptr = uv_basis;
  }
  return 0;
}
