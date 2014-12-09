#include "../common/IO.h"
#include "../tetmesh/hex_io.h"
#include "../common/vtk.h"
#include <iostream>

#include <zjucad/matrix/matrix.h>

using namespace std;
using namespace zjucad::matrix;

int zyz_mapping(int argc, char *argv[])
{
  if(argc != 4){
      cerr << "# [usage] zyz_mapping src_zyz mapping output_zyz" << endl;
      return __LINE__;
    }

  matrix<uint32_t> mapping;
  if(jtf::mesh::read_matrix(argv[2], mapping)){
      cerr << "# [error] can not open mapping matrix." << endl;
      return __LINE__;
    }

  matrixd src_zyz;

  if(jtf::mesh::read_matrix(argv[1], src_zyz)){
      cerr << "# [error] can not open src zyz matrix. " << endl;
      return __LINE__;
    }

  if(*max_element(mapping.begin(), mapping.end()) >= src_zyz.size(2)){
      cerr << "# [error] error mapping:" << endl;
      cerr << "# src_zyz size: " << src_zyz.size(1) << " " << src_zyz.size(2) << endl;
      cerr << "# mapping size " << mapping.size(1) << " " << mapping.size(2)
           << " max " << *max_element(mapping.begin(), mapping.end())
           << " min " << *min_element(mapping.begin(), mapping.end()) << endl;
      return __LINE__;
    }

  matrixd des_zyz = zeros<double>(3, mapping.size());
  for(size_t ti = 0; ti < des_zyz.size(2); ++ti){
      des_zyz(colon(), ti) = src_zyz(colon(),mapping[ti]);
    }

  if(jtf::mesh::write_matrix(argv[3], des_zyz)){
      cerr << "# [error] can not open output_zyz file." << endl;
      return __LINE__;
    }

  return 0;
}


int fix_type(int argc, char * argv[])
{
  if(argc != 7){
      cerr << "# [error] fix_type tet 0 1 orig face 2" << endl;
      return __LINE__;
    }

  jtf::mesh::meshes tm;
  if(jtf::mesh::tet_mesh_read_from_zjumat(argv[1], &tm.node_, &tm.mesh_))
    return __LINE__;

  unique_ptr<jtf::mesh::face2tet_adjacent> fa(jtf::mesh::face2tet_adjacent::create(tm.mesh_));
  size_t face_num = 0;
  {
    matrix<size_t> outside_face;
    get_outside_face(*fa, outside_face);
    face_num = outside_face.size(2);
  }
    matrix<size_t> outside_face_orig(3, face_num);
  matrix<size_t> outside_face_type(face_num,1);
  ifstream ifs_face(argv[5]);
  ifs_face >> face_num;
  for(size_t fi = 0; fi < outside_face_orig.size(2); ++fi){
      ifs_face >> outside_face_orig(0,fi)
               >> outside_face_orig(1,fi)
               >> outside_face_orig(2,fi);
    }

  vector<size_t> idx_map(outside_face_orig.size(2));
  for(size_t fi = 0; fi < outside_face_orig.size(2); ++fi){
      const size_t face_idx = fa->get_face_idx(&outside_face_orig(0,fi));
      idx_map[fi] = face_idx;
    }


  map<size_t,size_t> face_idx_type;

  ifstream ifs(argv[4]);
  size_t type;
  for(size_t fi = 0; fi < outside_face_type.size(); ++fi){
      ifs >> type;
      face_idx_type[idx_map[fi]] = type;
    }


  size_t num = 0;
  size_t idx = 0;

  ifstream ifs_0(argv[2]);
  ifs_0 >> num;

  for(size_t i = 0; i < num; ++i){
      ifs_0 >> idx;
      face_idx_type[idx_map[idx]] = 0;
    }

//  ifstream ifs_1(argv[3]);
//  ifs_1 >> num;

//    for(size_t i = 0; i < num; ++i){
//        ifs_1 >> idx;
//        face_idx_type[idx_map[idx]] = 1;
//      }

    ifstream ifs_2(argv[6]);
    ifs_2 >> num;

      for(size_t i = 0; i < num; ++i){
          ifs_2 >> idx;
          face_idx_type[idx_map[idx]] = 2;
        }

  matrix<size_t> face_type(outside_face_orig.size(2),1);
  for(size_t fi = 0; fi < outside_face_orig.size(2); ++fi){
      face_type[fi] = face_idx_type[idx_map[fi]];
    }

  ofstream ofs("fixed_surface_patch.vtk");
  tri2vtk(ofs, &tm.node_[0], tm.node_.size(2), &outside_face_orig[0], outside_face_orig.size(2));
  cell_data(ofs, &face_type[0], face_type.size(), "type");
  return 0;
}
