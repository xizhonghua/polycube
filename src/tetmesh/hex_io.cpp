#include "hex_io.h"

#include <iostream>
#include <fstream>
#include <memory>

#include <zjucad/matrix/io.h>
#include <jtflib/mesh/util.h>
#include "../tetmesh/tetmesh.h"
#include "../common/util.h"
#include "../common/IO.h"
#include "../common/zyz.h"

using namespace zjucad::matrix;
using namespace std;

int surface_uv2zyz_frame(const char *surf_node_uv_file, matrixd &zyz_frame)
{
  ifstream uv_ifs(surf_node_uv_file);
  if(uv_ifs.fail()) {
      cerr << "# open " << surf_node_uv_file << " fail." << endl;
      return __LINE__;
    }
  size_t face_num;
  uv_ifs >> face_num;
  cerr << "# get " << face_num << " in suv." << endl;
  if(zyz_frame.size(2) != face_num)
    return __LINE__;

  matrixd frame(3, 3);
  for(size_t i = 0; i < face_num; ++i) {
      for(int j = 0; j < 6; ++j) uv_ifs >> frame[3+j];
      frame(colon(), 0) = cross(frame(colon(), 1), frame(colon(), 2));
      if(norm(frame*trans(frame) - eye<double>(3)) > 1e-1)
        cerr << "# not a frame?" << frame << frame*trans(frame) << endl;
      rotation_matrix_2_zyz_angle(&frame[0], &zyz_frame(0, i), 0);
    }
  if(uv_ifs.fail()) {
      cerr << "# read " << surf_node_uv_file << " fail." << endl;
      return __LINE__;
    }
  return 0;
}

int surface_frame2volme_frame(
    const char *surf_node_uv_file, const char *surf_node2tet_node_file,
    const char *obj_file, const jtf::mesh::face2tet_adjacent &fa,
    matrixd &surface_zyz_frame, matrixst &idx)
{
  if(surface_uv2zyz_frame(surf_node_uv_file, surface_zyz_frame))
    return __LINE__;
  const size_t face_num = surface_zyz_frame.size(2);

  matrix<int32_t> surf_node2tet_node; {
    ifstream s2v_ifs(surf_node2tet_node_file, ifstream::binary);
    if(s2v_ifs.fail()) {
        cerr << "# open " << surf_node2tet_node_file << " fail." << endl;
        return __LINE__;
      }
    jtf::mesh::read_matrix(s2v_ifs, surf_node2tet_node);
  }
  matrixst outside_face;
  get_outside_face(fa, outside_face);
  if(face_num != outside_face.size(2)) {
      cerr << "# incompatible face uv file: " << face_num << " " << outside_face.size(2) << endl;
      return __LINE__;
    }

  matrixst surf_faces(3, face_num); {
    ifstream ifs(obj_file);
    if(ifs.fail()) {
        cerr << "# open the obj file fail." << endl;
        return __LINE__;
      }
    size_t count = 0;
    while(1) {
        string line;
        getline(ifs, line);
        if(ifs.fail())
          break;
        char c;
        istringstream iss(line.c_str());
        iss >> c;
        if(c == 'f') {
            iss >> surf_faces(0, count) >> surf_faces(1, count) >> surf_faces(2, count);
            ++count;
          }
        if(count == surf_faces.size(2))
          break;
      }
    if(count != surf_faces.size(2)) {
        cerr << "# incompatible surface obj." << endl;
      }
    surf_faces -= 1;
  }

  for(size_t i = 0; i < face_num; ++i) {
      matrixst tet_node_idx = surf_node2tet_node(surf_faces(colon(), i));
      size_t tet_f_idx = fa.get_face_idx(&tet_node_idx[0]);
      if(tet_f_idx >= fa.faces_.size()) {
          cerr << "# no such surface face: " << i
               << surf_faces(colon(), i) << endl;
          return __LINE__;
        }
      if(!fa.is_outside_face(fa.face2tet_[tet_f_idx])) {
          cerr << "# is not a surface triangle for tri: " << i
               << surf_faces(colon(), i) << tet_node_idx << endl;
          return __LINE__;
        }
      idx[i] = tet_f_idx;
    }
  return 0;
}

int surface_frame2volme_frame_inner(
    const char *surf_node_uv_file, const char *surf_node2tet_node_file,
    const char *obj_file, const jtf::mesh::face2tet_adjacent &fa,
    matrixd &surface_zyz_frame,
    matrixst &idx)
{
  if(surface_uv2zyz_frame(surf_node_uv_file, surface_zyz_frame))
    return __LINE__;
  const size_t face_num = surface_zyz_frame.size(2);

  matrix<int32_t> surf_node2tet_node; {
    ifstream s2v_ifs(surf_node2tet_node_file, ifstream::binary);
    if(s2v_ifs.fail()) {
        cerr << "# open " << surf_node2tet_node_file << " fail." << endl;
        return __LINE__;
      }
    jtf::mesh::read_matrix(s2v_ifs, surf_node2tet_node);
  }
  matrixst outside_face;
  get_outside_face(fa, outside_face);
  if(face_num != outside_face.size(2)) {
      cerr << "# incompatible face uv file: " << face_num << " " << outside_face.size(2) << endl;
      return __LINE__;
    }

  matrixst surf_faces(3, face_num); {
    matrix<double> surf_node;
    jtf::mesh::load_obj(obj_file, surf_faces, surf_node);
  }

  for(size_t i = 0; i < face_num; ++i) {
      matrixst tet_node_idx = surf_node2tet_node(surf_faces(colon(), i));
      size_t tet_f_idx = fa.get_face_idx(&tet_node_idx[0]);
      if(tet_f_idx >= fa.faces_.size()) {
          cerr << "# no such surface face: " << i
               << surf_faces(colon(), i) << endl;
          return __LINE__;
        }
      if(!fa.is_outside_face(fa.face2tet_[tet_f_idx])) {
          cerr << "# is not a surface triangle for tri: " << i
               << surf_faces(colon(), i) << tet_node_idx << endl;
          return __LINE__;
        }
      const pair<size_t,size_t> & tet_pair= fa.face2tet_[tet_f_idx];
      idx[i] = (tet_pair.first == -1?tet_pair.second:tet_pair.first);
    }
  return 0;
}

int load_from_tet_inner(
    const matrixd &node, const matrixst &tet,
    matrixd &fixed_frame, matrixst &fixed_frame_idx,
    matrixd &aligned, matrixst &aligned_idx,
    const char *surf_node_uv_file,
    const char *surf_node2tet_node_file, const char *obj_file)
{
  unique_ptr<jtf::mesh::face2tet_adjacent> fa(jtf::mesh::face2tet_adjacent::create(tet));
  if(!fa.get()) {
      cerr << "createjtf::mesh::face2tet_adjacent fail." << endl;
    }

  //  if(validate_face2tet_adjacent(*fa, tet, &node)) {
  //      cerr << "face2tet_adjacent is incorrect." << endl;
  //      return __LINE__;
  //    }

  // align by surface normal
  matrixst outside_face;
  get_outside_face(*fa, outside_face);
  cerr << "# number of surface node: " << outside_face.size(2) << endl;

  get_outside_face_idx(*fa, aligned_idx);
  aligned.resize(3, outside_face.size(2));
  for(size_t i = 0; i < aligned_idx.size(); ++i) {
      size_t face_id = fa->get_face_idx(&outside_face(0, i));
      if(aligned_idx[i] != face_id) {
          cerr << "load_from_tet strange error." << endl;
          return 1;
        }
      if(aligned_idx[i] == -1) {
          cerr << "surf tri error." << endl;
        }
      //    if(normal(node(colon(), outside_face(colon(), i)), &aligned(0, i)) < 1e-8) {
      //      cerr << "degenerated surface triangle for tet: " <<  aligned_idx[i] << endl;
      //    }
    }
  jtf::mesh::cal_face_normal(outside_face, node,aligned);

  // load fix frame from file, it contains a pair of vectors in the
  // face for each surface triangle (in the format of Chunfeng).
  if(surf_node_uv_file == 0 || surf_node2tet_node_file == 0 || obj_file == 0) {
      cerr << "# provide uv, s2v and the obj file for surface fixed frame." << endl;
    }
  else {
      const size_t face_num = outside_face.size(2);
      fixed_frame.resize(3, face_num);
      fixed_frame_idx.resize(face_num);
      if(surface_frame2volme_frame_inner(
           surf_node_uv_file, surf_node2tet_node_file, obj_file,
           *fa, fixed_frame, fixed_frame_idx))
        return __LINE__;
    }

  cerr << "# load_from_tet success." << endl;
  return 0;
}
