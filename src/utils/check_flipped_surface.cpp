#include <iostream>
#include <fstream>

#include <jtflib/mesh/mesh.h>
#include <jtflib/mesh/io.h>
#include <jtflib/mesh/util.h>

#include "../tetmesh/tetmesh.h"
#include "../tetmesh/util.h"

#include "../common/vtk.h"

using namespace std;
using namespace zjucad::matrix;

int check_flipped_surface(int argc, char *argv[])
{
  if(argc != 2){
    cerr << "# [usage] check_flipped_surface tet" << endl;
    return __LINE__;
  }

  jtf::mesh::meshes tm;
  if(jtf::mesh::tet_mesh_read_from_zjumat(argv[1], &tm.node_, &tm.mesh_))
    return __LINE__;

  unique_ptr<jtf::mesh::face2tet_adjacent> fa(jtf::mesh::face2tet_adjacent::create(tm.mesh_));
  if(!fa.get()){
    cerr << "# [error] can not buildjtf::mesh::face2tet_adjacent." << endl;
    return __LINE__;
  }

  matrixst outside_face;
  get_outside_face(*fa, outside_face);

  unique_ptr<jtf::mesh::edge2cell_adjacent> ea(
        jtf::mesh::edge2cell_adjacent::create(outside_face));
  if(!ea.get()){
    cerr << "# [error] can not build edge2cell_adjacent." << endl;
    return __LINE__;
  }

  vector<pair<size_t,size_t> > face_pairs;
  vector<size_t> face_pair_idx;

  matrixd face_normal = zjucad::matrix::zeros<double>(3, outside_face.size(2));
  jtf::mesh::cal_face_normal(outside_face, tm.node_, face_normal);

  size_t fi = 0;
  for(size_t ei = 0; ei < ea->edge2cell_.size(); ++ei){
    const pair<size_t,size_t> & tri_pair = ea->edge2cell_[ei];
    const double dot_val = dot(face_normal(colon(), tri_pair.first),
                               face_normal(colon(), tri_pair.second));
    if(fabs(dot_val + 1) < 1e-5){
      face_pairs.push_back(tri_pair);
      face_pair_idx.push_back(fi);
      face_pair_idx.push_back(fi);
      ++fi;
    }
  }

  {
    ofstream ofs("flipped_face.vtk");
    vector<size_t> faces;
    for(size_t fi = 0; fi < face_pairs.size(); ++fi){
      const pair<size_t,size_t> & one_pair = face_pairs[fi];
      faces.insert(faces.end(), outside_face(colon(), one_pair.first).begin(),
                   outside_face(colon(), one_pair.first).end());
      faces.insert(faces.end(), outside_face(colon(), one_pair.second).begin(),
                   outside_face(colon(), one_pair.second).end());
    }
    tri2vtk(ofs, &tm.node_[0], tm.node_.size(2), &faces[0], faces.size()/3) ;
    cell_data(ofs, &face_pair_idx[0], face_pair_idx.size(), "idx");
  }

  return 0;
}
