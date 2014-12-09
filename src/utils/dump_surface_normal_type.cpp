#include <fstream>
#include <jtflib/mesh/io.h>
#include <jtflib/mesh/mesh.h>
#include <jtflib/mesh/util.h>

#include "../common/vtk.h"

using namespace std;
using namespace jtf::mesh;
using namespace zjucad::matrix;

int dump_surface_normal_type(int argc, char * argv[])
{
  if(argc != 3){
      cerr << "# [usage] dump_surface_normal_type input_obj output_vtk" << endl;
      return __LINE__;
    }

  meshes tm;
  if(jtf::mesh::load_obj(argv[1], tm.mesh_, tm.node_))
    return __LINE__;

  matrix<double> face_normal(3, tm.mesh_.size(2));
  jtf::mesh::cal_face_normal(tm.mesh_, tm.node_, face_normal);

  matrix<size_t> face_type(tm.mesh_.size(2),1);
  const matrix<double> eye_mat = eye<double>(3);
  vector<pair<double,int> > dir_error(6);
  for(size_t fi = 0; fi < face_type.size(); ++fi){
      for(size_t i = 0; i < 3; ++i){
          dir_error[2*i+0] = make_pair(
                dot(face_normal(colon(),fi), eye_mat(colon(),i)),
                2*i+0);
          dir_error[2*i+1] = make_pair(
                -1*dot(face_normal(colon(),fi), eye_mat(colon(),i)),
                2*i+1);
        }
      sort(dir_error.begin(), dir_error.end());
      face_type[fi] = dir_error.back().second/2;
    }

  ofstream ofs(argv[2]);
  if(tm.mesh_.size(1) == 3)
    tri2vtk(ofs, &tm.node_[0], tm.node_.size(2), &tm.mesh_[0], tm.mesh_.size(2));
  if(tm.mesh_.size(1) == 4)
    quad2vtk(ofs, &tm.node_[0], tm.node_.size(2), &tm.mesh_[0], tm.mesh_.size(2));
  cell_data(ofs, &face_type[0], face_type.size(), "face_type");

  return 0;
}
