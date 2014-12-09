#include <jtflib/mesh/io.h>
#include <jtflib/mesh/mesh.h>
#include <jtflib/mesh/util.h>
#include <iostream>
#include "../common/vtk.h"
#include "../numeric/util.h"
using namespace std;
using namespace zjucad::matrix;

int angle_defect_fabs(int argc, char *argv[])
{
  if(argc != 3){
      cerr << "# [usage] angle_defect_fabs input_obj output_size_field" << endl;
      return __LINE__;
    }
  jtf::mesh::meshes trimesh;
  if(jtf::mesh::load_obj(argv[1], trimesh.mesh_, trimesh.node_))
    return __LINE__;

  shared_ptr<jtf::mesh::edge2cell_adjacent> ea(jtf::mesh::edge2cell_adjacent::create(trimesh.mesh_));
  if(!ea.get()){
      cerr << "# [error] can not build edge2cell_adjacent." << endl;
      return __LINE__;
    }

  jtf::mesh::one_ring_face_at_point orfap;
  orfap.add_all_faces(trimesh.mesh_,*ea);
  orfap.sort_int_loop(trimesh.mesh_, trimesh.node_);
  matrix<double> size = zeros<double>(trimesh.node_.size(2),1);

  for(jtf::mesh::one_ring_face_at_point::p2f_type::const_iterator it = orfap.p2f_.begin();
      it != orfap.p2f_.end(); ++it){
      const vector<size_t> & one_ring_face = it->second;
      if(one_ring_face.front() == -1 || one_ring_face.back() == -1) {
          size[it->first] = 0;
          continue;
        }
      double angle = 0;
      vector<double> angle_vec(3);
      for(const auto & one_face : one_ring_face){
          jtf::mesh::cal_face_angle(trimesh.mesh_(colon(), one_face), trimesh.node_, angle_vec);

          size_t idx = -1;
          for(size_t i = 0; i < trimesh.mesh_.size(1); ++i)
            if(trimesh.mesh_(i, one_face) == it->first)
              idx = i;
          angle += fabs(angle_vec[idx])/180.0*My_PI();
        }
      size[it->first] = 2*My_PI() - angle;
    }

  ofstream ofs(argv[2]);
  ofs << size.size() << endl;
  for(size_t i = 0; i < size.size(); ++i)
    ofs << size[i] << endl;

  ofstream ofs_vtk("size_field.vtk");
  tri2vtk(ofs_vtk, &trimesh.node_[0], trimesh.node_.size(2), &trimesh.mesh_[0], trimesh.mesh_.size(2));
  point_data(ofs_vtk, &size[0], size.size(), "size_field");
  return 0;
}
