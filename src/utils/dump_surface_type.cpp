#include <iostream>
#include <fstream>

#include <jtflib/mesh/io.h>
#include <jtflib/mesh/util.h>
#include "../tetmesh/tetmesh.h"
#include "../tetmesh/hex_io.h"
#include "../common/vtk.h"
#include "../hex_param/io.h"
#include "../common/zyz.h"
#include "../common/transition_type.h"

using namespace std;
using namespace zjucad::matrix;


int dump_surface_type(int argc, char * argv[])
{
  if(argc != 4){
      cerr << "# [info] useage: dump_surface_type tet surface_type output_vtk" << endl;
      return __LINE__;
    }

  jtf::mesh::meshes tm;
  if(jtf::mesh::tet_mesh_read_from_zjumat(argv[1],&tm.node_,&tm.mesh_)){
      cerr << "# [error] can not read tet file." << endl;
      return __LINE__;
    }

  unique_ptr<jtf::mesh::face2tet_adjacent> fa(jtf::mesh::face2tet_adjacent::create(tm.mesh_));
  if(!fa.get()){
      cerr << "# [error] can not buildjtf::mesh::face2tet_adjacent." << endl;
      return __LINE__;
    }

  boost::unordered_map<size_t,size_t> surface_type;
  int rtn = load_surface_type(argv[2], surface_type, fa.get());
  if(rtn != 0 && rtn != 1)
    return __LINE__;

  matrix<size_t> outside_face, outside_face_idx;
  get_outside_face(*fa, outside_face);
  get_outside_face_idx(*fa, outside_face_idx);


  for(size_t fi = 0; fi < outside_face_idx.size(); ++fi){
      boost::unordered_map<size_t,size_t>::const_iterator cit =
          surface_type.find(outside_face_idx[fi]);
      if(cit == surface_type.end()){
          cerr << "# [error] can not find surface_type of " << outside_face_idx[fi] << endl;
          return __LINE__;
        }
      outside_face_idx[fi] = cit->second;
    }

  ofstream ofs(argv[3]);
  tri2vtk(ofs, &tm.node_[0],tm.node_.size(2),&outside_face[0],outside_face.size(2));
  cell_data(ofs, &outside_face_idx[0],outside_face_idx.size(),"surface_type");
  return 0;
}

int dump_surface_zyz_type(int argc, char * argv[])
{
  if(argc != 4){
      cerr << "# [info] useage: dump_surface_zyz_type tet zyz output_vtk" << endl;
      return __LINE__;
    }

  jtf::mesh::meshes tm;
  if(jtf::mesh::tet_mesh_read_from_zjumat(argv[1],&tm.node_,&tm.mesh_)){
      cerr << "# [error] can not read tet file." << endl;
      return __LINE__;
    }

  unique_ptr<jtf::mesh::face2tet_adjacent> fa(jtf::mesh::face2tet_adjacent::create(tm.mesh_));
  if(!fa.get()){
      cerr << "# [error] can not buildjtf::mesh::face2tet_adjacent." << endl;
      return __LINE__;
    }

  matrix<double> zyz(3,tm.mesh_.size(2));
  matrix<matrix<double> > frame(tm.mesh_.size(2),1);

  if(jtf::mesh::read_matrix(argv[2], zyz))
    return __LINE__;

  for(size_t ti = 0; ti < tm.mesh_.size(2); ++ti){
      frame[ti].resize(3,3);
      zyz_angle_2_rotation_matrix1(&zyz(0,ti), &frame[ti][0]);
    }

  matrix<size_t> outside_face, outside_face_idx;
  get_outside_face(*fa, outside_face,true);
  get_outside_face_idx(*fa, outside_face_idx);

  vector<pair<double,size_t> > resi(24);

  matrix<double> normal(3,1);
  matrix<double> temp(3,3);
  vector<size_t> surface_type;
  for(size_t t = 0; t < outside_face_idx.size(); ++t){
      const pair<size_t,size_t> &tet_pair = fa->face2tet_[outside_face_idx[t]];
      assert(tet_pair.second == -1 || tet_pair.first == -1);
      const size_t &tet_idx = (tet_pair.first == -1)?tet_pair.second:tet_pair.first;

      jtf::mesh::cal_face_normal(tm.node_(colon(), outside_face(colon(),t)),normal);

      for(size_t i = 0; i < 24; ++i){
          temp = frame[tet_idx]*type_transition2(i);
          const matrixd U_vec = temp(colon(),0);
          resi[i] = make_pair(dot(U_vec,normal),i);
        }
      sort(resi.begin(),resi.end());

      for(size_t i = 23; i != -1; --i){
          if(resi[i].second < 10){ // choose one axis rotaion
              surface_type.push_back(resi[i].second);
              break;
            }
        }
    }

  ofstream ofs(argv[3]);
  tri2vtk(ofs, &tm.node_[0],tm.node_.size(2),&outside_face[0],outside_face.size(2));
  cell_data(ofs, &surface_type[0],surface_type.size(),"surface_type");
  return 0;
}
