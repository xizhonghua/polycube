#include <iostream>

#include <jtflib/mesh/mesh.h>
#include <jtflib/mesh/io.h>
#include <jtflib/mesh/util.h>

#include "../tetmesh/util.h"
#include "../tetmesh/hex_io.h"
#include "../common/IO.h"
#include "../common/zyz.h"
#include "../common/def.h"

using namespace std;
using namespace zjucad::matrix;

int dump_surf_zyz2fv(int argc, char * argv[])
{
  if(argc != 6 && argc != 4){
      cerr << "# [usage1] dump_surf_zyz2fv tet zyz obj s2v fv" << endl;
      cerr << "# [usage2] dump_surf_zyz2fv zyz obj fv" << endl;
      return __LINE__;
    }

  if(argc == 6){
      matrixst tet;
      matrixd node, zyz;
      if(jtf::mesh::tet_mesh_read_from_zjumat(argv[1], &node, & tet)){
          cerr << "# [error] can not read tet file." << endl;
          return __LINE__;
        }

      jtf::mesh::read_matrix(argv[2], zyz);
      if(zyz.size(2) != tet.size(2)){
          cerr << "# [error] wrong zyz file" << endl;
          return __LINE__;
        }

      unique_ptr<jtf::mesh::face2tet_adjacent> fa(jtf::mesh::face2tet_adjacent::create(tet));
      if(!fa.get()){
          cerr << "# [error] can not buildjtf::mesh::face2tet_adjacent" << endl;
          return __LINE__;
        }


      map<size_t, matrixd> face_to_normal;
      {
        matrixd face_normal;
        matrixst outside_face_idx, outside_face;
        get_outside_face_idx(*fa, outside_face_idx);
        get_outside_face(*fa, outside_face);
        jtf::mesh::cal_face_normal(outside_face, node, face_normal);
        jtf::tetmesh::orient_face_normal_outside_tetmesh(
              tet, node, outside_face, outside_face_idx, *fa, face_normal);

        for(size_t fi = 0 ; fi < outside_face_idx.size(); ++fi){
            face_to_normal[outside_face_idx[fi]] = face_normal(colon(), fi);
          }
      }


      zjucad::matrix::matrix<matrixd> frame(zyz.size(2));
      for(size_t fi = 0; fi < frame.size(); ++fi){
          frame[fi].resize(3,3);
          zyz_angle_2_rotation_matrix1(&zyz(0,fi), &frame[fi][0]);
        }

      matrix<int32_t> surf_node2tet_node; {
        ifstream s2v_ifs(argv[4], ifstream::binary);
        if(s2v_ifs.fail()) {
            cerr << "# open " << argv[4] << " fail." << endl;
            return __LINE__;
          }
        jtf::mesh::read_matrix(s2v_ifs, surf_node2tet_node);
      }

      jtf::mesh::meshes trim;
      if(jtf::mesh::load_obj(argv[3], trim.mesh_, trim.node_)){
          cerr << "# [error] can not load obj file." << endl;
          return __LINE__;
        }


      vector<matrixd> face_vector(trim.mesh_.size(2));
      vector<pair<double, size_t> > face_normal_alignment(6);
      for(size_t fi = 0; fi < trim.mesh_.size(2); ++fi){
          matrixst orig_surface = surf_node2tet_node(trim.mesh_(colon(),fi));
          const size_t & face_idx = fa->get_face_idx(&orig_surface[0]);
          const pair<size_t,size_t> & tet_pair = fa->face2tet_[face_idx];
          const size_t & tet_idx =
              (tet_pair.first == -1?tet_pair.second:tet_pair.first);
          for(size_t ri = 0; ri < 3; ++ri){
              face_normal_alignment[2 * ri + 0] =
                  make_pair(dot(face_to_normal[face_idx],frame[tet_idx](colon(), ri)),
                  2*ri+0);
              face_normal_alignment[2 * ri + 1] =
                  make_pair(dot(face_to_normal[face_idx],-1* frame[tet_idx](colon(), ri)),
                  2*ri+1);
            }
          sort(face_normal_alignment.begin(), face_normal_alignment.end());
          const size_t aligned_axis = face_normal_alignment.back().second / 2;

          face_vector[fi].resize(3,2);
          face_vector[fi](colon(),0) = frame[tet_idx](colon(), (aligned_axis+1)%3);
          face_vector[fi](colon(),1) = frame[tet_idx](colon(), (aligned_axis+2)%3);

        }

      cerr << "# [info] dump out fv." << endl;
      {
        ofstream ofs(argv[5]);
        if(ofs.fail()){
            cerr << "# [error] can not open output fv file" << endl;
            return __LINE__;
          }
        ofs << face_vector.size() << endl;
        for(size_t fi = 0; fi < face_vector.size(); ++fi){
            ofs << endl;
            for(size_t vi = 0; vi < face_vector[fi].size(2); ++vi){
                ofs << face_vector[fi](0,vi) << " "
                                             << face_vector[fi](1,vi) << " "
                                                                      <<  face_vector[fi](2,vi) << endl;
              }
          }
      }
    }else if(argc == 4){
      jtf::mesh::meshes trim;
      matrixd zyz;

      if(jtf::mesh::load_obj(argv[2], trim.mesh_, trim.node_))
        return __LINE__;

      jtf::mesh::read_matrix(argv[1], zyz);
      if(zyz.size(2) != trim.mesh_.size(2)){
          cerr << "# [error] wrong zyz file" << endl;
          return __LINE__;
        }

      matrix<double> face_normal(3, trim.mesh_.size(2));
      jtf::mesh::cal_face_normal(trim.mesh_ ,trim.node_,face_normal);

      zjucad::matrix::matrix<matrixd> frame(zyz.size(2));
      for(size_t fi = 0; fi < frame.size(); ++fi){
          frame[fi].resize(3,3);
          zyz_angle_2_rotation_matrix1(&zyz(0,fi), &frame[fi][0]);
        }

      vector<matrixd> face_vector(trim.mesh_.size(2));
      vector<pair<double, size_t> > face_normal_alignment(6);
      for(size_t fi = 0; fi < trim.mesh_.size(2); ++fi){
          for(size_t ri = 0; ri < 3; ++ri){
              face_normal_alignment[2 * ri + 0] =
                  make_pair(dot(face_normal(colon(),fi),frame[fi](colon(), ri)),
                  2*ri+0);
              face_normal_alignment[2 * ri + 1] =
                  make_pair(dot(face_normal(colon(),fi),-1.0* frame[fi](colon(), ri)),
                  2*ri+1);
            }
          sort(face_normal_alignment.begin(), face_normal_alignment.end());
          const size_t aligned_axis = face_normal_alignment.back().second / 2;

          face_vector[fi].resize(3,2);
          face_vector[fi](colon(),0) = frame[fi](colon(), (aligned_axis+1)%3);
          face_vector[fi](colon(),1) = frame[fi](colon(), (aligned_axis+2)%3);
        }

      cerr << "# [info] dump out fv." << endl;
      {
        ofstream ofs(argv[3]);
        if(ofs.fail()){
            cerr << "# [error] can not open output fv file" << endl;
            return __LINE__;
          }
        ofs << face_vector.size() << endl;
        for(size_t fi = 0; fi < face_vector.size(); ++fi){
            ofs << endl;

            for(size_t vi = 0; vi < face_vector[fi].size(2); ++vi){
                ofs << face_vector[fi](0,vi)
                       << " "
                       << face_vector[fi](1,vi)
                          << " "
                          <<  face_vector[fi](2,vi) << endl;
              }
          }
      }
    }

  return 0;
}
