#include <iostream>
#include "../tetmesh/tetmesh.h"
#include "../common/IO.h"
#include "../common/zyz.h"
#include "../common/visualize_tool.h"
#include "../common/vtk.h"
using namespace std;
using namespace zjucad::matrix;

int draw_frame_axis(int argc, char * argv[])
{
  if(argc != 4){
      cerr << "# [usage] draw_frame_axis tet zyz arrow_obj" << endl;
      cerr << "# [warning] input frame should be globally aligned" << endl;
      return __LINE__;
    }

  jtf::tet_mesh tm(argv[1]);
  matrix<double> zyz;
  if(read_zyz(argv[2], zyz))
    return __LINE__;

  matrix<matrix<double> > frame(zyz.size(2),1);
  for(size_t ti = 0; ti < zyz.size(2); ++ti){
      frame[ti].resize(3,3);
      zyz_angle_2_rotation_matrix1(&zyz(0,ti), &frame[ti][0]);
    }

  const double avg_len = 0.5*jtf::mesh::cal_average_edge(tm.tetmesh_.mesh_, tm.tetmesh_.node_);
  matrix<double> center_point = zeros<double>(3, tm.tetmesh_.mesh_.size(2));

  for(size_t ti = 0; ti < tm.tetmesh_.mesh_.size(2); ++ti){
      for(size_t pi = 0; pi < tm.tetmesh_.mesh_.size(1); ++pi){
          center_point(colon(), ti) += tm.tetmesh_.node_(colon(), tm.tetmesh_.mesh_(pi,ti));
        }
    }
  center_point /= 4.0;

  matrix<double> dir_node0(3,tm.tetmesh_.mesh_.size(2)*2),
      dir_node1(3, tm.tetmesh_.mesh_.size(2)*2),
      dir_node2(3, tm.tetmesh_.mesh_.size(2)*2);

  matrix<size_t>  dir_edge(2, tm.tetmesh_.mesh_.size(2));
  for(size_t ti = 0; ti < tm.tetmesh_.mesh_.size(2); ++ti){
      dir_edge(0,ti) = 2*ti;
      dir_edge(1,ti) = 2*ti+1;
      dir_node0(colon(), 2*ti) = center_point(colon(),ti);
      dir_node0(colon(), 2*ti+1) = center_point(colon(),ti) + avg_len * frame[ti](colon(),0);

      dir_node1(colon(), 2*ti) = center_point(colon(),ti);
      dir_node1(colon(), 2*ti+1) = center_point(colon(),ti) + avg_len * frame[ti](colon(),1);
      dir_node2(colon(), 2*ti) = center_point(colon(),ti);
      dir_node2(colon(), 2*ti+1) = center_point(colon(),ti) + avg_len * frame[ti](colon(),2);

    }

  directional_edge2arrow(&dir_edge[0], dir_edge.size(2), &dir_node0[0] , dir_node0.size(2), argv[3], "u_dir.obj");
  directional_edge2arrow(&dir_edge[0], dir_edge.size(2), &dir_node1[0] , dir_node1.size(2), argv[3], "v_dir.obj");
  directional_edge2arrow(&dir_edge[0], dir_edge.size(2), &dir_node2[0] , dir_node2.size(2), argv[3], "w_dir.obj");
  return 0;
}


//int draw_frame_axis_around_singularity(int argc, char * argv[])
//{
//  if(argc != 5){
//      cerr << "# [usage] draw_frame_axis_around_singularity tet zyz arrow_obj [arround_singularity?1:yes;0:no]" << endl;
//      cerr << "# [warning] input frame should be globally aligned" << endl;
//      return __LINE__;
//    }

//  jtf::tet_mesh tm(argv[1]);
//  matrix<double> zyz;
//  if(read_zyz(argv[2], zyz))
//    return __LINE__;

//  matrix<matrix<double> > frame(zyz.size(2),1);
//  for(size_t ti = 0; ti < zyz.size(2); ++ti){
//      frame[ti].resize(3,3);
//      zyz_angle_2_rotation_matrix1(&zyz(0,ti), &frame[ti][0]);
//    }

//  vector<deque<pair<size_t,size_t> > > chain_list;
//  vector<deque<size_t> > singularities_type;
//  vector<vector<size_t> > singularities_tet_loop;
//  singularity_extractor se(tm);

//  find_singularities_with_frame(tm.ortae_, frame, tm.outside_face_, chain_list, singularities_type, singularities_tet_loop);

//  set<size_t> tets_need_to_show;
//  {
//    for(size_t i = 0; i < singularities_tet_loop.size(); ++i){
//        for(size_t p = 0; p < singularities_tet_loop[i].size(); ++p){
//            tets_need_to_show.insert(singularities_tet_loop[i][p]);
//          }
//      }
//  }

//  const double avg_len = 0.5*jtf::mesh::cal_average_edge(tm.tetmesh_.mesh_, tm.tetmesh_.node_);
//  matrix<double> center_point = zeros<double>(3, tets_need_to_show.size());

//  size_t ti = 0;
//  for(const auto & one_tet : tets_need_to_show){
//      for(size_t pi = 0; pi < tm.tetmesh_.mesh_.size(1); ++pi){
//          center_point(colon(), ti) += tm.tetmesh_.node_(colon(), tm.tetmesh_.mesh_(pi,one_tet));
//        }
//      ++ti;
//    }
//  center_point /= 4.0;

//  matrix<double> dir_node0(3,tm.tetmesh_.mesh_.size(2)*2),
//      dir_node1(3, tm.tetmesh_.mesh_.size(2)*2),
//      dir_node2(3, tm.tetmesh_.mesh_.size(2)*2);

//  matrix<size_t>  dir_edge(2, tm.tetmesh_.mesh_.size(2));
//  ti = 0;
//  for(const auto & one_tet : tets_need_to_show){
//      dir_edge(0,ti) = 2*ti;
//      dir_edge(1,ti) = 2*ti+1;
//      dir_node0(colon(), 2*ti) = center_point(colon(),ti);
//      dir_node0(colon(), 2*ti+1) = center_point(colon(),ti) + avg_len * frame[one_tet](colon(),0);

//      dir_node1(colon(), 2*ti) = center_point(colon(),ti);
//      dir_node1(colon(), 2*ti+1) = center_point(colon(),ti) + avg_len * frame[one_tet](colon(),1);
//      dir_node2(colon(), 2*ti) = center_point(colon(),ti);
//      dir_node2(colon(), 2*ti+1) = center_point(colon(),ti) + avg_len * frame[one_tet](colon(),2);
//      ++ti;
//    }

//  directional_edge2arrow(&dir_edge[0], dir_edge.size(2), &dir_node0[0] , dir_node0.size(2), argv[3], "u_dir.obj");
//  directional_edge2arrow(&dir_edge[0], dir_edge.size(2), &dir_node1[0] , dir_node1.size(2), argv[3], "v_dir.obj");
//  directional_edge2arrow(&dir_edge[0], dir_edge.size(2), &dir_node2[0] , dir_node2.size(2), argv[3], "w_dir.obj");
//  return 0;
//}


int draw_frame_axis_on_surface(int argc, char * argv[])
{
  if(argc != 4){
      cerr << "# [usage] draw_frame_axis_on_surface tet zyz arrow_obj" << endl;
      cerr << "# [warning] input frame should be globally aligned" << endl;
      return __LINE__;
    }

  jtf::tet_mesh tm(argv[1]);
  matrix<double> zyz;
  if(read_zyz(argv[2], zyz))
    return __LINE__;

  matrix<matrix<double> > frame(zyz.size(2),1);
  for(size_t ti = 0; ti < zyz.size(2); ++ti){
      frame[ti].resize(3,3);
      zyz_angle_2_rotation_matrix1(&zyz(0,ti), &frame[ti][0]);
    }

  const double avg_len = 0.5*jtf::mesh::cal_average_edge(tm.tetmesh_.mesh_, tm.tetmesh_.node_);
  matrix<double> center_point = zeros<double>(3, tm.outside_face_.size(2));

  for(size_t fi = 0; fi < tm.outside_face_.size(2); ++fi){
      for(size_t pi = 0; pi < tm.outside_face_.size(1); ++pi){
          center_point(colon(),fi) += tm.tetmesh_.node_(colon(), tm.outside_face_(pi,fi));
        }
    }
  center_point /= 3.0;

  matrix<size_t> surface_tet_idx(tm.outside_face_.size(2),1);
  for(size_t i = 0; i < tm.outside_face_idx_.size(); ++i){
      const pair<size_t,size_t> & tet_pair = tm.fa_->face2tet_[tm.outside_face_idx_[i]];
      if(tet_pair.first == -1) surface_tet_idx[i]=tet_pair.second;
      if(tet_pair.second == -1) surface_tet_idx[i]=tet_pair.first;
    }

  matrix<double> dir_node0(3,tm.tetmesh_.mesh_.size(2)*2),
      dir_node1(3, tm.tetmesh_.mesh_.size(2)*2),
      dir_node2(3, tm.tetmesh_.mesh_.size(2)*2);

  matrix<size_t>  dir_edge(2, tm.tetmesh_.mesh_.size(2));
  size_t ti = 0;
  for(const auto & one_tet : surface_tet_idx){
      dir_edge(0,ti) = 2*ti;
      dir_edge(1,ti) = 2*ti+1;
      dir_node0(colon(), 2*ti) = center_point(colon(),ti);
      dir_node0(colon(), 2*ti+1) = center_point(colon(),ti) + avg_len * frame[one_tet](colon(),0);

      dir_node1(colon(), 2*ti) = center_point(colon(),ti);
      dir_node1(colon(), 2*ti+1) = center_point(colon(),ti) + avg_len * frame[one_tet](colon(),1);
      dir_node2(colon(), 2*ti) = center_point(colon(),ti);
      dir_node2(colon(), 2*ti+1) = center_point(colon(),ti) + avg_len * frame[one_tet](colon(),2);
      ++ti;
    }

  directional_edge2arrow(&dir_edge[0], dir_edge.size(2), &dir_node0[0] , dir_node0.size(2), argv[3], "u_dir.obj");
  directional_edge2arrow(&dir_edge[0], dir_edge.size(2), &dir_node1[0] , dir_node1.size(2), argv[3], "v_dir.obj");
  directional_edge2arrow(&dir_edge[0], dir_edge.size(2), &dir_node2[0] , dir_node2.size(2), argv[3], "w_dir.obj");
  return 0;
}
