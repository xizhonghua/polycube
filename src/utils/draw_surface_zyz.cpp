#include <jtflib/mesh/mesh.h>
#include <jtflib/mesh/io.h>
#include <jtflib/mesh/util.h>
#include "../common/zyz.h"
#include "../common/visualize_tool.h"
using namespace std;
using namespace jtf::mesh;
using namespace zjucad::matrix;

int draw_surface_zyz(int argc, char * argv[])
{
  if(argc != 4){
      cerr << "# [usage] draw_surface_zyz tet zyz arrow_obj" << endl;
      return __LINE__;
    }
  jtf::mesh::meshes tetmesh;
  if(jtf::mesh::tet_mesh_read_from_zjumat(argv[1], &tetmesh.node_, &tetmesh.mesh_))
    return __LINE__;
  matrix<double> zyz;
  if(jtf::mesh::read_matrix(argv[2], zyz))
    return __LINE__;

  matrix<matrix<double> > frame(zyz.size(2),1);
  for(size_t ti = 0; ti < zyz.size(2); ++ti){
      frame[ti].resize(3,3);
      zyz_angle_2_rotation_matrix1(&zyz(0,ti), &frame[ti][0]);
    }

  matrix<size_t> outside_face,outside_face_idx;
  unique_ptr<jtf::mesh::face2tet_adjacent> fa(jtf::mesh::face2tet_adjacent::create(tetmesh.mesh_));
  if(!fa.get()){
      cerr << "# [error] can not build face2tet_adjacent." << endl;
      return __LINE__;
    }
  get_outside_face(*fa, outside_face, true);
  get_outside_face_idx(*fa, outside_face_idx);
  matrix<double> outside_face_normal(3,outside_face.size(2));
  matrix<double> center_points = zeros<double>(3, outside_face.size(2)) ;
  vector<pair<double,size_t> > normal_axis_erro(6);
  matrix<double> next_direction1(3, outside_face.size(2));
  matrix<double> next_direction2(3, outside_face.size(2));

  double total_edge_len = 0;
  size_t edge_num = 0;
  matrix<double> one_face_normal(3,1);
  for(size_t fi = 0; fi < outside_face.size(2); ++fi){
      jtf::mesh::cal_face_normal(
            tetmesh.node_(colon(), outside_face(colon(),fi)),one_face_normal);

      outside_face_normal(colon(),fi) = one_face_normal;
      for(size_t pi = 0; pi < outside_face.size(1); ++pi){
          center_points(colon(),fi) += tetmesh.node_(colon(), outside_face(pi,fi));
          total_edge_len += norm(tetmesh.node_(colon(), outside_face(pi,fi))
                                 - tetmesh.node_(colon(),outside_face((pi+1)%3,fi)));
        }
      edge_num += 3;
      const pair<size_t,size_t> & tet_pair = fa->face2tet_[outside_face_idx[fi]];
      assert(tet_pair.first == -1 || tet_pair.second == -1);
      const size_t this_tet = (tet_pair.first==-1?tet_pair.second:tet_pair.first);
      for(size_t i = 0; i < 3; ++i){
          normal_axis_erro[2 * i + 0] =
              make_pair(dot(frame[this_tet](colon(), i), outside_face_normal(colon(),fi)),i);
          normal_axis_erro[2 * i + 1] =
              make_pair(dot(-1*frame[this_tet](colon(), i), outside_face_normal(colon(),fi)),i);
        }
      sort(normal_axis_erro.begin(), normal_axis_erro.end());

      const size_t choosed_axis = normal_axis_erro.back().second;
      next_direction1(colon(), fi) = frame[this_tet](colon(), (choosed_axis+1)%3);
      next_direction2(colon(), fi) = frame[this_tet](colon(), (choosed_axis+2)%3);
    }
  center_points /= 3.0;

  {
    vector<double> center_point_vec;
    center_point_vec.insert(center_point_vec.end(),
                            center_points.begin(), center_points.end());
    vector<size_t> lines;
    matrix<double> new_point;
    const double step = 0.5*total_edge_len / edge_num ;
    for(size_t fi = 0 ; fi < next_direction1.size(2); ++fi){
        lines.push_back(fi);
        lines.push_back(fi+next_direction1.size(2));
        new_point = center_points(colon(),fi) + step * next_direction1(colon(),fi);
        center_point_vec.insert(center_point_vec.end(), new_point.begin(), new_point.end());
      }
    directional_edge2arrow(&lines[0], lines.size()/2,
        &center_point_vec[0], center_point_vec.size()/3, argv[3], "axis_cross1.obj");
    cerr << "# lines " << lines.size()/2 << endl;
    cerr << "# points " << center_point_vec.size()/3 << endl;
  }

  {
    vector<double> center_point_vec;
    center_point_vec.insert(center_point_vec.end(),
                            center_points.begin(), center_points.end());
    vector<size_t> lines;
    matrix<double> new_point;
    const double step = 0.5*total_edge_len / edge_num ;
    for(size_t fi = 0 ; fi < next_direction1.size(2); ++fi){
        lines.push_back(fi);
        lines.push_back(fi+next_direction1.size(2));
        new_point = center_points(colon(),fi) + step * next_direction2(colon(),fi);
        center_point_vec.insert(center_point_vec.end(), new_point.begin(), new_point.end());
      }
    directional_edge2arrow(&lines[0], lines.size()/2,
         &center_point_vec[0], center_point_vec.size()/3, argv[3], "axis_cross2.obj");
    cerr << "# lines " << lines.size()/2 << endl;
    cerr << "# points " << center_point_vec.size()/3 << endl;
  }





  return 0;
}
