/*=============================================================================*\
 *                                                                             *
 *                               Volume Frame                                  *
 *         Copyright (C) 2011-2021 by State Key Lab of CAD&CG 410              *
 *			Zhejiang University, shizeyun                          *
 *                        http://www.cad.zju.edu.cn/                           *
 *                                                                             *
 *-----------------------------------------------------------------------------*
 *  This file is part of Volume Frame                                          *
 *  Created by shizeyun                                                        *
 *                                                                             *
 \*=============================================================================*/

//== INCLUDES ===================================================================

//STANDARD
#include<algorithm>
#include<cassert>
#include<iostream>
#include<fstream>
#include <memory>
#include<vector>
#include<zjucad/matrix/io.h>
#include<boost/unordered_set.hpp>

#include <jtflib/mesh/mesh.h>
#include <jtflib/util/util.h>

//LOCAL
#include "../common/util.h"
#include "../numeric/util.h"
#include "hex_process.h"
#include "../common/vtk.h"
#include <jtflib/util/container_operation.h>

//== NAMESPACES =================================================================

using zjucad::matrix::colon;
using zjucad::matrix::cross;
using zjucad::matrix::dot;
using zjucad::matrix::matrix;
using zjucad::matrix::max;
using zjucad::matrix::min;
using zjucad::matrix::sum;
using zjucad::matrix::norm;
using zjucad::matrix::zeros;
using jtf::mesh::edge2cell_adjacent;
using namespace std;

//===============================================================================

bool is_point_in_triangle(const matrixd &point,
                          const matrixd &triangle,
                          matrixd &weight_coord) {
  //area method
  matrixd tmp_matrix31(3, 1);
  double abc = norm(cross(&tmp_matrix31,
                          triangle(colon(), 0) - triangle(colon(), 1),
                          triangle(colon(), 0) - triangle(colon(), 2)));
  assert(abc > 1e-8);
  weight_coord = zeros<double>(3, 1);
  for (size_t i = 0; i < 3; ++i) {
    weight_coord[i] = norm(cross(&tmp_matrix31,
                                 point - triangle(colon(), (i+1)%3),
                                 point - triangle(colon(), (i+2)%3))) / abc;
  }
  return (fabs(sum(weight_coord) - 1.0) < 1e-8);
}

//-------------------------------------------------------------------------------

int create_map_from_hex_to_tet(
    const matrixd &tet_nodes,
    const matrixst &tet_faces,
    const matrixd &hex_nodes,
    const matrixst &hex_faces,
    const matrixst &hex,
    const matrixst & tet_surface_points,
    const jtf::mesh::face2hex_adjacent &fa,
    const jtf::mesh::edge2cell_adjacent &ea,
    const jtf::mesh::one_ring_face_at_point & ring_faces,
    const matrixst & hex_face_nodes,
    matrixd & nearest_nodes,
    matrixst & nearest_node_faces,
    matrixd & nearest_node_weights)
{
  //check input data
  if (0 == tet_nodes.size(2) || 0 == tet_faces.size(2) ||
      0 == hex_nodes.size(2) || 0 == hex_faces.size(2))
    return -1;

  std::cout << "# [info] node number of tet: " << tet_nodes.size(2) << std::endl;
  std::cout << "# [info] face number of tet: " << tet_faces.size(2) << std::endl;
  std::cout << "# [info] node number of hex: " << hex_nodes.size(2) << std::endl;
  std::cout << "# [info] face number of hex: " << hex_faces.size(2) << std::endl;
  
  std::cout << "# [info] face edge number of tet: "
            << ea.edges_.size() << std::endl;
  std::cout << "# [info] face node number of tet: "
            << tet_surface_points.size() << std::endl;
  std::cout << "# [info] face node number of hex: "
            << hex_face_nodes.size() << std::endl;

  //initial datas
  const size_t &tet_node_num = tet_surface_points.size();
  const size_t &hex_node_num = hex_face_nodes.size();
  const size_t &total_hex_node_num = hex_nodes.size(2);
  matrixd distances = zeros<double>(hex_node_num, 1);
  nearest_nodes = zeros<double>(3, hex_node_num);
  nearest_node_faces = zeros<size_t>(hex_node_num, 1);
  nearest_node_weights = zeros<double>(3, hex_node_num);

  //cal point to point
  matrixd tmp_matrix31(3, 1);
  for (size_t hid = 0; hid < hex_node_num; ++hid) {
    double dis = norm(hex_nodes(colon(), hex_face_nodes[hid])
                      - tet_nodes(colon(), tet_surface_points[0]));
    size_t choosed = tet_surface_points[0];
    for (size_t tid = 1; tid < tet_node_num; ++tid) {
      double tmp = norm(hex_nodes(colon(), hex_face_nodes[hid])
                        - tet_nodes(colon(), tet_surface_points[tid]));
      if (tmp < dis) {
        dis = tmp;
        choosed = tet_surface_points[tid];
      }
    }
    distances[hid] = dis;
    nearest_nodes(colon(), hid) = tet_nodes(colon(), choosed);
    //update face and weight info
    auto face_it = ring_faces.p2f_.find(choosed);
    size_t face_num = 0;
    if(face_it != ring_faces.p2f_.end())
      face_num = face_it->second.size();
    size_t i;
    for (i = 0; i < face_num; ++i)
      if (-1 != face_it->second[i]) {
        nearest_node_faces[hid] = face_it->second[i];
        for (size_t j = 0; j < 3; ++j)
          if (tet_faces(j, nearest_node_faces[hid]) == choosed) {
            nearest_node_weights(j, hid) = 1.0;
            break;
          }
        break;
      }
    if (i == face_num) {
      std::cerr << "can't find face for vertex: " << choosed << std::endl;
      return -1;
    }
  }
  
  matrixd result, p;
  
  //pre cal point to segment values
  const size_t tet_edge_num = ea.edges_.size();
  matrixd edge_directions = zeros<double>(3, tet_edge_num);
  matrixd a, b;
  for (size_t tid = 0; tid < tet_edge_num; ++tid) {
    a = tet_nodes(colon(), ea.edges_[tid].first);
    b = tet_nodes(colon(), ea.edges_[tid].second);
    double dis_ab = norm(a-b);
    assert(dis_ab > 1e-6);
    edge_directions(colon(), tid) = (a-b) / dis_ab;
  }
  
  //cal point to segment
  for (size_t hid = 0; hid < hex_node_num; ++hid) {
    p = hex_nodes(colon(), hex_face_nodes[hid]);
    double &dis = distances[hid];
    for (size_t tid = 0; tid < tet_edge_num; ++tid) {
      a = tet_nodes(colon(), ea.edges_[tid].first);
      b = tet_nodes(colon(), ea.edges_[tid].second);
      int d;
      for(d = 0; d < 3; ++d) {
        if((std::max(a[d], b[d]) + dis < p[d])
           || (std::min(a[d], b[d]) - dis > p[d]))
          break;
      }
      if(d < 3) continue;
      double t = dot((a-p), edge_directions(colon(), tid));
      if (t < -1e-6 || t > 1+1e-6)
        continue;
      result = a + t*(b-a);
      double tmp = norm(p-result);
      if (tmp < dis) {
        dis = tmp;
        nearest_nodes(colon(), hid) = result;

        //update face and weight info
        size_t fid = (-1 == ea.edge2cell_[tid].first)?
                       ea.edge2cell_[tid].first :
                       ea.edge2cell_[tid].second;
        if (-1 == fid) {
          std::cerr << "can't find face for edge: " << tid << std::endl;
          return -1;
        }
        nearest_node_faces[hid] = fid;
        for (size_t i = 0; i < 3; ++i) {
          nearest_node_weights(i, hid) = 0;
          if (tet_faces(i, fid) == ea.edges_[tid].first)
            nearest_node_weights(i, hid) = 1.0 - t;
          if (tet_faces(i, fid) == ea.edges_[tid].second)
            nearest_node_weights(i, hid) = t;
        }
      }
    }
  }
  
  //pre cal point to face values
  const size_t tet_face_num = tet_faces.size(2);
  matrixd normals = zeros<double>(3, tet_face_num);
  matrixd bounding_box = zeros<double>(6, tet_face_num);
  matrixd triangle, normal;
  for (size_t tid = 0; tid < tet_face_num; ++tid) {
    triangle = tet_nodes(colon(), tet_faces(colon(), tid));
    normal = cross(&tmp_matrix31,
                   triangle(colon(), 1) - triangle(colon(), 0),
                   triangle(colon(), 2) - triangle(colon(), 1));
    double norm_len = norm(normal);
    assert(norm_len > 1e-6);
    normals(colon(), tid) = normal / norm_len;

    for (size_t i = 0; i < 3; ++i) {
      bounding_box(i*2, tid) = max(triangle(i, colon()));
      bounding_box(i*2+1, tid) = min(triangle(i, colon()));
    }
  }

  //cal point to face
  matrixd weight_coord;
  for (size_t hid = 0; hid < hex_node_num; ++hid) {
    p = hex_nodes(colon(), hex_face_nodes[hid]);
    double &dis = distances[hid];
    for (size_t tid = 0; tid < tet_face_num; ++tid) {
      triangle = tet_nodes(colon(), tet_faces(colon(), tid));
      size_t d = 0;
      for (d = 0; d < 3; ++d) {
        if (bounding_box(d*2, tid) + dis < p[d]
            || bounding_box(d*2+1, tid) - dis > p[d])
          break;
      }
      if (d < 3) continue;
      double t = dot(normals(colon(), tid), triangle(colon(), 0)-p);
      result = p + t*normals(colon(), tid);
      double tmp = norm(p-result);
      if (tmp < dis) {
        if (is_point_in_triangle(result, triangle, weight_coord)) {
          dis = tmp;
          nearest_nodes(colon(), hid) = result;
          //update face and weight info
          nearest_node_faces[hid] = tid;
          nearest_node_weights(colon(), hid) = weight_coord;
        }
      }
    }
  }

  //test for code
  // for (size_t hid = 0; hid < hex_node_num; ++hid) {
  //   p = nearest_nodes(colon(), hid);
  //   triangle = tet_nodes(colon(), tet_faces(colon(), nearest_node_faces[hid]));
  //   weight_coord = nearest_node_weights(colon(), hid);
  //   if (!(max(weight_coord) > 1e-6))
  //     std::cerr << "ERRROR!!!!!!!!!!!!!!!!" << std::endl;
  //   if (!(norm(p-triangle*weight_coord) < 1e-6))
  //     std::cerr << "ERROR!!!!!!!!!!!!!!!ERROR" << std::endl;
  // }
  
  //check valid by write to vtk
#if 1
  {
    std::ofstream ofs("new_surface_node.vtk");
    std::vector<size_t> nearest_points;
    generate_N_number(nearest_points, nearest_nodes.size(2));
    point2vtk(ofs, &nearest_nodes[0], nearest_nodes.size(2),
              &nearest_points[0], nearest_points.size());

    //check valid by write to vtk
    std::ofstream outfile("quad.vtk");
    quad2vtk(outfile,
             &hex_nodes[0], hex_nodes.size(2),
             &hex_faces[0], hex_faces.size(2));

    //check valid by write to vtk
    std::ofstream outfile2("tri.vtk");
    tri2vtk(outfile2,
            &tet_nodes[0], tet_nodes.size(2),
            &tet_faces[0], tet_faces.size(2));
  }
#endif
  return 0;
}

static void cal_near_project_node(const matrixd & point,
                                  const matrixst & tri_feature,
                                  const matrixd & tri_node,
                                  matrixd & nearest_node)
{
  assert(point.size(1) == 3 && point.size(2) == 1);
  double min_dis = std::numeric_limits<double>::max();
  const matrixd &p = point;
  matrixd near_point = zeros<double>(3,1);
  matrixd foot = zeros<double>(3,1);

  for(size_t ti = 0; ti < tri_feature.size(2); ++ti){
    const matrixd p0 = tri_node(colon(),tri_feature(0,ti));
    const matrixd p1 = tri_node(colon(),tri_feature(1,ti));
    const double len_p0p1 = norm(p0-p1);

    if(norm(p0 - p) < norm(p1 - p) && norm(p0 - p) < min_dis){
      near_point = p0;
      min_dis = norm(p0 - p);
    }else if (norm(p1 - p) < norm(p0 - p) && norm(p1 - p) < min_dis){
      near_point = p1;
      min_dis = norm(p1 - p);
    }

    foot = p1 + dot((p1-p0)/norm(p1-p0), p-p1);
    if(norm(foot - p0) < len_p0p1 && norm(foot - p1) < len_p0p1){
      if(norm(foot-p) < min_dis){
        min_dis = norm(foot-p);
        near_point = foot;
      }
    }
  }
  nearest_node = near_point;
}

int project_quad_feature_to_tri_feature_new(
    const matrixst & quad_feature,
    const matrixd & quad_node,
    const matrixst & tri_feature,
    const matrixd & tri_node,
    boost::unordered_map<size_t, matrixd > & q2t)
{
  set<size_t> quad_fl_joints;
  {
    vector<pair<size_t,size_t> > quad_fl_edges(quad_feature.size(2));
    for(size_t ei = 0; ei < quad_feature.size(2); ++ei){
      quad_fl_edges[ei].first = quad_feature(0,ei);
      quad_fl_edges[ei].second = quad_feature(1,ei);
    }

    vector<deque<pair<size_t,size_t> > > chains;
    jtf::util::extract_chain_from_edges(quad_fl_edges, chains);
    for(size_t ci = 0; ci < chains.size(); ++ci){
      const deque<pair<size_t,size_t> > & one_deque = chains[ci];
      quad_fl_joints.insert(one_deque.front().first);
      quad_fl_joints.insert(one_deque.back().second);
    }
  }

  matrixd project_for_quad_points(3,quad_feature.size());
  matrixd nearest_point = zeros<double>(3,1);
  for(size_t qi = 0; qi < quad_feature.size(); ++qi){
    cal_near_project_node(quad_node(colon(),quad_feature[qi]), tri_feature, tri_node,
                          nearest_point);
    project_for_quad_points(colon(),qi) = nearest_point;
  }

  // create a threshold to determ which quad point should align tri feature
  // which is actually necessless points
  double avg_edge_length = 0;
  for(size_t qfli = 0; qfli < quad_feature.size(2); ++qfli){

    avg_edge_length += norm(quad_node(colon(), quad_feature(0,qfli))
                            -quad_node(colon(), quad_feature(1,qfli)));
  }
  avg_edge_length /= quad_feature.size(2);
  // for each quad feature line segment, is one point beyond the ave_edge_length,
  // then this line segment is no need.
  vector<bool> is_quad_fl_point_needed(quad_feature.size(),true);
  for(size_t qi = 0; qi < quad_feature.size(); ++qi){
    if(norm(project_for_quad_points(colon(),qi)
            - quad_node(colon(), quad_feature[qi]))
       > 2*avg_edge_length)
      is_quad_fl_point_needed[qi] = false;
  }

  // determin whether the quad feature line segment is parallel with
  // projected_edge segment
  const double cos_threshold = cos(60.0/180.0*My_PI());
  for(size_t qfi = 0; qfi < quad_feature.size(2); ++qfi){
    if(is_quad_fl_point_needed[2*qfi+0] &&
       is_quad_fl_point_needed[2*qfi+1]){
      const double proj_len =
          norm(project_for_quad_points(colon(),2*qfi + 0)
               -project_for_quad_points(colon(),2*qfi + 1));
      if(proj_len < 1e-6){
        is_quad_fl_point_needed[2*qfi+0] = false;
        is_quad_fl_point_needed[2*qfi+1] = false;
      }else{
        const matrixd proj_edge =
            (project_for_quad_points(colon(), 2*qfi+0)
             -project_for_quad_points(colon(), 2*qfi+1))/proj_len ;
        matrixd orig_edge =
            (quad_node(colon(), quad_feature(0,qfi))
             -quad_node(colon(), quad_feature(1,qfi)));
        orig_edge /= norm(orig_edge);
        const double cos_theta = dot(proj_edge, orig_edge);
        if(fabs(cos_theta) < cos_threshold){
          is_quad_fl_point_needed[2*qfi+0] = false;
          is_quad_fl_point_needed[2*qfi+1] = false;
        }
      }
    }
  }

  for(size_t qfi = 0; qfi < quad_feature.size(2); ++qfi){
    if(is_quad_fl_point_needed[2*qfi + 0] &&
       is_quad_fl_point_needed[2*qfi + 1]){
      q2t[quad_feature(0,qfi)] = project_for_quad_points(colon(), 2*qfi + 0);
      q2t[quad_feature(1,qfi)] = project_for_quad_points(colon(), 2*qfi + 1);
    }
  }

  return 0;
}

int project_quad_feature_to_tri_feature(
    const matrixst & quad_feature,
    const matrixd & quad_node,
    const matrixst & tri_feature,
    const matrixd & tri_node,
    boost::unordered_map<size_t, matrixd > & q2t)
{
  // TODO: left for shenxinxin
  matrixd node1(3, 1), node2(3, 1), node3(3, 1);
  int flag, temp_flag;           //flag=1,the nearest is foot point;
  //flag=2,the nearest is node2;
  //flag=3,the nearest is node3;
  double length, min_length;
  size_t feature_line_index;
  double dot_sign1, dot_sign2, lamda;
  std::vector<bool> is_point_visit(quad_node.size(2), false);
  matrixd foot_point;
  for(size_t i = 0; i < quad_feature.size(2); ++i){
    for(size_t j = 0; j < 2; ++j){
      if(!is_point_visit[quad_feature(j, i)]) {
        node1 = quad_node(colon(), quad_feature(j, i));
        min_length = 1000000;
        for(size_t k = 0; k < tri_feature.size(2); ++k){
          node2 = tri_node(colon(), tri_feature(0, k));
          node3 = tri_node(colon(), tri_feature(1, k));
          dot_sign1 = dot(node2 - node1, node2 - node3);
          dot_sign2 = dot(node3 - node1, node3 - node2);
          if(dot_sign1 > 0 && dot_sign2 > 0){
            temp_flag = 1;
            if(norm(node3 - node2) < 1e-8){
              //cout << node2 << endl;
              //cout << node3 << endl;
              cerr << "the length of the feature line is zero" <<endl;
              return 1;
            }
            length = norm(cross(node1 - node2, node1 - node3)) /
                     norm(node3 - node2);
          }else{
            if(dot_sign1 <= 0){
              temp_flag = 2;
              length = norm(node1 - node2);
            }else{
              temp_flag = 3;
              length = norm(node1 - node3);
            }
          }
          if(length < min_length){
            feature_line_index = k;
            flag = temp_flag;
            min_length = length;
          }
        }
        if(flag == 1) {
          node2 = tri_node(colon(), tri_feature(0, feature_line_index));
          node3 = tri_node(colon(), tri_feature(1, feature_line_index));
          if(fabs(dot(node2 - node3, node2 - node3)) < 1e-8) {
            // cout << "node2: " << node2 <<endl;
            // cout << "node3: " << node3 <<endl;
            cerr << "the length of the feature line is zero" << endl;
            return 1;
          }
          lamda = dot(node2 - node1, node2 - node3) /
                  dot(node2 - node3, node2 - node3);
          foot_point = node2 + lamda * (node3 - node2);
          q2t.insert(make_pair(quad_feature(j, i), foot_point));
        }else if(flag == 2){
          q2t.insert(make_pair(
                       quad_feature(j, i),
                       tri_node(colon(), tri_feature(0, feature_line_index))));
        }else{
          q2t.insert(make_pair(
                       quad_feature(j, i),
                       tri_node(colon(),
                                tri_feature(1, feature_line_index))));
        }
        is_point_visit[quad_feature(j,i)] = true;
      }
    }
  }

  std::cerr << "# [info] finish project quad line to triangle line."
            << std::endl;
  return 0;
}
//===============================================================================

//===============================================================================
