#include <iostream>
#include <numeric>
#include <map>
#include <fstream>
#include <memory>

#include <boost/tuple/tuple_comparison.hpp>
#include <boost/unordered_set.hpp>
#include <zjucad/matrix/matrix_expression.h>
#include <jtflib/util/container_operation.h>
#include <jtflib/mesh/util.h>
#include <mesquite/Mesquite_all_headers.hpp>

#include <hjlib/function/function.h>
#include "util.h"

#include "../common/vtk.h"
#include <zjucad/optimizer/optimizer.h>
#include <jtflib/util/container_operation.h>
#include <jtflib/optimizer/optimizer.h>
#include "../vol_param/descriptor/def.h"
#include "../vol_param/descriptor/func_terms/degree_smooth.h"
#include "../vol_param/descriptor/func_terms/arap.h"
#include "../vol_param/descriptor/func_terms/linear_equation.h"
#include "../vol_param/descriptor/func_terms/vol-anti-flip.h"
using namespace std;
using namespace zjucad::matrix;

namespace jtf{
  namespace tetmesh{

    int smooth_one_ring_point_normal(
        matrix<double> &point_normal,
        const boost::unordered_map<size_t, boost::unordered_set<size_t> > &one_ring_point_of_p,
        const size_t iter_num)
    {
      matrix<double> avg_normal = zeros<double>(3,1);
      for(size_t it = 0; it < iter_num; ++it){
          for(boost::unordered_map<size_t, boost::unordered_set<size_t> >::const_iterator cit = one_ring_point_of_p.begin();
              cit != one_ring_point_of_p.end(); ++cit){

              const boost::unordered_set<size_t> & linked_points = cit->second;
              avg_normal = zeros<double>(3,1);
              for(boost::unordered_set<size_t>::const_iterator scit = linked_points.begin();
                  scit != linked_points.end(); ++scit){
                  avg_normal += point_normal(colon(), *scit);
                }
              assert(linked_points.size() != 0);
              avg_normal /= linked_points.size();
              point_normal(colon(),cit->first) =
                  (point_normal(colon(),cit->first) + avg_normal)/2.0;
              const double len = norm(point_normal(colon(),cit->first));
              if(len > 1e-6) point_normal(colon(),cit->first) /= len;
            }
        }
      return 0;
    }

    void one_ring_tet_at_point::add_tets(const matrixst & tets)
    {
      const size_t max_idx = *max_element(tets.begin(), tets.end());
      p2t_.resize(max_idx + 1);
      for(size_t ti = 0; ti < tets.size(2); ++ti){
          for(size_t pi = 0; pi < tets.size(1); ++pi){
              p2t_[tets(pi,ti)].push_back(ti);
            }
        }
    }

    void one_ring_edge_at_point::add_tets(const matrixst &tets)
    {
      const size_t max_idx = *max_element(tets.begin(), tets.end());
      p2p_.resize(max_idx + 1);
      for(size_t ti = 0; ti < tets.size(2); ++ti){
          for(size_t pi = 0; pi < tets.size(1); ++pi){
              for(size_t j = 1; j < 4; ++j)
                p2p_[tets(pi,ti)].insert(tets((pi + j)%tets.size(1), ti));
            }
        }
    }

    void one_ring_point_at_point::build(
        const std::vector<std::pair<size_t,size_t> > & edges)
    {
      for(size_t ei = 0; ei < edges.size(); ++ei){
          const pair<size_t,size_t> & one_edge = edges[ei];
          p2p_[one_edge.first].insert(one_edge.second);
          p2p_[one_edge.second].insert(one_edge.first);
        }
    }

    void one_ring_point_at_point::build(const matrix<size_t> &tet)
    {
      for(size_t ti = 0; ti < tet.size(2); ++ti){
          for(size_t p = 0; p < tet.size(1); ++p){
              for(size_t np = p + 1; np < tet.size(1); ++np){
                  p2p_[tet(p, ti)].insert(tet(np, ti));
                  p2p_[tet(np, ti)].insert(tet(p, ti));
                }
            }
        }
    }

    //    int cal_one_face_normal(const matrixst &outside_face,
    //                            const matrixd &node,
    //                            matrixd &face_normal,
    //                            double *face_normal_weight)
    //    {
    //      assert(outside_face.size() == 3);
    //      face_normal.resize(3,1);

    //      matrixd normal = zeros<double>(3,1);
    //      matrixd edge0 = zeros<double>(3,1), edge1 = zeros<double>(3,1);

    //      edge0 = node(colon(),outside_face[1]) - node(colon(),outside_face[0]);
    //      edge1 = node(colon(),outside_face[2]) - node(colon(),outside_face[1]);

    //      normal = cross(edge0,edge1);
    //      double length = norm(normal);
    //      if(length < 1e-8){
    //          cerr << "# [error] this triangle is degenerated." << endl;
    //          length = 1;
    //        }
    //      if(face_normal_weight)
    //        *face_normal_weight = length;
    //      normal /= length;
    //      face_normal = normal;
    //      return 0;
    //    }

    //    int cal_face_normal(const matrixst &outside_face,
    //                        const matrixd &node,
    //                        matrixd &face_normal,
    //                        matrixd *face_normal_weight)
    //    {
    //      assert(outside_face.size(1) == 3);
    //      face_normal.resize(3,outside_face.size(2));

    //      if(face_normal_weight)
    //        face_normal_weight->resize(outside_face.size(2));

    //      matrixd normal = zeros<double>(3,1);
    //      matrixd edge0 = zeros<double>(3,1), edge1 = zeros<double>(3,1);

    //      for(size_t t = 0; t < outside_face.size(2); ++t){
    //          edge0 = node(colon(),outside_face(1,t)) - node(colon(),outside_face(0,t));
    //          edge1 = node(colon(),outside_face(2,t)) - node(colon(),outside_face(1,t));

    //          normal = cross(edge0,edge1);
    //          double length = norm(normal);
    //          if(length < 1e-8){
    //              cerr << "# [error] this triangle is degenerated." << endl;
    //              length = 1;
    //            }
    //          if(face_normal_weight)
    //            (*face_normal_weight)[t] = length;
    //          normal /= length;
    //          face_normal(colon(),t) = normal;
    //        }
    //      return 0;
    //    }

    int orient_one_face_normal_outside_tetmesh(
        const matrixst &tet,
        const matrixd &node,
        const matrixst &outside_face,
        const size_t &outside_face_idx,
        const jtf::mesh::face2tet_adjacent &fa,
        matrixd &face_normal)
    {
      assert(outside_face.size() == 3);
      assert(face_normal.size() == 3);
      matrixd outpoint_vec = zeros<double>(3,1);
      const pair<size_t,size_t> &two_tets = fa.face2tet_[outside_face_idx];
      assert(fa.is_outside_face(two_tets));
      const size_t &tet_idx = (two_tets.first == -1)?two_tets.second:two_tets.first;
      //TODO: this step may introduce overflow
      const size_t other_point_idx =
          accumulate(tet(colon(),tet_idx).begin(),tet(colon(),tet_idx).end(),
                     static_cast<size_t>(0)) -
          accumulate(outside_face.begin(), outside_face.end(),static_cast<size_t>(0));

      assert(find(tet(colon(),tet_idx).begin(), tet(colon(),tet_idx).end(),
                  other_point_idx) != tet(colon(),tet_idx).end());
      assert(find(outside_face.begin(),outside_face.end(),
                  other_point_idx) == outside_face.end());

      const size_t out_point = outside_face[0];
      outpoint_vec  = node(colon(),out_point) - node(colon(), other_point_idx);
      if(dot(outpoint_vec,face_normal) < 0){// this face_normal needs to be flipped
          face_normal *= -1;
        }
      return 0;
    }

    int orient_face_normal_outside_tetmesh(const matrixst &tet,
                                           const matrixd &node,
                                           const matrixst &outside_face,
                                           const matrixst &outside_face_idx,
                                           const jtf::mesh::face2tet_adjacent &fa,
                                           matrixd &face_normal)
    {
      assert(outside_face.size(2) == outside_face_idx.size());
      assert(face_normal.size(2) == outside_face.size(2));

      matrixd outpoint_vec = zeros<double>(3,1);
      for(size_t t = 0; t < outside_face.size(2); ++t){
          const pair<size_t,size_t> &two_tets = fa.face2tet_[outside_face_idx[t]];
          assert(fa.is_outside_face(two_tets));
          const size_t &tet_idx = (two_tets.first == -1)?two_tets.second:two_tets.first;
          // TODO: this step may introduce overflow
          const size_t other_point_idx = accumulate(tet(colon(),tet_idx).begin(),tet(colon(),tet_idx).end(),0)
              - accumulate(outside_face(colon(),t).begin(),outside_face(colon(),t).end(),0);

          assert(find(tet(colon(),tet_idx).begin(), tet(colon(),tet_idx).end(),other_point_idx) != tet(colon(),tet_idx).end());
          assert(find(outside_face(colon(),t).begin(),outside_face(colon(),t).end(),other_point_idx) == outside_face(colon(),t).end());

          const size_t out_point = outside_face(0,t);//(tet(0,tet_idx) == other_point_idx)?tet(1,tet_idx):tet(0,tet_idx);
          outpoint_vec  = node(colon(),out_point) - node(colon(), other_point_idx);
          if(dot(outpoint_vec,face_normal(colon(),t)) < 0)// this face_normal needs to be flipped
            {
              face_normal(colon(),t) *= -1;
              //swap(outside_face(0,t),outside_face(1,t));
            }

        }

      return 0;
    }


    int trans_face_normal_to_point_normal(const matrixd &face_normal,
                                          const matrixst &outside_face,
                                          matrixd & point_normal,
                                          matrixst & point_idx,
                                          map<size_t,size_t> &point_map,
                                          const matrixd *face_normal_weight)
    {
      assert(face_normal.size(2) == outside_face.size(2));
      if(face_normal_weight)
        assert(face_normal.size(2) == face_normal_weight->size());
      map<size_t,vector<size_t> > point_adj_face;

      for(size_t t = 0; t < outside_face.size(2); ++t){
          for(size_t i = 0; i < outside_face.size(1); ++i)
            point_adj_face[outside_face(i,t)].push_back(t);
        }

      point_normal.resize(3, point_adj_face.size());
      typedef map<size_t,vector<size_t> >::const_iterator mcit;
      matrixd normal = zeros<double>(3,1);
      //cerr << "# [info] outside point size = " << point_adj_face.size() << endl;
      point_idx.resize(point_adj_face.size());

      size_t point_idx_ = 0;
      for(mcit it = point_adj_face.begin(); it != point_adj_face.end(); ++it){
          normal = zeros<double>(3,1);
          vector<double> each_point_weight;
          for(size_t t = 0; t < it->second.size(); ++t){
              if(face_normal_weight)
                each_point_weight.push_back((*face_normal_weight)[it->second[t]]);
              else
                each_point_weight.push_back(1.0);
              normal += face_normal(colon(),it->second[t]) * each_point_weight.back();
            }
          double total_weight = accumulate(each_point_weight.begin(),each_point_weight.end(),0.0);

          if(norm(normal) < 1e-8 || fabs(total_weight) < 1e-8){
              cerr << "# [error] this normal is degenerated." << endl;
              return __LINE__;
            }
          normal /= total_weight;
          normal /= norm(normal);
          //cerr << "# [info] " << normal << endl;
          point_normal(colon(),point_idx_) = normal;
          point_idx[point_idx_] = it->first;
          point_map[it->first] = point_idx_;
          ++point_idx_;
        }

      return 0;
    }


    int laplace_smoothing_triangle_mesh(const matrixst &tri_mesh,
                                        matrixd &node)
    {
      map<size_t,set<size_t> > point_adj_point;
      for(size_t t = 0; t < tri_mesh.size(2); ++t){
          for(size_t i = 0; i < tri_mesh.size(1); ++i){
              point_adj_point[tri_mesh(i,t)].insert(tri_mesh((i+1)%3,t));
              point_adj_point[tri_mesh((i+1)%3,t)].insert(tri_mesh(i,t));
            }
        }

      const size_t point_num = point_adj_point.size();
      matrixd new_surface_node = zeros<double>(3,point_num);

      size_t t = 0;
      typedef map<size_t,set<size_t> >::const_iterator mcit_;
      typedef set<size_t>::const_iterator scit;
      for(mcit_ mcit = point_adj_point.begin(); mcit != point_adj_point.end(); ++mcit){
          // for(size_t i = 0; i < mcit->second.size(); ++i)
          for(scit it = mcit->second.begin(); it != mcit->second.end(); ++it)
            new_surface_node(colon(),t) += node(colon(),*it);
          new_surface_node(colon(),t) /= mcit->second.size();
          ++t;
        }

      t = 0;
      for(mcit_ mcit = point_adj_point.begin(); mcit != point_adj_point.end(); ++mcit){
          node(colon(),mcit->first) = new_surface_node(colon(),t);
          ++t;
        }
      return 0;
    }

    int subdivide_prism_into_tets(const matrixst &prism,
                                  matrixst &tets)
    {
#ifdef DEBUG
      set<size_t> validate_prism(prism.begin(),prism.end());
      if(validate_prism.size() != 6) return __LINE__;
#endif
      assert(prism.size(1) == 4 && prism.size(2) == 3);
      tets.resize(4,3);

      //edges: <p0,p1,face_idx>
      vector<tuple<size_t,size_t,size_t> > three_edges;
      three_edges.reserve(3);
      // each quad face of prism is:
      // 0----2
      // |    |
      // 1----3
      // the diagonal line is edge<idx,3-idx>
      typedef matrixst::const_iterator mcit;
      for(size_t t = 0; t < 3; ++t){
          mcit it_min = min_element(prism(colon(),t).begin(),prism(colon(),t).end());
          size_t min_idx = it_min - prism(colon(),t).begin();
          three_edges.push_back(make_tuple(*it_min,prism(3-min_idx,t),t));
        }
      sort(three_edges.begin(),three_edges.end());

      assert(get<0>(three_edges[0]) == get<0>(three_edges[1]));

      // construct tet0
      tets(0,0) = get<0>(three_edges[0]);
      tets(1,0) = get<1>(three_edges[0]);
      tets(2,0) = get<1>(three_edges[1]);
      tets(3,0) = get<0>(three_edges[2]);

      // construct tet1
      size_t face_idx_2 = 0;
      for(size_t t = 0; t < 3; ++t){ // find the face which contain point: three_edges[2].get<0>()
          if(t != get<2>(three_edges[2])){
              if(find(prism(colon(),t).begin(),prism(colon(),t).end(),get<0>(three_edges[2])) != prism(colon(),t).end())
                face_idx_2 = t;
            }
        }
      //  size_t a = accumulate(prism(colon(),face_idx_2).begin(),prism(colon(),face_idx_2).end(),0);
      //  size_t b = three_edges[2].get<0>() + three_edges[0].get<0>() + ((three_edges[0].get<2>() == face_idx_2)?three_edges[0].get<1>():three_edges[1].get<1>());
      //  cerr << a - b << endl;
      size_t other_point = accumulate(prism(colon(),face_idx_2).begin(),prism(colon(),face_idx_2).end(),0)
          - (get<0>(three_edges[2]) + get<0>(three_edges[0]) + ((get<2>(three_edges[0]) == face_idx_2)?get<1>(three_edges[0]):get<1>(three_edges[1])));

      tets(0,1) = get<0>(three_edges[0]);
      tets(1,1) = get<1>(three_edges[0]);
      tets(2,1) = get<1>(three_edges[1]);
      tets(3,1) = other_point;

      // construct tet2
      size_t face_idx_3 = get<2>(three_edges[2]);
      //  cerr << accumulate(prism(colon(),face_idx_3).begin(),prism(colon(),face_idx_3).end(),0) << " "
      //       << three_edges[0].get<1>() << " " << three_edges[1].get<1>() << " " <<  three_edges[2].get<0>();
      size_t last_point = accumulate(prism(colon(),face_idx_3).begin(),prism(colon(),face_idx_3).end(),0)
          - get<1>(three_edges[0]) - get<1>(three_edges[1]) - get<0>(three_edges[2]);

      tets(0,2) = get<0>(three_edges[2]);
      tets(1,2) = get<1>(three_edges[2]);
      tets(2,2) = get<0>(three_edges[0]);
      tets(3,2) = last_point;
      return 0;
    }


    int get_shared_face(const std::pair<size_t,size_t> &tet_pair,
                        const matrixst &tet,
                        size_t *face)
    {
      if(tet_pair.first == -1 || tet_pair.second == -1) return __LINE__;
      assert(face);

      matrixst tet0 = tet(colon(),tet_pair.first);
      matrixst tet1 = tet(colon(),tet_pair.second);
      vector<size_t> face_v(tet0.size());
      vector<size_t>::iterator it =
          find_intersection_set(tet0.begin(),tet0.end(),tet1.begin(),tet1.end(),face_v.begin());
      if(static_cast<size_t>(it - face_v.begin()) != 3) return __LINE__;
      //face.resize(3,1);
      //copy(face_v.begin(),it,face.begin());
      copy(face_v.begin(),it,face);
      return 0;
    }

    int extend_tetmesh(const matrix<size_t> & tri_face,
                       const matrix<double> & node,
                       const matrix<double> & point_normal,
                       jtf::mesh::meshes & tet_layer)
    {
      assert(point_normal.size(2) == node.size(2));
      const double average_len = jtf::mesh::cal_average_edge(tri_face,node)/5;
      set<size_t> all_used_node(tri_face.begin(), tri_face.end());

      map<size_t,size_t> surface_point_map;
      matrixd new_point_layer(3,all_used_node.size());
      size_t t = 0;
      for(const auto & one_p : all_used_node){
          new_point_layer(colon(),t) = node(colon(),one_p) + point_normal(colon(),one_p) * average_len;
          surface_point_map[one_p] = t + node.size(2);
          ++t;
        }

      tet_layer.node_.resize(3, node.size(2) + new_point_layer.size(2));
      copy(node.begin(),node.end(),tet_layer.node_.begin());
      copy(new_point_layer.begin(),new_point_layer.end(),tet_layer.node_.begin() + node.size());

      // construct new tets
      //       a
      //      /|\
      //     b-|-c
      //     | 0 |
      //     |/ \|
      //     1---2

      matrixst new_tet_layer(4,tri_face.size(2) * 3);
      // each outside face introduce a prism which contains 3 tets
      matrixst prism(4,3);
      matrixst tets(4,3);
      for(size_t t = 0; t < tri_face.size(2); ++t){
          { //contruct the prism
            prism(0,0) = surface_point_map[tri_face(0,t)];
            prism(1,0) = tri_face(0,t);
            prism(2,0) = surface_point_map[tri_face(1,t)];
            prism(3,0) = tri_face(1,t);

            prism(0,1) = surface_point_map[tri_face(1,t)];
            prism(1,1) = tri_face(1,t);
            prism(2,1) = surface_point_map[tri_face(2,t)];
            prism(3,1) = tri_face(2,t);

            prism(0,2) = surface_point_map[tri_face(2,t)];
            prism(1,2) = tri_face(2,t);
            prism(2,2) = surface_point_map[tri_face(0,t)];
            prism(3,2) = tri_face(0,t);

            if(jtf::tetmesh::subdivide_prism_into_tets(prism,tets)){
                cerr << "# [error] can not subdivide this prism into tets." << endl;
                copy(prism.begin(),prism.end(),ostream_iterator<size_t>(cerr," "));
                return __LINE__;
              }
            copy(tets.begin(),tets.end(),new_tet_layer.begin() + t * 4 * 3);
          }
        }

      tet_layer.mesh_.resize(4, new_tet_layer.size(2));
      copy(new_tet_layer.begin(),new_tet_layer.end(),tet_layer.mesh_.begin());

      orient_tet(tet_layer.node_,tet_layer.mesh_);
    }

    int extend_tetmesh(const matrixst & tet,
                       const matrixd & node,
                       matrixst &new_tet,
                       matrixd &new_node,
                       map<size_t,size_t> & surface_point_map,
                       matrix<double> * zyz)
    {
      unique_ptr<jtf::mesh::face2tet_adjacent> fa(jtf::mesh::face2tet_adjacent::create(tet));
      if(!fa.get())  return __LINE__;
      matrixst outside_face;
      matrixst outside_face_idx;
      jtf::mesh::get_outside_face(*fa,outside_face,true);
      get_outside_face_idx(*fa,outside_face_idx);

      matrixd face_normal,point_normal;//(3,outside_face.size(2));
      matrixst point_idx;
      // map<size_t,size_t> point_map;
      matrixd face_normal_weight = ones<double>(outside_face_idx.size(),1);
      jtf::mesh::cal_face_normal(outside_face,node,face_normal);
      jtf::tetmesh::orient_face_normal_outside_tetmesh(
            tet,node,outside_face,outside_face_idx,*fa,face_normal);
      jtf::tetmesh::trans_face_normal_to_point_normal(
            face_normal,outside_face,point_normal,point_idx,surface_point_map,&face_normal_weight);

      //const double average_len = jtf::mesh::cal_average_edge(tet,node)/5;
      const double average_len = jtf::mesh::cal_average_edge(tet,node)/5;
      assert(average_len > 0);

      // construct new points
      matrixd new_point_layer(3,point_idx.size());
      for(size_t t = 0; t < new_point_layer.size(2); ++t){
          new_point_layer(colon(),t) = node(colon(),point_idx[t])
              + point_normal(colon(),t) * average_len;
        }

      new_node.resize(3, node.size(2) + point_idx.size());
      copy(node.begin(),node.end(),new_node.begin());
      copy(new_point_layer.begin(),new_point_layer.end(),new_node.begin() + node.size());

      for(map<size_t,size_t>::iterator mit = surface_point_map.begin();
          mit != surface_point_map.end(); ++mit) mit->second += node.size(2);

      // construct new tets
      //       a
      //      /|\
      //     b-|-c
      //     | 0 |
      //     |/ \|
      //     1---2

      matrixst new_tet_layer(4,outside_face.size(2) * 3);
      // each outside face introduce a prism which contains 3 tets
      matrixst prism(4,3);
      matrixst tets(4,3);
      for(size_t t = 0; t < outside_face.size(2); ++t){
          { //contruct the prism
            prism(0,0) = surface_point_map[outside_face(0,t)];
            prism(1,0) = outside_face(0,t);
            prism(2,0) = surface_point_map[outside_face(1,t)];
            prism(3,0) = outside_face(1,t);

            prism(0,1) = surface_point_map[outside_face(1,t)];
            prism(1,1) = outside_face(1,t);
            prism(2,1) = surface_point_map[outside_face(2,t)];
            prism(3,1) = outside_face(2,t);

            prism(0,2) = surface_point_map[outside_face(2,t)];
            prism(1,2) = outside_face(2,t);
            prism(2,2) = surface_point_map[outside_face(0,t)];
            prism(3,2) = outside_face(0,t);

            if(jtf::tetmesh::subdivide_prism_into_tets(prism,tets)){
                cerr << "# [error] can not subdivide this prism into tets." << endl;
                copy(prism.begin(),prism.end(),ostream_iterator<size_t>(cerr," "));
                return __LINE__;
              }
            copy(tets.begin(),tets.end(),new_tet_layer.begin() + t * 4 * 3);
          }
        }

      new_tet.resize(4,tet.size(2) + new_tet_layer.size(2));
      copy(tet.begin(),tet.end(),new_tet.begin());
      copy(new_tet_layer.begin(),new_tet_layer.end(),new_tet.begin() + tet.size());

      orient_tet(new_node,new_tet);

      if(zyz){
          matrix<double> new_zyz(3, new_tet.size(2));
          std::copy(zyz->begin(), zyz->end(), new_zyz.begin());
          for(size_t i = 0; i < outside_face_idx.size(); ++i){
              const size_t & face_idx = outside_face_idx[i];
              assert(face_idx != -1);
              const pair<size_t,size_t> & tet_pair = fa->face2tet_[face_idx];
              const size_t other_tet = tet_pair.first==-1?tet_pair.second:tet_pair.first;

              for(size_t j = 0; j < 3; ++j){
                  new_zyz(colon(), 3 * i + tet.size(2) + j) = (*zyz)(colon(),other_tet);
                }
            }
          *zyz = new_zyz;
        }

#if 1
      for(size_t ti = 0; ti < new_tet.size(2); ++ti) {
          matrixd ele(3, 3);
          for(size_t ni = 0; ni < 3; ++ni)
            ele(zjucad::matrix::colon(), ni) = new_node(zjucad::matrix::colon(), new_tet(ni+1, ti))
                - new_node(zjucad::matrix::colon(), new_tet(0, ti));
          if(zjucad::matrix::dot(zjucad::matrix::cross(ele(zjucad::matrix::colon(), 0),
                                                       ele(zjucad::matrix::colon(), 1)),
                                 ele(zjucad::matrix::colon(), 2)) < 0) {
              cerr << "# [error] tet idx " << ti << " volume is not correct: " <<
                      new_tet(0,ti) << " " << new_tet(1,ti) << " " << new_tet(2,ti)
                   << " " << new_tet(3,ti) << endl;
            }
        }
#endif

      return 0;
    }

    int calculate_face_normal(const matrixst &cut_tet,
                              const matrixd &cut_tet_node,
                              const size_t &tet_idx,
                              const matrixst &one_face,
                              matrixd &one_normal)
    {
      const matrixd edge[2] = {
        cut_tet_node(colon(), one_face[1]) - cut_tet_node(colon(), one_face[0]),
        cut_tet_node(colon(), one_face[2]) - cut_tet_node(colon(), one_face[0]),
      };
      one_normal = cross(edge[0], edge[1]);
      const double len = norm(one_normal);
      if(len > 1e-8)
        {
          one_normal /= len;
          {// test whether the nomral is right or oppoiste
            matrixd tet_bary_center = zeros<double>(3,1);
            matrixd face_bary_center = zeros<double>(3,1);
            for(size_t t = 0; t < 4; ++t)
              tet_bary_center += cut_tet_node(colon(),cut_tet(t,tet_idx));
            tet_bary_center /= 4;
            for(size_t fi = 0; fi < 3; ++fi) {
                face_bary_center += cut_tet_node(colon(), one_face[fi]);
              }
            face_bary_center /= 3;
            const matrixd face_to_tet_bary = face_bary_center - tet_bary_center;

            if(dot(face_to_tet_bary,one_normal) < 0.0)
              one_normal = -one_normal;
          }
          return 0;
        }
      else
        return 1;
    }

    int calculate_face_normal_wrapper(const matrixst &tet,
                                      const matrixd &tet_node,
                                      const size_t &tet_idx,
                                      const size_t *face,
                                      matrixd &normal)
    {
      itr_matrix<const size_t*> face0(3,1,&face[0]);
      return calculate_face_normal(tet,tet_node, tet_idx, face0,normal);
    }

    int tetmesh_quality_improver::init()
    {
      fa_.reset(jtf::mesh::face2tet_adjacent::create(tet_));
      if(!fa_.get()){
          cerr << "# [error] can not buildjtf::mesh::face2tet_adjacent." << endl;
          return __LINE__;
        }

      get_outside_face(*fa_, outside_face_);
      get_outside_face_idx(*fa_, outside_face_idx_);
      return 0;
    }

    int tetmesh_quality_improver::improve_suface_by_one_ring()
    {
      matrixd face_normal(3, outside_face_.size(2));
      jtf::mesh::cal_face_normal(outside_face_, node_, face_normal);
      orient_face_normal_outside_tetmesh(tet_, node_, outside_face_, outside_face_idx_,
                                         *fa_, face_normal);

      boost::unordered_map<size_t, boost::unordered_set<size_t> > p2f, p2adj_p;
      for(size_t fi = 0; fi < outside_face_.size(2); ++fi){
          for(size_t pi = 0; pi < outside_face_.size(1); ++pi){
              p2f[outside_face_(pi,fi)].insert(fi);
              p2adj_p[outside_face_(pi,fi)].insert(
                    outside_face_((pi+1)%outside_face_.size(1),fi));
              p2adj_p[outside_face_(pi,fi)].insert(
                    outside_face_((pi+2)%outside_face_.size(1),fi));
            }
        }
      assert(p2f.size() == p2adj_p.size());

      matrixd normal = zeros<double>(3,1), average_node = zeros<double>(3,1);
      for(boost::unordered_map<size_t, boost::unordered_set<size_t> >::const_iterator
          cit = p2f.begin(); cit != p2f.end(); ++cit){
          normal = zeros<double>(3,1);
          average_node = zeros<double>(3,1);
          for(boost::unordered_set<size_t>::const_iterator ccit = cit->second.begin();
              ccit != cit->second.end(); ++ccit){
              normal += face_normal(colon(), *ccit);
            }
          const double len = norm(normal);
          if(len > 1e-6)
            normal /= len;
          const boost::unordered_set<size_t> & adj_points = p2adj_p[cit->first];
          if(adj_points.empty()) continue;
          for(boost::unordered_set<size_t>::const_iterator ccit = adj_points.begin();
              ccit != adj_points.end(); ++ccit){
              average_node += node_(colon(), *ccit);
            }
          average_node /= adj_points.size();

          // to calculate the project point on normal plane
          // for a point p with normal n, and another point p'
          // the project point from p' to the normal plane will be P = p' + d*n
          // then (P-p)n=0 => d = (p-p')n

          const double d = dot(node_(colon(),cit->first) - average_node, normal);
          node_(colon(), cit->first) = average_node + d * normal;
        }
      return 0;
    }


    int tetmesh_quality_improver::improve(const string surface_strategy)
    {
      //      using namespace Mesquite;

      //      bool fix_surface = true;
      //      vector<int> fix_node(node_.size(2), 0);

      //      if(surface_strategy == "one_ring"){
      //          improve_suface_by_one_ring();
      //          fix_surface = true;
      //        }

      //      if(fix_surface){
      //          for(size_t pi = 0; pi < outside_face_.size(); ++pi)
      //            fix_node[outside_face_[pi]] = 1;
      //        }

      //      MsqError err;

      //      ArrayMesh mesh(3, node_.size(2), &node_[0], &fix_node[0], tet_.size(2),
      //          TETRAHEDRON, &tet_[0]);

      //      Mesquite::ShapeImprover optimizer;
      //      //SizeAdaptShapeWrapper optimizer(1e-3);
      //      //LaplaceWrapper optimizer;
      //      optimizer.run_instructions(&mesh, err);

      //      if (err) {
      //          cerr << "# [error] jtf::mesh::meshes quality improve fail " << err << endl;
      //          return __LINE__;
      //        }
      //      return 0;
    }

    void deform_tet_accordint_to_surface_constraints(
        const zjucad::matrix::matrix<size_t> & def_mesh,
        zjucad::matrix::matrix<double> & def_node,
        const zjucad::matrix::matrix<size_t> & orig_mesh,
        const zjucad::matrix::matrix<double> & orig_node,
        const std::map<size_t,size_t> & old_surface_point_to_new_point_map,
        boost::property_tree::ptree & pt)
    {
      double w = 0.05;

      for(size_t j = 0; j < 5; ++j){
          vector<jtf_func_ptr> all_funcs;
          {
            SIMPLEX sim = TET;
            hj_func_ptr fptr(build_arap_func(def_mesh, def_node,sim));
            all_funcs.push_back(jtf_func_ptr(jtf::function::least_square_warpper(fptr)));
          }

          {
            w *= 1.5;
            // point constraints
            for(const auto & one_p2p : old_surface_point_to_new_point_map){
                const size_t point_idx = one_p2p.second;
                all_funcs.push_back(
                      jtf_func_ptr(
                        jtf::function::least_square_warpper(
                          shared_ptr<hj::function::function_t<double,int32_t> >(
                            new point_fix_hj(def_node.size(2), point_idx,
                                             orig_node(colon(), one_p2p.first),
                                             w)))));
              }
          }

          jtf_func_ptr sum_func_(new jtf::function::sum_function<double,int32_t,jtf::function::SMART_STD>(all_funcs));
          jtf::optimize(*sum_func_, def_node,pt, nullptr, nullptr, nullptr);
        }
    }
  }}
