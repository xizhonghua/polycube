#include <vector>
#include <numeric>
#include <stack>
#include <jtflib/mesh/util.h>
#include <jtflib/math/math.h>
#include <jtflib/function/function.h>
#include <jtflib/function/func_aux.h>
#include <hjlib/function/operation.h>

#include "../tetmesh/tetmesh.h"
#include "../tetmesh/util.h"
#include "frame_function.h"
#include "function_term.h"
#include "hex_function_term.h"
#include "../numeric/util.h"
#include "angle_defect.h"

using namespace std;
using namespace hj::function;
using namespace zjucad::matrix;

std::shared_ptr<sym_frame_opt_sys> sfptr;

hj::function::function_t<double, int32_t> *
build_frame_smooth_func(const matrixst & tet,
                        const matrixd & node,
                        std::shared_ptr<jtf::mesh::face2tet_adjacent> &fa,
                        const string &func_type,
                        const string smooth_strategy)
{
  if(!fa.get()){
      fa.reset(jtf::mesh::face2tet_adjacent::create(tet));
    }

  if(!sfptr.get())
    sfptr.reset(new sym_frame_opt_sys(tet.size(2)));

  std::shared_ptr<vector<std::shared_ptr<function_t<double,int32_t> > > > funcs(
        new vector<std::shared_ptr<function_t<double, int32_t> > >);

  std::vector<double> volume_of_tets(tet.size(2));
  for(size_t ti = 0; ti < volume_of_tets.size(); ++ti)
    volume_of_tets[ti] = fabs(jtf::mesh::cal_tet_vol(node(colon(), tet(colon(), ti))));

  const double total_v = std::accumulate(volume_of_tets.begin(), volume_of_tets.end(), 0.0);

  if(smooth_strategy == "face"){
      for(size_t fi = 0; fi < fa->face2tet_.size(); ++fi) {
          const pair<size_t,size_t> & tet_pair = fa->face2tet_[fi];
          if(fa->is_outside_face(tet_pair)) continue;

          matrixd bary_i = zeros<double>(3,1);
          matrixd bary_j = zeros<double>(3,1);
          for(size_t v = 0; v < 4; ++v) {
              bary_i += node(colon(), tet(v,tet_pair.first));
              bary_j += node(colon(), tet(v,tet_pair.second));
            }

          bary_i /= 4;
          bary_j /= 4;

          const double dij = norm(bary_i - bary_j);

          const double weight = sqrt((volume_of_tets[tet_pair.first] + volume_of_tets[tet_pair.second] )/4)
              * pow(total_v,-1.0/6.0) /dij  ;

          //            const double weight =
          //                    sqrt((volume_of_tets[tet_pair.first] + volume_of_tets[tet_pair.second] )/(4.0 * total_v));

          if(func_type == "sh2zyz" )
            funcs->push_back(
                  std::shared_ptr<inner_smooth_sh>(
                    new inner_smooth_sh(sfptr,tet_pair.first, tet_pair.second, weight)));
          else if(func_type == "sh")
            funcs->push_back(
                  std::shared_ptr<inner_smooth>(
                    new inner_smooth(tet.size(2), tet_pair.first, tet_pair.second, weight)));
          else {
              cerr << "# [error] unsupported func type: " << func_type << endl;
              return 0;
            }

        }
    }else if(smooth_strategy == "edge"){
      jtf::mesh::one_ring_tet_at_edge ortae;
      ortae.add_tets(tet, *fa);
      ortae.sort_into_loop(tet, node);

      for(jtf::mesh::one_ring_tet_at_edge::e2tet_type::const_iterator eit =
          ortae.e2t_.begin(); eit != ortae.e2t_.end(); ++eit){
          const vector<size_t> & loop = eit->second;
          size_t l_begin, l_end;
          if(ortae.is_inner_edge(loop)){
              l_begin = 0;
              l_end = loop.size();
            }else{
              l_begin = 1;
              l_end = loop.size() - 1;
            }

          //            for(size_t t = l_begin; t < l_end-1; ++t){
          //                const double weight =
          //                        sqrt((volume_of_tets[loop[t]] + volume_of_tets[loop[t+1]])/(6.0 * total_v));
          for(size_t t = l_begin; t < l_end-1; ++t){
              for(size_t tj = t + 1; tj < l_end-1; ++tj){
                  const double weight =
                      sqrt((volume_of_tets[loop[t]] + volume_of_tets[loop[tj]])
                      /(6.0 * (l_end - t - 2) * total_v));

                  if(func_type == "sh2zyz")
                    funcs->push_back(
                          std::shared_ptr<inner_smooth_sh>(
                            new inner_smooth_sh(sfptr,loop[t], loop[t+1], weight)));
                  else if(func_type == "sh"){
                      funcs->push_back(std::shared_ptr<inner_smooth>(
                                         new inner_smooth(tet.size(2), loop[t], loop[t+1], weight)));
                    }else {
                      cerr << "# [error] unsupported func type." << endl;
                      return 0;
                    }
                }
            }
        }
    }else{
      cerr << "# [error] unsupported smooth_strategy: " << smooth_strategy << endl;
      return 0;
    }

  return new_catenated_function<double, int32_t>(funcs);
}


hj::function::function_t<double, int32_t> *
build_normal_align_func(const matrixst &tet,
                        const matrixd &node,
                        const matrix<size_t> & aligned_face_idx,
                        const matrix<double> & aligned_normal,
                        std::shared_ptr<jtf::mesh::face2tet_adjacent> &fa,
                        const string & func_type,
                        const double weight)
{
  if(!fa.get()){
      fa.reset(jtf::mesh::face2tet_adjacent::create(tet));
    }

  if(!sfptr.get())
    sfptr.reset(new sym_frame_opt_sys(tet.size(2)));

  std::shared_ptr<vector<std::shared_ptr<function_t<double,int32_t> > > > funcs(
        new vector<std::shared_ptr<function_t<double, int32_t> > >);

  vector<double> area_of_surface_tri(aligned_face_idx.size());
  for(size_t fi = 0; fi < area_of_surface_tri.size(); ++fi)
    area_of_surface_tri[fi] = fabs(
          jtf::mesh::cal_face_area(&fa->faces_[aligned_face_idx[fi]][0],
        fa->faces_[aligned_face_idx[fi]].size(), node));
  const double total_area = accumulate(area_of_surface_tri.begin(), area_of_surface_tri.end(), 0.0);


  for(size_t fi = 0; fi < aligned_face_idx.size(); ++fi) {

      const pair<size_t,size_t> & tet_pair = fa->face2tet_[aligned_face_idx[fi]];
      if(!fa->is_outside_face(tet_pair)) {
          cerr << "# [error] face " << aligned_face_idx[fi] << "should be outside face." << endl;
          return 0;
        }
      const size_t & tet_idx = (tet_pair.first == -1? tet_pair.second:tet_pair.first);

      if(func_type == "sh2zyz")
        funcs->push_back(
              std::shared_ptr<align_function>(
                new align_function(sfptr, tet_idx, &aligned_normal(0,fi),
                                   sqrt(area_of_surface_tri[fi]/total_area * weight))));
      else if(func_type == "sh"){
          funcs->push_back(
                std::shared_ptr<align_sh_inner>(
                  new align_sh_inner(tet.size(2), tet_idx, &aligned_normal(0,fi),
                                     sqrt(area_of_surface_tri[fi]/total_area*weight))));
        }else {
          cerr << "# [error] unsupported func_type." << endl;
          return 0;
        }
    }

  return new_catenated_function<double, int32_t>(funcs);
}

hj::function::function_t<double, int32_t> *
build_normal_align_func_tri(
    const matrix<size_t> &tri,
    const matrix<double> &node,
    const matrix<double> & normal,
    const std::string strategy,
    const double w)
{
  std::shared_ptr<vector<std::shared_ptr<function_t<double,int32_t> > > > funcs(
        new vector<std::shared_ptr<function_t<double, int32_t> > >);

  vector<double> area_of_surface_tri(tri.size(2));
  for(size_t fi = 0; fi < area_of_surface_tri.size(); ++fi)
    area_of_surface_tri[fi] = fabs(
          jtf::mesh::cal_face_area(&tri(0,fi), tri.size(1), node));
  const double total_area = accumulate(area_of_surface_tri.begin(), area_of_surface_tri.end(), 0.0);


  for(size_t fi = 0; fi < tri.size(2); ++fi) {
      if(strategy == "sh2zyz")
        funcs->push_back(
              std::shared_ptr<align_function>(
                new align_function(sfptr, fi, &normal(0,fi),
                                   sqrt(area_of_surface_tri[fi]/total_area * w))));
      else if(strategy == "sh"){
          funcs->push_back(
                std::shared_ptr<align_sh_inner>(
                  new align_sh_inner(tri.size(2), fi, &normal(0,fi),
                                     sqrt(area_of_surface_tri[fi]/total_area*w))));
        }else {
          cerr << "# [error] unsupported func_type." << endl;
          return 0;
        }
    }

  return new_catenated_function<double, int32_t>(funcs);
}

//void
//build_normal_align_func_tri_no_quadratic(
//    const matrix<size_t> &tri,
//    const matrix<double> &node,
//    const matrix<double> & normal,
//    const std::string strategy,
//    shared_ptr<vector<shared_ptr<jtf::function::functionN1_t<double,int32_t> > > > constraints)
//{
//  if(strategy != "sh"){
//      throw std::invalid_argument("this version not support.");
//    }

//  for(size_t fi = 0; fi < tri.size(2); ++fi) {
//      constraints->push_back(
//            std::shared_ptr<jtf::function::functionN1_t<double,int32_t> >(
//              new align_sh_inner_no_quadratic(tri.size(2), fi, &normal(0,fi))));
//    }
//}

hj::function::function_t<double, int32_t> *
build_frame_fix_func(const matrixst & tet, const matrixd & node,
                     const matrixst & aligned_face_idx,
                     const matrixd & aligned_frame,
                     std::shared_ptr<jtf::mesh::face2tet_adjacent> & fa,
                     const double weight)
{
  if(!fa.get())
    fa.reset(jtf::mesh::face2tet_adjacent::create(tet));

  std::shared_ptr<vector<std::shared_ptr<function_t<double, int32_t> > > > funcs(
        new vector<std::shared_ptr<function_t<double, int32_t> > > );

  vector<double> area_of_surface_tri(aligned_face_idx.size());
  for(size_t fi = 0; fi < area_of_surface_tri.size(); ++fi)
    area_of_surface_tri[fi] = fabs(
          jtf::mesh::cal_face_area(&fa->faces_[aligned_face_idx[fi]][0] , fa->faces_[aligned_face_idx[fi]].size() , node));
  const double total_area = accumulate(area_of_surface_tri.begin(), area_of_surface_tri.end(), 0.0);

  for(size_t ai = 0; ai < aligned_face_idx.size(); ++ai){
      const double w = sqrt(area_of_surface_tri[ai]/total_area * weight);

      const pair<size_t,size_t> & tet_pair = fa->face2tet_[aligned_face_idx[ai]];
      assert(fa->is_outside_face(tet_pair));

      const size_t & tet_idx = (tet_pair.first == -1?tet_pair.second:tet_pair.first);

      funcs->push_back(
            std::shared_ptr<frame_fix_func>(
              new frame_fix_func(tet.size(2), aligned_frame(colon(),ai), tet_idx, w)));
    }
  return new_catenated_function<double, int32_t>(funcs);
}

hj::function::function_t<double, int32_t> *
add_box_aligned_function(const matrixst & tet, const control_unit::UNIT_TYPE & unit,
                         const double weight)
{
  std::shared_ptr<vector<std::shared_ptr<function_t<double, int32_t> > > > funcs(
        new vector<std::shared_ptr<function_t<double, int32_t> > > );
  matrixd box_sh(9,1);

  for(size_t bi = 0 ; bi < unit.size(); ++bi){
      calc_rot_cubic_f_sh_(&box_sh[0], &get<0>(unit[bi])[0]);
      for(boost::unordered_set<size_t>::const_iterator bucit =
          get<1>(unit[bi]).begin(); bucit != get<1>(unit[bi]).end();
          ++bucit){
          funcs->push_back(
                std::shared_ptr<function_t<double, int32_t> >(
                  new frame_fix_func(tet.size(2), *bucit, box_sh, weight)));
        }
    }
  return new_catenated_function<double, int32_t>(funcs);
}


hj::function::function_t<double, int32_t> *
add_plane_aligned_function(const matrixst & tet, const control_unit::UNIT_TYPE & unit,
                           const double weight)
{
  std::shared_ptr<vector<std::shared_ptr<function_t<double, int32_t> > > > funcs(
        new vector<std::shared_ptr<function_t<double, int32_t> > > );

  for(size_t bi = 0; bi < unit.size(); ++bi){
      //plane_sh[bi].resize(9,1);
      //calc_rot_cubic_f_sh_(&plane_sh[bi][0], &plane_vec[bi].get<0>()[0]);
      for(boost::unordered_set<size_t>::const_iterator bucit =
          get<1>( unit[bi]).begin(); bucit != get<1>(unit[bi]).end();
          ++bucit){
          funcs->push_back(
                std::shared_ptr<function_t<double, int32_t> >(
                  new align_sh_inner(tet.size(2), *bucit, &get<0>(unit[bi])[0], weight)));
        }
    }
  return new_catenated_function<double, int32_t>(funcs);
}

hj::function::function_t<double, int32_t> *
add_line_aligned_function(const matrixst & tet, const control_unit::UNIT_TYPE & unit,
                          const double weight)
{
  std::shared_ptr<vector<std::shared_ptr<function_t<double, int32_t> > > > funcs(
        new vector<std::shared_ptr<function_t<double, int32_t> > > );

  for(size_t bi = 0; bi < unit.size(); ++bi){
      // line_sh[bi].resize(3,1);
      //calc_rot_cubic_f_sh_(&line_sh[bi][0], &line_vec[bi].get<0>()[0]);
      for(boost::unordered_set<size_t>::const_iterator bucit =
          get<1>(unit[bi]).begin(); bucit != get<1>(unit[bi]).end();
          ++bucit){
          funcs->push_back(
                std::shared_ptr<function_t<double, int32_t> >(
                  new align_sh_inner(tet.size(2), *bucit, &get<0>(unit[bi])[0],1.0)));
        }
    }
  return new_catenated_function<double, int32_t>(funcs);
}

hj::function::function_t<double, int32_t> *
build_control_function(const matrixst & tet, const matrixd &node,
                       const control_unit &cu, const double control_weight)
{
  std::shared_ptr<vector<std::shared_ptr<function_t<double, int32_t> > > > funcs(
        new vector<std::shared_ptr<function_t<double, int32_t> > > );

  for(size_t cui = 0; cui < cu.units_.size(); ++cui){
      if(cu.units_[cui].first == "box" && cu.units_[cui].second.get()){
          funcs->push_back(std::shared_ptr<function_t<double, int32_t> >(
                             add_box_aligned_function(tet, *(cu.units_[cui].second), control_weight)));
          continue;
        }
      if(cu.units_[cui].first == "plane" && cu.units_[cui].second.get()){
          funcs->push_back(std::shared_ptr<function_t<double, int32_t> >(
                             add_plane_aligned_function(tet, *(cu.units_[cui].second), control_weight)));
          continue;
        }
      if(cu.units_[cui].first == "line" && cu.units_[cui].second.get()){
          funcs->push_back(std::shared_ptr<function_t<double, int32_t> >(
                             add_line_aligned_function(tet, *(cu.units_[cui].second), control_weight)));
          continue;
        }
    }
  if(funcs->empty()) return 0;
  return new_catenated_function<double, int32_t>(funcs);
}


int cal_rotation(const matrixd &d1, const matrixd &d2, double *rot)
{
  if(d1.size() != d2.size()){
      cerr << "# [error] two direction do not have the same length." << std::endl;
      return __LINE__;
    }
  itr_matrix<double*> rot_mat(3,3,rot);
  matrixd cross_v = zjucad::matrix::cross(d1,d2);
  matrixd d2_temp = d2;
  double cross_norm = zjucad::matrix::norm(cross_v);
  double dot_norm = zjucad::matrix::dot(d1,d2);
  if(cross_norm < 1e-6 && dot_norm > 0){ // parallel, same direction
      rot_mat = zjucad::matrix::eye<double>(3);
      return 0;
    }else if(cross_norm < 1e-6 && dot_norm < 0){
      const static matrixd eye_static = eye<double>(3);
      for(size_t i = 0; i < 3; ++i){
          cross_v = zjucad::matrix::cross(eye_static(colon(),i),d1);
          if(zjucad::matrix::norm(cross_v) > 1e-6){
              d2_temp = eye_static(colon(),i);
              dot_norm = zjucad::matrix::dot(d1,d2_temp);
              break;
            }
        }
    }

  const double len_d1 = norm(d1);
  const double len_d2 = norm(d2_temp);
  assert(len_d1 > 1e-6);
  assert(len_d2 > 1e-6);
  const double angle = jtf::math::safe_acos(dot_norm/(len_d1*len_d2));
  static zjucad::matrix::matrix<double> M33 = eye<double>(3);
  from_angle_to_rotation_matrix(angle, cross_v, M33);
  rot_mat =  M33;
  return 0;
}

inline bool is_acute_triangle(const zjucad::matrix::matrix<size_t> & one_tri,
                              const zjucad::matrix::matrix<double> & node)
{
  assert(one_tri.size() == 3);
  double edge_len[] = {
    zjucad::matrix::norm(node(colon(), one_tri[0]) - node(colon(), one_tri[1])),
    zjucad::matrix::norm(node(colon(), one_tri[1]) - node(colon(), one_tri[2])),
    zjucad::matrix::norm(node(colon(), one_tri[2]) - node(colon(), one_tri[0]))
  };
  sort(&edge_len[0], &edge_len[0]+3);
  if(fabs(edge_len[0]*edge_len[0]+edge_len[1]*edge_len[1]-edge_len[2]*edge_len[2]) < 1e-6)
    return false;
  if(edge_len[0]*edge_len[0]+edge_len[1]*edge_len[1]>edge_len[2]*edge_len[2])
    return true;
  else return false;
}

inline double cal_tri_area(const double e1, const double e2, const double e3)
{
  const double p = (e1+e2+e3)/2;
  //assert(p > e1 && p > e2 && p > e3);
  return sqrt(p*(p-e1)*(p-e2)*(p-e3));
}

inline double cal_circum_radius(const zjucad::matrix::matrix<size_t> &one_tri,
                                const zjucad::matrix::matrix<double> & node)
{
  assert(one_tri.size() == 3);
  const double edge_len[] = {
    zjucad::matrix::norm(node(colon(), one_tri[0]) - node(colon(), one_tri[1])),
    zjucad::matrix::norm(node(colon(), one_tri[1]) - node(colon(), one_tri[2])),
    zjucad::matrix::norm(node(colon(), one_tri[2]) - node(colon(), one_tri[0]))
  };
  const double p = (edge_len[0] + edge_len[1] + edge_len[2])/2.0;
  const double area = sqrt(p*(p-edge_len[0])*(p-edge_len[1])*(p-edge_len[2]));
  return edge_len[0]*edge_len[1]*edge_len[2]/(4*area);
}

int generate_cross_field_according_to_angle_defet(
    const jtf::mesh::edge2cell_adjacent &ea,
    const matrix<size_t> & tri,
    const matrix<double> & node,
    const matrix<double> & angle_defect,
    const matrix<double> & normal,
    const char * fv_file)
{
  jtf::mesh::one_ring_face_at_point orfap;
  orfap.add_all_faces(tri,ea);
  orfap.sort_int_loop_with_normal_info(tri, node, ea, normal);

  matrix<matrix<double> > cross_v(tri.size(2),1);
  cross_v[0].resize(3,2);
  cross_v[0](colon(),0) = node(colon(),tri(1,0)) - node(colon(),tri(0,0));
  cross_v[0](colon(),0) /= norm(cross_v[0](colon(),0));

  manually_set_angle_defect msae;
  vector<bool> visited_face(tri.size(2), false);

  stack<size_t> need_to_assign_face;
  need_to_assign_face.push(0);
  size_t common_edge[2];
  matrix<double> rot1(3,3), rot2(3,3);
  while(find(visited_face.begin(), visited_face.end(), false) != visited_face.end()){
      size_t one_face_idx = need_to_assign_face.top();
      need_to_assign_face.pop();

      visited_face[one_face_idx] = true;

      for(size_t pi = 0; pi < tri.size(1); ++pi){
          const size_t edge_idx = ea.get_edge_idx(tri(pi,one_face_idx),
                                                  tri((pi+1)%tri.size(1), one_face_idx));
          const pair<size_t,size_t> & face_pair = ea.edge2cell_[edge_idx];
          const size_t other_face_idx = face_pair.first + face_pair.second - one_face_idx;
          if(visited_face[other_face_idx]) continue;

          common_edge[0] = ea.edges_[edge_idx].first;
          common_edge[1] = ea.edges_[edge_idx].second;
          msae.orient_edge(common_edge, tri(colon(), one_face_idx));

          // two rotations
          // first, rotate around normal of one_face_idx, with angle_defect
          double angle_this_edge = angle_defect[edge_idx];
          if(face_pair.first == other_face_idx && face_pair.second == one_face_idx)
            angle_this_edge *= -1;
          from_angle_to_rotation_matrix(angle_this_edge, normal(colon(), one_face_idx), rot1);

          double dia_angle =
              jtf::math::safe_acos(dot(normal(colon(), one_face_idx),normal(colon(), other_face_idx))/
                                   (norm(normal(colon(),one_face_idx))*norm(normal(colon(),other_face_idx))));
          matrix<double> n1xn2 = cross(normal(colon(), one_face_idx), normal(colon(), other_face_idx));
          matrix<double> common_edge_ref = node(colon(), common_edge[1]) - node(colon(), common_edge[0]);
          if(dot(n1xn2, common_edge_ref) < 0) dia_angle *= -1;
          from_angle_to_rotation_matrix(dia_angle, common_edge_ref, rot2);

          cross_v[other_face_idx].resize(3,2);
          cross_v[other_face_idx](colon(),0) = rot2 * temp(rot1 * cross_v[one_face_idx](colon(),0));

          visited_face[other_face_idx] = true;
          need_to_assign_face.push(other_face_idx);
        }
    }

  for(size_t fi = 0; fi < tri.size(2); ++fi){
      cross_v[fi](colon(),1) = cross(normal(colon(), fi), cross_v[fi](colon(),0));
    }

  ofstream ofs(fv_file);
  if(ofs.fail()){
      cerr << "# [error] can not open fv file." << endl;
      return __LINE__;
    }

  ofs << tri.size(2) << endl;
  ofs << endl;
  for(size_t fi = 0; fi < tri.size(2) ; ++fi){
      for(size_t j = 0; j < 2; ++j){
          for(size_t i = 0;  i < tri.size(1); ++i){
              ofs << cross_v[fi](i,j) << " ";
            }
          ofs << endl;
        }
      ofs << endl;
    }

  return 0;
}

shared_ptr<jtf::function::functionN1_t<double, int32_t> >
build_frame_surface_smooth_function_jtf_tri(
    const size_t cell_num,
    const matrix<size_t> & tri,
    const matrix<double> & node,
    const matrix<double> & normal,
    const jtf::mesh::edge2cell_adjacent &ea,
    const std::string strategy,
    const double w,
    boost::property_tree::ptree &pt,
    const zjucad::matrix::matrix<size_t> * tri_idx,
    const jtf::mesh::face2tet_adjacent * fa)
{
  if(tri_idx != 0 && fa == 0){
      throw std::invalid_argument("for tet mesh, tri_idx and fa should be provided togather.");
    }
  ///init
  // \|R_z*R_n1_n2*f1-f2\|^2
  zjucad::matrix::matrix<double> R12(9, ea.edges_.size());
  for(size_t ei = 0; ei < ea.edge2cell_.size(); ++ei){
      const pair<size_t,size_t> & p2p = ea.edge2cell_[ei];
      if(ea.is_boundary_edge(p2p)) continue;
      cal_rotation(normal(colon(), p2p.first), normal(colon(), p2p.second), &R12(0,ei));
    }

  matrix<double> area(tri.size(2),1);
  for(size_t fi = 0; fi < tri.size(2); ++fi){
      area[fi] = jtf::mesh::cal_face_area(&tri(0,fi), tri.size(1), node);
    }
  const double total_area = std::accumulate(area.begin(), area.end(), 0.0);

  /////////////////////////////////////////////////////////////////////////

  shared_ptr<vector<std::shared_ptr<jtf::function::functionN1_t<double,int32_t> > > > funcs(
        new vector<std::shared_ptr<jtf::function::functionN1_t<double,int32_t> > >);

#define angle_DEFECT 0
#if angle_DEFECT == 0 // trivial connection
  matrix<double> angle = zeros<double>(ea.edges_.size(),1);
  matrix<double> rot(3,3),zyz(3,1);
  std::unique_ptr<angle_defect> ad(new manually_set_angle_defect);
  ad->opt(ea,tri, node, angle,pt);

  {
    ofstream ofs("line_angle.vtk");
    matrix<size_t> edges(2, ea.edge2cell_.size());
    for(size_t ei = 0; ei < ea.edges_.size(); ++ei){
        edges(0,ei) = ea.edges_[ei].first;
        edges(1,ei) = ea.edges_[ei].second;
      }
    line2vtk(ofs, &node[0], node.size(2), &edges[0], edges.size(2));
    cell_data(ofs, &angle[0], angle.size(), "angle");
    generate_cross_field_according_to_angle_defet(ea, tri, node, angle, normal, "test_fv");
  }

  {
    ofstream ofs("point_angle.vtk");
    matrix<double> angle_point(node.size(2),1);
    size_t common_edge[2];
    jtf::mesh::one_ring_face_at_point orfap;
    orfap.add_all_faces(tri, ea);
    orfap.sort_int_loop_with_normal_info(tri, node, ea, normal);

    for(jtf::mesh::one_ring_face_at_point::p2f_type::const_iterator cit = orfap.p2f_.begin();
        cit != orfap.p2f_.end(); ++cit){
        const vector<size_t> & faces_loop = cit->second;
        assert(faces_loop.front() == faces_loop.back());
        double angle_this_point = 0;
        for(size_t fi = 0; fi < faces_loop.size()-1;++fi){
            jtf::mesh::find_common_edge(tri(zjucad::matrix::colon(), faces_loop[fi]),
                                        tri(zjucad::matrix::colon(), faces_loop[fi+1]),
                common_edge);
            const size_t edge_idx = ea.get_edge_idx(common_edge);
            const pair<size_t,size_t> & face_pair = ea.edge2cell_[edge_idx];
            if(face_pair.first == faces_loop[fi+1] && face_pair.second == faces_loop[fi]){
                angle_this_point += angle[edge_idx] * -1.0;
              }else
              angle_this_point += angle[edge_idx];
          }
        angle_point[cit->first] = angle_this_point;
      }

    tri2vtk(ofs, &node[0], node.size(2), &tri[0], tri.size(2));
    point_data(ofs, &angle_point[0], angle_point.size(), "point_angle");
  }

  jtf::mesh::one_ring_face_at_point orfap;
  orfap.add_all_faces(tri,ea);
  orfap.sort_int_loop_with_normal_info(tri, node, ea, normal);
  matrix<double> rn1n2(3,3);

  for(size_t ei = 0 ; ei < ea.edges_.size(); ++ei){
      std::copy(R12(colon(),ei).begin(), R12(colon(),ei).end(), rn1n2.begin());
      const pair<size_t,size_t> & two_face_pair = ea.edge2cell_[ei];
      from_angle_to_rotation_matrix(angle[ei],normal(colon(),two_face_pair.second), rot);
      rot = temp(rot * rn1n2);
      rotation_matrix_2_zyz_angle(&rot[0], &zyz[0], 0);
      double weight = (area[two_face_pair.first] + area[two_face_pair.second])/(3.0*total_area);
      pair<size_t,size_t> idx_pair = two_face_pair;
      if(tri_idx != 0 && fa != 0){
          const pair<size_t,size_t> &tet_pair_0 = fa->face2tet_.at((*tri_idx)[idx_pair.first]);
          const pair<size_t,size_t> &tet_pair_1 = fa->face2tet_.at((*tri_idx)[idx_pair.second]);
          assert(fa->is_outside_face(tet_pair_0) && fa->is_outside_face(tet_pair_1));
          idx_pair.first = (tet_pair_0.first==-1?tet_pair_0.second:tet_pair_0.first);
          idx_pair.second = (tet_pair_1.first==-1?tet_pair_1.second:tet_pair_1.first);
        }
      funcs->push_back(shared_ptr<jtf::function::functionN1_t<double,int32_t> >(
                         new surface_smooth_jtf_func(zyz, idx_pair, cell_num, w*weight)));
    }
#elif angle_DEFECT == 1

  matrix<double> dual_edge_area = zeros<double>(ea.edges_.size(),1);
  map<size_t,double> angle_defect_map;
  {
    jtf::mesh::one_ring_face_at_point orfap;
    orfap.add_all_faces(tri, ea);
    naive_angle_defect nad;
    nad.cal_angle_defect(orfap, tri, node, angle_defect_map);
    if(angle_defect_map.size() != node.size(2)){
        throw std::logic_error("angle_defect_map should have the same size with nodes");
      }
    double total_angle = 0;
    for(const auto & one_point : angle_defect_map){
        total_angle += one_point.second;
      }
    cerr << "# [info] total angle_defect " << total_angle << endl;
  }

  matrix<double> Radius(tri.size(2),1);
  matrix<bool> is_acute_tri(tri.size(),1);
  for(size_t i = 0 ; i < tri.size(2); ++i){
      Radius[i] = cal_circum_radius(tri(colon(),i), node);
      is_acute_tri[i] = is_acute_triangle(tri(colon(),i), node);
    }

  map<size_t,double> cell_area4prime_vertex;
  for(size_t ei = 0 ; ei < ea.edges_.size(); ++ei){
      const pair<size_t,size_t> & face_pair = ea.edge2cell_[ei];
      assert(!ea.is_boundary_edge(face_pair));
      const double edge_len = norm(node(colon(), ea.edges_[ei].first) - node(colon(), ea.edges_[ei].second));
      if(is_acute_tri[face_pair.first]){
          dual_edge_area[ei] += cal_tri_area(Radius[face_pair.first], Radius[face_pair.first], edge_len);
        }else{
          dual_edge_area[ei] += area[face_pair.first]/3.0;
        }
      if(is_acute_tri[face_pair.second]){
          dual_edge_area[ei] += cal_tri_area(Radius[face_pair.second], Radius[face_pair.second], edge_len);
        }else
        dual_edge_area[ei] += area[face_pair.second]/3.0;

      auto it = cell_area4prime_vertex.find(ea.edges_[ei].first);
      if(it != cell_area4prime_vertex.end()){
          it->second += dual_edge_area[ei]/2.0;
        }else
        cell_area4prime_vertex[ea.edges_[ei].first] = dual_edge_area[ei]/2.0;
      it = cell_area4prime_vertex.find(ea.edges_[ei].second);
      if(it != cell_area4prime_vertex.end()){
          it->second += dual_edge_area[ei]/2.0;
        }else
        cell_area4prime_vertex[ea.edges_[ei].second] = dual_edge_area[ei]/2.0;
    }

  double weight = 1;
  double angle_ = 0;
  matrix<double> rot = eye<double>(3);
  matrix<double> zyz(3,1);
  double total_angle = 0;
  for(size_t ei = 0; ei < ea.edge2cell_.size(); ++ei){
      pair<size_t,size_t> face_pair = ea.edge2cell_[ei];
      pair<size_t,size_t> &idx_pair = face_pair;
      const pair<size_t,size_t> & one_edge = ea.edges_[ei];
      if(ea.is_boundary_edge(face_pair)) continue;

      itr_matrix<double*> rn1n2(3,3, &R12(0,ei));
      weight = dual_edge_area[ei]/(2*total_area);

      if(tri_idx != 0 && fa != 0){
          const pair<size_t,size_t> &tet_pair_0 = fa->face2tet_.at((*tri_idx)[idx_pair.first]);
          const pair<size_t,size_t> &tet_pair_1 = fa->face2tet_.at((*tri_idx)[idx_pair.second]);
          assert(fa->is_outside_face(tet_pair_0) && fa->is_outside_face(tet_pair_1));
          idx_pair.first = (tet_pair_0.first==-1?tet_pair_0.second:tet_pair_0.first);
          idx_pair.second = (tet_pair_1.first==-1?tet_pair_1.second:tet_pair_1.first);
        }
      {
        angle_ = angle_defect_map.at(one_edge.first)*dual_edge_area[ei]/(cell_area4prime_vertex.at(one_edge.first)*2.0);
        total_angle += angle_;
        from_angle_to_rotation_matrix(angle_,normal(colon(),face_pair.second), rot);
        rot = temp(rot * rn1n2);
        rotation_matrix_2_zyz_angle(&rot[0], &zyz[0], 0);
        funcs->push_back(shared_ptr<hj::function::function_t<double,int32_t> >(
                           new surface_smooth_func(zyz, idx_pair, tri.size(2), sqrt(w*weight))));
      }

      {
        angle_ = angle_defect_map.at(one_edge.second)*dual_edge_area[ei]/(cell_area4prime_vertex.at(one_edge.second)*2.0);
        total_angle += angle_;
        from_angle_to_rotation_matrix(angle_,normal(colon(),face_pair.second), rot);
        rot = temp(rot * rn1n2);
        rotation_matrix_2_zyz_angle(&rot[0], &zyz[0], 0);
        //        funcs->push_back(shared_ptr<jtf::function::functionN1_t<double,int32_t> >(
        //                           new surface_smooth_jtf_func(zyz, face_pair, tri.size(2), w*weight)));
        funcs->push_back(shared_ptr<hj::function::function_t<double,int32_t> >(
                           new surface_smooth_func(zyz, idx_pair, cell_num, sqrt(w*weight))));
      }
    }

  cerr << "total angle defect " << total_angle/My_PI() << "Pi" << endl;

#endif

  //return shared_ptr<hj::function::function_t<double,int32_t> >(new_catenated_function<double,int32_t>(funcs));
  return shared_ptr<jtf::function::functionN1_t<double,int32_t> >(
        new jtf::function::sum_function<double,int32_t, jtf::function::SMART_STD>(*funcs));
}
