#include <iostream>
#include <fstream>
#include <string>
#include <memory>
#include <boost/algorithm/string.hpp>
#include <hjlib/sparse/fast_AAT.h>
#include <hjlib/sparse/format.h>
#include <hjlib/sparse/operation.h>

#include <zjucad/matrix/matrix.h>
#include <hjlib/sparse/sparse.h>
#include <jtflib/function/func_aux.h>
#include <jtflib/function/operation.h>
#include <jtflib/mesh/io.h>

#include "../common/util.h"
#include "../common/io.h"
#include "../../common/IO.h"
#include "../../hex_param/io.h"
#include "../../common/zyz.h"
#include "../../common/vtk.h"

#include "../common/mesh_topology.h"

#include "descriptor_vol.h"
#include "func_terms/arap.h"
#include "func_terms/frame_align.h"
#include "func_terms/linear_equation.h"
#include "func_terms/degree_smooth.h"
#include "func_terms/normal_diff.h"
#include "func_terms/l1-normal.h"
#include "func_terms/area_sum.h"
#include "func_terms/normal-anti-flip.h"
#include "func_terms/vol-anti-flip.h"
#include "func_terms/quality_check.h"

#include "../../equation_graph/equation_graph.h"
using namespace std;
using namespace zjucad::matrix;
using namespace hj::sparse;

double linear_equation_hj_func::weight_ = 1.0;
double smooth_L1::scalar_ = 0.05;
double smooth_L1::L1_sqrt_eps_ = 0.5;

////////////////////////////////////////////////////////////////////////////////
/// \brief least_square_Gauss_Newton::least_square_Gauss_Newton
/// \param f
///
//! assume the dim and pattern of f is constant


////////////////////////////////////////////////////////////////////////////////
/// descriptor
////////////////////////////////////////////////////////////////////////////////

void construct_node_mapping(
    const zjucad::matrix::matrix<size_t> & cut_tet,
    const zjucad::matrix::matrix<size_t> & uncut_tet,
    const zjucad::matrix::matrix<size_t> & cut_tet2tet,
    hj::sparse::csc<double,int32_t> & eqn,
    std::set<size_t> & zero_set,
    const boost::unordered_map<pair<size_t,size_t>,size_t> &inner_type,
    const boost::unordered_map<size_t,size_t> & surface_type,
    int is_restricted_type,
    matrix<size_t> & cut_face,
    std::vector<std::vector<std::pair<size_t,size_t> > > & cut_patches,
    std::vector<size_t> & rot_type,
    vector<vector<size_t> > & node_group_vec,
    vector<bool> &integer_group_flag ,
    const zjucad::matrix::matrix<double> &node)
{
  unique_ptr<transition_elimination> te(
        transition_elimination::create_with_gap(
          uncut_tet, cut_tet, cut_tet2tet, inner_type,
          surface_type,(is_restricted_type==0?true:false), node));

  te->get_node_comp(eqn);
  te->get_zero_idx(zero_set);
  te->get_ordered_cut_face_patches(cut_face, cut_patches, rot_type);
  vector<group<size_t> > node_group = te->out();
  vector<size_t> integer_group_idx = te->get_integer_group_idx();

  vector<size_t> one_group_vec;
  //for(const auto & one_group : node_group){
  for(size_t gi = 0; gi < node_group.size(); ++gi){
      const auto & one_group = node_group[gi];
      if(one_group.empty()) continue;
      one_group_vec.resize(one_group.size());
      copy(one_group.begin(), one_group.end(), one_group_vec.begin());
      node_group_vec.push_back(one_group_vec);

      if(find(integer_group_idx.begin(), integer_group_idx.end(), gi)
         != integer_group_idx.end())
        integer_group_flag.push_back(true);
      else
        integer_group_flag.push_back(false);
    }


  // write out
  {// save node_groups
    ofstream ofs("node_group");
    assert(node_group_vec.size() == integer_group_flag.size());
    for(size_t i = 0 ; i < node_group_vec.size(); ++i){
        ofs << "g " << i << " " << node_group_vec[i].size() << " " << integer_group_flag[i] << endl;
        for(size_t vi = 0; vi < node_group_vec[i].size(); ++vi){
            ofs << node_group_vec[i][vi] << " ";
          }
        ofs << endl;
      }
  }
}

int descriptor_vol::init(const zjucad::matrix::matrix<size_t> &mesh,
                         const zjucad::matrix::matrix<double> &node,
                         boost::property_tree::ptree &pt)
{
  bi_.cut_tet = mesh;
  bi_.cut_node = node;
  cb_ = nullptr;
  wf_.clear();

  bi_.fa.reset(jtf::mesh::face2tet_adjacent::create(mesh));
  if(!bi_.fa.get()){
      cerr << "# [error] can not build face2tet_adjacent." << endl;
      return __LINE__;
    }

  sim_ = TET;
  jtf::mesh::get_outside_face(*(bi_.fa), bi_.outside_face,true);
  jtf::mesh::get_outside_face_idx(*(bi_.fa), bi_.outside_face_idx);

  bi_.vol_weight = zjucad::matrix::ones<double>(mesh.size(2),1);

  obj_vec_.clear();
  obj_type_.clear();

  eqn_cons_.clear();
  eqn_cons_type_.clear();

  ineqn_cons_.clear();
  ineqn_cons_type_.clear();

  bi_.vol_func_.clear();

  if(bi_.uncut_tet.size()  ==  0){
      if(jtf::mesh::tet_mesh_read_from_zjumat(pt.get<string>("input/uncut_tet.value").c_str(),
                                              &bi_.uncut_node, &bi_.uncut_tet))
        return __LINE__;
    }
  bi_.fa_uncut.reset(jtf::mesh::face2tet_adjacent::create(bi_.uncut_tet));
  jtf::mesh::get_outside_face(*bi_.fa_uncut, bi_.outside_face_uncut);
  if(!bi_.fa_uncut.get()){
      throw std::logic_error("can not build face2tet adjacent.");
    }

  bi_.cut_tet2tet.resize(max(mesh)+1);
  bi_.cut_tet2tet(mesh) = bi_.uncut_tet(colon());

  construct_original_face_in_cut_mesh(bi_.uncut_tet, bi_.uncut_node,
                                      *bi_.fa_uncut,mesh, *(bi_.fa),
                                      bi_.orig_face_in_cut, bi_.cut_tet2tet);

  bi_.ea_orig_in_cut.reset(jtf::mesh::edge2cell_adjacent::create(bi_.orig_face_in_cut));
  const string need_node_compression =
      pt.get<string>("input/need_node_compression.value").c_str();
  if(need_node_compression == "yes" || need_node_compression == "YES" ||
     need_node_compression == "y" || need_node_compression == "Y"){
      if(!zjucad::has("input/node_compression.value",pt)){
          if(bi_.inner_type.empty()){
              if(load_inner_face_jump_type(pt.get<string>("input/inner_type.value").c_str(), bi_.inner_type))
                return __LINE__;
            }
          int is_restricted_type = 0;
          if(bi_.surface_type.empty()){
              if(zjucad::has("input/surface_type.value", pt)){
                  is_restricted_type = load_surface_type(pt.get<string>("input/surface_type.value").c_str(),
                                                         bi_.surface_type, bi_.fa_uncut.get());
                  if(is_restricted_type == 1) { // normal align type
                      convert_surface_normal_type2_restricted_type(bi_.surface_type);
                      is_restricted_type = 0;
                    }
                  cerr << "# [info] find surface type." << endl;
                }else{
                  cerr << "# [info] no surface type." << endl;
                }
            }
          construct_node_mapping(mesh, bi_.uncut_tet, bi_.cut_tet2tet,
                                 bi_.NM.ZT, bi_.zero_index_,
                                 bi_.inner_type, bi_.surface_type,is_restricted_type,
                                 bi_.gap_faces, bi_.cut_patches,
                                 bi_.rot_type, bi_.node_group, bi_.integer_group_flag,
                                 node);
          bi_.NM.q = zeros<double>(bi_.NM.ZT.size(2),1); // original q is zero
          bi_.NM.is_fixed = zeros<int>(bi_.NM.ZT.size(2),1);
        }else{
          //          if(bi_.load_node_comp_eqn(
          //               pt.get<string>("input/node_compression.value").c_str(), bi_.node_mapping))
          //            throw std::invalid_argument("invalide node_compression_equation");

          //          if(bi_.load_zero_index(
          //               pt.get<string>("input/zero_index.value").c_str(), bi_.zero_index_))
          //            throw std::invalid_argument("invalide zero_index");

        }
    }

  return 0;
}

void descriptor_vol::set_objective(
    const zjucad::matrix::matrix<size_t> & mesh,
    const zjucad::matrix::matrix<double> & node,
    boost::property_tree::ptree &pt)
{
  pt.put("object/type.desc", "[arap/asap/frame/integer/len/degree/equation/l1-normal]");
  const string obj_type = pt.get<string>("object/type.value");
  vector<string> types;
  boost::split(types, obj_type, boost::is_any_of(","));
  obj_vec_.clear();
  obj_type_.clear();
  for(const auto &obj_type_term: types)
    switch(str2int(obj_type_term.c_str())) {
      case str2int("arap"):
        {
          hj_func_cons_ptr hfptr(
                build_arap_func(mesh, node, sim_, has_node_mapping()?&get_node_mapping():0));

          obj_vec_.push_back(jtf_func_ptr(jtf::function::least_square_warpper(hfptr)));
          obj_type_["arap"] = obj_vec_.back();
          cerr << "# [info] add object/type: arap" << endl;
          wf_.push_back(make_pair(obj_vec_.back().get(),1));

          break;
        }
      case str2int("frame"):
        {
          if(bi_.frames.size() == 0){
              zjucad::matrix::matrix<double> zyz;
              if(read_zyz(pt.get<string>("input/zyz.value").c_str(), zyz)){
                  throw std::invalid_argument("invalid zyz.");
                }
              bi_.frames.resize(zyz.size(2),1);
              for(size_t ti = 0; ti < zyz.size(2); ++ti){
                  bi_.frames[ti].resize(3,3);
                  zyz_angle_2_rotation_matrix1(&zyz(0,ti), &bi_.frames[ti][0]);
                }
            }
          assert(bi_.frames.size() == mesh.size(2));
          build_frame_align_func(mesh, node, bi_.frames, bi_, bi_.vol_func_, 0, TET,
                                 (has_node_mapping()?&(get_node_mapping()):0));
          hj_func_cons_ptr all_func(
                hj::function::new_catenated_function<double,int32_t>(bi_.vol_func_));

          obj_vec_.push_back(jtf_func_ptr(jtf::function::least_square_warpper(all_func)));

          cerr << "# [info] add object/type: frame" << endl;
          obj_type_["frame"] = obj_vec_.back();
          break;
        }
      case str2int("vol-anti-flip"):
        {
          const double w = pt.get<double>("weight/vol-anti-flip.value",0.01);
          if(w < 1e-6) break;
          matrix<double> weight;
          hj_func_ptr hfptr(
                build_vol_anti_flip_hj_func(mesh, node, weight,
                                            has_node_mapping()?&get_node_mapping():0));

          weight *= w;
          obj_vec_.push_back(jtf_func_ptr(jtf::function::neg_log_warpper(hfptr, weight)));

          obj_type_["vol-anti-flip"] = obj_vec_.back();
          cerr << "# [info] add object/type: vol-anti-flip" << endl;

          break;
        }
      case str2int("surface-type"):
        {
          assert(bi_.NM.ZT.size(1) == 0 && bi_.NM.ZT.size(2) == 0);
          const double surface_type_w = pt.get<double>("weight/surface-type-align.value");
          if(bi_.surface_type.empty()){
              int is_restricted_type = load_surface_type(pt.get<string>("input/surface_type.value").c_str(),
                                                         bi_.surface_type, bi_.fa_uncut.get());
              if(is_restricted_type == 1) { // normal align type
                  convert_surface_normal_type2_restricted_type(bi_.surface_type);
                }
              cerr << "# [info] find surface type." << endl;
            }
          shared_ptr<vector<hj_func_cons_ptr> > funcs(new vector<hj_func_cons_ptr>);
          vector<size_t> idx(2); vector<double> coeff(2);
          for(size_t fi = 0; fi < bi_.orig_face_in_cut.size(2); ++fi){
              const size_t orig_face_idx =
                  bi_.fa_uncut->get_face_idx(
                    bi_.cut_tet2tet[bi_.orig_face_in_cut(0,fi)],
                  bi_.cut_tet2tet[bi_.orig_face_in_cut(1,fi)],
                  bi_.cut_tet2tet[bi_.orig_face_in_cut(2,fi)]);
              assert(orig_face_idx != -1);
              const auto it = bi_.surface_type.find(orig_face_idx);
              assert(it != bi_.surface_type.end());
              const size_t di = it->second;
              idx[0] = 3 * bi_.orig_face_in_cut(0, fi) + di;
              coeff[0] = 1.0; coeff[1] = -1.0;
              for(size_t i = 1; i < bi_.orig_face_in_cut.size(1); ++i){
                  idx[1] = 3 * bi_.orig_face_in_cut(i, fi) + di;
                  funcs->push_back(hj_func_cons_ptr(
                                     new linear_equation_hj_func(
                                       bi_.cut_node.size(),
                                       idx, coeff, 0.0)));
                }
            }
          set_linear_equation_hj_func_weight(surface_type_w);
          hj_func_cons_ptr all_func(hj::function::new_catenated_function<double,int32_t>(funcs));
          obj_vec_.push_back(jtf_func_ptr(jtf::function::least_square_warpper(all_func)));
          obj_type_["surface-type"] = obj_vec_.back();
          cerr << "# [info] add object/type: surface-type" << endl;

          break;
        }
        //      case str2int("degree"):
        //        {
        //          hj_func_cons_ptr hfptr(
        //                build_degree_smooth_func(mesh,has_node_mapping()?&get_node_mapping():0));

        //          obj_vec_.push_back(jtf_func_cons_ptr(jtf::function::least_square_warpper(hfptr)));

        //          cerr << "# [info] add objectype/type: degree" << endl;
        //          obj_type_["degree"] = obj_vec_.back();
        //          break;
        //        }
        //      case str2int("equation"):
        //        {
        //          if(bi_.variable_idx.empty()){
        //              if(load_equation_file(pt.get<string>("input/equation.value").c_str(),
        //                                    bi_.variable_idx, bi_.coefficient,
        //                                    (has_node_mapping()?&(get_node_mapping()):0)))
        //                throw std::invalid_argument("load equation file fail.");
        //            }
        //          build_linear_equation_funcs(
        //                (has_node_mapping()?max(get_node_mapping())+1:node.size()),
        //                bi_.variable_idx, bi_.coefficient, obj_vec_, true);
        //          cerr << "# [info] add object/type: equation" << endl;
        //          break;
        //        }
      case str2int("normal-diff"):{
          const double w = pt.get<double>("weight/normal-diff.value");
          if(w < 1e-6) break;

          if(bi_.orig_face_in_cut.size() == 0){
              construct_original_face_in_cut_mesh(
                    bi_.uncut_tet, bi_.uncut_node, *bi_.fa_uncut,mesh, *(bi_.fa),
                    bi_.orig_face_in_cut, bi_.cut_tet2tet);
            }

          hj_func_cons_ptr hfptr(
                build_normal_diff_constraint(
                  bi_.orig_face_in_cut, node, w,
                  (has_node_mapping()?&(get_node_mapping()):0)));

          obj_vec_.push_back(jtf_func_ptr(jtf::function::least_square_warpper(hfptr)));
          obj_type_["normal-diff"] = obj_vec_.back();
          cerr << "# [info] add object/type: normal-diff" << endl;

          break;
        }
      case str2int("l1-normal"):{
          if(bi_.orig_face_in_cut.size() == 0){
              construct_original_face_in_cut_mesh(
                    bi_.uncut_tet, bi_.uncut_node, *bi_.fa_uncut,mesh, *(bi_.fa),
                    bi_.orig_face_in_cut, bi_.cut_tet2tet);
            }

          matrix<double> face_area(bi_.orig_face_in_cut.size(2),1);
          for(size_t fi = 0; fi < bi_.orig_face_in_cut.size(2); ++fi){
              face_area[fi] = jtf::mesh::cal_face_area(bi_.orig_face_in_cut(colon(),fi), node);
            }

          jtf_func_ptr l1_normal_func(
                build_smooth_L1_area_normal(node, bi_.orig_face_in_cut, face_area,
                                            (has_node_mapping()?&get_node_mapping():0)));

          cerr << "# [info] add object/type: l1_normal" << endl;
          obj_vec_.push_back(l1_normal_func);
          obj_type_["l1-normal"] = obj_vec_.back();

          wf_.push_back(make_pair(l1_normal_func.get(),get_l1_normal_align_weight()));
          break;
        }
      case str2int("normal-anti-flip"):{
          const double w = pt.get<double>("weight/normal-anti-flip.value");
          if(w < 1e-6) break;

          if(bi_.orig_face_in_cut.size() == 0){
              construct_original_face_in_cut_mesh(
                    bi_.uncut_tet, bi_.uncut_node, *bi_.fa_uncut,mesh, *(bi_.fa),
                    bi_.orig_face_in_cut, bi_.cut_tet2tet);
            }

          matrix<double> weight;
          hj_func_ptr hfptr(
                build_adj_normal_func(
                  node, bi_.orig_face_in_cut, weight,
                  (has_node_mapping()?&(get_node_mapping()):0)));

          weight *= w;

          obj_vec_.push_back(jtf_func_ptr(jtf::function::neg_log_warpper(hfptr, weight)));

          obj_type_["normal-diff"] = obj_vec_.back();
          cerr << "# [info] add object/type: normal-anti-flip" << endl;
          break;
        }
      case str2int("zero-fix"):{
          if(has_node_mapping()){
              const hj::sparse::csc<double,int32_t> &node_mapping = get_node_mapping().ZT;
              for(size_t di = 0 ; di < 3; ++di){
                  jtf_func_ptr fptr(
                        jtf_func_ptr(
                          jtf::function::least_square_warpper(
                            hj_func_ptr(new variant_fix_hj(
                                          has_node_mapping()?node_mapping.size(1):node.size(), di, 0.0, &node_mapping, 1.0)))));
                  obj_vec_.push_back(fptr);
                }
            }else{
              for(size_t di = 0 ; di < 3; ++di){
                  jtf_func_ptr fptr(
                        jtf_func_ptr(
                          jtf::function::least_square_warpper(
                            hj_func_ptr(new variant_fix_hj(
                                          node.size(), di, 0.0, 0, 1.0)))));
                  obj_vec_.push_back(fptr);
                }
            }
          obj_type_["zero-fix"] = obj_vec_.back();
          cerr << "# [info] add object/type: zero-fix" << endl;
          break;
        }
      default: cerr << "# [error] unsupported object/type." << endl; break;
      }


  if(!bi_.zero_index_.empty()){ // after elimination, I will find thoese gap variants with
      assert(has_node_mapping());
      const node_mapping &node_mapping_ = get_node_mapping();
      add_zero_index_equation(obj_vec_, bi_.zero_index_.begin(),
                              bi_.zero_index_.end(),node_mapping_.ZT.size(1), &node_mapping_);
    }

  obj_.reset(new jtf::function::sum_function<double,int32_t, jtf::function::SMART_STD>(obj_vec_));
  wf_.push_back(make_pair(obj_.get(), (1+get_l1_normal_align_weight())));
}

void descriptor_vol::set_constraint(
    const zjucad::matrix::matrix<size_t> & mesh,
    const zjucad::matrix::matrix<double> & node,
    boost::property_tree::ptree &pt)
{
  pt.put("cons/type.desc", "[equation/group]");
  if(!zjucad::has("cons/type.value", pt)){
      cerr << "# [info] no cons/type."  << endl;
      return;
    }

  // here I use "," to separate the cons/type
  const string cons_type = pt.get<string>("cons/type.value");
  vector<string> types;
  boost::split(types, cons_type, boost::is_any_of(","));

  for(const string & type : types)
    switch(str2int(type.c_str())){
//      case str2int("equation"):
//        {
//          const size_t cons_num = eqn_cons_.size();
//          if(bi_.variable_idx.empty()){
//              if(load_equation_file(pt.get<string>("input/equation.value").c_str(),
//                                    bi_.variable_idx, bi_.coefficient,
//                                    (has_node_mapping()?&(get_node_mapping()):0)))
//                throw std::invalid_argument("load equation file fail.");
//            }
//          build_linear_equation_funcs(
//                (has_node_mapping()?max(get_node_mapping())+1:node.size()),
//                bi_.variable_idx, bi_.coefficient, eqn_cons_);
//          cerr << "# [info] equation constraints range: [" << cons_num << ","
//               << eqn_cons_.size() << ")." << endl;
//          eqn_cons_type_["equation"] = make_pair(cons_num, eqn_cons_.size());
//          break;
//        }
      case str2int("area-keep"):{
          if(bi_.orig_face_in_cut.size() == 0){
              construct_original_face_in_cut_mesh(
                    bi_.uncut_tet, bi_.uncut_node,
                    *bi_.fa_uncut,mesh, *(bi_.fa),
                    bi_.orig_face_in_cut, bi_.cut_tet2tet);
            }

          matrix<double> face_area(bi_.orig_face_in_cut.size(2),1);
          for(size_t fi = 0; fi < bi_.orig_face_in_cut.size(2); ++fi){
              face_area[fi] = jtf::mesh::cal_face_area(bi_.orig_face_in_cut(colon(),fi), node);
            }

          const double total_area = std::accumulate(face_area.begin(), face_area.end(), 0.0);

          eqn_cons_.push_back(
                jtf_func_ptr(new area_sum(
                               node.size(2),bi_.orig_face_in_cut,total_area,
                               (has_node_mapping()?&get_node_mapping():0))));
          cerr << "# [info] area-keep constraints range: [" << eqn_cons_.size()-1 << ","
               << eqn_cons_.size() << ")." << endl;
          eqn_cons_type_["area-keep"]= make_pair(eqn_cons_.size()-1, eqn_cons_.size());

          cb_.reset(new unnormalized_normal_quality_checker(
                      node, mesh, bi_.orig_face_in_cut, wf_, *(eqn_cons_.back())));
          break;
        }
      default: break;
      }
}

void descriptor_vol::add_eqn_constraint(jtf_func_ptr fc,
                                        int hard_or_soft,
                                        const char * eqn_type)
{
  if(hard_or_soft == 0){
      obj_vec_.push_back(fc);
      obj_.reset(new jtf::function::sum_function<double,int32_t,jtf::function::SMART_STD>(obj_vec_));
    }else{
      const size_t eqn_num = eqn_cons_.size();
      eqn_cons_.push_back(fc);
      const string type = (eqn_type == 0?"null":eqn_type);
      eqn_cons_type_[type] = make_pair(eqn_num, eqn_cons_.size());
    }
}


jtf_func_ptr descriptor_vol::get_objective(const std::string obj_type) const
{
  const auto & it = obj_type_.find(obj_type);
  if(it == obj_type_.end()) {
      return nullptr;
    }
  return it->second;
}

void descriptor_vol::frame_local_stiffening(
    const zjucad::matrix::matrix<size_t> & mesh,
    const zjucad::matrix::matrix<double> & orig_node,
    const zjucad::matrix::matrix<double> & real_node)
{
  jtf_func_cons_ptr frame_align_fptr = get_objective("frame");

  vector<pair<double,int> > diff_idx(mesh.size(2));

  assert(bi_.vol_func_.size() == mesh.size(2));
  zjucad::matrix::matrix<double> f(3,3);
  for(size_t ti = 0; ti < bi_.vol_func_.size(); ++ti){
      const auto & it = bi_.vol_func_[ti];
      it->val(&real_node[0], &f[0]);
      diff_idx[ti] = make_pair(norm(f), ti);
    }
  sort(diff_idx.begin(), diff_idx.end(), std::greater<std::pair<double,int> >());

  const double diff_max = diff_idx[0].first;
  const double diff_min = diff_idx[diff_idx.size()/10].first;

  // here, to tune the local stiffen, I use a smooth quadratic function
  // weight = a * (x-diff_min)^2+weight_min
  const double weight_min = 2;
  const double weight_max = 10;

  const double a = weight_max/ pow(diff_max-diff_min,2);

  for(size_t i = 0; i < diff_idx.size()/10; ++i){
      bi_.vol_weight[diff_idx[i].second] *=  a * pow(diff_idx[i].first-diff_min,2) + diff_min;
    }

  build_frame_align_func(mesh, orig_node, bi_.frames, bi_, bi_.vol_func_, &bi_.vol_weight,
                         TET,has_node_mapping()?&(get_node_mapping()):0);

  hj_func_ptr all_func(
        hj::function::new_catenated_function<double,int32_t>(bi_.vol_func_));

  frame_align_fptr.reset(jtf::function::least_square_warpper(all_func));

  obj_.reset(new jtf::function::sum_function<double,int32_t,jtf::function::SMART_STD>(obj_vec_));
}

int descriptor_vol::basic_infor::load_node_mapping(
    const char* node_mapping_file,
    zjucad::matrix::matrix<size_t> & node_mapping,
    const size_t variable_number)
{
  ifstream ifs(node_mapping_file);
  if(ifs.fail()){
      cerr << "# [error] can not open node_mapping_file." << endl;
      return __LINE__;
    }
  node_mapping.resize(variable_number,1);

  for(size_t i = 0; i < node_mapping.size(); ++i)
    node_mapping[i] = i;

  string temp;
  vector<string> idx_str;
  size_t st_v;
  while(!ifs.eof()){
      getline(ifs, temp); // group name or other things
      if(ifs.eof()) break;
      getline(ifs, temp); // variables in this group
      idx_str.clear();
      boost::split(idx_str, temp, boost::is_any_of(" "));
      if(idx_str.back() == "") idx_str.pop_back();
      if(idx_str.size() < 2) continue;
      const size_t first_v = std::stoi(idx_str.front());
      for(size_t i = 1; i < idx_str.size() ; ++i){
          const size_t v_idx = std::stoi(idx_str[i]);
          node_mapping[v_idx] = node_mapping[first_v];
        }
    }

  set<size_t> independent_variables(node_mapping.begin(), node_mapping.end());
  map<size_t,size_t> to_iv; // to independent variables
  size_t idx_ = 0;
  for(const auto & idx : independent_variables){
      to_iv[idx] = idx_++;
    }

  for(size_t i = 0; i < node_mapping.size(); ++i){
      node_mapping[i] = to_iv[node_mapping[i]];
    }

  return 0;
}

int descriptor_vol::basic_infor::load_zero_index(const char *zero_idx_file,
                                                 std::set<size_t> &zero_idx_set)
{
  zero_idx_set.clear();
  ifstream ifs(zero_idx_file);
  if(ifs.fail()){
      cerr << "# [error] can not load zero index file." << endl;
      return __LINE__;
    }
  size_t idx;
  while(!ifs.eof()){
      ifs >> idx;
      zero_idx_set.insert(idx);
    }
  return 0;
}


int descriptor_vol::basic_infor::load_node_comp_eqn(
    const char *node_mapping_file,
    hj::sparse::csc<double,int32_t> & node_mapping)
{
  ifstream ifs(node_mapping_file);
  if(ifs.fail()){
      cerr << "# [error] can not load node compression equation." << endl;
      return __LINE__;
    }
  size_t vn, ivn;
  ifs >> vn >> ivn;

  vector<vector<pair<double,size_t> > > eqn;
  eqn.resize(vn);

  string temp;
  size_t trash, index_num;
  size_t nnz = 0;
  for(size_t i = 0; i < vn; ++i){
      ifs >> temp >> trash >> index_num;
      nnz += index_num;
      eqn[i].resize(index_num);
      for(size_t j = 0; j < index_num; ++j){
          ifs >> eqn[i][j].first;
        }
      for(size_t j = 0; j < index_num; ++j){
          ifs >> eqn[i][j].second;
        }
    }

  node_mapping.resize(vn, ivn, nnz);
  for(size_t i = 0; i < vn; ++i){
      node_mapping.ptr()[i+1] = node_mapping.ptr()[i] + eqn[i].size();
      for(size_t j = 0; j < eqn[i].size(); ++j){
          node_mapping.val()[node_mapping.ptr()[i] + j] = eqn[i][j].second;
          node_mapping.idx()[node_mapping.ptr()[i] + j] = eqn[i][j].first;
        }
    }

  return 0;
}

void descriptor_vol::reset(
    const zjucad::matrix::matrix<size_t> &uncut_tet,
    const zjucad::matrix::matrix<double> &uncut_node,
    const zjucad::matrix::matrix<size_t> &cut_tet,
    const zjucad::matrix::matrix<double> &cut_node,
    boost::property_tree::ptree &pt)
{
  bi_.uncut_tet = uncut_tet;
  bi_.uncut_node = uncut_node;

  init(cut_tet,cut_node,pt);
  set_objective(cut_tet, cut_node,pt);
  set_constraint(cut_tet, cut_node, pt);

}

int descriptor_vol::update_Z_q(const std::pair<size_t, double> &fix_integer)
{
  // warning fix integer should use orig variables
  // in this function, I should update the equations based on this fix_integer
  // all real_variables v are reindexed as v+orig_variable_number

  const hj::sparse::csc<double,int32_t> & ZT = bi_.NM.ZT;
  const size_t orig_node_number = bi_.NM.ZT.size(2);

  jtf::algorithm::equation<double> one_eqn;
  assert(fix_integer.first < ZT.size(2));
  for(size_t i = ZT.ptr()[fix_integer.first]; i != ZT.ptr()[fix_integer.first+1]; ++i){
      one_eqn.add_expression(jtf::algorithm::make_expression(
                               ZT.idx()[i] + orig_node_number, ZT.val()[i]));
    }
  one_eqn.value() = fix_integer.second - bi_.NM.q[fix_integer.first];
  if(one_eqn.e_vec_.empty()){
      if(bi_.NM.is_fixed[fix_integer.first]){
          if(fabs(bi_.NM.q[fix_integer.first] - fix_integer.second) > 1e-6){
              cerr << "# [error] this integer has alread beed fixed to " << bi_.NM.q[fix_integer.first]
                   << ", however it's required to be " << fix_integer.second << " again." << endl;
              cerr << "# [error] I will keep original values, and ignore the new one." << endl;
            }
          return 1;
        }else{
          cerr << "# [error] strange, this variable do not associate with any real variables, however, it's value is not fixed." << endl;
          return __LINE__;
        }
    }

  vector<double> node(bi_.NM.ZT.size(1) + bi_.NM.ZT.size(2), 0);
  boost::dynamic_bitset<> node_flag(bi_.NM.ZT.size(1) + bi_.NM.ZT.size(2));
  jtf::algorithm::gauss_eliminator<double> ge(node, node_flag);
  ge.add_equation(one_eqn);

  for(size_t i = 0; i < bi_.NM.ZT.ptr().size()-1; ++i){ // each orig_variable stands for a equation
      jtf::algorithm::equation<double> each_eqn;

      for(size_t j = bi_.NM.ZT.ptr()[i]; j != bi_.NM.ZT.ptr()[i+1]; ++j){
          each_eqn.add_expression(jtf::algorithm::make_expression(
                                    bi_.NM.ZT.idx()[j] + orig_node_number,
                                    -1* bi_.NM.ZT.val()[j]));
        }
      if(each_eqn.e_vec_.empty()) continue; // means original variable vi do not associate with any real variables
      each_eqn.add_expression(jtf::algorithm::make_expression(
                                i,1.0));
      each_eqn.value() = bi_.NM.q[i];
      ge.add_equation(each_eqn);
    }

  // after elimination, reassemble them into bi_.NM

  vector<jtf::algorithm::equation<double> > all_eqn(orig_node_number);
  set<size_t> real_variables;
  size_t nnz = 0;
  for(const auto & each_eqn : ge){
      const size_t prime_idx = each_eqn.get_prime_idx();
      if(prime_idx < orig_node_number){
          jtf::algorithm::equation<double> reduced_eqn; // move prime idx to left side of eqn, and move
          reduced_eqn.value() = each_eqn.value();
          for(const auto & one_exp : each_eqn){
              if(one_exp.index == prime_idx) continue; // ignore prime idx
              assert(one_exp.index+1 > orig_node_number);
              reduced_eqn.add_expression(jtf::algorithm::make_expression(
                                           one_exp.index-orig_node_number, -1 * one_exp.coefficient));
              real_variables.insert(one_exp.index-orig_node_number);
              ++nnz;
            }
          all_eqn[prime_idx] = reduced_eqn;
        }
    }

  map<size_t,size_t> prev_real_variable2new_idx;
  size_t new_idx = 0;
  for(const auto & idx : real_variables){
      prev_real_variable2new_idx[idx] = new_idx++;
    }
  bi_.NM.ZT.resize(real_variables.size(), orig_node_number, nnz);
  for(size_t i = 0; i < all_eqn.size(); ++i){

      const auto & each_eqn = all_eqn[i];
      bi_.NM.ZT.ptr()[i+1] = bi_.NM.ZT.ptr()[i] + each_eqn.e_vec_.size();
      size_t idx = bi_.NM.ZT.ptr()[i];
      for(const auto & one_exp : each_eqn){
          bi_.NM.ZT.idx()[idx] = prev_real_variable2new_idx[one_exp.index];
          bi_.NM.ZT.val()[idx] = one_exp.coefficient;
          ++idx;
        }
      if(each_eqn.e_vec_.empty() && bi_.NM.is_fixed[i]) continue; // if already fixed, ignore.
      if(each_eqn.e_vec_.empty()) bi_.NM.is_fixed[i] = 1; // fix it
      bi_.NM.q[i] = each_eqn.value();
    }
  assert(max(bi_.NM.ZT.idx()) < bi_.NM.ZT.size(1));
  return 0;
}

void descriptor_vol::recover_gap_node(
    zjucad::matrix::matrix<double> & node) const
{
  assert(node.size() % 3 == 0);
  itr_matrix<double*> node_mat(3, node.size()/3, &node[0]);
  assert(has_node_mapping());
  const hj::sparse::csc<double,int32_t> & NMT = get_node_mapping().ZT;
  assert(node.size() == NMT.size(2));
  assert(node.size() == (bi_.cut_tet2tet.size() + bi_.cut_patches.size()) * 3);
  matrix<double> gap(3,1);
  vector<matrix<double> > each_gaps;
  for(size_t pi = 0 ; pi < bi_.cut_patches.size(); ++pi){
      //      {
      //        stringstream ss;
      //        ss << "cut_patch_face_" << pi << ".vtk";
      //        ofstream ofs(ss.str().c_str());
      //        vector<size_t> cut_patch;
      //        vector<size_t> cut_patch_faces;
      //        vector<size_t> cut_patch_rot;
      //        for(const auto & face_pair : bi_.cut_patches[pi]){
      //            cut_patch_faces.insert(cut_patch_faces.end(),
      //                                   bi_.gap_faces(colon(),face_pair.first).begin(),
      //                                   bi_.gap_faces(colon(),face_pair.first).end());
      //            cut_patch_faces.insert(cut_patch_faces.end(),
      //                                   bi_.gap_faces(colon(),face_pair.second).begin(),
      //                                   bi_.gap_faces(colon(),face_pair.second).end());
      //            cut_patch.push_back(0);
      //            cut_patch.push_back(1);
      //            cut_patch_rot.push_back(bi_.rot_type[pi]);
      //            cut_patch_rot.push_back(bi_.rot_type[pi]);
      //          }
      //        tri2vtk(ofs, &node[0], node.size(2), &cut_patch_faces[0], cut_patch_faces.size()/3);
      //        cell_data(ofs, &cut_patch[0], cut_patch.size(), "patch");
      //        vtk_data(ofs, &cut_patch_rot[0], cut_patch_rot.size(), "rot");
      //      }
      gap *= 0;
      const vector<pair<size_t,size_t> > & one_patch = bi_.cut_patches[pi];
      const size_t & type = bi_.rot_type[pi];
      const zjucad::matrix::matrix<double> rot_mat = trans(type_transition2(type));
      set<pair<size_t,size_t> > point_mapping;
      for(size_t fi = 0; fi < one_patch.size(); ++fi){
          const pair<size_t,size_t> & face_pair = one_patch[fi];
          for(size_t pj = 0; pj < 3; ++pj){
              point_mapping.insert(make_pair(bi_.gap_faces(pj, face_pair.first),
                                             bi_.gap_faces(pj, face_pair.second)));
            }
        }

      ////////////////
      //      each_gaps.clear();
      //      stringstream ss;
      //      ss  << "gaps_error_" << pi;
      //      ofstream ofs(ss.str().c_str());
      //////////////////

      for(const auto & one_point_pair : point_mapping){
          gap += -1*rot_mat * node_mat(colon(), one_point_pair.first)
              + node_mat(colon(), one_point_pair.second);
          //          {
          //            matrix<double> one_gap = -1*rot_mat * node_mat(colon(), one_point_pair.first)
          //                + node_mat(colon(), one_point_pair.second);
          //            each_gaps.push_back(one_gap);
          //            ofs << one_point_pair.first << " to " << one_point_pair.second << endl;
          //            ofs << one_gap << endl;
          //          }
        }
      gap /= (point_mapping.size());
      //      ofs << " avg gap."<< endl;
      //      ofs << gap << endl;
      //      {
      //        double v = 0;
      //        for(size_t i = 0; i < each_gaps.size(); ++i){
      //            v += norm(each_gaps[i] - gap);
      //          }
      //        v /= each_gaps.size();
      //        cerr << "# [info] patch " << pi << " gap error " << v << endl;
      //      }
      node_mat(colon(), bi_.cut_tet2tet.size() + pi) = gap;
      for(size_t di = 0; di < 3; ++di){
          if(bi_.NM.is_fixed[3*(bi_.cut_tet2tet.size() + pi)+di])
            node_mat(di, bi_.cut_tet2tet.size() + pi) = bi_.NM.q[3*(bi_.cut_tet2tet.size() + pi)+di];
        }
    }
}

void descriptor_vol::basic_infor::get_independent_node(
    const zjucad::matrix::matrix<size_t> & node_mapping,
    const zjucad::matrix::matrix<double> & node,
    zjucad::matrix::matrix<double> & independent_node)const
{
  set<size_t> independent_variable(node_mapping.begin(), node_mapping.end());
  assert(*independent_variable.rbegin() == independent_variable.size() -1);
  independent_node.resize(independent_variable.size());
  vector<size_t> v2orig_v(independent_variable.size());
  for(size_t i = 0; i < node_mapping.size(); ++i)
    v2orig_v[node_mapping[i]] = i;

  for(size_t i = 0; i < v2orig_v.size(); ++i){
      independent_node[i] = node[v2orig_v[i]];
    }
}
