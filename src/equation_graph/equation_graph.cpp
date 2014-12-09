#include <numeric>
#include <tuple>
#include <vector>

#include <jtflib/util/vertex_connection.h>
#include <jtflib/util/container_operation.h>
#include <jtflib/mesh/mesh.h>
#include <jtflib/math/math.h>


#include "cycle_detection.h"
#include "equation_graph.h"
#include "util.h"
#include "../common/timer.h"
#include "../numeric/util.h"
#include "../common/transition_type.h"
#include "../common/extract_loop_from_undirected_edges.h"
#include "../tetmesh/util.h"
#include "../common/util.h"
#include "../common/vtk.h"
#include "../common/visualize_tool.h"
#include "../common/transition_type.h"

using namespace std;
using namespace zjucad::matrix;

#define debug

size_t get_edge_type_with_part_face_type_map(
    const std::vector<size_t> & tet_seq,
    const boost::unordered_map<std::pair<size_t,size_t>,size_t> & face_type )
{
  typedef  boost::unordered_map<pair<size_t,size_t>,size_t>::const_iterator mpscit;
  matrixd rot = eye<double>(3);
  for(size_t t = 0; t < tet_seq.size()-1; ++t){
      mpscit it = face_type.find(make_pair(tet_seq[t], tet_seq[t+1]));
      if(it == face_type.end()) continue;
      rot = temp(rot * type_transition2(it->second));
    }
  return type_transition1(rot);
}

int transition_elimination::get_trivial_cut_face_pair(
    const boost::unordered_map<std::pair<size_t,size_t>,size_t> & inner_type,
    const zjucad::matrix::matrix<size_t> & cut_tet,
    const zjucad::matrix::matrix<size_t> & cut_tet2tet,
    vector<matrix<size_t> > & trivial_cut_face_pairs) const
{
  trivial_cut_face_pairs.clear();
  map<size_t,vector<size_t> >  point_mapping;
  matrix<size_t> one_pair_trivial_cut_face(3,2);
  for(const auto & tet_pair_idx : bi_.tet_pair2g_idx){
      const pair<size_t,size_t> & tet_pair = tet_pair_idx.first;

      if(gnode_flag_[3 * tet_pair_idx.second + 0] &&
         gnode_flag_[3 * tet_pair_idx.second + 1] &&
         gnode_flag_[3 * tet_pair_idx.second + 2]){
          assert(fabs(gnodes_[3 * tet_pair_idx.second + 0]) < 1e-6);
          assert(fabs(gnodes_[3 * tet_pair_idx.second + 1]) < 1e-6);
          assert(fabs(gnodes_[3 * tet_pair_idx.second + 2]) < 1e-6);
          const auto it = inner_type.find(tet_pair);
          if(it != inner_type.end() && it->second != TRIVIAL_TYPE) continue;// not trivial cut face
          point_mapping.clear();
          for(size_t pi = 0; pi < cut_tet.size(1); ++pi){
              point_mapping[cut_tet2tet[cut_tet(pi, tet_pair.first)]].push_back(
                    cut_tet(pi, tet_pair.first));
            }
          for(size_t pi = 0; pi < cut_tet.size(1); ++pi){
              point_mapping[cut_tet2tet[cut_tet(pi, tet_pair.second)]].push_back(
                    cut_tet(pi, tet_pair.second));
            }

          size_t j = 0;
          for(const auto & one_point : point_mapping){
              if(one_point.second.size() == 2){
                  one_pair_trivial_cut_face(j,0) = one_point.second.front();
                  one_pair_trivial_cut_face(j,1) = one_point.second.back();
                  ++j;
                }
            }
          if(j != 3){
              cerr << "# [error] strange, cut face pair error." << endl;
              return __LINE__;
            }
        }
      trivial_cut_face_pairs.push_back(one_pair_trivial_cut_face);
    }
}
transition_elimination * transition_elimination::create(
    const zjucad::matrix::matrix<size_t> & tetmesh,
    const zjucad::matrix::matrix<size_t> & cut_mesh,
    const zjucad::matrix::matrix<size_t> & cut_tet2tet,
    const zjucad::matrix::matrix<double> & tet_node,
    const boost::unordered_map<std::pair<size_t,size_t>,size_t> & inner_face_type,
    const boost::unordered_map<size_t,size_t> & surface_type,
    const bool is_restricted_type)
{
  unique_ptr<transition_elimination> te(new transition_elimination);
  if(te->init(tetmesh,cut_mesh, cut_tet2tet, tet_node, inner_face_type, surface_type,
              is_restricted_type))
    return 0;
  return te.release();
}

transition_elimination * transition_elimination::create_with_gap(
    const zjucad::matrix::matrix<size_t> & tetmesh,
    const zjucad::matrix::matrix<size_t> & cut_mesh,
    const zjucad::matrix::matrix<size_t> & cut_tet2tet,
    const boost::unordered_map<std::pair<size_t,size_t>,size_t> & inner_face_type,
    const boost::unordered_map<size_t,size_t> & surface_type,
    const bool is_restricted_type,
    const zjucad::matrix::matrix<double> &node)
{
  unique_ptr<transition_elimination> te(new transition_elimination);
  if(te->init_with_gap(tetmesh,cut_mesh, cut_tet2tet, inner_face_type, surface_type,
                       is_restricted_type, node))
    return 0;
  return te.release();
}

int transition_elimination::init(
    const matrix<size_t> &tetmesh,
    const matrix<size_t> &cut_mesh,
    const matrix<size_t> &cut_tet2tet,
    const matrix<double> & tet_node,
    const boost::unordered_map<pair<size_t,size_t>,size_t > &inner_face_type,
    const boost::unordered_map<size_t,size_t> &surface_type,
    const bool is_restricted_type)
{
  vector<double> nodes;
  boost::dynamic_bitset<> node_flag;
  construct_basic_info(tetmesh, cut_mesh, cut_tet2tet, tet_node,
                       nodes, node_flag, gnodes_, gnode_flag_);
  get_gaps_stupid(tetmesh, cut_mesh, cut_tet2tet, inner_face_type, tet_node);

  get_gaps_elimination(tetmesh, cut_mesh, cut_tet2tet, inner_face_type, tet_node);
  add_inner_transition(tetmesh, cut_mesh, cut_tet2tet, inner_face_type);

  if(!surface_type.empty())
    add_surface_transition(tetmesh, cut_mesh, cut_tet2tet, surface_type,is_restricted_type);

  collect_group_variants_from_equations();
  return 0;
}

int transition_elimination::construct_basic_info_with_gaps(
    const zjucad::matrix::matrix<size_t> & tetmesh,
    const zjucad::matrix::matrix<size_t> & cut_mesh,
    const zjucad::matrix::matrix<size_t> & cut_tet2tet,
    const boost::unordered_map<std::pair<size_t,size_t>,size_t> & inner_face_type,
    const zjucad::matrix::matrix<double> & node)
{
  bi_.fa_cut.reset(jtf::mesh::face2tet_adjacent::create(cut_mesh));
  bi_.fa.reset(jtf::mesh::face2tet_adjacent::create(tetmesh));
  if(!bi_.fa_cut.get() || !bi_.fa.get()){
      throw std::logic_error("can not build face2tet_adjacent.");
    }

  jtf::mesh::get_outside_face(*bi_.fa_cut, bi_.outside_face_cut);

  vector<size_t> face(3), sort_cut_face(3);
  for(size_t fi = 0; fi < bi_.outside_face_cut.size(2); ++fi){
      for(size_t pi = 0; pi < bi_.outside_face_cut.size(1); ++pi)
        face[pi] = cut_tet2tet[bi_.outside_face_cut(pi,fi)];
      const size_t face_idx_orig = bi_.fa->get_face_idx(&face[0]);
      if(face_idx_orig == -1){
          cerr << "# [error] can not find face idx." << endl;
          return __LINE__;
        }
      if(bi_.fa->is_outside_face(bi_.fa->face2tet_[face_idx_orig])) continue;
      sort(face.begin(), face.end());
      for(size_t pi = 0; pi < bi_.outside_face_cut.size(1); ++pi)
        for(size_t pj = 0; pj < bi_.outside_face_cut.size(1); ++pj){
            if(face[pi] == cut_tet2tet[bi_.outside_face_cut(pj, fi)]){
                sort_cut_face[pi] = bi_.outside_face_cut(pj, fi);
                continue;
              }
          }
      bi_.orig_face2_cut_faces[face].push_back(sort_cut_face);
    }

  bi_.cut_face_pairs.resize(3, bi_.orig_face2_cut_faces.size()*2);

  size_t idx = 0;
  for(const auto & orig_face2cut_face: bi_.orig_face2_cut_faces){
      const vector<basic_info::one_face> & cut_face = orig_face2cut_face.second;
      pair<size_t,size_t> tet_pair =
          bi_.get_tet_pair_from_cut_face(cut_face.front(),  cut_face.back());
      bi_.tet_pair2g_idx[tet_pair] = idx;
      bi_.orig_face2g_idx[orig_face2cut_face.first] = idx;

      const vector<basic_info::one_face> & two_faces = orig_face2cut_face.second;
      std::copy(two_faces.front().begin(), two_faces.front().end(),
                bi_.cut_face_pairs(colon(),2 * idx).begin());
      std::copy(two_faces.back().begin(), two_faces.back().end(),
                bi_.cut_face_pairs(colon(),2 * idx + 1).begin());
      ++idx;
    }

  bi_.ea_cut_face_pair.reset(jtf::mesh::edge2cell_adjacent::create(bi_.cut_face_pairs));
  if(!bi_.ea_cut_face_pair.get()){
      throw std::logic_error("can not build edge2tet_adjacent.");
    }

  bi_.ortae.add_tets(tetmesh, *(bi_.fa.get()));
  bi_.ortae.sort_into_loop(tetmesh, node, 1);

  return 0;
}

size_t transition_elimination::get_face_rot_type(
    const zjucad::matrix::matrix<size_t> & one_face_0,
    const zjucad::matrix::matrix<size_t> & one_face_1,
    const jtf::mesh::face2tet_adjacent &fa_cut,
    const boost::unordered_map<pair<size_t,size_t>,size_t> & inner_type)
{
  const auto it_0 = fa_cut.get_face_idx(&one_face_0[0]);
  const auto it_1 = fa_cut.get_face_idx(&one_face_1[0]);
  assert(it_0 != -1 && it_1 != -1);
  const pair<size_t,size_t> & tet_pair_0 = fa_cut.face2tet_[it_0];
  const pair<size_t,size_t> & tet_pair_1 = fa_cut.face2tet_[it_1];
  pair<size_t,size_t> tet_pair(
        tet_pair_0.first==-1?tet_pair_0.second:tet_pair_0.first,
        tet_pair_1.first==-1?tet_pair_1.second:tet_pair_1.first);

  const auto it_type = inner_type.find(tet_pair);
  if(it_type == inner_type.end())
    return TRIVIAL_TYPE;
  else
    return it_type->second;
}

void transition_elimination::reorder_nodes(
    const zjucad::matrix::matrix<size_t> & cut_tet,
    matrix<size_t> & M_new2old, // recorded mapping from new index to new one
    matrix<size_t> & M_old2new,
    size_t &P_out_end_idx)
{
  assert(bi_.outside_face_cut.size());
  set<size_t> outside_face_point_cut(bi_.outside_face_cut.begin(), bi_.outside_face_cut.end());
  set<size_t> all_points(cut_tet.begin(), cut_tet.end());

  vector<size_t> other_points_cut(all_points.size() - outside_face_point_cut.size());

  set_difference(all_points.begin(), all_points.end(), outside_face_point_cut.begin(),
                 outside_face_point_cut.end(), other_points_cut.begin());

  const size_t total_points_num = outside_face_point_cut.size() + other_points_cut.size();
  assert(max(cut_tet) + 1 == total_points_num);
  M_new2old.resize(3 *(total_points_num + bi_.cut_face_patches_.size()),1);
  size_t pi = 0;
  // fill with outside_face_points
  for(const auto & point : outside_face_point_cut){
      for(size_t di = 0; di < 3; ++di)
        M_new2old[3 * pi + di] = 3 * point + di;
      ++pi;
    }

  // fill with gaps
  for(size_t gi = 0; gi < bi_.cut_face_patches_.size(); ++gi){
      for(size_t di = 0; di < 3; ++di){
          M_new2old[3 * pi +di] = 3 * (total_points_num + gi) + di;
        }
      ++pi;
    }
  P_out_end_idx =  pi;

  // fill with inner points
  for(const auto & point : other_points_cut){
      for(size_t di = 0; di < 3; ++di){
          M_new2old[3 * pi + di] = 3 * point + di;
        }
      ++pi;
    }

  M_old2new.resize(M_new2old.size(1), M_new2old.size(2));
  for(size_t i = 0; i < M_new2old.size(); ++i){
      M_old2new[M_new2old[i]] = i;
    }

  for(const auto & idx : outside_face_point_cut){
      for(size_t di = 0; di < 3; ++di){
          if(M_old2new[3 * idx + di] >= 3 * P_out_end_idx)
            throw std::logic_error("invalide node mapping.");
        }
    }
  for(size_t gi = 3 * total_points_num;
      gi < 3 * (total_points_num + bi_.cut_face_patches_.size()); ++gi){
      if(M_old2new[gi] >= 3 * P_out_end_idx)
        throw std::logic_error("invalid node mapping.");

    }
}

void transition_elimination::assemble_inner_transition_and_surface_type(
    const zjucad::matrix::matrix<size_t> & tetmesh,
    const zjucad::matrix::matrix<size_t> & cut_mesh,
    const zjucad::matrix::matrix<size_t> & cut_tet2tet,
    const boost::unordered_map<pair<size_t,size_t>,size_t> & inner_type,
    const boost::unordered_map<size_t,size_t> & surface_type,
    const bool is_restricted_type,
    const size_t original_node_number_with_gaps)
{
  add_surface_transition(tetmesh, cut_mesh, cut_tet2tet,
                         surface_type, is_restricted_type);

  add_inner_transition_edge2(tetmesh, cut_mesh, cut_tet2tet, inner_type);

  add_inner_transition_gap2(tetmesh, cut_mesh, cut_tet2tet,original_node_number_with_gaps);

  add_inner_transition_gap3(tetmesh, cut_mesh, cut_tet2tet, inner_type, original_node_number_with_gaps);
}

void transition_elimination::get_node_comp(
    hj::sparse::csc<double,int32_t> & MT) const
{
  MT = node_comp_mapping;
}

void transition_elimination::get_zero_idx(std::set<size_t> & zero_set) const
{
  zero_set.clear();
  for(size_t gi = 0; gi < gnode_flag_.size(); ++gi){
      if(gnode_flag_[gi] == true){
          zero_set.insert(gi);
        }
    }
}

void transition_elimination::get_ordered_cut_face_patches(
    zjucad::matrix::matrix<size_t> & cut_faces,
    std::vector<std::vector<std::pair<size_t,size_t> > > & patches,
    std::vector<size_t> & rot_type)const
{
  cut_faces = bi_.cut_face_pairs;
  patches = bi_.cut_face_patches_;
  rot_type = bi_.cut_face_patch_rot_type_;
}

void transition_elimination::compress_equation(
    const zjucad::matrix::matrix<size_t> &M_old2new,
    const size_t P_out_end_index,
    vector<jtf::algorithm::equation<double> > & all_eqn,
    size_t & independent_variable_number)
{
  // ge should be already gauss eliminated
  assert(ge.get());

  // equation convertion from old index to new index
  vector<double> node(3*P_out_end_index);
  boost::dynamic_bitset<> node_flag(3*P_out_end_index);

  unique_ptr<jtf::algorithm::gauss_eliminator<double> > ge_new(
        new jtf::algorithm::gauss_eliminator<double>(node, node_flag));
  if(!ge_new.get()){
      throw std::logic_error("can not build gauss_elimination");
    }

  for(jtf::algorithm::gauss_eliminator<double>::equation_ptr one_eqn_it = ge->begin();
      one_eqn_it != ge->end(); ++one_eqn_it){
      if(one_eqn_it->e_vec_.size() == 1) continue;
      jtf::algorithm::equation<double> new_eqn = *one_eqn_it;

      for(auto & exp : new_eqn){
          assert(M_old2new[exp.index] < 3 * P_out_end_index);
          exp.index = M_old2new[exp.index];
          assert(exp.index < 3*P_out_end_index);
        }
      ge_new->add_equation(new_eqn);
    }


  vector<bool> is_independent_variables(3*P_out_end_index, true);
  vector<jtf::algorithm::equation<double> > v2eqn(3 * P_out_end_index);
  // x0 - x2 + x3 = 0; is rewrited as x0 = x2-x3
  //map<size_t, jtf::algorithm::equation<double> > dependent_v_2_eqn;

  for(const auto & one_eqn : *ge_new){
      //const jtf::algorithm::equation<double> & one_eqn = *one_eqn_it;
      assert(one_eqn.e_vec_.size() != 1);
      size_t pidx = one_eqn.get_prime_idx();

      jtf::algorithm::equation<double> reduced_eqn;
      jtf::algorithm::equation<double>::eq_const_iterator it = one_eqn.begin();
      ++it;
      for(;it != one_eqn.end(); ++it){
          const auto & one_exp = *it;
          reduced_eqn.add_expression(
                jtf::algorithm::make_expression(
                  one_exp.index,-1 *one_exp.coefficient));
        }
      v2eqn[pidx] = reduced_eqn;
      is_independent_variables[pidx] = false;
    }
  // after this step, I find all dependent variables

  // next step is to reorder them and keep each equation on up of variable list
  vector<size_t> independent_varible;
  vector<size_t> dependent_vaiable;

  for(size_t i = 0; i < is_independent_variables.size(); ++i){
      if(is_independent_variables[i])
        independent_varible.push_back(i);
      else
        dependent_vaiable.push_back(i);
    }
  matrix<size_t> variable_mapping_new2old(3*P_out_end_index,1),
      variable_mapping_old2new(3*P_out_end_index,1);
  size_t i = 0;
  for(; i < dependent_vaiable.size(); ++i){
      variable_mapping_new2old[i] = dependent_vaiable[i];
    }
  for(; i < 3*P_out_end_index; ++i){
      variable_mapping_new2old[i] = independent_varible[i - dependent_vaiable.size()];
    }
  for(size_t i = 0; i < variable_mapping_new2old.size(); ++i){
      variable_mapping_old2new[variable_mapping_new2old[i]] = i;
    }

  vector<jtf::algorithm::equation<double> > eqns;
  eqns.resize(M_old2new.size());


  // assemble all node mapping
  size_t eqi = 0;
  size_t offset = dependent_vaiable.size();
  for(; eqi < dependent_vaiable.size(); ++eqi){
      jtf::algorithm::equation<double> & one_eqn = v2eqn[variable_mapping_new2old[eqi]];
      for(auto & one_exp : one_eqn){
          one_exp.index = variable_mapping_old2new[one_exp.index]-offset;
        }
      eqns[variable_mapping_new2old[eqi]] = one_eqn;

    }
  for(; eqi < 3*P_out_end_index; ++eqi){
      jtf::algorithm::equation<double> one_eqn;
      one_eqn.add_expression(
            jtf::algorithm::make_expression(
              variable_mapping_old2new[independent_varible[eqi - dependent_vaiable.size()]]-offset,1.0));
      eqns[variable_mapping_new2old[eqi]] = one_eqn;
    }
  for(; eqi < M_old2new.size(); ++eqi){
      jtf::algorithm::equation<double> one_eqn;
      one_eqn.add_expression(
            jtf::algorithm::make_expression(eqi-offset,1.0));
      eqns[eqi] = one_eqn;
    }

  //reorder node mapping according to M_old2new
  all_eqn.resize(eqns.size());
  size_t nnz = 0;
  for(size_t i = 0; i < M_old2new.size(); ++i){
      all_eqn[i] = eqns[M_old2new[i]];
      nnz += all_eqn[i].e_vec_.size();
    }

  independent_variable_number = independent_varible.size() + M_old2new.size() - 3*P_out_end_index;

  // X = Z*u
  // eliminate col of matrix

  {//convert equations to hj::space::csc
    // get Z^T
    hj::sparse::csc<double,int32_t> ZT(independent_variable_number, all_eqn.size(), nnz);

    for(size_t eqi = 0; eqi < all_eqn.size(); ++eqi){
        ZT.ptr()[eqi+1] = ZT.ptr()[eqi] + all_eqn[eqi].e_vec_.size();
        size_t idx = ZT.ptr()[eqi];
        for(const auto & one_exp : all_eqn[eqi]){
            ZT.idx()[idx] = one_exp.index;
            ZT.val()[idx] = one_exp.coefficient;
            ++idx;
          }
      }
    // I want to get each row of ZT, which means get each col of Z
    hj::sparse::csc<double,int32_t> Z;
    hj::sparse::trans(ZT,Z);

    // eliminate each col of Z
    vector<double> gnode_Z(all_eqn.size());
    boost::dynamic_bitset<> gnode_Z_flag(all_eqn.size());

    jtf::algorithm::gauss_eliminator<double> ge_Z(gnode_Z, gnode_Z_flag);
    for(size_t i = 0; i < Z.ptr().size()-1; ++i){
        jtf::algorithm::equation<double> one_eqn;
        for(size_t idx = Z.ptr()[i]; idx < Z.ptr()[i+1]; ++idx){
            one_eqn.add_expression(jtf::algorithm::make_expression(
                                     Z.idx()[idx], Z.val()[idx]));
          }
        ge_Z.add_equation(one_eqn);
      }

    // get finial ge to ZT_compressed
    size_t nnz_finial = 0;
    for(const auto & one_eqn : ge_Z){
        nnz_finial += one_eqn.e_vec_.size();
      }
    hj::sparse::csc<double,int32_t> ZT_compressed(all_eqn.size(), ge_Z.get_eqn_number(), nnz_finial);
    size_t eqi = 0;
    for(const auto & one_eqn : ge_Z){
        ZT_compressed.ptr()[eqi+1] = ZT_compressed.ptr()[eqi] + one_eqn.e_vec_.size();
        size_t idx = ZT_compressed.ptr()[eqi];
        for(const auto & one_exp : one_eqn){
            ZT_compressed.idx()[idx] = one_exp.index;
            ZT_compressed.val()[idx] = one_exp.coefficient;
            ++idx;
          }
        ++eqi;
      }

    // get finial Z_compressed
    hj::sparse::csc<double,int32_t> Z_compressed;
    hj::sparse::trans(ZT_compressed, Z_compressed);

    for(size_t i = 0; i < Z_compressed.ptr().size()-1; ++i){
        jtf::algorithm::equation<double> one_eqn;
        for(size_t idx = Z_compressed.ptr()[i]; idx != Z_compressed.ptr()[i+1]; ++idx){
            one_eqn.add_expression(jtf::algorithm::make_expression(Z_compressed.idx()[idx],
                                                                   Z_compressed.val()[idx]));
          }
        all_eqn[i] = one_eqn;
      }
    independent_variable_number = Z_compressed.size(1);
  }
}

void transition_elimination::build_transitions(
    const zjucad::matrix::matrix<size_t> &uncut_mesh,
    const zjucad::matrix::matrix<size_t> &cut_mesh,
    const zjucad::matrix::matrix<size_t> &cut_tet2tet,
    const boost::unordered_map<std::pair<size_t,size_t>,size_t> &inner_face_type,
    const boost::unordered_map<size_t,size_t> & surface_type,
    bool is_restricted_type)
{//in this step, I will finish a null space basis extraction
  // denote (P:G)^T as a colon vector which stores all node variables and gap variables
  // First, (P:G)^T can be reorder to gather all node on cut faces and gap node ahead,
  // by (P:G)^T = M*(P_out:P_inner)^T, note that gaps are inside P_out
  // Second, We can apply nullspace basis extraction on P_out by P_out^T = Z*u^T
  // u is independ variables, thus
  // | P_out  |    | Z  | |  u    |
  // |        | =  |    | |       |
  // | P_inner|    |   I| |P_inner|
  // Third, combine them togather:
  // (P:G)^T = M * (Z, I) * (u: P_inner)^T

  // M_inv recorded new idx of each original variable
  matrix<size_t> M_new2old;
  matrix<size_t> M_old2new;
  size_t P_out_end_index  = -1; // end of P_out: [0,P_out_end_index)
  reorder_nodes(cut_mesh, M_new2old, M_old2new, P_out_end_index);

  {//
    node2group_.clear();
    node_group_.clear();
    node2group_.resize(M_old2new.size());
    node_group_.resize(M_old2new.size());
    for(size_t i = 0; i < M_old2new.size(); ++i){
        node2group_[i] = i;
        node_group_[i].insert(i);
      }
  }

  gnodes_.resize(M_old2new.size());
  gnode_flag_.resize(M_old2new.size());
  ge.reset(
        new jtf::algorithm::gauss_eliminator<double>(gnodes_, gnode_flag_));

  assert(M_old2new.size() %3 == 0);
  assemble_inner_transition_and_surface_type(
        uncut_mesh, cut_mesh, cut_tet2tet, inner_face_type,
        surface_type, is_restricted_type, M_old2new.size()/3);

  {
    ofstream ofs("debug_eqn");
    for(const auto & one_eqn : *ge){
        ofs << "eqi" << endl;
        for(const auto & one_exp : one_eqn){
            ofs << one_exp.index << " ";
          }
        ofs << endl;
        for(const auto & one_exp : one_eqn){
            ofs << one_exp.coefficient << " ";
          }
        ofs << endl;
      }
  }

  vector<jtf::algorithm::equation<double> > eqn;

  size_t independent_variable_num;
  compress_equation(M_old2new,P_out_end_index, eqn,independent_variable_num);

  cerr << "# [info] original variables " << M_old2new.size()
       << " independent variables " << independent_variable_num << endl;

  size_t nnz = 0;
  vector<pair<size_t,double> > one_eqn;
  for(size_t i = 0; i < eqn.size(); ++i){
      nnz += eqn[i].e_vec_.size();
    }
  node_comp_mapping.resize(independent_variable_num, M_old2new.size(), nnz);

  for(size_t i = 0; i < eqn.size(); ++i){
      const jtf::algorithm::equation<double> & one_eqn = eqn[i];
      node_comp_mapping.ptr()[i+1] = node_comp_mapping.ptr()[i] + one_eqn.e_vec_.size();
      size_t vi = 0;
      for(const auto & exp : one_eqn){
          node_comp_mapping.idx()[node_comp_mapping.ptr()[i] + vi] = exp.index;
          node_comp_mapping.val()[node_comp_mapping.ptr()[i] + vi] = exp.coefficient;
          ++vi;
        }
    }

  {
    ofstream ofs("node_compression"), ofs_m("node_compress.mat");
    ofs << M_old2new.size() << " " << independent_variable_num << endl;
    for(size_t eqi = 0; eqi < eqn.size(); ++eqi){
        ofs << "eqi " << eqi << " " << eqn[eqi].e_vec_.size() << endl;
        const jtf::algorithm::equation<double> & one_eqn = eqn[eqi];
        for(const auto & one_exp : one_eqn){
            ofs << one_exp.index << " ";
            ofs_m << eqi << " " << one_exp.index << " " << one_exp.coefficient << endl;
          }
        ofs << endl;
        for(const auto & one_exp : one_eqn){
            ofs << one_exp.coefficient << " ";
          }
        ofs << endl;
      }
  }
  {
    ofstream ofs("zero_gap");
    for(size_t gi = 0; gi < gnode_flag_.size(); ++gi){
        if(gnode_flag_[gi] == true){
            if(fabs(gnodes_[gi]) > 1e-6)
              throw std::logic_error("Found one known node not to be zero after elimination.");
            ofs << gi << " ";
          }
      }
  }
}

void extract_groups_from_equations(const jtf::algorithm::gauss_eliminator<double> &ge_gap,
                                   const vector<double> &gnode,
                                   const boost::dynamic_bitset<> &gnode_flag,
                                   vector<group<size_t>> & gnode_groups,
                                   vector<size_t> &gnode_to_group)
{// this function assume all equations are well eliminated, each equation contains only
  // one dependent variable
  map<vector<pair<size_t,double> > , vector<size_t> >  eqn2variable;
  vector<size_t> zero_node;
  vector<pair<size_t,double> > reduced_eqn;
  // add basic groups
  for(size_t i = 0; i < gnode_to_group.size(); ++i){
      reduced_eqn.clear();
      reduced_eqn.push_back(make_pair(i,1));
      eqn2variable[reduced_eqn].push_back(i);
    }

  for(const auto & one_eqn : ge_gap){
      if(one_eqn.e_vec_.size() == 1){
          const size_t prime_idx = one_eqn.get_prime_idx();
          assert(fabs(gnode[prime_idx]) < 1e-6);
          zero_node.push_back(prime_idx);
        }else{
          const size_t prime_idx = one_eqn.get_prime_idx();
          reduced_eqn.clear();
          for(const auto & one_exp : one_eqn){
              if(one_exp.index == prime_idx) continue;
              reduced_eqn.push_back(make_pair(one_exp.index, -1*one_exp.coefficient));
            }
          eqn2variable[reduced_eqn].push_back(prime_idx);
        }
    }

  gnode_to_group.resize(gnode.size());
  gnode_groups.resize(gnode.size());
  for(size_t i = 0; i < gnode.size(); ++i){
      gnode_to_group[i] = i;
      gnode_groups[i].insert(i);
    }

  for(size_t i = 1; i < zero_node.size(); ++i){
      const size_t &from = zero_node[i];
      if(gnode_to_group[from] == gnode_to_group[zero_node[0]]) continue;
      group<size_t> & from_group = gnode_groups[gnode_to_group[from]];
      for(const auto & idx : from_group) gnode_to_group[idx] = gnode_to_group[zero_node[0]];
      gnode_groups[gnode_to_group[zero_node[0]]].merge(from_group);
    }

  for(const auto & one_group : eqn2variable){
      const vector<size_t> & group_of_index = one_group.second;
      for(size_t i = 1; i < group_of_index.size(); ++i){
          const size_t &from = group_of_index[i];
          if(gnode_to_group[from] == gnode_to_group[group_of_index[0]]) continue;
          group<size_t> & from_group = gnode_groups[gnode_to_group[from]];
          for(const auto & idx : from_group) gnode_to_group[idx] = gnode_to_group[group_of_index[0]];
          gnode_groups[gnode_to_group[group_of_index[0]]].merge(from_group);
        }
    }
}
void transition_elimination::cluster_transition_patch_by_elimination(
    const zjucad::matrix::matrix<size_t> & tetmesh,
    const zjucad::matrix::matrix<size_t> & cut_mesh,
    const zjucad::matrix::matrix<size_t> & cut_tet2tet,
    const zjucad::matrix::matrix<double> & tet_node,
    const boost::unordered_map<std::pair<size_t,size_t>,size_t> & inner_type)
{
  // assume function construct_basic_info_with_gaps is called.
  // to avoid impact of order of each gap, I add double variables:
  // g_st, and g_ts, and add equation g_st = -Pi_st*g_ts
  vector<double> gnode_(bi_.tet_pair2g_idx.size() * 2 * 3);
  boost::dynamic_bitset<> gnode_flag(gnode_.size());
  unique_ptr<jtf::algorithm::gauss_eliminator<double> > ge_gap(
        new jtf::algorithm::gauss_eliminator<double>(gnode_, gnode_flag));
  double tt = 0;
  double tt_idx = 0;

  { // assemble gap equations
    timer tim, tm_eqn;
    tim.start();
    vector<size_t> zero_variant;

    for(const auto & one_edge2tet_loop : bi_.ortae.e2t_){
        const pair<size_t,size_t> & raw_edge = one_edge2tet_loop.first;

        const vector<size_t> & tet_loop = one_edge2tet_loop.second;
        if(!bi_.ortae.is_inner_edge(tet_loop)) continue;

        vector<matrix<double> > gnode;
        matrix<double> total_rot = eye<double>(3);
        matrix<double> one_gnode;
        size_t edge_type = get_edge_type_with_part_face_type_map(tet_loop, inner_type);

        for(size_t i = tet_loop.size() - 2; i != -1; --i){
            size_t forward_type = bi_.get_inner_face_type(
                  make_pair(tet_loop[i+1], tet_loop[(i+2)%tet_loop.size()]),
                inner_type);
            total_rot = temp(total_rot * trans(type_transition2(forward_type)));

            int rtn = bi_.get_gnode_idx(make_pair(tet_loop[i], tet_loop[i+1]),one_gnode);
            one_gnode += 1.0; // to avoid jtf::math::get_sign(0)
            if(rtn != 0 && rtn != 1) continue;
            if(rtn == 1){
                // f_t = Pi_st*f_s+g_st ==> f_s = Pi_ts*f_t-Pi_ts*g_st
                // Goal: g_st; Get: g_ts
                // -Pi_ts*g_st=g_ts ==> g_st = -Pi_st*g_ts
                size_t type_ts = bi_.get_inner_face_type(
                      make_pair(tet_loop[i],tet_loop[i+1]), inner_type);
                assert(type_ts < 24); // only 24 roation typeswe
                matrix<double> rot_mat = -1*trans(type_transition2(type_ts));
                one_gnode = temp(rot_mat * one_gnode);
              }

            one_gnode = temp(total_rot * one_gnode);
            gnode.push_back(one_gnode);
          }

        zero_variant.clear();
        if(edge_type == TRIVIAL_TYPE){
            zero_variant.push_back(0); // three variants are zeros
            zero_variant.push_back(1);
            zero_variant.push_back(2);
          }else{
            size_t axis_type = axis_to_around(edge_type);
            if(axis_type == -1){ // compound edge
                zero_variant.push_back(0);
                zero_variant.push_back(1);
              }else
              zero_variant.push_back(axis_type);
          }

        for(size_t di = 0; di < zero_variant.size(); ++di){
            jtf::algorithm::equation<double> eqn;
            for(size_t gi = 0; gi < gnode.size(); ++gi){
                eqn.add_expression(jtf::algorithm::make_expression(
                                     abs(number_rounding(gnode[gi][zero_variant[di]]))-1,
                                   jtf::math::get_sign(gnode[gi][zero_variant[di]])));
              }
            tm_eqn.start();
            ge_gap->add_equation(eqn);
            tm_eqn.finish();
            tt += tm_eqn.result();
            tt_idx += 1;
          }
      }
    tim.finish();
    cerr << "# [info] +++ gauss elimination around edge : time " << tim.result() << endl;

    cerr << "# [info] gauss eliminator eqn number " << tt_idx << endl;
    cerr << "# [info] gauss eliminator sec/eqn " << tt/tt_idx << endl;
    // add inverse order gaps
    // g_ts + Pi_ts*g_st = 0
    tim.start();
    tt= 0;
    tt_idx = 0;
    matrix<double> gap_node_pair(3,2);
    matrix<double> rot_mat(3,3);
    for(size_t fi = 0; fi < bi_.cut_face_pairs.size(2)/2; ++fi){
        const size_t rot_type_st =
            get_face_rot_type(bi_.cut_face_pairs(colon(),fi * 2 + 0),
                              bi_.cut_face_pairs(colon(),fi * 2 + 1),
                              *bi_.fa_cut, inner_type);
        gap_node_pair(colon(),0) = 3 * fi;
        gap_node_pair(colon(),1) = 3 * (fi + bi_.cut_face_pairs.size(2)/2);
        gap_node_pair(0,colon()) += 1;
        gap_node_pair(1,colon()) += 2;
        gap_node_pair(2,colon()) += 3;
        rot_mat = type_transition2(rot_type_st); // rot_ts = trans(rot_st), and rot_ts should be transed if it * g
        gap_node_pair(colon(),0) = temp(rot_mat * gap_node_pair(colon(),0));
        for(size_t di = 0; di < 3; ++di){

            jtf::algorithm::equation<double> eqn;
            eqn.add_expression(jtf::algorithm::make_expression(
                                 abs(gap_node_pair(di,1))-1,
                                 jtf::math::get_sign(gap_node_pair(di,1))));
            eqn.add_expression(jtf::algorithm::make_expression(
                                 abs(gap_node_pair(di,0))-1,
                                 jtf::math::get_sign(gap_node_pair(di,0))));
            tm_eqn.start();
            ge_gap->add_equation(eqn);
            tm_eqn.finish();
            tt += tm_eqn.result();
            tt_idx += 1;
          }
      }
    tim.finish();
    cerr << "# [info] gauss eliminator time " << tim.result() << endl;
  }

  cerr << "# [info] gauss eliminator eqn number " << tt_idx << endl;
  cerr << "# [info] gauss eliminator sec/eqn " << tt/tt_idx << endl;
  {
    ofstream ofs("debug_gap_eqn");
    for(const auto & one_eqn : *ge_gap){
        ofs << "eqn" << endl;
        for(const auto & one_exp : one_eqn){
            ofs << one_exp.index << " ";
          }
        ofs <<  endl;
        for(const auto & one_exp : one_eqn){
            ofs << one_exp.coefficient << " ";
          }
        ofs << endl;
      }
  }

  vector<group<size_t> > gnode_group(gnode_.size());
  vector<size_t> gnode_to_group(gnode_.size());

  timer tim;
  tim.start();
  extract_groups_from_equations(*ge_gap, gnode_, gnode_flag, gnode_group, gnode_to_group);

  itr_matrix<const size_t*> g_mat(3, gnode_to_group.size()/3, &gnode_to_group[0]);

  map<vector<size_t>, vector<size_t> > g_coord_2_groups;
  vector<size_t> one_coord(3);
  for(size_t gi = 0; gi < g_mat.size(2); ++gi){
      std::copy(g_mat(colon(),gi).begin(), g_mat(colon(),gi).end(), one_coord.begin());
      g_coord_2_groups[one_coord].push_back(gi);
    }

  vector<size_t> g2group(g_mat.size(2));
  size_t group_idx = 0;
  for(const auto & one_group : g_coord_2_groups){
      for(const auto & idx: one_group.second){
          g2group[idx] = group_idx;
        }
      ++group_idx;
    }

  tim.finish();
  cerr << "# [info] gap groups extraction time " << tim.result() << endl;

  {
    tim.start();
    map<pair<size_t,size_t>, vector<pair<size_t,size_t> > > rot_gap2_face_pairs;
    for(size_t fi = 0; fi < bi_.cut_face_pairs.size(2)/2; ++fi){
        const size_t r_st = get_face_rot_type
            (bi_.cut_face_pairs(colon(),2*fi+0),
             bi_.cut_face_pairs(colon(),2*fi+1), *(bi_.fa_cut), inner_type);
        const size_t r_ts = get_trans_type(r_st);
        const size_t g_st = g2group[fi];
        const size_t g_ts = g2group[fi + bi_.cut_face_pairs.size(2)/2];

        auto it = rot_gap2_face_pairs.find(make_pair(r_st, g_st));
        if(it == rot_gap2_face_pairs.end()){
            it = rot_gap2_face_pairs.find(make_pair(r_ts, g_ts));
            if(it == rot_gap2_face_pairs.end())
              rot_gap2_face_pairs[make_pair(r_st, g_st)].push_back(make_pair(2*fi+0,2*fi+1));
            else{
                it->second.push_back(make_pair(2*fi+1,2*fi+0));
              }
          }else{
            it->second.push_back(make_pair(2*fi+0,2*fi+1));
          }
      }
    for(const auto & one_patch : rot_gap2_face_pairs){
        bi_.cut_face_patches_.push_back(one_patch.second);
        bi_.cut_face_patch_rot_type_.push_back(one_patch.first.first);
      }
    tim.finish();
    cerr << "# [info] patch extraction time " << tim.result() << endl;
  }
  {
    ofstream ofs("cut_patch_face.vtk");
    tri2vtk(ofs, &tet_node[0], tet_node.size(2), &bi_.cut_face_pairs[0], bi_.cut_face_pairs.size(2));
    vector<size_t> cut_patch(bi_.cut_face_pairs.size(2));
    for(size_t pi = 0; pi < bi_.cut_face_patches_.size(); ++pi){
        for(const auto & face_pair : bi_.cut_face_patches_[pi]){
            cut_patch[face_pair.first] = pi;
            cut_patch[face_pair.second] = pi;
          }
      }
    cell_data(ofs, &cut_patch[0], cut_patch.size(), "patch");
  }
}

void transition_elimination::cluster_transition_patch_each(
    const boost::unordered_map<pair<size_t,size_t>,size_t> & inner_type,
    const zjucad::matrix::matrix<double> * node_ptr)
{// this function takes such assumption, in cut_face_pairs,
  // all 2i,2i+1 faces are paris
  //face_stack.push(0); // seed face idx in cut_face_pairs_
  std::vector<pair<size_t,size_t> > one_patch; // recorded one pair of faces

  for(size_t fi = 0; fi < bi_.cut_face_pairs.size(2)/2; ++fi){
      one_patch.clear();
      one_patch.push_back(make_pair(2*fi,2*fi+1));

      size_t rot_type_of_face_pair = get_face_rot_type
          (bi_.cut_face_pairs(colon(),2*fi),
           bi_.cut_face_pairs(colon(),2*fi+1), *(bi_.fa_cut), inner_type);
      bi_.cut_face_patches_.push_back(one_patch);
      bi_.cut_face_patch_rot_type_.push_back(rot_type_of_face_pair);
    }

  if(node_ptr){
      {// debug
        ofstream ofs("cut_patch.vtk");
        matrix<size_t> cut_face_pairs_type(bi_.cut_face_pairs.size(2),1);
        for(size_t pi = 0; pi < bi_.cut_face_patches_.size(); ++pi){
            const vector<pair<size_t,size_t> > & one_patch = bi_.cut_face_patches_[pi];
            for(const auto & one_pair : one_patch){
                cut_face_pairs_type[one_pair.first] = pi;
                cut_face_pairs_type[one_pair.second] = pi;
              }
          }

        tri2vtk(ofs, &(*node_ptr)[0], node_ptr->size(2),
            &bi_.cut_face_pairs[0], bi_.cut_face_pairs.size(2));
        cell_data(ofs, &cut_face_pairs_type[0], cut_face_pairs_type.size(), "type");
      }

      {// debug
        ofstream ofs("cut_patch_type.vtk");
        matrix<size_t> cut_face_pairs_type(bi_.cut_face_pairs.size(2),1);
        for(size_t pi = 0; pi < bi_.cut_face_patches_.size(); ++pi){
            const vector<pair<size_t,size_t> > & one_patch = bi_.cut_face_patches_[pi];
            for(const auto & one_pair : one_patch){
                cut_face_pairs_type[one_pair.first] = bi_.cut_face_patch_rot_type_[pi];
                cut_face_pairs_type[one_pair.second] = bi_.cut_face_patch_rot_type_[pi];
              }
          }

        tri2vtk(ofs, &(*node_ptr)[0], node_ptr->size(2),
            &bi_.cut_face_pairs[0], bi_.cut_face_pairs.size(2));
        cell_data(ofs, &cut_face_pairs_type[0], cut_face_pairs_type.size(), "type");
      }

      //      {// debug
      //        // ofstream ofs("cut_patch_type_i.vtk");
      //        matrix<size_t> cut_face_pairs_type(bi_.cut_face_pairs.size(2),1);
      //        for(size_t pi = 0; pi < bi_.cut_face_patches_.size(); ++pi){
      //            const vector<pair<size_t,size_t> > & one_patch = bi_.cut_face_patches_[pi];
      //            stringstream ss;
      //            ss << "cut_patch_" << pi << ".vtk";
      //            ofstream ofs(ss.str().c_str());
      //            vector<size_t> cut_face_of_one_patch;
      //            vector<size_t> type;
      //            size_t ffi = 0;
      //            for(const auto & one_pair : one_patch){
      //                cut_face_of_one_patch.insert(cut_face_of_one_patch.end(),
      //                                             bi_.cut_face_pairs(colon(),one_pair.first).begin(),
      //                                             bi_.cut_face_pairs(colon(),one_pair.first).end());
      //                cut_face_of_one_patch.insert(cut_face_of_one_patch.end(),
      //                                             bi_.cut_face_pairs(colon(),one_pair.second).begin(),
      //                                             bi_.cut_face_pairs(colon(),one_pair.second).end());
      //                type.push_back(ffi);
      //                type.push_back(ffi);
      //                ++ffi;
      //              }
      //            tri2vtk(ofs, &(*node_ptr)[0], node_ptr->size(2),
      //                &cut_face_of_one_patch[0], cut_face_of_one_patch.size()/3);
      //            cell_data(ofs, &type[0], type.size(), "type");
      //          }


      //      }
    }
}

void transition_elimination::cluster_transition_patch(
    const boost::unordered_map<pair<size_t,size_t>,size_t> & inner_type,
    const zjucad::matrix::matrix<double> * node_ptr)
{// this function takes such assumption, in cut_face_pairs,
  // all 2i,2i+1 faces are paris
  std::vector<bool> face_visited(bi_.cut_face_pairs.size(2), false);
  stack<size_t> face_stack;
  //face_stack.push(0); // seed face idx in cut_face_pairs_
  std::vector<pair<size_t,size_t> > one_patch; // recorded one pair of faces
  size_t rot_type;
  while(1){
      auto it = find(face_visited.begin(), face_visited.end(), false);
      if(it == face_visited.end()) break;
      assert(face_stack.empty());
      face_stack.push(it-face_visited.begin());
      //it->second = true; // label is as visited

      one_patch.clear();
      rot_type = -1;
      while(!face_stack.empty()){
          const size_t face_idx = face_stack.top();
          face_stack.pop();
          if(face_visited[face_idx]) continue;

          const size_t face_idx_pair = face_idx%2==0?face_idx+1:face_idx-1;

          size_t rot_type_of_face_pair = get_face_rot_type
              (bi_.cut_face_pairs(colon(),face_idx),
               bi_.cut_face_pairs(colon(),face_idx_pair), *(bi_.fa_cut), inner_type);
          if(rot_type == -1)
            rot_type = rot_type_of_face_pair;

          if(rot_type != rot_type_of_face_pair){ // meet different rot type
              continue;
            }else{
              // check each adjacent face
              for(size_t pi = 0; pi < bi_.cut_face_pairs.size(1); ++pi){
                  const size_t edge_idx =
                      bi_.ea_cut_face_pair->get_edge_idx(
                        bi_.cut_face_pairs(pi,face_idx),
                        bi_.cut_face_pairs((pi+1)%bi_.cut_face_pairs.size(1), face_idx));
                  assert(edge_idx != -1);
                  const pair<size_t,size_t> &two_cells = bi_.ea_cut_face_pair->edge2cell_[edge_idx];
                  const pair<size_t,size_t> & one_edge = bi_.ea_cut_face_pair->edges_[edge_idx];

                  assert(two_cells.first == face_idx || two_cells.second == face_idx);
                  const size_t other_face_idx = two_cells.first + two_cells.second - face_idx;
                  if(other_face_idx == -1) continue;
                  face_stack.push(other_face_idx);
                }
              //              // check each adjacent face of pairded face
              //              // adjacent face of "adjacent face of paired face" also need to be checked
              //              for(size_t pi = 0; pi < bi_.cut_face_pairs.size(1); ++pi){
              //                  const size_t edge_idx =
              //                      bi_.ea_cut_face_pair->get_edge_idx(
              //                        bi_.cut_face_pairs(pi,face_idx_pair),
              //                        bi_.cut_face_pairs((pi+1)%bi_.cut_face_pairs.size(1), face_idx_pair));
              //                  assert(edge_idx != -1);
              //                  const pair<size_t,size_t> &two_cells = bi_.ea_cut_face_pair->edge2cell_[edge_idx];
              //                  const pair<size_t,size_t> & one_edge = bi_.ea_cut_face_pair->edges_[edge_idx];

              //                  assert(two_cells.first == face_idx_pair || two_cells.second == face_idx_pair);
              //                  const size_t other_face_idx = two_cells.first + two_cells.second - face_idx_pair;
              //                  if(other_face_idx == -1) continue;

              //                  const size_t face_idx_pair_pair = other_face_idx%2==0?other_face_idx+1:other_face_idx-1;

              //                  face_stack.push(face_idx_pair_pair);
              //                }
            }

          face_visited[face_idx] = true;
          face_visited[face_idx_pair] = true;

          one_patch.push_back(make_pair(face_idx, face_idx_pair));
        }
      bi_.cut_face_patches_.push_back(one_patch);
      bi_.cut_face_patch_rot_type_.push_back(rot_type);
    }

  //  if(node_ptr){
  //      {// debug
  //        ofstream ofs("cut_patch.vtk");
  //        matrix<size_t> cut_face_pairs_type(bi_.cut_face_pairs.size(2),1);
  //        for(size_t pi = 0; pi < bi_.cut_face_patches_.size(); ++pi){
  //            const vector<pair<size_t,size_t> > & one_patch = bi_.cut_face_patches_[pi];
  //            for(const auto & one_pair : one_patch){
  //                cut_face_pairs_type[one_pair.first] = pi;
  //                cut_face_pairs_type[one_pair.second] = pi;
  //              }
  //          }

  //        tri2vtk(ofs, &(*node_ptr)[0], node_ptr->size(2),
  //            &bi_.cut_face_pairs[0], bi_.cut_face_pairs.size(2));
  //        cell_data(ofs, &cut_face_pairs_type[0], cut_face_pairs_type.size(), "type");
  //      }

  //      {// debug
  //        ofstream ofs("cut_patch_type.vtk");
  //        matrix<size_t> cut_face_pairs_type(bi_.cut_face_pairs.size(2),1);
  //        for(size_t pi = 0; pi < bi_.cut_face_patches_.size(); ++pi){
  //            const vector<pair<size_t,size_t> > & one_patch = bi_.cut_face_patches_[pi];
  //            for(const auto & one_pair : one_patch){
  //                cut_face_pairs_type[one_pair.first] = bi_.cut_face_patch_rot_type_[pi];
  //                cut_face_pairs_type[one_pair.second] = bi_.cut_face_patch_rot_type_[pi];
  //              }
  //          }

  //        tri2vtk(ofs, &(*node_ptr)[0], node_ptr->size(2),
  //            &bi_.cut_face_pairs[0], bi_.cut_face_pairs.size(2));
  //        cell_data(ofs, &cut_face_pairs_type[0], cut_face_pairs_type.size(), "type");
  //      }

  //      {// debug
  //        // ofstream ofs("cut_patch_type_i.vtk");
  //        matrix<size_t> cut_face_pairs_type(bi_.cut_face_pairs.size(2),1);
  //        for(size_t pi = 0; pi < bi_.cut_face_patches_.size(); ++pi){
  //            const vector<pair<size_t,size_t> > & one_patch = bi_.cut_face_patches_[pi];
  //            stringstream ss;
  //            ss << "cut_patch_" << pi << ".vtk";
  //            ofstream ofs(ss.str().c_str());
  //            vector<size_t> cut_face_of_one_patch;
  //            vector<size_t> type;
  //            size_t ffi = 0;
  //            for(const auto & one_pair : one_patch){
  //                cut_face_of_one_patch.insert(cut_face_of_one_patch.end(),
  //                                             bi_.cut_face_pairs(colon(),one_pair.first).begin(),
  //                                             bi_.cut_face_pairs(colon(),one_pair.first).end());
  //                cut_face_of_one_patch.insert(cut_face_of_one_patch.end(),
  //                                             bi_.cut_face_pairs(colon(),one_pair.second).begin(),
  //                                             bi_.cut_face_pairs(colon(),one_pair.second).end());
  //                type.push_back(ffi);
  //                type.push_back(ffi);
  //                ++ffi;
  //              }
  //            tri2vtk(ofs, &(*node_ptr)[0], node_ptr->size(2),
  //                &cut_face_of_one_patch[0], cut_face_of_one_patch.size()/3);
  //            cell_data(ofs, &type[0], type.size(), "type");
  //          }
  //      }
  //    }
}

int  transition_elimination::init_with_gap(
    const zjucad::matrix::matrix<size_t> & uncut_mesh,
    const zjucad::matrix::matrix<size_t> & cut_mesh,
    const zjucad::matrix::matrix<size_t> & cut_tet2tet,
    const boost::unordered_map<std::pair<size_t,size_t>,size_t> & inner_face_type,
    const boost::unordered_map<size_t,size_t> & surface_type,
    const bool is_restricted_type,
    const zjucad::matrix::matrix<double> &node)
{
  construct_basic_info_with_gaps(uncut_mesh, cut_mesh, cut_tet2tet, inner_face_type, node);

  //cluster_transition_patch(inner_face_type,&node);
  cluster_transition_patch_by_elimination(uncut_mesh, cut_mesh, cut_tet2tet, node, inner_face_type);

  build_transitions(uncut_mesh, cut_mesh, cut_tet2tet, inner_face_type,
                    surface_type, is_restricted_type);
}

void transition_elimination::get_equation(
    std::vector<std::vector<size_t> > & eqn_idx,
    std::vector<std::vector<double> > & eqn_coeff)const
{
  map<size_t, vector<jtf::algorithm::gauss_eliminator<double>::equation_ptr> > A;
  for(jtf::algorithm::gauss_eliminator<double>::equation_ptr it = ge->begin(); it != ge->end(); ++it){
      A[it->get_prime_idx()].push_back(it);
    }

  eqn_idx.clear();
  eqn_coeff.clear();

  eqn_idx.reserve(A.size());
  eqn_coeff.reserve(A.size());


  vector<size_t> one_eqn_idx;
  vector<double> one_eqn_coeff;

  for(const auto &it : A){
      const vector<jtf::algorithm::gauss_eliminator<double>::equation_ptr> &one_eqn_vec = it.second;
      one_eqn_coeff.clear();
      one_eqn_idx.clear();

      for(size_t ei = 0; ei < one_eqn_vec.size(); ++ei){
          const jtf::algorithm::equation<double> & one_eqn = *(one_eqn_vec[ei]);

          for(const auto & one_exp : one_eqn){
              one_eqn_idx.push_back(one_exp.index);
              one_eqn_coeff.push_back(one_exp.coefficient);
            }
        }
      eqn_idx.push_back(one_eqn_idx);
      eqn_coeff.push_back(one_eqn_coeff);
    }
}

int transition_elimination::get_gaps_stupid(
    const zjucad::matrix::matrix<size_t> &tetmesh,
    const zjucad::matrix::matrix<size_t> &cut_mesh,
    const zjucad::matrix::matrix<size_t> &cut_tet2tet,
    const boost::unordered_map<std::pair<size_t, size_t>, size_t> &inner_face_type,
    const zjucad::matrix::matrix<double> & node)
{
  typedef vector<size_t> one_face;

  deque<pair<size_t,size_t> > edges;
  for(const auto & edge_tet_loop: bi_.ortae.e2t_){
      edges.push_back(edge_tet_loop.first);
    }

  while(1){
      const size_t edge_size = edges.size();
      for(size_t i = 0; i < edge_size; ++i){
          const pair<size_t,size_t>  &one_edge = edges[i];
          const auto it = bi_.ortae.e2t_.find(one_edge);
          assert(it != bi_.ortae.e2t_.end());
          const vector<size_t> & tet_loop = it->second;
          if(!bi_.ortae.is_inner_edge(tet_loop)) continue;

          vector<matrix<double> > gnode;
          matrix<double> total_rot = eye<double>(3);
          matrix<double> one_gnode;
          size_t edge_type = get_edge_type_with_part_face_type_map(tet_loop, inner_face_type);

          for(size_t i = tet_loop.size() - 2; i != -1; --i){
              int rtn = bi_.get_gnode_idx(make_pair(tet_loop[i], tet_loop[i+1]),one_gnode);
              one_gnode += 1.0; // to avoid jtf::math::get_sign(0)
              if(rtn != 0 && rtn != 1) continue;
              if(rtn == 1){
                  // f_t = Pi_st*f_s+g_st ==> f_s = Pi_ts*f_t-Pi_ts*g_st
                  // Goal: g_st; Get: g_ts
                  // -Pi_ts*g_st=g_ts ==> g_st = -Pi_st*g_ts
                  size_t type = bi_.get_inner_face_type(
                        make_pair(tet_loop[i],tet_loop[i+1]), inner_face_type);
                  assert(type < 24); // only 24 roation typeswe
                  matrix<double> rot_mat = -1*trans(type_transition2(type));
                  one_gnode = temp(rot_mat * one_gnode);
                }
              size_t forward_type = bi_.get_inner_face_type(
                    make_pair(tet_loop[i+1], tet_loop[(i+2)%tet_loop.size()]),
                  inner_face_type);
              total_rot = temp(total_rot * trans(type_transition2(forward_type)));
              one_gnode = temp(total_rot * one_gnode);
              gnode.push_back(one_gnode);
            }

          vector<size_t> zero_variant;
          if(edge_type == TRIVIAL_TYPE){
              zero_variant.push_back(0); // three variants are zeros
              zero_variant.push_back(1);
              zero_variant.push_back(2);
            }else{
              size_t axis_type = axis_to_around(edge_type);
              if(axis_type == -1){ // compound edge
                  zero_variant.push_back(0);
                  zero_variant.push_back(1);
                }else
                zero_variant.push_back(axis_type);
            }

          bool is_equation_solved = true;
          for(size_t di = 0; di < zero_variant.size(); ++di){
              jtf::algorithm::equation<double> eqn;
              for(size_t gi = 0; gi < gnode.size(); ++gi){
                  if(gnode_flag_[abs(number_rounding(gnode[gi][zero_variant[di]]))-1]){
                      assert(fabs(gnodes_[abs(number_rounding(gnode[gi][zero_variant[di]]))-1])<1e-6);
                      continue;
                    }
                  eqn.add_expression(jtf::algorithm::make_expression(
                                       abs(number_rounding(gnode[gi][zero_variant[di]]))-1,
                                     jtf::math::get_sign(gnode[gi][zero_variant[di]])));
                  is_equation_solved = false;
                }
              eqn.standardization();
              if(eqn.state() == 0) continue;
              if(eqn.state() == 1) {
                  gnodes_[eqn.get_prime_idx()] = 0; // assert
                  gnode_flag_[eqn.get_prime_idx()] = true;
                  continue;
                }
            }
          if(!is_equation_solved) edges.push_back(one_edge);
        }
      edges.erase(edges.begin(), edges.begin() + edge_size);
      if(edges.size() == edge_size) break;
    }

  {
    ofstream ofs("gap_unknown_face.vtk");
    vector<size_t> faces;
    for(const auto & tet_pair_idx : bi_.tet_pair2g_idx){
        if(!gnode_flag_[3 * tet_pair_idx.second + 0] ||
           !gnode_flag_[3 * tet_pair_idx.second + 1] ||
           !gnode_flag_[3 * tet_pair_idx.second + 2]){
            one_face of(3);
            jtf::mesh::find_common_face(tetmesh(colon(), tet_pair_idx.first.first),
                                        tetmesh(colon(), tet_pair_idx.first.second),
                                        &of[0]);
            faces.insert(faces.end(), of.begin(), of.end());
          }
      }
    tri2vtk(ofs, &node[0], node.size(2), &faces[0], faces.size()/3);
  }
  return 0;
}

int transition_elimination::get_gaps_elimination(
    const zjucad::matrix::matrix<size_t> &tetmesh,
    const zjucad::matrix::matrix<size_t> &cut_mesh,
    const zjucad::matrix::matrix<size_t> &cut_tet2tet,
    const boost::unordered_map<std::pair<size_t, size_t>, size_t> &inner_face_type,
    const zjucad::matrix::matrix<double> & node)
{
  for(const auto & one_edge2tet_loop : bi_.ortae.e2t_){

      const vector<size_t> & tet_loop = one_edge2tet_loop.second;
      if(!bi_.ortae.is_inner_edge(tet_loop)) continue;

      vector<matrix<double> > gnode;
      matrix<double> total_rot = eye<double>(3);
      matrix<double> one_gnode;
      size_t edge_type = get_edge_type_with_part_face_type_map(tet_loop, inner_face_type);

      for(size_t i = tet_loop.size() - 2; i != -1; --i){
          size_t forward_type = bi_.get_inner_face_type(
                make_pair(tet_loop[i+1], tet_loop[(i+2)%tet_loop.size()]),
              inner_face_type);
          total_rot = temp(total_rot * trans(type_transition2(forward_type)));

          int rtn = bi_.get_gnode_idx(make_pair(tet_loop[i], tet_loop[i+1]),one_gnode);
          one_gnode += 1.0; // to avoid jtf::math::get_sign(0)
          if(rtn != 0 && rtn != 1) continue;
          if(rtn == 1){
              // f_t = Pi_st*f_s+g_st ==> f_s = Pi_ts*f_t-Pi_ts*g_st
              // Goal: g_st; Get: g_ts
              // -Pi_ts*g_st=g_ts ==> g_st = -Pi_st*g_ts
              size_t type_ts = bi_.get_inner_face_type(
                    make_pair(tet_loop[i],tet_loop[i+1]), inner_face_type);
              assert(type_ts < 24); // only 24 roation typeswe
              matrix<double> rot_mat = -1*trans(type_transition2(type_ts));
              one_gnode = temp(rot_mat * one_gnode);
            }

          one_gnode = temp(total_rot * one_gnode);
          gnode.push_back(one_gnode);
        }

      vector<size_t> zero_variant;
      if(edge_type == TRIVIAL_TYPE){
          zero_variant.push_back(0); // three variants are zeros
          zero_variant.push_back(1);
          zero_variant.push_back(2);
        }else{
          size_t axis_type = axis_to_around(edge_type);
          if(axis_type == -1){ // compound edge
              zero_variant.push_back(0);
              zero_variant.push_back(1);
            }else
            zero_variant.push_back(axis_type);
        }

      for(size_t di = 0; di < zero_variant.size(); ++di){
          jtf::algorithm::equation<double> eqn;
          for(size_t gi = 0; gi < gnode.size(); ++gi){
              eqn.add_expression(jtf::algorithm::make_expression(
                                   abs(number_rounding(gnode[gi][zero_variant[di]]))-1,
                                 jtf::math::get_sign(gnode[gi][zero_variant[di]])));
            }
          ge_gap->add_equation(eqn);
        }
    }

  {
    typedef vector<size_t> one_face;
    ofstream ofs("gap_unknown_face.vtk");
    vector<size_t> faces;
    for(const auto & tet_pair_idx : bi_.tet_pair2g_idx){
        if(!gnode_flag_[3 * tet_pair_idx.second + 0] ||
           !gnode_flag_[3 * tet_pair_idx.second + 1] ||
           !gnode_flag_[3 * tet_pair_idx.second + 2]){
            one_face of(3);
            jtf::mesh::find_common_face(tetmesh(colon(), tet_pair_idx.first.first),
                                        tetmesh(colon(), tet_pair_idx.first.second),
                                        &of[0]);
            faces.insert(faces.end(), of.begin(), of.end());
          }
      }
    tri2vtk(ofs, &node[0], node.size(2), &faces[0], faces.size()/3);
  }
  {
    ofstream ofs("unknown_gap");
    for(const auto & one_eqn : *ge_gap){
        ofs << "eqi" << endl;
        for(const auto & one_exp : one_eqn){
            ofs << one_exp.index << " " ;
          }
        ofs << endl;
        for(const auto & one_exp : one_eqn){
            ofs << one_exp.coefficient << " " ;
          }
        ofs << endl;
      }
  }
  {
    typedef vector<size_t> one_face;
    ofstream ofs("zero_gap_rot_face.vtk");
    vector<size_t> faces;
    for(const auto & tet_pair_idx : bi_.tet_pair2g_idx){
        if(gnode_flag_[3 * tet_pair_idx.second + 0] &&
           gnode_flag_[3 * tet_pair_idx.second + 1] &&
           gnode_flag_[3 * tet_pair_idx.second + 2]){
            const auto it = inner_face_type.find(
                  make_pair(tet_pair_idx.first.first,
                            tet_pair_idx.first.second));
            if(it == inner_face_type.end()) continue;
            if(it->second == TRIVIAL_TYPE) continue;

            one_face of(3);
            jtf::mesh::find_common_face(tetmesh(colon(), tet_pair_idx.first.first),
                                        tetmesh(colon(), tet_pair_idx.first.second),
                                        &of[0]);
            faces.insert(faces.end(), of.begin(), of.end());
          }
      }
    tri2vtk(ofs, &node[0], node.size(2), &faces[0], faces.size()/3);
  }

  {
    typedef vector<size_t> one_face;
    ofstream ofs("zero_gap_identity_face.vtk");
    vector<size_t> faces;
    for(const auto & tet_pair_idx : bi_.tet_pair2g_idx){
        if(gnode_flag_[3 * tet_pair_idx.second + 0] &&
           gnode_flag_[3 * tet_pair_idx.second + 1] &&
           gnode_flag_[3 * tet_pair_idx.second + 2]){
            const auto it = inner_face_type.find(
                  make_pair(tet_pair_idx.first.first,
                            tet_pair_idx.first.second));
            if(it != inner_face_type.end() &&
               it->second != TRIVIAL_TYPE) continue;

            one_face of(3);
            jtf::mesh::find_common_face(tetmesh(colon(), tet_pair_idx.first.first),
                                        tetmesh(colon(), tet_pair_idx.first.second),
                                        &of[0]);
            faces.insert(faces.end(), of.begin(), of.end());
          }
      }
    tri2vtk(ofs, &node[0], node.size(2), &faces[0], faces.size()/3);
  }

  return 0;
}

int transition_elimination::construct_basic_info(
    const matrix<size_t> & tetmesh,
    const matrix<size_t> & cut_mesh,
    const matrix<size_t> & cut_tet2tet,
    const matrix<double> & tet_node,
    std::vector<double> & node,
    boost::dynamic_bitset<> & node_flag,
    std::vector<double> & gnode,
    boost::dynamic_bitset<> & gnode_flag)
{
  bi_.fa.reset(jtf::mesh::face2tet_adjacent::create(tetmesh));
  bi_.fa_cut.reset(jtf::mesh::face2tet_adjacent::create(cut_mesh));
  if(!bi_.fa.get() || !bi_.fa_cut.get()){
      cerr << "# [error] can not build face2tet_adjacent." << endl;
      return __LINE__;
    }

  jtf::mesh::get_outside_face(*bi_.fa_cut, bi_.outside_face_cut);
  node_group_.resize(3*(max(cut_mesh)+1));
  node2group_.resize(3*(max(cut_mesh)+1));
  for(size_t ni = 0; ni < node_group_.size(); ++ni) {
      node_group_[ni] << ni;
      node2group_[ni] = ni;
    }

  node.resize(node_group_.size());
  node_flag.resize(node.size());
  ge.reset(new jtf::algorithm::gauss_eliminator<double>(node, node_flag));

  bi_.ortae.add_tets(tetmesh, *(bi_.fa.get()));
  bi_.ortae.sort_into_loop(tetmesh, tet_node, 1);

  vector<size_t> face(3), sort_cut_face(3);
  for(size_t fi = 0; fi < bi_.outside_face_cut.size(2); ++fi){
      for(size_t pi = 0; pi < bi_.outside_face_cut.size(1); ++pi)
        face[pi] = cut_tet2tet[bi_.outside_face_cut(pi,fi)];
      const size_t face_idx_orig = bi_.fa->get_face_idx(&face[0]);
      if(face_idx_orig == -1){
          cerr << "# [error] can not find face idx." << endl;
          return __LINE__;
        }
      if(bi_.fa->is_outside_face(bi_.fa->face2tet_[face_idx_orig])) continue;
      sort(face.begin(), face.end());
      for(size_t pi = 0; pi < bi_.outside_face_cut.size(1); ++pi)
        for(size_t pj = 0; pj < bi_.outside_face_cut.size(1); ++pj){
            if(face[pi] == cut_tet2tet[bi_.outside_face_cut(pj, fi)]){
                sort_cut_face[pi] = bi_.outside_face_cut(pj, fi);
                continue;
              }
          }
      bi_.orig_face2_cut_faces[face].push_back(sort_cut_face);
    }

  gnode.resize(bi_.orig_face2_cut_faces.size() * 3);
  gnode_flag.resize(gnode.size());
  ge_gap.reset(new jtf::algorithm::gauss_eliminator<double>(gnode,gnode_flag));

  size_t idx = 0;
  for(const auto & orig_face2cut_face: bi_.orig_face2_cut_faces){
      const vector<basic_info::one_face> & cut_face = orig_face2cut_face.second;
      pair<size_t,size_t> tet_pair =
          bi_.get_tet_pair_from_cut_face(cut_face.front(),  cut_face.back());
      bi_.tet_pair2g_idx[tet_pair] = idx;
      bi_.orig_face2g_idx[orig_face2cut_face.first] = idx;
      ++idx;
    }
  return 0;
}

int transition_elimination::add_inner_transition(
    const matrix<size_t> &tetmesh,
    const matrix<size_t> &cut_mesh,
    const matrix<size_t> &cut_tet2tet,
    const boost::unordered_map<pair<size_t,size_t>,size_t> &inner_face_type)
{
  int rtn0 = add_inner_transition_edge2(tetmesh, cut_mesh, cut_tet2tet, inner_face_type);
  int rtn1 = add_inner_transition_gap(tetmesh, cut_mesh, cut_tet2tet, inner_face_type);

  if(!rtn0 && !rtn1)
    return 0;
  else return __LINE__;
}

void transition_elimination::add_inner_transition_gap2(
    const zjucad::matrix::matrix<size_t> & tetmesh,
    const zjucad::matrix::matrix<size_t> & cut_mesh,
    const zjucad::matrix::matrix<size_t> & cut_tet2tet,
    const size_t original_node_number_with_gaps)
{
  assert(original_node_number_with_gaps > bi_.cut_face_patches_.size());
  const size_t original_points_number = original_node_number_with_gaps - bi_.cut_face_patches_.size();
  matrix<double> variable_s(3,1), variable_t(3,1), gap(3,1);
  matrix<double> rot(3,3);
  for(size_t pi = 0; pi < bi_.cut_face_patches_.size(); ++pi){
      const vector<pair<size_t,size_t> > & one_patch = bi_.cut_face_patches_[pi];
      const size_t rot_type = bi_.cut_face_patch_rot_type_[pi];
      gap[0] = 3 * (original_points_number + pi);
      gap[1] = gap[0] + 1;
      gap[2] = gap[1] + 1;
      rot = trans(type_transition2(rot_type));

      for(size_t fi = 0; fi < one_patch.size(); ++fi){
          const pair<size_t,size_t> & one_face_pair = one_patch[fi];
          for(size_t pj = 0; pj < 3; ++pj){
              for(size_t di = 0; di < 3; ++di){
                  variable_s[di] = 3*bi_.cut_face_pairs(pj, one_face_pair.first)+di;
                  variable_t[di] = 3*bi_.cut_face_pairs(pj, one_face_pair.second)+di;
                }
              variable_s += 1; // to avoid 0 index
              variable_t += 1;
              variable_s = temp(rot * variable_s);
              for(size_t di = 0; di < 3; ++di){
                  jtf::algorithm::equation<double> eqn;

                  eqn.add_expression(jtf::algorithm::make_expression(
                                       abs(variable_s[di])-1,
                                       1.0*jtf::math::get_sign(variable_s[di])));
                  eqn.add_expression(jtf::algorithm::make_expression(
                                       abs(variable_t[di])-1,
                                       -1.0*jtf::math::get_sign(variable_t[di])));
                  eqn.add_expression(jtf::algorithm::make_expression(
                                       gap[di],1.0));
                  ge->add_equation(eqn);
                }
            }
        }
    }
}

void transition_elimination::add_inner_transition_gap3(
    const zjucad::matrix::matrix<size_t> & tetmesh,
    const zjucad::matrix::matrix<size_t> & cut_mesh,
    const zjucad::matrix::matrix<size_t> & cut_tet2tet,
    const boost::unordered_map<std::pair<size_t,size_t>,size_t> & inner_face_type,
    const size_t original_node_number)
{
  map<pair<size_t,size_t>,size_t> tet_pair2cut_patch_idx;
  const size_t original_points_number = original_node_number - bi_.cut_face_patches_.size();
  {
    for(size_t i = 0 ; i < bi_.cut_face_patches_.size(); ++i){
        const vector<pair<size_t,size_t> > & one_patch = bi_.cut_face_patches_[i];
        ///warning!!! this face pair index is for bi_.cut_face_pais
        for(const pair<size_t,size_t> & one_face_pair : one_patch){
            const size_t face_idx_0 = bi_.fa_cut->get_face_idx(&bi_.cut_face_pairs(0, one_face_pair.first));
            const size_t face_idx_1 = bi_.fa_cut->get_face_idx(&bi_.cut_face_pairs(0, one_face_pair.second));
            assert(face_idx_0 != -1 && face_idx_1 != -1);
            const pair<size_t,size_t> & tet_pair_0 = bi_.fa_cut->face2tet_[face_idx_0];
            const pair<size_t,size_t> & tet_pair_1 = bi_.fa_cut->face2tet_[face_idx_1];
            pair<size_t,size_t> tet_pair(tet_pair_0.first == -1?tet_pair_0.second:tet_pair_0.first,
                                         tet_pair_1.first == -1?tet_pair_1.second:tet_pair_1.first);
            tet_pair2cut_patch_idx[tet_pair] = i;
          }
      }
  }
  
  for(const auto & one_edge2tet_loop : bi_.ortae.e2t_){
      const pair<size_t,size_t> & one_edge = one_edge2tet_loop.first;
      const auto it = bi_.ortae.e2t_.find(one_edge);
      assert(it != bi_.ortae.e2t_.end());
      const vector<size_t> & tet_loop = it->second;
      if(!bi_.ortae.is_inner_edge(tet_loop)) continue;

      vector<matrix<double> > gnode;
      matrix<double> total_rot = eye<double>(3);
      matrix<double> one_gnode;
      //     matrix<double> one_fnode;
      size_t edge_type = get_edge_type_with_part_face_type_map(tet_loop, inner_face_type);

      for(size_t i = tet_loop.size() - 2; i != -1; --i){
          size_t forward_type = bi_.get_inner_face_type(
                make_pair(tet_loop[i+1], tet_loop[(i+2)%tet_loop.size()]),
              inner_face_type);
          total_rot = temp(total_rot * trans(type_transition2(forward_type)));

          int rtn = bi_.get_gnode_idx_of_patch(
                make_pair(tet_loop[i], tet_loop[i+1]), tet_pair2cut_patch_idx,
              one_gnode);
          one_gnode += 1.0 + 3 * original_points_number; // to avoid jtf::math::get_sign(0)
          if(rtn != 0 && rtn != 1) continue;
          if(rtn == 1){
              // f_t = Pi_st*f_s+g_st ==> f_s = Pi_ts*f_t-Pi_ts*g_st
              // Goal: g_st; Get: g_ts
              // -Pi_ts*g_st=g_ts ==> g_st = -Pi_st*g_ts
              size_t type = bi_.get_inner_face_type(
                    make_pair(tet_loop[i],tet_loop[i+1]), inner_face_type);
              assert(type < 24); // only 24 roation typeswe
              matrix<double> rot_mat = -1*trans(type_transition2(type));
              one_gnode = temp(rot_mat * one_gnode);
            }

          one_gnode = temp(total_rot * one_gnode);
          gnode.push_back(one_gnode);
        }

      vector<size_t> zero_variant;
      if(edge_type == TRIVIAL_TYPE){
          zero_variant.push_back(0); // three variants are zeros
          zero_variant.push_back(1);
          zero_variant.push_back(2);
        }else{
          size_t axis_type = axis_to_around(edge_type);
          if(axis_type == -1){ // compound edge
              zero_variant.push_back(0);
              zero_variant.push_back(1);
            }else
            zero_variant.push_back(axis_type);
        }

      for(size_t di = 0; di < zero_variant.size(); ++di){
          jtf::algorithm::equation<double> eqn;
          for(size_t gi = 0; gi < gnode.size(); ++gi){
              eqn.add_expression(jtf::algorithm::make_expression(
                                   abs(number_rounding(gnode[gi][zero_variant[di]]))-1,
                                 jtf::math::get_sign(gnode[gi][zero_variant[di]])));
            }
          ge->add_equation(eqn);
        }
    }
}

int transition_elimination::add_inner_transition_gap(
    const matrix<size_t> &tetmesh,
    const matrix<size_t> &cut_mesh,
    const matrix<size_t> &cut_tet2tet,
    const boost::unordered_map<pair<size_t,size_t>,size_t> &inner_face_type)
{
  vector<matrix<double> > face_point(6);
  for(const auto & orig_face2cut_face: bi_.orig_face2_cut_faces){
      const vector<basic_info::one_face> & cut_faces = orig_face2cut_face.second;
      const auto it = bi_.orig_face2g_idx.find(orig_face2cut_face.first);
      if(it == bi_.orig_face2g_idx.end()) {
          cerr << "# [error] can not find gnode." << endl;
          return __LINE__;
        }
      const size_t g_idx = it->second;

      pair<size_t,size_t> tet_pair =
          bi_.get_tet_pair_from_cut_face(cut_faces.front(), cut_faces.back());
      if(tet_pair.first == -1 || tet_pair.second == -1) return __LINE__;
      const size_t rot_type = bi_.get_inner_face_type(tet_pair, inner_face_type);

      const matrix<double> rot_st = -1.0*trans(type_transition2(rot_type));
      for(size_t pi = 0; pi < 3; ++pi){
          bi_.get_fnode_idx(cut_faces.front()[pi], face_point[pi]);
          face_point[pi] += 1.0;
          face_point[pi] = temp(rot_st * face_point[pi]);
          bi_.get_fnode_idx(cut_faces.back()[pi], face_point[pi+3]);
          face_point[pi+3] += 1.0;
        }

      for(size_t di = 0; di < 3; ++di){
          if(gnode_flag_[3 * g_idx + di]){
              for(size_t pi = 0; pi < 3; ++pi){
                  jtf::algorithm::equation<double> eqn;
                  eqn.add_expression(
                        jtf::algorithm::make_expression(abs(number_rounding(face_point[pi][di]))-1,
                                                        jtf::math::get_sign(face_point[pi][di])));
                  eqn.add_expression(
                        jtf::algorithm::make_expression(abs(number_rounding(face_point[pi+3][di]))-1,
                        jtf::math::get_sign(face_point[pi+3][di])));
                  ge->add_equation(eqn);

                  eqn.standardization();
                  if(eqn.e_vec_.size() == 2 &&
                     fabs(eqn.e_vec_.front().coefficient +
                          eqn.e_vec_.back().coefficient) < 1e-6)
                    merge_node(eqn.e_vec_.front().index, eqn.e_vec_.back().index);
                }
            }else{
              for(size_t pi = 0; pi < 3; ++pi){
                  jtf::algorithm::equation<double> eqn;
                  eqn.add_expression(
                        jtf::algorithm::make_expression(abs(number_rounding(face_point[pi][di]))-1,
                                                        jtf::math::get_sign(face_point[pi][di])));
                  eqn.add_expression(
                        jtf::algorithm::make_expression(abs(number_rounding(face_point[(pi+1)%3][di]))-1,
                        -1*jtf::math::get_sign(face_point[(pi+1)%3][di])));

                  eqn.add_expression(
                        jtf::algorithm::make_expression(abs(number_rounding(face_point[pi+3][di]))-1,
                        jtf::math::get_sign(face_point[pi+3][di])));

                  eqn.add_expression(
                        jtf::algorithm::make_expression(abs(number_rounding(face_point[(pi+1)%3+3][di]))-1,
                        -1*jtf::math::get_sign(face_point[(pi+1)%3+3][di])));
                  ge->add_equation(eqn);

                  eqn.standardization();
                  if(eqn.e_vec_.size() == 2 &&
                     fabs(eqn.e_vec_.front().coefficient +
                          eqn.e_vec_.back().coefficient) < 1e-6)
                    merge_node(eqn.e_vec_.front().index, eqn.e_vec_.back().index);
                }
            }
        }
    }
  return 0;
}

int transition_elimination::add_inner_transition_edge(
    const matrix<size_t> &tetmesh,
    const matrix<size_t> &cut_mesh,
    const matrix<size_t> &cut_tet2tet,
    const boost::unordered_map<pair<size_t,size_t>,size_t> &inner_face_type)
{
  typedef vector<size_t> one_face;
  matrix<double> rot_mat = eye<double>(3);
  vector<matrix<double> > face_point(6);
  // for each inner face, build transition equations
  for(const auto & face2cut_face : bi_.orig_face2_cut_faces){

      const vector<one_face> & cut_faces = face2cut_face.second;
      if(cut_faces.size() != 2){
          cerr << "# [error] can not find two cut faces" << endl;
          return __LINE__;
        }

      pair<size_t,size_t> tet_pair =
          bi_.get_tet_pair_from_cut_face(cut_faces.front(),cut_faces.back());
      if(tet_pair.first == -1 || tet_pair.second == -1) return __LINE__;

      // regard face point variants as float variant, so that I can multiple rot
      for(size_t pi = 0; pi < 3; ++pi){
          bi_.get_fnode_idx(cut_faces.front()[pi], face_point[pi]);
          bi_.get_fnode_idx(cut_faces.back()[pi], face_point[pi+3]);
          face_point[pi] += 1.0;
          face_point[pi+3] += 1.0; // to avoid jtf::math::get_sign(0)
        }

      for(size_t pi = 0; pi < cut_faces[0].size(); ++pi){
          int point_order = 1;
          const pair<size_t,size_t> one_edge(
                cut_faces[0][pi],
              cut_faces[0][(pi+1)%cut_faces[0].size()]);

          auto ortae_it = bi_.ortae.e2t_.find(make_pair(cut_tet2tet[one_edge.first],
                                              cut_tet2tet[one_edge.second]));
          if(ortae_it == bi_.ortae.e2t_.end()){
              ortae_it = bi_.ortae.e2t_.find(make_pair(cut_tet2tet[one_edge.second],
                                             cut_tet2tet[one_edge.first]));
              point_order = -1;
            }
          if(ortae_it == bi_.ortae.e2t_.end()){
              cerr << "# [error] can not find edge in ortae.e2t_" << endl;
              return __LINE__;
            }
          vector<size_t> tet_loop = ortae_it->second;
          if(!bi_.ortae.is_inner_edge(tet_loop)) continue;

          {  // adjust tet_loop order
            tet_loop.pop_back();
            auto tet_loop_mid = find(tet_loop.begin(), tet_loop.end(), tet_pair.first);
            if(tet_loop_mid == tet_loop.end()) {
                cerr << "# [error] can not find tet in tet_loop: " << endl;
                return __LINE__;
              }
            std::rotate(tet_loop.begin(), tet_loop_mid, tet_loop.end());
            tet_loop.push_back(tet_loop.front());
          }

          const size_t edge_type =
              get_edge_type_with_part_face_type_map(tet_loop,inner_face_type);

          if(edge_type == TRIVIAL_TYPE) continue;

          rot_mat = type_transition2(edge_type);

          if(point_order == -1) rot_mat = temp(trans(rot_mat));

          {
            matrix<double> f0_p0, f0_p1, f1_p0, f1_p1;
            matrix<double> rot_f0_p0, rot_f0_p1, rot_f1_p0,rot_f1_p1;
            {
              rot_f0_p1 = face_point[(pi+1)%cut_faces.front().size()];
              rot_f0_p0 = face_point[pi];

              rot_f1_p1 = face_point[3+(pi+1)%cut_faces.front().size()];
              rot_f1_p0 = face_point[3+pi];

              f0_p0 = rot_f0_p0;
              f0_p1 = rot_f0_p1;
              f1_p0 = rot_f1_p0;
              f1_p1 = rot_f1_p1;
            }

            rot_f0_p0 = temp(rot_mat * rot_f0_p0);
            rot_f0_p1 = temp(rot_mat * rot_f0_p1);

            matrix<double> trans_mat = trans(rot_mat);
            rot_f1_p0 = temp(trans_mat * rot_f1_p0);
            rot_f1_p1 = temp(trans_mat * rot_f1_p1);

            for(size_t di = 0; di < 3; ++di){
                jtf::algorithm::equation<double> eq0,eq1;
                eq0.add_expression(
                      jtf::algorithm::make_expression(
                        abs(number_rounding(f0_p1[di]))-1, jtf::math::get_sign(f0_p1[di])));
                eq0.add_expression(
                      jtf::algorithm::make_expression(
                        abs(number_rounding(f0_p0[di]))-1, -1.0*jtf::math::get_sign(f0_p0[di])));

                eq0.add_expression(
                      jtf::algorithm::make_expression(
                        abs(number_rounding(rot_f0_p1[di]))-1, -1.0*jtf::math::get_sign(rot_f0_p1[di])));
                eq0.add_expression(
                      jtf::algorithm::make_expression(
                        abs(number_rounding(rot_f0_p0[di]))-1, 1.0*jtf::math::get_sign(rot_f0_p0[di])));

                ge->add_equation(eq0);

                eq1.add_expression(
                      jtf::algorithm::make_expression(
                        abs(number_rounding(f1_p1[di]))-1, jtf::math::get_sign(f1_p1[di])));
                eq1.add_expression(
                      jtf::algorithm::make_expression(
                        abs(number_rounding(f1_p1[di]))-1, -1.0*jtf::math::get_sign(f1_p1[di])));

                eq1.add_expression(
                      jtf::algorithm::make_expression(
                        abs(number_rounding(rot_f1_p1[di]))-1, -1.0*jtf::math::get_sign(rot_f1_p1[di])));
                eq1.add_expression(
                      jtf::algorithm::make_expression(
                        abs(number_rounding(rot_f1_p1[di]))-1, 1.0*jtf::math::get_sign(rot_f1_p1[di])));

                ge->add_equation(eq1);
              }
          }

          { // add merge info
            const size_t axis = axis_to_around(edge_type);
            merge_node(3*cut_faces.front()[pi]+(axis+1)%3,
                       3*cut_faces.front()[(pi+1)%3]+(axis+1)%3);
            integer_variable_.insert(3*cut_faces.front()[pi]+(axis+1)%3);
            integer_variable_.insert(3*cut_faces.front()[(pi+1)%3]+(axis+1)%3);

            merge_node(3*cut_faces.front()[pi]+(axis+2)%3,
                       3*cut_faces.front()[(pi+1)%3]+(axis+2)%3);
            integer_variable_.insert(3*cut_faces.front()[pi]+(axis+2)%3);
            integer_variable_.insert(3*cut_faces.front()[(pi+1)%3]+(axis+2)%3);

            merge_node(3*cut_faces.back()[pi]+(axis+1)%3,
                       3*cut_faces.back()[(pi+1)%3]+(axis+1)%3);
            integer_variable_.insert(3*cut_faces.back()[pi]+(axis+1)%3);
            integer_variable_.insert(3*cut_faces.back()[(pi+1)%3]+(axis+1)%3);


            merge_node(3*cut_faces.back()[pi]+(axis+2)%3,
                       3*cut_faces.back()[(pi+1)%3]+(axis+2)%3);
            integer_variable_.insert(3*cut_faces.back()[pi]+(axis+2)%3);
            integer_variable_.insert(3*cut_faces.back()[(pi+1)%3]+(axis+2)%3);
          }
        }
    }

  return 0;
}

int transition_elimination::add_inner_transition_edge2(
    const matrix<size_t> &tetmesh,
    const matrix<size_t> &cut_mesh,
    const matrix<size_t> &cut_tet2tet,
    const boost::unordered_map<pair<size_t,size_t>,size_t> &inner_face_type)
{
  typedef vector<size_t> one_face;

  vector<matrix<double> > face_point(6);
  // for each inner face, build transition equations
  for(const auto & face2cut_face : bi_.orig_face2_cut_faces){

      const vector<one_face> & cut_faces = face2cut_face.second;
      if(cut_faces.size() != 2){
          cerr << "# [error] can not find two cut faces" << endl;
          return __LINE__;
        }

      pair<size_t,size_t> tet_pair =
          bi_.get_tet_pair_from_cut_face(cut_faces.front(),cut_faces.back());
      if(tet_pair.first == -1 || tet_pair.second == -1) return __LINE__;

      // regard face point variants as float variant, so that I can multiple rot
      for(size_t pi = 0; pi < 3; ++pi){
          bi_.get_fnode_idx(cut_faces.front()[pi], face_point[pi]);
          bi_.get_fnode_idx(cut_faces.back()[pi], face_point[pi+3]);
          face_point[pi] += 1.0;
          face_point[pi+3] += 1.0; // to avoid jtf::math::get_sign(0)
        }

      for(size_t pi = 0; pi < cut_faces[0].size(); ++pi){
          int point_order = 1;
          const pair<size_t,size_t> one_edge(
                cut_faces[0][pi],
              cut_faces[0][(pi+1)%cut_faces[0].size()]);

          auto ortae_it = bi_.ortae.e2t_.find(make_pair(cut_tet2tet[one_edge.first],
                                              cut_tet2tet[one_edge.second]));
          if(ortae_it == bi_.ortae.e2t_.end()){
              ortae_it = bi_.ortae.e2t_.find(make_pair(cut_tet2tet[one_edge.second],
                                             cut_tet2tet[one_edge.first]));
              point_order = -1;
            }
          if(ortae_it == bi_.ortae.e2t_.end()){
              cerr << "# [error] can not find edge in ortae.e2t_" << endl;
              cerr << bi_.ortae.e2t_.size() << endl;
              return __LINE__;
            }
          vector<size_t> tet_loop = ortae_it->second;
          if(!bi_.ortae.is_inner_edge(tet_loop)) continue;

          {  // adjust tet_loop order
            tet_loop.pop_back();
            auto tet_loop_mid = find(tet_loop.begin(), tet_loop.end(), tet_pair.first);
            if(tet_loop_mid == tet_loop.end()) {
                cerr << "# [error] can not find tet in tet_loop: " << endl;
                return __LINE__;
              }
            std::rotate(tet_loop.begin(), tet_loop_mid, tet_loop.end());
            tet_loop.push_back(tet_loop.front());
          }

          const size_t edge_type =
              get_edge_type_with_part_face_type_map(tet_loop,inner_face_type);

          if(edge_type == TRIVIAL_TYPE) continue;

          const size_t axis_type = axis_to_around(edge_type);

          if(axis_type == -1){
              for(size_t axis_type_i = 0; axis_type_i < 3; ++axis_type_i){
                  jtf::algorithm::equation<double> eq0;
                  eq0.add_expression(
                        jtf::algorithm::make_expression(
                          3 *one_edge.first + axis_type_i, 1.0));
                  eq0.add_expression(
                        jtf::algorithm::make_expression(
                          3 *one_edge.second + axis_type_i, -1.0));
                  ge->add_equation(eq0);

                  integer_variable_.insert(3 *one_edge.first + axis_type_i);
                  integer_variable_.insert(3 *one_edge.second + axis_type_i);
                  merge_node(3 *one_edge.first + axis_type_i,
                             3 *one_edge.second + axis_type_i);
                }
              cerr << "# [error] meet compound edge." << endl;
            }else{

              for(size_t axis_type_i = 1; axis_type_i < 3; ++axis_type_i){
                  const size_t axis_type_idx = (axis_type + axis_type_i)%3;
                  jtf::algorithm::equation<double> eq0;
                  eq0.add_expression(
                        jtf::algorithm::make_expression(
                          3 *one_edge.first + axis_type_idx, 1.0));
                  eq0.add_expression(
                        jtf::algorithm::make_expression(
                          3 *one_edge.second + axis_type_idx, -1.0));
                  ge->add_equation(eq0);

                  integer_variable_.insert(3 *one_edge.first + axis_type_idx);
                  integer_variable_.insert(3 *one_edge.second + axis_type_idx);
                  merge_node(3 *one_edge.first + axis_type_idx,
                             3 *one_edge.second + axis_type_idx);
                }
            }


        }
    }

  return 0;
}

int transition_elimination::add_surface_transition(
    const zjucad::matrix::matrix<size_t> & tetmesh,
    const zjucad::matrix::matrix<size_t> & cut_mesh,
    const zjucad::matrix::matrix<size_t> & cut_tet2tet,
    const boost::unordered_map<size_t,size_t> & surface_type,
    const bool is_restricted_type)
{
  for(size_t fi = 0; fi < bi_.outside_face_cut.size(2); ++fi){
      const size_t face_idx = bi_.fa->get_face_idx(
            cut_tet2tet[bi_.outside_face_cut(0,fi)],
          cut_tet2tet[bi_.outside_face_cut(1,fi)],
          cut_tet2tet[bi_.outside_face_cut(2,fi)]);
      if(face_idx == -1){
          cerr << "# [error] can not find face idx." << endl;
          return __LINE__;
        }
      auto it = surface_type.find(face_idx);
      if(it == surface_type.end()) continue; // this face is not on surface.
      const size_t restricted_type = it->second;
      // combine face restricted type together
      if(is_restricted_type){
          for(size_t pi = 0; pi < bi_.outside_face_cut.size(1) -1; ++pi){
              merge_node(3 *bi_.outside_face_cut(pi, fi) + restricted_type,
                         3 * bi_.outside_face_cut(pi+1, fi) + restricted_type);
              integer_variable_.insert(3 *bi_.outside_face_cut(pi, fi) + restricted_type);
              integer_variable_.insert( 3 * bi_.outside_face_cut(pi+1, fi) + restricted_type);

              jtf::algorithm::equation<double> eq0;
              eq0.add_expression(
                    jtf::algorithm::make_expression(
                      3 *bi_.outside_face_cut(pi, fi) + restricted_type, 1.0));
              eq0.add_expression(
                    jtf::algorithm::make_expression(
                      3 * bi_.outside_face_cut(pi+1, fi) + restricted_type, -1.0));

              ge->add_equation(eq0);
            }
        }else{
          const size_t trans_type = get_trans_type(restricted_type);
          size_t rt = -1;
          matrix<double> ttm = trans(type_transition2(trans_type));
          for(size_t i = 0; i < 3; ++i)
            if(fabs(fabs(ttm(0,i)) -1) < 1e-6){
                rt = i; break;
              }
          for(size_t pi = 0; pi < bi_.outside_face_cut.size(1) -1; ++pi){
              merge_node(3 *bi_.outside_face_cut(pi, fi) + rt,
                         3 * bi_.outside_face_cut(pi+1, fi) + rt);
              integer_variable_.insert(3 *bi_.outside_face_cut(pi, fi) + rt);
              integer_variable_.insert(3 * bi_.outside_face_cut(pi+1, fi) + rt);

              jtf::algorithm::equation<double> eq0;
              eq0.add_expression(
                    jtf::algorithm::make_expression(
                      3 *bi_.outside_face_cut(pi, fi) + rt, 1.0));
              eq0.add_expression(
                    jtf::algorithm::make_expression(
                      3 * bi_.outside_face_cut(pi+1, fi) + rt, -1.0));
              ge->add_equation(eq0);
            }
        }
    }
  return 0;
}


//int transition_elimination::add_surface_transition_with_gaps(
//    const zjucad::matrix::matrix<size_t> & tetmesh,
//    const zjucad::matrix::matrix<size_t> & cut_mesh,
//    const zjucad::matrix::matrix<size_t> & cut_tet2tet,
//    const boost::unordered_map<size_t,size_t> & surface_type,
//    const bool is_restricted_type)
//{
//  for(size_t fi = 0; fi < bi_.outside_face_cut.size(2); ++fi){
//      const size_t face_idx = bi_.fa->get_face_idx(
//            cut_tet2tet[bi_.outside_face_cut(0,fi)],
//          cut_tet2tet[bi_.outside_face_cut(1,fi)],
//          cut_tet2tet[bi_.outside_face_cut(2,fi)]);
//      if(face_idx == -1){
//          cerr << "# [error] can not find face idx." << endl;
//          return __LINE__;
//        }
//      auto it = surface_type.find(face_idx);
//      if(it == surface_type.end()) continue; // this face is not on surface.
//      const size_t restricted_type = it->second;
//      // combine face restricted type together
//      if(is_restricted_type){
//          for(size_t pi = 0; pi < bi_.outside_face_cut.size(1) -1; ++pi){
//              merge_node(3 *bi_.outside_face_cut(pi, fi) + restricted_type,
//                         3 * bi_.outside_face_cut(pi+1, fi) + restricted_type);
//              integer_variable_.insert(3 *bi_.outside_face_cut(pi, fi) + restricted_type);
//              integer_variable_.insert( 3 * bi_.outside_face_cut(pi+1, fi) + restricted_type);
//            }
//        }else{
//          const size_t trans_type = get_trans_type(restricted_type);
//          size_t rt = -1;
//          matrix<double> ttm = type_transition2(trans_type);
//          for(size_t i = 0; i < 3; ++i)
//            if(fabs(fabs(ttm(0,i)) -1) < 1e-6){
//                rt = i; break;
//              }
//          for(size_t pi = 0; pi < bi_.outside_face_cut.size(1) -1; ++pi){
//              merge_node(3 *bi_.outside_face_cut(pi, fi) + rt,
//                         3 * bi_.outside_face_cut(pi+1, fi) + rt);
//              integer_variable_.insert(3 *bi_.outside_face_cut(pi, fi) + rt);
//              integer_variable_.insert(3 * bi_.outside_face_cut(pi+1, fi) + rt);
//            }
//        }
//    }
//  return 0;
//}


void transition_elimination::collect_group_variants_from_equations()
{
  // step0 : gather all equation with two variant
  ofstream ofs("node_equation");
  size_t eqi = 0;
  //  for(const auto & one_eqn : *ge){
  //      //ofs << one_eqn << endl;
  //      ofs << "eqn " << eqi++ << " " << one_eqn.e_vec_.size() << endl;
  //      for(const auto & one_exp : one_eqn)
  //        ofs << one_exp.index << " ";
  //      ofs << endl;
  //      for(const auto & one_exp : one_eqn)
  //        ofs << one_exp.coefficient << " ";
  //      ofs << endl;

  //      if(one_eqn.e_vec_.size() != 2) {
  //          continue;
  //        }

  //      if(fabs(one_eqn.e_vec_.front().coefficient +
  //              one_eqn.e_vec_.back().coefficient) < 1e-6){
  //          merge_node(one_eqn.e_vec_.front().index,
  //                     one_eqn.e_vec_.back().index);
  //        }
  //    }

  map<size_t, vector<jtf::algorithm::gauss_eliminator<double>::equation_ptr> > A;
  for( jtf::algorithm::gauss_eliminator<double>::equation_ptr it = ge->begin(); it != ge->end(); ++it){
      A[it->get_prime_idx()].push_back(it);
    }

  for(const auto &it: A){
      const vector<jtf::algorithm::gauss_eliminator<double>::equation_ptr> &one_eqn_vec = it.second;
      for(size_t ei = 0; ei < one_eqn_vec.size(); ++ei){
          const jtf::algorithm::equation<double> & one_eqn = *(one_eqn_vec[ei]);
          ofs << "eqn " << eqi++ << " " << one_eqn.e_vec_.size() << endl;
          for(const auto & one_exp : one_eqn)
            ofs << one_exp.index << " ";
          ofs << endl;
          for(const auto & one_exp : one_eqn)
            ofs << one_exp.coefficient << " ";
          ofs << endl;

          if(one_eqn.e_vec_.size() != 2) {
              continue;
            }

          if(fabs(one_eqn.e_vec_.front().coefficient +
                  one_eqn.e_vec_.back().coefficient) < 1e-6){
              merge_node(one_eqn.e_vec_.front().index,
                         one_eqn.e_vec_.back().index);
            }
        }
    }

  //step1 : gather all equations who imply equation constraints
  //  typedef pair<size_t,double> exp_pair;
  //  typedef set<exp_pair> expression_set;
  //  map<expression_set,set<size_t> > expression_map;
  //  for(const auto & one_eqn : *ge){
  //      const size_t major_idx = one_eqn.get_prime_idx();
  //      set<exp_pair> one_set;
  //      for(const auto & exp : one_eqn){
  //          if(exp.index == major_idx) continue;
  //          one_set.insert(make_pair(exp.index, exp.coefficient));
  //        }
  //      expression_map[one_set].insert(major_idx);
  //    }

  //  for(const auto & one_map_set : expression_map){
  //      if(one_map_set.second.size() == 1) continue;
  //      for(const auto & v_idx: one_map_set.second){
  //          merge_node(v_idx, *one_map_set.second.begin());
  //        }
  //    }

  //  ofstream ofs2("node_group");
  //  size_t gi = 0;
  //  for(const auto & g: node_group_){
  //      if(g.empty() || g.size() == 1) continue;
  //      ofs2 << "g " << gi++ << " " << g.size() << " " << 0 << endl; // TODO: 0: integer constraint, 1: no-interger constraints
  //      for(const auto & idx : g)
  //        ofs2 << idx << " ";
  //      ofs2 << endl << endl;
  //    }
}

void transition_elimination::merge_node(const size_t &from, const size_t &to)
{
  if(node2group_[from] == node2group_[to]) return;
  group<size_t> & from_group = node_group_[node2group_[from]];
  for(const auto & idx : from_group) node2group_[idx] = node2group_[to];

  node_group_[node2group_[to]].merge(from_group);
}


///////////////////////////////////////////////////////////////////////////////
singularity_graph * singularity_graph::create(
    const jtf::mesh::one_ring_tet_at_edge & ortae_cut,
    const jtf::mesh::one_ring_tet_at_edge & ortae,
    const matrixst & cut_tet2tet,
    const matrixst & outside_face_cut,
    const matrixst & cut_face_pair,
    const matrixst & outside_face_idx_cut,
    const boost::unordered_map<size_t,size_t> & surface_idx_to_rot_idx,
    const boost::unordered_map<pair<size_t,size_t>,size_t> & jump_face_to_rot_idx,
    const bool no_surface)
{
  unique_ptr<singularity_graph> sg(new singularity_graph);
  if(sg->init(ortae_cut,ortae, cut_tet2tet,outside_face_cut,cut_face_pair,
              outside_face_idx_cut, surface_idx_to_rot_idx,
              jump_face_to_rot_idx, no_surface)){
      return 0;
    }
  return sg.release();
}

int singularity_graph::init(
    const jtf::mesh::one_ring_tet_at_edge & ortae_cut,
    const jtf::mesh::one_ring_tet_at_edge & ortae,
    const matrixst & cut_tet2tet,
    const matrixst & outside_face_cut,
    const matrixst & cut_face_pair,
    const matrixst & outside_face_idx_cut,
    const boost::unordered_map<size_t,size_t> & surface_idx_to_rot_idx,
    const boost::unordered_map<pair<size_t,size_t>,size_t> & jump_face_to_rot_idx,
    const bool no_surface)
{
  assert(outside_face_cut.size(2) == cut_face_pair.size());

  boost::unordered_set<size_t> points_set(outside_face_cut.begin(),
                                          outside_face_cut.end());
  outside_points_.resize(points_set.size());
  copy(points_set.begin(), points_set.end(), outside_points_.begin());

  for(size_t t = 0; t < outside_points_.size(); ++t)
    outside_point_idx_.insert(make_pair(outside_points_[t], t));

  fnode_.resize(3 * outside_points_.size()); // init all fnode to separated groups
  groups_.resize(fnode_.size());

  //jump_faces_ = jump_faces;
  jump_faces_.reserve(jump_face_to_rot_idx.size());
  gnode_.resize(3 * jump_face_to_rot_idx.size(),0); // initial the gnodes between jump face
  gnode_flag_.resize(3 * jump_face_to_rot_idx.size()); // this means all gnode is not ready
  gnode2rot_idx_.resize(3 * jump_face_to_rot_idx.size());

  for(boost::unordered_map<pair<size_t,size_t>,size_t>::const_iterator it =
      jump_face_to_rot_idx.begin(); it != jump_face_to_rot_idx.end(); ++it){
      rot_idx2g_idx_[it->second] = jump_faces_.size();
      gnode2rot_idx_[jump_faces_.size()] = it->second;
      jump_face_to_gnode_idx_[it->first] = jump_faces_.size();
      //gnode_idx_to_jump_face_[jump_faces_.size()] = it->first;
      jump_faces_.push_back(it->first);
    }

  for(size_t t = 0; t < fnode_.size(); ++t) {
      fnode_[t] = t;
      groups_[t] << t;
    }

  // test: do not eliminate equations
  ge_ptr.reset(
        new jtf::algorithm::gauss_eliminator<double>(gnode_, gnode_flag_));

  jtf::mesh::one_ring_face_at_point orfap;
  unique_ptr<jtf::mesh::edge2cell_adjacent> ea(
        jtf::mesh::edge2cell_adjacent::create(outside_face_cut));
  if(!ea.get()){
      cerr << "# [error] can not build edge2cell_adjacent." << endl;
      return __LINE__;
    }
  orfap.add_all_faces(outside_face_cut, *ea);

  if(!no_surface){
      for(jtf::mesh::one_ring_face_at_point::p2f_type::const_iterator
          p2f_it = orfap.p2f_.begin(); p2f_it != orfap.p2f_.end(); ++p2f_it){
          const vector<size_t> & linked_faces = p2f_it->second;
          //boost::unordered_set<size_t> linked_face_idx;
          boost::unordered_set<size_t> linked_face_idx_with_rot_idx;
          for(size_t li = 0; li < linked_faces.size(); ++li){
              if(cut_face_pair[linked_faces[li]] == -1){
                  boost::unordered_map<size_t,size_t>::const_iterator bcit =
                      surface_idx_to_rot_idx.find(outside_face_idx_cut[linked_faces[li]]);
                  if(bcit == surface_idx_to_rot_idx.end()){
                      cerr << "# [error] strange can not find surface idx to rot_idx "
                           << outside_face_idx_cut[linked_faces[li]] << endl;
                      return __LINE__;
                    }
                  linked_face_idx_with_rot_idx.insert(bcit->second);
                }
            }
          if(!linked_face_idx_with_rot_idx.empty())
            point_around_faces_with_rot_idx[p2f_it->first] = linked_face_idx_with_rot_idx;
        }
    }

  unique_ptr<jtf::mesh::edge2cell_adjacent> ea_cut(
        jtf::mesh::edge2cell_adjacent::create(outside_face_cut));
  if(!ea_cut.get()){
      cerr << "# [error] can not build edge2cell_adjacent." << endl;
      return __LINE__;
    }

  //for(jtf::mesh::one_ring_tet_at_edge::e2tet_type::const_iterator
  //  cit = ortae_cut.e2t_.begin(); cit != ortae_cut.e2t_.end(); ++cit){

  for(size_t ei = 0; ei < ea_cut->edges_.size(); ++ei){
      const pair<size_t,size_t> & cut_edge = ea_cut->edges_[ei];
      pair<size_t,size_t> orig_edge(cut_tet2tet[cut_edge.first],
          cut_tet2tet[cut_edge.second]);
      jtf::mesh::one_ring_tet_at_edge::e2tet_type::const_iterator ecit
          = ortae.e2t_.find(orig_edge);
      if(ecit == ortae.e2t_.end())
        ecit = ortae.e2t_.find(make_pair(orig_edge.second, orig_edge.first));
      if(ecit == ortae.e2t_.end()){
          cerr << "# [error] strange can not find edge " << orig_edge.first
               << " " << orig_edge.second << endl;
          return __LINE__;
        }
      const vector<size_t> & tet_loop = ecit->second;
      if(!no_surface && ortae.is_inner_edge(tet_loop)){
          // to check whether this edge is near surface
          boost::unordered_map<size_t, boost::unordered_set<size_t> >::const_iterator
              bsscit = point_around_faces_with_rot_idx.find(cut_edge.first);
          if(bsscit != point_around_faces_with_rot_idx.end())
            possible_edge_[cut_edge] = true;
          else{
              bsscit = point_around_faces_with_rot_idx.find(cut_edge.second);
              if(bsscit != point_around_faces_with_rot_idx.end())
                possible_edge_[cut_edge] = true;
              else
                possible_edge_[cut_edge] = false;
            }
        }
      else
        possible_edge_[cut_edge] = false;
    }

  singularity_edges_of_cut_tet_.resize(1);
  singularity_type_of_cut_tet_.resize(1);
  compound_edges_.resize(1);

  return 0;
}

int singularity_graph::insert_transition_group(
    const std::vector<std::pair<size_t,size_t> > & face_groups,
    const matrixst & cut_tet,
    const matrixst & cut_tet2tet,
    const jtf::mesh::face2tet_adjacent & fa_cut,
    const jtf::mesh::one_ring_tet_at_edge & ortae_original,
    const std::vector<size_t> & rot_type,
    const boost::unordered_map<std::pair<size_t,size_t>,size_t> & face_pair_to_rot_idx,
    const boost::unordered_map<std::pair<size_t,size_t>,size_t> & tet_pair2rot_idx,
    const size_t mode )
{
  for(size_t fi = 0; fi < face_groups.size(); ++fi){
      insert_transition(
            face_groups[fi], cut_tet, cut_tet2tet, fa_cut, ortae_original,
            rot_type, face_pair_to_rot_idx, tet_pair2rot_idx, mode);
    }

  return 0;
}

int singularity_graph::insert_transition(
    const std::pair<size_t,size_t> & face,
    const matrixst & cut_tet,
    const matrixst & cut_tet2tet,
    const jtf::mesh::face2tet_adjacent & fa_cut,
    const jtf::mesh::one_ring_tet_at_edge & ortae_original,
    const vector<size_t> & rot_type,
    const boost::unordered_map<pair<size_t,size_t>,size_t> & face_pair_to_rot_idx,
    const boost::unordered_map<pair<size_t,size_t>,size_t> & tet_pair2rot_idx,
    const size_t  mode)
{
  map<pair<size_t,size_t>, vector<size_t> > edges;
  const vector<size_t> & face_first_vec = fa_cut.faces_[face.first];
  const vector<size_t> & face_second_vec = fa_cut.faces_[face.second];

  assert(face_first_vec.size() == face_second_vec.size());

  for(size_t t = 0; t < face_first_vec.size() ;++t){
      pair<size_t,size_t> one_edge(face_first_vec[t],
                                   face_first_vec[(t+1)%face_first_vec.size()]);
      if(one_edge.first > one_edge.second)
        swap(one_edge.first, one_edge.second);
      edges[one_edge].push_back(face.first);
    }

  for(size_t t = 0; t < face_second_vec.size() ;++t){
      pair<size_t,size_t> one_edge(
            face_second_vec[t], face_second_vec[(t+1)%face_second_vec.size()]);
      if(one_edge.first > one_edge.second)
        swap(one_edge.first, one_edge.second);
      edges[one_edge].push_back(face.second);
    }

  // to record the current gnode flag
  boost::dynamic_bitset<> gnode_flag_prev = gnode_flag_;
#if 0
  { // here need a determination, if this face pair has intersection
    // which means f_s and f_t share one or more points,
    // the translation will be determined, for p \in f_s \cap f_t
    // (I - \Pi_st) * f(p) = g_st
    // if row_x of (I - \Pi_st) is (0,0,0), then g_st[x] = 0

    boost::unordered_set<size_t> points;
    points.insert(face_first_vec.begin(), face_first_vec.end());
    points.insert(face_second_vec.begin(), face_second_vec.end());

    if(points.size() < 6){ // there are shared points
        boost::unordered_map<pair<size_t,size_t>,size_t>::const_iterator
            bumpcit = face_pair_to_rot_idx.find(face);

        if(bumpcit == face_pair_to_rot_idx.end()){
            cerr << "# [error] can not find face pair <" << face.first << " "
                 << face.second << endl;
            return __LINE__;
          }
        assert(rot_type[bumpcit->second] != -1);
        const size_t &rot_type_ = rot_type[bumpcit->second];
        matrixd rot_st = trans(type_transition2(rot_type_));
        matrixd IPST = eye<double>(3) - rot_st;

        for(size_t ri = 0; ri < IPST.size(1); ++ri){
            if(norm(IPST(ri,colon())) < 1e-6){ // zero row
                // WARNING: here is an assumption, gnode idx is 0 ~ 3 * n -1,
                // n is the jump face number
                const size_t gnode_idx = get_gnode_idx_from_rot_idx(bumpcit->second,ri);
                if(gnode_idx == -1){
                    cerr << "# [error] can not find gnode_idx of " << bumpcit->second
                         << " " << ri << endl;
                    return __LINE__;
                  }
                jtf::algorithm::equation<double> eq;
                eq.add_expression(jtf::algorithm::make_expression(gnode_idx, 1.0));
                ge_ptr->add_equation(eq);
              }
          }
      }
  }
#endif
  typedef jtf::mesh::one_ring_tet_at_edge::e2tet_type::const_iterator oecit;

  vector<size_t> reordered_tet_loop;
  for(map<pair<size_t,size_t>,vector<size_t> >::const_iterator
      buscit = edges.begin(); buscit != edges.end(); ++buscit){
      const pair<size_t,size_t> &edge = buscit->first;

      const vector<size_t> & edge_begin_face = buscit->second;

      const pair<size_t,size_t> & tet_pair_0 =
          fa_cut.face2tet_[edge_begin_face.front()];

      assert(fa_cut.is_outside_face(tet_pair_0));

      const size_t begin_tet_0 =
          (tet_pair_0.first == -1? tet_pair_0.second:tet_pair_0.first);

      pair<size_t,size_t> orig_edge(
            cut_tet2tet[edge.first], cut_tet2tet[edge.second]);

      oecit it = ortae_original.e2t_.find(orig_edge);
      if(it == ortae_original.e2t_.end()){
          it = ortae_original.e2t_.find(make_pair(orig_edge.second,orig_edge.first));

          if(it == ortae_original.e2t_.end()){
              cerr << "# [error] can not find edge <" << edge.first << ","
                   << edge.second << "> in ortae." << endl;
              return __LINE__;
            }

          //swap(orig_edge.first, orig_edge.second);
          //swap(edge.first, edge.second);
        }
      const vector<size_t> & loop = it->second;

      if(!ortae_original.is_inner_edge(loop)) continue;

      int rtn = cal_singularity_type_at_given_tet(
            loop,begin_tet_0,rot_type, tet_pair2rot_idx,&reordered_tet_loop);
      if(rtn == -1){
          cerr << "# [error] shit happens while calculating edge type." << endl;
          return __LINE__;
        }else if(rtn == -2) {// this loop is not ready
          unready_edges_.push_back(
                std::make_tuple(edge.first, edge.second,begin_tet_0));
          continue;
        }
      else if(is_trivial_type(rtn)){
          add_trivial_edge(reordered_tet_loop, cut_tet, cut_tet2tet,
                           rot_type, fa_cut, tet_pair2rot_idx);

        }else if(is_black_line_new(rtn)){
          //cerr << "# [info] meet compound edge." << endl;

          pair<size_t,size_t> compound_type = get_compound_axis(rtn);
          assert(compound_type.first != -1 && compound_type.second != -1);
          edge_with_type edge0(edge.first, edge.second, compound_type.first);
          edge_with_type edge1(edge.first, edge.second, compound_type.second);

#ifdef debug
          {
            singularity_edges_of_cut_tet_[0].push_back(edge);
            singularity_type_of_cut_tet_[0].push_back(rtn);
          }
#endif
          add_singularity(reordered_tet_loop,cut_tet,cut_tet2tet,
                          tet_pair2rot_idx,rot_type,fa_cut, edge0);
          add_singularity(reordered_tet_loop,cut_tet,cut_tet2tet,
                          tet_pair2rot_idx,rot_type,fa_cut, edge1);
          if(mode == 1) // fast return out
            return __LINE__;

        }else {
          assert(is_regular_type(rtn));
#ifdef debug
          {
            singularity_edges_of_cut_tet_[0].push_back(edge);
            singularity_type_of_cut_tet_[0].push_back(rtn);
          }
#endif
          edge_with_type edge0(edge.first, edge.second, axis_to_around(rtn));
          add_singularity(reordered_tet_loop,cut_tet,cut_tet2tet,
                          tet_pair2rot_idx,rot_type,fa_cut, edge0);
        }
    }

  {
    // after the singularity and trivial edge setting, link the fnodes
    // according to gnode
    gnode_flag_prev ^= gnode_flag_;
    boost::dynamic_bitset<>::size_type modified_gnode_idx =
        gnode_flag_prev.find_first();
    if(modified_gnode_idx != boost::dynamic_bitset<>::npos) {
        while(modified_gnode_idx != boost::dynamic_bitset<>::npos){
            assert(modified_gnode_idx < 3 * jump_faces_.size());
            ///const size_t jump_face_pair_idx = modified_gnode_idx / 3;
            // if gnode \neq 0, ignore it.
            if(fabs(gnode_[static_cast<size_t>(modified_gnode_idx)]) > 1e-6)
              continue;
            connect_fnode_according_to_gnode_at_modified_face_pair(
                  cut_tet2tet, rot_type, tet_pair2rot_idx,fa_cut,
                  get_jump_face_from_gnode_idx(modified_gnode_idx),
                  static_cast<size_t>(modified_gnode_idx));
            modified_gnode_idx = gnode_flag_prev.find_next(modified_gnode_idx);
          }
      }
  }

  if(mode != 2){
      int rtn = check_unready_edges(ortae_original, fa_cut, cut_tet, cut_tet2tet,
                                    rot_type, tet_pair2rot_idx);
      if(mode == 1)
        return rtn;
    }
  return 0;
}

int singularity_graph::insert_surface_transition_group(
    const size_t & rot_type,
    const jtf::mesh::face2tet_adjacent & fa_cut,
    const matrixst & cut_tet2tet,
    const std::vector<size_t> & face_group,
    const std::vector<size_t> & rot_type_vec,
    const boost::unordered_map<std::pair<size_t,size_t>,size_t> & tet_pair2rot_idx)
{
  for(size_t fi = 0; fi < face_group.size(); ++fi){
      insert_surface_transition(rot_type, fa_cut, cut_tet2tet, face_group[fi],
                                rot_type_vec, tet_pair2rot_idx);
    }
  return 0;
}

int singularity_graph::insert_surface_transition(
    const size_t & uvw_type,
    const jtf::mesh::face2tet_adjacent & fa_cut,
    const matrixst & cut_tet2tet,
    const size_t & face_idx,
    const vector<size_t> & rot_type,
    const boost::unordered_map<pair<size_t,size_t>,size_t> & tet_pair2rot_idx)
{
  assert(face_idx < fa_cut.faces_.size());
  const vector<size_t> & face = fa_cut.faces_[face_idx];

  // means not ready
  if(uvw_type == -1)
    return 0;

  assert(uvw_type == 0 || uvw_type == 1 || uvw_type == 2);

  // surface node should be group according to the type,
  // for example, u type face should make their u coordinates of the same group
  for(size_t pi = 1; pi < face.size(); ++pi){
      const size_t gnode_0 = get_node_idx_from_point_idx(face[pi-1],uvw_type);
      const size_t gnode_1 = get_node_idx_from_point_idx(face[pi],uvw_type);
      combine_fnode(gnode_0, gnode_1);
    }

  return 0;
}

int singularity_graph::add_trivial_edge(
    const std::vector<size_t> &tet_loop_vec,
    const matrixst & cut_tet,
    const matrixst & cut_tet2tet,
    const vector<size_t> &rot_type,
    const jtf::mesh::face2tet_adjacent & fa_cut,
    const boost::unordered_map<pair<size_t,size_t>,size_t> & tet_pair2rot_idx)
{
  assert(tet_loop_vec.front() == tet_loop_vec.back());

  vector<jtf::algorithm::equation<double> > eqs;

  add_trivial_edge_equation_for_gnode(tet_loop_vec, tet_pair2rot_idx, rot_type,
                                      cut_tet, cut_tet2tet, fa_cut, eqs);

  //  vector<size_t> rev_tet_loop_vec = tet_loop_vec;
  //  reverse(rev_tet_loop_vec.begin(), rev_tet_loop_vec.end());

  //  if(vec_norm(tet_loop_vec,rev_tet_loop_vec) > 0)
  //    add_trivial_edge_equation_for_gnode(
  //          rev_tet_loop_vec, tet_pair2rot_idx, rot_type,cut_tet,
  //          cut_tet2tet, fa_cut, eqs);

  for(size_t eqi = 0; eqi < eqs.size(); ++eqi){
      ge_ptr->add_equation(eqs[eqi]);
    }

  return 0;
}

int singularity_graph::connect_fnode_according_to_gnode_at_modified_face_pair(
    const matrixst & cut_tet2tet,
    const std::vector<size_t> & rot_type,
    const boost::unordered_map<std::pair<size_t,size_t>,size_t> & rot_type_map_idx,
    const jtf::mesh::face2tet_adjacent &fa_cut,
    const std::pair<size_t,size_t> & face_pair,
    const size_t & gnode_idx)
{
  assert(fabs(gnode_[gnode_idx]) < 1e-6);
  assert(gnode_flag_[gnode_idx]);

  const pair<size_t,size_t> &tet_pair_0 = fa_cut.face2tet_[face_pair.first];

  const pair<size_t,size_t> &tet_pair_1 = fa_cut.face2tet_[face_pair.second];

  const size_t & tet0 =
      (tet_pair_0.first == -1? tet_pair_0.second:tet_pair_0.first);
  const size_t & tet1 =
      (tet_pair_1.first == -1? tet_pair_1.second:tet_pair_1.first);

  boost::unordered_map<pair<size_t,size_t>,size_t>::const_iterator bumcit =
      rot_type_map_idx.find(make_pair(tet0, tet1));
  if(bumcit == rot_type_map_idx.end()){
      cerr << "# [error] strange can not find rot type between tet pair <"
           << tet0 << "," << tet1 << ">." << endl;
      return __LINE__;
    }

  assert(rot_type[bumcit->second] != -1); // this rot type shoud be calculated

  const matrixd rot = trans(type_transition2(rot_type[bumcit->second]));
  const size_t uvw_idx_t = gnode_idx % 3;

  //  assert(rot_type[bumcit->second] <= TRIVIAL_TYPE);
  // f_t(uvw_idx) should be grouped with  f_s(x)
  // where rot * f_s(x) = f_t(uvw_idx)
  matrixd sf = ones<double>(3,1);
  sf[1] = 2; sf[2] = 3;
  sf = temp(rot * sf);
  int uvw_idx_s = number_rounding(sf[uvw_idx_t]);//static_cast<size_t>(fabs(sf[uvw_idx_t]) + 0.5);

  if(uvw_idx_s < 0) {
      //cerr << "# [info] different sign should not be glued " << endl;
      return 0;
    }

  if(is_trivial_type(rot_type[bumcit->second])){
      if(uvw_idx_t != uvw_idx_s-1){
          cerr << "# [error] strange: two uvw_idx should be the same." << endl;
          return __LINE__;
        }
    }

  --uvw_idx_s;

  const vector<size_t> & face_s_vec = fa_cut.faces_[face_pair.first];
  const vector<size_t> & face_t_vec = fa_cut.faces_[face_pair.second];
  assert(face_s_vec.size() == face_t_vec.size());

  for(size_t pi = 0; pi < face_s_vec.size(); ++pi){
      const size_t & point_i = face_s_vec[pi];
      for(size_t pj = 0; pj < face_t_vec.size(); ++pj){
          const size_t & point_j = face_t_vec[pj];
          if(cut_tet2tet[point_i] == cut_tet2tet[point_j]){
              combine_fnode( get_node_idx_from_point_idx(point_i, uvw_idx_s),
                             get_node_idx_from_point_idx(point_j, uvw_idx_t));
              break;
            }
        }
    }
  return 0;
}


int singularity_graph::convert_g_eq_to_f_eq(
    const vector<size_t> & rot_type,
    const matrixst & cut_tet2tet,
    const jtf::mesh::face2tet_adjacent & fa_cut,
    const boost::unordered_map<std::pair<size_t,size_t>,size_t> & jump_face_to_rot_idx,
    std::shared_ptr<jtf::algorithm::gauss_eliminator<double> > &fe,
    std::vector<double> &fnode_vec,
    boost::dynamic_bitset<> & fnode_flag)const
{
  fnode_vec.resize(fnode_.size());
  fnode_flag.resize(fnode_.size());
  fe.reset(new jtf::algorithm::gauss_eliminator<double>(fnode_vec,fnode_flag,false));

  typedef list<jtf::algorithm::equation<double> > equation_list_type;
  const equation_list_type &eqs = ge_ptr->get_equations();

  // for each equation of gnode, replace g_st with f_t(p) - \Pi_st*f_s(p)
  matrixd rot_st = eye<double>(3);
  for(equation_list_type::const_iterator eqcit = eqs.begin(); eqcit != eqs.end();
      ++eqcit){
      const jtf::algorithm::equation<double> & eq = *eqcit;
      if(eq.state() == 0) continue; // empty equation
      if(eq.state() == -1) {
          cerr << "# [error] conflict equations." << endl;
          continue;
        }

      vector<jtf::algorithm::equation<double> > eq_fnode(3);

      vector<vector<size_t> > face_st(2);
      for(jtf::algorithm::equation<double>::eq_const_iterator ecit = eq.begin();
          ecit != eq.end(); ++ecit){
          const jtf::algorithm::expression<double> & exp = *ecit;
          const size_t & gnode_idx = exp.index;
          const double & gnode_cof = exp.coefficient;
          const size_t rot_idx = get_rot_idx_from_gnode_idx(gnode_idx).first;
          const size_t uvw = get_rot_idx_from_gnode_idx(gnode_idx).second;

          rot_st = trans(type_transition2(rot_type[rot_idx]));

          //assert(jump_face_idx < jump_faces_.size());
          const pair<size_t,size_t>  & jump_face_pair = jump_faces_[gnode_idx/3];
          const vector<size_t> &face_s = fa_cut.faces_[jump_face_pair.first];
          const vector<size_t> &face_t = fa_cut.faces_[jump_face_pair.second];
          assert(face_s.size() == face_t.size());
          for(size_t pis = 0; pis < face_s.size(); ++pis){
              for(size_t pit = 0; pit < face_t.size(); ++pit){
                  if(cut_tet2tet[face_s[pis]]== cut_tet2tet[face_t[pit]]){
                      face_st[0].push_back(face_s[pis]);
                      face_st[1].push_back(face_t[pit]);
                      break;
                    }
                }
            }
          assert(face_st[0].size() == face_st[1].size());
          matrixd face_s_p_mat = ones<double>(3,1);
          face_s_p_mat[0] = get_node_idx_from_point_idx(face_st[0][0],0);
          face_s_p_mat[1] = get_node_idx_from_point_idx(face_st[0][0],1);
          face_s_p_mat[2] = get_node_idx_from_point_idx(face_st[0][0],2);

          face_s_p_mat = temp(rot_st * face_s_p_mat);
          for(size_t ai = 0; ai < 3; ++ai){
              eq_fnode[ai].add_expression(
                    jtf::algorithm::make_expression(
                      get_node_idx_from_point_idx(face_st[1][0],ai),
                  gnode_cof));
              eq_fnode[ai].add_expression(
                    jtf::algorithm::make_expression(
                      //abs(number_rounding(face_s_p_mat[ai])) - 1,
                      static_cast<size_t>(fabs(face_s_p_mat[ai])),
                      -1.0*gnode_cof*jtf::math::get_sign(face_s_p_mat[ai])));
            }
        }
      fe->add_equation(eq_fnode[0]);
      fe->add_equation(eq_fnode[1]);
      fe->add_equation(eq_fnode[2]);
    }

  // for each gnode: g_st = f_t(p)  - \Pi_st*f_s(p)
  // add equations: f_t(p) - \Pi_st*f_s(p) = f_t(q) - \Pi_st*f_s(q)
  // add equations: f_t(q) - \Pi_st*f_s(q) = f_t(r) - \Pi_st*f_s(r)
  vector<vector<size_t> > face_st_2(2);
  vector<matrixd> s_pqr(3);
  for(boost::unordered_map<pair<size_t,size_t>,size_t>::const_iterator bucit
      = jump_face_to_rot_idx.begin(); bucit != jump_face_to_rot_idx.end(); ++bucit){

      //for(size_t ji = 0; ji < jump_faces_.size(); ++ji){
      const pair<size_t,size_t> & jump_face_st = bucit->first;

      const vector<size_t> &face_s = fa_cut.faces_[jump_face_st.first];
      const vector<size_t> &face_t = fa_cut.faces_[jump_face_st.second];
      assert(face_s.size() == face_t.size());
      for(size_t pis = 0; pis < face_s.size(); ++pis){
          for(size_t pit = 0; pit < face_t.size(); ++pit){
              if(cut_tet2tet[face_s[pis]]== cut_tet2tet[face_t[pit]]){
                  face_st_2[0].push_back(face_s[pis]);
                  face_st_2[1].push_back(face_t[pit]);
                  break;
                }
            }
        }

      rot_st = trans(type_transition2(rot_type[bucit->second]));
      for(size_t pi = 0; pi < 3; ++pi){
          s_pqr[pi] = zeros<double>(3,1);
          for(size_t ai = 0; ai < 3; ++ai){
              s_pqr[pi][ai] = get_node_idx_from_point_idx(face_st_2[0][pi],ai);
            }
          s_pqr[pi] = temp(rot_st * s_pqr[pi]);
        }

      for(size_t pi = 0; pi < 2; ++pi){
          vector<jtf::algorithm::equation<double> > eqs0(3);
          for(size_t ai = 0; ai < 3; ++ai){
              eqs0[ai].add_expression(
                    jtf::algorithm::make_expression(
                      get_node_idx_from_point_idx(face_st_2[1][pi],ai), 1.0));
              eqs0[ai].add_expression(
                    jtf::algorithm::make_expression(
                      static_cast<size_t>(fabs(s_pqr[pi][ai])),
                      -1.0* jtf::math::get_sign(s_pqr[pi][ai])));

              eqs0[ai].add_expression(
                    jtf::algorithm::make_expression(
                      get_node_idx_from_point_idx(face_st_2[1][pi+1],ai), -1.0));

              eqs0[ai].add_expression(
                    jtf::algorithm::make_expression(
                      static_cast<size_t>(fabs(s_pqr[pi+1][ai])),
                    1.0*jtf::math::get_sign(s_pqr[pi+1][ai])));
              fe->add_equation(eqs0[ai]);
            }
        }
    }

  return 0;
}

int singularity_graph::add_singularity(
    const std::vector<size_t> &tet_loop_vec,
    const matrixst & cut_tet,
    const matrixst & cut_tet2tet,
    const boost::unordered_map<pair<size_t,size_t>,size_t> & rot_type_map_idx,
    const vector<size_t> &rot_type,
    const jtf::mesh::face2tet_adjacent & fa_cut,
    const edge_with_type & se)
{
  const boost::unordered_map<size_t,size_t>::const_iterator msscit_first =
      outside_point_idx_.find(se.point0());
  const boost::unordered_map<size_t,size_t>::const_iterator msscit_second =
      outside_point_idx_.find(se.point1());

  if(msscit_first == outside_point_idx_.end()){
      cerr << "# [error] strange, can not find edge point "<< se.point0() << endl;
      return __LINE__;
    }

  if(msscit_second == outside_point_idx_.end()){
      cerr << "# [error] strange, can not find edge point "<< se.point1() << endl;
      return __LINE__;
    }

  assert(se.type() < 3);

  vector<size_t> type_equal;
  type_equal << (se.type()+1)%3 << (se.type()+2)%3;

  // add euqal constraints
  for(size_t t = 0; t < type_equal.size(); ++t){
      const size_t gnode_0 = get_node_idx_from_point_idx(msscit_first->first,type_equal[t]);
      const size_t gnode_1 = get_node_idx_from_point_idx(msscit_second->first,type_equal[t]);
      combine_fnode(gnode_0, gnode_1);

      assert(fnode_[gnode_1] == fnode_[gnode_0]);

      //    jtf::algorithm::equation<double> eq;
      //    eq.add_expression(jtf::algorithm::make_expression(gnode_0,1.0));
      //    eq.add_expression(jtf::algorithm::make_expression(gnode_1,-1.0));
      //    ge_ptr->add_equation(eq);
    }

  jtf::algorithm::equation<double> eq;
  add_singularity_equation_for_gnode(
        tet_loop_vec, rot_type_map_idx,rot_type, cut_tet,cut_tet2tet,
        fa_cut, se.type(),eq);

  ge_ptr->add_equation(eq);

  return 0;
}

int singularity_graph::add_trivial_edge_equation_for_gnode(
    const std::vector<size_t> & tet_loop_vec,
    const boost::unordered_map<pair<size_t,size_t>,size_t> & tet_pair2rot_idx,
    const std::vector<size_t> & rot_type,
    const matrixst & cut_tet,
    const matrixst & cut_tet2tet,
    const jtf::mesh::face2tet_adjacent &fa_cut,
    vector<jtf::algorithm::equation<double> > &eqs)
{
  vector<size_t> rot_type_of_loop;
  rot_type_of_loop.reserve(tet_loop_vec.size());

  typedef boost::unordered_map<pair<size_t,size_t>,size_t>::const_iterator
      bupcit;

  for(size_t fi = 0; fi < tet_loop_vec.size() - 1; ++fi){
      bupcit it =
          tet_pair2rot_idx.find(make_pair(tet_loop_vec[fi], tet_loop_vec[fi+1]));
      if(it != tet_pair2rot_idx.end()){
          if(rot_type[it->second] == -1){
              cerr << "# [error] strange, this rot type should not be unkonwn." << endl;
              return __LINE__;
            }
          rot_type_of_loop.push_back(rot_type[it->second]);
        }else{
          it = tet_pair2rot_idx.find(make_pair(tet_loop_vec[fi+1], tet_loop_vec[fi]));
          if(it != tet_pair2rot_idx.end())
            rot_type_of_loop.push_back(get_trans_type(rot_type[it->second]));
          else{
              rot_type_of_loop.push_back(TRIVIAL_TYPE); //it's glued inside the cut mesh
            }
        }
    }

  {
    matrixd rot_test = eye<double>(3);
    for(size_t ri = 0; ri < rot_type_of_loop.size(); ++ri){
        rot_test = temp(rot_test * type_transition2(rot_type_of_loop[ri]));
      }
    if(norm(rot_test - eye<double>(3)) > 1e-6){
        cerr << "# [error] trivial edge has non-trivial type "
             << type_transition1(rot_test) << endl;
        return __LINE__;
      }
  }

  vector<jtf::algorithm::equation<double> > eq3(3);

  vector<size_t> common_face(3,0);
  for(size_t ti = 0; ti < tet_loop_vec.size() - 1; ++ti){
      // if can not find a common face of these two tet in cut_tet meshn
      // that means there need a transition, or the transition is zeros
      if(jtf::mesh::find_common_face(cut_tet(colon(), tet_loop_vec[ti]),
                                     cut_tet(colon(),tet_loop_vec[ti+1]), &common_face[0]))
        {// can not find common face, which means there must be a transistion
          vector<size_t> face0, face1;
          for(size_t pi = 0; pi < cut_tet.size(1); ++pi){
              const size_t & point_i = cut_tet(pi, tet_loop_vec[ti]);
              for(size_t pj = 0; pj < cut_tet.size(1); ++pj){
                  const size_t & point_j = cut_tet(pj, tet_loop_vec[ti+1]);
                  if(cut_tet2tet[point_i] == cut_tet2tet[point_j]){
                      face0.push_back(point_i);
                      face1.push_back(point_j);
                      break;
                    }
                }
            }
          if(face0.size() != 3 || face1.size() != 3){
              cerr << "# [error] can not find jump face pair of tets <" << tet_loop_vec[ti]
                      << "," << tet_loop_vec[ti+1] << ">." << endl;
              return __LINE__;
            }
          const size_t face_idx0 = fa_cut.get_face_idx(&face0[0]);
          const size_t face_idx1 = fa_cut.get_face_idx(&face1[0]);
          if(face_idx0 == -1 || face_idx1 == -1){
              cerr << "# [error] can not find faces " << endl;
              return __LINE__;
            }

          double direction = 1.0;
          bupcit fit =jump_face_to_gnode_idx_.find(make_pair(face_idx0, face_idx1));

          if(fit == jump_face_to_gnode_idx_.end()){
              fit = jump_face_to_gnode_idx_.find(make_pair(face_idx1, face_idx0));

              if(fit != jump_face_to_gnode_idx_.end())
                direction = -1.0;
              else{
                  cerr << "# [error] can not find face pair <" << face_idx0 << ","
                       << face_idx1 << "> in jump_face_to_gnode_idx." << endl;
                  return __LINE__;
                }
            }
          // each gnode has u,v,w components
          matrixd node_idx = ones<double>(3,1);//+ 3 * fit->second;
          node_idx[1] += 1;
          node_idx[2] += 2;

          if(direction < 0){
              // if direction is -1, g_ts = - R_ts * g_st
              node_idx = -1.0 * temp(trans(type_transition2(rot_type_of_loop[ti])) * node_idx);
            }

          for(size_t ri = ti +1 ; ri < rot_type_of_loop.size(); ++ri)
            node_idx = temp(trans(type_transition2(rot_type_of_loop[ri])) * node_idx);

          for(size_t ni = 0; ni < node_idx.size(); ++ni){
              node_idx[ni] = jtf::math::get_sign(node_idx[ni]) *
                  (fabs(node_idx[ni]) + 3 * fit->second-1);
            }

          for(size_t ui = 0; ui < node_idx.size(); ++ui){
              eq3[ui].add_expression(jtf::algorithm::make_expression(
                                       abs(number_rounding(node_idx[ui])),
                                       //static_cast<size_t>(fabs(node_idx[ui])+0.5),
                                       jtf::math::get_sign(node_idx[ui])));
            }
        }
    }

  for(size_t eqi = 0; eqi < eq3.size(); ++eqi)
    eqs.push_back(eq3[eqi]);

  return 0;
}

int singularity_graph::add_singularity_equation_for_gnode(
    const std::vector<size_t> & tet_loop_vec,
    const boost::unordered_map<pair<size_t,size_t>,size_t> & rot_type_map_idx,
    const std::vector<size_t> & rot_type,
    const matrixst & cut_tet,
    const matrixst & cut_tet2tet,
    const jtf::mesh::face2tet_adjacent &fa_cut,
    const size_t & singularity_axis_type,
    jtf::algorithm::equation<double> &eq)
{
  typedef boost::unordered_map<pair<size_t,size_t>, size_t>
      ::const_iterator bupcit;

  vector<size_t> rot_type_of_loop;
  rot_type_of_loop.reserve(tet_loop_vec.size());

  for(size_t fi = 0; fi < tet_loop_vec.size() - 1; ++fi){
      bupcit it = rot_type_map_idx.find(
            make_pair(tet_loop_vec[fi], tet_loop_vec[fi+1]));
      if(it != rot_type_map_idx.end()){
          rot_type_of_loop.push_back(rot_type[it->second]);
        }else{
          it = rot_type_map_idx.find(make_pair(tet_loop_vec[fi+1], tet_loop_vec[fi]));
          if(it != rot_type_map_idx.end())
            rot_type_of_loop.push_back(
                  type_transition1(trans(type_transition2(rot_type[it->second]))));
          else
            rot_type_of_loop.push_back(TRIVIAL_TYPE);
        }
    }

  vector<size_t> common_face(3,0);

  for(size_t ti = 0; ti < tet_loop_vec.size() - 1; ++ti){
      // if can not find a common face of these two tet in cut_tet meshn
      // that means there need a transition, or the transition is zeros
      if(jtf::mesh::find_common_face(cut_tet(colon(), tet_loop_vec[ti]),
                                     cut_tet(colon(),tet_loop_vec[ti+1]),
                                     &common_face[0]) == 1)
        {// can not find common face, which means there must be a transistion
          vector<size_t> face0, face1;
          for(size_t pi = 0; pi < cut_tet.size(1); ++pi){
              const size_t & point_i = cut_tet(pi, tet_loop_vec[ti]);
              for(size_t pj = 0; pj < cut_tet.size(1); ++pj){
                  const size_t & point_j = cut_tet(pj, tet_loop_vec[ti+1]);
                  if(cut_tet2tet[point_i] == cut_tet2tet[point_j]){
                      face0.push_back(point_i);
                      face1.push_back(point_j);
                      break;
                    }
                }
            }
          if(face0.size() != 3 || face1.size() != 3){
              cerr << "# [error] can not find jump face pair of tets <" << tet_loop_vec[ti]
                      << "," << tet_loop_vec[ti+1] << ">." << endl;
              return __LINE__;
            }
          const size_t face_idx0= fa_cut.get_face_idx(&face0[0]);
          const size_t face_idx1= fa_cut.get_face_idx(&face1[0]);
          if(face_idx0 == -1 || face_idx1 == -1){
              cerr << "# [error] can not find faces " << endl;
              return __LINE__;
            }

          double direction = 1.0;
          bupcit fit =
              jump_face_to_gnode_idx_.find(make_pair(face_idx0, face_idx1));

          if(fit == jump_face_to_gnode_idx_.end()){
              fit = jump_face_to_gnode_idx_.find(make_pair(face_idx1, face_idx0));

              if(fit != jump_face_to_gnode_idx_.end())
                direction = -1.0;
              else{
                  cerr << "# [error] can not find face pair <" << face_idx0 << ","
                       << face_idx1 << "> in jump_face_to_gnode_idx." << endl;
                  return __LINE__;
                }
            }

          // each gnode has u,v,w components
          matrixd node_idx = ones<double>(3,1);
          node_idx += 3 * fit->second;
          node_idx[1] += 1;
          node_idx[2] += 2;

          if(direction < 0){
              // if direction is -1, g_ts = - R_ts * g_st
              node_idx = -1*temp(trans(type_transition2(rot_type_of_loop[ti])) * node_idx);
            }

          for(size_t ri = ti +1 ; ri < rot_type_of_loop.size(); ++ri)
            node_idx = temp(trans(type_transition2(rot_type_of_loop[ri])) * node_idx);

          eq.add_expression(
                jtf::algorithm::make_expression(
                  abs(number_rounding(node_idx[singularity_axis_type]))-1,
                  //static_cast<size_t>(fabs(node_idx[singularity_axis_type])+0.5),
                  jtf::math::get_sign(node_idx[singularity_axis_type])));
        }
    }

  return 0;
}

std::pair<size_t,size_t> singularity_graph::get_point_idx_from_node_idx(
    const size_t &node_idx) const
{
  if(node_idx >= fnode_.size()) return make_pair(-1,-1);

  const size_t point_idx = node_idx / 3;
  const size_t uvw_type = node_idx % 3;

  return make_pair(outside_points_[point_idx], uvw_type);
}

size_t singularity_graph::get_node_idx_from_point_idx(
    const size_t & point_idx, const size_t & type)const
{
  assert(type < 3);
  boost::unordered_map<size_t,size_t>::const_iterator bumcit =
      outside_point_idx_.find(point_idx);
  if(bumcit == outside_point_idx_.end()){
      cerr << "# [error] can not find point " << point_idx << endl;
      return -1;
    }

  return 3 * bumcit->second + type;
}

int singularity_graph::get_step_state(step_state &ss) const
{
  ss.gnode_ = gnode_;
  ss.gnode_flag_ = gnode_flag_;
  ss.fnode_ = fnode_;
  ss.es_ = ge_ptr->get_equations();
  ss.unready_edges_ = unready_edges_;
  return 0;
}

int singularity_graph::set_step_state(step_state &ss)
{
  assert(gnode_.size() == ss.gnode_.size());
  assert(fnode_.size() == ss.fnode_.size());
  assert(gnode_flag_.size() == ss.gnode_flag_.size());
  copy(ss.gnode_.begin(), ss.gnode_.end(), gnode_.begin());
  ss.gnode_flag_ = gnode_flag_;
  copy(ss.fnode_.begin(), ss.fnode_.end(), fnode_.begin());

  unready_edges_ = ss.unready_edges_;

  groups_.clear();
  groups_.resize(fnode_.size());
  for(size_t fi = 0; fi < fnode_.size(); ++fi){
      //groups_[fnode_[fi]].insert(fi);
      groups_[fnode_[fi]] << fi;
    }

  ge_ptr->clear();
  for(list<jtf::algorithm::equation<double> >::const_iterator eqcit =
      ss.es_.begin(); eqcit != ss.es_.end(); ++eqcit){
      if(eqcit->state() == 2)
        ge_ptr->add_equation(*eqcit);
    }
  //ss.es_ = ge_ptr->get_equations();

  return 0;
}

int singularity_graph::save_state(const std::string & file_name )const
{
  std::ofstream ofs(file_name.c_str(), ios::binary);
  if(ofs.fail()){
      cerr << "# [error] can not open state file." << endl;
      return __LINE__;
    }
  ge_ptr->save_equations(ofs);
  const size_t fnode_num = fnode_.size();
  ofs.write((char*)&fnode_num, sizeof(size_t));
  ofs.write((char*)&fnode_[0], fnode_num * sizeof(size_t));

  const size_t unready_edges_num = unready_edges_.size();
  ofs.write((char*)&unready_edges_num, sizeof(size_t));
  for(std::list<std::tuple<size_t,size_t,size_t> >::const_iterator
      cit = unready_edges_.begin(); cit != unready_edges_.end(); ++cit){
      const std::tuple<size_t,size_t,size_t> & ed = *cit;
      ofs.write((char*)&ed,
                sizeof(list<std::tuple<size_t,size_t,size_t> >::value_type));
    }
  return 0;
}

int singularity_graph::load_state(const string &file_name)
{
  std::ifstream ifs(file_name.c_str(), ios::binary);
  if(ifs.fail()){
      cerr << "# [error] can not open state file." << endl;
      return __LINE__;
    }
  ge_ptr->load_equations(ifs);
  size_t fnode_num = 0;
  ifs.read((char*)&fnode_num, sizeof(size_t));
  if(fnode_num != fnode_.size()){
      cerr << "# [error] wrong state file." << endl;
      return __LINE__;
    }
  ifs.read((char*)&fnode_[0], fnode_num * sizeof(size_t));
  groups_.clear();
  groups_.resize(fnode_num);
  for(size_t fi = 0; fi < fnode_num; ++fi){
      //groups_[fnode_[fi]].insert(fi);
      groups_[fnode_[fi]] << fi;
    }

  size_t unready_edges_num = 0;
  ifs.read((char*)&unready_edges_num, sizeof(size_t));
  std::tuple<size_t,size_t,size_t> ue;
  unready_edges_.clear();
  for(size_t uei = 0; uei < unready_edges_num; ++uei){
      ifs.read((char*)&ue,
               sizeof(list<std::tuple<size_t,size_t,size_t> >::value_type));
      unready_edges_.push_back(ue);
    }
  return 0;
}

int singularity_graph::save_state_mem(state_each_step & ses)const
{
  ge_ptr->save_equations_mem(ses.nodes_, ses.node_flag_, ses.eq_vec_);

  ses.fnode_ = fnode_;
  ses.unready_edges_ = unready_edges_;

  return 0;
}

int singularity_graph::clear_state_mem()
{
  vector<double> nodes_(gnode_.size(),0);
  boost::dynamic_bitset<> node_flag_(gnode_.size());
  std::vector<jtf::algorithm::equation<double> > eq_vec_;

  ge_ptr->load_equations_mem(nodes_, node_flag_, eq_vec_);

  unready_edges_.clear();
  groups_.clear();
  groups_.resize(fnode_.size());
  for(size_t i = 0; i < fnode_.size(); ++i) {
      fnode_[i] = i;
      //groups_[i].insert(i);
      groups_[i] << i;
    }

  return 0;
}

int singularity_graph::load_state_mem(const state_each_step & ses)
{
  ge_ptr->load_equations_mem(ses.nodes_, ses.node_flag_, ses.eq_vec_);
  fnode_ = ses.fnode_;
  unready_edges_ = ses.unready_edges_;

  groups_.clear();
  groups_.resize(fnode_.size());
  for(size_t fi = 0; fi < fnode_.size(); ++fi){
      groups_[fnode_[fi]] << fi;
    }
  return 0;
}

bool singularity_graph::is_group_info_broken()const
{
  group<size_t> one_group;
  for(size_t gi = 0; gi < groups_.size(); ++gi){
      if(groups_[gi].empty()){
          if(fnode_[gi] == gi)
            return true;
        }
      one_group << groups_[gi];
    }

  if(one_group.size() != fnode_.size())
    return true;

  for(size_t fi = 0; fi < fnode_.size(); ++fi){
      const size_t &gi = fnode_[fi];
      boost::unordered_set<size_t>::const_iterator bucit = groups_[gi].find(fi);
      if(bucit == groups_[gi].end()){
          return true;
        }
    }
  return false;
}

int singularity_graph::check_unready_edges(
    const jtf::mesh::one_ring_tet_at_edge & ortae_original,
    const jtf::mesh::face2tet_adjacent & fa_cut,
    const matrixst & cut_tet,
    const matrixst & cut_tet2tet,
    const vector<size_t> & rot_type,
    const boost::unordered_map<pair<size_t,size_t>,size_t> & tet_pari2rot_idx,
    const size_t mode)
{
  typedef std::tuple<size_t,size_t,size_t> unready_edges_type;
  typedef jtf::mesh::one_ring_tet_at_edge::e2tet_type::const_iterator oecit;

  boost::dynamic_bitset<> node_flag_prev = gnode_flag_;
  vector<size_t> reordered_tet_loop;
  for(list<unready_edges_type>::iterator lucit = unready_edges_.begin();
      lucit != unready_edges_.end();){
      const unready_edges_type & uet = *lucit;
      pair<size_t,size_t> orig_edge(cut_tet2tet[get<0>(uet)],
          cut_tet2tet[get<1>(uet)]);
      oecit it = ortae_original.e2t_.find(orig_edge);
      if(it == ortae_original.e2t_.end())
        it = ortae_original.e2t_.find(make_pair(orig_edge.second, orig_edge.first));
      assert(it != ortae_original.e2t_.end());
      const vector<size_t> &loop = it->second;
      const size_t &begin_tet = get<2>(uet);
      int rtn =
          cal_singularity_type_at_given_tet(
            loop,begin_tet,rot_type, tet_pari2rot_idx, &reordered_tet_loop);
      if(rtn == -2) {
          ++lucit;
        }else {
          if(is_trivial_type(rtn)){
              add_trivial_edge(reordered_tet_loop, cut_tet, cut_tet2tet,
                               rot_type, fa_cut, tet_pari2rot_idx);
            }else if(is_black_line_new(rtn)){

              pair<size_t,size_t> compound_type = get_compound_axis(rtn);
              assert(compound_type.first != -1 && compound_type.second != -1);
              edge_with_type edge0(get<0>(uet), get<1>(uet), compound_type.first);
              edge_with_type edge1(get<0>(uet), get<1>(uet), compound_type.second);

              add_singularity(reordered_tet_loop,cut_tet,cut_tet2tet,
                              tet_pari2rot_idx,rot_type,fa_cut, edge0);
              add_singularity(reordered_tet_loop,cut_tet,cut_tet2tet,
                              tet_pari2rot_idx,rot_type,fa_cut, edge1);
              if(mode == 1)
                return __LINE__;
            }else {
              assert(is_regular_type(rtn));
              edge_with_type edge0(get<0>(uet), get<1>(uet), axis_to_around(rtn));
              add_singularity(reordered_tet_loop,cut_tet,cut_tet2tet,
                              tet_pari2rot_idx,rot_type,fa_cut, edge0);
            }

          unready_edges_.erase(lucit++);
        }
    }

  {
    // after the singularity and trivial edge setting, link the fnodes
    // according to gnode
    node_flag_prev ^= gnode_flag_;
    boost::dynamic_bitset<>::size_type modified_gnode_idx =
        node_flag_prev.find_first();
    if(modified_gnode_idx != boost::dynamic_bitset<>::npos) {
        while(modified_gnode_idx != boost::dynamic_bitset<>::npos){
            assert(modified_gnode_idx < 3 * jump_faces_.size());
            //const size_t jump_face_pair_idx = modified_gnode_idx / 3;
            // if gnode \neq 0, ignore it.
            if(fabs(gnode_[static_cast<size_t>(modified_gnode_idx)]) > 1e-6)
              continue;
            connect_fnode_according_to_gnode_at_modified_face_pair(
                  cut_tet2tet, rot_type, tet_pari2rot_idx,fa_cut,
                  get_jump_face_from_gnode_idx(modified_gnode_idx),
                  static_cast<size_t>(modified_gnode_idx));
            modified_gnode_idx = node_flag_prev.find_next(modified_gnode_idx);
          }
      }
  }
  return 0;
}

pair<size_t,size_t> singularity_graph::get_rot_idx_from_gnode_idx(
    const size_t & g_idx)const
{
  assert(g_idx < gnode_.size());
  const size_t rot_idx = gnode2rot_idx_[g_idx / 3];
  const size_t uvw = g_idx % 3;
  return make_pair(rot_idx,uvw);
}

const std::pair<size_t,size_t>& singularity_graph::get_jump_face_from_gnode_idx(
    const size_t gnode_idx)const
{
  assert(gnode_idx < gnode_.size());
  const size_t g_idx = gnode_idx/3;
  return jump_faces_[g_idx];
}

size_t singularity_graph::get_gnode_idx_from_rot_idx( const size_t & rot_idx,
                                                      const size_t & type)const
{
  boost::unordered_map<size_t,size_t>::const_iterator it =
      rot_idx2g_idx_.find(rot_idx);
  if(it != rot_idx2g_idx_.end()){
      return 3 * it->second + type;
    }
  return -1;
}

int singularity_graph::combine_fnode(const size_t & fnode_idx_0,
                                     const size_t & fnode_idx_1)
{
  if(fnode_[fnode_idx_0] == fnode_[fnode_idx_1]) return 0;

  if(fnode_[fnode_idx_0]  != fnode_[fnode_idx_1]){
      assert(!groups_[fnode_[fnode_idx_1]].empty());

      groups_[fnode_[fnode_idx_0]] << groups_[fnode_[fnode_idx_1]];
      const size_t group_need_to_be_clean = fnode_[fnode_idx_1];
      for(boost::unordered_set<size_t>::const_iterator cit =
          groups_[fnode_[fnode_idx_1]].begin();
          cit != groups_[fnode_[fnode_idx_1]].end(); ++cit){
          assert(fnode_[*cit] != fnode_[fnode_idx_0]);
          fnode_[*cit] = fnode_[fnode_idx_0];
        }
      groups_[group_need_to_be_clean].clear();
    }

  assert(!is_group_info_broken());
  assert(fnode_[fnode_idx_0] == fnode_[fnode_idx_1]);
  return 0;
}

size_t singularity_graph::extract_free_axis_of_edge(const size_t &p0,
                                                    const size_t &p1)
{
  vector<size_t> fixed_components;
  fixed_components.reserve(3);
  for(size_t ai = 0; ai < 3; ++ai){
      const size_t fnode_0 = get_node_idx_from_point_idx(p0,ai);
      const size_t fnode_1 = get_node_idx_from_point_idx(p1,ai);
      assert(fnode_0 != -1 && fnode_1 != -1);
      if(fnode_[fnode_0] == fnode_[fnode_1]){
          fixed_components.push_back(ai);
        }
    }

  if(fixed_components.size() != 2){
      cerr << "# [error] this two points have " << 3 - fixed_components.size()
           << " free degrees." << endl;
      return -1;
    }

  return (3 - std::accumulate(fixed_components.begin(), fixed_components.end(),0));
  //return -1;
}

int singularity_graph::extract_directed_chain_from_node_edges(
    const boost::unordered_set<std::pair<size_t,size_t> > & node_edges,
    std::vector<std::deque<std::pair<size_t,size_t> > > & loops,
    const matrixst &  cut_tet2tet,
    const vector<size_t> & rot_type,
    const jtf::mesh::one_ring_tet_at_edge & ortae,
    const boost::unordered_map<size_t,size_t >& surface_face2rot_idx,
    const size_t mode)
{
  // I want to build connections between node edges, here nodes of each edge
  // are components of points, some node can be grouped if and only if three
  // components of the two points are the same.

  typedef boost::unordered_set<pair<size_t,size_t> > node_edge_type;
  typedef vector<size_t> uvw_type;
  typedef boost::unordered_map<uvw_type, boost::unordered_set<size_t> >
      group2node_type;
  group2node_type group2node;

  // first extract chains from these node_edges
  vector<pair<size_t,size_t> > candidate_edges(node_edges.begin(),
                                               node_edges.end());
  vector<deque<pair<size_t,size_t> > > chains;
  jtf::util::extract_chain_from_edges(candidate_edges, chains);
  matrixst uvw = zeros<double>(3,2);
  vector<deque<pair<size_t,size_t> > > trivial_loops;
  for(size_t ci = 0; ci < chains.size(); ++ci){
      const deque<pair<size_t,size_t> > & deq = chains[ci];
      if(deq.front().first == deq.back().second){
          if(mode == 1)
            return __LINE__;
          else
            trivial_loops.push_back(deq);
        }else{
          pair<size_t,size_t> end_nodes(deq.front().first, deq.back().second);
          end_nodes.first = 3 * (end_nodes.first/3);
          end_nodes.second = 3 * (end_nodes.second/3);
          uvw = zeros<double>(3,2);
          for(size_t pi = 0; pi < 3; ++pi){
              uvw(pi,0) = fnode_[end_nodes.second + pi];
              uvw(pi,1) = fnode_[end_nodes.first + pi];
            }
          if(norm(uvw(colon(),0) - uvw(colon(),1)) < 1e-6){
              pair<size_t,size_t> cut_edge(
                    get_point_idx_from_node_idx(end_nodes.first).first,
                    get_point_idx_from_node_idx(end_nodes.second).first);
              if(cut_tet2tet[cut_edge.first] == cut_tet2tet[cut_edge.second])
                continue;
              if(mode == 1)
                return __LINE__;
              else
                trivial_loops.push_back(deq);
            }
        }
    }

  node_edge_type node_edge_after_simple_check = node_edges;
  for(size_t ci = 0; ci < trivial_loops.size(); ++ci){
      const deque<pair<size_t,size_t> > & deq = trivial_loops[ci];
      for(size_t ei = 0; ei < deq.size(); ++ei){
          pair<size_t,size_t> one_edge = deq[ei];
          if(one_edge.first > one_edge.second)
            swap(one_edge.first, one_edge.second);
          node_edge_type::iterator it = node_edge_after_simple_check.find(one_edge);
          if(it != node_edge_after_simple_check.end()){
              node_edge_after_simple_check.erase(it);
            }
        }
    }

  for(node_edge_type::const_iterator ncit = node_edge_after_simple_check.begin();
      ncit != node_edge_after_simple_check.end(); ++ncit){
      const pair<size_t,size_t> & node_edge = *ncit;

      const pair<size_t,size_t> fnode_edge(
            get_point_idx_from_node_idx(node_edge.first).first,
            get_point_idx_from_node_idx(node_edge.second).first);

      vector<uvw_type> uvw(2);
      for(size_t pi = 0;  pi < 3; ++pi){
          const pair<size_t,size_t> gnode_edge(
                get_node_idx_from_point_idx(fnode_edge.first, pi),
                get_node_idx_from_point_idx(fnode_edge.second, pi));

          uvw[0].push_back(fnode_[gnode_edge.first]);
          uvw[1].push_back(fnode_[gnode_edge.second]);
        }
      assert(node_edge.first % 3 == node_edge.second % 3);
      uvw[0].push_back(node_edge.first % 3);
      uvw[1].push_back(node_edge.second % 3);
      group2node[uvw[0]].insert(node_edge.first);
      group2node[uvw[1]].insert(node_edge.second);
    }

  map<size_t,size_t> p2p;
  map<size_t, boost::unordered_set<size_t> > degenerated_points;
  for(group2node_type::const_iterator gtcit = group2node.begin();
      gtcit != group2node.end(); ++gtcit){
      if(gtcit->second.size() > 1){
          // those points have the same uvw axis, then they must be the same point
          // in original tetmesh, or they are degenerated
          size_t original_node_idx = -1;
          const boost::unordered_set<size_t> & linked_points = gtcit->second;
          for(boost::unordered_set<size_t>::const_iterator buscit =
              linked_points.begin(); buscit != linked_points.end(); ++buscit){
              if(original_node_idx == -1){
                  original_node_idx =
                      cut_tet2tet[get_point_idx_from_node_idx(*buscit).first];
                }else if(cut_tet2tet[get_point_idx_from_node_idx(*buscit).first]
                         != original_node_idx)
                {
                  if(mode == 0){
                      cerr << "# [error] meet degenerated points "
                           << "cut_point "
                           <<  get_point_idx_from_node_idx(*buscit).first << " , "
                            << get_point_idx_from_node_idx(*(linked_points.begin())).first
                            << endl;

                      deque<pair<size_t,size_t> > one_chain;
                      one_chain.push_back(
                            make_pair(get_point_idx_from_node_idx(*buscit).first,
                                      get_point_idx_from_node_idx(*(linked_points.begin())).first));
                      loop_edges_.push_back(one_chain);

                    }
                  if(mode == 1)
                    return __LINE__;
                }
              // to make the index valid
              p2p[*buscit] = ((*buscit)%3 + 1)*fnode_.size() + *(linked_points.begin()) ;
              degenerated_points[((*buscit)%3 + 1)*fnode_.size() + *(linked_points.begin())]
                  .insert(*buscit);
              degenerated_points[((*buscit)%3 + 1)*fnode_.size() + *(linked_points.begin())]
                  .insert(*(linked_points.begin()));
            }
        }
    }

  // replace some node of edges with grouped new index
  boost::unordered_set<pair<size_t,size_t> > node_edge_new;
  map<pair<size_t,size_t>, vector<pair<size_t,size_t> > > new_edge2old_edge;
  // vector<deque<pair<size_t,size_t> > >  self_loops;
  for(node_edge_type::const_iterator ncit = node_edge_after_simple_check.begin();
      ncit != node_edge_after_simple_check.end(); ++ncit){
      pair<size_t,size_t> one_edge = *ncit;
      pair<size_t,size_t> old_edge(one_edge);
      map<size_t,size_t>::const_iterator bumcit_first =
          p2p.find(one_edge.first);
      map<size_t,size_t>::const_iterator bumcit_second =
          p2p.find(one_edge.second);

      if(bumcit_first != p2p.end()){
          one_edge.first = bumcit_first->second;
        }
      if(bumcit_second != p2p.end()){
          one_edge.second = bumcit_second->second;
        }

      new_edge2old_edge[one_edge].push_back(old_edge);

      if(one_edge.first > one_edge.second)
        swap(one_edge.first, one_edge.second);
      node_edge_new.insert(one_edge);
    }

  int rtn = extract_chains_from_undirected_edges(node_edge_new, chains_,
                                                 loops, mode);
  if(mode == 1 && rtn)
    return __LINE__; // fast return

  loops.insert(loops.end(), trivial_loops.begin(), trivial_loops.end());

  {
    typedef map<pair<size_t,size_t>, vector<pair<size_t,size_t> > >::const_iterator  new2old_it;
    for(size_t ci = 0; ci < chains_.size(); ++ci){
        for(size_t ei = 0; ei < chains_[ci].size(); ++ei){
            pair<size_t,size_t> & one_edge = chains_[ci][ei];
            //        new2old_it citt = new_edge2old_edge.find(one_edge);
            //        bool reverse = false;
            //        if(citt == new_edge2old_edge.end()){
            //          citt = new_edge2old_edge.find(make_pair(one_edge.second,
            //                                                  one_edge.first));
            //          reverse = true;
            //        }
            //        if(citt == new_edge2old_edge.end()){
            //          cerr << "# [error] strange, can not find " << one_edge.first << ","
            //               << one_edge.second << " to original edge mapping." << endl;
            //          return __LINE__;
            //        }
            //        if(!reverse)
            //          one_edge = citt->second.front();
            //        else
            //          one_edge = make_pair(citt->second.front().second,
            //                               citt->second.front().first);
            one_edge.first = one_edge.first % fnode_.size();
            one_edge.second = one_edge.second % fnode_.size();
          }
      }
  }

#if 0
  // check each chain to find near miss degeneration
  for(size_t ci = 0; ci < chains_.size(); ++ci) {
      const size_t &chain_end_0 = chains_[ci].front().first;
      const size_t &chain_end_1 = chains_[ci].back().second;

      const size_t point_0 = get_point_idx_from_node_idx(chain_end_0).first;
      const size_t point_1 = get_point_idx_from_node_idx(chain_end_1).first;
      const size_t point_2 = get_point_idx_from_node_idx(chains_[ci].front().second).first;

      // to detect the chain is inner or outside
      // we only handle the inner edges
      pair<size_t,size_t> orig_edge(cut_tet2tet[point_0], cut_tet2tet[point_2]);

      jtf::mesh::one_ring_tet_at_edge::e2tet_type::const_iterator ecit =
          ortae.e2t_.find(orig_edge);
      if(ecit == ortae.e2t_.end()){
          ecit = ortae.e2t_.find(make_pair(orig_edge.second, orig_edge.first));
        }
      if(ecit == ortae.e2t_.end()){
          cerr << "# [error] can not find edge < " << orig_edge.first << " "
               << orig_edge.second << " > cut_edge < "
               << point_0 << "," <<  point_1 << ">." << endl;
          return __LINE__;
        }
      const vector<size_t> & tet_loop = ecit->second;
      if(!ortae.is_inner_edge(tet_loop)) continue;

      boost::unordered_map<size_t, boost::unordered_set<size_t> >::const_iterator
          bcit_0 = point_around_faces_.find(point_0);
      boost::unordered_map<size_t, boost::unordered_set<size_t> >::const_iterator
          bcit_1 = point_around_faces_.find(point_1);
      if(bcit_0 == point_around_faces_.end() ||
         bcit_1 == point_around_faces_.end())
        continue;

      // to extract chain type
      const size_t free_axis = extract_free_axis_of_edge(point_0, point_1);
      assert(free_axis == 0 || free_axis == 1 || free_axis == 2);

      for(boost::unordered_set<size_t>::const_iterator bucit =
          bcit_0->second.begin(); bucit != bcit_0->second.end(); ++bucit){
          boost::unordered_map<size_t,size_t>::const_iterator sur2rot_it =
              surface_face2rot_idx.find(*bucit);
          if(sur2rot_it == surface_face2rot_idx.end()) {
              cerr << "# [error] can not find surface " << *bucit << " around point "
                   << point_0 << endl;
              return __LINE__;
            }
          assert(sur2rot_it->second < rot_type.size());
          if(rot_type[sur2rot_it->second] == -1) continue;
          if(rot_type[sur2rot_it->second] != free_axis)
            return __LINE__;
        }

      for(boost::unordered_set<size_t>::const_iterator bucit =
          bcit_1->second.begin(); bucit != bcit_1->second.end(); ++bucit){
          boost::unordered_map<size_t,size_t>::const_iterator sur2rot_it =
              surface_face2rot_idx.find(*bucit);
          if(sur2rot_it == surface_face2rot_idx.end()) {
              cerr << "# [error] can not find surface " << *bucit << " around point "
                   << point_1 << endl;
              return __LINE__;
            }
          assert(sur2rot_it->second < rot_type.size());
          if(rot_type[sur2rot_it->second] == -1) continue;
          if(rot_type[sur2rot_it->second] != free_axis)
            return __LINE__;
        }
    }
#endif
  return 0;
}

bool singularity_graph::is_valid_with_info(
    const matrixst & cut_tet2tet,
    const std::vector<size_t> & rot_type,
    const jtf::mesh::one_ring_tet_at_edge & ortae,
    const boost::unordered_map<size_t,size_t> & surface_idx_to_rot_idx,
    const size_t mode)
{
  loop_edges_.clear();

  bool is_valid_graph = true;
  boost::unordered_set<pair<size_t,size_t> > edges;
  boost::unordered_map<pair<size_t,size_t>,size_t> cut_edges_with_type;

  for(boost::unordered_map<pair<size_t,size_t>,bool >::const_iterator
      bucit = possible_edge_.begin(); bucit != possible_edge_.end(); ++bucit){
      const pair<size_t,size_t> & one_edge_cut = bucit->first;

      vector<size_t> grouped_axes;
      for(size_t ai = 0; ai < 3; ++ai){
          pair<size_t,size_t> node_edge(
                get_node_idx_from_point_idx(one_edge_cut.first, ai),
                get_node_idx_from_point_idx(one_edge_cut.second, ai));
          if(fnode_[node_edge.first] == fnode_[node_edge.second])
            grouped_axes.push_back(ai);
        }
      if(grouped_axes.size() < 2) continue;
      if(grouped_axes.size() == 2){
          const size_t axis_need_to_build_path =
              (0+1+2) - (grouped_axes[0] + grouped_axes[1]);//(std::accumulate(grouped_axes.begin(), grouped_axes.end(),0));
          pair<size_t,size_t> path(
                get_node_idx_from_point_idx(one_edge_cut.first, axis_need_to_build_path),
                get_node_idx_from_point_idx(one_edge_cut.second, axis_need_to_build_path));
          edges.insert(path);
          if(bucit->second)
            cut_edges_with_type[one_edge_cut] = axis_need_to_build_path;
        }
      if(grouped_axes.size() == 3){
#ifdef debug
          cerr << "# [info] meet compound edge: " << endl;
          cerr << "# ------ node "
               << get_node_idx_from_point_idx(one_edge_cut.first, 0) << ","
               << get_node_idx_from_point_idx(one_edge_cut.second, 0)
               << " cut_point: "<< one_edge_cut.first <<  ","
               << one_edge_cut.second << " orig edg <"
               << cut_tet2tet[one_edge_cut.first] << ","
               << cut_tet2tet[one_edge_cut.second] << ">." << endl;

          deque<pair<size_t,size_t> > compound_chain;
          compound_chain.push_back(one_edge_cut);
          loop_edges_.push_back(compound_chain);

          //      compound_edges[0].push_back(make_pair(one_edge_cut.first,
          //                                            one_edge_cut.second));
          is_valid_graph = false;
#endif
          if(mode == 1)
            return false;

        }
    }

  // check each restricted edge, if it's an inner edge in original tet mesh
  // and has surface points, then check the arounding patches around it, if
  // those patches' normal align type is not suitable with such edge point,
  // return false
  //#define near_miss_limit_weak
#ifdef near_miss_limit_strong
  {
    typedef boost::unordered_map<size_t, boost::unordered_set<size_t> > point2adj_face_type;
    for(boost::unordered_map<pair<size_t,size_t>,size_t>::const_iterator ecit =
        cut_edges_with_type.begin(); ecit != cut_edges_with_type.end(); ++ecit){
        const size_t &edge_type = ecit->second;

        vector<size_t> cut_edge_vec(2);
        cut_edge_vec[0] = ecit->first.first;
        cut_edge_vec[1] = ecit->first.second;

        for(size_t pi = 0; pi < cut_edge_vec.size(); ++pi){
            point2adj_face_type::const_iterator pcit =
                point_around_faces_with_rot_idx.find(cut_edge_vec[pi]);
            if(pcit == point_around_faces_with_rot_idx.end()) continue; // not surface point
            const boost::unordered_set<size_t> & linked_faces_idx = pcit->second;
            for(boost::unordered_set<size_t>::const_iterator bcit =
                linked_faces_idx.begin(); bcit != linked_faces_idx.end(); ++bcit){
                if(rot_type[*bcit] == -1) continue;
                if(rot_type[*bcit] != edge_type)
                  return false;
              }
          }
      }
  }
#endif

#ifdef near_miss_limit_weak
  {
    typedef boost::unordered_map<size_t, boost::unordered_set<size_t> > point2adj_face_type;
    for(boost::unordered_map<pair<size_t,size_t>,size_t>::const_iterator ecit =
        cut_edges_with_type.begin(); ecit != cut_edges_with_type.end(); ++ecit){
        const size_t &edge_type = ecit->second;

        vector<size_t> cut_edge_vec(2);
        cut_edge_vec[0] = ecit->first.first;
        cut_edge_vec[1] = ecit->first.second;
        // -1: not surfacep oint
        //  0: has surface type align the singularity type
        //  1: no surface type align the singularity type
        vector<size_t> meet_align_surface(2,-1);
        for(size_t pi = 0; pi < cut_edge_vec.size(); ++pi){
            point2adj_face_type::const_iterator pcit =
                point_around_faces_with_rot_idx.find(cut_edge_vec[pi]);
            if(pcit == point_around_faces_with_rot_idx.end()) continue; // not surface point
            const boost::unordered_set<size_t> & linked_faces_idx = pcit->second;
            for(boost::unordered_set<size_t>::const_iterator bcit =
                linked_faces_idx.begin(); bcit != linked_faces_idx.end(); ++bcit){
                if(rot_type[*bcit] == -1) continue;
                if(rot_type[*bcit] == edge_type){
                    meet_align_surface[pi] = 0;
                    break;
                  }
              }
            if(meet_align_surface[pi] == -1)
              meet_align_surface[pi] = 1; // flag it as no surface type aligned
          }
        if(meet_align_surface[0] == 1 || meet_align_surface[1] == 1)
          return false;
      }
  }
#endif

  vector<deque<pair<size_t,size_t> > > chains;

  int rtn = extract_directed_chain_from_node_edges(
        edges, chains,cut_tet2tet, rot_type, ortae,
        surface_idx_to_rot_idx, mode);
  if(mode == 1 && rtn ){
      return false;
    }
  if(rtn) is_valid_graph = false;

  //jtf::util::extract_chain_from_edges(edges_vec, chains);

  for(size_t t = 0; t < chains.size(); ++t){
      deque<pair<size_t,size_t> > & one_chain = chains[t];
      if(one_chain.front().first == one_chain.back().second) // find a loop
        {
#ifdef debug
          for(size_t ei = 0; ei < one_chain.size(); ++ei){
              pair<size_t,size_t> & one_edge = one_chain[ei];
              if(one_edge.first > fnode_.size()-1){
                  one_edge.first -= fnode_.size();
                }
              if(one_edge.second > fnode_.size()-1){
                  one_edge.second -= fnode_.size();
                }
            }


          deque<pair<size_t,size_t> > one_chain_loop;

          cerr << "# [info] meet loop: " << endl;
          for(size_t ei = 0; ei < one_chain.size(); ++ei){
              const pair<size_t,size_t> & one_edge = one_chain[ei];
              pair<size_t,size_t> cut_edge(
                    get_point_idx_from_node_idx(one_edge.first%fnode_.size()).first,
                    get_point_idx_from_node_idx(one_edge.second%fnode_.size()).first);

              cerr << "# ------ cut point "
                   << cut_edge.first << "," << cut_edge.second;
              cerr << "  orig point "
                   << cut_tet2tet[cut_edge.first] << ","
                   << cut_tet2tet[cut_edge.second] << " type "
                   << one_edge.first % 3 << ", " << one_edge.second %3 << endl;
              one_chain_loop.push_back(make_pair(cut_edge.first,
                                                 cut_edge.second));
            }
          loop_edges_.push_back(one_chain_loop);
          is_valid_graph = false;
#endif
          if(mode == 1)
            return false;
        }
    }
  return is_valid_graph;
}

//bool singularity_graph::is_valid_with_info_accelerate(
//    const matrixst & cut_tet2tet,
//    const std::vector<size_t> & rot_type,
//    const jtf::mesh::one_ring_tet_at_edge & ortae,
//    const boost::unordered_map<size_t,size_t> & surface_idx_to_rot_idx,
//    const size_t mode = 0 )
//{
//  cerr << "# "
//  boost::unordered_map<pair<size_t,size_t>,size_t> cut_edges_with_type;

//  return true;
//}

//! @brief return -1 : shit happen
//         return -2 : this edge is not ready
//         return 0-23: edge type
int singularity_graph::cal_singularity_type_at_given_tet(
    const vector<size_t> & tet_loop,
    const size_t &begin_tet,
    const vector<size_t> &rot_type,
    const boost::unordered_map<pair<size_t,size_t>,size_t> & tet_pair2rot_idx,
    std::vector<size_t> * reordered_tet_loop_ptr)
{
  vector<size_t> tet_loop_new = tet_loop;
  vector<size_t>::iterator it =
      find(tet_loop_new.begin(), tet_loop_new.end(), begin_tet);

  if(it == tet_loop_new.end())
    return -1;
  assert(tet_loop_new.front() == tet_loop_new.back());

  tet_loop_new.pop_back();
  std::rotate(tet_loop_new.begin(), it, tet_loop_new.end());
  assert(tet_loop_new.front() == begin_tet);
  tet_loop_new.push_back(tet_loop_new.front());

  if(reordered_tet_loop_ptr)
    *reordered_tet_loop_ptr = tet_loop_new;

  typedef boost::unordered_map<pair<size_t,size_t>,size_t>::const_iterator
      bupscit;
  matrixd rot = eye<double>(3);
  for(size_t t = 0; t < tet_loop_new.size() - 1; ++t){
      bupscit bit = tet_pair2rot_idx.find(make_pair(tet_loop_new[t],
                                                    tet_loop_new[t+1]));
      if(bit != tet_pair2rot_idx.end()){
          if(rot_type[bit->second] == -1){
              return -2; // this loop is not ready
            }
          rot = temp(rot * type_transition2(rot_type[bit->second]));
          continue;
        }else{
          bit = tet_pair2rot_idx.find(make_pair(tet_loop_new[t+1],
                                      tet_loop_new[t]));
          if(bit == tet_pair2rot_idx.end())
            continue; // it's a trivial face which is inside the cut mesh
          if(rot_type[bit->second] == -1)
            return -2; // this loop is not ready
          // continue;
          const size_t type = get_trans_type(rot_type[bit->second]);
          rot = temp(rot * type_transition2(type));
          continue;
        }
    }

  return type_transition1(rot);
}

/////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////
void equation_graph::assemble_equation_graph_generally(
    const std::vector<std::vector<size_t> > & node_group,
    const vector<pair<size_t,size_t> > &edges)
{
  assert(!node_group.empty());

  for(const auto & one_group: node_group)
    for(const auto & t: one_group)
      add_equal_constraint(t, *one_group.begin());

  convert2virtual_point(v2g_idx_,vp_);

  auto_add_inequal_cons_based_on_equal_cons(edges);
}

void equation_graph::assemble_equation_graph_generally(
    const std::vector<std::vector<std::vector<size_t> > > & node_group_dim,
    const vector<pair<size_t,size_t> > & edges)
{
  assert(!node_group_dim.empty());
  for(const auto & one_dim : node_group_dim)
    for(const auto & one_group: one_dim)
      for(const auto & t: one_group)
        add_equal_constraint(t, *one_group.begin());

  convert2virtual_point(v2g_idx_,vp_);

  auto_add_inequal_cons_based_on_equal_cons(edges);
}


void equation_graph::add_equal_constraint(const size_t &from, const size_t &to)
{
  assert(from < v2g_idx_.size() && to < v2g_idx_.size());
  if(v2g_idx_[from] == v2g_idx_[to]) return;
  const size_t group_need_be_cleared = v2g_idx_[from];
  for(group<size_t>::const_iterator cit = groups_[v2g_idx_[from]].begin();
      cit != groups_[v2g_idx_[from]].end(); ++cit)
    v2g_idx_[*cit] = v2g_idx_[to];
  groups_[v2g_idx_[to]].merge(groups_[group_need_be_cleared]);
}

void equation_graph::add_inequal_constraint(const tuple<size_t,size_t,size_t> & edge)
{
  if(get<0>(edge) < get<1>(edge))
    edges_.insert(edge);
  else
    edges_.insert(make_tuple(get<1>(edge), get<0>(edge), get<2>(edge)));
}

void equation_graph::auto_add_inequal_cons_based_on_equal_cons(
    const std::vector<std::pair<size_t,size_t> > &edges)
{
  vector<size_t> free_type;
  for(size_t ei = 0; ei < edges.size(); ++ei){
      const pair<size_t,size_t> & one_edge = edges[ei];
      free_type.clear();
      for(size_t di = 0; di < dim_; ++di){
          if(v2g_idx_[get_variant_idx(vp_[p2vp_[one_edge.first]], di)] !=
             v2g_idx_[get_variant_idx(vp_[p2vp_[one_edge.second]], di)])
            free_type.push_back(di);
        }
      if(free_type.size() == 1 ){ // one restricted edge
          edges_.insert(make_tuple(one_edge.first, one_edge.second,free_type.front()));
        }
      if(free_type.size() == 0){
          for(size_t di = 0; di < 2; ++di){
              edges_.insert(make_tuple(one_edge.first, one_edge.second,di));
            }
        }
    }
  //  {
  //    vector<size_t> edges;
  //    vector<size_t> type;
  //    for(const auto & edge : edges_){
  //        edges.push_back(get<0>(edge));
  //        edges.push_back(get<1>(edge));
  //        type.push_back(get<2>(edge));
  //      }
  //    ofstream ofs("edge.vtk");
  //    line2vtk(ofs, &(*node_)[0], node_->size(2), &edges[0], edges.size()/2);
  //    cell_data(ofs, &type[0], type.size(), "type");
  //  }
}

bool equation_graph::next_assignment(boost::dynamic_bitset<> & order,
                                     const vector<path_type> & orig_path,
                                     vector<path_type> & output_path,
                                     const vector<group<int> > & dir_v_idx)const
{
  size_t n = 0;
  while(n < order.size()){
      order.flip(n);
      if(order[n]==1) break;
      ++n;
    }
  if(n == order.size())
    return false;

  if(!order[0] && !order[1] && !order[2] && !order[3])
    cout << bitset_to_hex_string(order) << endl;

  output_path = orig_path;
  for(size_t i = 0; i < order.size(); ++i){
      if(order[i] == 1){
          const group<int> & g = dir_v_idx[i];
          for(const auto & path_idx : g)
            reverse_container(output_path[path_idx]);
        }
    }
  return true;
}


pair<size_t,size_t> equation_graph::find_path_intersection(
    const std::vector<size_t> & path_idx_0,
    const std::vector<size_t> & path_idx_1,
    const std::shared_ptr<vertex_connection<DIRECT> >  vc0,
    const std::shared_ptr<vertex_connection<DIRECT> >  vc1,
    const std::vector<path_type> &dir_path)const
{
  if(path_idx_0.empty() || path_idx_1.empty()) return make_pair(-1,-1);
  // for each path in path_idx_1, check whether it intersect with pathes
  // in path_idx_0, if yes return true, or return false

  // first find pathes on the same plane
  const size_t free_type_0 = get<2>(dir_path[path_idx_0.front()].front());
  const size_t free_type_1 = get<2>(dir_path[path_idx_1.front()].front());
  const size_t common_plane_type = (dim_-1)*dim_/2 - free_type_0 - free_type_1;
  assert(common_plane_type < dim_);

  vector<size_t> reachable_path;
  for(size_t i = 0; i < path_idx_0.size(); ++i){
      const path_type & path_i = dir_path[path_idx_0[i]];
      for(size_t j = 0; j < path_idx_1.size(); ++j){
          const path_type & path_j = dir_path[path_idx_1[j]];
          if(v2g_idx_[dim_*(get<0>(path_i.front())) + common_plane_type]
             != v2g_idx_[dim_*(get<0>(path_j.front())) + common_plane_type])
            continue;
          const pair<size_t,size_t> free_type_0_dir(v2g_idx_[get<0>(path_i.front())],
              v2g_idx_[get<1>(path_i.back())]);
          const pair<size_t,size_t> free_type_1_dir(v2g_idx_[get<0>(path_j.front())],
              v2g_idx_[get<1>(path_j.back())]);
          const size_t fix_type_0 = v2g_idx_[get_variant_idx(get<0>(path_i.front()), free_type_1)];
          const size_t fix_type_1 = v2g_idx_[get_variant_idx(get<0>(path_i.front()), free_type_0)];
          // here i check fix_type_0\in (free_type_0_dir.first,free_type_0_dir.second)
          // fix_type_1\in (free_type_1_dir.first,free_type_1_dir.second)
          if(fix_type_0 == free_type_0_dir.first || fix_type_0 == free_type_0_dir.second)
            continue;
          if(fix_type_1 == free_type_1_dir.first || fix_type_1 == free_type_1_dir.second)
            continue;
          int rtn = vc0->get_shortest_path(free_type_0_dir.first, fix_type_0, reachable_path);
          if(reachable_path.empty()) continue; // can not find a path
          rtn = vc0->get_shortest_path(fix_type_0, free_type_0_dir.second, reachable_path);
          if(reachable_path.empty()) continue; // can not find a path
          rtn = vc1->get_shortest_path(free_type_1_dir.first, fix_type_1, reachable_path);
          if(reachable_path.empty()) continue;
          rtn = vc1->get_shortest_path(fix_type_1, free_type_1_dir.second, reachable_path);
          if(reachable_path.empty()) continue;
          return make_pair(i,j); // find intersection
        }
    }
  return make_pair(-1,-1);
}

bool equation_graph::check_loop(
    const std::vector<path_type> & dir_path,
    bool show_info)
{
  vector<pair<size_t,size_t> > path_end_edges;
  for(size_t pi = 0; pi < dir_path.size(); ++pi){
      const path_type & one_path = dir_path[pi];
      path_end_edges.push_back(make_pair(get<0>(one_path.front()),
                                         get<1>(one_path.back())));
    }
  return cycle_detection::has_cycle<DIRECT2>(path_end_edges);
}


bool equation_graph::check_edge_intersection(
    const std::vector<path_type> & dir_path,
    bool show_info) const
{
  // separated dir_path into three variant graph which
  // helps to specify direction
  vector<map<pair<size_t,size_t>,double> > uvw_dir_edge_w(dim_);
  vector<vector<size_t> > uvw_path_idx(dim_);
  size_t path_idx = 0;
  for(const auto & path : dir_path){
      const size_t type = get<2>(path.front());
      uvw_dir_edge_w[type][
          make_pair(
            v2g_idx_[get_variant_idx(vp_[get<0>(path.front())], type)],
          v2g_idx_[get_variant_idx(vp_[get<1>(path.back())], type)])] = 1.0;
      uvw_path_idx[type].push_back(path_idx++);
    }
  vector<std::shared_ptr<vertex_connection<DIRECT> > > vc(dim_);
  for(size_t i = 0; i < dim_; ++i)
    vc[i].reset(vertex_connection<DIRECT>::create(uvw_dir_edge_w[i]));

  // for every two edge belong to different variant graph,
  // check whether they intersect to each other

  for(size_t i = 0; i < dim_; ++i)
    for(size_t j = i + 1; j < dim_; ++j){
        pair<size_t,size_t> has_found = find_path_intersection(
              uvw_path_idx[i], uvw_path_idx[j], vc[i],vc[j], dir_path);
        if(has_found.first != -1 && has_found.second != -1) {
            if(has_node() && show_info){
                vector<size_t> edges;
                ofstream ofs("edge_intersection.vtk");
                for(const auto & edge: dir_path[has_found.first]){
                    const auto it = vpedge2rawedge_.find(
                          make_pair(get<0>(edge),get<1>(edge)));
                    assert(it != vpedge2rawedge_.end());
                    for(const auto &raw_edge : it->second){
                        edges.push_back(raw_edge.first);
                        edges.push_back(raw_edge.second);
                      }
                  }
                for(const auto & edge: dir_path[has_found.second]){
                    const auto it = vpedge2rawedge_.find(
                          make_pair(get<0>(edge),get<1>(edge)));
                    assert(it != vpedge2rawedge_.end());
                    for(const auto &raw_edge : it->second){
                        edges.push_back(raw_edge.first);
                        edges.push_back(raw_edge.second);
                      }
                  }
                line2vtk(ofs, &(*node_)[0], node_->size(2), &edges[0], edges.size()/2);
              }
            return true;
          }
      }
  return false;
}
bool equation_graph::check_path_direction(
    const std::vector<path_type> & dir_path,
    bool show_info) const
{
  // uvw_point_degree is used to separate dir_path according to its free type
  // and for each ending, i record this path direction: p--> 1; -->p -1.
  // so that:
  // if ending degree > 2, error
  // if ending degree = 2, dir should be 1,-1, otherwise error
  vector<map<size_t,vector<int> > > uvw_point_degree(dim_);
  for(const auto & path : dir_path){
      const size_t type = get<2>(path.front());
      uvw_point_degree[type][get<0>(path.front())].push_back(1);
      uvw_point_degree[type][get<1>(path.back())].push_back(-1);
    }

  for(const auto & point_degree : uvw_point_degree){
      for(const auto & point2degree : point_degree){
          if(point2degree.second.size() > 2){
              if(has_node() && show_info){
                  ofstream ofs("path_degree_and_direction_of_point.vtk");
                  vector<size_t> points;
                  points.push_back(point2degree.first / dim_);
                  point2vtk(ofs, &(*node_)[0], node_->size(2), &points[0], points.size());
                }
              return false;
            }
          if(point2degree.second.size() == 2){
              const int acc =  std::accumulate(point2degree.second.begin(),
                                               point2degree.second.end(),0);
              if(acc != 0){
                  if(has_node() && show_info){
                      ofstream ofs("path_degree_and_direction_of_point.vtk");
                      vector<size_t> points;
                      points.push_back(point2degree.first / dim_);
                      point2vtk(ofs, &(*node_)[0], node_->size(2),&points[0], points.size());
                    }
                  return false;
                }
            }
        }
    }
  return true;
}

bool equation_graph::check_duplicate_type_path(
    const std::vector<path_type> & dir_path,
    bool show_info)const
{
  map<pair<size_t,size_t>, vector<size_t> > endings_to_idx;

  // for each path
  for(size_t pi = 0; pi < dir_path.size(); ++pi){
      const path_type & path = dir_path[pi];

      pair<size_t,size_t> both_ends(get<0>(path.front()), get<1>(path.back()));
      vector<size_t> & to_idx = endings_to_idx[both_ends];
      to_idx.push_back(pi);
      if(to_idx.size() != 1){
          if(has_node() && show_info){
              ofstream ofs("duplicated_type_path.vtk");
              vector<size_t> edges;
              for(size_t pj = 0; pj < to_idx.size(); ++pj){
                  const path_type &path_pj = dir_path[to_idx[pj]];
                  for(const auto & edge: path_pj){
                      const auto it = vpedge2rawedge_.find(make_pair(get<0>(edge),
                                                                     get<1>(edge)));
                      assert(it != vpedge2rawedge_.end());
                      for(const auto & raw_edge: it->second){
                          edges.push_back(raw_edge.first);
                          edges.push_back(raw_edge.second);
                        }
                    }
                }
              line2vtk(ofs, &(*node_)[0], node_->size(2), &edges[0], edges.size()/2);
            }
          return true;
        }
    }

  // for each edge
  endings_to_idx.clear();
  for(size_t pi = 0; pi < dir_path.size(); ++pi){
      const path_type & path = dir_path[pi];

      for(const auto & one_edge: path){
          pair<size_t,size_t> both_ends(get<0>(one_edge), get<1>(one_edge));
          vector<size_t> & to_idx = endings_to_idx[both_ends];
          to_idx.push_back(pi);
          if(to_idx.size() != 1){
              if(has_node() && show_info){
                  ofstream ofs("duplicated_type_path.vtk");
                  vector<size_t> edges;
                  for(size_t pj = 0; pj < to_idx.size(); ++pj){
                      const path_type &path_pj = dir_path[to_idx[pj]];
                      for(const auto & edge: path_pj){
                          const auto it = vpedge2rawedge_.find(make_pair(get<0>(edge),
                                                                         get<1>(edge)));
                          assert(it != vpedge2rawedge_.end());
                          for(const auto & raw_edge: it->second){
                              edges.push_back(raw_edge.first);
                              edges.push_back(raw_edge.second);
                            }
                        }
                    }
                  line2vtk(ofs, &(*node_)[0], node_->size(2), &edges[0], edges.size()/2);
                }
              return true;
            }
        }
    }

  return false;
}

bool equation_graph::check_self_loop(const std::vector<path_type> & dir_path,
                                     bool show_info)const
{
  for(const auto & g: dir_path){
      if(get<0>(g.front()) == get<1>(g.back())){
          if(has_node() && show_info){
              ofstream ofs("self_loop.vtk");
              vector<size_t> edges;
              for(const auto & one_edge : g){
                  const auto it = vpedge2rawedge_.find(make_pair(get<0>(one_edge), get<1>(one_edge)));
                  for(const auto & vpedge2rawedge_one: it->second){
                      edges.push_back(vpedge2rawedge_one.first);
                      edges.push_back(vpedge2rawedge_one.second);
                    }
                }
              line2vtk(ofs, &(*node_)[0], node_->size(2),&edges[0], edges.size()/2);
            }

          return true;
        }
    }
  return false;
}

void equation_graph::convert2virtual_point(
    const std::vector<size_t> & variant2group,
    std::vector<virtual_point> & vp)
{
  map<vector<size_t>,vector<size_t> > f2p;
  for(size_t i = 0 ; i < variant2group.size()/dim_; ++i){
      f2p[make_vector(&variant2group[i*dim_],dim_)].push_back(i);
    }
  vp.resize(f2p.size());
  p2vp_.resize(variant2group.size()/dim_);
  size_t idx = 0;
  for(const auto & f2o : f2p){
      vp[idx].orig_point = f2o.second;
      vp[idx].uvw_group = f2o.first;
      for(const auto & t : f2o.second)
        p2vp_[t] = idx;
      ++idx;
    }
}


equation_graph::graph_state equation_graph::check_directional_edge_graph_step0(
    const vector<path_type> & dir_path,
    bool show_info)
{
  if(check_self_loop(dir_path, show_info)){
      if(show_info)
        cerr << "# [error] find self loop." << endl;
      return graph_state::ABSOLUTE_FAIL;
    }
  if(check_duplicate_type_path(dir_path, show_info)){
      if(show_info)
        cerr << "# [error] find duplicate_type_path." << endl;
      return graph_state::ABSOLUTE_FAIL;
    }

  return graph_state::SUCCEED;
}

equation_graph::graph_state equation_graph::check_directional_edge_graph_step1(
    const vector<path_type> & dir_path,
    bool show_info)
{
  if(check_edge_intersection(dir_path,show_info))
    return graph_state::COMMON_FAIL;
  if(check_loop(dir_path, show_info))
    return graph_state::COMMON_FAIL;

  if(show_info)
    cerr << "# [info] graph is ok." << endl;
  return graph_state::SUCCEED;
}


bool equation_graph::check_valid(bool show_info)
{
  vector<path_type> dir_path;
  vector<group<int> > dir_v;

  edge2path(dir_path);

  {
    if(has_node() && show_info){
        dump_singularity_to_vtk("all_singularity.vtk", *node_ , edges_);
      }

    // Notice that, dir_path store virtual edges, not real edges
    ofstream ofs("restricted_edges");
    ofs << dir_path.size() << endl;
    for(size_t pi = 0; pi < dir_path.size(); ++pi){
        const auto & one_path = dir_path[pi];
        ofs << one_path.size() << endl;
        for(const auto & one_edge: one_path){
            const auto it = vpedge2rawedge_.find(make_pair(get<0>(one_edge), get<1>(one_edge)));
            assert(it != vpedge2rawedge_.end());

            ofs << it->second.front().first << " " << it->second.front().second
                <<  " " <<  std::get<2>(one_edge) << endl;
          }
      }
  }

  collect_direction_variant(dir_path, dir_v);

  size_t dir_v_num = 0;
  vector<group<int> > order2v;
  order2v.reserve(dir_v.size());
  std::for_each(dir_v.begin(), dir_v.end(), [&dir_v_num, &order2v](const group<int> & g)
  {if(!g.empty()){ ++dir_v_num; order2v.push_back(g);}});

  boost::dynamic_bitset<> order(dir_v_num);
  vector<path_type> dir_path_temp = dir_path;

  const graph_state st = check_directional_edge_graph_step0(dir_path_temp, show_info);
  if(st == graph_state::ABSOLUTE_FAIL) return false;

  const graph_state st2 = check_directional_edge_graph_step1(dir_path_temp, show_info);
  if(st2 == graph_state::ABSOLUTE_FAIL) return false;
  //  do{
  //      const graph_state st = check_directional_edge_graph_step1(dir_path_temp, show_info);
  //      if(st == graph_state::SUCCEED) return true;
  //      if(st == graph_state::ABSOLUTE_FAIL) return false;
  //    }while(next_assignment(order, dir_path, dir_path_temp, order2v));

  return true;
}

void equation_graph::edge2path(vector<path_type> & path)
{
  vector<set<pair<size_t,size_t> > > edges(dim_);
  vector<deque<pair<size_t,size_t> > > one_path;
  deque<tuple<size_t,size_t,size_t> > path_with_free_type;
  path.clear();

  for(const auto& one_edge: edges_){
      pair<size_t,size_t> vp_edge(p2vp_[get<0>(one_edge)], p2vp_[get<1>(one_edge)]);
      if(vp_edge.first > vp_edge.second){
          vpedge2rawedge_[vp_edge].push_back(make_pair(get<0>(one_edge),get<1>(one_edge)));
          swap(vp_edge.first, vp_edge.second);
          vpedge2rawedge_[vp_edge].push_back(make_pair(get<1>(one_edge),get<0>(one_edge)));
        }else{
          vpedge2rawedge_[vp_edge].push_back(make_pair(get<0>(one_edge),get<1>(one_edge)));
          vpedge2rawedge_[make_pair(vp_edge.second, vp_edge.first)].push_back(make_pair(get<1>(one_edge),get<0>(one_edge)));
        }
      edges[get<2>(one_edge)].insert(make_pair(get<0>(one_edge), get<1>(one_edge)));
    }

  for(size_t di = 0; di < dim_; ++di){
      vector<pair<size_t,size_t> > edge_vec(edges[di].begin(), edges[di].end());
      jtf::util::extract_chain_from_edges(edge_vec, one_path);
      for(size_t pi = 0; pi < one_path.size(); ++pi){
          path_with_free_type.resize(one_path[pi].size());
          for(size_t ei = 0; ei < one_path[pi].size(); ++ei){
              get<0>(path_with_free_type[ei]) = p2vp_[one_path[pi][ei].first];
              get<1>(path_with_free_type[ei]) = p2vp_[one_path[pi][ei].second];
              get<2>(path_with_free_type[ei]) = di;
            }
          path.push_back(path_with_free_type);
        }
    }
}

void equation_graph::collect_direction_variant(vector<path_type> & directional_path,
                                               vector<group<int> > & direction_v) const
{
  direction_v.clear();
  // first step:
  vector<vector<size_t> > uvw_pathes(3);
  for(size_t pi = 0; pi < directional_path.size(); ++pi){
      const path_type &one_path = directional_path[pi];
      const tuple<size_t,size_t,size_t> & one_edge = one_path.front();
      uvw_pathes[get<2>(one_edge)].push_back(pi);
    }

  direction_v.resize(directional_path.size());
  for(size_t i = 0; i < direction_v.size(); ++i) direction_v[i] << i;

  for(const auto &p: uvw_pathes){
      if(p.empty()) continue;
      unify_axis_order(p, directional_path, direction_v);
    }

#if 0 //        // visualize the grouped pathes
  {
    if(node_){
        ofstream ofs("grouped_edges.vtk");
        vector<size_t> edges;
        vector<size_t> edge_type;
        size_t type = -1;
        for(const auto & g : direction_v){
            if(g.empty()) continue;
            ++type;
            for(const size_t & path_idx : g){
                const path_type & path = directional_path[path_idx];
                for(const auto & one_edge : path){
                    edges.push_back(one_edge.first/dim_);
                    edges.push_back(one_edge.second/dim_);
                    edge_type.push_back(type);
                  }
              }
          }
        line2vtk(ofs, &((*node_)[0]), node_->size(2), &edges[0], edges.size()/2);
        cell_data(ofs, &edge_type[0], edge_type.size(), "group_edge_idx");

        directional_edge2arrow(&edges[0], edges.size()/2,
            &(*node_)[0], node_->size(2), "arrow.obj", "edge_arrow.obj");
      }
  }
#endif
}

void equation_graph::unify_axis_order(
    const std::vector<size_t> & path_idx_of_one_axis,
    std::vector<path_type> & directional_path,
    vector<group<int> > & direction_v)const
{
  const size_t free_type = get<2>(directional_path[path_idx_of_one_axis.front()].front());
  vector<pair<size_t,size_t> > path_end_pairs;
  // ends of one path, and its direction variant and order
  map<pair<size_t,size_t>,pair<size_t,int> > path_end_idx_with_direction;
  for(size_t i = 0; i < path_idx_of_one_axis.size(); ++i){
      const path_type & path = directional_path[path_idx_of_one_axis[i]];
      auto it = path_end_idx_with_direction.find(
            make_pair(get<0>(path.front()), get<1>(path.back())));

      if(it == path_end_idx_with_direction.end()) {
          path_end_pairs.push_back(
                make_pair(get<0>(path.front()), get<1>(path.back())));
          path_end_idx_with_direction[
              make_pair(get<0>(path.front()), get<1>(path.back()))]
              = make_pair(path_idx_of_one_axis[i],1);
          path_end_idx_with_direction[
              make_pair(get<1>(path.back()),get<0>(path.front()))]
              = make_pair(path_idx_of_one_axis[i], -1);
        }else{
          if(it->second.second == -1)
            reverse_container(directional_path[path_idx_of_one_axis[i]]);
          direction_v[it->second.first].merge(direction_v[path_idx_of_one_axis[i]]);
        }
    }

  vector<path_type> path_end_path;
  {
    vector<deque<pair<size_t,size_t> > > path_end_path_deq;
    jtf::util::extract_chain_from_edges(path_end_pairs, path_end_path_deq);
    {
      path_end_path.resize(path_end_path_deq.size());
      for(size_t di = 0; di < path_end_path_deq.size(); ++di){
          const auto & edge_path = path_end_path_deq[di];
          path_type one_path;
          const pair<size_t,size_t> first_edge_pair = path_end_path_deq[di].front();
          auto it = path_end_idx_with_direction.find(first_edge_pair);
          if(it == path_end_idx_with_direction.end())
            it = path_end_idx_with_direction.find(
                  make_pair(first_edge_pair.second, first_edge_pair.first));
          assert(it != path_end_idx_with_direction.end());

          for(const auto & one_edge: edge_path)
            one_path.push_back(make_tuple(one_edge.first, one_edge.second, it->second.first));

          path_end_path[di] = one_path;
        }
    }

    class path_reverser{
    public:
      static void func(const path_type & one_path_path,
                       map<pair<size_t,size_t>,pair<size_t,int> > &path_end_idx_with_direction,
                       vector<group<int> > &direction_v,
                       vector<path_type> & directional_path,
                       size_t & first_path_idx){
        for(auto &edge : one_path_path){
            auto it = path_end_idx_with_direction.find(
                  make_pair(get<0>(edge), get<1>(edge)));

            assert_iterator(path_end_idx_with_direction, it,
                            "#[error] path_end_pair should be found.");
            if(first_path_idx == -1)
              first_path_idx = it->second.first;
            if(it->second.second == -1){
                reverse_container(directional_path[it->second.first]);
                it->second.second *= -1;
                auto itt = path_end_idx_with_direction.find(
                      make_pair(get<1>(edge), get<0>(edge)));
                if(itt == path_end_idx_with_direction.end())
                  throw "#[error] stange can not find path end.";
                itt->second.second *= -1;
              }
            direction_v[first_path_idx].merge(direction_v[it->second.first]);
          }
      }
    };
    // for each path, unify their direction
    // 1. for these pathes which have the same virtual points, they can be connected ,
    // and direction variants be unified. This addressing is used for arbitrarily cutting
    // 2. unify all direction of pathes which connect by a shared point


    for(size_t i = 0; i < path_end_path.size(); ++i){
        const path_type & one_path_path = path_end_path[i]; // a path of path-endings
        size_t first_path_idx = -1;
        path_reverser::func(one_path_path, path_end_idx_with_direction,
                            direction_v, directional_path, first_path_idx);
      }
  }

  // for different path, if they have the same endings, unify their direction
  // first: path endings, second: directional_variant(), direction
  map<pair<size_t,size_t>, pair<size_t,int> > path_end_g_dir_v;
  for(size_t i = 0; i < path_end_path.size(); ++i){

      size_t dir_v = -1;
      {// get direction variant of this path
        const path_type & one_path_path = path_end_path[i];
        auto itt = path_end_idx_with_direction.find(
              make_pair(get<0>(one_path_path.front()),get<1>(one_path_path.front())));
        assert_iterator(path_end_idx_with_direction, itt,
                        "#[error] can not find direction variant of path end.");
        dir_v = itt->second.first;
      }

      const pair<size_t,size_t> ends =
          make_pair(v2g_idx_[get_variant_idx(vp_[get<0>(path_end_path[i].front())],free_type)],
          v2g_idx_[get_variant_idx(vp_[get<1>(path_end_path[i].back())], free_type)]);
      auto it = path_end_g_dir_v.find(ends);
      if(it == path_end_g_dir_v.end()) {
          path_end_g_dir_v[ends] = make_pair(dir_v,1);
          path_end_g_dir_v[make_pair(ends.second, ends.first)] = make_pair(dir_v,-1);
        }else{
          if(it->second.second == -1){
              direction_v[it->second.first].merge(direction_v[dir_v]);

              const path_type & one_path_path = path_end_path[i]; // a path of path-endings

              for(const auto & edge: one_path_path){
                  auto itt = path_end_idx_with_direction.find(
                        make_pair(get<0>(edge), get<1>(edge)));
                  assert_iterator(path_end_idx_with_direction, itt,
                                  "#[error] path_end_pair should be found.");

                  if(itt->second.second == 1){
                      reverse_container(directional_path[itt->second.first]);
                      itt->second.second *= -1;
                      auto ittt = path_end_idx_with_direction.find(
                            make_pair(get<1>(edge), get<0>(edge)));
                      assert_iterator(path_end_idx_with_direction, ittt,
                                      "#[error] stange can not find path end.");
                      ittt->second.second *= -1;
                    }
                }
            }else if (it->second.second == 1){ // this path direction is the same as current
              direction_v[it->second.first].merge(direction_v[dir_v]);
            }
        }
    }
}
