#include <iostream>
#include <fstream>
#include <map>
#include <deque>

#include "../../common/vtk.h"
#include <jtflib/function/operation.h>
#include <jtflib/algorithm/equation.h>
#include "../solver/solver_jtf.h"
#include "../solver/solver_ipopt.h"
#include "../descriptor/descriptor_base.h"
#include "../descriptor/descriptor_vol.h"
#include "../descriptor/func_terms/linear_equation.h"
#include "../common/util.h"
#include "../../common/transition_type.h"
#include "../../hex_param/io.h"
#include "util.h"

using namespace std;
using namespace zjucad::matrix;

class rounding_check
{
public:
  rounding_check(const vector<vector<size_t> > &iv,
                 const vector<vector<pair<size_t,size_t> > > & ure)
    :integer_variants_(iv), uvw_restricted_edges_(ure){
    for(size_t i = 0; i < integer_variants_.size(); ++i){
        for(size_t vi = 0 ; vi < integer_variants_[i].size(); ++vi){
            integer_v2g_[integer_variants_[i][vi]] = i;
          }
      }
  }

  ///
  /// \brief add_integer_value
  /// \param v_idx
  /// \param value
  /// \return
  ///
  int add_integer_value(const size_t v_idx, const int value,
                        const matrix<double> & node)
  {
    const auto it = integer_val_.find(v_idx);
    if(it != integer_val_.end()) {
        cerr << "# [error] " << endl;
        return __LINE__;
      }
    const auto g_it = integer_v2g_.find(v_idx);
    assert(g_it != integer_v2g_.end());
    map<size_t,int> temp_integer_val = integer_val_;
    for(size_t vi = 0; vi < integer_variants_[g_it->second].size(); ++vi){
        const size_t &one_v = integer_variants_[g_it->second][vi];
        integer_val_[one_v] = value;
      }
    bool rtn = has_degenerated_restriced_edges();
    if(rtn == true) {
        integer_val_ = temp_integer_val;
        return 1;
      }
    rtn = has_coupled_restriced_edges(node);
    if(rtn == true) {
        integer_val_ = temp_integer_val;
        return 1;
      }
    return 0;
  }


  bool has_degenerated_restriced_edges()const
  {
    for(size_t di = 0; di < uvw_restricted_edges_.size(); ++di){
        const vector<pair<size_t,size_t> > & one_group_edges = uvw_restricted_edges_[di];
        for(const auto & one_edge : one_group_edges){
            const auto it_0 = integer_val_.find(3*one_edge.first+di);
            const auto it_1 = integer_val_.find(3*one_edge.second+di);
            if(it_0 == integer_val_.end() || it_1 == integer_val_.end()) continue;
            if(it_0->second == it_1->second) return true;
          }
      }
    return false;
  }

  bool has_coupled_restriced_edges(const matrix<double> & node)const
  {
    for(size_t di = 0; di < uvw_restricted_edges_.size(); ++di){
        const vector<pair<size_t,size_t> > & one_group_edges = uvw_restricted_edges_[di];

        // step 1: check whether two restricted edges share the same variables
        // it store each restricted edge,
        // for example, if edge is u edge, the key of map is <v,w>, value is u range of edge
        map<pair<int,int>,vector<pair<int,int> > > plane2edge;
        for(size_t ei = 0; ei < one_group_edges.size(); ++ei){
            const pair<size_t,size_t> & one_edge = one_group_edges[ei];
            const auto it_0 = integer_val_.find(3*one_edge.first+di);
            const auto it_1 = integer_val_.find(3*one_edge.second+di);
            if(it_0 == integer_val_.end() || it_1 == integer_val_.end()) continue;
            assert(it_0->second != it_1->second);

            const auto it_0_1 = integer_val_.find(3*one_edge.first+(di+1)%3);
            const auto it_1_1 = integer_val_.find(3*one_edge.second+(di+1)%3);
            if(it_0_1 == integer_val_.end() || it_1_1 == integer_val_.end()) continue;
            assert(it_0_1->second == it_1_1->second);

            const auto it_0_2 = integer_val_.find(3*one_edge.first+(di+2)%3);
            const auto it_1_2 = integer_val_.find(3*one_edge.second+(di+2)%3);
            if(it_0_2 == integer_val_.end() || it_1_2 == integer_val_.end()) continue;
            assert(it_0_2->second == it_1_2->second);

            pair<int,int> v_range(it_0->second, it_1->second);
            if(v_range.first > v_range.second)
              swap(v_range.first, v_range.second);

            plane2edge[make_pair(it_0_1->second, it_0_2->second)].push_back(v_range);
          }

        if(plane2edge.empty()) continue;
        for(auto & one_plane : plane2edge){
            if(one_plane.second.empty()) continue;
            vector<pair<int,int> > & edges = one_plane.second;
            sort(edges.begin(), edges.end());
            for(size_t i = 1; i < edges.size(); ++i){
                if(edges[i].first < edges[i-1].second) return true;
              }
          }

        // step 2: check under current nodes, whether there are edges overlaped
        // this function assume input node has already satisify all linear constraints (inner transition and surface align)
        {

        }

      }
    return false;
  }
private:
  rounding_check();
  rounding_check(const rounding_check &);

  std::map<size_t, int> integer_val_;
  std::map<size_t, size_t> integer_v2g_;
  const vector<vector<size_t> > &integer_variants_;
  const vector<vector<pair<size_t,size_t> > > &uvw_restricted_edges_;
};


///
/// @brief update_integer_variants
/// @param integer_variants, it stores all integer variants, each group share the same integer variants
/// @param uvw_v2g           it contains three maps (u/v/w), in each map, it stores variant to group in integer_variants
/// @param integer_fixed     stores bool flags for integer_variants, means whether an integer  group is fixed
/// @param rc                rounding checker
/// @param node              input node
/// @param integer_value     fixed integer value
/// @return
///
int update_integer_variants(const vector<vector<size_t> > &integer_variants,
                            const vector<map<size_t,size_t> > &uvw_v2g,
                            vector<bool> &integer_fixed,
                            rounding_check & rc,
                            const matrix<double> & node,
                            std::vector<pair<size_t, int> > & integer_group_value)
{
  integer_group_value.clear();
  // each tuple stores min_diff to one side integer, max_diff_to_other_side_integer,
  // one_side_integer, other_side_integer, group
  std::deque<tuple<double,double, int,int,size_t> > integer_diff;
  tuple<double, double, int, int, size_t> one_tuple;
  for(size_t i = 0; i < integer_fixed.size(); ++i){
      if(!integer_fixed[i]){
          const size_t & variant = integer_variants[i].front();
          const double &value = node[variant];
          int left = std::floor(value);
          int right = std::ceil(value);
          const double left_diff = fabs(value-left);
          const double right_diff = fabs(value-right);
          if(left_diff < right_diff){
              get<0>(one_tuple) = left_diff;
              get<1>(one_tuple) = right_diff;
              get<2>(one_tuple) = left;
              get<3>(one_tuple) = right;
              get<4>(one_tuple) = i;
            }else{
              get<0>(one_tuple) = right_diff;
              get<1>(one_tuple) = left_diff;
              get<2>(one_tuple) = right;
              get<3>(one_tuple) = left;
              get<4>(one_tuple) = i;
            }
          integer_diff.push_back(one_tuple);
        }
    }

  if(integer_diff.empty()) return 0;

  sort(integer_diff.begin(), integer_diff.end());

  while(1){
      size_t left_number = integer_diff.size();

      for(size_t i = 0; i < left_number; ++i){
          const auto & diff_tuple = integer_diff.front();
          int rtn =
              rc.add_integer_value(
                integer_variants[get<4>(diff_tuple)].front(), get<2>(diff_tuple),node);
          if(rtn == 0){
              integer_group_value.push_back(
                    make_pair(get<4>(diff_tuple), get<2>(diff_tuple)));
              integer_fixed[get<4>(diff_tuple)] = true;
              integer_diff.pop_front();
              break;
            }else{
              integer_diff.push_back(std::make_tuple(
                                       get<1>(diff_tuple), get<0>(diff_tuple),
                                       get<3>(diff_tuple), get<2>(diff_tuple),
                                       get<4>(diff_tuple)));
              integer_diff.pop_front();
            }
        }
      if(integer_diff.size() != left_number) break;
      break;
      cerr << "# [invalid rounding.] " << endl;
    }

  return 0;
}

int update_integer_variants2(const vector<vector<size_t> > &integer_variants,
                             const vector<vector<size_t> > &uvw_g,
                             vector<bool> &integer_fixed,
                             rounding_check & rc,
                             const matrix<double> & node,
                             std::vector<pair<size_t, int> > & integer_group_value,
                             std::vector<vector<int> > &order_integer_value_uvw )
{
  integer_group_value.clear();
  size_t windows = 0;
  for(size_t i = 0; i != uvw_g.size(); ++i){
      const vector<size_t> & one_g = uvw_g[i];
      for(size_t vi = 0; vi < one_g.size(); ++vi){
          const size_t & g_idx = one_g[vi];
          if(!integer_fixed[g_idx]){
              const double &value = node[integer_variants[g_idx].front()];
              int left = std::floor(value);
              int right = std::ceil(value);

              int near_integer = fabs(value-left) < fabs(value-right)?left:right;
              int far_integer = left + right - near_integer;

              int integer = near_integer;

              int rtn = rc.add_integer_value(integer_variants[g_idx].front(), near_integer,node);
              if(rtn != 0){
                  rtn =  rc.add_integer_value(integer_variants[g_idx].front(), far_integer,node);
                  integer = far_integer;
                }
              if(rtn == 0){ // valid rounding
                  integer_group_value.push_back(make_pair(g_idx, integer));
                  order_integer_value_uvw[i].push_back(integer);
                  integer_fixed[g_idx] = true;
                  break;
                }else{ // totally invalid
                  int prev_integer;
                  if(vi == 0){
                      cerr << "# [error] strange, why will the first variant of axis " << i
                           << " can not find an integer value." << endl;
                      //return __LINE__;
                    }
                  if(order_integer_value_uvw[i].empty())
                    prev_integer = near_integer;
                  else
                    prev_integer = order_integer_value_uvw[i].back();
                  if(value > prev_integer)
                    integer = std::ceil(value-prev_integer) + prev_integer;
                  else
                    integer = prev_integer + 1;
                  rtn = rc.add_integer_value(integer_variants[g_idx].front(), integer,node);
                  if(rtn != 0) {
                      cerr << "# [error] strange, this integer should be valid." << endl;
                      cerr << "# [info] begin to try integers 2 times: current value " << value << endl;
                      for(size_t try_i = 0; try_i < 2; ++try_i){
                          integer += try_i+1;
                          rtn = rc.add_integer_value(integer_variants[g_idx].front(), integer,node);
                          cerr << "# [info] +++++++ try v = "  << integer ;
                          if(rtn == 0) {
                              cerr << " seccess." << endl;
                              break;
                            }
                          cerr << " fail." << endl;
                        }
                      if(rtn != 0){
                          cerr << "# [error] can not find an integer, I simply choose " << std::ceil(value)+1
                               << "to enforce the integer rounding."<< endl;
                          integer = std::ceil(value)+1;
                        }
                    }
                  integer_group_value.push_back(make_pair(g_idx, integer));
                  order_integer_value_uvw[i].push_back(integer);
                  integer_fixed[g_idx] = true;
                  return 0;
                }
            }
        }
    }
  return 0;

}

void init_group_order(const matrix<double> & node,
                      const vector<vector<size_t> > & integer_variants,
                      vector<vector<size_t> > & uvw_g)
{
  if(uvw_g.size() != 3) uvw_g.resize(3);
  // this function assume that parameterization is already applied, thus
  // in each group, variants should share almost the same value

  vector<vector<pair<double,size_t> > > diff_variant(3);
  for(size_t i = 0; i < integer_variants.size(); ++i){
      const vector<size_t> & one_group = integer_variants[i];
      diff_variant[one_group.front()%3].push_back(make_pair(node[one_group.front()], i));
    }

  for(size_t di = 0; di < 3; ++di){
      sort(diff_variant[di].begin(), diff_variant[di].end());
      uvw_g[di].resize(diff_variant[di].size());
      for(size_t i = 0; i < diff_variant[di].size(); ++i){
          uvw_g[di][i] = diff_variant[di][i].second;
        }
    }
}

///
/// \brief eliminate_one_eqn_and_remove_variables, original node mapping will be
///
/// \param des_vol
/// \param fix_integer
/// \return 0: Z/q has been changed, 1: nothing changed, others: error
///
int eliminate_one_eqn_and_remove_variables(std::shared_ptr<descriptor_vol> des_vol,
                                           pair<size_t,double> fix_integer)
{
  vector<jtf_func_ptr> &eqns = des_vol->get_eqn_constraint();
  if(eqns.size() > 1){
      cerr << "# [error] only support one equation." << endl;
      return __LINE__;
    }

  eqns.clear();
  return des_vol->update_Z_q(fix_integer);
}

int eliminate_eqn_and_remove_variables(std::shared_ptr<descriptor_vol> des_vol,
                                       const vector<pair<size_t,int> > & fix_integer,
                                       size_t begin_idx)
{
  vector<jtf_func_ptr> &eqns = des_vol->get_eqn_constraint();

  eqns.clear();
  for(size_t t = begin_idx; t != fix_integer.size(); ++t)
    des_vol->update_Z_q(make_pair(fix_integer[t].first, 1.0*fix_integer[t].second));
  return 0;
}

int integer_rounding_integer_variable(
    std::shared_ptr<descriptor_base> des,
    const matrix<size_t> &mesh,
    matrix<double> &node,
    boost::property_tree::ptree &pt,
    vector<vector<size_t> > &integer_variants,
    vector<vector<pair<size_t,size_t> > > &uvw_restricted_edges)
{
  shared_ptr<descriptor_vol> des_vol = dynamic_pointer_cast<descriptor_vol>(des);
  cerr << "# [info] tetmesh to this function should be parameterizated. " << endl;
  if(!des->has_node_mapping()){
      throw std::logic_error("I can not round integer variables without nodemapping.");
    }

  //store group idx and integer value
  std::vector<pair<size_t, int> >  integer_group_value;

  shared_ptr<solver_base> sv;

  switch(str2int(pt.get<string>("package.value").c_str())) {
    case str2int("ipopt"): {
        sv.reset(new solver_ipopt); break;
      }
    case str2int("jtf"): {
        sv.reset(new solver_jtf); break;
      }
dafault:
      throw std::invalid_argument("unknown solver.");
    }

  rounding_check rc(integer_variants, uvw_restricted_edges);
  vector<map<size_t,size_t> > uvw_v2g(3); // uvw

  for(size_t gi = 0; gi < integer_variants.size(); ++gi){
      const vector<size_t> & one_integer_g = integer_variants[gi];
      for(const auto & one_v : one_integer_g){
          const size_t uvw = one_v % 3;
          uvw_v2g[uvw][one_v] = gi;
        }
    }

  vector<vector<size_t> > uvw_g(3); // groups in u/v/w it should be sort later.

  size_t i = 0;
  vector<pair<size_t,int> > variant2_ptr;
  size_t flag_v = 0;
  init_group_order(node,integer_variants, uvw_g);
  std::vector<vector<int> > order_integer_value_uvw(3) ;
  static vector<size_t> restricted_edges_vtk;

  const size_t f_node_num = max(mesh) + 1;
  matrix<double> init_node = node;
  vector<bool> integer_fixed(integer_variants.size(), false);

  size_t rounded_variables = 0;
  while(find(integer_fixed.begin(), integer_fixed.end(), false) != integer_fixed.end()){
      cerr << "# [info] rounded_variables " <<  rounded_variables << "/" << integer_variants.size() << endl;

      flag_v = variant2_ptr.size();
      update_integer_variants2(integer_variants, uvw_g, integer_fixed,rc, init_node,
                               integer_group_value, order_integer_value_uvw);
      if(integer_group_value.empty()){
          break;
          cerr << "# [error] invalid."<< endl;
        }

      bool need_to_round = false;
      for(size_t j = 0; j < integer_group_value.size(); ++j){

          const size_t idx = integer_variants[integer_group_value[j].first].front();
          variant2_ptr.push_back(
                make_pair(integer_variants[integer_group_value[j].first].front(),
                integer_group_value[j].second));

          //if(des_vol->bi_.NM.is_fixed[idx]) continue;
          //          des->add_eqn_constraint(
          //                jtf_func_cons_ptr(
          //                  jtf_func_cons_ptr(
          //                    new variant_fix(node.size(),idx,integer_group_value[j].second,
          //                                    (des->has_node_mapping()?&des->get_node_mapping():0)))),1);

          des->add_eqn_constraint(
                jtf_func_ptr(
                  jtf::function::least_square_warpper(
                    shared_ptr<const hj::function::function_t<double,int32_t> >(
                      new variant_fix_hj2(node.size(),idx,integer_group_value[j].second,
                                          (des->has_node_mapping()?&des->get_node_mapping():0),10)))),0);
          need_to_round = true;
        }

      if(need_to_round == false) continue;

      rounded_variables += integer_group_value.size();
      cerr << "# before integer fix: " << endl;
      for(size_t j = 0; j < variant2_ptr.size(); ++j){
          cerr << "\t v" << variant2_ptr[j].first << "\t"
               << variant2_ptr[j].second << "\t"
               << init_node[variant2_ptr[j].first] << endl;
        }

      if(init_node.size(2) == f_node_num){
          assert(des->has_node_mapping());
          const hj::sparse::csc<double,int32_t> &NMT = des->get_node_mapping().ZT;
          const size_t vnum = NMT.size(2);
          assert(vnum %3 == 0);
          init_node = zeros<double>(3, vnum/3);
          std::copy(node.begin(), node.end(), init_node.begin());
          des_vol->recover_gap_node(init_node);
        }

      sv->solve(init_node, des, pt);

      cerr << "# after integer fix: " << endl;
      for(size_t j = 0; j < variant2_ptr.size(); ++j){
          cerr << "\t v" << variant2_ptr[j].first << "\t"
               << variant2_ptr[j].second << "\t"
               << init_node[variant2_ptr[j].first] << endl;
        }

      //      {// eliminate variables
      //        int rtn = eliminate_one_eqn_and_remove_variables(
      //              des_vol, std::make_pair(variant2_ptr.back().first, 1.0*variant2_ptr.back().second));
      //        //TODO: I hack it!
      //        //        if(rtn != 0 && rtn != 1) {
      //        //            //throw std::logic_error("conflict integer rounding");
      //        //          }
      //        //        if(rtn == 0) { // Z/q has been changed
      //        des->set_objective(mesh, node, pt);
      //        //          }
      //      }

      {
        for(size_t i = flag_v; i < variant2_ptr.size(); ++i){
            eliminate_eqn_and_remove_variables(des_vol, variant2_ptr, flag_v);
          }
        des->set_objective(mesh, node, pt);
      }

      /////////////////////////////////////////////////////////////////////////
      //// visualize

      {
        stringstream ss ;
        ss << "output_integer_" << i++ <<".vtk";
        ofstream ofs(ss.str().c_str());
        tet2vtk(ofs, &init_node[0], init_node.size(2), &mesh[0], mesh.size(2));
      }

      {
        if(restricted_edges_vtk.size() == 0){
            for(size_t di = 0; di < uvw_restricted_edges.size(); ++di){
                const vector<pair<size_t,size_t> > & one_group = uvw_restricted_edges[di];
                for(const auto & one_edge : one_group){
                    restricted_edges_vtk.push_back(one_edge.first);
                    restricted_edges_vtk.push_back(one_edge.second);
                  }
              }
          }
        stringstream ss;
        ss << "restricted_edges_" << i-1 << ".vtk";
        ofstream ofs(ss.str().c_str());
        line2vtk(ofs, &init_node[0], init_node.size(2), &restricted_edges_vtk[0], restricted_edges_vtk.size()/2);
      }
    }

  sv->solve(init_node, des, pt); //

  cerr << "# after integer fix: " << endl;
  for(size_t j = 0; j < variant2_ptr.size(); ++j){
      cerr << "\t v" << variant2_ptr[j].first << "\t"
           << variant2_ptr[j].second << "\t"
           << init_node[variant2_ptr[j].first] << endl;
    }

  if(init_node.size() == node.size()) node = init_node;
  else
    {
      itr_matrix<const double *> node_m(node.size(1), node.size(2), &init_node[0]);
      node = node_m;
    }
  return 0;
}

void construct_mapping_face(
    const matrix<size_t> & uncut_tet,
    const matrix<size_t> & cut_tet,
    const boost::unordered_map<pair<size_t,size_t>,size_t> & inner_type,
    map<pair<vector<size_t>, vector<size_t> >, size_t> & ff2t)
{
  matrix<size_t> uncut_outside_faces;
  matrix<size_t> cut_outside_faces;

  unique_ptr<jtf::mesh::face2tet_adjacent> fa_cut(jtf::mesh::face2tet_adjacent::create(cut_tet));
  unique_ptr<jtf::mesh::face2tet_adjacent> fa_uncut(jtf::mesh::face2tet_adjacent::create(uncut_tet));

  if(!fa_cut.get() && !fa_uncut.get()){
      throw std::logic_error("can not build face2tet adjacent.");
    }

  jtf::mesh::get_outside_face(*fa_uncut, uncut_outside_faces);
  jtf::mesh::get_outside_face(*fa_cut, cut_outside_faces);

  map<vector<size_t>, vector<vector<size_t> > > orig_face2cut_faces;

  matrix<size_t> cut_tet2tet(max(cut_tet)+1);
  cut_tet2tet(cut_tet) = uncut_tet(colon());

  vector<size_t> one_face(3), one_face_cut(3);
  for(size_t fi = 0; fi < cut_outside_faces.size(2); ++fi){
      for(size_t pi = 0; pi < cut_outside_faces.size(1); ++pi){
          one_face[pi] = cut_tet2tet[cut_outside_faces(pi,fi)];
        }
      sort(one_face.begin(), one_face.end());

      for(size_t pi = 0; pi < cut_outside_faces.size(1); ++pi){
          for(size_t pj = 0; pj < cut_outside_faces.size(1); ++pj){
              if(cut_tet2tet[cut_outside_faces(pj, fi)] == one_face[pi]){
                  one_face_cut[pi] = cut_outside_faces(pj,fi);
                  break;
                }
            }
        }

      orig_face2cut_faces[one_face].push_back(one_face_cut);
    }

  for(const auto & one_face_uncut : orig_face2cut_faces){
      if(one_face_uncut.second.size() != 2) continue;
      const vector<vector<size_t> > & two_faces = one_face_uncut.second;
      const size_t face_idx_left = fa_cut->get_face_idx(&two_faces[0][0]);
      const size_t face_idx_right = fa_cut->get_face_idx(&two_faces[1][0]);
      assert(face_idx_left != -1 && face_idx_right != -1);
      const pair<size_t,size_t> & tet_pair_left = fa_cut->face2tet_[face_idx_left];
      const pair<size_t,size_t> & tet_pair_right = fa_cut->face2tet_[face_idx_right];

      assert(fa_cut->is_outside_face(tet_pair_left));
      assert(fa_cut->is_outside_face(tet_pair_right));

      pair<size_t,size_t> tet_pair(tet_pair_left.first==-1?tet_pair_left.second:tet_pair_left.first,
                                   tet_pair_right.first==-1?tet_pair_right.second:tet_pair_right.first);
      const auto type_it = inner_type.find(tet_pair);
      if(type_it == inner_type.end())
        ff2t[make_pair(two_faces.front(), two_faces.back())] = TRIVIAL_TYPE;
      else
        ff2t[make_pair(two_faces.front(), two_faces.back())] = type_it->second;
    }
}

template <typename T>
double jtf::math::get_sign(const T & v){
  return v > 0?1:-1;
}

void group_gap_variable(vector<size_t> & gap_variable2group,
                        const size_t f_node_num,
                        const hj::sparse::csc<double,int32_t> &NMT)
{
  const size_t total_variable_num = NMT.size(2);
  // form [f_node_num, total_variable_num) is gap variables
  typedef vector<pair<size_t,double> > equation_vec;
  map<equation_vec, vector<size_t> > eqn2variables;
  for(size_t i = 3*f_node_num; i < total_variable_num; ++i){
      equation_vec ev;
      for(size_t idx = NMT.ptr()[i]; idx != NMT.ptr()[i+1]; ++idx){
          ev.push_back(make_pair(NMT.idx()[idx], NMT.val()[idx]));
        }
      for(size_t pi = 0; pi < ev.size(); ++pi){
          ev[pi].second /= ev[0].second;
        }
      eqn2variables[ev].push_back(i-3*f_node_num);
    }

  if(gap_variable2group.size() != total_variable_num - 3*f_node_num)
    gap_variable2group.resize(total_variable_num - 3*f_node_num);
  for(size_t i = 0; i < gap_variable2group.size(); ++i)
    gap_variable2group[i] = i;

  for(const auto & one_eqn : eqn2variables){
      const auto & group_of_eqn = one_eqn.second;
      for(size_t idx_i = 1; idx_i < group_of_eqn.size(); ++idx_i){
          gap_variable2group.at(group_of_eqn[idx_i]) = gap_variable2group.at(group_of_eqn[0]);
        }
    }
}

// check each gap of face pair, to ensure that gap should be integer
int integer_rounding_gap(std::shared_ptr<descriptor_base> des,
                         const matrix<size_t> &mesh,
                         matrix<double> &node,
                         boost::property_tree::ptree &pt)
{
  shared_ptr<descriptor_vol> des_vol = dynamic_pointer_cast<descriptor_vol>(des);
  const size_t f_node_num = max(mesh)+1;

  const hj::sparse::csc<double,int32_t> & NMT = des_vol->get_node_mapping().ZT;

  matrix<double> init_node;
  if(node.size(2) == f_node_num){
      const size_t vnum = NMT.size(2);
      assert(vnum %3 == 0);
      init_node = zeros<double>(3, vnum/3);
      std::copy(node.begin(), node.end(), init_node.begin());
      des_vol->recover_gap_node(init_node);
    }

  const size_t gap_number  = init_node.size(2) - f_node_num;
  cerr << "# [info] " << gap_number << " gaps need to be rounded." << endl;

  const set<size_t> & zero_index = des_vol->bi_.zero_index_; // it stores all zero index which is assigned in objective function
  shared_ptr<solver_base> sv;

  switch(str2int(pt.get<string>("package.value").c_str())) {
    case str2int("ipopt"): {
        sv.reset(new solver_ipopt); break;
      }
    case str2int("jtf"): {
        sv.reset(new solver_jtf); break;
      }
dafault:
      throw std::invalid_argument("unknown solver.");
    }

  //  // to avoid KKT matrix become singulary, I should test each gap variable
  //  // if two gap variabe g_i = a*g_j, then g_j should not be added if g_i is already there
  //  vector<bool> gap_added(gap_number * 3,false);
  //  vector<size_t> gap_variable2group(gap_number * 3);

  //group_gap_variable(gap_variable2group, f_node_num, NMT);

  matrix<int> targets(3);
  pair<size_t,double> fix_integer;

  ofstream ofs_gaps("gaps");

  for(size_t gap_idx = 0; gap_idx < gap_number; ++gap_idx){
      ofs_gaps << init_node(colon(), gap_idx + f_node_num);
      targets *= 0;
      bool is_add_eqn = false;
      for(size_t di = 0; di < 3 ; ++di){
          if(find(zero_index.begin(), zero_index.end(), 3*(f_node_num+gap_idx)+di)
             != zero_index.end()) continue;
          //if(gap_added[gap_variable2group[3*gap_idx+di]]) continue;

          if(des_vol->bi_.NM.q[3*(f_node_num+gap_idx)+di]) continue;

          if(des_vol->bi_.NM.is_irrelevant_variable(3*(f_node_num+gap_idx)+di))
            continue;

          const int val_ceil = std::ceil(init_node(di, gap_idx + f_node_num));
          const int val_floor = std::floor(init_node(di, gap_idx + f_node_num));
          const double dis_ceil = fabs(init_node(di, gap_idx + f_node_num) - val_ceil);
          const double dis_floor = fabs(init_node(di, gap_idx + f_node_num) - val_floor);

          int target = 0;
          if(dis_ceil < dis_floor)
            target = val_ceil;
          else
            target = val_floor;
          targets[di] = target;

          {
            des->add_eqn_constraint(
                  jtf_func_ptr (
                    jtf::function::least_square_warpper(
                      shared_ptr<const hj::function::function_t<double,int32_t> >(
                        new variant_fix_hj2(
                          node.size(),3*(f_node_num+gap_idx)+di, target,
                          (des->has_node_mapping()?&des->get_node_mapping():0),100)))),0);

            fix_integer.first = 3*(f_node_num+gap_idx)+di;
            fix_integer.second = target;
          }

          is_add_eqn = true;
          //gap_added[gap_variable2group[3*gap_idx+di]] = true;
          break;
        }
      if(is_add_eqn == false) continue;
      cerr << "# [info] before rounding gap " << gap_idx << " " << init_node(colon(), gap_idx + f_node_num) << endl;
      cerr << "# [info] gap target " << targets << endl;

      sv->solve(init_node, des, pt);
      cerr << "# [info] after rounding gap " << gap_idx << " " << init_node(colon(), gap_idx + f_node_num) << endl;

      {// eliminate variables
        int rtn = eliminate_one_eqn_and_remove_variables(
              des_vol,fix_integer);
        if(rtn != 0 && rtn != 1) {
            throw std::logic_error("conflict integer rounding");
          }
        if(rtn == 0) { // Z/q has been changed
            des->set_objective(mesh, node, pt);
          }
      }
    }

  sv->solve(init_node, des, pt);
  if(init_node.size() == node.size()) node = init_node;
  else{
      itr_matrix<const double *> node_m(node.size(1), node.size(2), &init_node[0]);
      node = node_m;
    }

  return 0;
}

int integer_rounding(std::shared_ptr<descriptor_base> des,
                     const matrix<size_t> &mesh,
                     matrix<double> &node,
                     boost::property_tree::ptree &pt,
                     vector<vector<size_t> > & integer_variants,
                     vector<vector<pair<size_t,size_t> > > & uvw_restricted_edges)
{
  cerr << "# [info] round integer variable." << endl;
  integer_rounding_integer_variable(des, mesh, node, pt, integer_variants, uvw_restricted_edges);
  cerr << "# [info] round gaps." << endl;
  integer_rounding_gap(des, mesh, node, pt);
  return 0;
}
