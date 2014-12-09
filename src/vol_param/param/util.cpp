#include "util.h"
#include <set>
#include <boost/unordered_map.hpp>
#include <jtflib/mesh/util.h>
#include <jtflib/mesh/mesh.h>
#include "../../tet_mesh_sxx/tet_mesh_sxx.h"
#include <iostream>
#include <fstream>
#include <stdexcept>
#include "../../common/transition_type.h"
#include "../../common/util.h"
#include "../../equation_graph/equation_graph.h"
#include "../../common/vtk.h"

using namespace std;
using namespace zjucad::matrix;

int fix_sxx_tet_mapping(const zjucad::matrix::matrix<size_t> & uncut_tet,
                        zjucad::matrix::matrix<size_t> & cut_tet,
                        const std::vector<size_t> & cut_tet2tet,
                        const zjucad::matrix::matrix<size_t> & new_tet2orig_tet_uncut,
                        const zjucad::matrix::matrix<size_t> & new_tet2orig_tet_cut)
{
  itr_matrix<const size_t*> cut_tet2tet_mat(cut_tet2tet.size(),1, &cut_tet2tet[0]);

  return fix_sxx_tet_mapping(uncut_tet, cut_tet, cut_tet2tet_mat,
                             new_tet2orig_tet_uncut, new_tet2orig_tet_cut);
}

int fix_sxx_tet_mapping(const zjucad::matrix::matrix<size_t> & uncut_tet,
                        zjucad::matrix::matrix<size_t> & cut_tet,
                        const zjucad::matrix::matrix<size_t> & cut_tet2tet,
                        const zjucad::matrix::matrix<size_t> & new_tet2orig_tet_uncut,
                        const zjucad::matrix::matrix<size_t> & new_tet2orig_tet_cut)
{
  assert(cut_tet.size(2) == uncut_tet.size(2));
  set<size_t> A(new_tet2orig_tet_uncut.begin(), new_tet2orig_tet_uncut.end());
  set<size_t> B(new_tet2orig_tet_cut.begin(), new_tet2orig_tet_cut.end());

  assert(A == B);

  // first step tet order may be not the same, so I should reoder all tets
  {
    map<size_t,size_t> orig_tet2new_tet_idx_cut;
    for(size_t ti = 0; ti < new_tet2orig_tet_cut.size(); ++ti)
      orig_tet2new_tet_idx_cut[new_tet2orig_tet_cut[ti]] = ti;

    matrix<size_t> cut_tet_new = cut_tet;
    for(size_t ti = 0; ti < new_tet2orig_tet_uncut.size(); ++ti){
        const size_t & orig_idx = new_tet2orig_tet_uncut[ti];
        const auto orig2new_cut_it = orig_tet2new_tet_idx_cut.find(orig_idx);
        if(orig2new_cut_it == orig_tet2new_tet_idx_cut.end()){
            throw std::logic_error("strange can not find orig tet idx to new in cut");
          }
        cut_tet_new(colon(), ti) = cut_tet(colon(), orig2new_cut_it->second);
      }
    cut_tet = cut_tet_new;
  }

  matrix<size_t> one_cut2orig_tet(cut_tet.size(1),1);
  for(size_t ti = 0; ti < cut_tet.size(2); ++ti){
      for(size_t pi = 0; pi < cut_tet.size(1); ++pi){
          for(size_t pj = pi; pj < cut_tet.size(1); ++pj){
              if(cut_tet2tet[cut_tet(pj,ti)] == uncut_tet(pi,ti)){
                  swap(cut_tet(pi,ti),cut_tet(pj,ti));
                  break;
                }
            }
          one_cut2orig_tet[pi] = cut_tet2tet[cut_tet(pi,ti)];
        }

      if(norm(one_cut2orig_tet - uncut_tet(colon(),ti)) > 1e-6){
          cerr << "# orig_tet " << uncut_tet(colon(),ti) << endl;
          cerr << "# cut_tet " << cut_tet(colon(),ti) << endl;
          cerr << "# cut2orig_tet " << one_cut2orig_tet << endl;
          throw std::logic_error("one cut tet convert to orig tet fail.");
        }
    }

  return 0;
}

inline void swap_edge(pair<size_t,size_t> & edge)
{
  if(edge.first > edge.second)
    swap(edge.first, edge.second);
}

int update_inner_type_by_edge_collapsing(
    const pair<size_t,size_t> & one_edge,
    boost::unordered_map<pair<size_t,size_t>,size_t> & inner_type,
    const matrix<size_t> & tet,
    vector<bool> & tet_flag, // true means tet is valid, false means tet is removed
    const jtf::mesh::one_ring_tet_at_edge & ortae,
    const jtf::mesh::face2tet_adjacent &fa) // assume ortae is already sorted
{
  auto edge_it = ortae.e2t_.find(one_edge);
  if(edge_it == ortae.e2t_.end())
    edge_it = ortae.e2t_.find(make_pair(one_edge.second, one_edge.first));
  if(edge_it == ortae.e2t_.end())
    throw std::logic_error("can not find edge in one ring tet at edge.") ;

  const vector<size_t> & arounding_tets = edge_it->second;

  {// clean tet pair type among arounding_tets
    //assert(arounding_tets.front() == -1 && arounding_tets.back() == -1);// surface edge
    for(size_t i = 0; i < arounding_tets.size()-1; ++i){
        auto it = inner_type.find(make_pair(arounding_tets[i], arounding_tets[i+1]));
        if(it == inner_type.end()) continue;
        inner_type.erase(it);
        it = inner_type.find(make_pair(arounding_tets[i+1], arounding_tets[i]));
        assert(it != inner_type.end());
        inner_type.erase(it);
      }

  }
  // this arounding tets should be sorted
  set<size_t> arounding_tets_set(arounding_tets.begin(), arounding_tets.end());
  for(const auto & idx : arounding_tets_set){
      if(idx != -1) { // inner tets
          if(tet_flag[idx]){// two face of this tet will be merged
              vector<size_t> other_edges;
              for(size_t pi = 0; pi < tet.size(1); ++pi){
                  if(tet(pi, idx) != one_edge.first &&
                     tet(pi, idx) != one_edge.second)
                    other_edges.push_back(tet(pi, idx));
                }
              if(other_edges.size() != 2){
                  throw std::logic_error("strange, this tet does not contain given edge. ");
                }
              const size_t face_idx_0 = fa.get_face_idx(other_edges.front(),
                                                        other_edges.back(),
                                                        one_edge.first);
              const size_t face_idx_1 = fa.get_face_idx(other_edges.front(),
                                                        other_edges.back(),
                                                        one_edge.second);
              assert(face_idx_0 != -1 && face_idx_1 != -1);
              const pair<size_t,size_t> &tet_pair_0 = fa.face2tet_[face_idx_0];
              const pair<size_t,size_t> &tet_pair_1 = fa.face2tet_[face_idx_1];
              if(fa.is_outside_face(tet_pair_0) || fa.is_outside_face(tet_pair_1)){ // has face outside on surface
                  // this collapsion will remove inner face type across these faces
                  if(!fa.is_outside_face(tet_pair_0)){
                      auto it = inner_type.find(tet_pair_0);
                      if(it != inner_type.end()) inner_type.erase(it);
                      it = inner_type.find(make_pair(tet_pair_0.second, tet_pair_0.first));
                      if(it != inner_type.end()) inner_type.erase(it);
                    }
                  if(!fa.is_outside_face(tet_pair_1)){
                      auto it = inner_type.find(tet_pair_1);
                      if(it != inner_type.end()) inner_type.erase(it);
                      it = inner_type.find(make_pair(tet_pair_1.second, tet_pair_1.first));
                      if(it != inner_type.end()) inner_type.erase(it);
                    }
                }else{ // both faces are inner face, and the tet pair of them should be merged together.
                  // this tet seq contains three tets,
                  // middle is the tet need to be collapsed.
                  vector<size_t> tet_seq(3);
                  tet_seq[1] = idx;
                  tet_seq[0] = tet_pair_0.first + tet_pair_0.second - idx;
                  tet_seq[2] = tet_pair_1.first + tet_pair_1.second - idx;

                  matrix<double> rot = eye<double>(3);
                  auto it_0 = inner_type.find(make_pair(tet_seq[0], tet_seq[1]));
                  if(it_0  != inner_type.end()){
                      rot = temp(rot * type_transition2(it_0->second));
                      inner_type.erase(it_0);
                      it_0 = inner_type.find(make_pair(tet_seq[1], tet_seq[0]));
                      assert(it_0 != inner_type.end());
                      inner_type.erase(it_0);
                    }
                  auto it_1 = inner_type.find(make_pair(tet_seq[1], tet_seq[2]));
                  if(it_1 != inner_type.end()){
                      rot = temp(rot * type_transition2(it_1->second));
                      inner_type.erase(it_1);
                      it_1 = inner_type.find(make_pair(tet_seq[2], tet_seq[1]));
                      assert(it_1 != inner_type.end());
                      inner_type.erase(it_1);
                    }
                  const size_t type = type_transition1(rot);
                  if(type != TRIVIAL_TYPE){
                      inner_type[make_pair(tet_seq[0], tet_seq[2])] = type;
                      inner_type[make_pair(tet_seq[2], tet_seq[0])] = get_trans_type(type);
                    }
                }
              tet_flag[idx] = false;
            }
        }
    }
  return 0;
}

inline void orig_edge2cut_edge(
    const matrix<size_t> & cut_tet,
    const vector<size_t> & cut_tet2tet,
    map<pair<size_t,size_t>, set<pair<size_t,size_t> > > &orig_edge2cut_edges)
{
  for(size_t ti = 0;  ti < cut_tet.size(2); ++ti){
      for(size_t pi = 0; pi < cut_tet.size(1); ++pi){
          for(size_t pj = pi; pj < cut_tet.size(1); ++pj){
              pair<size_t,size_t> cut_edge(
                    cut_tet(pi, ti),
                    cut_tet((pj+1)%cut_tet.size(1),ti));// cut_edge
              pair<size_t,size_t> orig_edge(cut_tet2tet[cut_edge.first],
                  cut_tet2tet[cut_edge.second]);
              if(orig_edge.first > orig_edge.second){
                  swap(orig_edge.first, orig_edge.second);
                  swap(cut_edge.first, cut_edge.second);
                }
              orig_edge2cut_edges[orig_edge].insert(cut_edge);
            }
        }
    }
}
int collapse_small_angle_edges(
    const set<pair<size_t,size_t> > & edges_in_cut,
    const jtf::mesh::face2tet_adjacent & fa_uncut,
    matrix<size_t> & uncut_tet,
    matrix<double> & uncut_node,
    matrix<size_t> & cut_tet,
    matrix<double> & cut_node,
    matrix<double> & param_node,
    boost::unordered_map<pair<size_t,size_t>,size_t> & inner_type,
    matrix<matrix<double> > & frame,
    std::vector<std::vector<size_t> > * eqn_idx ,
    std::vector<std::vector<double> > * eqn_coeff)
{
  assert(uncut_tet.size(2) == cut_tet.size(2));

  vector<size_t> cut_tet2tet;
  map<size_t, set<size_t> > orig_point2new_point;
  {
    cut_tet2tet.resize(max(cut_tet) + 1);
    itr_matrix<size_t*> cut_tet2tet_mat(cut_tet2tet.size(), 1, &cut_tet2tet[0]);
    cut_tet2tet_mat(cut_tet) = uncut_tet(colon());
    for(size_t pi = 0; pi < cut_tet2tet_mat.size(); ++pi){
        orig_point2new_point[cut_tet2tet[pi]].insert(pi);
      }
  }

  assert(max(cut_tet) + 1 == cut_node.size(2)); // cut_node number shoule be max(cut_tet) +1

  vector<bool> tet_flag(uncut_tet.size(2), true); // flag to indicate whether a tet is valid not removed

  jtf::mesh::one_ring_tet_at_edge ortae;
  ortae.add_tets(uncut_tet, fa_uncut);
  ortae.sort_into_loop(uncut_tet, uncut_node);

  // construct orig edge to cut_edge mapping
  map<pair<size_t,size_t>, set<pair<size_t,size_t> > > orig_edge2cut_edges;

  orig_edge2cut_edge(cut_tet, cut_tet2tet, orig_edge2cut_edges);

  assert(cut_tet.size(2) == uncut_tet.size(2));
  //TODO: need to speed up
  sxx::tet_mesh stm_cut,stm_uncut,stm_param;
  matrix<size_t> param_tet = cut_tet;
  stm_cut.create_tetmesh(cut_node, cut_tet);
  stm_param.create_tetmesh(param_node, param_tet);
  stm_uncut.create_tetmesh(uncut_node, uncut_tet);

  size_t new_point_idx_in_orig_tet = uncut_node.size(2) ;

  bool real_collapsed = false;
  jtf::mesh::meshes cut_tet_mesh, param_tet_mesh; // for backup
  cut_tet_mesh.mesh_ = cut_tet;
  cut_tet_mesh.node_ = cut_node;

  param_tet_mesh.mesh_ = cut_tet;
  param_tet_mesh.node_ = param_node;

  // it store the mapping in a sequence
  map<size_t,size_t> new_tet_idx2orig_tet_idx_seq;
  for(size_t i = 0; i < cut_tet.size(2); ++i) new_tet_idx2orig_tet_idx_seq[i] = i;

  bool is_reloaed = false;
  //  pair<bool,bool> is_orig_cut_tet_updated(false,false); // first: is cut_tet_mesh modified, second: is modified_cut_tet_mesh reloaded
  for(const auto &one_edge : edges_in_cut){
      pair<size_t,size_t> orig_one_edge(cut_tet2tet[one_edge.first],
          cut_tet2tet[one_edge.second]);
      //cerr << "collapse uncut_edge " << orig_one_edge.first << " " << orig_one_edge.second << endl;

      swap_edge(orig_one_edge);

      vector<size_t> cut_tet2tet_bkp = cut_tet2tet;
      const auto it = orig_edge2cut_edges.find(orig_one_edge);
      if(it == orig_edge2cut_edges.end()) {
          // can not find this orig_edge, that's because cut_tet2tet is updated
          // and point mapping is changed to new inserted one, thus I
          // ignore this edge just now.
          continue;
          //throw std::logic_error("strange can not find orig edge in orig2cut edge mapping.");
        };

      const set<pair<size_t,size_t> > & cut_edges = it->second;
      bool collapsed = true;

      for(const auto & one_cut_edge : cut_edges){

          int rtn_cut = stm_cut.collapse_edge(one_cut_edge);
          int rtn_param = stm_param.collapse_edge(one_cut_edge);
          assert(rtn_cut == rtn_param);

          if(rtn_cut == -1) {
              collapsed = false;
              break;
            }

          cut_tet2tet.push_back(new_point_idx_in_orig_tet);

          const auto pit_first = orig_point2new_point.find(cut_tet2tet[one_cut_edge.first]);
          const auto pit_second = orig_point2new_point.find(cut_tet2tet[one_cut_edge.second]);

          if(pit_first != orig_point2new_point.end()){ // this point has not been modified before
              for(const auto & idx_first : pit_first->second){
                  cut_tet2tet[idx_first] = new_point_idx_in_orig_tet;
                }
            }

          if(pit_second != orig_point2new_point.end()){
              for(const auto & idx_second : pit_second->second){
                  cut_tet2tet[idx_second] = new_point_idx_in_orig_tet;
                }
            }
        }

      if(collapsed){
          int rtn = stm_uncut.collapse_edge(orig_one_edge); // this collapse will add new vertex

          if(rtn != -1){ // can collapse in uncut mesh
              update_inner_type_by_edge_collapsing(orig_one_edge, inner_type, uncut_tet,
                                                   tet_flag, ortae, fa_uncut);
              real_collapsed = true;
              ++new_point_idx_in_orig_tet;

              stm_cut.write_tetmesh_to_matrix(cut_tet_mesh.node_, cut_tet_mesh.mesh_);
              stm_param.write_tetmesh_to_matrix(param_tet_mesh.node_, param_tet_mesh.mesh_);

              matrix<size_t> n2o_mapping;
              stm_cut.get_tet2orginal_index(cut_tet_mesh.mesh_, n2o_mapping);

              map<size_t,size_t> n2o_mapping_map;
              for(size_t noi = 0 ; noi < n2o_mapping.size(); ++noi){
                  const auto it = new_tet_idx2orig_tet_idx_seq.find(n2o_mapping[noi]);
                  assert(it != new_tet_idx2orig_tet_idx_seq.end());
                  n2o_mapping_map[noi] = it->second;
                }

              new_tet_idx2orig_tet_idx_seq = n2o_mapping_map;
            }else{
              cut_tet2tet = cut_tet2tet_bkp;
            }

        }else{
          cut_tet2tet = cut_tet2tet_bkp;
        }

      stm_cut.clear();
      stm_cut.create_tetmesh(cut_tet_mesh.node_, cut_tet_mesh.mesh_);

      stm_param.clear();
      stm_param.create_tetmesh(param_tet_mesh.node_, param_tet_mesh.mesh_);
    }

  if(!real_collapsed) return 0;

  const size_t orig_cut_tet_number = cut_tet.size(2);
  const size_t orig_uncut_tet_number = uncut_tet.size(2);

  stm_cut.write_tetmesh_to_matrix(cut_node, cut_tet);
  stm_param.write_tetmesh_to_matrix(param_node, param_tet);
  stm_uncut.write_tetmesh_to_matrix(uncut_node, uncut_tet);

  matrix<size_t> new_tet2orig_tet_uncut, new_tet2orig_tet_cut;
  stm_uncut.get_tet2orginal_index(uncut_tet, new_tet2orig_tet_uncut);
  stm_cut.get_tet2orginal_index(cut_tet, new_tet2orig_tet_cut);


  for(size_t noi = 0 ; noi < new_tet2orig_tet_cut.size(); ++noi){
      const auto it = new_tet_idx2orig_tet_idx_seq.find(noi);
      assert(it != new_tet_idx2orig_tet_idx_seq.end());
      new_tet2orig_tet_cut[noi] = it->second;
    }

  assert(max(new_tet2orig_tet_cut) < orig_cut_tet_number);
  assert(max(new_tet2orig_tet_uncut) < orig_uncut_tet_number);

  fix_sxx_tet_mapping(uncut_tet, param_tet, cut_tet2tet,
                      new_tet2orig_tet_uncut, new_tet2orig_tet_cut );

  cut_tet = param_tet;


  if(frame.size() > 0) { // in case polycube configuration
      matrix<matrix<double> > new_frame(uncut_tet.size(2),1);
      for(size_t ti = 0; ti < uncut_tet.size(2); ++ti){
          new_frame[ti] = frame[new_tet2orig_tet_uncut[ti]];
        }
      frame = new_frame;
    }

  remove_extra_node(cut_tet, cut_node);
  remove_extra_node(uncut_tet, uncut_node);
  remove_extra_node(param_tet, param_node);


  assert(max(cut_tet)+1 == cut_node.size(2));
  assert(max(uncut_tet)+1 == uncut_node.size(2));
  assert(max(param_tet)+1 == param_node.size(2));

  map<size_t,size_t> orig_tet2new_tet;
  for(size_t i = 0; i < new_tet2orig_tet_uncut.size(); ++i)
    orig_tet2new_tet[new_tet2orig_tet_uncut[i]] = i;

  boost::unordered_map<pair<size_t,size_t>,size_t> updated_inner_type;
  for(const auto & one_face : inner_type){
      const auto t0 = orig_tet2new_tet.find(one_face.first.first);
      const auto t1 = orig_tet2new_tet.find(one_face.first.second);

      if(t0 == orig_tet2new_tet.end() ||
         t1 == orig_tet2new_tet.end()){
          cerr << endl;
        }else
        updated_inner_type[make_pair(t0->second, t1->second)] = one_face.second;
    }

  inner_type = updated_inner_type;

  return 0;
}

int split_large_angle_edges(
    const std::set<std::pair<size_t,size_t> > & edges_in_cut,
    const jtf::mesh::face2tet_adjacent & fa_uncut,
    zjucad::matrix::matrix<size_t> & uncut_tet,\
    zjucad::matrix::matrix<double> & uncut_node,
    zjucad::matrix::matrix<size_t> & cut_tet,
    zjucad::matrix::matrix<double> & cut_node,
    zjucad::matrix::matrix<double> & param_node,
    boost::unordered_map<std::pair<size_t,size_t>,size_t> & inner_type,
    std::vector<std::vector<size_t> > * eqn_idx,
    std::vector<std::vector<double> > * eqn_coeff)
{
  assert(cut_tet.size(2) == uncut_tet.size(2));
  vector<size_t> cut_tet2tet(max(cut_tet)+1);
  {
    itr_matrix<size_t*> cut_tet2tet_mat(cut_tet2tet.size(),1, &cut_tet2tet[0]);
    cut_tet2tet_mat(cut_tet) = uncut_tet(colon());
  }

  // construct orig edge to cut_edge mapping
  map<pair<size_t,size_t>, set<pair<size_t,size_t> > > orig_edge2cut_edges;

  orig_edge2cut_edge(cut_tet, cut_tet2tet, orig_edge2cut_edges);

  sxx::tet_mesh stm_uncut, stm_cut, stm_param;
  matrix<size_t> param_tet = cut_tet;
  stm_uncut.create_tetmesh(uncut_node, uncut_tet);
  stm_cut.create_tetmesh(cut_tet, cut_node);
  stm_param.create_tetmesh(param_tet, param_node);

  size_t point_number_in_uncut_tet = uncut_node.size(2);
  for(const auto & one_cut_edge : edges_in_cut){
      pair<size_t,size_t> orig_edge(cut_tet2tet[one_cut_edge.first],
          cut_tet2tet[one_cut_edge.second]);
      swap_edge(orig_edge);

      stm_uncut.split_edge(orig_edge);

      const auto it  = orig_edge2cut_edges.find(orig_edge);
      if(it == orig_edge2cut_edges.end()){
          throw std::logic_error("strange can not find orig edge in orig_edge2cut_edges");
        }

      const set<pair<size_t,size_t> > & cut_edges = it->second;
      for(const auto & one_cut_edge : cut_edges){
          stm_cut.split_edge(one_cut_edge);
          stm_param.split_edge(one_cut_edge);
          cut_tet2tet.push_back(point_number_in_uncut_tet);
        }
      ++point_number_in_uncut_tet;

      //      update_inner_type_by_edge_splitting(orig_edge, inner_type, uncut_tet,
      //                                           ortae, fa_uncut);
    }


  return 0;
}

int remove_degenerated_face_of_tet(
    const matrix<size_t> & orig_face_in_cut,
    const jtf::mesh::face2tet_adjacent & fa_uncut,
    zjucad::matrix::matrix<size_t> & uncut_tet,
    zjucad::matrix::matrix<double> & uncut_node,
    zjucad::matrix::matrix<size_t> & cut_tet,
    zjucad::matrix::matrix<double> & cut_node,
    zjucad::matrix::matrix<double> & param_node,
    boost::unordered_map<std::pair<size_t,size_t>,size_t> & inner_type,
    zjucad::matrix::matrix<zjucad::matrix::matrix<double> > & frame,
    const jtf::mesh::edge2cell_adjacent & ea_cut,
    const string remesh_strategy)
{
  set<pair<size_t,size_t> > edges_to_be_collapsed;
  vector<double> angle;

  matrix<size_t> cut_tet2tet(max(cut_tet)+1);
  cut_tet2tet(cut_tet) = uncut_tet(colon());
  //////////////////////////////////////////////////////////////////////////////
  /// collapse edges which make volume of a tet to be zero
  ///

  if(remesh_strategy == "vol"){
      matrix<double> vol(cut_tet.size(2),1);

      for(size_t ti = 0; ti < cut_tet.size(2); ++ti){
          vol[ti] = jtf::mesh::cal_tet_vol(cut_node(colon(),cut_tet(colon(),ti)));
        }

      const double avg_vol = std::accumulate(vol.begin(), vol.end(), 0.0)/cut_tet.size(2);

      vector<pair<double,size_t> > edges;
      for(size_t ti = 0; ti < cut_tet.size(2); ++ti){
          if(fabs(vol[ti]) < fabs(avg_vol)/100){ // this tet is already degenerated
              edges.clear();
              for(size_t pi = 0; pi < cut_tet.size(1); ++pi){
                  for(size_t pj = pi + 1; pj < cut_tet.size(1); ++pj){
                      edges.push_back(
                            make_pair(norm(cut_node(colon(), cut_tet(pi, ti))
                                           - cut_node(colon(), cut_tet(pj%cut_tet.size(1), ti))),pi));
                    }
                }
              sort(edges.begin(), edges.end());
              pair<size_t,size_t> one_edge(cut_tet(edges.front().second,ti),
                                           cut_tet((edges.front().second+1)%cut_tet.size(1),ti));

              one_edge.first = cut_tet2tet[one_edge.first];
              one_edge.second = cut_tet2tet[one_edge.second];
              if(one_edge.first > one_edge.second) swap(one_edge.first, one_edge.second);
              edges_to_be_collapsed.insert(one_edge);
            }
        }
    }

  //////////////////////////////////////////////////////////////////////////////
  /// collapse edge based on surface triangle angle
  //////////////////////////////////////////////////////////////////////////////
  if(remesh_strategy == "surface"){
      const double  small_angle = 5;
      for(size_t fi = 0; fi < orig_face_in_cut.size(2); ++fi){
          jtf::mesh::cal_face_angle(orig_face_in_cut(colon(),fi), param_node, angle);
          if(*min_element(angle.begin(), angle.end()) < small_angle){
              size_t min_angle_point = min_element(angle.begin(), angle.end()) - angle.begin();
              pair<size_t,size_t> one_edge(orig_face_in_cut((min_angle_point+1)%3, fi),
                                           orig_face_in_cut((min_angle_point+2)%3, fi));
              const size_t edge_idx = ea_cut.get_edge_idx(one_edge.first, one_edge.second);
              assert(edge_idx != -1);
              if(ea_cut.is_boundary_edge(ea_cut.edge2cell_[edge_idx])) continue; // if it's a boundary edge, I do not collapse it.
              swap_edge(one_edge);

              edges_to_be_collapsed.insert(one_edge);
            }
        }
    }

  cerr << "# [info] edges need to be collapsed: " << edges_to_be_collapsed.size() << endl;

  if(edges_to_be_collapsed.size()){
      collapse_small_angle_edges(edges_to_be_collapsed, fa_uncut, uncut_tet, uncut_node, cut_tet, cut_node,
                                 param_node, inner_type, frame, 0, 0);
    }

  //  if(eqn_idx != 0 && eqn_coeff != 0){ // update node_equation
  //      matrix<size_t> cut_tet2tet(max(cut_tet)+1);
  //      cut_tet2tet(cut_tet) = uncut_tet(colon());
  //      boost::unordered_map<size_t,size_t> surface_type; // empty surface type
  //      unique_ptr<transition_elimination> te(
  //            transition_elimination::create(
  //              uncut_tet, cut_tet, cut_tet2tet, uncut_node, inner_type, surface_type, true));
  //      if(!te.get()){
  //          throw std::logic_error("can not build transition elimination");
  //        }
  //      te->get_equation(*eqn_idx, *eqn_coeff);


  //    }
  return edges_to_be_collapsed.size();
}


int load_integer_groups(const char * filename,
                        vector<vector<size_t> > & integer_v_group)
{
  ifstream ifs(filename);
  if(ifs.fail()){
      cerr << "# [error] can not open integer group file." << endl;
      return __LINE__;
    }
  vector<size_t> one_group;
  int integer_or_not; // 1: integer, 0: not integer
  size_t temp, v_size;
  string str_temp;
  while(!ifs.eof()){
      one_group.clear();
      ifs >> str_temp >> temp >> v_size >> integer_or_not;
      if(ifs.eof()) break;
      for(size_t i = 0; i < v_size; ++i){
          ifs >> temp;
          if(integer_or_not)
            one_group.push_back(temp);
        }
      if(integer_or_not)
        integer_v_group.push_back(one_group);
    }
  return 0;
}


int load_restricted_path(
    const char * filename,
    std::vector<std::vector<std::pair<size_t,size_t> > > & restricted_edges)
{
  ifstream ifs(filename);
  if(ifs.fail()){
      cerr << "# [error] can not load restricted_edges." << endl;
      return __LINE__;
    }
  restricted_edges.clear();

  restricted_edges.resize(3);

  size_t edge_num;
  ifs >> edge_num;
  std::tuple<size_t,size_t,size_t> one_edge;
  size_t segments_num;
  for(size_t ei = 0; ei < edge_num; ++ei){
      ifs >> segments_num;
      for(size_t ei = 0; ei < segments_num; ++ei){
          // only record first and last point of each pathl
          if(ei == 0)
            ifs >> std::get<0>(one_edge) >> std::get<1>(one_edge) >> std::get<2>(one_edge) ;
          else
            ifs >> std::get<1>(one_edge) >> std::get<1>(one_edge) >> std::get<2>(one_edge) ;
        }
      restricted_edges[std::get<2>(one_edge)].push_back(
            std::make_pair( std::get<0>(one_edge), std::get<1>(one_edge)));
    }

  return 0;
}


double get_minimal_distance(const vector<vector<size_t> > & integer_variants,
                            const matrix<double> & node,
                            const double threshold)
{
  vector<vector<double> > uvw(3);
  for(size_t i = 0; i < integer_variants.size(); ++i){
      const size_t di = integer_variants[i].front()%3;
      const size_t vi = integer_variants[i].front();
      uvw[di].push_back(node[vi]);
    }
  vector<double> dis;
  for(size_t i = 0; i < uvw.size(); ++i){
      if(uvw[i].empty()) continue;
      sort(uvw[i].begin(), uvw[i].end());
      for(size_t j = 0 ;j < uvw[i].size()-1; ++j){
          const double distance = fabs(uvw[i][j]-uvw[i][j+1]);
          if( distance > threshold) dis.push_back(distance);
        }
    }
  sort(dis.begin(), dis.end());
  return dis.front();
}
