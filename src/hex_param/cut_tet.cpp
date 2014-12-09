#include "topology_operation.h"
#include "hex_param.h"
#include "cut_tet.h"
#include "find_singularities.h"
#include "normalize_cut.h"
#include "topology_analysis.h"
#include "../common/graph.h"
#include "../common/zyz.h"
#include "../common/transition.h"
#include "../common/transition_type.h"
#include "../common/vtk.h"
#include "../common/IO.h"
#include "../common/visualize_tool.h"
#include <jtflib/util/container_operation.h>
#include "common.h"

extern "C" {
#include "../spherical_harmonics/rot_cubic_f_SH.h"
}

#include <iostream>
#include <stack>
#include <fstream>
#include <ctime>
#include <limits>
#include <numeric>

#include <zjucad/matrix/io.h>
#include <zjucad/matrix/itr_matrix.h>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/kruskal_min_spanning_tree.hpp>
#include <boost/graph/dijkstra_shortest_paths.hpp>
#include <boost/graph/graph_traits.hpp>

using namespace std;
using namespace zjucad::matrix;
using namespace jtf::mesh;

int tetmesh_cutter::glue_tet(size_t a, size_t b, bool is_span_tree_edge)
{
  if(a == b) return 0;
  size_t face[3], ab[2] = {a, b};
  common_tet_face(&tm_.tetmesh_.mesh_(0, a), &tm_.tetmesh_.mesh_(0, b), face);
  vector<pair<size_t, size_t> > merge_node; // merge second in b to first in a;
  merge_node.reserve(4);
  for(int i = 0; i < 3; ++i) { // unify the node idx in b to a
      int in_ab[2];
      for(int j = 0; j < 2; ++j) {
          for(in_ab[j] = 0; in_ab[j] < 4; ++in_ab[j]) {
              if(tm_.tetmesh_.mesh_(in_ab[j], ab[j]) == face[i]) break;
            }
        }
      assert(in_ab[0] < 4 && in_ab[1] < 4);
      const size_t unify_to = cut_tm_.mesh_(in_ab[0], a);
      const size_t unify_from = cut_tm_.mesh_(in_ab[1], b);
      if(unify_from == unify_to)
        continue;
      merge_node.push_back(make_pair(unify_to, unify_from));
    }
  if(merge_node.empty()) // already glued
    return 1;

  if(is_span_tree_edge && merge_node.size() != 3) // error spanning_tree edge
    cerr << "# error: merge_node.size() != 3 in spanning_tree." << endl;
  if(!is_span_tree_edge && merge_node.size() != 1) // cause loop
    return 2;

  for(size_t ai = 0; ai < merge_node.size(); ++ai) {
      set<matrixst::value_type *> &node_set = acc_table_[merge_node[ai].second];
      assert(!node_set.empty());
      for(set<matrixst::value_type *>::const_iterator iter = node_set.begin();
          iter != node_set.end(); ++iter) {
          assert(*(*iter) == merge_node[ai].second);
          *(*iter) = merge_node[ai].first;
          acc_table_[merge_node[ai].first].insert(*iter);
        }
      node_set.clear();
    }

  ec_.node_num -= merge_node.size();
  ec_.edge_num -= 3;
  ec_.face_num -= 2;
  return 0;
}

void tetmesh_cutter::init()
{
  tet2dual_graph(tm_.tetmesh_.mesh_, fg_, *tm_.fa_);
  cut_tm_.mesh_.resize(4, tm_.tetmesh_.mesh_.size(2));

  cut_tm_.mesh_(colon()) = colon(0, cut_tm_.mesh_.size()-1);
  acc_table_.clear();
  acc_table_.resize(cut_tm_.mesh_.size());

  for(size_t i = 0; i < acc_table_.size(); ++i)
    acc_table_[i].insert(&cut_tm_.mesh_[i]);

  ec_.node_num = cut_tm_.mesh_.size(2)*4;
  ec_.edge_num = cut_tm_.mesh_.size(2)*6;
  ec_.face_num = cut_tm_.mesh_.size(2)*4;

  cerr << "# Euler_characteristic before glue: " << ec_() << " "
       << ec_.node_num << " " << ec_.edge_num << " " << ec_.face_num << endl;
}

void tetmesh_cutter::update_cut_node()
{
  map<size_t,size_t> mapping;
  remove_extra_node(cut_tm_.mesh_, mapping);
  matrix<size_t> cut_tet2tet(max(cut_tm_.mesh_)+1,1);
  cut_tet2tet(cut_tm_.mesh_) = tm_.tetmesh_.mesh_(colon());
  cut_tm_.node_.resize(3, cut_tet2tet.size());;
  for(size_t pi = 0; pi < cut_tet2tet.size(); ++pi){
      cut_tm_.node_(colon(),pi) = tm_.tetmesh_.node_(colon(), cut_tet2tet[pi]);
    }
}

void tetmesh_cutter::cut(matrix<double> &zyz, bool only_mst)
{
  init();

  matrixd tet_cf_sh(9, zyz.size(2));
  for(size_t fi = 0; fi < tet_cf_sh.size(2); ++fi){
      calc_rot_cubic_f_sh_(&tet_cf_sh(0, fi), &zyz(0,fi));
    }

  std::vector<double> err(fg_.edge_num());
  for(size_t ni = 0; ni < fg_.node_num(); ++ni) {
      for(size_t ei = fg_.ptr_[ni]; ei < fg_.ptr_[ni+1]; ++ei)
        err[ei] = norm(tet_cf_sh(colon(), ni) - tet_cf_sh(colon(), fg_.idx_[ei]));
    }

  std::vector<size_t> spanning_tree;
  matrix<matrixd > frame;
  {
    fix_graph spanning_tree_graph;
    build_minimum_spanning_tree(fg_, err, spanning_tree);
    subgraph(fg_, spanning_tree, spanning_tree_graph);
    cerr << "# [info] graph node: " << spanning_tree_graph.node_num()
         << ", edges in spanning_tree: " << spanning_tree_graph.edge_num()
         << endl;

    matrixd orient_err(spanning_tree_graph.edge_num());
    orient_frame(zyz, frame, spanning_tree_graph, orient_err);
    frame2zyz(frame, zyz);
  }

  for(size_t ei = 0; ei < spanning_tree.size(); ++ei) {
      std::pair<node_idx, node_idx> edge = fg_.get_edge(spanning_tree[ei]);
      if(glue_tet(edge.first, edge.second, true) > 1) {
          cerr << "# [error] glue according spanning_tree error." << endl;
        }
    }
  cerr << "# [info] Euler_characteristic after glue MST: " << ec_() << " "
       << ec_.node_num << " " << ec_.edge_num << " " << ec_.face_num << endl;
  update_cut_node();

  if(only_mst) return;

  vector<deque<pair<size_t,size_t> > > singularity_chain;
  vector<deque<size_t> > singularity_type;

  singularity_extractor se(tm_);
  se.extract(frame, singularity_chain, singularity_type);

  jtf::mesh::dyn_one_ring_face_at_point dorfap;

  {
    unique_ptr<jtf::mesh::face2tet_adjacent> fa_test(jtf::mesh::face2tet_adjacent::create(cut_tm_.mesh_));
    matrixst outside_face_cut_tet;
    get_outside_face(*fa_test, outside_face_cut_tet);
    unique_ptr<jtf::mesh::edge2cell_adjacent> ea(
          jtf::mesh::edge2cell_adjacent::create(outside_face_cut_tet));
    if(!ea.get()){
        cerr << "# [error] non-manifold surface after MST." << endl;
        return ;
      }
    dorfap.add_all_faces(outside_face_cut_tet);
  }

  cerr << "# [info] begin to glue tets avoid singularity crease." << endl;
  glue_tet_avoid_singularity_crease(tm_.tetmesh_.mesh_, frame, *tm_.fa_,
                                    singularity_chain, tm_.ortae_, fg_,
                                    dorfap,cut_tm_.mesh_,ec_,acc_table_);
  update_cut_node();
}

void tetmesh_cutter::cut(const boost::unordered_map<std::pair<size_t,size_t>,size_t> & inner_type)
{
  init();

  typedef boost::unordered_map<pair<size_t,size_t>,size_t>::const_iterator mcit;

  std::vector<double> err(fg_.edge_num());
  for(size_t ni = 0; ni < fg_.node_num(); ++ni) {
      for(size_t ei = fg_.ptr_[ni]; ei < fg_.ptr_[ni+1]; ++ei){
          mcit it = inner_type.find(make_pair(ni,fg_.idx_[ei]));
          if(it == inner_type.end()) err[ei] = 0;
          else {
              if(is_trivial_type(it->second)) err[ei] = 0.;
              else err[ei] = 1.;
            }
        }
    }

  std::vector<size_t> spanning_tree;
  {
    fix_graph spanning_tree_graph;
    build_minimum_spanning_tree(fg_, err, spanning_tree);
    subgraph(fg_, spanning_tree, spanning_tree_graph);
    cerr << "# graph node: " << spanning_tree_graph.node_num()
         << ", edges in spanning_tree: " << spanning_tree_graph.edge_num() << endl;
  }

  for(size_t ei = 0; ei < spanning_tree.size(); ++ei) {
      std::pair<node_idx, node_idx> edge = fg_.get_edge(spanning_tree[ei]);
      if(glue_tet(edge.first, edge.second,true) > 1) {
          cerr << "glue according spanning_tree error." << endl;
        }
    }
  cerr << "# Euler_characteristic after glue MST: " << ec_() << " "
       << ec_.node_num << " " << ec_.edge_num << " " << ec_.face_num << endl;
}

bool is_singularity_edge_glue_inside_cut_tet(
    const vector<pair<size_t,size_t> > & merge_node,
    const matrixst & tet,
    const matrixst & cut_tet,
    const jtf::mesh::one_ring_tet_at_edge & ortae,
    const pair<size_t,size_t> & singularity_edge);

bool is_contain_singularity_edge_and_surface_edge(
    const vector<size_t> &face,
    const deque<pair<size_t,size_t> > &chain,
    const jtf::mesh::edge2cell_adjacent &ea,
    set<pair<size_t,size_t> > *singularity_edge = 0);

void build_minimum_spanning_tree(
    const fix_graph &fg, std::vector<double> &weights,
    std::vector<size_t> &edge_list)
{
  using namespace boost;
  typedef adjacency_list < vecS, vecS, undirectedS,
      no_property, property < edge_weight_t, double > > Graph;
  typedef graph_traits < Graph >::edge_descriptor Edge;
  typedef graph_traits < Graph >::vertex_descriptor Vertex;
  typedef std::pair<int, int> E;

  std::vector<E> edge_array(fg.edge_num());
  for(int ni = 0, ei = 0; ni < fg.node_num(); ++ni) {
      for(vector<node_idx>::const_iterator iter = fg.begin(ni);
          iter != fg.end(ni); ++iter, ++ei) {
          if(ni < *iter)
            edge_array[ei] = E(ni, *iter);
        }
    }

  Graph g(edge_array.begin(), edge_array.end(), &weights[0], fg.node_num());

  property_map < Graph, edge_weight_t >::type weight = get(edge_weight, g);

  std::vector < Edge > spanning_tree0;
  kruskal_minimum_spanning_tree(g, std::back_inserter(spanning_tree0));

  edge_list.reserve(spanning_tree0.size()*2);
  edge_list.clear();
  for (std::vector < Edge >::iterator ei = spanning_tree0.begin();
       ei != spanning_tree0.end(); ++ei) {
      edge_list.push_back(fg.find_edge(source(*ei, g), target(*ei, g)));
      edge_list.push_back(fg.find_edge(target(*ei, g), source(*ei, g)));
    }
}

int tet2dual_graph(const matrixst &tet, fix_graph &fg, const jtf::mesh::face2tet_adjacent &fa)
{
  matrixst ptr, idx;
  if(tet2tet_adj_matrix(tet, ptr, idx, &fa))
    return __LINE__;
  fg.ptr_.resize(ptr.size());
  fg.idx_.resize(idx.size());
  std::copy(ptr.begin(), ptr.end(), fg.ptr_.begin());
  std::copy(idx.begin(), idx.end(), fg.idx_.begin());
  return 0;
}

void subgraph(const fix_graph &graph,
              const std::vector<size_t> &edge_list,
              fix_graph &sub)
{
  cerr << "# subgraph: edge: " << edge_list.size() << endl;
  dyn_graph dg(graph.node_num());
  for(size_t ei = 0; ei < edge_list.size(); ++ei) {
      if(!fix_graph::is_valid_idx(edge_list[ei])) {
          cerr << "error in edge_list." << endl;
        }
      std::pair<node_idx, node_idx> edge = graph.get_edge(edge_list[ei]);
      if(!fix_graph::is_valid_idx(edge.first)
         || !fix_graph::is_valid_idx(edge.second)) {
          cerr << "error in edge." << endl;
        }
      dg[edge.first].insert(edge.second);
      dg[edge.second].insert(edge.first);
    }
  convert(dg, sub);
}


/**
 * @return 0: successfully glue tet a and b, 1: already glued, 2:
 * topology changed
 *
 * NOTICE, when there're 2 common nodes, the edge of decrease may be 2
 * or 3, depends on whether the edge has been fully surrounded.
 * However, since we restrict the surface of the jtf::mesh::meshes is manifold,
 * so edge_num alwasy decreases by 3.
 */
int glue_tet(const matrixst &input_tet,
             matrixst &output_tet,
             size_t a, size_t b, Euler_characteristic &ec,
             vector<set<matrixst::value_type *> > &acc_table,
             bool is_span_tree_edge)
{
  assert(a != b);
  size_t face[3], ab[2] = {a, b};
  common_tet_face(&input_tet(0, a), &input_tet(0, b), face);
  vector<pair<size_t, size_t> > merge_node; // merge second in b to first in a;
  merge_node.reserve(4);
  for(int i = 0; i < 3; ++i) { // unify the node idx in b to a
      int in_ab[2];
      for(int j = 0; j < 2; ++j) {
          for(in_ab[j] = 0; in_ab[j] < 4; ++in_ab[j]) {
              if(input_tet(in_ab[j], ab[j]) == face[i]) break;
            }
        }
      assert(in_ab[0] < 4 && in_ab[1] < 4);
      const size_t unify_to = output_tet(in_ab[0], a);
      const size_t unify_from = output_tet(in_ab[1], b);
      if(unify_from == unify_to)
        continue;
      merge_node.push_back(make_pair(unify_to, unify_from));
    }
  if(merge_node.empty()) // already glued
    return 1;

  if(is_span_tree_edge && merge_node.size() != 3) // error spanning_tree edge
    cerr << "# error: merge_node.size() != 3 in spanning_tree." << endl;
  if(!is_span_tree_edge && merge_node.size() != 1) // cause loop
    return 2;

  for(size_t ai = 0; ai < merge_node.size(); ++ai) {
      set<matrixst::value_type *> &node_set = acc_table[merge_node[ai].second];
      assert(!node_set.empty());
      for(set<matrixst::value_type *>::const_iterator iter = node_set.begin();
          iter != node_set.end(); ++iter) {
          assert(*(*iter) == merge_node[ai].second);
          *(*iter) = merge_node[ai].first;
          acc_table[merge_node[ai].first].insert(*iter);
        }
      node_set.clear();
    }

  ec.node_num -= merge_node.size();
  ec.edge_num -= 3;
  ec.face_num -= 2;
  return 0;
}

// this function glue two faces only if they are not on span_tree
int glue_tet_specail(
    const matrixst &input_tet,
    matrixst &output_tet,
    size_t a, size_t b, Euler_characteristic &ec,
    vector<set<matrixst::value_type *> > &acc_table,
    bool is_span_tree_edge)
{
  assert(a != b);
  size_t face[3], ab[2] = {a, b};
  common_tet_face(&input_tet(0, a), &input_tet(0, b), face);
  vector<pair<size_t, size_t> > merge_node; // merge second in b to first in a;
  merge_node.reserve(4);
  for(int i = 0; i < 3; ++i) { // unify the node idx in b to a
      int in_ab[2];
      for(int j = 0; j < 2; ++j) {
          for(in_ab[j] = 0; in_ab[j] < 4; ++in_ab[j]) {
              if(input_tet(in_ab[j], ab[j]) == face[i]) break;
            }
        }
      assert(in_ab[0] < 4 && in_ab[1] < 4);
      const size_t unify_to = output_tet(in_ab[0], a);
      const size_t unify_from = output_tet(in_ab[1], b);
      if(unify_from == unify_to)
        continue;
      merge_node.push_back(make_pair(unify_to, unify_from));
    }
  if(merge_node.empty()) // already glued
    return 1;

  if(is_span_tree_edge && merge_node.size() != 3) // error spanning_tree edge
    cerr << "# error: merge_node.size() != 3 in spanning_tree." << endl;

  for(size_t ai = 0; ai < merge_node.size(); ++ai) {
      set<matrixst::value_type *> &node_set = acc_table[merge_node[ai].second];
      assert(!node_set.empty());
      for(set<matrixst::value_type *>::const_iterator iter = node_set.begin();
          iter != node_set.end(); ++iter) {
          assert(*(*iter) == merge_node[ai].second);
          *(*iter) = merge_node[ai].first;
          acc_table[merge_node[ai].first].insert(*iter);
        }
      node_set.clear();
    }

  ec.node_num -= merge_node.size();
  ec.edge_num -= 3;
  ec.face_num -= 2;
  return 0;
}

int glue_each_tet_pair_avoid_singularity_edge(
    const matrixst & original_tet,
    matrixst & output_cut_tet,
    const size_t & a, const size_t & b,
    Euler_characteristic & ec,
    vector<set<size_t* > > & acc_table,
    bool is_span_tree_edge,
    const set<pair<size_t,size_t> > & singularity_edges,
    const jtf::mesh::one_ring_tet_at_edge & ortae,
    jtf::mesh::dyn_one_ring_face_at_point & dorfap)
{
  assert(a != b);
  size_t face[3], ab[2] = {a, b};
  common_tet_face(&original_tet(0, a), &original_tet(0, b), face);
  vector<pair<size_t, size_t> > merge_node; // merge second in b to first in a;
  merge_node.reserve(4);
  vector<size_t> same_vertex_in_cut_tet;
  for(int i = 0; i < 3; ++i) { // unify the node idx in b to a
      int in_ab[2];
      for(int j = 0; j < 2; ++j) {
          for(in_ab[j] = 0; in_ab[j] < 4; ++in_ab[j]) {
              if(original_tet(in_ab[j], ab[j]) == face[i]) break;
            }
        }
      assert(in_ab[0] < 4 && in_ab[1] < 4);
      const size_t unify_to = output_cut_tet(in_ab[0], a);
      const size_t unify_from = output_cut_tet(in_ab[1], b);
      if(unify_from == unify_to){
          same_vertex_in_cut_tet.push_back(unify_from);
          continue;
        }
      merge_node.push_back(make_pair(unify_to, unify_from));
    }
  if(merge_node.empty()) // already glued
    return 1;

  if(is_span_tree_edge && merge_node.size() != 3) // error spanning_tree edge
    cerr << "# error: merge_node.size() != 3 in spanning_tree." << endl;
  if(!is_span_tree_edge && merge_node.size() != 1) // cause loop
    return 2;

  // if is not span tree edge, need to avoid glue singulairty edge
  if(!is_span_tree_edge){
      // check whether the singularity edge is glued inside
      assert(merge_node.size() == 1);
      for(set<pair<size_t,size_t> >::const_iterator spcit
          = singularity_edges.begin(); spcit != singularity_edges.end(); ++spcit){
          if(is_singularity_edge_glue_inside_cut_tet(
               merge_node, original_tet, output_cut_tet, ortae,*spcit))
            return 3;
        }
      // check  whether the merging will result in non-maniofld outside surface
#if 1 // grep the local area field , ana detect whether its manifold
      {
        typedef dyn_one_ring_face_at_point::point2face_type::const_iterator pit;
        typedef dyn_one_ring_face_at_point::linked_face_ptr_types  lfpt;

        vector<size_t> needs_to_handle_points;
        needs_to_handle_points << same_vertex_in_cut_tet.front()
                               << same_vertex_in_cut_tet.back()
                               << merge_node.front().first
                               << merge_node.front().second ;

        set<vector<size_t> > local_faces_set; // here each face is already sorted
        for(size_t t = 0; t < needs_to_handle_points.size(); ++t){
            pit it = dorfap.p2f_.find(needs_to_handle_points[t]);
            if(it == dorfap.p2f_.end()) {
                cerr << "# [error] can not find point " << needs_to_handle_points[t]
                        << " in dyn_one_ring_face_at_point." << endl;
                return __LINE__;
              }
            const lfpt &current_point_linking_faces_ptr = it->second;
            for(lfpt::const_iterator lit = current_point_linking_faces_ptr.begin();
                lit != current_point_linking_faces_ptr.end(); ++lit){
                const dyn_one_ring_face_at_point::face_type &face_ = *(*lit);
                local_faces_set.insert(face_);
              }
          }

        vector<size_t> local_faces_vec;
        local_faces_vec.reserve(local_faces_set.size() * 3); // for triangle face
        for(set<vector<size_t> >::const_iterator svcit = local_faces_set.begin();
            svcit != local_faces_set.end(); ++svcit){
            const vector<size_t> & one_face = *svcit;
            local_faces_vec.insert(local_faces_vec.end(),
                                   one_face.begin(), one_face.end());
          }

        itr_matrix<size_t*> local_faces_mat(3, local_faces_vec.size()/3,
                                            &local_faces_vec[0]);

        dyn_one_ring_face_at_point local_area_dorfap;
        local_area_dorfap.add_all_faces(local_faces_mat);
        pair<size_t,size_t> crease(same_vertex_in_cut_tet.front(),
                                   same_vertex_in_cut_tet.back());
#if 0// debug
        {
          vector<size_t> face0_test, face1_test;
          face0_test << same_vertex_in_cut_tet.front() << same_vertex_in_cut_tet.back()
                     << merge_node.front().first;
          face1_test << same_vertex_in_cut_tet.front() << same_vertex_in_cut_tet.back()
                     << merge_node.front().second;
          sort(face0_test.begin(), face0_test.end());
          sort(face1_test.begin(), face1_test.end());

          if(!local_area_dorfap.is_face_inside(face0_test)){
              cerr << "# [error] strange can not find face before glue ";
              copy(face0_test.begin(), face0_test.end(),
                   ostream_iterator<size_t>(cerr, ","));
              cerr << endl;
            }
          if(!local_area_dorfap.is_face_inside(face1_test)){
              cerr << "# [error] strange can not find face before glue ";
              copy(face1_test.begin(), face1_test.end(),
                   ostream_iterator<size_t>(cerr, ","));
              cerr << endl;
            }
        }
#endif
#if 0  // debug
        {
          matrixst local_face_before_glue;
          local_area_dorfap.convert_face_to_matrix(local_face_before_glue);
          if(!is_surface_mesh_manifold(local_face_before_glue)){
              cerr << "# [error] !!!!!! non-manifold before glueing" << endl;
              return 4;
            }
        }
#endif
        local_area_dorfap.glue_two_faces_at_shared_edge(
              crease, merge_node.front().second, merge_node.front().first);
        matrixst local_face_after_glue;
        local_area_dorfap.convert_face_to_matrix(local_face_after_glue);
        if(is_surface_mesh_manifold(local_face_after_glue))
          dorfap.glue_two_faces_at_shared_edge(
                crease, merge_node.front().second, merge_node.front().first);
        else{
#if 1 // debug
            cerr << "# [error] meet non-manifold while gluing tet "
                 << a << "," << b << endl;
#endif
            return 2;
          }
      }
#endif

#if 0// too slow
      {
        matrixst temp_cut_tet = output_cut_tet;
        const pair<size_t, size_t> & one_merge_node = merge_node.front();
        for(size_t t = 0; t < temp_cut_tet.size(); ++t){
            if(temp_cut_tet[t] == one_merge_node.second)
              temp_cut_tet[t] = one_merge_node.first;
          }
        matrixst outside_face_temp;
        unique_ptr<jtf::mesh::face2tet_adjacent> fa_test(jtf::mesh::face2tet_adjacent::create(temp_cut_tet));
        get_outside_face(*fa_test, outside_face_temp);
        unique_ptr<jtf::mesh::edge2cell_adjacent> ea(jtf::mesh::edge2cell_adjacent::create(outside_face_temp));
        if(!ea.get()) {
            cerr << "# [error] non-manifold happes while gluing tets " << a << "," << b << endl;
            return 2;
          }
      }
#endif
    }

  for(size_t ai = 0; ai < merge_node.size(); ++ai) {
      set<matrixst::value_type *> &node_set = acc_table[merge_node[ai].second];
      assert(!node_set.empty());
      for(set<matrixst::value_type *>::const_iterator iter = node_set.begin();
          iter != node_set.end(); ++iter) {
          assert(*(*iter) == merge_node[ai].second);
          *(*iter) = merge_node[ai].first;
          acc_table[merge_node[ai].first].insert(*iter);
        }
      node_set.clear();
    }

  ec.node_num -= merge_node.size();
  ec.edge_num -= 3;
  ec.face_num -= 2;

  return 0;
}

int  glue_tet_check_inner_singularity(
    const matrixst &input_tet,
    matrixst &output_tet,
    size_t a, size_t b, Euler_characteristic &ec,
    const jtf::mesh::edge2cell_adjacent ea,
    vector<set<matrixst::value_type *> > &acc_table,
    bool is_span_tree_edge,
    const deque<pair<size_t,size_t> > & chain,
    const jtf::mesh::one_ring_tet_at_edge & ortae)
{
  assert(a != b);
  vector<size_t> face(3);
  size_t ab[2] = {a, b};
  common_tet_face(&input_tet(0, a), &input_tet(0, b), &face[0]);
  vector<pair<size_t, size_t> > merge_node; // merge second in b to first in a;
  merge_node.reserve(4);
  //  vector<size_t> merge_node_idx_in_face;
  //  merge_node_idx_in_face.reserve(3);
  for(int i = 0; i < 3; ++i) { // unify the node idx in b to a
      int in_ab[2];
      for(int j = 0; j < 2; ++j) {
          for(in_ab[j] = 0; in_ab[j] < 4; ++in_ab[j]) {
              if(input_tet(in_ab[j], ab[j]) == face[i]) break;
            }
        }
      assert(in_ab[0] < 4 && in_ab[1] < 4);
      const size_t unify_to = output_tet(in_ab[0], a);
      const size_t unify_from = output_tet(in_ab[1], b);
      if(unify_from == unify_to)
        continue;
      merge_node.push_back(make_pair(unify_to, unify_from));
      // merge_node_idx_in_face.push_back(i);
    }
  if(merge_node.empty()) // already glued
    return 1;

  if(is_span_tree_edge && merge_node.size() != 3) // error spanning_tree edge
    cerr << "# error: merge_node.size() != 3 in spanning_tree." << endl;
  //    if(!is_span_tree_edge && merge_node.size() != 1) // cause loop
  //      return 2;

#if 1
  {// detech whether singularity edge will be glued inside
    // when not span_tree edge,

    if(!is_span_tree_edge ){
        // this means the two tets will merge on one node, while our gluing strategy
        // will keep the topology of original mode, so it's impossible to merge just one node
        // which will result in non-manifold. Then this merge operation is based on that
        // one edge has already been merged. After this glue operation, that edge will be
        // glued inside. So checking whether this edge is singularity becomes important
        //      for(size_t t = 0; t < 3; ++t){
        //        //const size_t & need_to_merge_node_idx = merge_node_idx_in_face.front();
        //        const size_t &edge_0 = face[t];
        //        const size_t &edge_1 = face[(t+1)%3];

        //        deque<pair<size_t,size_t> >::const_iterator dit =
        //            find(chain.begin(),chain.end(),make_pair(edge_0,edge_1));
        //        if(dit != chain.end()) return 3; // this singularity edge will be glued inside
        //        dit = find(chain.begin(),chain.end(),make_pair(edge_1,edge_0));
        //        if(dit != chain.end()) return 3;
        //      }
        //if(merge_node.size() == 1){

        for(size_t t = 0; t < chain.size(); ++t){
            if(is_singularity_edge_glue_inside_cut_tet(merge_node,
                                                       input_tet,
                                                       output_tet,
                                                       ortae,chain[t]))
              return 3;
          }
        static set<pair<size_t,size_t> > reach_surface_edge;
        size_t rse = reach_surface_edge.size();
        if(is_contain_singularity_edge_and_surface_edge(face,chain,ea,&reach_surface_edge)){

            if(reach_surface_edge.size() > rse)
              return 4;
          }
        // }
      }
  }
#endif

  for(size_t ai = 0; ai < merge_node.size(); ++ai) {
      set<matrixst::value_type *> &node_set = acc_table[merge_node[ai].second];
      assert(!node_set.empty());
      for(set<matrixst::value_type *>::const_iterator iter = node_set.begin();
          iter != node_set.end(); ++iter) {
          assert(*(*iter) == merge_node[ai].second);
          *(*iter) = merge_node[ai].first;
          acc_table[merge_node[ai].first].insert(*iter);
        }
      node_set.clear();
    }

  ec.node_num -= merge_node.size();
  ec.edge_num -= 3;
  ec.face_num -= 2;
  return 0;
}

void shrink_boundary(const matrixst &input_tet,
                     const fix_graph &fg,
                     matrix<matrixd > &frame_in_tet,
                     const vector<size_t> &spanning_tree,
                     matrixst &output_tet)
{
  output_tet.resize(4, input_tet.size(2));
  output_tet(colon()) = colon(0, output_tet.size()-1);
  vector<set<matrixst::value_type *> > acc_table(output_tet.size());
  for(size_t i = 0; i < acc_table.size(); ++i)
    acc_table[i].insert(&output_tet[i]);

  Euler_characteristic ec;
  ec.node_num = output_tet.size(2)*4,
      ec.edge_num = output_tet.size(2)*6,
      ec.face_num = output_tet.size(2)*4;

  cerr << "# Euler_characteristic before glue: " << ec() << " "
       << ec.node_num << " " << ec.edge_num << " " << ec.face_num << endl;
  for(size_t ei = 0; ei < spanning_tree.size(); ++ei) {
      std::pair<node_idx, node_idx> edge = fg.get_edge(spanning_tree[ei]);
      if(glue_tet(input_tet, output_tet, edge.first, edge.second, ec, acc_table, true) > 1) {
          cerr << "glue according spanning_tree error." << endl;
        }
    }
  cerr << "# Euler_characteristic after glue: " << ec() << " "
       << ec.node_num << " " << ec.edge_num << " " << ec.face_num << endl;
#if 1 // debug
  vector<pair<double, pair<size_t, size_t> > > missing_edge;
  vector<pair<size_t, size_t> > rot_jump;
  matrixd shuffle_rot(3, 3);
  missing_edge.reserve(fg.edge_num());
  for(size_t ai = 0; ai < fg.node_num(); ++ai) {
      for(size_t ei = fg.ptr_[ai]; ei < fg.ptr_[ai+1]; ++ei) {
          const size_t bi = fg.idx_[ei];
          if(ai == bi) continue;

          get_best_alignment(&frame_in_tet[ai][0], &frame_in_tet[bi][0],
              &shuffle_rot[0]);
          if(norm(shuffle_rot-eye<double>(3)) > 1e-6) { // need rotation transition
              rot_jump.push_back(make_pair(ai, bi));
              continue;
            }
          // if(ai < 100)
          // 	continue;
          const double err = norm(frame_in_tet[ai]-frame_in_tet[bi]);
          missing_edge.push_back(make_pair(err, make_pair(ai, bi)));
        }
    }
  cerr << "# possible additional edges: " << missing_edge.size() << endl;
  sort(missing_edge.begin(), missing_edge.end());

  cout << "# max edge error: " << max_element(missing_edge.begin(),
                                              missing_edge.end())->first << endl;

  while(!missing_edge.empty()) {
      cout << "# missing_edge: " << missing_edge.size() << endl;
      vector<pair<double, pair<size_t, size_t> > > tmp;
      for(size_t i = 0; i < missing_edge.size(); ++i) {
          const int rtn =
              glue_tet(input_tet, output_tet,
                       missing_edge[i].second.first,
                       missing_edge[i].second.second,
                       ec, acc_table, false);
          assert(ec() == 2);
          if(rtn < 2) // sucess or has been added
            continue;
          if(rtn == 2)
            tmp.push_back(missing_edge[i]);
        }
      if(tmp.size() == missing_edge.size()) // no more edge can be added
        break;
      swap(tmp, missing_edge);
    }
#endif
  cerr << "# Euler_characteristic after glue: " << ec() << " "
       << ec.node_num << " " << ec.edge_num << " " << ec.face_num << endl;
  cerr << "# rotation jump: " << rot_jump.size() << endl;
  cerr << "# integration jump: " << missing_edge.size() << endl;

#if 1
  {
    size_t count = 0;
    for(size_t ai = 0; ai < acc_table.size(); ++ai) {
        if(acc_table[ai].empty()) continue;
        ++count;
      }
    if(count != ec.node_num) {
        cerr << "error in glue." << endl;
      }
  }
#endif
}

//int modify_face_jump_type_at_given_tets_edge(const size_t tet_idx_a,
//                                             const size_t tet_idx_b,
//                                             map<pair<size_t,size_t>,size_t> &inner_face_jump_type,
//                                             const vector<size_t> &edge0_tet);

bool is_singularity_edge(const size_t e0,
                         const size_t e1,
                         const vector<pair<size_t,size_t> > &singularity_edges)
{
  if(find(singularity_edges.begin(),singularity_edges.end(),make_pair(e0,e1)) == singularity_edges.end()
     && find(singularity_edges.begin(),singularity_edges.end(),make_pair(e1,e0)) == singularity_edges.end())
    return false;
  return true;
}

bool is_edge_glue_inside(const size_t e0,
                         const size_t e1,
                         const size_t other_v,
                         const jtf::mesh::one_ring_tet_at_edge &ortae,
                         const matrixst &input_tet,
                         const matrixst &output_tet)
{
  typedef jtf::mesh::one_ring_tet_at_edge::e2tet_type::const_iterator ci;
  ci i = ortae.e2t_.find(make_pair(e0,e1));
  if(i == ortae.e2t_.end()) i = ortae.e2t_.find(make_pair(e1,e0));
  if(i == ortae.e2t_.end()) {
      cerr << "# strange, can not find edge:" << e0 << "," << e1 << endl;
      return false;
    }else{
      map<size_t,vector<size_t> > tet_v_idx; // stroe the vertex idx of tet_i which is different with e0 and e1
      const vector<size_t> &tet_loop = i->second;

      if(tet_loop.front() != tet_loop.back()) return false; // sharp edge

      size_t common_face_num = tet_loop.size() - 1;

      size_t t = 0;
      if(tet_loop[0] == -1) {
          t = 1;
        }
      for(; t < tet_loop.size() - 1; ++t) // the last tet is the same with the first
        {
          for(size_t j = 0; j < 4; ++j) {
              if(input_tet(j,tet_loop[t]) != e0
                 && input_tet(j,tet_loop[t]) != e1
                 && input_tet(j,tet_loop[t]) != other_v) // other_v may be glue to one vertex
                tet_v_idx[tet_loop[t]].push_back(j);
            }
        }
      set<size_t> output_vertex;
      //vector<size_t> one_ring_vertex;
      typedef map<size_t,vector<size_t> >::const_iterator mci;
      for(mci ci = tet_v_idx.begin(); ci != tet_v_idx.end(); ++ci){
          for(size_t t = 0; t < ci->second.size(); ++t)
            {
              output_vertex.insert(output_tet(ci->second[t],ci->first));
              //      one_ring_vertex.push_back(output_tet(ci->second[t],ci->first));
            }
        }
      if(output_vertex.size() + 1 == common_face_num )
        return true;
      else
        return false;
    }
}

bool is_contain_singularity_edge_and_surface_edge(
    const vector<size_t> &face,
    const deque<pair<size_t,size_t> > &chain,
    const jtf::mesh::edge2cell_adjacent &ea,
    set<pair<size_t,size_t> > *singularity_edge)
{
  bool is_singularity = false;
  bool is_on_surface = false;
  for(size_t t = 0; t < face.size(); ++t){
      pair<size_t,size_t> edge(face[t],face[(t+1)%face.size()]);
      if(!is_singularity){
          if(find(chain.begin(),chain.end(),edge) != chain.end()) is_singularity = true;
          else if(find(chain.begin(),chain.end(),make_pair(edge.second,edge.first)) != chain.end())
            is_singularity = true;
          if(singularity_edge != 0 && is_singularity){
              if(edge.first > edge.second)
                singularity_edge->insert(make_pair(edge.second,edge.first));
              else
                singularity_edge->insert(edge);
            }
        }
      if(!is_on_surface){
          if(ea.get_edge_idx(edge.first,edge.second) != -1
             || ea.get_edge_idx(edge.second,edge.first) != -1){
              is_on_surface = true;
            }
        }
      if(is_singularity && is_on_surface) return true;
    }
  return false;
}

bool is_singularity_edge_glue_inside_cut_tet(
    const vector<pair<size_t,size_t> > & merge_node,
    const matrixst & tet,
    const matrixst & cut_tet,
    const jtf::mesh::one_ring_tet_at_edge & ortae,
    const pair<size_t,size_t> & singularity_edge)
{
  // assert(merge_node.first != merge_node.second);
  jtf::mesh::one_ring_tet_at_edge::e2tet_type::const_iterator it
      = ortae.e2t_.find(singularity_edge);
  if(it == ortae.e2t_.end())
    it = ortae.e2t_.find(make_pair(singularity_edge.second,
                                   singularity_edge.first));
  if(it == ortae.e2t_.end()) {
      cerr << "# [error] can not find singularity edge in such ortae." << endl;
      return true;
    }

  const vector<size_t> & loop = it->second;
  assert(loop.front() == loop.back());
  assert(find(loop.begin(), loop.end(), -1) == loop.end());

  set<size_t> cut_tet_out_vertex;

  size_t check_num = 0;
  for(size_t t = 0; t < loop.size()-1; ++t){
      check_num = 0;
      assert(find(tet(colon(), loop[t]).begin(), tet(colon(), loop[t]).end(),
                  singularity_edge.first) != tet(colon(), loop[t]).end());
      assert(find(tet(colon(), loop[t]).begin(), tet(colon(), loop[t]).end(),
                  singularity_edge.second) != tet(colon(), loop[t]).end());
      for(size_t i = 0; i < 4; ++i){
          if(tet(i,loop[t]) != singularity_edge.first &&
             tet(i,loop[t]) != singularity_edge.second) {
              ++check_num;
              cut_tet_out_vertex.insert(cut_tet(i,loop[t]));
            }
        }
      if(check_num != 2){
          cerr << "# [error] out vertex num is not 2: " << check_num << endl;
        }
    }
  size_t out_vertex_size = cut_tet_out_vertex.size();
  // outside vertex num > around tet +1, because loop.size() = around tet + 1
  if(out_vertex_size > loop.size() ) return false;

  for(size_t t = 0; t < merge_node.size(); ++t){
      const pair<size_t,size_t> & merge_pair = merge_node[t];
      assert(merge_pair.first != merge_pair.second);
      if(find(cut_tet_out_vertex.begin(),
              cut_tet_out_vertex.end(),
              merge_pair.first) != cut_tet_out_vertex.end() &&
         find(cut_tet_out_vertex.begin(),
              cut_tet_out_vertex.end(),
              merge_pair.second) != cut_tet_out_vertex.end()){
          --out_vertex_size;
        }
    }
  assert(out_vertex_size >= loop.size() - 1);
  if(out_vertex_size == loop.size() - 1) return true;
  return false;
}

bool is_singularity_edge_glue_inside(const size_t *face,
                                     const matrixst &input_tet,
                                     const matrixst &cut_tet,
                                     const jtf::mesh::one_ring_tet_at_edge &ortae,
                                     const vector<pair<size_t,size_t> > &singularity_edges)
{
  typedef jtf::mesh::one_ring_tet_at_edge::e2tet_type::const_iterator ci;
  for(size_t t = 0; t < singularity_edges.size(); ++t){
      ci it = ortae.e2t_.find(singularity_edges[t]);
      if(it == ortae.e2t_.end()) it = ortae.e2t_.find(make_pair(singularity_edges[t].second,singularity_edges[t].first));
      if(it == ortae.e2t_.end()) {
          cerr << "# strange can not find this singularity edge." << endl;
          return false;
        }

      const pair<size_t,size_t> & edge = singularity_edges[t];
      const vector<size_t>& tet_loop = it->second;
      if(find(face,face+3, edge.first) != face+3
         && find(face,face+3, edge.second) != face+3) // this edge is singularity edge
        {
          size_t other_vertex;
          for(size_t i = 0; i < 3; ++i) {
              if(face[i] != edge.first && face[i] != edge.second){
                  other_vertex = face[i];
                  break;
                }
            }
          if(is_edge_glue_inside(edge.first,edge.second,other_vertex,ortae,input_tet,cut_tet))
            return true;
        }else // this singularity edge is not on the face, but may have one vertex on face.
        //need to check whether this singularity edge will be glue inside if three vertex of this face is merged
        {
          size_t common_face_num = tet_loop.size() - 1; // record the common face num
          map<size_t,list<size_t> > tet_v_idx;
          set<size_t> one_ring_vertex_around_edge;
          for(size_t i = 0; i < tet_loop.size() - 1; ++i) // the last tet is the same with the first
            {
              for(size_t j = 0; j < 4; ++j) {
                  if(input_tet(j,tet_loop[i]) != edge.first
                     && input_tet(j,tet_loop[i]) != edge.second)
                    {
                      one_ring_vertex_around_edge.insert(input_tet(j,tet_loop[i]));
                    }

                  if(input_tet(j,tet_loop[i]) != edge.first
                     && input_tet(j,tet_loop[i]) != edge.second
                     /*                            && input_tet(j,tet_loop[i]) != face[0]
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                        && input_tet(j,tet_loop[i]) != face[1]
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                        && input_tet(j,tet_loop[i]) != face[2]*/){
                      tet_v_idx[tet_loop[i]].push_back(j);
                    }
                }
            }

          size_t find_face_vertex_num = 0;
          vector<size_t> find_face_vertex;
          for(size_t i = 0; i < 3; ++i){
              if(find(one_ring_vertex_around_edge.begin(),one_ring_vertex_around_edge.end(),face[i])
                 != one_ring_vertex_around_edge.end())
                // ++find_face_vertex_num;
                find_face_vertex.push_back(face[i]);
            }

          //            if(find_face_vertex.size() == 0) continue; // the face glue can not impact this singularity edge.

          typedef map<size_t,list<size_t> >::const_iterator mci;
          set<size_t> merged_vertex;
          if(find_face_vertex.size() == 1){
              for(mci ci = tet_v_idx.begin(); ci != tet_v_idx.end(); ++ci){
                  const list<size_t> &list_ = ci->second;
                  for(list<size_t>::const_iterator lci = list_.begin(); lci != list_.end(); ++lci){
                      if(input_tet(*lci,ci->first) == find_face_vertex[0]){
                          merged_vertex.insert(cut_tet(*lci,ci->first));
                        }
                    }
                }
              if(merged_vertex.size() == 1) continue; //this ver
              if(is_edge_glue_inside(edge.first,edge.second,find_face_vertex[0],ortae,input_tet,cut_tet))
                return true;
              continue;
            }

          //            set<size_t> output_vertex;

          //            //vector<size_t> one_ring_vertex;
          //            typedef map<size_t,list<size_t> >::const_iterator mci;
          //            for(mci ci = tet_v_idx.begin(); ci != tet_v_idx.end(); ++ci){
          //                const list<size_t> &list_ = ci->second;
          //                if(list_.empty()) continue;
          //                for(list<size_t>::const_iterator lci = list_.begin(); lci != list_.end(); ++lci)
          //                    output_vertex.insert(cut_tet(*lci,ci->first));
          //            }

          //            if(output_vertex.size() + find_face_vertex_num == common_face_num ) // this edge is glued inside
          //                return true;
          //            else
          //                return false;
          //        }
        }
    }
  return false; // the face do not impact any singularity edge
}

void shrink_boundary_new(const matrixst &input_tet,
                         const matrixd &node,
                         const fix_graph &fg,
                         const matrix<matrixd > &frame_in_tet,
                         const vector<size_t> &spanning_tree,
                         matrixst &output_tet,
                         const jtf::mesh::one_ring_tet_at_edge &ortae,
                         const vector<pair<size_t,size_t> > &singularity_edges)
{
  output_tet.resize(4, input_tet.size(2));
  output_tet(colon()) = colon(0, output_tet.size()-1);
  vector<set<matrixst::value_type *> > acc_table(output_tet.size());
  for(size_t i = 0; i < acc_table.size(); ++i)
    acc_table[i].insert(&output_tet[i]);

  Euler_characteristic ec;
  ec.node_num = output_tet.size(2)*4,
      ec.edge_num = output_tet.size(2)*6,
      ec.face_num = output_tet.size(2)*4;

  cerr << "# Euler_characteristic before glue: " << ec() << " "
       << ec.node_num << " " << ec.edge_num << " " << ec.face_num << endl;
  for(size_t ei = 0; ei < spanning_tree.size(); ++ei) {
      std::pair<node_idx, node_idx> edge = fg.get_edge(spanning_tree[ei]);
      if(glue_tet(input_tet, output_tet, edge.first, edge.second, ec, acc_table, true) > 1) {
          cerr << "glue according spanning_tree error." << endl;
        }
    }
  cerr << "# Euler_characteristic after glue: " << ec() << " "
       << ec.node_num << " " << ec.edge_num << " " << ec.face_num << endl;
#if 1
  {
    size_t count = 0;
    for(size_t ai = 0; ai < acc_table.size(); ++ai) {
        if(acc_table[ai].empty()) continue;
        ++count;
      }
    if(count != ec.node_num) {
        cerr << "error in glue." << endl;
      }
  }
#endif
}

void orient_type(
    const fix_graph &graph,
    const size_t & tet_num,
    const jtf::mesh::face2tet_adjacent & fa,
    boost::unordered_map<std::pair<size_t,size_t>,size_t> & inner_face_jump_type,
    boost::unordered_map<size_t,size_t> & surface_type,
    const size_t surface_type_flag,
    zjucad::matrix::matrix<matrixd> * frame_ptr)
{
  // cerr << "# [degenerou] not used" << endl;
  matrixst  frame_type_to_seed(tet_num,1);
  stack<size_t> s;
  s.push(0);
  vector<size_t> father_tet_idx(tet_num,-1);
  father_tet_idx[0]=0;
  frame_type_to_seed[0] = TRIVIAL_TYPE;
  vector<bool> visted_flag(tet_num, false);
  visted_flag[0] = true;

  while(!s.empty()) {
      size_t current = s.top();
      visted_flag[current] = true;
      s.pop();
      for(size_t ni = graph.ptr_[current]; ni < graph.ptr_[current+1]; ++ni) {
          size_t next = graph.idx_[ni];

          if(visted_flag[next]) continue;
          father_tet_idx[next] = current;
          // F_next * \Pi_{next, current} = F_current
          boost::unordered_map<pair<size_t,size_t>,size_t>::const_iterator bumpcit
              = inner_face_jump_type.find(make_pair(next,current));
          if(bumpcit == inner_face_jump_type.end()){
              frame_type_to_seed[next] = frame_type_to_seed[current];
            }
          else{
              matrixd rot_next = type_transition2(bumpcit->second) *
                  type_transition2(frame_type_to_seed[current]);
              frame_type_to_seed[next] = type_transition1(rot_next);
            }

          if(frame_ptr){
              (*frame_ptr)[next] = temp((*frame_ptr)[next] *
                                        type_transition2(frame_type_to_seed[next]));
            }

          s.push(next);
        }
    }

#if 1
  // replace all jump_type
  //inner_face_jump_type.clear();
  bool need_to_orient_surface = true;
  if(surface_type.empty())
    need_to_orient_surface = false;
  if(need_to_orient_surface){
      cerr << "# [WARNING!!!] only surface type from frame field can be global aligned." << endl;
    }

  const matrixd I = eye<double>(3);
  vector<pair<double,int> > res;

  boost::unordered_map<pair<size_t,size_t>,size_t> inner_face_jump_type_temp
      = inner_face_jump_type;
  inner_face_jump_type.clear();
  matrixd frame_type_rot = eye<double>(3);
  for(size_t fi = 0; fi < fa.face2tet_.size(); ++fi){
      const pair<size_t,size_t> & tet_pair = fa.face2tet_[fi];
      if(fa.is_outside_face(tet_pair)){
          if( need_to_orient_surface && surface_type_flag == 1){
              const size_t tet_idx =
                  (tet_pair.first == -1? tet_pair.second: tet_pair.first);
              frame_type_rot = type_transition2(frame_type_to_seed[tet_idx]);
              boost::unordered_map<size_t, size_t>::iterator bumit =
                  surface_type.find(fi);
              if(bumit == surface_type.end()){
                  surface_type.insert(
                        make_pair(fi, type_transition1(trans(frame_type_rot))));
                  //make_pair(fi, type_transition1((frame_type_rot))));
                }else{
                  size_t & surface_type_idx  = bumit->second;
                  surface_type_idx =
                      type_transition1(trans(frame_type_rot) *
                                       type_transition2(surface_type_idx));
                  //                      type_transition1(frame_type_rot*
                  //                                       type_transition2(surface_type_idx));
                }
            }else if( need_to_orient_surface && surface_type_flag == 0){
              cerr << "# [error] can not orient restricted surface type." << endl;
              //        const vector<size_t> & face_vec = fa.faces_[fi];
              ////        if(find(face_vec.begin(), face_vec.end(), 46) != face_vec.end() &&
              ////           find(face_vec.begin(), face_vec.end(), 47) != face_vec.end() &&
              ////           find(face_vec.begin(), face_vec.end(), 1189) != face_vec.end())
              ////          cerr << endl;
              //        const size_t tet_idx =
              //            (tet_pair.first == -1? tet_pair.second: tet_pair.first);
              //        frame_type_rot = type_transition2(frame_type_to_seed[tet_idx]);
              //        boost::unordered_map<size_t, size_t>::iterator bumit =
              //            surface_type.find(fi);
              //        if(bumit == surface_type.end()){
              //          cerr << "# [error] strange can not find surface type of face " << fi << endl;
              //          return __LINE__;
              //        }

              //        assert(bumit->second < 3);
              //        res.clear();
              //        for(size_t i = 0; i < 3; ++i){
              //          res.push_back(make_pair(dot(frame_type_rot(colon(), bumit->second),
              //                                      I(colon(),i)), i * 2 + 0));
              //          res.push_back(make_pair(dot(frame_type_rot(colon(), bumit->second),
              //                                      -1*I(colon(),i)), i * 2 + 1));
              //        }
              //        sort(res.begin(), res.end());
              //        bumit->second = res.back().second/2;
            }
        }else{ // inner face
          boost::unordered_map<pair<size_t,size_t>,size_t>::iterator bumit =
              inner_face_jump_type_temp.find(tet_pair);
          matrixd rot = eye<double>(3);
          if(bumit != inner_face_jump_type_temp.end()){
              rot = type_transition2(bumit->second);
            }

          const size_t rot_type =
              type_transition1(
                trans(type_transition2(frame_type_to_seed[tet_pair.first]))
              * rot * type_transition2(frame_type_to_seed[tet_pair.second]));

          //      if(father_tet_idx[tet_pair.first] == tet_pair.second
          //         || father_tet_idx[tet_pair.second] == tet_pair.first)
          //      {
          //        if(rot_type != TRIVIAL_TYPE)
          //          cerr << endl;
          //      }
          if(is_trivial_type(rot_type)) continue;
          inner_face_jump_type[tet_pair] = rot_type;
          inner_face_jump_type[make_pair(tet_pair.second, tet_pair.first)] =
              type_transition1( trans(type_transition2(rot_type)));
        }
    }
#endif
}

void orient_frame(matrix<matrixd > &frame_in_tet,
                  const fix_graph &graph,
                  matrixd &orient_err)
{
  assert(frame_in_tet[0].size() != 0);
  stack<size_t> s;
  s.push(0);
  frame_in_tet[0].resize(3, 3);
  matrixd R(3, 3);
  vector<bool> visited_tets(frame_in_tet.size(), false);
  visited_tets[0] = true;
  while(!s.empty()) {
      size_t current = s.top();
      s.pop();
      for(size_t ni = graph.ptr_[current]; ni < graph.ptr_[current+1]; ++ni) {
          size_t next = graph.idx_[ni];
          if(visited_tets[next]) continue;
          get_best_alignment(&frame_in_tet[next][0], &frame_in_tet[current][0], &R[0]);
          frame_in_tet[next] = temp(frame_in_tet[next])*R;

          visited_tets[next] = true;
#if 0
          if(norm(frame_in_tet[next]*trans(frame_in_tet[next]) - eye<double>(3)) > 1e-8) {
              cerr << "# bad frame: " << frame_in_tet[next] << endl;
            }
#endif

          orient_err[ni] = norm(frame_in_tet[next]-frame_in_tet[current]);

#if 0 // check
          if(norm(F-frame_in_tet[current]) < orient_err[ni])
            cerr << "# error in orient_frame." << endl;
#endif
          s.push(next);
        }
    }
}


void orient_frame(const matrixd &zyz_frame_in_tet,
                  matrix<matrixd > &frame_in_tet,
                  const fix_graph &graph,
                  matrixd &orient_err)
{
  if(frame_in_tet.size() != zyz_frame_in_tet.size(2))
    frame_in_tet.resize(zyz_frame_in_tet.size(2));

  const size_t seed = 0;
  stack<size_t> s;
  s.push(seed);
  frame_in_tet[seed].resize(3, 3);
  zyz_angle_2_rotation_matrix1(&zyz_frame_in_tet(seed, 0), &frame_in_tet[seed][0]);
  matrixd F(3, 3), R(3, 3);
  while(!s.empty()) {
      size_t current = s.top();
      s.pop();
      for(size_t ni = graph.ptr_[current]; ni < graph.ptr_[current+1]; ++ni) {
          size_t next = graph.idx_[ni];
          if(frame_in_tet[next].size() == 9) continue;
          zyz_angle_2_rotation_matrix1(&zyz_frame_in_tet(0, next), &F[0]);
          get_best_alignment(&F[0], &frame_in_tet[current][0], &R[0]);
          frame_in_tet[next] = F*R;
#if 0
          if(norm(frame_in_tet[next]*trans(frame_in_tet[next]) - eye<double>(3)) > 1e-8) {
              cerr << "# bad frame: " << frame_in_tet[next] << endl;
            }
#endif

          orient_err[ni] = norm(frame_in_tet[next]-frame_in_tet[current]);

#if 0 // check
          if(norm(F-frame_in_tet[current]) < orient_err[ni])
            cerr << "# error in orient_frame." << endl;
#endif			
          s.push(next);
        }
    }
}

bool is_edge_singularity_in_original_tet(
    const pair<size_t,size_t> &edge_in_cut_tet,
    const matrixst & cut_tet,
    const matrixst & cut_tet2tet,
    const deque<pair<size_t,size_t> > &chain)
{
  pair<size_t,size_t> edge_in_ori_tet(cut_tet2tet[edge_in_cut_tet.first],
      cut_tet2tet[edge_in_cut_tet.second]);

  typedef deque<pair<size_t,size_t> >::const_iterator dpcit;
  dpcit it = find(chain.begin(),chain.end(),edge_in_ori_tet);
  if(it != chain.end()) return true;
  it = find(chain.begin(),chain.end(),
            make_pair(edge_in_ori_tet.second,edge_in_ori_tet.first));
  if(it != chain.end()) return true;
  return false;
}

static bool is_original_singularity(
    const pair<size_t,size_t> &each_edge_in_cut_tet,
    const matrixst & cut_tet2tet,
    const set<pair<size_t,size_t> > &singularity_edges)
{
  pair<size_t,size_t>  edge_in_original(
        cut_tet2tet[each_edge_in_cut_tet.first],
      cut_tet2tet[each_edge_in_cut_tet.second]);
  set<pair<size_t,size_t> >::const_iterator
      spcit = singularity_edges.find(edge_in_original);
  if(spcit == singularity_edges.end())
    spcit = singularity_edges.find(make_pair(edge_in_original.second,
                                             edge_in_original.first));
  if(spcit == singularity_edges.end())
    return false;
  return true;
}

int analysis_face_type(
    const matrixst & input_tet,
    const matrixd & input_node,
    const matrixst & cut_tet,
    const jtf::mesh::face2tet_adjacent &input_fa,
    const boost::unordered_map<pair<size_t,size_t>,size_t> & inner_face_jump_type,
    const set<pair<size_t,size_t> > & singularity_edges,
    matrixst & face_pair,
    matrixst & face_type,
    const matrixd * param)
{
  typedef boost::unordered_map<pair<size_t,size_t>,size_t>::const_iterator mpscit;
  unique_ptr<jtf::mesh::face2tet_adjacent> output_fa(jtf::mesh::face2tet_adjacent::create(cut_tet));
  const jtf::mesh::face2tet_adjacent *fa[2] = { &input_fa, output_fa.get()};
  // compact node index, setup the mapping of cut_tet2tet
  matrixst cut_tet2tet(max(cut_tet)+1);
  cut_tet2tet(cut_tet) = input_tet(colon());
  matrixst outside_face_original;
  get_outside_face(*fa[0], outside_face_original);

  matrixst outside_face_idx_in_output;
  get_outside_face_idx(*fa[1], outside_face_idx_in_output);

  matrixst outside_face_in_cut_tet;
  get_outside_face(*fa[1], outside_face_in_cut_tet);

  unique_ptr<jtf::mesh::edge2cell_adjacent> ea_original_surface(jtf::mesh::edge2cell_adjacent::create(outside_face_original));
  if(!ea_original_surface.get()) {
      cerr << "# [error] build edge2cell_adjacent fail." << endl;
      return __LINE__;
    }

  unique_ptr<jtf::mesh::edge2cell_adjacent_general>
      ea(jtf::mesh::edge2cell_adjacent_general::create(outside_face_in_cut_tet));
  if(!ea.get()) {
      cerr << "# [error] build edge2cell_adjacent fail." << endl;
      return __LINE__;
    }

  //! map key of original tet face to the index in outside_face_idx_in_output
  typedef boost::unordered_map<vector<size_t>, pair<size_t, size_t> > face_pair_type;
  face_pair_type outside_face;
  {
    matrixst face(3);
    for(size_t fi = 0; fi < outside_face_idx_in_output.size(); ++fi) {
        const size_t face_id = outside_face_idx_in_output[fi];
        const pair<size_t, size_t> &nb_tet_id = fa[1]->face2tet_[face_id];
        VERIFY(fa[1]->is_outside_face(nb_tet_id), "# error in nb_tet_id.");
        const size_t tet_id = (nb_tet_id.first == -1)?
              nb_tet_id.second:nb_tet_id.first;

        vector<size_t> ori_face = to_vec(cut_tet2tet(fa[1]->get_face(face_id)));
        sort(ori_face.begin(), ori_face.end());

        face_pair_type::iterator i = outside_face.find(ori_face);
        if(i == outside_face.end()) {
            outside_face[ori_face] = make_pair(fi, -1);
          }
        else {
            VERIFY(outside_face[ori_face].second == -1, "# error in outside_face.");
            outside_face[ori_face].second = fi;
          }
      }
    // face_pair stores the index to outside_face_idx_in_output
    face_pair.resize(outside_face_idx_in_output.size());
    face_type.resize(outside_face_idx_in_output.size());
    for(size_t fi = 0; fi < outside_face_idx_in_output.size(); ++fi) {
        const size_t face_id = outside_face_idx_in_output[fi];
        vector<size_t> ori_face = to_vec(cut_tet2tet(fa[1]->get_face(face_id)));
        sort(ori_face.begin(), ori_face.end());
        const pair<size_t, size_t> &id = outside_face[ori_face];
        VERIFY(id.first == fi || id.second == fi, "# outside_face map error.");
        //			face_pair[fi] = id.first+id.second-fi;
        face_pair[fi] = (id.first == fi)?id.second:id.first;

#if 1 // add face type
        if(face_pair[fi] == -1)
          face_type[fi] = -1;
        else{
            const pair<size_t,size_t> & tet_pair_from =
                fa[1]->face2tet_[outside_face_idx_in_output[fi]];
            const pair<size_t,size_t> & tet_pair_to =
                fa[1]->face2tet_[outside_face_idx_in_output[face_pair[fi]]];
            assert(tet_pair_from.first == -1 || tet_pair_from.second == -1);
            assert(tet_pair_to.first == -1 || tet_pair_to.second == -1);
            const size_t & tet_from_idx = (tet_pair_from.first == -1?
                                             tet_pair_from.second:
                                             tet_pair_from.first);
            const size_t & tet_to_idx = (tet_pair_to.first == -1?
                                           tet_pair_to.second:
                                           tet_pair_to.first);
            mpscit it = inner_face_jump_type.find(make_pair(tet_from_idx, tet_to_idx));
            if(it == inner_face_jump_type.end()){
                face_type[fi] = 9; // trivial type
              }else
              face_type[fi] = it->second;
          }

#endif

#if 1
        if(face_pair[fi] != -1) {
            matrixst oriA = cut_tet2tet(fa[1]->get_face(
                  outside_face_idx_in_output[fi]));
            matrixst oriB =
                cut_tet2tet(fa[1]->get_face(
                  outside_face_idx_in_output[face_pair[fi]]));
            sort(oriA.begin(), oriA.end());
            sort(oriB.begin(), oriB.end());
            if(norm(oriA-oriB) > 1e-8) {
                cerr << "# error in face_pair." << endl;
              }
          }
#endif
      }
  }

#if 1 // visual
  {
    vector<size_t> cut_faces;
    typedef face_pair_type::const_iterator mvpcit;
    for(mvpcit it = outside_face.begin(); it != outside_face.end(); ++it){
        const pair<size_t,size_t> & each_face_pair = it->second;
        if(each_face_pair.first == -1 || each_face_pair.second == -1) continue;
        const vector<size_t> & face_in_cut_tet = it->first;
        cut_faces.insert(cut_faces.end(), face_in_cut_tet.begin(), face_in_cut_tet.end());
      }
    ofstream ofs("cut_face.vtk");
    tri2vtk(ofs, &input_node[0], input_node.size(2), &cut_faces[0], cut_faces.size() /3);

  }
#endif

  vector< set<size_t> > face_groups; // stores the face groups
#if 1 // group face accronding the rotation jump type
  {
    const size_t face_seeding = 0;

    stack<size_t> one_group_stack;
    //one_group.push(face_seeding);

#if 1 // group faces
    set<size_t> not_visited_face_idx;
    assert(face_pair.size() == outside_face_idx_in_output.size());
    for(size_t i = 0; i < outside_face_idx_in_output.size(); ++i) {
        if(face_pair[i] == -1) continue;
        not_visited_face_idx.insert(i);
      }

    vector<size_t> each_face(3); // store the face with vertex idx in original tet
    //size_t parent_face_type = -1;
    typedef map<pair<size_t,size_t>,size_t>::const_iterator mpscit;
    vector<size_t> common_edge(2);
    while(!not_visited_face_idx.empty()) // some faces are still not visited
      {
        assert(one_group_stack.empty());

        cerr << "# [info] left " << not_visited_face_idx.size()
             << " face not visited" << endl;
        one_group_stack.push(*(not_visited_face_idx.begin()));

        set<size_t> one_grounp_vec;
        // grow to be a group
        while(!one_group_stack.empty()){
            const size_t current_face_idx = one_group_stack.top();
            one_group_stack.pop();

            one_grounp_vec.insert(outside_face_idx_in_output[current_face_idx]);

            set<size_t>::iterator sit =
                find(not_visited_face_idx.begin(), not_visited_face_idx.end(),
                     current_face_idx);
            if(sit == not_visited_face_idx.end())  // this face has been visited
              continue;

            not_visited_face_idx.erase(sit);

            const size_t oppsite_face_idx = face_pair[current_face_idx];
            each_face = to_vec(fa[1]->get_face(
                  outside_face_idx_in_output[current_face_idx]));

            for(size_t t = 0; t < each_face.size(); ++t){
                const pair<size_t,size_t> each_edge(each_face[t],
                                                    each_face[(t+1)%each_face.size()]);

                if(is_original_singularity(each_edge, cut_tet2tet,
                                           singularity_edges)){
                    continue;
                  }
                const size_t edge_idx =
                    ea->get_edge_idx(each_edge.first, each_edge.second);
                if(edge_idx == -1) {
                    cerr << "# [error] can not find such edge <" << each_edge.first
                         << each_edge.second  << ">."<< endl;
                    return __LINE__;
                  }
                const vector<size_t> & adj_faces = ea->edge2cell_[edge_idx];
                if(adj_faces.size() == 1) continue;
                if(adj_faces.size() > 2){
                    cerr << "# [error] edge meets " << adj_faces.size() << " faces"
                         << " it's non-manifold." << endl;
                    cerr << "# [error] edge: <" << cut_tet2tet[each_edge.first] << ","
                         << cut_tet2tet[each_edge.second] << ">." << endl;
                    return __LINE__;
                  } // if adj_faces > 2 it's non-manifold
                assert(adj_faces.size() == 2);
                assert(adj_faces.front()== current_face_idx ||
                       adj_faces.back() == current_face_idx);
                const size_t adjacent_face_idx = adj_faces.front() +
                    adj_faces.back() - current_face_idx;
                if(find(not_visited_face_idx.begin(),  not_visited_face_idx.end(),
                        adjacent_face_idx) == not_visited_face_idx.end())
                  continue;  // this adjacent face has been visited

                const size_t oppo_adjacent_face_idx = face_pair[adjacent_face_idx];
                if(oppo_adjacent_face_idx == -1) // adjacent_face_idx is a  outside surface triangle
                  continue;
#if 0
                {
                  const matrixst adj_face_mat =
                      fa[1]->get_face(outside_face_idx_in_output[adjacent_face_idx]);
                  cerr << "pause" << endl;
                }
#endif
                const matrixst oppo_adj_face_mat =
                    fa[1]->get_face(outside_face_idx_in_output[oppo_adjacent_face_idx]);
                const matrixst oppsite_face_mat =
                    fa[1]->get_face(outside_face_idx_in_output[oppsite_face_idx]);

                if(!find_common_edge(oppo_adj_face_mat, oppsite_face_mat,&common_edge[0]))
                  one_group_stack.push(adjacent_face_idx);
              }
          }
        face_groups.push_back(one_grounp_vec);
      }
#endif

    // visuallize each group
    { // visual
      ofstream ofs("cut_tet_surface_patch.vtk");
      vector<size_t> out_faces;
      vector<size_t> out_faces_type;
      vector<size_t> each_face(3);
      vector<size_t> out_faces_in_param;
      vector<size_t> each_face_in_param(3);
      for(size_t t = 0; t < face_groups.size(); ++t){
          ostringstream os;
          os << t;
          ofstream ofs_group((os.str() + "group.vtk").c_str());
          vector<size_t> each_group;
          const set<size_t> & one_set = face_groups[t];
          for(set<size_t>::const_iterator scit = one_set.begin();
              scit != one_set.end(); ++scit){
              each_face = to_vec(cut_tet2tet(fa[1]->get_face(*scit)));
              each_face_in_param = to_vec(fa[1]->get_face(*scit));

              out_faces.insert(out_faces.end(), each_face.begin(), each_face.end());
              out_faces_in_param.insert(out_faces_in_param.end(),
                                        each_face_in_param.begin(),
                                        each_face_in_param.end());
              out_faces_type.push_back(t);
              each_group.insert(each_group.end(), each_face.begin(), each_face.end());
            }
          tri2vtk(ofs_group, &input_node[0], input_node.size(2),
              &each_group[0], each_group.size()/3);
        }
      tri2vtk(ofs, &input_node[0], input_node.size(2), &out_faces[0], out_faces.size()/3);
      cell_data(ofs, &out_faces_type[0], out_faces_type.size(), "group");

      if(param){
          ofstream ofs_p("cut_tet_surface_path_param.vtk");
          tri2vtk(ofs_p, &(*param)[0], param->size(2), &out_faces_in_param[0],
              out_faces_in_param.size()/3);
          cell_data(ofs_p, &out_faces_type[0], out_faces_type.size(), "group");
        }
    }
  }
#endif

#if 1 // analysis each face group
  cerr << "# [info] left " << face_groups.size()
       << "surface groups to analysis" << endl;
  for(size_t t = 0; t < face_groups.size(); ++t){
      const set<size_t> & each_group = face_groups[t];
      analysis_each_surface_face_group(t, each_group, input_node,
                                       singularity_edges, *fa[1],
          *ea_original_surface, cut_tet2tet, param);
    }
#endif


  return 0;
}

int analysis_each_surface_face_group(
    const size_t & group_idx,
    const std::set<size_t> & one_group,
    const matrixd & node,
    const std::set<std::pair<size_t,size_t> > & singularity_edges,
    const jtf::mesh::face2tet_adjacent & fa_cut_tet,
    const jtf::mesh::edge2cell_adjacent &ea_original_surface,
    const matrixst & cut_tet2tet,
    const matrixd * param)
{
  matrixst face_group_mat(3, one_group.size());
  size_t face_idx = 0;
  for(set<size_t>::const_iterator scit = one_group.begin();
      scit != one_group.end(); ++scit, ++face_idx){
      assert(*scit < fa_cut_tet.faces_.size());
      const vector<size_t> & face = fa_cut_tet.faces_[*scit];
      copy(face.begin(), face.end(), face_group_mat(colon(), face_idx).begin());
    }

  assert(is_surface_mesh_manifold(face_group_mat));
  unique_ptr<edge2cell_adjacent> ea(edge2cell_adjacent::create(face_group_mat));

  vector<pair<size_t,size_t> > boundary_edges;
  for(size_t t = 0; t < ea->edge2cell_.size(); ++t){
      const pair<size_t,size_t> & adj_faces = ea->edge2cell_[t];
      if(ea->is_boundary_edge(adj_faces)){
          //  boundary_edges << ea->edges_[t];
          if(ea->edges_[t].first > ea->edges_[t].second)
            boundary_edges << make_pair(ea->edges_[t].second, ea->edges_[t].first);
          else
            boundary_edges << ea->edges_[t];
        }
    }

  const size_t singularity_edge_flag = 0;
  const size_t surface_edge_flag = 1;
  const size_t crease_flag = 2;

  vector<size_t> edge_on_boundary;
  vector<size_t> edge_on_boundary_in_param;
  vector<size_t> boundary_edge_type;
  for(size_t t = 0; t < boundary_edges.size(); ++t){
      const pair<size_t,size_t> & one_edge = boundary_edges[t];
      const pair<size_t,size_t> edge_original(cut_tet2tet[one_edge.first],
          cut_tet2tet[one_edge.second]);

      if(find(singularity_edges.begin(), singularity_edges.end(),
              edge_original) != singularity_edges.end() ||
         find(singularity_edges.begin(), singularity_edges.end(),
              make_pair(edge_original.second, edge_original.first))
         != singularity_edges.end()){
          edge_on_boundary << edge_original.first << edge_original.second;
          edge_on_boundary_in_param << one_edge.first << one_edge.second;
          boundary_edge_type << singularity_edge_flag;
        }else if(ea_original_surface.get_edge_idx(edge_original.first,
                                                  edge_original.second) != -1)
        { // find this edge on original surface
          edge_on_boundary << edge_original.first << edge_original.second;
          edge_on_boundary_in_param << one_edge.first << one_edge.second;
          boundary_edge_type << surface_edge_flag;
        }else{
          edge_on_boundary << edge_original.first << edge_original.second;
          edge_on_boundary_in_param << one_edge.first << one_edge.second;
          boundary_edge_type << crease_flag;
        }
    }

  if(find(boundary_edge_type.begin(), boundary_edge_type.end(), crease_flag)
     == boundary_edge_type.end()){
      cerr << "# [info] face group: " << group_idx << " is regulary." << endl;
      return 0;
    }

  ostringstream os;
  os << group_idx;
  ofstream ofs((os.str() + "_face_group_boundary.vtk").c_str());
  line2vtk(ofs, &node[0], node.size(2), &edge_on_boundary[0],
      edge_on_boundary.size()/2);
  cell_data(ofs, &boundary_edge_type[0], boundary_edge_type.size(), "flag");

  if(param){
      ofstream ofs_p((os.str() + "_face_group_boundary_param.vtk").c_str());
      line2vtk(ofs_p, &(*param)[0], param->size(2), &edge_on_boundary_in_param[0],
          edge_on_boundary_in_param.size()/2);
      cell_data(ofs_p, &boundary_edge_type[0], boundary_edge_type.size(), "flag");
    }
  return 0;
}

void analysis_transition(const jtf::mesh::meshes &input,
                         const jtf::mesh::face2tet_adjacent &input_fa,
                         const matrix<matrixd > &frame_in_tet,
                         const matrixst &cut_tet,
                         matrixst &face_pair)
{
  unique_ptr<jtf::mesh::face2tet_adjacent> output_fa(jtf::mesh::face2tet_adjacent::create(cut_tet));
  analysis_transition_raw(input.mesh_, input.node_, input_fa, *output_fa,
                          cut_tet, face_pair);
}

void analysis_transition_raw(
    const matrixst &input_tet,
    const matrixd &input_node,
    const jtf::mesh::face2tet_adjacent &input_fa,
    const jtf::mesh::face2tet_adjacent & fa_cut,
    const matrixst &cut_tet,
    matrixst &face_pair)
{
  //unique_ptr<jtf::mesh::face2tet_adjacent> output_fa(jtf::mesh::face2tet_adjacent::create(cut_tet));
  const jtf::mesh::face2tet_adjacent *fa[2] = {
    &input_fa, & fa_cut
  };
  // compact node index, setup the mapping of cut_tet2tet
  matrixst cut_tet2tet(max(cut_tet)+1);
  cut_tet2tet(cut_tet) = input_tet(colon());
  matrixst outside_face_idx_in_output;
  get_outside_face_idx(*fa[1], outside_face_idx_in_output);

  matrixd node = input_node(colon(), cut_tet2tet);

  vector<size_t> jump_face_idx_in_input;
  matrixst jump_faces_in_input;
  //! map key of original tet face to the index in outside_face_idx_in_output
  typedef map<vector<size_t>, pair<size_t, size_t> > face_pair_type;
  face_pair_type outside_face;
  {
    matrixst face(3);
    for(size_t fi = 0; fi < outside_face_idx_in_output.size(); ++fi) {
        const size_t face_id = outside_face_idx_in_output[fi];
        const pair<size_t, size_t> &nb_tet_id = fa[1]->face2tet_[face_id];
        VERIFY(fa[1]->is_outside_face(nb_tet_id), "# error in nb_tet_id.");
        const size_t tet_id = (nb_tet_id.first == -1)?nb_tet_id.second:nb_tet_id.first;

        vector<size_t> ori_face = to_vec(cut_tet2tet(fa[1]->get_face(face_id)));
        sort(ori_face.begin(), ori_face.end());

        face_pair_type::iterator i = outside_face.find(ori_face);
        if(i == outside_face.end()) {
            outside_face[ori_face] = make_pair(fi, -1);
          }
        else {
            VERIFY(outside_face[ori_face].second == -1, "# error in outside_face.");
            outside_face[ori_face].second = fi;
          }
      }
    // face_pair stores the index to outside_face_idx_in_output
    face_pair.resize(outside_face_idx_in_output.size());
    for(size_t fi = 0; fi < outside_face_idx_in_output.size(); ++fi) {
        const size_t face_id = outside_face_idx_in_output[fi];
        vector<size_t> ori_face = to_vec(cut_tet2tet(fa[1]->get_face(face_id)));
        sort(ori_face.begin(), ori_face.end());
        const pair<size_t, size_t> &id = outside_face[ori_face];
        VERIFY(id.first == fi || id.second == fi, "# outside_face map error.");
        //			face_pair[fi] = id.first+id.second-fi;
        face_pair[fi] = (id.first == fi)?id.second:id.first;
#if 1
        if(face_pair[fi] != -1) {
            matrixst oriA = cut_tet2tet(fa[1]->get_face(outside_face_idx_in_output[fi]));
            matrixst oriB = cut_tet2tet(fa[1]->get_face(outside_face_idx_in_output[face_pair[fi]]));
            sort(oriA.begin(), oriA.end());
            sort(oriB.begin(), oriB.end());
            if(norm(oriA-oriB) > 1e-8) {
                cerr << "# error in face_pair." << endl;
              }
          }
#endif
      }

    // the following code is for visualization
    set<size_t> jump_face_set_in_input;
    for(size_t fi = 0; fi < outside_face_idx_in_output.size(); ++fi) {
        if(face_pair[fi] == -1) continue;
        const size_t face_id = outside_face_idx_in_output[fi];
        const vector<size_t> &vface = fa[1]->faces_[face_id];
        vector<size_t> ori_face(3);
        for(int i = 0; i < 3; ++i)
          ori_face[i] = cut_tet2tet[vface[i]];
        jump_face_set_in_input.insert(fa[0]->get_face_idx(&ori_face[0]));
      }
    jump_face_idx_in_input.resize(jump_face_set_in_input.size());
    copy(jump_face_set_in_input.begin(), jump_face_set_in_input.end(), jump_face_idx_in_input.begin());
    jump_faces_in_input.resize(3, jump_face_idx_in_input.size());
    for(size_t fi = 0; fi < jump_face_idx_in_input.size(); ++fi) {
        copy(fa[0]->faces_[jump_face_idx_in_input[fi]].begin(),
            fa[0]->faces_[jump_face_idx_in_input[fi]].end(),
            &jump_faces_in_input(0, fi));
      }
  }

  // group jump_faces[0] by manifold part.
  map<set<size_t>, vector<size_t> > edge2face; // for the jump_faces_in_output
  for(size_t fi = 0; fi < jump_faces_in_input.size(2); ++fi) {
      for(int ei = 0; ei < 3; ++ei) {
          set<size_t> edge;
          edge.insert(jump_faces_in_input(ei, fi));
          edge.insert(jump_faces_in_input((ei+1)%3, fi));
          edge2face[edge].push_back(fi);
        }
    }
  matrixst face2group = zeros<size_t>(jump_face_idx_in_input.size(), 1)-1;
  set<size_t> ungrouped_face;
  for(size_t i = 0; i < face2group.size(); ++i)
    ungrouped_face.insert(i);
  size_t group_id = 0;
  while(!ungrouped_face.empty()) {
      size_t seed = *ungrouped_face.begin();
      ungrouped_face.erase(seed);
      stack<size_t> s;
      s.push(seed);
      while(!s.empty()) {
          size_t current = s.top();
          s.pop();
          face2group[current] = group_id;
          ungrouped_face.erase(current);
          // cerr << "# group face: " << current << " -> " << group_id << endl;
          for(int ei = 0; ei < 3; ++ei) {
              set<size_t> edge;
              edge.insert(jump_faces_in_input(ei, current));
              edge.insert(jump_faces_in_input((ei+1)%3, current));
              const vector<size_t> &face_idx = edge2face[edge];
              if(face_idx.size() == 0) {
                  cerr << "# incorrect edge2face." << endl;
                }
              if(face_idx.size() != 2)
                continue; // non-manifold or boundary
              const size_t next = face_idx[0]+face_idx[1]-current;
              if(face2group[next] != -1) continue;
              if(next == current) {
                  cerr << "# strange error." << endl;
                }
              s.push(next);
            }
        }
      ++group_id;
    }
  cerr << "# get " << group_id << " groups." << endl;

  // output for visualization
  {
    {
      if(jump_faces_in_input.size()) {
          ofstream ofs("jump_in_input.vtk");
          tri2vtk(ofs, &input_node[0], input_node.size(2), &jump_faces_in_input[0], jump_faces_in_input.size(2));
          ofs << "CELL_DATA " << face2group.size() << "\n";
          ofs << "SCALARS " << "group_id" << " float\nLOOKUP_TABLE my_table\n";
          copy(face2group.begin(), face2group.end(), ostream_iterator<size_t>(ofs, "\n"));
        }
    }
  }
}

void cut_tetmesh(
    const matrixst & tet,
    const matrixd & node,
    const jtf::mesh::face2tet_adjacent & fa,
    boost::unordered_map<pair<size_t,size_t>,size_t> & inner_face_jump_type,
    boost::unordered_map<size_t,size_t> & surface_type,
    const bool is_restrcited_type,
    matrixst & cut_tet,
    bool only_mst,
    zjucad::matrix::matrix<matrixd> *frame_ptr)
{
  Euler_characteristic ec;
  vector<set<matrixst::value_type*> > acc_table;
  fix_graph fg;
  set<pair<size_t,size_t> > spanning_edges;

  tet2dual_graph(tet,fg,fa);

  cerr << "# [info] begin to glue tets accronding a MST." << endl;
  glue_tet_mst_with_type(tet,fa, spanning_edges, fg,cut_tet,
                         ec,acc_table, inner_face_jump_type, surface_type,
                         (is_restrcited_type?0:1), frame_ptr);

  if(!only_mst){
      unique_ptr<jtf::mesh::face2tet_adjacent> fa_first(jtf::mesh::face2tet_adjacent::create(cut_tet));
      {
        matrix<matrix<double> > frame_fake(tet.size(2));
        if(frame_ptr){
            frame_fake = *frame_ptr;
          }else{
            for(size_t ti = 0; ti < frame_fake.size(); ++ti){
                frame_fake[ti] = eye<double>(3);
              }
          }
        jtf::mesh::one_ring_tet_at_edge ortae;
        ortae.add_tets(tet,fa);
        ortae.sort_into_loop(tet, node);
        jtf::mesh::dyn_one_ring_face_at_point dorfap;

        {
          if(!fa_first.get()){
              cerr << "# [error] can not build face2tet_adjacent." << endl;
              return ;
            }
          matrixst outside_face_cut_tet;
          get_outside_face(*fa_first, outside_face_cut_tet);
          unique_ptr<jtf::mesh::edge2cell_adjacent> ea(
                jtf::mesh::edge2cell_adjacent::create(outside_face_cut_tet));
          if(!ea.get()){
              cerr << "# [error] non-manifold surface after MST." << endl;
              return ;
            }
          dorfap.add_all_faces(outside_face_cut_tet);
        }

        jtf::mesh::meshes tm_mesh;
        tm_mesh.mesh_ = tet;
        tm_mesh.node_ = node;
        jtf::tet_mesh tm(tm_mesh);
        singularity_extractor se(tm);
        vector<deque<pair<size_t,size_t> > > singularity_edges;
        std::vector<deque<size_t> > singularities_type;
        se.extract(inner_face_jump_type, singularity_edges, singularities_type);

        if(dump_singularity_chain_to_vtk_2("singularity_inner_face_jump_type.vtk",
                                           node,
                                           singularity_edges, singularities_type))
          return ;
        glue_tet_avoid_singularity_crease(tet, frame_fake,fa,
                                          singularity_edges, ortae, fg, dorfap, cut_tet, ec, acc_table);
      }
    }
  std::map<size_t,size_t> mapping;
  remove_extra_node(cut_tet,mapping);
}

int glue_tet_mst_with_type(
    const matrixst &original_tet,
    const jtf::mesh::face2tet_adjacent & fa,
    std::set<std::pair<size_t,size_t> > & spanning_edge,
    fix_graph &fg,
    matrixst &output_cut_tet,
    Euler_characteristic &ec,
    std::vector<std::set<size_t *> > &acc_table,
    boost::unordered_map<pair<size_t,size_t>, size_t> & inner_face_jump_type,
    boost::unordered_map<size_t,size_t> & surface_type,
    const size_t surface_type_flag,
    zjucad::matrix::matrix<matrixd> * frame_ptr)
{
  if(output_cut_tet.size(1) != original_tet.size(1) ||
     output_cut_tet.size(2) != original_tet.size(2)){
      output_cut_tet.resize(4, original_tet.size(2));

      output_cut_tet(colon()) = colon(0, output_cut_tet.size()-1);
      acc_table.resize(output_cut_tet.size());

      for(size_t i = 0; i < acc_table.size(); ++i)
        acc_table[i].insert(&output_cut_tet[i]);
    }

  assert(original_tet.size() == output_cut_tet.size());
  assert(original_tet.size() == acc_table.size());

  ec.node_num = output_cut_tet.size(2)*4;
  ec.edge_num = output_cut_tet.size(2)*6;
  ec.face_num = output_cut_tet.size(2)*4;

  cerr << "# Euler_characteristic before glue: " << ec() << " "
       << ec.node_num << " " << ec.edge_num << " " << ec.face_num << endl;

  std::vector<double> err(fg.edge_num());
  const matrixd Identity = eye<double>(3);
  for(size_t ni = 0; ni < fg.node_num(); ++ni) {
      for(size_t ei = fg.ptr_[ni]; ei < fg.ptr_[ni+1]; ++ei){
          boost::unordered_map<pair<size_t,size_t>,size_t>::const_iterator bumcit
              = inner_face_jump_type.find(make_pair(ni, fg.idx_[ei]));
          if(bumcit == inner_face_jump_type.end())
            err[ei] = 0; // tivial type
          else
            err[ei] = norm(type_transition2(bumcit->second) - Identity);
        }
    }

  std::vector<size_t> spanning_tree;
  {
    fix_graph spanning_tree_graph;
    build_minimum_spanning_tree(fg, err, spanning_tree);
    subgraph(fg, spanning_tree, spanning_tree_graph);
    cerr << "# graph node: " << spanning_tree_graph.node_num()
         << ", edges in spanning_tree: " << spanning_tree_graph.edge_num()
         << endl;

    //orient_frame(zyz_frame_in_tet, frame, spanning_tree_graph, orient_err);

    orient_type(spanning_tree_graph, original_tet.size(2), fa,
                inner_face_jump_type, surface_type, surface_type_flag,frame_ptr);
  }

  for(size_t ei = 0; ei < spanning_tree.size(); ++ei) {
      std::pair<node_idx, node_idx> edge = fg.get_edge(spanning_tree[ei]);
      if(glue_tet(original_tet, output_cut_tet,
                  edge.first, edge.second, ec, acc_table, true) > 1) {
          cerr << "# [error] glue according spanning_tree error." << endl;
        }
    }
  cerr << "# [info] Euler_characteristic after glue MST: " << ec() << " "
       << ec.node_num << " " << ec.edge_num << " " << ec.face_num << endl;
  return 0;
}

int glue_tet_avoid_singularity_crease(
    const matrixst &original_tet,
    const matrix<matrixd > &frame,
    const jtf::mesh::face2tet_adjacent &fa,
    const vector<deque<pair<size_t,size_t> > > & singularity_edges,
    const jtf::mesh::one_ring_tet_at_edge & ortae,
    fix_graph &fg,
    jtf::mesh::dyn_one_ring_face_at_point & dorfap,
    matrixst &output_cut_tet,
    Euler_characteristic &ec,
    vector<set<size_t *> > &acc_table)
{
  assert(acc_table.size() == original_tet.size());
  assert(acc_table.size() == output_cut_tet.size());

  vector<pair<double, pair <size_t,size_t > > > possible_tets_pair;
  typedef map<pair<size_t,size_t>,size_t>::const_iterator mpscit;
  possible_tets_pair.reserve(fg.edge_num());

  for(size_t t = 0; t < fa.face2tet_.size(); ++t){
      const pair<size_t,size_t> & tet_pair = fa.face2tet_[t];
      if(fa.is_outside_face(tet_pair)) continue;
      possible_tets_pair.push_back(
            make_pair(norm(frame[tet_pair.first] - frame[tet_pair.second]),
          make_pair(tet_pair.first, tet_pair.second)));
    }

  cerr << "# possible additional edges: " << possible_tets_pair.size() << endl;
  sort(possible_tets_pair.begin(), possible_tets_pair.end());

  set<pair<size_t,size_t> > singularity_edges_set;
  for(size_t t = 0; t < singularity_edges.size(); ++t){
      const deque<pair<size_t,size_t> > & one_chain = singularity_edges[t];
      singularity_edges_set.insert(one_chain.begin(), one_chain.end());
    }

  while(1){
      const size_t possible_tets_num = possible_tets_pair.size();
      cerr << "# [info] left " << possible_tets_num << " tet pair to be glued. " << endl;
      vector<pair<double, pair<size_t,size_t> > > temp_possible_tet_pair;
      const size_t per_100 = possible_tets_num / 100;
      for(size_t t = 0; t < possible_tets_num; ++t){

          if(possible_tets_num > 1000 && t % per_100 == 0){
              cerr << "# [info] ------ " << t / per_100 << "%." << endl;
            }
          const pair<size_t,size_t> & tet_pair = possible_tets_pair[t].second;
          const int rtn = glue_each_tet_pair_avoid_singularity_edge(
                original_tet, output_cut_tet, tet_pair.first,
                tet_pair.second, ec, acc_table, false, singularity_edges_set,
                ortae, dorfap);
          //      const int rtn =
          //          glue_tet(original_tet, output_cut_tet,
          //                   tet_pair.first,tet_pair.second,
          //                   ec, acc_table, false);
          if(rtn == 2) // 0 or 1 means glued, 3 means singularity glued inside
            temp_possible_tet_pair.push_back(possible_tets_pair[t]);
        }
      if(temp_possible_tet_pair.size() == possible_tets_num)
        break;
      swap(possible_tets_pair, temp_possible_tet_pair);

#if 0// debug
      {
        matrixst outside_face;
        unique_ptr<jtf::mesh::face2tet_adjacent> fa_test(jtf::mesh::face2tet_adjacent::create(output_cut_tet));
        get_outside_face(*fa_test, outside_face);
        unique_ptr<jtf::mesh::edge2cell_adjacent> ea(jtf::mesh::edge2cell_adjacent::create(outside_face));
        if(!ea.get()){
            cerr << "# [error] non-mainfold surface of cut tet." << endl;
            return __LINE__;
          }
      }
#endif

    }

  cerr << "# [info] left " << possible_tets_pair.size()
       << " tet pair not glued." << endl;
  cerr << "# [info] Euler_characteristic after glue avoid singularity: " << ec()
       << " " << ec.node_num << " " << ec.edge_num << " " << ec.face_num << endl;
  return 0;
}

//int minimal_cut_tet(
//    const jtf::tet_mesh &tm,
//    const matrixd &zyz,
//    matrixst & cut_tet,
//    matrixst & face_pair,
//    matrixst & face_type,
//    const matrixd * param)
//{
//  tetmesh_cutter tmc(tm);
//  matrixd new_zyz = zyz;
//  tmc.cut(new_zyz, false);

//  cut_tet = tmc.cut_tm_.mesh_;

//  singularity_extractor se(tm);
//  vector<deque<pair<size_t,size_t> > > singularity_chain;
//  vector<deque<size_t> > singularity_type;
//  se.extract(new_zyz, singularity_chain, singularity_type);

//  set<pair<size_t,size_t> > singularity_edges_set;
//  for(size_t t = 0; t < singularity_chain.size(); ++t){
//      const deque<pair<size_t, size_t> > & one_chain = singularity_chain[t];
//      singularity_edges_set.insert(one_chain.begin(), one_chain.end());
//    }

//  // normalize cut tet
//  if(!is_normalized_cut_tet(tm.tetmesh_.mesh_, cut_tet)){
//      normalize_cut_tet_node_idx(tm.tetmesh_.mesh_, cut_tet);

//      if(!is_normalized_cut_tet(tm.tetmesh_.mesh_, cut_tet)){
//          cerr << "# [error] fail to normalize cut tet." << endl;
//        }else
//        cerr << "# [info] success normalize cut tet." << endl;
//    }

//  boost::unordered_map<pair<size_t,size_t>,size_t> inner_face_jump_type_after_aligned;
//  matrix<matrixd> new_frame;
//  zyz2frame(new_zyz, new_frame);
//  matrixd rot(3,3);
//  for(size_t fi = 0; fi < tm.fa_->face2tet_.size(); ++fi){
//      if(tm.fa_->is_outside_face(tm.fa_->face2tet_[fi])) continue;
//      const pair<size_t,size_t> & tet_pair = tm.fa_->face2tet_[fi];
//      get_best_alignment(&new_frame[tet_pair.first][0], &new_frame[tet_pair.second][0], &rot[0]);
//      inner_face_jump_type_after_aligned[tet_pair] = type_transition1(rot);
//      inner_face_jump_type_after_aligned[make_pair(tet_pair.second, tet_pair.first)] = type_transition1(trans(rot));
//    }

//  analysis_face_type(tm.tetmesh_.mesh_, tm.tetmesh_.node_, cut_tet, *tm.fa_,
//                     inner_face_jump_type_after_aligned,
//                     singularity_edges_set, face_pair, face_type, param);

//  return 0;
//}

int glue_tet_with_given_cut_faces(
    const matrixst &tet,
    const matrixd &node,
    const std::vector<std::pair<size_t,size_t> > & cut_tet_pair,
    matrixst &new_cut_tet)
{
  cerr << "# [info] cut_tet pair number " << cut_tet_pair.size() << endl;
  new_cut_tet.resize(tet.size(1), tet.size(2));

  new_cut_tet(colon()) = colon(0, new_cut_tet.size()-1);

  unique_ptr<jtf::mesh::face2tet_adjacent> fa(jtf::mesh::face2tet_adjacent::create(tet));
  if(!fa.get()){
      cerr << "# [error] can not buildjtf::mesh::face2tet_adjacent." << endl;
      return __LINE__;
    }

  Euler_characteristic ec;

  ec.node_num = new_cut_tet.size(2)*4;
  ec.edge_num = new_cut_tet.size(2)*6;
  ec.face_num = new_cut_tet.size(2)*4;

  vector<set<matrixst::value_type*> > acc_table;
  fix_graph fg;
  tet2dual_graph(tet,fg,*fa);

  acc_table.resize(new_cut_tet.size());

  for(size_t i = 0; i < acc_table.size(); ++i)
    acc_table[i].insert(&new_cut_tet[i]);


  cerr << "# Euler_characteristic before glue: " << ec() << " "
       << ec.node_num << " " << ec.edge_num << " " << ec.face_num << endl;


  std::vector<double> err(fg.edge_num());
  for(size_t ni = 0; ni < fg.node_num(); ++ni) {
      for(size_t ei = fg.ptr_[ni]; ei < fg.ptr_[ni+1]; ++ei){
          if(find(cut_tet_pair.begin(), cut_tet_pair.end(), make_pair(ni, fg.idx_[ei]))
             == cut_tet_pair.end() &&
             find(cut_tet_pair.begin(), cut_tet_pair.end(), make_pair(fg.idx_[ei], ni))
             == cut_tet_pair.end()){
              err[ei] = 0;// norm(tet_cf_sh(colon(), ni) - tet_cf_sh(colon(), fg.idx()[ei]));
            }else
            err[ei] = 1;
        }
    }

  std::vector<size_t> spanning_tree;
  {
    fix_graph spanning_tree_graph;
    build_minimum_spanning_tree(fg, err, spanning_tree);
    subgraph(fg, spanning_tree, spanning_tree_graph);
    cerr << "# graph node: " << spanning_tree_graph.node_num()
         << ", edges in spanning_tree: " << spanning_tree_graph.edge_num()
         << endl;
  }

  for(size_t ei = 0; ei < spanning_tree.size(); ++ei) {
      std::pair<node_idx, node_idx> edge = fg.get_edge(spanning_tree[ei]);
      if(glue_tet(tet, new_cut_tet,
                  edge.first, edge.second, ec, acc_table, true) > 1) {
          cerr << "# [error] glue according spanning_tree error." << endl;
        }
    }

  cerr << "# [info] Euler_characteristic after glue MST: " << ec() << " "
       << ec.node_num << " " << ec.edge_num << " " << ec.face_num << endl;
#if 1
  vector<pair<double, pair<size_t,size_t> > > possible_tets_pair;

  //  for(size_t t = 0; t < fa->face2tet_.size(); ++t){
  //    const pair<size_t,size_t> & tet_pair = fa->face2tet_[t];
  //    if(fa->is_outside_face(tet_pair)) continue;
  //    if(find(cut_tet_pair.begin(), cut_tet_pair.end(),
  //            make_pair(tet_pair.first, tet_pair.second))
  //       == cut_tet_pair.end() &&
  //       find(cut_tet_pair.begin(), cut_tet_pair.end(),
  //            make_pair(tet_pair.second, tet_pair.first))
  //       == cut_tet_pair.end()){
  //      possible_tets_pair.push_back(make_pair(1,tet_pair));
  //    }
  //  }

  //! TODO: why is the face number does not equal the edge number in tet_dual_graph???
  for(size_t ni = 0; ni < fg.node_num(); ++ni) {
      for(size_t ei = fg.ptr_[ni]; ei < fg.ptr_[ni+1]; ++ei){
          const size_t bi = fg.idx_[ei];
          if(ni == bi) continue;
          if(find(cut_tet_pair.begin(), cut_tet_pair.end(), make_pair(ni, bi))
             != cut_tet_pair.end()) continue;
          if(find(cut_tet_pair.begin(), cut_tet_pair.end(), make_pair(bi, ni))
             != cut_tet_pair.end()) continue;
          possible_tets_pair.push_back(make_pair(0,make_pair(ni,bi)));
        }
    }
  cerr << "# possible additional edges: " << possible_tets_pair.size() << endl;
  sort(possible_tets_pair.begin(), possible_tets_pair.end());


  //while(!possible_tets_pair.empty()){
  //const size_t possible_tets_num = possible_tets_pair.size();
  while(1){
      cerr << "# [info] left " << possible_tets_pair.size()
           << " tet pair to be glued. " << endl;
      vector<pair<double, pair<size_t,size_t> > > temp_possible_tet_pair;
      for(size_t t = 0; t < possible_tets_pair.size(); ++t){
          const pair<size_t,size_t> & tet_pair = possible_tets_pair[t].second;

          const int rtn =  glue_tet_specail(tet, new_cut_tet,
                                            tet_pair.first,tet_pair.second,
                                            ec, acc_table, false);
          if(rtn == 2)
            temp_possible_tet_pair.push_back(possible_tets_pair[t]);
        }
      if(temp_possible_tet_pair.size() == possible_tets_pair.size())
        break;
      swap(possible_tets_pair, temp_possible_tet_pair);
    }

  cerr << "# [info] left " << possible_tets_pair.size()
       << " tet pair not glued." << endl;
#endif
  cerr << "# [info] Euler_characteristic after glue avoid singularity: " << ec()
       << " " << ec.node_num << " " << ec.edge_num << " " << ec.face_num << endl;

  {
    // vis the left tet pairs
    ofstream ofs("left_faces.vtk");
    vector<size_t> faces;
    size_t face[3];

    for(size_t ti = 0; ti < possible_tets_pair.size(); ++ti){
        const pair<size_t,size_t> & tet_pair = possible_tets_pair[ti].second;
        common_tet_face(&tet(0, tet_pair.first),
                        &tet(0, tet_pair.second), face);
        faces.insert(faces.end(), face, face+3);
      }
    tri2vtk(ofs, &node[0], node.size(2), &faces[0], faces.size()/3);
  }

  return 0;
}
