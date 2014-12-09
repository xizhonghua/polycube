#ifndef SOLUTION_SEARCHING_H
#define SOLUTION_SEARCHING_H

#include <vector>
#include <map>
#include <boost/dynamic_bitset.hpp>


#include <jtflib/algorithm/gauss_elimination.h>
#include "../common/def.h"
#include "../tetmesh/tetmesh.h"
#include "../tetmesh/util.h"
//#include <jtflib/mesh/mesh.h>
#include <jtflib/mesh/mesh.h>

class singularity_graph;

//class state_each_step{
//public:
//  std::vector<double> nodes_;
//  boost::dynamic_bitset<> node_flag_;
//  std::vector<jtf::algorithm::equation<double> > eq_vec_;
//  std::vector<size_t> fnode_;
//  std::list<std::tuple<size_t,size_t,size_t> > unready_edges_;
//};

class searching_strategy{
public:
  static int prev_process(
      const matrixst & orig_tet,
      const matrixd & orig_node,
      const matrixst &cut_tet,
      const jtf::mesh::face2tet_adjacent & fa,
      const jtf::mesh::face2tet_adjacent & fa_cut,
      const matrixst & face_pair_cut,
	  const boost::unordered_map<size_t,size_t> & surface_type,
      matrixd &cut_node,
      matrixst & cut_tet2tet,
      matrixst & outside_face,
      matrixst & outside_face_cut,
      matrixst & outside_face_idx_in_cut,
      matrixst & outside_face_idx,
      std::vector<std::pair<size_t,size_t> > &jump_face_vec,
      boost::unordered_map<size_t,size_t> &surface_type_cut,
      jtf::mesh::one_ring_tet_at_edge  &ortae,
      jtf::mesh::one_ring_tet_at_edge  &ortae_cut);

  // store states on hard disk
  static int searching_solutions_by_using_disk(
	const matrixst & orig_tet,
    const matrixst & cut_tet,
    const matrixd & node,
	const jtf::mesh::face2tet_adjacent & fa,
    const jtf::mesh::face2tet_adjacent & fa_cut,
    const matrixst & face_pair_cut,
	const boost::unordered_map<size_t,size_t> & surface_type,
    boost::unordered_map<std::pair<size_t,size_t>,size_t>  &inner_face_jump_type,
    const zjucad::matrix::matrix<matrixd> * frame_ptr = 0);

  // store states on memory which restricted models to small size
  static int searching_solutions_by_using_memory(
	const matrixst & orig_tet,
    const matrixst & cut_tet,
    const matrixd & node,
	const jtf::mesh::face2tet_adjacent & fa,
    const jtf::mesh::face2tet_adjacent & fa_cut,
    const matrixst & face_pair_cut,
	const boost::unordered_map<size_t,size_t> & surface_type,
    boost::unordered_map<std::pair<size_t,size_t>,size_t>  &inner_face_jump_type,
    const zjucad::matrix::matrix<matrixd> * frame_ptr = 0);

  static int searching_solutions_by_using_memory_binary_search(
	const matrixst & orig_tet,
    const matrixst & cut_tet,
    const matrixd & node,
	const jtf::mesh::face2tet_adjacent & fa,
    const jtf::mesh::face2tet_adjacent & fa_cut,
    const matrixst & face_pair_cut,
	const boost::unordered_map<size_t,size_t> & surface_type,
    boost::unordered_map<std::pair<size_t,size_t>,size_t>  &inner_face_jump_type,
    const zjucad::matrix::matrix<matrixd> * frame_ptr = 0);
  static int adjust_face_and_type_order_to_speed_up(
      const std::vector<std::pair<size_t,size_t> > & jump_face_vec,
      const boost::unordered_map<std::pair<size_t,size_t>,size_t>  &inner_face_jump_type,
      const boost::unordered_map<size_t,size_t> & surface_type_cut,
      const jtf::mesh::face2tet_adjacent & fa_cut,
      const matrixst & cut_tet2tet,
      matrixst & type_candidates,
      boost::unordered_map<size_t,size_t> & surface_idx_to_rot_idx,
      boost::unordered_map<size_t,size_t> & surface_idx_with_rot_idx,
      boost::unordered_map<std::pair<size_t,size_t>,size_t > & face_pair_to_rot_idx,
      boost::unordered_map<size_t, std::pair<size_t,size_t> > & face_pair_with_rot_idx,
      boost::unordered_map<std::pair<size_t,size_t>,size_t> & tet_pair_to_rot_idx,
      std::vector<bool> & is_surface,
      const zjucad::matrix::matrix<matrixd> * frame_ptr = 0);

  static int searching_solutions_by_minimal_error(
    const jtf::mesh::face2tet_adjacent & fa,
    const jtf::mesh::face2tet_adjacent & fa_cut,
	const matrixst & orig_tet,
    const matrixst & cut_tet,
    const matrixd & node,
    const matrixst & face_pair_cut,
    const zjucad::matrix::matrix<matrixd> & frame,
    boost::unordered_map<std::pair<size_t,size_t>,size_t> & inner_face_jump_type);

  static int searching_solutions_by_minimal_error_binary_search(
    const jtf::mesh::face2tet_adjacent & fa,
    const jtf::mesh::face2tet_adjacent & fa_cut,
	const matrixst & orig_tet,
    const matrixst & cut_tet,
    const matrixd & node,
    const matrixst & face_pair_cut,
    const zjucad::matrix::matrix<matrixd> & frame,
    boost::unordered_map<std::pair<size_t,size_t>,size_t> & inner_face_jump_type);

  static int post_process(
      const std::vector<size_t> & rot_type,
      const matrixd & cut_node,
      const matrixd & node,
      const matrixst &outside_face,
      const jtf::mesh::face2tet_adjacent & fa_cut,
      const jtf::mesh::one_ring_tet_at_edge & ortae,
      const boost::unordered_map<size_t,size_t> & surface_idx_to_rot_idx,
      const boost::unordered_map<std::pair<size_t,size_t>,size_t> & tet_pair_to_rot_idx,
      const boost::unordered_map<std::pair<size_t,size_t>,size_t> & face_pair_to_rot_idx,
      boost::unordered_map<std::pair<size_t,size_t>,size_t> & inner_jump_type);

  static int simple_check(
      const std::vector<size_t> & rot_type,
      const std::vector<bool> & is_surface,
      const std::vector<std::pair<size_t,size_t> > & jump_face_vec,
      const matrixst & cut_tet,
      const matrixd & cut_node,
      const matrixd & node,
      const jtf::mesh::face2tet_adjacent & fa_cut,
      const matrixst &cut_tet2tet,
      const jtf::mesh::one_ring_tet_at_edge & ortae,
      const jtf::mesh::one_ring_tet_at_edge & ortae_cut,
      const matrixst & outside_face_cut,
      const matrixst & outside_face_idx_in_cut,
      const matrixst & face_pair_cut,
      const boost::unordered_map<size_t,size_t> &surface_idx_to_rot_idx,
      const boost::unordered_map<std::pair<size_t,size_t>,size_t> &face_pair_to_rot_idx,
      const boost::unordered_map<std::pair<size_t,size_t>,size_t> &tet_pair_to_rot_idx,
      bool no_surface = false);
  /// cluster the inner faces into several groups to decrease the searching space
  static int group_faces(
      const matrixst & cut_tet,
      const matrixd & cut_node,
      const jtf::mesh::face2tet_adjacent & fa_cut,
      const matrixst & outside_face_cut,
      const matrixst & outside_face_cut_idx,
      const matrixst & face_pair_cut,
      const zjucad::matrix::matrix<matrixd> & frame,
      std::vector<std::vector<std::pair<size_t,size_t> > > &jump_face_groups,
      std::vector<std::vector<size_t> > &surface_groups,
      std::vector<std::pair<size_t, double> > &group_type_with_frame_reliablity,
      std::vector<std::pair<size_t, double> > &group_type_with_surface_reliablity,
      const double frame_reliablity_threshold = 1);

  static int group_faces_according_face_reliablity(
      const matrixst & cut_tet,
      const matrixd & cut_node,
      const jtf::mesh::face2tet_adjacent & fa_cut,
      const matrixst & outside_face_cut,
      const matrixst & outside_face_cut_idx,
      const matrixst & face_pair_cut,
      const zjucad::matrix::matrix<matrixd> & frame,
      std::vector<std::vector<std::pair<size_t,size_t> > > &jump_face_groups,
      std::vector<std::vector<size_t> > &surface_groups,
      boost::unordered_map<std::pair<size_t,size_t>,size_t> & face_pair_to_rot_idx,
      boost::unordered_map<std::pair<size_t,size_t>,size_t> & tet_pair_to_rot_idx,
      boost::unordered_map<size_t,size_t> & surface_to_rot_idx,
      matrixst & type_candidates,
      std::vector<bool> & is_surface,
      const double frame_reliablity_threshold = 1);

  static int group_faces_according_face_reliablity_for_minimal_error(
      const matrixst & cut_tet,
      const matrixd & cut_node,
      const jtf::mesh::face2tet_adjacent & fa_cut,
      const matrixst & outside_face_cut,
      const matrixst & outside_face_cut_idx,
      const matrixst & face_pair_cut,
      const zjucad::matrix::matrix<matrixd> & frame,
      std::vector<std::vector<std::pair<size_t,size_t> > > &jump_face_groups,
      std::vector<std::vector<size_t> > &surface_groups,
      boost::unordered_map<std::pair<size_t,size_t>,size_t> & face_pair_to_rot_idx,
      boost::unordered_map<std::pair<size_t,size_t>,size_t> & tet_pair_to_rot_idx,
      boost::unordered_map<size_t,size_t> & surface_to_rot_idx,
      std::vector<std::deque<std::pair<double,size_t> > > & diff,
      std::vector<bool> & is_surface,
      const double frame_reliablity_threshold = 1);

  static int grow_mst_for_jump_faces(
      const jtf::mesh::face2tet_adjacent & fa_cut,
      const std::pair<size_t,size_t> & inner_face_pair,
      const boost::unordered_map<size_t,size_t> & outside_face_idx2idx_in_vec,
      const double &frame_reliablity_threshold,
      const jtf::mesh::edge2cell_adjacent & ea,
      const zjucad::matrix::matrix<matrixd> & frame,
      const matrixst & outside_face_idx_cut,
      const boost::unordered_map<size_t,size_t> & jump_face_pair,
      std::vector<std::pair<size_t,size_t> > &one_jump_face_group,
      std::vector<bool> &visited_face_flag,
      size_t & type,
      double & frame_reliablity);
  static int grow_mst_for_surface_faces(
      const size_t & surface_idx,
      const double & frame_reliablity_threshold,
      const jtf::mesh::edge2cell_adjacent &ea,
      const zjucad::matrix::matrix<matrixd> & frame,
      const matrixst & cut_tet,
      const matrixd & cut_node,
      const jtf::mesh::face2tet_adjacent & fa_cut,
      const matrixst & outside_face_cut_idx,
      const boost::unordered_map<size_t,size_t> &surface_to_idx_in_vec,
      std::vector<bool> & visited_face_flag,
      std::vector<size_t> &one_surface_group,
      size_t &type,
      double &frame_reliablity);

  static int cal_frame_difference_from_face_idx(
      const size_t & face_from,
      const size_t & face_to,
      const jtf::mesh::face2tet_adjacent & fa_cut,
      const zjucad::matrix::matrix<matrixd> & frame,
      std::vector<std::pair<double,size_t> > & frame_difference);

  static int cal_frame_difference(
      const size_t & tet_from,
      const size_t & tet_to,
      const zjucad::matrix::matrix<matrixd> & frame,
      std::vector<std::pair<double,size_t> > & frame_difference);

  static int cal_frame_to_normal_difference(
      const size_t & face_idx,
      const matrixst & cut_tet,
      const matrixd & cut_node,
      const jtf::mesh::face2tet_adjacent & fa_cut,
      const zjucad::matrix::matrix<matrixd> & frame,
      std::vector<std::pair<double, size_t> > & frame_to_normal_difference);

  static int searching_solutions_by_using_memory_group(
    const matrixst & orig_tet,
    const matrixst & cut_tet,
    const matrixd & node,
    const jtf::mesh::face2tet_adjacent & fa,
    const jtf::mesh::face2tet_adjacent & fa_cut,
    const matrixst & face_pair_cut,
    const boost::unordered_map<size_t,size_t> & surface_type,
    boost::unordered_map<std::pair<size_t,size_t>,size_t>  &inner_face_jump_type,
    const zjucad::matrix::matrix<matrixd> * frame_ptr = 0);

  static int searching_solutions_by_using_memory_binary_search_group(
    const matrixst & orig_tet,
    const matrixst & cut_tet,
    const matrixd & node,
    const jtf::mesh::face2tet_adjacent & fa,
    const jtf::mesh::face2tet_adjacent & fa_cut,
    const matrixst & face_pair_cut,
    const boost::unordered_map<size_t,size_t> & surface_type,
    boost::unordered_map<std::pair<size_t,size_t>,size_t>  &inner_face_jump_type,
    const zjucad::matrix::matrix<matrixd> * frame_ptr = 0);

  static int searching_solutions_by_using_memory_binary_search_group_polycube(
    const matrixst & orig_tet,
    const matrixst & cut_tet,
    const matrixd & node,
    const jtf::mesh::face2tet_adjacent & fa,
    const jtf::mesh::face2tet_adjacent & fa_cut,
    const matrixst & face_pair_cut,
    const boost::unordered_map<size_t,size_t> & surface_type,
    boost::unordered_map<std::pair<size_t,size_t>,size_t>  &inner_face_jump_type,
    const zjucad::matrix::matrix<matrixd> * frame_ptr = 0);

  static int searching_solutions_by_minimal_error_search_group(
      const jtf::mesh::face2tet_adjacent & fa,
      const jtf::mesh::face2tet_adjacent & fa_cut,
      const matrixst & orig_tet,
      const matrixst & cut_tet,
      const matrixd & node,
      const matrixst & face_pair_cut,
      const zjucad::matrix::matrix<matrixd> & frame,
      boost::unordered_map<std::pair<size_t,size_t>,size_t> & inner_face_jump_type);
  static int orig_check_graph(
      const matrixst & orig_tet,
      const matrixst & cut_tet,
      const matrixd & node,
      const jtf::mesh::face2tet_adjacent & fa,
      const jtf::mesh::face2tet_adjacent & fa_cut,
      const matrixst & face_pair_cut,
      const boost::unordered_map<size_t,size_t> & surface_type, // orig surface type
      boost::unordered_map<std::pair<size_t,size_t>,size_t>  &inner_face_jump_type,
      const bool no_surface,
      const zjucad::matrix::matrix<matrixd> * frame_ptr = 0);
};

#endif
