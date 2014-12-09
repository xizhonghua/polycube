#ifndef TOPOLOGY_ANALYSIS_H
#define TOPOLOGY_ANALYSIS_H

#include <deque>
#include <set>
#include <stack>
#include <zjucad/matrix/matrix.h>
#include <boost/property_tree/ptree.hpp>
#include <boost/unordered_set.hpp>

#include <jtflib/mesh/mesh.h>
#include <jtflib/algorithm/gauss_elimination.h>

#include "../tetmesh/tetmesh.h"
#include "../tetmesh/util.h"
#include "../hex_param/find_singularities.h"

class singularity_graph;

/**
 * @brief extract jump type from a spanning tree, and compute a suitable type
 *        as gluing. Each type should be checked whether it will result in
 *        compound singularity or global conflicts, if no, choose it, or choose
 *        other options. If can not find an appropriate type, go back and
 *        research.
 *
 * @param orig_mesh original tet mesh
 * @param node      original tet node
 * @param frame     cross frame
 * @param face_type output face type setting
 * @return int      if works fine return 0, or return non-zeros
 */
int aggressive_extract_jump_type(
    const matrixst & orig_tet,
    const matrixd & node,
    const matrixd & zyz,
    boost::unordered_map<std::pair<size_t,size_t>,size_t> & face_type);


int aggressive_assign_transition(
    const matrixst & orig_tet,
    const matrixd & node,
    const matrixd & zyz,
    const bool no_surface, // whether need surface constraints
    boost::unordered_map<std::pair<size_t,size_t>,size_t> & face_type);

int aggressive_assign_transition_with_origin_type(
    const matrixst & orig_tet,
    const matrixd & node,
    boost::unordered_map<std::pair<size_t,size_t>,size_t> & face_type,
    boost::property_tree::ptree & pt,
    const bool no_surface, // whether need surface constraints
    matrixd * zyz_ptr = 0);

int aggressive_assign_transition_with_type(
    const matrixst & orig_tet,
    const matrixst & cut_tet,
    const matrixd & node,
    boost::unordered_map<std::pair<size_t,size_t>,size_t> & face_type,
    boost::unordered_map<size_t,size_t> & surface_type,
    const bool no_surface,
    const zjucad::matrix::matrix<matrixd> * frame_ptr = 0 );

int detect_topology_graph_valid(
    const std::vector<size_t> &rot_type,
    const matrixst & cut_tet,
    const matrixst & cut_outside_face,
    const matrixst & cut_outside_face_idx,
    const matrixst & cut_tet2tet,
    const jtf::mesh::face2tet_adjacent & fa_cut,
    const matrixd & node,
    const matrixd & cut_node,
    const std::vector<std::pair<size_t,size_t> > & jump_face_vec,
    const matrixst &  face_pair_cut,
    const jtf::mesh::one_ring_tet_at_edge & ortae,
    singularity_graph & g,
    const std::vector<bool> & is_surface,
    const boost::unordered_map<std::pair<size_t,size_t>,size_t> &tet_pair2rot_idx,
    const boost::unordered_map<std::pair<size_t,size_t>,size_t> & face_pair_to_rot_type_idx,
    const boost::unordered_map<size_t,size_t> & surface_idx_to_rot_idx,
    bool no_surface = false);

/**
 * @brief convert surface rotation type to axis align type
 *        surface align type stands for a rotation which align the u axis to
 *        normal
 *
 * @param type input rotation type 0-23
 * @return int
 */
int convert_surface_rotation_to_axis_type(const size_t &type);

/**
 * @brief set up integer constraint per edge
 *
 * @param pt
 * @param fa
 * @param tet
 * @param node
 * @param outside_face
 * @param outside_face_idx
 * @param std::vector<std::deque<std::pair<size_t
 * @param singularities_chain
 * @param singularities_type
 * @return int
 */
int calc_singularities_integer_constraint_per_edge(
    boost::property_tree::ptree &pt,
    boost::unordered_map<std::pair<size_t,size_t>,size_t> inner_face_jump_type,
    const jtf::mesh::face2tet_adjacent &fa,
    const jtf::mesh::one_ring_tet_at_edge & ortae,
    const matrixst &tet,
    const matrixd &node,
    const matrixst & outside_face,
    const matrixst & outside_face_idx,
    const std::vector<std::deque<std::pair<size_t,size_t> > > &singularities_chain,
    const std::vector<std::deque<size_t> > &singularities_type);

/**
 * @brief
 *
 * @param pt
 * @param fa
 * @param ea
 * @param orfap
 * @param tet
 * @param node
 * @param outside_face
 * @param outside_face_idx
 * @param std::vector<std::deque<std::pair<size_t
 * @param singularities_chain
 * @param singularities_type
 * @param outside_face_normal
 * @return int
 */
int calc_surface_singularities_integer_constraint(
    boost::property_tree::ptree &pt,
    const jtf::mesh::face2tet_adjacent &fa,
    const jtf::mesh::edge2cell_adjacent &ea,
    const jtf::mesh::one_ring_face_at_point & orfap,
    const matrixst &tet,
    const matrixd &node,
    const matrixst & outside_face,
    const matrixst & outside_face_idx,
    const std::vector<std::deque<std::pair<size_t,size_t> > > &singularities_chain,
    const std::vector<std::deque<size_t> > &singularities_type,
    const matrixd & outside_face_normal);

/**
 * @brief cal jump type between two frame, and sort the geometry error from small
 *        to large
 *
 * @param tet_pair  input tet pair
 * @param frame     input cross frame
 * @param priority_queue output queue sorted by geometry error
 * @return int  if works fine return 0, else return non-zeros
 */
int cal_frame_geometry_priority_queue(
    const std::pair<size_t,size_t>& tet_pair,
    const zjucad::matrix::matrix<matrixd > & frame,
    std::vector<std::pair<double,size_t> > & priority_queue);

/**
 * @brief find a closed face patch  on surface with two closed boundary, return
 * the smaller one, which is a simple connected one limitions: can only work on
 * closed manifold surfaec
 *
 * @param node
 * @param boundary0
 * @param boundary1
 * @param ea
 * @param outside_face
 * @param smaller_face_patch_with_boundary
 * @return int
 */
int find_face_patch_with_two_closed_boundary(
    const matrixd & node,
    const std::vector<size_t> & boundary0,
    const std::vector<size_t> & boundary1,
    const jtf::mesh::edge2cell_adjacent &ea,
    const matrixst &outside_face,
    std::vector<size_t> &smaller_face_patch_with_boundary);


//
bool is_edge_type_accetable(const size_t & edge_type);


/**
 * @brief This function is used to check whether current configuration with
 *        attempt modification will result in conflicts
 *
 * @param face_type         current
 * @param attempt_tet_pair
 * @param attempt_face_type
 * @return int
 */
bool is_singularity_configuration_accetable(
    const std::map<std::pair<size_t,size_t>,size_t> & face_type,
    const std::pair<size_t,size_t> & attempt_tet_pair,
    const size_t  & attempt_face_type);
/////////////////////////////////////////////////////////////////////
///////////////////////////////// IO ////////////////////////////////

/**
 * @brief  this function is used to dump out singulairty fatal degenerated points
 *
 * @param filename output filename
 * @param degenerated_points  store the degenerated points, each vector stands
 *                            for the same type of degeneration
 * @return int if works well, return 0, or non-zeros
 */
int dump_singularity_degenerated_points(
    const char * filename,
    const std::vector<std::vector<size_t> > & degenerated_points);

/**
 * @brief  this function is used to load singulairty fatal degenerated points
 *
 * @param filename input filename
 * @param degenerated_points  store the degenerated points, each vector stands
 *                            for the same type of degeneration
 * @return int if works well, return 0, or non-zeros
 */
int load_singularity_degenerated_points(
    const char * filename,
    std::vector<std::vector<size_t> > & degenerated_points);

/**
 * @brief
 *
 * @param filename
 * @param std::map<std::pair<size_t
 * @param size_t>
 * @param integer_constraints
 * @return int
 */
int load_integer_constraint(
    const char* filename,
    std::map<std::pair<size_t,size_t>,std::set<size_t> > & integer_constraints);


/**
 * @brief
 *
 * @param filename
 * @param std::map<std::pair<size_t
 * @param size_t>
 * @param integer_cons
 * @return int
 */
int load_integer_constraint_linking_info(
    const char* filename,
    std::map<std::pair<size_t,size_t>,std::set<size_t> > & integer_cons);


/**
 * @brief
 *
 * @param filename
 * @param std::map<std::pair<size_t
 * @param size_t>
 * @param integer_cons
 * @param group
 * @return int
 */
int load_integer_constraint_group_info(
    const char* filename,
    boost::unordered_map<std::pair<size_t,size_t>,size_t> & integer_cons,
    std::vector<std::set<size_t> > & group);


/**
 * @brief
 *
 * @param filename
 * @param std::map<std::pair<size_t
 * @param size_t>
 * @param inner_face_jump_type
 * @return int
 */
int dump_inner_face_jump_type(
    const char * filename,
    const boost::unordered_map<std::pair<size_t,size_t>,size_t> & inner_face_jump_type);

template <typename OS,typename FACE_Iterator, typename TYPE_Iterator, typename INT>
int dump_surface_axis_type(OS &os, FACE_Iterator face_begin,
                           INT face_num, TYPE_Iterator type_begin,INT face_points = 3)
{
  os << face_num << std::endl;
  for(size_t fi = 0; fi < face_num; ++fi, ++type_begin){
    for(size_t pi = 0; pi < face_points; ++pi, ++face_begin)
      os << *face_begin << " ";
    os << *type_begin << std::endl;
  }
  return 0;
}
/**
 * @brief
 *
 * @param filename
 * @param tet_rot
 * @return int
 */
int dump_tet_rotation_from_array(
    const char* filename,
    const std::vector<size_t> & tet_rot);


/**
 * @brief dump out tet rot type, base is tet_0
 *
 * @param filename
 * @param tet_rot
 * @return int return 0 if works ok, or non-zero
 */
int dump_tet_rotation_from_raw(
    const char* filename,
    const size_t *tet_rot,
    const size_t tet_num);

/**
 * @brief
 *
 * @param filename
 * @param singularity_point_uvw_ptr
 * @param idx_to_real_point
 * @return int
 */
int dump_inner_singularity_integer_cons(
    const char * filename,
    const std::vector<std::set<size_t*> > &singularity_point_uvw_ptr,
    const matrixst & idx_to_real_point);


/**
 * @brief
 *
 * @param filename
 * @param uvw_of_point_on_surface_ptr
 * @param idx_to_real_point
 * @return int
 */
int dump_surface_integer_cons(
    const char * filename,
    const zjucad::matrix::matrix<size_t*> & uvw_of_point_on_surface_ptr,
    const matrixst & idx_to_real_point);


/**
 * @brief
 *
 * @param filename
 * @param uvw_of_point
 * @param idx_to_real_point
 * @return int
 */
int dump_union_integer_constraints(
    const char * filename,
    const matrixst & uvw_of_point,
    const matrixst & idx_to_real_point);


/**
 * @brief
 *
 * @param filename
 * @param tet_rot
 * @return int
 */
int load_tet_rotation(const char * filename,
                      matrixst & tet_rot);


/**
 * @brief
 *
 * @param filename
 * @param tet_rot
 * @return int
 */
int load_tet_rotation_array(const char* filename,
                            std::vector<size_t> & tet_rot);
/////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////


/**
 * @brief
 *
 * @param tet
 * @param std::map<std::pair<size_t
 * @param size_t>
 * @param inner_face_jump_type
 * @param tet_rot
 * @return int
 */
int calc_tet_rotation(
    const jtf::tet_mesh & tm,
    const boost::unordered_map<std::pair<size_t,size_t>,size_t> & inner_face_jump_type,
    std::vector<size_t> & tet_rot);




/**
 * @brief
 *
 * @param cut_tet
 * @param fa_cut
 * @param std::map<std::pair<size_t
 * @param size_t>
 * @param inner_face_jump_type
 * @param tet_rot
 * @return int
 */
int get_tet_rot_type(
    const matrixst & cut_tet,
    const jtf::mesh::face2tet_adjacent &fa_cut,
    const boost::unordered_map<std::pair<size_t,size_t>,size_t> & inner_face_jump_type,
    matrixst & tet_rot);

int get_tet_rot_type_vec(
    const matrixst & cut_tet,
    const jtf::mesh::face2tet_adjacent &fa_cut,
    const boost::unordered_map<std::pair<size_t,size_t>,size_t> & inner_face_jump_type,
    std::vector<size_t> & tet_rot);
///**
// * @brief
// *
// * @param tet
// * @param std::map<std::pair<size_t
// * @param size_t>
// * @param inner_face_jump_type
// * @param fa
// * @param cut_tet
// * @return int
// */
//int cut_tet_accronding_to_face_type(
//    const matrixst & tet,
//    const boost::unordered_map<std::pair<size_t,size_t> ,size_t> & inner_face_jump_type,
//    const jtf::mesh::face2tet_adjacent &fa,
//    matrixst & cut_tet);


/**
 * @brief
 *
 * @param pt
 * @param fa
 * @param tet
 * @param node
 * @param ortae
 * @param std::map<std::pair<size_t
 * @param size_t>
 * @param inner_face_jump_type
 * @param ea
 * @param outside_face
 * @param outside_face_idx
 * @param std::vector<std::deque<std::pair<size_t
 * @param singularity_chain
 * @param singularity_type
 * @return int
 */
int check_degeneration_under_inner_singularity_and_surface(
    boost::property_tree::ptree & pt,
   jtf::mesh::face2tet_adjacent & fa,
    std::vector<size_t> & tet,
    std::vector<double> & node,
    jtf::mesh::one_ring_tet_at_edge & ortae,
    boost::unordered_map<std::pair<size_t,size_t>,size_t> inner_face_jump_type,
    const jtf::mesh::edge2cell_adjacent & ea,
    const matrixst & outside_face,
    const matrixst & outside_face_idx,
    const std::vector<std::deque<std::pair<size_t,size_t> > > &singularity_chain,
    const std::vector<std::deque<size_t> > & singularity_type);

/**
 * @brief this function is used to remove inner degeneration loop by modify inner
 *        face jump type
 *
 * @param loops   each deque make a loop chain, changing the type of faces
 *                arounded by the loop is able to remove the loop
 * @param inner_face_jump_type  store the inner face jump type
 * @param tet     tet matrix
 * @param node    node matrix
 * @return int    if works fine reture 0, or return non-zeros
 */
int remove_singularity_degenerated_edge(
    const std::vector<std::deque<std::pair<size_t,size_t> > > &loops,
    boost::unordered_map<std::pair<size_t,size_t>,size_t> & inner_face_jump_type,
    const matrixst &tet,
    const matrixd &node,
    const jtf::mesh::one_ring_tet_at_edge & ortae,
    const jtf::mesh::face2tet_adjacent &fa);

/**
 * @brief this function is used to remove inner degeneration point by modify inner
 *        face jump type, the main method is to break the loop which means to find
 *        two linkeds edge(whose fixed axes are the same) at given degenerated point,
 *        and by modifying the face type to move thoes two linked edges avoiding
 *        given degenerated point
 *
 * @param point_idx   degenerated point
 * @param linked_edges   edges linked at given degenerated point
 * @param singularity_edge_with_type  a map which contains all singularity edge and
 *                                    its type
 * @param inner_face_jump_type     inner_face_jump_type
 * @param tet    tet
 * @return fa    if works fine reture 0, or return non-zeros
 */
int remove_singularity_degenerated_point(
    const size_t &point_idx,
    const std::set<std::pair<size_t,size_t> > & linked_edges,
    const std::map<std::pair<size_t,size_t>,size_t> &singularity_edge_with_type,
    boost::unordered_map<std::pair<size_t,size_t>, size_t> &inner_face_jump_type,
    const matrixst &tet,
    const jtf::mesh::face2tet_adjacent &fa);
/**
 * @brief this function is used to remove one inner loop degeneration loop by modify inner
 *        face jump type
 *
 * @param loops   each deque make a loop chain, changing the type of faces
 *                arounded by the loop is able to remove the loop
 * @param inner_face_jump_type  store the inner face jump type
 * @param tet     tet matrix
 * @param node    node matrix
 * @return int    if works fine reture 0, or return non-zeros
 */

int remove_one_inner_loop_by_modify_face_type(
    const std::deque<std::pair<size_t,size_t> > & one_loop,
    boost::unordered_map<std::pair<size_t,size_t>,size_t> &inner_face_jump_type,
    const jtf::mesh::face2tet_adjacent & fa,
    const jtf::mesh::one_ring_tet_at_edge & ortae,
    const matrixst &tet,
    const matrixd &node);

/**
 * @brief this function is used to remove degeneration points by modify inner face jump type and surface type
 *
 * @param degenerated_points
 * @param outside_face
 * @param std::vector<std::deque<std::pair<size_t
 * @param singularity_chain
 * @param singularity_type
 * @param outside_face_normal_align_type
 * @param std::map<std::pair<size_t
 * @param size_t>
 * @param inner_face_jump_type
 * @return int
 */

int remove_union_degeneration_points(
    const std::vector<std::vector<size_t> > & degenerated_points,
    const matrixst &outside_face,
    const std::vector<std::deque<std::pair<size_t,size_t> > > &singularity_chain,
    const std::vector<std::deque<size_t> > &singularity_type,
    matrixst &outside_face_normal_align_type,
    boost::unordered_map<std::pair<size_t,size_t>,size_t> & inner_face_jump_type);

/**
 * @brief
 *
 * @param degenerated_points
 * @param outside_face
 * @param std::vector<std::deque<std::pair<size_t
 * @param singularity_chain
 * @param singularity_type
 * @param uvw_integer
 * @param ea
 * @param fa
 * @param new_tet
 * @param new_node
 * @param ortae
 * @param outside_face_normal_align_type
 * @param std::map<std::pair<size_t
 * @param size_t>
 * @param inner_face_jump_type
 * @return int
 */

int remove_union_degeneration_edge(
    const std::vector<std::vector<size_t> > & degenerated_points,
    const matrixst &outside_face,
    const std::vector<std::deque<std::pair<size_t,size_t> > > &singularity_chain,
    const std::vector<std::deque<size_t> > &singularity_type,
    const matrixst & uvw_integer,
    const jtf::mesh::edge2cell_adjacent & ea,
   jtf::mesh::face2tet_adjacent & fa,
    std::vector<size_t> & new_tet,
    std::vector<double> & new_node,
    jtf::mesh::one_ring_tet_at_edge & ortae,
    matrixst &outside_face_normal_align_type,
    boost::unordered_map<std::pair<size_t,size_t>,size_t> & inner_face_jump_type);

/**
 * @brief this function is used to combine singularity chain
 *
 * @param std::vector<std::deque<std::pair<size_t
 * @param singularities_chain
 * @param singularities_type
 * @return int
 */
int combine_singularity_chain(
    std::vector<std::deque<std::pair<size_t,size_t> > >  &singularities_chain,
    std::vector<std::deque<size_t> > &singularities_type);


int check_cut_tet_inner_after_aligned(
    const jtf::mesh::face2tet_adjacent &fa_cut,
    const boost::unordered_map<std::pair<size_t,size_t>,size_t> &inner_face_jump_type);

int check_tet_rotation_type_compatible_with_jump_type(
    const std::vector<size_t> & tet_rot_type,
    const jtf::mesh::face2tet_adjacent & fa_cut,
    const boost::unordered_map<std::pair<size_t,size_t>,size_t> & inner_face_jump_type);


///**
// * @brief this function is used to align cut_tet to keep the inner face jump type
// *        become identity type.
// *
// * @param std::vector<std::deque<std::pair<size_t
// * @param singularities_chain
// * @param singularities_type
// * @return int
// */
//int global_align_cut_tet(
//    const jtf::mesh::face2tet_adjacent & fa_original,
//    const jtf::mesh::face2tet_adjacent & fa_cut,
//    std::map<std::pair<size_t,size_t>, size_t> & inner_face_jump_type,
//    zjucad::matrix::matrix<matrixd > & frame );
/**
 * @brief
 *
 * @param faces
 * @return bool
 */

bool is_simple_connected_patch(const matrixst & faces);

/**
 * @brief
 *
 * @param faces
 * @param outside_face_type
 * @return bool
 */

bool is_patch_with_same_axis_fix(
    const std::vector<size_t> & faces,
    const matrixst &outside_face_type);

/**
 * @brief
 *
 * @param degenerated_points
 * @param linked_near_surface_chain
 * @return bool
 */

bool is_near_surface_degenerated_case(
    const std::vector<size_t> & degenerated_points,
    std::vector<size_t> & linked_near_surface_chain);

/**
 * @brief
 *
 * @param same_coord_points_vec
 * @param uvw_integer
 * @param std::map<size_t
 * @param real_point_to_idx
 * @param singularity_points
 * @param std::vector<std::deque<std::pair<size_t
 * @param singularity_chain
 * @param singularity_type
 * @param outside_face
 * @return bool
 */

bool is_union_edge_degenerated_case(
    const std::vector<size_t> &same_coord_points_vec,
    const matrixst &uvw_integer,
    const boost::unordered_map<size_t,size_t> & real_point_to_idx,
    const std::set<size_t> &singularity_points,
    const std::vector<std::deque<std::pair<size_t,size_t> > > &singularity_chain,
    const std::vector<std::deque<size_t> > &singularity_type,
    const matrixst & outside_face);


/**
 * @brief check whether these points can be linked as a singularity loop
 *
 * @param same_coord_points                   stores the points with the integer coordinates fixed
 * @param singularity_edge_map_to_chain_idx   all inner singularity edge map to chain idx, edge is stored as <1,2>,
 *                                            the smaller is first, the bigger the second
 * @param loop_chains                         if there points can be linked as a singularity loop,
 *                                            loop_chains store the loop edges, each deque stands for a chain,
 *                                            otherwise it's empty
 * @return bool                               if points can be linked as a singularity loop,
 *                                            return true, or false
 */

bool is_exist_inner_loop(
    const std::vector<size_t> &same_coord_points,
    const boost::unordered_map<std::pair<size_t,size_t>,size_t> &singularity_edge_map_to_chain_idx,
    std::vector<std::deque<std::pair<size_t,size_t> > > & loop_chains);

int extract_inner_face_type_from_zyz(
    const jtf::mesh::face2tet_adjacent & fa,
    const matrixd & zyz,
    boost::unordered_map<std::pair<size_t,size_t>, size_t> & inner_face_jump_type);

int extract_inner_face_type_from_frame(
    const jtf::mesh::face2tet_adjacent & fa,
    const zjucad::matrix::matrix<matrixd > & frame,
    boost::unordered_map<std::pair<size_t,size_t>, size_t> & inner_face_jump_type);

/**
 * @brief         get surface integer fix coordinate axis
 *
 * @param type    input surface face normal aligned type
 * @return size_t axis (0,1,2) stands for u,v,w axis, or return -1 if input type in invalid
 */
size_t get_surface_fix_axis(const size_t type);

/**
 * @brief This function is used to gradually fix reliable constraints according
 *  to the inner face jump type
 * @param cut_tet the input cut tet mesh
 * @param frame input oriented frame
 * @param node the input node, it may be modified if parameterization is used
 * @param inner_face_jump_type input inner face jump type
 * @return int  return 0 if works well, or return non-zero
 */
int gradually_fix_reliable_constraints(
    const matrixst & cut_tet,
    matrixd &node,
    boost::unordered_map<std::pair<size_t,size_t>,size_t> &inner_face_jump_type,
    boost::unordered_map<size_t,size_t> & surface_type);

/**
 * @brief This function is used to remove near miss, it will first find a patch
 *        along each chain, and allow the type flipping operation applied on
 *        faces, but it can not remove all near miss, since:
 *         1. lack of clear definition of near miss
 *         2. can not find a simple patch along the singularity: these g_unknown
 *            faces are not manifold
 * @param cut_tet the input cut tet mesh
 * @param tet   input tet
 * @param cut_node input cut_node
 * @param node input node
 * @param inner_face_jump_type  input inner_face_jump_type, will be modified
 * @param surface_type  input surface type
 * @param g_unknown_face_pair input g_unknown_face_pair
 * @param fa  inputjtf::mesh::face2tet_adjacent for original tet
 * @param fa_cut  inputjtf::mesh::face2tet_adjacent for original cut tet
 * @param fnode input fnode coming from equation graph, fnode is used to detect near miss
 * @return int  return 0 if works well, or return non-zero
 */
int remove_near_miss_by_minimal_cut_new(
    const matrixst &cut_tet,
    const matrixst &tet,
    const matrixd & cut_node,
    const matrixd & node,
    boost::unordered_map<std::pair<size_t,size_t>,size_t> &inner_face_jump_type,
    boost::unordered_map<size_t,size_t> &surface_type,
    const std::vector<std::pair<size_t,size_t> > & g_unknown_face_pair,
    const jtf::mesh::face2tet_adjacent &fa,
    const jtf::mesh::face2tet_adjacent &fa_cut,
    const std::vector<size_t> & fnode);


//! @brief this function is used to remove near miss, but it can not remove all
//  near miss
int remove_near_miss_by_minimal_cut_old(
    const matrixst &cut_tet,
    const matrixst &tet,
    const matrixd & cut_node,
    const matrixd & node,
    boost::unordered_map<std::pair<size_t,size_t>,size_t> &inner_face_jump_type,
    boost::unordered_map<size_t,size_t> &surface_type,
    const std::vector<std::pair<size_t,size_t> > & g_unknown_face_pair,
    const jtf::mesh::face2tet_adjacent &fa,
    const jtf::mesh::face2tet_adjacent &fa_cut,
    const std::vector<size_t> & fnode);

/**
 * @brief This function is used to remove near miss for each singularity chain
 *        by given patch, it will go through the patch and flip each face to remove
 *        singulariy, but this process is not able to remove all singularity
 * @param one_chain the input one_chain
 * @param fa input face2tet adjacent
 * @param ortae input one_ring_tet_around_edge
 * @param orfae input one_ring_face_around_edge
 * @param ea input edge2triangle adjacent
 * @param all_singularity_edges input  all singularity edge gathered in set
 * @param g_unknown_face_in_original_without_patch
 * @param g_unknown_face_in_original
 * @param inner_face_jump_type
 * @param outside_face
 * @param node
 * @return int  return 0 if works well, or return non-zero
 */
int remove_near_miss_chain_by_modifying_g_unknown_faces_new(
    const std::deque<std::pair<size_t,size_t> > & one_chain,
    const jtf::mesh::face2tet_adjacent &fa,
    const jtf::mesh::one_ring_tet_at_edge & ortae,
    const jtf::mesh::one_ring_face_at_edge & orfae,
    const jtf::mesh::edge2cell_adjacent & ea, // surface edge
    const boost::unordered_set<std::pair<size_t,size_t> > & all_singularity_edges,
    const boost::unordered_set<size_t> & g_unknown_face_in_original_without_patch,
    boost::unordered_set<size_t> & g_unknown_face_in_original,
    boost::unordered_map<std::pair<size_t,size_t>,size_t> & inner_face_jump_type,
    const matrixst & outside_face,
    const matrixd & node);

//! @brief: this function will left some chains, sinece these ones can not be
//          removed only by type flipping on candidate faces
int remove_near_miss_chain_by_modifying_g_unknown_faces_old(
    const std::deque<std::pair<size_t,size_t> > & one_chain,
    const jtf::mesh::face2tet_adjacent &fa,
    const jtf::mesh::one_ring_tet_at_edge & ortae,
    const jtf::mesh::one_ring_face_at_edge & orfae,
    const jtf::mesh::edge2cell_adjacent & ea, // surface edge
    boost::unordered_set<size_t> & g_unknown_face_in_original,
    boost::unordered_map<std::pair<size_t,size_t>,size_t> & inner_face_jump_type,
    const matrixst & outside_face,
    const matrixd & node);

bool is_near_miss_chain(
    const std::deque<std::pair<size_t,size_t> > & chain_orig,
    const boost::unordered_map<std::pair<size_t,size_t>, std::vector<
    std::tuple<size_t,size_t,size_t> > > &orig_edges2cut_edges,
    const boost::unordered_map<size_t,size_t> & surface_cut_type,
    const jtf::mesh::one_ring_face_at_point &orfap,
    const matrixst & outside_face_cut_idx,
    const matrixst & outside_face);

int iterately_remove_surface_degenerated_patch(
    matrixst &tet, matrixd &node,
    boost::unordered_map<std::pair<size_t,size_t>,size_t> & inner_face_jump_type,
    boost::unordered_map<size_t,size_t> & surface_type);
#endif // TOPOLOGY_ANALYSIS_H
