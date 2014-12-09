#ifndef REMOVE_SURFACE_WEDGE_H
#define REMOVE_SURFACE_WEDGE_H

#include "../common/def.h"

#include <jtflib/mesh/mesh.h>
#include "topology_analysis.h"
#include <boost/unordered_set.hpp>
#include <boost/unordered_map.hpp>


/**
 * @brief   this function is used to remove surface wedge, surface degeneration
 *          consists of many complex cases, currently I can only handle simple
 *          case: one patch which is inside another patch
 *
 * @param tet  input original tet mesh
 * @param node  input node
 * @param cut_tet  input cut_tet
 * @param cut_node input cut_node
 * @param inner_face_jump_type input inner_face_jump_type
 * @param surface_type the initial surface typ, and will be modified
 * @param loop_points_in_cut all degenerated points detected by graph
 * @return int  return 0 if works fine, or return non-zeros
 */
int remove_surface_wedge(
    const matrixst & tet,
    const matrixd &node,
    const matrixst &cut_tet,
    const matrixd & cut_node,
    const std::vector<std::pair<size_t,size_t> > &g_unknown_face_pair,
    const boost::unordered_map<std::pair<size_t,size_t>,size_t > & restricted_edge_from_graph,
    const bool use_uvw_restricted_surface_type,
    boost::unordered_map<std::pair<size_t,size_t>,size_t> & inner_face_jump_type,
    matrixst & new_tet,
    matrixd & new_node,
    boost::unordered_map<size_t,size_t> & surface_type_cut,
    matrixd * zyz_ptr = 0);

int remove_face_degeneration_by_splitting(
   jtf::mesh::face2tet_adjacent &fa,
    boost::unordered_map<size_t,size_t> &surface_type_original,
    zjucad::matrix::matrix<double> &new_node,
    zjucad::matrix::matrix<size_t> &new_tet,
    boost::unordered_map<std::pair<size_t,size_t>,size_t> &inner_face_jump_type);

///**
// *  NOT USED ANY MORE !!!!!!!!!
// * @brief   this function is used to remove each degenerated patch, in this
// *          function, it will only handle those local degenerated ones, here I
// *          defined a local one as: a patch which is detected by checking its
// *          arounding degenerated points, if all boundary points of this patch
// *          is degenerated and all boundary points are not critical, flip the type
// *
// * @param patches input patches which will be degenerated
// * @param orfap   input one_ring_face_at_point for surface points
// * @param outside_face_idx input outside_face_idx
// * @param boundary_points boundary points for each patch
// * @param surface type inital surface type which will be modified.
// * @return int  return 0 if works fine, or return non-zeros
// */
//int remove_surface_degenerated_patch_first_filter(
//    const std::vector<boost::unordered_set<size_t> > & patchs,
//    const jtf::mesh::one_ring_face_at_point & orfap,
//    const matrixst & outside_face_idx,
//    const std::vector<boost::unordered_set<size_t> > & boundary_points,
//    boost::unordered_map<size_t,size_t> & surface_type);


/**
 * @brief   this function is used to detect whether a surface points is critical
 *          points, which means the arounding surface type of such point is fractional
 *          For example:
 *           1.  0 0 1 1 2 2 is not critical
 *           2.  0 1 0 1 is critical
 *
 * @param arounding_faces input arounding faces, !!!!!! MUST BE ORDERED !!!!!
 * @param outside_face_idx   input outside face idx
 * @param surface type input surface type, 0/1/2
 * @return bool  true for critical point, false for not
 */
bool is_critical_point(
    const std::vector<size_t> & arounding_faces,
    const matrixst & outside_face_idx,
    const boost::unordered_map<size_t,size_t> & surface_type);

int remove_isolated_patch(
    const matrixst & tet,
    const matrixd &node,
    const jtf::mesh::face2tet_adjacent & fa,
    const std::vector<boost::unordered_set<size_t> > patches,
    const boost::unordered_map<size_t, boost::unordered_set<size_t> > &patch_links_info,
    boost::unordered_map<size_t,size_t> & surface_type_original);

int remove_fractional_degree_two_patches(
    const matrixst & tet,
    const matrixd &node,
    const jtf::mesh::face2tet_adjacent & fa,
    const std::vector<boost::unordered_set<size_t> > patches,
    const boost::unordered_map<size_t, boost::unordered_set<size_t> > &patch_links_info,
    boost::unordered_map<size_t,size_t> & surface_type_original,
    const boost::unordered_map<std::pair<size_t,size_t>, double>  & patch2patch_len,
    const size_t fractional_cell_num = -1);

///////////////////////////////////////////////////////////////////////

int remove_multi_orientation_group(
    const matrixd & node,
    const jtf::mesh::edge2cell_adjacent & ea,
    const jtf::mesh::face2tet_adjacent & fa,
    const matrixst & outside_face,
    const matrixst & outside_face_idx,
    const boost::unordered_map<size_t, matrixd> & surface_normal,
    const boost::unordered_map<size_t,size_t> &surface_type_orient_orig,
    const std::vector<boost::unordered_set<size_t> > &patches,
    const boost::unordered_map<size_t, boost::unordered_set<size_t> > &patch_links_info,
    const size_t N_ring_to_modify,
    boost::unordered_map<size_t,size_t> &surface_type_original);

int modify_surface_type_to_remove_separated_edges(
    const std::deque<std::pair<size_t,size_t> > & one_chain,
    const boost::unordered_set<size_t> &one_patch,
    const matrixst & outside_face_idx,
    const boost::unordered_map<size_t, matrixd> & surface_normal,
    boost::unordered_map<size_t,size_t> & surface_type_original,
    const jtf::mesh::edge2cell_adjacent &ea,
    const jtf::mesh::N_ring_face_at_point &nrfap,
    const size_t N_ring_to_modify);
/**
 * @brief   this function is used to extrace surface patch, it uses surface
 *          restricted type and cut edge to depart surface patch
 *
 * @param tet input tet mesh
 * @param cut_tet input cut tetmesh
 * @param node input node of tet mesh
 * @param cut_node input node of cut tet mesh
 * @param cut_tet2tet input cut_tet2tet relations
 * @param outside_face input outside
 * @param outside_face_idx input outside face idx mapping
 * @param fa  inputjtf::mesh::face2tet_adjacent of original tetmesh
 * @param fa_cut inputjtf::mesh::face2tet_adjacent of cut tetmesh
 * @param ea  input edge2face adjacent of outside face
 * @param g_unknown_face_pair  input g_unknown_faces, also the cutting faces
 * @param surface_type_original input surface type
 * @param patches output patches
 * @param group_linking_info output group linking info
 * @return 0 if work ok, or return non-zeros
 */
int extract_surface_patch_graph(
    const matrixst & tet,
    const matrixst & cut_tet,
    const matrixd & node,
    const matrixst & cut_tet2tet,
    const matrixst & outside_face,
    const matrixst & outside_face_idx,
    const jtf::mesh::face2tet_adjacent & fa,
    const jtf::mesh::face2tet_adjacent & fa_cut,
    const jtf::mesh::edge2cell_adjacent & ea,
    const std::vector<std::pair<size_t,size_t> > & g_unknown_face_pair,
    const boost::unordered_map<size_t,size_t> & surface_type_original,
    std::vector<boost::unordered_set<size_t> > & patches,
    boost::unordered_map<size_t,boost::unordered_set<size_t> > & group_linking_info,
    boost::unordered_map<std::pair<size_t,size_t>,double> *patch_to_patch_boundary_len = 0);

int remove_fractional_patches_by_degenerated_edges(
    const matrixst & tet,
    const matrixd & node,
    const jtf::mesh::face2tet_adjacent &fa,
    const jtf::mesh::edge2cell_adjacent & ea,
    const boost::unordered_set<std::pair<size_t,size_t> > & degenerated_edges,
    boost::unordered_map<size_t,size_t> &surface_type_original);
////////////////////////////////////////////////////////////////////////

int gather_restricted_edges_on_polycube(
    const matrixst & tet,
    const jtf::mesh::edge2cell_adjacent & ea,
    const matrixst & outside_face_idx,
    const boost::unordered_map<size_t,size_t> & surface_type_original,
    boost::unordered_map<std::pair<size_t,size_t>,size_t > & restricted_edges);

int remove_one_face_degeneration(
    const matrixst &tet,
    const matrixst & outside_face,
    const matrixst & outside_face_idx,
    const boost::unordered_map<std::pair<size_t,size_t>,size_t > & restricted_edges,
    boost::unordered_map<size_t,size_t> &surface_type_original);

int remove_surface_zigzag_by_restricted_edges(
    const matrixst &tet,
    const matrixd & node,
    const matrixst & outside_face_idx,
    const jtf::mesh::edge2cell_adjacent & ea,
    const boost::unordered_map<std::pair<size_t,size_t>,size_t> &restricted_edges,
    const jtf::mesh::face2tet_adjacent &fa,
    matrixst &new_tet,
    matrixd &new_node,
    boost::unordered_map<size_t,size_t> & surface_type,
    boost::unordered_map<std::pair<size_t,size_t>,size_t> & inne_face_jump_type,
    matrixd * zyz_ptr = 0);

int remove_surface_zigzag_by_restricted_edges_from_graph(
    const matrixst &tet,
    const matrixd &node,
    const matrixst &outside_face_idx,
   jtf::mesh::face2tet_adjacent &fa,
    const jtf::mesh::edge2cell_adjacent &ea,
    const boost::unordered_map<std::pair<size_t,size_t>,size_t> &restricted_edge_from_graph_orig,
    boost::unordered_map<size_t,size_t> &surface_type_original,
    matrixst &ew_tet,
    matrixd &new_node,
    boost::unordered_map<std::pair<size_t,size_t>,size_t> & inner_face_jump_type,
    matrixd * zyz_ptr = 0);

//! @brief This function assume all local degenerations has been cleaned, all left
// is global degeneration which we should remove one whole patch or break through
// between two patches
int remove_degenerated_patch_and_break_through(
    const matrixst &tet,
    const matrixd &node,
    const matrixst &cut_tet,
    const matrixst &cut_tet2tet,
    const jtf::mesh::face2tet_adjacent & fa,
    const jtf::mesh::face2tet_adjacent & fa_cut,
    const jtf::mesh::edge2cell_adjacent & ea,
    const std::vector<std::pair<size_t,size_t> > & g_unknown_face_pair,
    const matrixst & outside_face,
    const matrixst & outside_face_idx,
    const matrixd & polycube_node,
    boost::unordered_map<size_t,size_t> &surface_type);

////////////////////////////////////////////////////////////////////////////

int remove_surface_critcial_points(
    const matrixst &tet,
    const matrixd &node,
    const matrixst &outside_face,
    const matrixst &outside_face_idx,
    const jtf::mesh::face2tet_adjacent &fa,
    const jtf::mesh::edge2cell_adjacent &ea,
    boost::unordered_map<size_t,size_t> & surface_type);

int extract_relax_surface(
    const matrixst & tet,
    const matrixd & node,
    const boost::unordered_map<size_t,size_t> &surface_type,
    const boost::unordered_set<size_t> &loop_points,
    boost::unordered_map<size_t,size_t> &restricted_surface_type);

//! @brief this function is used to remove compound edge on surface by collapsing
//         edges
// @return 0: remove compound edge succeed
//         1: there is no compound edge need to collapse
//         other: shit happens
int collapse_degenerated_edges(
    matrixst & tet,
    matrixd & node,
   jtf::mesh::face2tet_adjacent & fa,
    boost::unordered_map<std::pair<size_t,size_t>,size_t> & inner_face_jump_type,
    boost::unordered_map<size_t,size_t> & surface_type,
    boost::unordered_map<std::pair<size_t,size_t>,size_t> & restricted_edges);

int remove_dgenerated_patch_by_modify_type(
    const matrixst &tet,
    const matrixst & cut_tet,
    const matrixd &node,
    const jtf::mesh::face2tet_adjacent & fa,
    const matrixst & outside_face,
    const matrixst & outside_face_idx,
    const jtf::mesh::edge2cell_adjacent & ea,
    const std::vector<std::pair<size_t,size_t> > & g_unknown_face_pair,
    boost::unordered_map<size_t,size_t> &surface_type);

int straighten_patch_boundary(
    const matrixst &outside_face,
    const matrixst &outside_face_idx,
    const jtf::mesh::edge2cell_adjacent &ea,
    boost::unordered_map<size_t,size_t> &surface_type_original);

int straighten_patch_boundary(
    const matrixst &outside_face,
    const matrixst &outside_face_idx,
    const jtf::mesh::edge2cell_adjacent &ea,
    const jtf::mesh::face2tet_adjacent &fa,
    boost::unordered_map<size_t,size_t> &surface_type_original,
    const boost::unordered_map<std::pair<size_t,size_t>,size_t> & rot_type);

int update_surface_type(
    const jtf::mesh::face2tet_adjacent & fa,
    const jtf::mesh::face2tet_adjacent & fa_new,
    const boost::unordered_map<std::vector<size_t>,std::vector<size_t> > & f2omap,
    boost::unordered_map<size_t,size_t> & surface_type);

int update_inner_face_jump_type(
    const jtf::mesh::face2tet_adjacent &fa_new,
    const matrixst & tet_idx_map,
    boost::unordered_map<std::pair<size_t,size_t>,size_t> & inner_face_jump_type);

#endif
