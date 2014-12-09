#ifndef HEX_PARAM_IO_H
#define HEX_PARAM_IO_H

#include <boost/unordered_set.hpp>
#include <boost/unordered_map.hpp>
#include "../common/def.h"
#include "../tetmesh/tetmesh.h"
class singularity_graph;

int dump_out_param_config(
    const char * group_file,
    const char * equation_file,
    const char * chain_file,
    const char * cut_tet_file,
    const matrixst & cut_tet,
    const matrixd & cut_node,
    const matrixst & cut_tet2tet,
    const std::vector<size_t> & rot_type,
    const jtf::mesh::face2tet_adjacent & fa_cut,
    const singularity_graph & sg,
    const boost::unordered_set<size_t> & restricted_nodes,
    const boost::unordered_map<std::pair<size_t,size_t>,size_t> & jump_face2rot_idx,
    const std::vector<std::pair<size_t,size_t> > &g_unknown_face_pair);

//! @brief: this function is used to load fnode_with_group info from file
//          notice that: fnode size equals 3 * vertex number in cut tet
//  @param filename   input file
//  @param fnode      output fnode with group info
int load_fnode_group(
    const char * filename,
    const matrixst & cut_tet,
    std::vector<size_t> & fnode);

int dump_surface_restricted_type_to_vtk(
    const char * filename,
    const std::string & type_name,
    const matrixd & node,
    const jtf::mesh::face2tet_adjacent & fa,
    const boost::unordered_map<size_t,size_t> & surface_type);

/**
 * @brief
 *
 * @param filname
 * @param std::map<std::pair<size_t
 * @param size_t>
 * @param inner_face_jump_type
 * @return int
 */
int load_inner_face_jump_type(
    const char * filname,
    boost::unordered_map<std::pair<size_t,size_t>,size_t> & inner_face_jump_type);

int load_g_unknown_face(
    const char * filename,
    std::vector<std::pair<size_t,size_t> > & face_pairs);

int load_loop_points(
    const char * filename,
    boost::unordered_set<size_t> & loop_points);

int load_restricted_edges(
    const char * filename,
    boost::unordered_set<std::pair<size_t,size_t> > &degenerated_edges);


///
/// @brief load_surface_type,
///        WARNING!!! since three surface type files are used:
///        (p0,p1,p2, normal_align rotation matrix type)
///        (face_idx, restricted type)
///        (face_idx, normal_align rotation matrix type)
///        Need to unify them
/// @param filename
/// @param surface_type
/// @param fa
/// @return 0 if read restricted type, 1 if read normal_align rotation matrix type,
///         other-non-zeros if error
///
int load_surface_type(
    const char * filename,
    boost::unordered_map<size_t,size_t> & surface_type,
    const jtf::mesh::face2tet_adjacent *fa = 0);

///
/// \brief input surface normal align type, output surface restricted type
///
/// \param surface_type
/// \return
///
int convert_surface_normal_type2_restricted_type(
    boost::unordered_map<size_t,size_t> & surface_type);
/**
 * @brief
 *
 * @param filename
 * @param fa
 * @param outside_face_idx
 * @param outside_face_type
 * @return int
 */
int load_surface_type(
    const char * filename,
    const jtf::mesh::face2tet_adjacent &fa,
    const matrixst & outside_face_idx,
    matrixst & outside_face_type);

//int load_surface_normal_align_type(
//    const char * filename,
//    const jtf::mesh::face2tet_adjacent &fa,
//    boost::unordered_map<size_t,size_t> &outside_face_type);


///**
// * @brief load surface restricted_type,
// *
// * @param filename  input surface type name
// * @param surface_type  output surface type
// * @return return 0, if load restricted type, return 1 if normal rotation type, return others if meet error
// */
int load_surface_restricted_type(
    const char * filename,
    boost::unordered_map<size_t,size_t> & surface_type);

int dump_surface_restricted_type(
    const char * filename,
    const boost::unordered_map<size_t,size_t> & surface_type);


/**
 * @brief
 *
 * @param filename
 * @param outside_face
 * @param outside_face_type
 * @return int
 */
int dump_surface_normal_align_type(
    const char * filename,
    const matrixst & outside_face,
    const matrixst & outside_face_type);

int dump_surface_normal_align_type_map(
    const char * filename,
    const matrixst & outside_face,
    const matrixst & outside_face_idx,
    const boost::unordered_map<size_t,size_t> & face_jump_type);

int load_restricted_edge(
    const char * filename,
    boost::unordered_map<std::pair<size_t,size_t>,size_t> & restrict_edge_type);

int dump_restricted_edge(
    const char *filename,
    const std::vector<std::pair<size_t,size_t> > & restriced_edges,
    const std::vector<size_t> & restriced_edges_type);

int load_chain_file(const char * chain_file,
                    std::vector<std::vector<size_t> > & chains);

///
/// @brief load_group_file, load variants group file
/// @param group_file
/// @param groups
/// @param integer_group_flag
/// @return
///
int load_group_file(const char * group_file,
		    std::vector<std::vector<size_t> > & groups,
		    std::vector<bool> * integer_group_flag = 0);

int dump_surface_normal_align_type2vtk(
    const char *filename,
    const zjucad::matrix::matrix<double> & node,
    const jtf::mesh::face2tet_adjacent & fa,
    const boost::unordered_map<size_t,size_t> &outside_face_type);

int dump_inner_face_jump_type(
    const char * filename,
    const boost::unordered_map<std::pair<size_t,size_t>,size_t> & inner_face_jump_type);


int dump_inner_face_jump_type2vtk(
    const char *filename,
    const zjucad::matrix::matrix<double> & node,
    const zjucad::matrix::matrix<size_t> & tet,
    const boost::unordered_map<std::pair<size_t,size_t>,size_t> &inner_type);

#endif // HEX_PARAM_IO_H
