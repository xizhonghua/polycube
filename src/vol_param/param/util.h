#ifndef UTIL_H
#define UTIL_H

#include <boost/unordered_map.hpp>
#include <zjucad/matrix/matrix.h>
#include <set>
#include <jtflib/mesh/mesh.h>

///
/// @brief remove_degenerated_face_of_tet, and update inner_type
/// @param orig_face_in_cut
/// @param tet
/// @param node
/// @param param_node
/// @param inner_type
/// @return
///
int remove_degenerated_face_of_tet(
    const zjucad::matrix::matrix<size_t> & orig_face_in_cut,
    const jtf::mesh::face2tet_adjacent & fa_uncut,
    zjucad::matrix::matrix<size_t> & uncut_tet,
    zjucad::matrix::matrix<double> & uncut_node,
    zjucad::matrix::matrix<size_t> & cut_tet,
    zjucad::matrix::matrix<double> & cut_node,
    zjucad::matrix::matrix<double> & param_node,
    boost::unordered_map<std::pair<size_t,size_t>,size_t> & inner_type,
    zjucad::matrix::matrix<zjucad::matrix::matrix<double> > & frame,
    const jtf::mesh::edge2cell_adjacent & ea_cut,
    const std::string remesh_strategy);


///
/// @brief collapse_small_angle_edges, defualt is 5 degree
/// @param edges_in_cut
/// @param tet
/// @param node
/// @param param_node
/// @param inner_type
/// @param eqn_idx
/// @param eqn_coeff
/// @return
///
int collapse_small_angle_edges(
    const std::set<std::pair<size_t,size_t> > & edges_in_cut,
    const jtf::mesh::face2tet_adjacent & fa_uncut,
    zjucad::matrix::matrix<size_t> & uncut_tet,
    zjucad::matrix::matrix<double> & uncut_node,
    zjucad::matrix::matrix<size_t> & cut_tet,
    zjucad::matrix::matrix<double> & cut_node,
    zjucad::matrix::matrix<double> & param_node,
    boost::unordered_map<std::pair<size_t,size_t>,size_t> & inner_type,
    zjucad::matrix::matrix< zjucad::matrix::matrix<double> > & frame,
    std::vector<std::vector<size_t> > * eqn_idx = 0,
    std::vector<std::vector<double> > * eqn_coeff = 0);

///
/// @brief split_large_angle_edges, default angle is 170
/// @param edges_in_cut
/// @param fa_uncut
/// @param uncut_tet
/// @param cut_tet
/// @param cut_node
/// @param param_node
/// @param inner_type
/// @param eqn_idx
/// @param eqn_coeff
/// @return
///
int split_large_angle_edges(
    const std::set<std::pair<size_t,size_t> > & edges_in_cut,
    const jtf::mesh::face2tet_adjacent & fa_uncut,
    zjucad::matrix::matrix<size_t> & uncut_tet,
    zjucad::matrix::matrix<double> & uncut_node,
    zjucad::matrix::matrix<size_t> & cut_tet,
    zjucad::matrix::matrix<double> & cut_node,
    zjucad::matrix::matrix<double> & param_node,
    boost::unordered_map<std::pair<size_t,size_t>,size_t> & inner_type,
    std::vector<std::vector<size_t> > * eqn_idx = 0,
    std::vector<std::vector<double> > * eqn_coeff = 0);

///
/// \brief fix_sxx_tet_mapping, in library tet_mesh_sxx, each tet is reorderd, thus
///        connection between two tets should be rebuilt.
/// \param uncut_tet
/// \param cut_tet
/// \param cut_tet2tet
/// \return
///
int fix_sxx_tet_mapping(const zjucad::matrix::matrix<size_t> & uncut_tet,
                        zjucad::matrix::matrix<size_t> & cut_tet,
                        const std::vector<size_t> & cut_tet2tet,
                        const zjucad::matrix::matrix<size_t> & new_tet2orig_tet_uncut,
                        const zjucad::matrix::matrix<size_t> & new_tet2orig_tet_cut);

///
/// \brief fix_sxx_tet_mapping, in library tet_mesh_sxx, each tet is reorderd, thus
///        connection between two tets should be rebuilt.
/// \param uncut_tet
/// \param cut_tet
/// \param cut_tet2tet
/// \return
///
int fix_sxx_tet_mapping(const zjucad::matrix::matrix<size_t> & uncut_tet,
                        zjucad::matrix::matrix<size_t> & cut_tet,
                        const zjucad::matrix::matrix<size_t> & cut_tet2tet,
                        const zjucad::matrix::matrix<size_t> & new_tet2orig_tet_uncut,
                        const zjucad::matrix::matrix<size_t> & new_tet2orig_tet_cut);


int load_integer_groups(const char * filename,
                        std::vector<std::vector<size_t> > & integer_v_group);

int load_restricted_path(
    const char * filename,
    std::vector<std::vector<std::pair<size_t,size_t> > > & restricted_edges);


double get_minimal_distance(const std::vector<std::vector<size_t> > & integer_variants,
                            const zjucad::matrix::matrix<double> & node,
                            const double threshold = 0);
#endif // UTIL_H
