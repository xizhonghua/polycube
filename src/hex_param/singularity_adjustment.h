#ifndef SINGULARITY_ADJUSTMENT_H
#define SINGULARITY_ADJUSTMENT_H
#include "../tetmesh/tetmesh.h"
#include "hex_param.h"
#include "topology_operation.h"
#include <jtflib/util/vertex_connection.h>
#include "find_singularities.h"
#include <zjucad/matrix/io.h>
#include <hjlib/function/func_aux.h>
#include <zjucad/matrix/itr_matrix.h>

#include <vector>

/**
 * @brief
 *
 * @param tet
 * @param node
 * @param frame
 * @param ortae
 * @param std::vector<std::deque<std::pair<size_t
 * @param chain_list
 * @param singularities_type_
 * @param std::map<std::pair<size_t
 * @param size_t>
 * @param inner_face_jump_type
 * @param fa
 * @param pt
 * @param tet_rot_to_root
 * @return int
 */
int relabel_singularity_chain_by_splitting(
    jtf::tet_mesh &tm,
    zjucad::matrix::matrix<matrixd > &frame,
    std::vector<std::deque<std::pair<size_t,size_t> > > &chain_list,
    std::vector<std::deque<size_t> > &singularities_type_,
    boost::unordered_map<std::pair<size_t,size_t>,size_t>  &inner_face_jump_type,
    boost::property_tree::ptree &pt,
    std::vector<size_t> &tet_rot_to_root);

/**
 * @brief
 *
 * @param tet
 * @param node
 * @param outside_face
 * @param fa
 * @param ortae
 * @param std::vector<std::deque<std::pair<size_t
 * @param singularity_edges
 * @param singularity_type
 * @param std::map<std::pair<size_t
 * @param size_t>
 * @param inner_face_jump_type
 * @return int
 */
int relabel_singularity_chain_by_modify_face_jump_type(
    matrixst &tet,
    matrixd &node,
    matrixst &outside_face,
    matrixst &outside_face_idx,
    zjucad::matrix::matrix<matrixd > & frame,
    jtf::mesh::face2tet_adjacent &fa,
    jtf::mesh::one_ring_tet_at_edge &ortae,
    std::vector<std::deque<std::pair<size_t,size_t> > > &singularity_edges,
    std::vector<std::deque<std::size_t> > &singularity_type,
    boost::unordered_map<std::pair<size_t,size_t>,size_t> &inner_face_jump_type);


/**
 * @brief This function is used to decrease singularities number, detect each
 *        face. If there is a face which contains
 *
 * @param ortae   input one_ring_tet_at_edge
 * @param fa      inputjtf::mesh::face2tet_adjacent
 * @param singularity_edges input singularity_edges, and will be modified
 * @param singularity_type  input singularity_type, also will be modified
 * @param inner_face_jump_type  input face_jump_type which will be modified too
 * @return int
 */
int relabel_to_converge_singularities(
    const jtf::mesh::one_ring_tet_at_edge & ortae,
    const jtf::mesh::face2tet_adjacent& fa,
    std::vector<std::deque<std::pair<size_t,size_t> > > &singularity_edges,
    std::vector<std::deque<size_t> > &singularity_type,
    std::map<std::pair<size_t,size_t>, size_t> &inner_face_jump_type);


/**
 * @brief
 *
 * @param tet
 * @param node
 * @param fa
 * @param outside_face
 * @param ortae
 * @param std::vector<std::deque<std::pair<size_t
 * @param singularity_edges
 * @param singularity_type
 * @param std::map<std::pair<size_t
 * @param size_t>
 * @param inner_face_jump_type
 * @return int
 */
int relabel_remove_black_lines_by_modify_face_jump_type(
    const matrixst &tet,
    const matrixd &node,
    const jtf::mesh::face2tet_adjacent &fa,
    const matrixst &outside_face,
    const jtf::mesh::one_ring_tet_at_edge &ortae,
    std::vector<std::deque<std::pair<size_t,size_t> > > &singularity_edges,
    std::vector<std::deque<size_t> > &singularity_type,
    std::map<std::pair<size_t,size_t>,size_t> &inner_face_jump_type);

int relabel_remove_black_lines_by_modify_face_jump_type_and_splitting_fans(
    std::vector<size_t> &tet_array,
    std::vector<double> &node_array,
    std::vector<matrixd > & frame_array,
    jtf::mesh::face2tet_adjacent &fa,
    const matrixst &outside_face,
    jtf::mesh::one_ring_tet_at_edge &ortae,
    std::vector<std::deque<std::pair<size_t,size_t> > > &singularity_edges,
    std::vector<std::deque<size_t> > &singularity_type,
    boost::unordered_map<std::pair<size_t,size_t>,size_t> &inner_face_jump_type);

/**
 * @brief
 *
 * @param tet
 * @param node
 * @param outside_face
 * @param fa
 * @param ortae
 * @param std::vector<std::deque<std::pair<size_t
 * @param singularity_edges
 * @param std::map<std::pair<size_t
 * @param size_t>
 * @param inner_face_jump_type
 * @return int
 */
int relabel_zigzag_by_modify_face_jump_type(
    const matrixst &tet,
    const matrixd &node,
    const matrixst &outside_face,
    const jtf::mesh::face2tet_adjacent &fa,
    const jtf::mesh::one_ring_tet_at_edge &ortae,
    std::vector<std::deque<std::pair<size_t,size_t> > > &singularity_edges,
    std::vector<std::deque<size_t> > &singularity_type,
    boost::unordered_map<std::pair<size_t,size_t>,size_t> &inner_face_jump_type);


/**
 * @brief
 *
 * @param tet
 * @param node
 * @param outside_face
 * @param fa
 * @param ortae
 * @param std::vector<std::deque<std::pair<size_t
 * @param singularity_edges
 * @param std::map<std::pair<size_t
 * @param size_t>
 * @param inner_face_jump_type
 * @return int
 */
int relabel_remove_near_surface_by_modify_face_jump_type(
    const matrixst &tet,
    const matrixd &node,
    const matrixst &outside_face,
    const jtf::mesh::face2tet_adjacent &fa,
    const jtf::mesh::one_ring_tet_at_edge &ortae,
    std::vector<std::deque<std::pair<size_t,size_t > > > &singularity_edges,
    boost::unordered_map<std::pair<size_t,size_t>,size_t> & inner_face_jump_type);


/**
 * @brief
 *
 * @param fa
 * @param black_chain
 * @param other_chain
 * @param ortae
 * @param std::map<std::pair<size_t
 * @param size_t>
 * @param inner_face_jump_type
 * @param remove_faces
 * @return int
 */
int iterative_modify_face_type(
    const jtf::mesh::face2tet_adjacent &fa,
    const std::deque<size_t> & black_chain,
    const std::deque<size_t> &other_chain,
    const jtf::mesh::one_ring_tet_at_edge & ortae,
    boost::unordered_map<std::pair<size_t,size_t>,size_t> &inner_face_jump_type,
    std::vector<size_t> &remove_faces);

int split_fans_between_two_path(
    const std::deque<size_t> &black_chain,
    std::deque<size_t> &path,
    std::vector<size_t> &tet_array,
    std::vector<double> &node_array,
    std::vector<matrixd > & frame,
    jtf::mesh::face2tet_adjacent &fa,
    jtf::mesh::one_ring_tet_at_edge &ortae,
    boost::unordered_map<std::pair<size_t,size_t>,size_t> &inner_face_jump_type);

// return -1 if is not zigzag, else return the common face idx

/**
 * @brief
 *
 * @param std::deque<std::pair<size_t
 * @param chain
 * @param chain_type
 * @return int
 */
int reverse_singularity_with_type(
    std::deque<std::pair<size_t,size_t> > & chain,
    std::deque<size_t> & chain_type);

/**
 * @brief
 *
 * @param std::deque<std::pair<size_t
 * @param chain
 * @return int
 */
int reverse_singularity(std::deque<std::pair<size_t,size_t> > & chain);

/**
 * @brief
 *
 * @param std::pair<size_t
 * @param edge0
 * @param std::pair<size_t
 * @param edge1
 * @param face_idx
 * @param ortae
 * @param fa
 * @param std::map<std::pair<size_t
 * @param size_t>
 * @param inner_face_jump_type
 * @param tet
 * @param node
 * @return int
 */
int remove_zigzag(
    const std::pair<size_t,size_t> &edge0,
    const std::pair<size_t,size_t> &edge1,
    const size_t face_idx,
    const jtf::mesh::one_ring_tet_at_edge & ortae,
    const jtf::mesh::face2tet_adjacent &fa,
    boost::unordered_map<std::pair<size_t,size_t>,size_t> &inner_face_jump_type,
    const matrixst &tet,
    const matrixd &node);

//// used new singularity definition under global reference
//int remove_zigzag_new(const std::pair<size_t,size_t> & edge0,
//                      const std::pair<size_t,size_t> & edge1,
//                      const size_t face_idx,
//                      const one_ring_tet_at_edge &ortae,
//                      const jtf::mesh::face2tet_adjacent & fa,
//                      std::map<std::pair<size_t,size_t>,size_t> &inner_face_jump_type,
//                      const matrixst &tet);

//int relabel_face_type_to_remove_zigzag(
//    std::vector<std::deque<std::pair<size_t,size_t> > > &chain_list,
//    std::vector<std::deque<size_t> >  &singularities_type_,
//    const matrixst &tet,
//    const jtf::mesh::face2tet_adjacent & fa,
//    const matrixst &outside_face,
//    std::map<std::pair<size_t,size_t>,size_t> &inner_face_jump_type,
//    matrixst &tet_rot_type,
//    const one_ring_tet_at_edge &ortae);

/**
 * @brief
 *
 * @param tet
 * @param node
 * @param outside_face
 * @param fa
 * @param ortae
 * @param std::deque<std::pair<size_t
 * @param chain
 * @param std::map<std::pair<size_t
 * @param size_t>
 * @param inner_face_jump_type
 * @param vc
 * @return int
 */
int remove_near_surface_chain(
    const matrixst &tet,
    const matrixd &node,
    const matrixst &outside_face,
    const jtf::mesh::face2tet_adjacent &fa,
    const jtf::mesh::one_ring_tet_at_edge &ortae,
    const std::deque<std::pair<size_t,size_t> > &chain,
    boost::unordered_map<std::pair<size_t,size_t>,size_t> &inner_face_jump_type,
    const vertex_connection<UNDIRECT> & vc);


/**
 * @brief
 *
 * @param tet_array
 * @param node_array
 * @param tet_idx
 * @param fa
 * @param ortae
 * @param frame
 * @param std::pair<size_t
 * @param black_edge
 * @param angle_zyx
 * @param std::map<std::pair<size_t
 * @param size_t>
 * @param inner_face_jump_type
 * @param std::map<std::pair<size_t
 * @param size_t>
 * @param std::vector<std::pair<size_t
 * @param split_edge_map
 * @param black_line_type
 * @return int
 */
int split_to_remove_black_line(
    std::vector<size_t> &tet_array,
    std::vector<double> &node_array,
    const size_t &tet_idx,
    const jtf::mesh::face2tet_adjacent &fa,
    jtf::mesh::one_ring_tet_at_edge & ortae,
    std::vector<matrixd > &frame,
    const std::pair<size_t,size_t> &black_edge,
    const matrixd &angle_zyx,
    boost::unordered_map<std::pair<size_t,size_t>,size_t> &inner_face_jump_type,
    boost::unordered_map<std::pair<size_t,size_t>,std::vector<std::pair<size_t,size_t> > > &split_edge_map,
    size_t black_line_type); // 1:compound of two axes, 0 means compound of three axes

//! @brief use splitting strategy to remove zigzag lines
/**
 * @brief
 *
 * @param tet_array
 * @param node_array
 * @param std::vector<std::deque<std::pair<size_t
 * @param chain_list
 * @param ortae
 * @param std::map<std::pair<size_t
 * @param size_t>
 * @param inner_face_jump_type
 * @param frame
 * @return int
 */
int split_to_remove_zigzag(
    std::vector<size_t> &tet_array,
    std::vector<double> &node_array,
    const std::vector<std::deque<std::pair<size_t,size_t> > > &chain_list,
    jtf::mesh::one_ring_tet_at_edge &ortae,
    boost::unordered_map<std::pair<size_t,size_t>,size_t > &inner_face_jump_type,
    std::vector<matrixd > &frame);

/**
 * @brief
 *
 * @param tet
 * @param node
 * @param outside_face
 * @param ea
 * @param fa
 * @param vc
 * @param ortae
 * @param std::deque<std::pair<size_t
 * @param chain
 * @param std::map<std::pair<size_t
 * @param size_t>
 * @param inner_face_jump_type
 * @return int
 */
int modify_face_type_to_remove_near_surface(
    const matrixst &tet,
    const matrixd &node,
    const matrixst outside_face,
    const jtf::mesh::edge2cell_adjacent &ea,
    const jtf::mesh::face2tet_adjacent &fa,
    const vertex_connection<UNDIRECT> &vc,
    const jtf::mesh::one_ring_tet_at_edge & ortae,
    const std::deque<std::pair<size_t,size_t> > &chain,
    boost::unordered_map<std::pair<size_t,size_t>,size_t> &inner_face_jump_type);

//! @brief split edge <first,second> at (1-lambda) * first + lambda * second
// WARNINGl: if you want to iteratorly subdivide the tets around edge, you should go along first to second
// because I used fa to get faces around point "edge_second" and I didn't update fa;
// only if you go along first--> second, can you get correct answer.
// after this function called, the node_array will add one point.
int split_tets_around_edge(
    std::vector<size_t> &tet_array,
    std::vector<double> &node_array,
    jtf::mesh::one_ring_tet_at_edge &ortae,
    std::vector<matrixd > &frame,
    const std::pair<size_t,size_t> &edge,
    const double lambda,
    boost::unordered_map<std::pair<size_t,size_t>,std::vector<std::pair<size_t,size_t> > > &split_edge_map,
    boost::unordered_map<size_t,size_t > & old_tet_map_to_new,
    boost::unordered_map<std::pair<size_t,size_t>,size_t> &inner_face_jump_type);

int apply_global_U_rotation(
    const matrixst &tet,
    const matrixd &node,
    std::map<size_t,size_t> &face_jump,
    const zjucad::matrix::matrix<matrixd > &frame_inner,
    const jtf::mesh::face2tet_adjacent & fa,
    const matrixst &outside_face_idx);

/**
 * @brief This function is used to extract surface normal align type according to
 *        frame,
 * TODO: This name is not much suitable, need to replace it
 *
 * @param tet         input tet mesh
 * @param node        input tet node
 * @param face_jump   output surface type
 * @param frame_inner input frame
 * @param fa          input face to tet adjacent relations of tet mesh
 * @param outside_face_idx  input outside face index
 * @return int        if work fine, return 0, or return non-zeros
 */
int apply_global_U_rotation_with_normal(
    const matrixst &tet,
    const matrixd &node,
    boost::unordered_map<size_t,size_t> &face_jump,
    const zjucad::matrix::matrix<matrixd > &frame_inner,
    const jtf::mesh::face2tet_adjacent & fa,
    const matrixst &outside_face_idx,
    const matrixd & outside_face_normal);
#endif // SINGULARITY_ADJUSTMENT_H
