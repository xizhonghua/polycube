#ifndef TOPOLOGY_OPERATION_H
#define TOPOLOGY_OPERATION_H


#include <memory>
#include <iostream>
#include <zjucad/matrix/io.h>
#include <zjucad/matrix/itr_matrix.h>
#include "../hex_param/find_singularities.h"
#include <jtflib/mesh/mesh.h>
#include <jtflib/mesh/util.h>


/**
 * @brief WARNING!!! Experimental dangeous! This function is used to remove one
 *        edge in original tet mesh, and then
 *        update all corresponding informations
 *
 * @param edge      input edge,
 * @param tet_vec   input tet mesh with vector,
 * @param node_vec  input node with vector
 * @param zyz_vec   input zyz with vector
 * @param ortae     input one_ring_tet_at_edge info
 * @return int      return 0 if works well, non-zeros for erros.
 */
int remove_edge_and_update_info(
    const std::pair<size_t,size_t> & edge,
    std::vector<size_t> & tet_vec,
    const std::vector<double> & node_vec,
    std::vector<matrixd > & zyz_vec,
    jtf::mesh::one_ring_tet_at_edge & ortae);


/**
 * @brief This function is used to remove one edge in original tet mesh,
 *
 * @param edge      input edge,
 * @param tet_vec   input tet mesh with vector,
 * @param node_vec  input node with vector
 * @param zyz_vec   input zyz with vector
 * @return int      return 0 if works well, non-zeros for erros.
 */
int remove_edge_in_tet_raw(
    const std::pair<size_t,size_t> & edge,
    std::vector<size_t> & tet_vec,
    const std::vector<double> & node_vec,
    std::vector<matrixd > & zyz_vec);
/**
 * @brief This function is used to remove a chain in original tet mesh, and then
 *        update all corresponding informations
 *
 * @param one_chain     input chain which need to be removed
 * @param tet_vec   input tet mesh with vector,
 * @param node_vec  input node with vector
 * @param zyz_vec   input zyz with vector
 * @param ortae     input one_ring_tet_at_edge info
 * @return int      return 0 if works well, non-zeros for erros.
 */
int remove_one_chain_and_update_info(
    const std::deque<std::pair<size_t,size_t> > & one_chain,
    std::vector<size_t> & tet_vec,
    const std::vector<double> & node_vec,
    std::vector<matrixd > & zyz_vec,
    jtf::mesh::one_ring_tet_at_edge & ortae);

/**
 * @brief
 *
 * @param std::pair<size_t
 * @param edge
 * @param std::pair<size_t
 * @param tet_pair
 * @param tet
 * @param ortae
 * @param std::map<std::pair<size_t
 * @param size_t>
 * @param inner_face_jump_type
 * @return int
 */
int modidy_edge_type_by_change_face_type(
    const std::pair<size_t,size_t> & edge,
    const std::pair<size_t,size_t> & tet_pair,
    const matrixst & tet,
    const jtf::mesh::one_ring_tet_at_edge & ortae,
    boost::unordered_map<std::pair<size_t,size_t>,size_t> & inner_face_jump_type);

/**
 * @brief
 *
 * @param ortae
 * @param std::map<std::pair<size_t
 * @param size_t>
 * @param inner_face_jump_type
 * @param singularity_edges
 * @return int
 */
int extract_singulairty_edges(
    const jtf::mesh::one_ring_tet_at_edge &ortae,
    const std::map<std::pair<size_t,size_t>,size_t> &inner_face_jump_type,
    std::vector<size_t> &singularity_edges);

/**
 * @brief
 *
 * @param tet_idx_a
 * @param tet_idx_b
 * @param std::map<std::pair<size_t
 * @param size_t>
 * @param inner_face_jump_type
 * @param edge0_tet
 * @return int
 */
int modify_face_jump_type_at_given_tets_edge(
    const size_t tet_idx_a,
    const size_t tet_idx_b,
   boost::unordered_map<std::pair<size_t,size_t>,size_t> &inner_face_jump_type,
    const std::vector<size_t> &edge0_tet);

/**
 * @brief
 *
 * @param tet_idx_a
 * @param tet_idx_b
 * @param std::map<std::pair<size_t
 * @param size_t>
 * @param inner_face_jump_type
 * @param edge0_tets
 * @return int
 */
int modify_face_jump_type_to_remove_singularity_edge(
    const size_t tet_idx_a,
    const size_t tet_idx_b,
    boost::unordered_map<std::pair<size_t,size_t>,size_t> &inner_face_jump_type,
    const std::vector<size_t> &edge0_tets);
#if 0
int find_shortest_path(const size_t &begin_,
                       const size_t &end_,
                       const matrixst & tet,
                       const matrixd &node,
                       const std::deque<size_t>  &black_chain,
                       std::map<std::pair<size_t,size_t>, bool > &map_changed_edge);
#endif
/**
 * @brief
 *
 * @param begin_
 * @param end_
 * @param tet
 * @param node
 * @param black_chain
 * @param path
 * @return int
 */
int find_shortest_path_new(
    const size_t begin_,
    const size_t end_,
    const matrixst &tet,
    const matrixd &node,
    const std::deque<size_t> &black_chain,
    std::deque<size_t> &path);

int find_triangle_fan_path(
    const size_t begin_point_of_path,
    const matrixst &tet,
    const matrixd &node,
    const std::deque<size_t> &black_chain,
    const  jtf::mesh::one_ring_tet_at_edge & ortae,
    const jtf::mesh::face2tet_adjacent & fa,
    std::deque<size_t> &path);
/**
 * @brief
 *
 * @param std::pair<size_t
 * @param edge0
 * @param std::pair<size_t
 * @param edge1
 * @param fa
 * @return size_t
 */
size_t is_zigzag(const std::pair<size_t,size_t> &edge0,
                 const std::pair<size_t,size_t> &edge1,
                 const jtf::mesh::face2tet_adjacent &fa);

//! @brief detect whether two edges are zigzag edges
/**
 * @brief
 *
 * @param std::pair<size_t
 * @param edge0
 * @param std::pair<size_t
 * @param edge1
 * @param tet_array
 * @param tet_num
 * @return bool
 */
bool is_zigzag_edge(const std::pair<size_t,size_t> &edge0,
                    const std::pair<size_t,size_t> &edge1,
                    const size_t * tet_array,
                    size_t tet_num);

/**
 * @brief
 *
 * @param tet
 * @param std::pair<size_t
 * @param edge
 * @param ortae
 * @param one_ring_vertex
 * @return int
 */
int find_one_ring_vertex_around_edge(
    const matrixst &tet,
    const std::pair<size_t,size_t>& edge,
    const jtf::mesh::one_ring_tet_at_edge & ortae,
    std::vector<size_t> &one_ring_vertex);

//! @brief find one ring edges around edge, the order of edges is in right hand order
/**
 * @brief
 *
 * @param tet
 * @param node
 * @param std::pair<size_t
 * @param around_edge
 * @param ortae
 * @param std::vector<std::pair<size_t
 * @param one_ring_edges
 * @return int
 */
int find_one_ring_edges_around_edge(
    const matrixst &tet,
    const matrixd &node,
    const std::pair<size_t,size_t>& around_edge,
    const jtf::mesh::one_ring_tet_at_edge & ortae,
    std::vector<std::pair<size_t,size_t> > &one_ring_edges);

//! @brief find adjacent tets of face <a,b,c>,
// if face is an outside face tets will be <-1,t>
// if face is not valid, tets will be <-1,-1>
/**
 * @brief
 *
 * @param point_a
 * @param point_b
 * @param point_c
 * @param tet_array
 * @param tet_num
 * @return std::pair<size_t, size_t>
 */
std::pair<size_t,size_t> get_adj_tets(const size_t point_a,
                                      const size_t point_b,
                                      const size_t point_c,
                                      const size_t * tet_array,
                                      const size_t tet_num);

//! @brief check if "other_edge" is in right hand order of "around_edge"
// around_edge : first-->second
// other_edge: [1] --> [0]
/**
 * @brief
 *
 * @param tet
 * @param node
 * @param std::pair<size_t
 * @param around_edge
 * @param other_edge
 * @return int
 */
int adjust_cross_edges_right_hand_order_of_tet(
    const matrixst &tet,
    const matrixd &node,
    const std::pair<size_t,size_t> &around_edge,
    std::vector<size_t> &other_edge);


// this function only define the chain which has all vertex one ring near surface as near surface chain.
/**
 * @brief
 *
 * @param std::deque<std::pair<size_t
 * @param chain
 * @param tet
 * @param outside_face
 * @return bool
 */
bool is_one_ring_near_surface_chain(
    const std::deque<std::pair<size_t,size_t> > &chain,
    const matrixst &tet,
    const matrixst &outside_face);

/**
 * @brief
 *
 * @param vertex
 * @param tet
 * @param outside_face
 * @return bool
 */
bool is_one_ring_near_surface(
    const size_t vertex,
    const matrixst &tet,
    const matrixst &outside_face);

//! @brief this function is used to get surface patch idx,
// for each surface triganle, define an unique patch idx
/**
 * @brief
 *
 * @param surface_path_idx
 * @param surface_normal_align_type
 * @param tet
 * @param node
 * @param outside_face
 * @param outside_face_idx
 * @param fa_new
 * @return int
 */
int get_surface_patch_type(
    std::vector<size_t> &surface_path_idx,
    const std::vector<size_t> &surface_normal_align_type,
    const matrixst & tet,
    const matrixd & node,
    const matrixst & outside_face,
    const matrixst & outside_face_idx,
    const jtf::mesh::face2tet_adjacent &fa_new);

/**
 * @brief
 *
 * @param surface_patch_idx
 * @param surface_normal_align_type
 * @param tet
 * @param node
 * @param outside_face
 * @param outside_face_idx
 * @param fa_new
 * @return int
 */
int get_surface_patch_type_mat(
    matrixst &surface_patch_idx,
    const matrixst &surface_normal_align_type,
    const matrixst & tet,
    const matrixd & node,
    const matrixst & outside_face,
    const matrixst & outside_face_idx,
    const jtf::mesh::face2tet_adjacent &fa_new);

/**
 * @brief
 *
 * @param surface_normal_align_type
 * @param tet
 * @param node
 * @param outside_face
 * @param outside_face_idx
 * @param fa
 * @param frame_array
 * @return int
 */
int get_surface_normal_align_type(
    std::vector<size_t> & surface_normal_align_type,
    const matrixst &tet,
    const matrixd &node,
    const matrixst &outside_face,
    const matrixst &outside_face_idx,
    const jtf::mesh::face2tet_adjacent &fa,
    const std::vector<matrixd > &frame_array);

int get_surface_normal_align_type(
    boost::unordered_map<size_t,size_t> &surface_normal_align_type,
    const matrixst &tet,
    const matrixd &node,
    const matrixst &outside_face,
    const matrixst &outside_face_idx,
    const jtf::mesh::face2tet_adjacent &fa,
    const zjucad::matrix::matrix<matrixd>& frame);

//! @brief this function is used to find near surface singularity
/**
 * @brief
 *
 * @param tet
 * @param node_array
 * @param outside_face
 * @param std::vector<std::deque<std::pair<size_t
 * @param chain_list
 * @param surface_path_idx
 * @param surface_path_align_uvw_idx
 * @param near_surface_singularity_list
 * @return int
 */
int find_near_surface_singularity(
    const matrixst &tet,
    const matrixd &node_array,
    const matrixst & outside_face,
    const std::vector<std::deque<std::pair<size_t,size_t> > > &chain_list,
    const std::vector<size_t> &surface_path_idx,
    const std::vector<size_t> &surface_path_align_uvw_idx,
    std::vector<size_t> &near_surface_singularity_list);

//! @brief this function is used to find adjacent face of one point on surface mesh
/**
 * @brief
 *
 * @param outside_face
 * @param ea
 * @param std::map<size_t
 * @param point_adjacent_tri
 * @return int
 */
int find_adjacent_face_of_point(
    const matrixst &outside_face,
    const jtf::mesh::edge2cell_adjacent &ea,
    std::map<size_t,std::set<size_t> > & point_adjacent_tri);

//! @brief this function is used to check whether two points are in the same situation
// if the set of adjacent patches of v0 and v1 intersect, they are compatible
// others are not
/**
 * @brief
 *
 * @param v0
 * @param v1
 * @param adj_tri_0
 * @param adj_tri_1
 * @param surface_path_idx
 * @param surface_path_align_uvw_idx
 * @return bool
 */
bool is_compatible_situation(
    const size_t v0,const size_t v1,
    const std::vector<size_t> &adj_tri_0,
    const std::vector<size_t> & adj_tri_1,
    const std::vector<size_t> &surface_path_idx,
    const std::vector<size_t> &surface_path_align_uvw_idx);

/**
 * @brief
 *
 * @param v0
 * @param v1
 * @param outside_face
 * @return bool
 */
inline bool is_outside_edge(const size_t v0, const size_t v1,
                            const matrixst &outside_face)
{
  for(size_t t = 0; t < outside_face.size(2); ++t){
    if(std::find(outside_face(zjucad::matrix::colon(),t).begin(),
                 outside_face(zjucad::matrix::colon(),t).end(),v0)
       != outside_face(zjucad::matrix::colon(),t).end()
       && std::find(outside_face(zjucad::matrix::colon(),t).begin(),
                    outside_face(zjucad::matrix::colon(),t).end(),v1)
       != outside_face(zjucad::matrix::colon(),t).end() )
      return true;
  }
  return false;
}


///
/// @brief get_face_patches_according_to_boundary
/// @param outside_face
/// @param ea
/// @param boundary_edges
/// @param patches
/// @return
///
int get_face_patches_according_to_boundary(
    const zjucad::matrix::matrix<size_t> &outside_face,
    const jtf::mesh::edge2cell_adjacent & ea,
    const boost::unordered_set<std::pair<size_t,size_t> > &boundary_edges,
    std::vector<std::vector<size_t> > &patches);

///
/// @brief get_face_patches_according_to_type
/// @param faces
/// @param face_type
/// @param orfap
/// @param patches
/// @return
///
int get_face_patches_according_to_type(
    const zjucad::matrix::matrix<size_t> & faces,
    const zjucad::matrix::matrix<size_t> & face_type,
    const jtf::mesh::one_ring_face_at_point &orfap,
    std::vector<std::vector<size_t> > &patches);

/**
 * @brief
 *
 * @param a
 * @param b
 * @param c
 * @return int
 */
template <typename T1, typename T2>
int common_tet_face(const T1 *a, const T1 *b, T2 *c)
{
  T1 ab[8];
  std::copy(a, a+4, &ab[0]);
  std::copy(b, b+4, &ab[4]);
  std::sort(ab, ab+8);

  int j = 0;
  for(size_t i = 0; i < 7; ++i) {
    if(ab[i+1] == ab[i]) {
      c[j++] = ab[i];
      ++i;
    }
  }
  if(j != 3)
    return 1;

  return 0;
  // the above code can improve the performance a little.
  std::map<T1, char> count;
  for(int i = 0 ; i < 4; ++i) {
    count[a[i]] += 1;
    count[b[i]] += 1;
  }
  if(count.size() != 5) {
    std::cerr << "error in common_tet_face." << std::endl;
    return __LINE__;
  }
  int ni = 0;
  for(typename std::map<T1, char>::const_iterator iter = count.begin();
      iter != count.end(); ++iter) {
    if(iter->second == 1) continue;
    assert(iter->second == 2);
    c[ni++] = iter->first;
  }
  if(ni != 3) {
    std::cerr << "error in common_tet_face." << std::endl;
    return __LINE__;
  }
  return 0;
}

#endif // TOPOLOGY_OPERATION_H
