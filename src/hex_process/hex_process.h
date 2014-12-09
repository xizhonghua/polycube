#ifndef HEX_PROCESS_H
#define HEX_PROCESS_H

#include <jtflib/mesh/mesh.h>
#include "../tetmesh/tetmesh.h"
#include "../tetmesh/util.h"
#include "../hexmesh/util.h"

#include <memory>
#include <boost/property_tree/ptree.hpp>
#include <boost/function.hpp>
#include <boost/unordered_set.hpp>
#include <zjucad/matrix/matrix.h>



// @brief pare one layer hexmesh from outside surface
int pare_hexmesh_from_surface(const matrixst &hex,
                              const matrixd &node,
                              matrixst &hex_new);

#if 0 // no used
int extrapolate_hexmesh_from_surface(const matrixst &hex,
                                     const matrixd &node,
                                     matrixst &hex_new,
                                     matrixd &node_new,
                                     boost::property_tree::ptree &pt);
#endif

//! @brief return the state of one hex and one face
// 0: outside along the triangle normal
// 1: outside along -1 * triangle normal
// 2: intersect
// not used!
size_t cal_adjacent_state(const matrixst &each_hex,
                          const matrixd &hex_node,
                          const matrixst &face,
                          const matrixd &face_normal,
                          const matrixd &tet_node);

//! @brief define the hit function
class hit_func{
public:
  virtual ~hit_func(){}
  virtual int hex_tri_hit_func(const matrixst &outside_tri_face,
                               const matrixd &tet_node,
                               const matrixst &hex,
                               const matrixd &hex_node,
                               const matrixd &face_normal,
                               matrixst &hex_hit) = 0;
};

class trivial_hit_func : public hit_func{
public:
  //enum dirty_flag {USE_NEW_BOUNDING_BOX,USE_CACHED_BOUNDING_BOX};
  virtual ~trivial_hit_func(){}
  virtual int hex_tri_hit_func(const matrixst &outside_tri_face,
                               const matrixd &tet_node,
                               const matrixst &hex,
                               const matrixd &hex_node,
                               const matrixd &face_normal,
                               matrixst &hex_hit);

  //! @brief use the ray casting method to determin a point is inside the tet or not
  // we
  // return: 0, outside; 1, inside
  size_t ray_casting_determin(const matrixd &point_,
                              const matrixd &direction,
                              const matrixst &outside_tri_face,
                              const matrixd &face_normal,
                              const matrixd &tet_node);

  //! @brief use the bary centering method to determin a point is inside or not.
  bool is_point_in_triangle(const matrixd & point_,
                            const matrixst &triangle,
                            const matrixd &node)const;
private:
  matrixd bounding_box;
};

// @brief pare hex which is outside tet outside surface
// accroding to the original tet surface, pare all hex which are out of tet surface
// for those hex which are on the boundary surface, we need an extra policy to judge
int pare_hex_outside_tet_surface(const matrixst &tet,
                                 const matrixd &tet_node,
                                 const matrixst &ori_hex,
                                 const matrixd &ori_hex_node,
                                 hit_func &func, // func is used to varify whether a hex around surface is in or out
                                 matrixst &new_hex,
                                 matrixd &new_hex_node);

/**
 * @brief  This function is used to deform a hex to make it looks like original
 *         tet. This deformation is based on constrained boundary arap.
 *
 * @param hex       original hex mesh
 * @param hex_node  original hex node
 * @param tet       original tet mesh
 * @param tet_node  original tet node
 * @return int      if works fine return 0, or return non-zeros
 */
int deform_hex_to_original_tet(matrixst & hex,
                               matrixd & hex_node,
                               const matrixst & tet,
                               const matrixd & tet_node,
                               boost::property_tree::ptree &pt);

int create_map_from_hex_to_tet(
    const matrixd &tet_nodes,
    const matrixst &tet_faces,
    const matrixd &hex_nodes,
    const matrixst &hex_faces,
    const matrixst &hex,
    const matrixst & tet_surface_points,
    const jtf::mesh::face2hex_adjacent &fa,
    const jtf::mesh::edge2cell_adjacent &ea,
    const jtf::mesh::one_ring_face_at_point & ring_faces,
    const matrixst & hex_face_points,
    matrixd & nearest_nodes,
    matrixst & nearest_node_faces,
    matrixd & nearest_node_weights);


int one_deformation(
    const matrixst & hex,
    matrixd & hex_node,
    const matrixst & tet,
    const matrixd & tet_node,
    const matrixst & tet_surface,
    const matrixst & hex_surface,
    const matrixst & tet_surface_points,
    const matrixst & hex_surface_points,
    const boost::unordered_set<size_t> & hex_surface_point_set,
    const jtf::mesh::face2hex_adjacent &fa_hex,
    const jtf::mesh::edge2cell_adjacent &ea,
    const jtf::mesh::one_ring_face_at_point & orfap,
    const boost::unordered_map<size_t,boost::unordered_set<size_t> > &point_adj_points,
    const matrixd & tri_face_normal,
    boost::property_tree::ptree &pt);

/**
 * @brief This function is used to cut hexahedral into six tets, no extra points
 *        are introduced
 *
 * @param hex   input hexahedral mesh
 * @param node  input hexahedral node
 * @param tet   output tet mesh
 * @return int  return 0 if works fine, or return non-zeros
 */
int create_tet_from_hex_cut_raw(
    const matrixst & hex,
    const matrixd & node,
    matrixst & tet);


//int merge_vertex(
//    matrixst & hex,
//    const std::vector<std::pair<size_t,size_t> > & merge_vertex_list);

int calc_point_normal(const size_t & tri_face_idx,
                      const matrixd & uvw,
                      const jtf::mesh::one_ring_face_at_point & orfap,
                      const matrixst & tri_faces,
                      const matrixd & tri_nodes,
                      const matrixd & tri_face_normal,
                      matrixd & normal);

/**
 * @brief This function is used to remove invalid hex.
 *        Here invalid hex is determined by a face, where
 *
 * @param hex   input hexahedral mesh
 * @param node  input hexahedral node
 * @param tet   output tet mesh
 * @return int  return 0 if works fine, or return non-zeros
 */
int remove_invalid_hex(matrixst & hex,
                       const matrixd & node);

/**
 * @brief This function is used to project each point of quad feature line to
 *        triangle feature line, and return the nearest point of triangle feature
 *        line for each point of quad feature line
 * TODO:  left for shenxinxin
 *
 * @param quad_feature   input hexahedral mesh
 * @param quad_node      input hexahedral node
 * @param tri_feature    input tet mesh
 * @param tri_node       input triangle node
 * @param q2t            output q2t store nearest node for each point of quad
 * @return int  return 0 if works fine, or return non-zeros
 */
int project_quad_feature_to_tri_feature(
    const matrixst & quad_feature,
    const matrixd & quad_node,
    const matrixst & tri_feature,
    const matrixd & tri_node,
    boost::unordered_map<size_t, matrixd > & q2t);

int project_quad_feature_to_tri_feature_new(
    const matrixst & quad_feature,
    const matrixd & quad_node,
    const matrixst & tri_feature,
    const matrixd & tri_node,
    boost::unordered_map<size_t, matrixd > & q2t);
#endif // HEX_PROCESS_H
