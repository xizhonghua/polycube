#ifndef TETMESH_REFINE_H
#define TETMESH_REFINE_H
#include <zjucad/matrix/matrix.h>
#include <vector>
#include <deque>
#include <map>
#include <memory>
#include <boost/property_tree/ptree.hpp>
#include <hjlib/function/function.h>
#include <hjlib/function/func_aux.h>
#include "../tetmesh/tetmesh.h"

/**
 * @brief
 *
 * @param tet
 * @param node
 * @param std::vector<std::deque<std::pair<size_t
 * @param singularities
 * @param fa
 * @param outside_face
 * @param pt
 * @return int
 */
int relax_singularities(
    const matrixst & tet,
    matrixd &node,
    const std::vector<std::deque<std::pair<size_t,size_t> > > &singularities,
    const jtf::mesh::face2tet_adjacent &fa,
    const matrixst &outside_face,
    boost::property_tree::ptree &pt);

/**
 * @brief design the singularities coord, it can be get from simply laplacing
 *        or get from input
 *
 * @param node
 * @param std::vector<std::deque<std::pair<size_t
 * @param singularities
 * @param std::map<size_t
 * @param design_singularities_coord
 * @param pt
 * @return int
 */
int design_singularities_coord(
    const matrixd &node,
    const std::vector<std::deque<std::pair<size_t,size_t> > > &singularities,
    std::map<size_t,matrixd > &design_singularities_coord,
    boost::property_tree::ptree &pt);

/**
 * @brief
 *
 */
typedef hj::function::function_t<double, int32_t> function_di;
/**
 * @brief
 *
 */
typedef std::shared_ptr<const function_di> func_ptr;
/**
 * @brief
 *
 */
typedef std::vector<func_ptr> vec_func;
/**
 * @brief
 *
 */
typedef std::shared_ptr<std::vector<func_ptr> > vec_func_ptr;

/**
 * @brief
 *
 * @param vec_all_funcs
 * @param std::map<size_t
 * @param desired_vertex_coord
 * @param tet
 * @param node
 * @param outside_face
 * @param pt
 * @return int
 */
int construct_singularity_relaxing_func(
    vec_func_ptr vec_all_funcs,
    const std::map<size_t,matrixd > &desired_vertex_coord,
    const matrixst &tet,
    const matrixd &node,
    const matrixst &outside_face,
    boost::property_tree::ptree &pt);

/**
 * @brief   arap deformation with a searies constraints that degenerated points
 *          should be concentrated as close as possible
 *
 * @param func store the functions
 * @param tet original tet mesh
 * @param degenerated_points points which have integer constraints conflicted
 * @param node nodes of tet mesh which will be changed after deformation
 * @param pt boost::property_tree::ptree
 * @return int  if works well, return 0, or return non-zeros
 */
int arap_with_degenerated_points_driven(
    vec_func_ptr func,
    const matrixst & tet,
    const std::vector<std::vector<size_t> > & degenerated_points,
    matrixd & node,
    boost::property_tree::ptree &pt);


/**
 * @brief   apply arap deformation with a searies constraints that degenerated
 *          points to concentrate points
 *
 * @param tet original tet mesh
 * @param node nodes of tet mesh which will be changed after deformation
 * @param pt boost::property_tree::ptree
 * @return int  if works well, return 0, or return non-zeros
 */
int arap_to_concentrate_points(
    const matrixst & tet,
    matrixd & node,
    boost::property_tree::ptree &pt);


/**
 * @brief   This function is used to deform a sphere like(genus = 2) model, so
 *          the input should be sphere like
 *
 * @param tet original tet mesh
 * @param node nodes of tet mesh which will be changed after deformation
 * @param pt boost::property_tree::ptree
 * @return int  if works well, return 0, or return non-zeros
 */
int tet_mesh_inflation_for_ball(
    const matrixst & tet,
    matrixd & ndoe,
    boost::property_tree::ptree &pt);

#endif // TETMESH_REFINE_H
