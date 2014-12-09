#ifndef CUT_TET_H
#define CUT_TET_H

#include <map>
#include <zjucad/matrix/matrix.h>
#include "../tetmesh/tetmesh.h"
#include "../common/graph.h"
#include "../hex_param/hex_param.h"


struct Euler_characteristic
{
  int node_num, edge_num, face_num;
  int operator()(void) const {
    return node_num - edge_num + face_num;
  }
};

class tetmesh_cutter
{
public:
  tetmesh_cutter(const jtf::tet_mesh &tm):tm_(tm){}
  void cut(zjucad::matrix::matrix<double> & zyz, bool only_mst = true);
  void cut(const boost::unordered_map<std::pair<size_t,size_t>,size_t> & inner_type);
  jtf::mesh::meshes cut_tm_;
private:
  void init();
  void update_cut_node();
  /**
   * @return 0: successfully glue tet a and b, 1: already glued, 2:
   * topology changed
   *
   * NOTICE, when there're 2 common nodes, the edge of decrease may be 2
   * or 3, depends on whether the edge has been fully surrounded.
   * However, since we restrict the surface of the jtf::mesh::meshes is manifold,
   * so edge_num alwasy decreases by 3.
   */
  int glue_tet(size_t t0, size_t t1, bool is_span_tree_edge);
private:
  const jtf::tet_mesh &tm_;
  Euler_characteristic ec_;
  vector<set<size_t*> > acc_table_;
  fix_graph fg_;
};
///
/// @brief cut_tet cut tet mesh open according to type settings
/// @param tet input tetmesh
/// @param fa  input face tet_adjacent
/// @param inner_face_jump_type input surface_type
/// @param surface_type input surface type
/// @param is_restrcited_type control surface type to be restrited type or not
/// @param cut_tet output cut tet mesh
/// @param frame_ptr output  frame field
///
void cut_tetmesh(
    const matrixst & tet,
    const matrixd & node,
    const jtf::mesh::face2tet_adjacent & fa,
    boost::unordered_map<std::pair<size_t,size_t>,size_t> &inner_face_jump_type,
    boost::unordered_map<size_t,size_t> & surface_type,
    const bool is_restrcited_type,
    matrixst & cut_tet,
    bool only_mst = false,
    zjucad::matrix::matrix<matrixd> *frame_ptr = 0);

///**
// * @brief   cut tet with minial cuts, which means :
// *          1. each singularity crease is outside on surface
// *          2. each non-singularity crease can not be glued anymore
// *
// * @param original_tet  input tet mesh
// * @param node  input node
// * @param inner_face_jump_type  input inner face jump type
// * @param fa inputjtf::mesh::face2tet_adjacent
// * @param ortae input one_ring_tet_at_edge info
// * @param singularity_chain input singularity chains
// * @param cut_tet output cut tet
// * @param face_pair store face pair for each original face
// * @param face_type output face type
// * @param param
// * @return int  return 0 if works fine, or return non-zeros
// */
//int minimal_cut_tet(
//    const jtf::tet_mesh &tm,
//    const zjucad::matrix::matrix<double > & zyz,
//    matrixst & cut_tet,
//    matrixst & face_pair,
//    matrixst & face_type,
//    const matrixd * param = 0);

int glue_tet_mst_with_type(
    const matrixst &original_tet,
    const jtf::mesh::face2tet_adjacent & fa,
    std::set<std::pair<size_t,size_t> > & spanning_edge,
    fix_graph &fg,
    matrixst &output_cut_tet,
    Euler_characteristic &ec,
    std::vector<std::set<size_t *> > &acc_table,
    boost::unordered_map<std::pair<size_t,size_t>, size_t> & inner_face_jump_type,
    boost::unordered_map<size_t,size_t> & surface_type,
    const size_t surface_type_flag, // 0: restricted type, 1 : normal aligned type
    zjucad::matrix::matrix<matrixd> * frame_ptr = 0);
/**
 * @brief   WARNING:  This function takes an assumption: cut tet has been glued
 *                    with a minimal spanning tree !!!!
 *          glue tet and avoid gluing singularity edge inner, that means all
 *          singularity crease exist on surface, and non-singularity crease will
 *          be reduced into the maximum extent.
 *
 * @param original_tet  input tet mesh
 * @param inner_face_jump_type  input inner face jump type
 * @param fa  inputjtf::mesh::face2tet_adjacent
 * @param spanning_edges  spanning tree edges
 * @param ortae input one_ring_tet_at_edge
 * @param fg  input fix graph
 * @param dorfap input dynamical one ring faec at point
 * @param ec            store euler characteristic to keep topology like a sphere
 * @param acc_table
 * @return int  return 0 if works fine, or return non-zeros
 */
int glue_tet_avoid_singularity_crease(
    const matrixst &original_tet,
    const zjucad::matrix::matrix<matrixd > & frame,
    const jtf::mesh::face2tet_adjacent &fa,
    const std::vector<std::deque<std::pair<size_t,size_t> > > & singularity_edges,
    const jtf::mesh::one_ring_tet_at_edge & ortae,
    fix_graph &fg,
    jtf::mesh::dyn_one_ring_face_at_point & dorfap,
    matrixst &output_cut_tet,
    Euler_characteristic &ec,
    std::vector<std::set<size_t *> > &acc_table);

/**
 * @brief   This function is used to glue all non-spanning tree edge(here edge
 *          means adjacent tet pair), and keep the toplogy like a sphere
 *
 * @param original_tet  input tet mesh
 * @param output_cut_tet output cut tet
 * @param tet_a  input tet a of adjacent tet pair
 * @param tet_b  input tet b of adjacent tet pair
 * @param ec     store euler characteristic to keep topology like a sphere
 * @param acc_table
 * @param is_span_tree_edge determin whether this tet pair is on spanning edge
 * @param singularity_edges all singularity edges
 * @param ortae  input one rign tet at edge info
 * @param dorfap input dynamical_one_ring_face_at_point relations
 * @return int  return 0 if works fine, or return non-zeros
 */
int glue_each_tet_pair_avoid_singularity_edge(
    const matrixst & original_tet,
    matrixst & output_cut_tet,
    const size_t & a, const size_t & b,
    Euler_characteristic & ec,
    std::vector<std::set<size_t* > > & acc_table,
    bool is_span_tree_edge,
    const std::set<std::pair<size_t,size_t> > & singularity_edges,
    const jtf::mesh::one_ring_tet_at_edge & ortae,
    jtf::mesh::dyn_one_ring_face_at_point & dorfap);



/**
 * @brief   Analysis the face type of cut tet surface
 *
 * @param input_tet  input tet mesh
 * @param input_node input node
 * @param cut_tet    input cut tet mesh
 * @param input_fa   inputjtf::mesh::face2tet_adjacent
 * @param inner_face_jump_type  input inner_face_jump_type
 * @param singularity_edges     input all singularity edges
 * @param face_pair  output face pair, if one face has -1 as its pair, that means
 *                   it's original surface face
 * @param face_type  possible output face type
 * @return int  return 0 if works fine, or return non-zeros
 */
int analysis_face_type(
    const matrixst & input_tet,
    const matrixd & input_node,
    const matrixst & cut_tet,
    const jtf::mesh::face2tet_adjacent &input_fa,
    const boost::unordered_map<std::pair<size_t,size_t>, size_t> & inner_face_jump_type,
    const std::set<std::pair<size_t,size_t> > & singularity_edges,
    matrixst & face_pair,
    matrixst & face_type,
    const matrixd * param = 0);



int analysis_each_surface_face_group(
    const size_t & group_idx,
    const std::set<size_t> & one_group,
    const matrixd & node,
    const std::set<std::pair<size_t,size_t> > & singularity_edges_original,
    const jtf::mesh::face2tet_adjacent & fa_cut_tet,
    const jtf::mesh::edge2cell_adjacent & ea_original_surface,
    const matrixst & cut_tet2tet,
    const matrixd * param = 0);

void analysis_transition_raw(
    const matrixst &input_tet,
    const matrixd &input_node,
    const jtf::mesh::face2tet_adjacent &input_fa,
    const jtf::mesh::face2tet_adjacent & fa_cut,
    const matrixst &cut_tet,
    matrixst &face_pair);

int glue_tet_with_given_cut_faces(
    const matrixst &tet,
    const matrixd & node,
    const std::vector<std::pair<size_t,size_t> > & cut_tet_pair,
    matrixst &new_cut_tet);
#endif // CUT_TET_H

