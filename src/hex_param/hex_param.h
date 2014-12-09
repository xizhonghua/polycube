#ifndef HJ_HEX_PARAM_H_
#define HJ_HEX_PARAM_H_

#include "find_singularities.h"
#include <vector>
#include <string>
#include <deque>
#include <zjucad/matrix/matrix.h>

#include <boost/property_tree/ptree.hpp>
#include "../tetmesh/tetmesh.h"
#include "../tetmesh/util.h"
#include "../hex_frame_opt/hex_frame_opt.h"
#include "../common/graph.h"
#include <jtflib/mesh/mesh.h>

int refine_frame_field_after_aligned_new(boost::property_tree::ptree &pt,
                                         jtf::tet_mesh &tm,
                                         const matrixd &zyz,
                                         matrixd & new_zyz);

//! @brief this function is used to remove black lines and zigzag by splitting tets
// this method can make sure to clean all black lines and zigzag,
// but the final singularity distribution is not good enough
int refine_frame_field_after_aligned_split(boost::property_tree::ptree &pt,
                                           jtf::tet_mesh &tm,
                                           const matrixd &zyz,
                                           matrixd & new_zyz);



//int dump_surface_singularity_to_vtk(const char * vtk_file,
//                                    const matrixd &node,
//                                    const std::vector<std::deque<std::pair<size_t,size_t> > > &singularity_chain,
//                                    const std::vector<std::deque<size_t> > &singularity_type);


int solve_jump_face(const matrixst &cut_tet,
                    const matrixd &cut_node,
                    const matrixst &face_pair,
                    matrixst &jump_face_type,
                    const matrixst &outside_face_cut,
                    const matrixst &outside_face_cut_idx,
                    const matrixst &cut_tet2tet,
                    const jtf::mesh::face2tet_adjacent &fa_cut,
                    const std::vector<std::deque<std::pair<size_t,size_t> > > &singularity_chain,
                    const std::vector<std::deque<size_t> > &singularity_type);

int solve_jump_face_using_zyz(
    const matrixst &cut_tet,
    const matrixd &cut_node,
    const matrixst &face_pair,
    matrixst &jump_face_type,
    const matrixst &outside_face_cut,
    const matrixst &outside_face_cut_idx,
    const matrixst &cut_tet2tet,
    const jtf::mesh::face2tet_adjacent &fa_cut,
    const std::vector<std::deque<std::pair<size_t, size_t> > > &singularity_chain,
    const std::vector<std::deque<size_t> > &singularity_type,
    const zjucad::matrix::matrix<matrixd > &frame_inner);

int solve_jump_face_using_zyz_new(
    const matrixst &cut_tet,
    const matrixd &cut_node,
    const matrixst &tet,
    const matrixd &node,
    const matrixst &face_pair,
    matrixst &jump_face_type,
    const matrixst &outside_face_cut,
    const matrixst &outside_face_cut_idx,
    const matrixst &cut_tet2tet,
    const jtf::mesh::face2tet_adjacent &fa_cut,
    const jtf::mesh::face2tet_adjacent &fa,
    std::map<std::pair<size_t,size_t>,size_t> &singularity_segment,
    const zjucad::matrix::matrix<matrixd > &frame_inner,
    const jtf::mesh::one_ring_tet_at_edge &ortae,
    const boost::property_tree::ptree& pt,
    std::map<std::pair<size_t,size_t>,size_t> &tet_pair_type);


int split_tet_at_zigzag(const matrixst &tet,
                        const matrixd &node,
                        const matrixd &zyz,
                        const std::vector<std::deque<std::pair<size_t,size_t> > > &singularity_chain,
                        matrixst &new_tet,
                        matrixd &new_node,
                        matrixd &new_zyz);

int remove_surface_zigzag(
    const matrixst &outside_face,
    const matrixst &tet,
    const matrixd &node,
    const jtf::mesh::face2tet_adjacent &fa,
    const boost::unordered_map<std::pair<size_t,size_t>,size_t> &jump_type_between_tets,
    const matrixst &outside_face_idx,
    const zjucad::matrix::matrix<matrixd > &frame_inner,
    const jtf::mesh::one_ring_tet_at_edge &ortae,
    const boost::property_tree::ptree &pt,
    boost::unordered_map<size_t,size_t> &face_jump);

//matrixd type_transition2(const size_t type);
//size_t type_transition1(const matrixd& transition);

#include "../common/graph.h"
void orient_frame(const matrixd &zyz_frame_in_tet,
                  zjucad::matrix::matrix<matrixd > &frame_in_tet,
                  const fix_graph &graph,
                  matrixd &orient_err);

void orient_frame(matrix<matrixd > &frame_in_tet,
                  const fix_graph &graph,
                  matrixd &orient_err);

void orient_type(
    const fix_graph &graph,
    const size_t & tet_num,
    const jtf::mesh::face2tet_adjacent & fa,
    boost::unordered_map<std::pair<size_t,size_t>,size_t> & inner_face_jump_type,
    boost::unordered_map<size_t,size_t> & surface_type,
    const size_t surface_type_flag, // 0: restricted type, 1: normal aligned type
    zjucad::matrix::matrix<matrixd> * frame_ptr = 0);

int tet2dual_graph(const matrixst &tet,
                   fix_graph &fg,
                   const jtf::mesh::face2tet_adjacent &fa);
void build_minimum_spanning_tree(const fix_graph &fg,
                                 std::vector<double> &weights,
                                 std::vector<size_t> &edge_list);
void subgraph(const fix_graph &graph,
              const std::vector<size_t> &edge_list,
              fix_graph &sub);
#endif
