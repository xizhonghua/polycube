#ifndef FIND_SINGULARITIES_H
#define FIND_SINGULARITIES_H

#include <zjucad/matrix/matrix.h>
#include <zjucad/matrix/itr_matrix.h>
#include <vector>
#include <deque>
#include <map>
#include <set>
#include <boost/unordered_set.hpp>
#include <boost/unordered_map.hpp>

#include <jtflib/mesh/mesh.h>
#include <jtflib/util/util.h>

#include "../tetmesh/tetmesh.h"
#include "../tetmesh/util.h"



class singularity_extractor
{
public:
  singularity_extractor(const jtf::tet_mesh &tm):tm_(tm) {}
  void extract(const zjucad::matrix::matrix<double> & zyz,
               std::vector<std::deque<std::pair<size_t,size_t> > > &chain,
               std::vector<std::deque<size_t> > &singularities_type)const;
  void extract(const zjucad::matrix::matrix<zjucad::matrix::matrix<double> > & frame,
               std::vector<std::deque<std::pair<size_t,size_t> > > &chain,
               std::vector<std::deque<size_t> > &singularities_type)const;
  void extract(const boost::unordered_map<std::pair<size_t,size_t>, size_t> & inner_type,
               std::vector<std::deque<std::pair<size_t,size_t> > > &chain,
               std::vector<std::deque<size_t> > &singularities_type,
               bool with_unknown_face_type = false)const;
private:
  const jtf::tet_mesh &tm_;
};


/**
 * @brief This function is used to find surface singularity edge by using
 *        surface normal alignment type
 * @param outside_face  input outside face
 * @param outside_face_align_type input outside face normal align type
 * @param e2t  input edge to triangle adjacnet releations
 * @param singularity_chain output singularity chain
 * @param singularity_type  output singularity chain type
 * @return int  if works fine return 0, or return non-zeros
 */
int find_surface_singularity_edge_using_normal_align_type(
    const matrixst & outside_face,
    const std::vector<size_t> & outside_face_align_type,
    const jtf::mesh::edge2cell_adjacent & e2t,
    std::vector<std::deque<std::pair<size_t,size_t> > > & singularity_chain,
    std::vector<std::deque<size_t> > & singularity_type);

/**
 * @brief
 *
 * @param sc_file
 * @param singularity_chain
 * @param singualrity_type
 * @return int
 */
int dump_out_singularity_chain_2(const char * sc_file,
                                 const std::vector<std::deque<size_t> > &singularity_chain,
                                 const std::vector<std::deque<size_t> > &singualrity_type);


int find_singularities_global_align_with_type(
    const jtf::mesh::one_ring_tet_at_edge & ortae,
    const matrixst & tet_rot_type,
    const boost::unordered_map<std::pair<size_t,size_t>,size_t > &inner_face_jump,
    const matrixst &outside_face,
    std::vector<std::deque<std::pair<size_t,size_t> > > &chain_list,
    std::vector<std::deque<size_t> > &singularities_type);

/**
 * @brief   Extract chains from edges which keep end of each chain either on surface
 *              or linked with other chains
 *
 * @param segments  input al edges
 * @param outside_face  input outside_face
 * @param chain_list    output chain_list
 * @return int
 */
int extract_chain_from_edges_with_outside_points(
    const boost::unordered_set<std::pair<size_t,size_t> > & segments,
    const matrixst &outside_face,
    std::vector<std::deque<std::pair<size_t,size_t> > > &chain_list);

size_t get_edge_type_with_part_face_type_map(
    const jtf::mesh::one_ring_tet_at_edge & ortae,
    const std::pair<size_t,size_t> & edge,
    const boost::unordered_map<std::pair<size_t,size_t>,size_t> & face_type);

size_t get_edge_type_with_part_face_type_map(
    const std::vector<size_t> & tet_seq,
    const boost::unordered_map<std::pair<size_t,size_t>,size_t> & face_type );

size_t get_edge_type_with_full_face_type_map(
    const jtf::mesh::one_ring_tet_at_edge & ortae,
    const std::pair<size_t,size_t> & edge,
    const boost::unordered_map<std::pair<size_t,size_t>,size_t> & face_type);

#endif // FIND_SINGULARITIES_H
