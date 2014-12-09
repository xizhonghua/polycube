#ifndef _POLYCUBE_UTIL_H_
#define _POLYCUBE_UTIL_H_

#include <vector>
#include <boost/unordered_map.hpp>
#include <boost/unordered_set.hpp>
#include <math.h>
#include <limits>
#include <zjucad/matrix/matrix.h>
#include <jtflib/mesh/mesh.h>
#include <jtflib/mesh/util.h>

#include "../common/def.h"

template <typename T>
bool non_zero(T & d){
    return fabs(d) > 1e-6;
}

bool is_degenerate(const zjucad::matrix::matrix<double> &tri);

int convert_surface_type_to_surface_patches(
        const zjucad::matrix::matrix<size_t> &tet,
        const boost::unordered_map<size_t,size_t> &surface_type,
        std::vector<zjucad::matrix::matrix<size_t> > &surface_patches,
        boost::unordered_set<std::pair<size_t,size_t> > &patch_boundary);

int smooth_boundary(
        const boost::unordered_set<std::pair<size_t,size_t> > & boundary,
        const zjucad::matrix::matrix<double> & orig_node,
        zjucad::matrix::matrix<double> & polycube_node,
        const size_t iter = 5);

int smooth_patch(const std::vector<zjucad::matrix::matrix<size_t> > & surface_patch,
                 const boost::unordered_set<std::pair<size_t,size_t> > &boundary,
                 zjucad::matrix::matrix<double> & node,
                 const size_t iter = 5);

int smooth_volume(
        const zjucad::matrix::matrix<size_t> & tet,
        const std::vector<zjucad::matrix::matrix<size_t> > &surface_patches,
        zjucad::matrix::matrix<double> &node,
        const size_t iter_v = 5);

int get_face_patches_according_to_boundary(
        const zjucad::matrix::matrix<size_t> &outside_face,
        const jtf::mesh::edge2cell_adjacent & ea,
        const boost::unordered_set<std::pair<size_t,size_t> > &boundary_edges,
        std::vector<std::vector<size_t> > &patches);

int remove_flipped_tet(
        matrixst &orig_tet,
        matrixd &node,
        matrixst &polycube_tet,
        matrixd &polycube_node);

int dump_out_normal_flipped_face(
        const zjucad::matrix::matrix<double> &node,
        const std::vector<zjucad::matrix::matrix<size_t> > &surface_patches,
        const std::vector<zjucad::matrix::matrix<double> > &patch_normal,
        const size_t i,
        const zjucad::matrix::matrix<double> &vtk_node);

double get_tet_arap_distortion(
        const zjucad::matrix::matrix<double> & orig_node,
        const zjucad::matrix::matrix<double> & new_node,
        const zjucad::matrix::matrix<size_t> & one_tet);

double tet_scaled_jacobian(const zjucad::matrix::matrix<double> & tet_node);

int local_remesh_tet(zjucad::matrix::matrix<size_t> & tet,
                     zjucad::matrix::matrix<double> & orig_node,
                     zjucad::matrix::matrix<double> & deform_node);

int split_tet_at_large_arap_distortion(
        zjucad::matrix::matrix<size_t> &orig_tet,
        zjucad::matrix::matrix<double> &orig_node,
        zjucad::matrix::matrix<double> &deform_node,
        zjucad::matrix::matrix<size_t> &new_node_parent_id,
        const double split_percent = 1.0);


int uniform_vertex_order(std::vector<size_t> &tri_a,
                         std::vector<size_t> &tri_b,
                         const zjucad::matrix::matrix<size_t> *cut_tet2tet = 0);

int load_chain_file(const char * chain_file,
                    std::vector<std::vector<size_t> > & chains);


int dump_out_group_file(const char * group_file,
                        const std::vector<std::vector<size_t> > &groups,
                        const std::vector<size_t> &group_integer);

int dump_out_chain_file(const char * chain_file,
                        const std::vector<std::vector<size_t> > & chain);

//! This function is used to remove surface low quality triangle and the corresponding
// tet, low quality is defined by min angle
// @param min angle threshold degrees
int remove_degenerated_face_of_tet(
        zjucad::matrix::matrix<size_t> &orig_tet,
        zjucad::matrix::matrix<double> &orig_node,
        zjucad::matrix::matrix<double> &def_node,
        zjucad::matrix::matrix<size_t> &tri_faces,
        const double min_angle_threshold = 10.0,
        const double max_angle_threshold = 170.0);
#endif
