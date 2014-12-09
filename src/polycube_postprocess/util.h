#ifndef POLYCUBE_UTIL_H
#define POLYCUBE_UTIL_H

#include <vector>
#include <boost/unordered_map.hpp>
#include <boost/unordered_set.hpp>
#include <math.h>
#include <limits>
#include <zjucad/matrix/matrix.h>

#include <jtflib/mesh/mesh.h>
#include <jtflib/mesh/util.h>

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

#endif
