#ifndef POLYCUBE_IO_H
#define POLYCUBE_IO_H

#include <boost/unordered_map.hpp>
#include <boost/unordered_set.hpp>
#include <vector>

#include <jtflib/mesh/mesh.h>
#include <zjucad/matrix/matrix.h>

#include "../tetmesh/tetmesh.h"

int load_surface_restricted_type_static(
    const char * filename,
    boost::unordered_map<size_t,size_t> & result_surface_type);

int dump_surface_restricted_type_to_vtk_static(
    const char * filename,
    const std::string & type_name,
    const zjucad::matrix::matrix<double>& node,
    const jtf::mesh::face2tet_adjacent & fa,
    const boost::unordered_map<size_t,size_t> & surface_type);

int load_fnode_group_static(
    const char * filename,
    const zjucad::matrix::matrix<size_t> & cut_tet,
    std::vector<boost::unordered_set<size_t> > & fnode_group);

int load_point_diffused_value(
    const char * filename,
    zjucad::matrix::matrix<double> & point_weight_diffused);

int load_surface_patch(
    const char * surface_patch_file,
    const char * obj_file,
    const char * s2v_file,
    std::vector<zjucad::matrix::matrix<size_t> > &patch_faces);

int save_to_uv(const char * filename,
               const std::vector<size_t> &uv_basis_point,
               const std::vector<size_t> &boundary_nodes,
               const zjucad::matrix::matrix<double> &polycube_node);

int load_from_uv(const char * filename,
                 const jtf::mesh::meshes & trm,
                 zjucad::matrix::matrix<double> & uv,
                 zjucad::matrix::matrix<double> * orig_node = 0,
                 zjucad::matrix::matrix<double> * uv_basis = 0);

int load_inner_face_jump_type_(
    const char * filename,
    boost::unordered_map<std::pair<size_t,size_t>,size_t> & inner_face_jump_type);
#endif
