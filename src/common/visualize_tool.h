#ifndef VISUALIZE_TOOL_H
#define VISUALIZE_TOOL_H


#include <boost/dynamic_bitset.hpp>
#include <tuple>

#include <zjucad/matrix/io.h>
#include <zjucad/matrix/itr_matrix.h>
#include <jtflib/mesh/io.h>

#include "../numeric/util.h"
#include "../hex_param/find_singularities.h"
#include "../tetmesh/util.h"
/**
 * @brief
 *
 * @param file_name
 * @param node
 * @param std::vector<std::deque<std::pair<size_t
 * @param singularity_edges
 * @return int
 */
int  dump_singularity_to_vtk(
    const char * file_name,
    const matrixd &node,
    const std::vector<std::deque<std::pair<size_t,size_t> > > &singularity_edges);

///
/// @brief dump_singularity_to_vtk this function requires input edge to be a set
///        of tuple <p0,p1,type>
/// @param vtk_file output vtk_file
/// @param node input node
/// @param edges
/// @return 0 if ok, or return non-zeros
///
int dump_singularity_to_vtk(
    const char *vtk_file,
    const matrixd & node,
    const std::set<std::tuple<size_t,size_t,size_t> > & edges);

/**
 * @brief
 *
 * @param vtk_file
 * @param node
 * @param singularity_edges
 * @param singularity_type
 * @return int
 */
int dump_singularity_to_vtk(
    const char * vtk_file,
    const matrixd &node,
    const std::vector<std::deque<size_t> > &singularity_edges,
    const std::vector<std::deque<size_t> > &singularity_type);

/**
 * @brief
 *
 * @param vtk_file
 * @param node
 * @param std::map<std::pair<size_t
 * @param size_t>
 * @param edge_info
 * @return int
 */
int dump_singularity_to_vtk(
    const char * vtk_file,
    const matrixd &node,
    const std::map<std::pair<size_t,size_t>,size_t> &edge_info);

/**
 * @brief
 *
 * @param vtk_file
 * @param node
 * @param std::vector<std::deque<std::pair<size_t
 * @param singularity_edges
 * @param singularity_type
 * @return int
 */
int dump_singularity_to_vtk(
    const char * vtk_file,
    const matrixd &node,
    const std::vector<std::deque<std::pair<size_t,size_t> > > &singularity_edges,
    const std::vector<size_t> &singularity_type);


int dump_singularity_to_cylinder(
    const char * obj_file,
    const matrixd &node,
    const std::vector<std::deque<std::pair<size_t,size_t> > > &singularity_edges,
    const double radious);
/**
 * @brief
 *
 * @param vtk_file
 * @param node
 * @param std::vector<std::deque<std::pair<size_t
 * @param singularity_edges
 * @param singularity_type
 * @return int
 */
int dump_singularity_chain_to_vtk_2(
    const char * vtk_file,
    const matrixd &node,
    const std::vector<std::deque<std::pair<size_t,size_t> > > &singularity_edges,
    const std::vector<std::deque<size_t> > &singularity_type);

int dump_singularity_chain_to_vtk_2_type(
    const char * vtk_file,
    const matrixd &node,
    const std::vector<std::deque<std::pair<size_t,size_t> > > &singularity_edges,
    const std::vector<std::deque<size_t> > &singularity_type);

int dump_singularity_chain_to_vtk_3(
    const char * vtk_file,
    const matrixd &node,
    const std::vector<std::deque<std::tuple<size_t,size_t,size_t> > > &singularity_edges);
/**
 * @brief
 *
 * @param vtk_file
 * @param tet
 * @param node
 * @param uvw
 * @return int
 */
int dump_uvw(const char* vtk_file,
             const matrixst &tet,
             const matrixd &node,
             const matrixd &uvw);

/**
 * @brief
 *
 * @param vtk_file
 * @param tet
 * @param node
 * @param uvw
 * @param freq
 * @return int
 */
int dump_uvw_wave(const char * vtk_file,
                  const matrixst &tet,
                  const matrixd &node,
                  const matrixd &uvw,
                  const double freq = 5);

/**
 * @brief
 *
 * @param vtk_file
 * @param tet_
 * @param node_
 * @param uvw
 * @return int
 */
int dump_uvw_vol(const char * vtk_file,
                 const matrixst &tet_,
                 const matrixd &node_,
                 const matrixd &uvw);

/**
 * @brief
 *
 * @param vtk_file
 * @param frame_in_tet
 * @param tet
 * @param node
 * @param uvw
 * @return int
 */
int dump_frame_align_error(
    const char * vtk_file,
    const zjucad::matrix::matrix<matrixd > &frame_in_tet,
    const matrixst &tet,
    const matrixd &node,
    const matrixd &uvw);


/**
 * @brief dump out frame difference to vtk
 *
 * @param vtk_file  output vtk file
 * @param tet input tet mesh
 * @param node  input node
 * @param fa  inputjtf::mesh::face2tet_adjacent
 * @param frame input frame
 * @return int
 */
int dump_frame_difference_to_vtk(
    const char * vtk_file,
    const matrixst & tet,
    const matrixd & node,
    const jtf::mesh::face2tet_adjacent & fa,
    const zjucad::matrix::matrix<matrixd >& frame);


int dump_sh_difference_to_vtk(
    const char * vtk_file,
    const matrixst & tet,
    const matrixd & node,
    const jtf::mesh::face2tet_adjacent & fa,
    const matrixd& sh);
/**
 * @brief
 *
 * @param file_name
 * @param std::vector<std::deque<std::pair<size_t
 * @param singularities_chain
 * @param singularities_type
 * @return int
 */
int dump_out_singularity_chain(
    const char * file_name,
    const std::vector<std::deque<std::pair<size_t,size_t> > > &singularities_chain,
    const std::vector<size_t> & singularities_type);

/**
 * @brief
 *
 * @param file_name
 * @param singularities_chain
 * @param singularities_type
 * @return int
 */
int dump_out_singularity_chain_2(
    const char * file_name,
    const std::vector<std::deque<size_t> > &singularities_chain,
    const std::vector<std::deque<size_t> > & singularities_type);

/**
 * @brief
 *
 * @param file_name
 * @param std::vector<std::deque<std::pair<size_t
 * @param singularities_chain
 * @param singularities_type
 * @return int
 */
int dump_out_singularity_chain_3(
    const char * file_name,
    const std::vector<std::deque<std::pair<size_t,size_t > > > &singularities_chain,
    const std::vector<std::deque<size_t> > & singularities_type);

/**
 * @brief dump out singularity chain and the beginning and ending tets around each edge
 *
 * @param file_name output file
 * @param ortae input one_ring_tet_at_edge info
 * @param singularities_chain
 * @param singularities_type
 * @return int  return 0 if works fine, or return non-zeros
 */
int dump_out_singularity_chain_with_tet_ends(
    const char * file_name,
    const jtf::mesh::one_ring_tet_at_edge & ortae,
    const std::vector<std::deque<std::pair<size_t,size_t > > > &singularities_chain,
    const std::vector<std::deque<size_t> > & singularities_type);

/**
 * @brief
 *
 * @param file_name
 * @param node
 * @param std::vector<std::deque<std::pair<size_t
 * @param singularities_chain
 * @param singularities_loop
 * @param singularities_type
 * @param frame_inner
 * @return int
 */
int dump_out_singularity_axis(
    const char *file_name,
    const matrixd &node,
    const std::vector<std::deque<std::pair<size_t,size_t> > > &singularities_chain,
    const std::vector<std::vector<size_t> > &singularities_loop,
    const std::vector<size_t> & singularities_type,
    const zjucad::matrix::matrix<matrixd > & frame_inner);

/**
 * @brief
 *
 * @param file
 * @param faces
 * @param fa
 * @param cut_tet2tet
 * @param node
 * @return int
 */
int dump_out_faces_in_cut_tet(const char * file,
                              const std::vector<size_t> faces,
                              const jtf::mesh::face2tet_adjacent &fa,
                              const matrixst &cut_tet2tet,
                              const matrixd &node);

/**
 * @brief
 *
 * @param file_name
 * @param std::vector<std::deque<std::pair<size_t
 * @param singularities_chain
 * @param singularities_type
 * @return int
 */
int load_singularity_chain(
    const char * file_name,
    std::vector<std::deque<std::pair<size_t,size_t> > > &singularities_chain,
    std::vector<size_t> & singularities_type);

/**
 * @brief
 *
 * @param file_name
 * @param std::vector<std::deque<std::pair<size_t
 * @param singularities_chain
 * @param singularities_type
 * @return int
 */
int load_singularity_chain_new(
    const char * file_name,
    std::vector<std::deque<std::pair<size_t,size_t> > > &singularities_chain,
    std::vector<std::deque<size_t> > & singularities_type);

//! @brief this function is used to convert directional edge into arrows
template <typename T1, typename T2>
void directional_edge2arrow(const T1 *edges_ptr,  const size_t edge_num,
                            const T2 *nodes_ptr, const size_t node_num,
                            const char * arrow_obj_file,
                            const char * output_file)
{
  // this function assume inout arrow_obj_file is inside x\in [0,1], and
  // the arrow is at x = 1
  zjucad::matrix::matrix<size_t> arrow_mesh;
  zjucad::matrix::matrix<double> arrow_node;
  if(jtf::mesh::load_obj(arrow_obj_file, arrow_mesh, arrow_node))
    throw std::invalid_argument("# [error] can not load arrow_obj_file.");

  for(size_t ni = 0; ni < arrow_node.size(2); ++ni)
    arrow_node(zjucad::matrix::colon(),ni) -= arrow_node(zjucad::matrix::colon(),0);
  std::vector<size_t> dyn_mesh;
  std::vector<double> dyn_node;

  zjucad::matrix::itr_matrix<const T1*> edges(2, edge_num, edges_ptr);
  zjucad::matrix::itr_matrix<const T2*> nodes(3, node_num, nodes_ptr);

  zjucad::matrix::matrix<double> temp_arrow_node;
  zjucad::matrix::matrix<size_t> temp_arrow_mesh;
  zjucad::matrix::matrix<double> x_axis = zjucad::matrix::zeros<double>(3,1);
  x_axis[0] = 1.0;
  zjucad::matrix::matrix<double> dir(3,1);
  zjucad::matrix::matrix<double> axis(3,1);
  zjucad::matrix::matrix<double> rot(3,3);

  for(size_t ei = 0; ei < edges.size(2); ++ei){
      dir = nodes(zjucad::matrix::colon(),edges(1,ei)) -
          nodes(zjucad::matrix::colon(),edges(0,ei));
      const double len = zjucad::matrix::norm(dir);
      dir /= (len > 1e-6?len:1.0);
      temp_arrow_node = arrow_node;
      temp_arrow_mesh = arrow_mesh;
      temp_arrow_node *= len;
      axis = zjucad::matrix::cross(x_axis, dir);
      if(zjucad::matrix::norm(axis) < 1e-5) {// x_axis is parallel with dir
          axis = zjucad::matrix::zeros<double>(3,1);
          axis[2] = 1.0;
        }
      const double angle =
          acos(std::min(1.0,std::max(-1.0, zjucad::matrix::dot(x_axis, dir))));

      from_angle_to_rotation_matrix(angle, axis, rot);

      temp_arrow_node = zjucad::matrix::temp(rot * temp_arrow_node);

      for(size_t ni = 0; ni < temp_arrow_node.size(2); ++ni)
        temp_arrow_node(zjucad::matrix::colon(),ni) +=
            nodes(zjucad::matrix::colon(), edges(0,ei)) ;
      temp_arrow_mesh += dyn_node.size()/3;
      dyn_mesh.insert(dyn_mesh.end(), temp_arrow_mesh.begin(), temp_arrow_mesh.end());
      dyn_node.insert(dyn_node.end(), temp_arrow_node.begin(), temp_arrow_node.end());
    }

  zjucad::matrix::itr_matrix<const size_t * > dyn_mesh_m(3, dyn_mesh.size()/3, &dyn_mesh[0]);
  zjucad::matrix::itr_matrix<const double * > dyn_node_m(3, dyn_node.size()/3, &dyn_node[0]);
  if(jtf::mesh::save_obj(output_file, dyn_mesh_m, dyn_node_m )){
      throw std::invalid_argument("# [error] can not open obj file.");
    }
}

std::string bitset_to_hex_string(const boost::dynamic_bitset<> & bit);
#endif // VISUALIZE_TOOL_H
