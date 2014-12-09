#ifndef JTF_TETMESH_UTIL_H
#define JTF_TETMESH_UTIL_H

#include <set>
#include <boost/unordered_map.hpp>
#include <boost/unordered_set.hpp>

#include <zjucad/matrix/matrix.h>
#include <jtflib/mesh/mesh.h>
#include <boost/property_tree/ptree.hpp>
#include "tetmesh.h"

namespace jtf{
  namespace tetmesh{

    inline const std::pair<size_t,size_t> face_pair2tet_pair(
        const std::pair<size_t,size_t> & face_pair,
        const jtf::mesh::face2tet_adjacent & fa){
        using namespace std;
      const pair<size_t,size_t> & t0_pair = fa.face2tet_.at(face_pair.first);
      const pair<size_t,size_t> & t1_pair = fa.face2tet_.at(face_pair.second);
      assert(fa.is_outside_face(t0_pair) && fa.is_outside_face(t1_pair));
      pair<size_t,size_t> tet_pair;
      tet_pair.first = (t0_pair.first==-1?t0_pair.second:t0_pair.first);
      tet_pair.second = (t1_pair.first==-1?t1_pair.second:t1_pair.first);
      return tet_pair;
    }
    ///
    /// @brief smooth_one_ring_point_normal
    /// @param input/output point_normal
    /// @param one_ring_point_of_p
    /// @param iter_num
    /// @return
    ///
    int smooth_one_ring_point_normal(
        zjucad::matrix::matrix<double> &point_normal,
        const boost::unordered_map<size_t, boost::unordered_set<size_t> > &one_ring_point_of_p,
        const size_t iter_num = 3);

    class one_ring_tet_at_point
    {
    public:
      void add_tets(const matrixst & tets);
      std::vector<std::vector<size_t> > p2t_;
    };

    class one_ring_edge_at_point
    {
    public:
      void add_tets(const matrixst & tets);
      std::vector<boost::unordered_set<size_t> > p2p_;
    };

    class one_ring_point_at_point
    {
    public:
      void build(const zjucad::matrix::matrix<size_t> & tet);
      void build(const std::vector<std::pair<size_t,size_t> > & edges);
      typedef boost::unordered_map<size_t, boost::unordered_set<size_t> > p2p_type;
      p2p_type p2p_;
    };

    int calculate_face_normal(const matrixst &tet,
                              const matrixd &tet_node,
                              const size_t &tet_idx,
                              const matrixst &face,
                              matrixd &normal);

    int calculate_face_normal_wrapper(const matrixst &tet,
                                      const matrixd &tet_node,
                                      const size_t &tet_idx,
                                      const size_t *face,
                                      matrixd &normal);

    int orient_one_face_normal_outside_tetmesh(
        const matrixst &tet,
        const matrixd &node,
        const matrixst &outside_face,
        const size_t &outside_face_idx,
        const jtf::mesh::face2tet_adjacent &fa,
        matrixd &face_normal);

    //! @brief orient the face normal to make it point outside,
    //  and will also adjust the three points order of each face
    int orient_face_normal_outside_tetmesh(const matrixst &tet,
                                           const matrixd &node,
                                           const matrixst &outside_face,
                                           const matrixst &outside_face_idx,
                                           const jtf::mesh::face2tet_adjacent &fa,
                                           matrixd &face_normal);

    int trans_face_normal_to_point_normal(const matrixd &face_normal,
                                          const matrixst &outside_face,
                                          matrixd & point_normal,
                                          matrixst & point_idx,
                                          std::map<size_t,size_t> &point_map,
                                          const matrixd *face_normal_weight = 0);

    int laplace_smoothing_triangle_mesh(const matrixst &tri_mesh,
                                        matrixd &node);

    //! @brief subdivide a prism into three tets, if subdividion failed, return non-zero
    // the way to subdivide a prism into tets is publised in:
    // "How to Subdivide Pyramids, Prisms and Hexahedra into Tetrahedra" [Julien Dompierre,etc]
    int subdivide_prism_into_tets(const matrixst &prism,
                                  matrixst &tets);

    //! @brief find a common face between tets, if can not find such a face, return non-zero
    // face should be 3*1
    int get_shared_face(const std::pair<size_t,size_t> &tet_pair,
                        const matrixst &tet,
                        size_t *face);

    // WARNING: experimental function which is used to extented a jtf::mesh::meshes with a layer of
    // jtf::mesh::meshes outside the original surface.
    // TODO: this function is not able to avoid intersection.
    int extend_tetmesh(const matrixst & tet,
                       const matrixd & node,
                       matrixst &new_tet,
                       matrixd &new_node,
                       std::map<size_t,size_t> & surface_point_map,
                       matrixd * zyz = 0);

    int extend_tetmesh(const zjucad::matrix::matrix<size_t> & tri_face,
                       const zjucad::matrix::matrix<double> & node,
                       const zjucad::matrix::matrix<double> & point_normal,
                       jtf::mesh::meshes & tet_layer);

    template <typename Iter>
    int check_manifold(Iter begin, Iter end,
                       const jtf::mesh::face2tet_adjacent & fa,
                       bool output_info = true)
    {
      boost::unordered_map<std::pair<size_t,size_t>, boost::unordered_set<size_t> >
          edge2faces;

      for(Iter i = begin; i != end; ++i){
          const size_t &face_idx = *i;
          if(face_idx > fa.faces_.size()){
              std::cerr << "# [error] wrong face index " << face_idx << std::endl;
              return __LINE__;
            }
          const std::vector<size_t> & face_vec = fa.faces_[face_idx];
          for(size_t i = 0; i < face_vec.size(); ++i){
              std::pair<size_t,size_t> edge(face_vec[i], face_vec[(i+1)%face_vec.size()]);
              if(edge.first > edge.second)
                std::swap(edge.first, edge.second);
              edge2faces[edge].insert(face_idx);
            }
        }
      for(boost::unordered_map<std::pair<size_t,size_t>,boost::unordered_set<size_t> >::const_iterator
          cit = edge2faces.begin(); cit != edge2faces.end(); ++cit){
          if(cit->second.size() > 2){
              if(output_info)
                std::cerr << "# [error] non-manifold: edge  " <<  cit->first.first
                          << " " << cit->first.second << " meet " << cit->second.size()
                          << " faces." << std::endl;
              return __LINE__;
            }
        }
      return 0;
    }

    class tetmesh_quality_improver
    {
    public:
      tetmesh_quality_improver(const matrixst& tet, matrixd & node)
        :tet_(tet), node_(node){
        init();
      }

      int improve(const std::string surface_strategy = "one_ring");
    private:
      int init();
      int improve_suface_by_one_ring();

      std::unique_ptr<jtf::mesh::face2tet_adjacent> fa_;
      const matrixst & tet_;
      matrixd & node_;
      matrixst outside_face_;
      matrixst outside_face_idx_;
    };

    ///
    /// \brief deform_tet_accordint_to_surface_constraints, this function apply
    ///        arap and surface point constraint to deform it.
    /// \param def_mesh
    /// \param def_node
    /// \param orig_mesh
    /// \param orig_node
    /// \param old_surface_point_to_new_point_map
    ///
    void deform_tet_accordint_to_surface_constraints(
        const zjucad::matrix::matrix<size_t> & def_mesh,
        zjucad::matrix::matrix<double> & def_node,
        const zjucad::matrix::matrix<size_t> & orig_mesh,
        const zjucad::matrix::matrix<double> & orig_node,
        const std::map<size_t,size_t> & old_surface_point_to_new_point_map,
        boost::property_tree::ptree & pt);

  }
}

#endif // JTF_TETMESH_UTIL_H

