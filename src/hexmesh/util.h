#ifndef HEXMESH_UTIL_H
#define HEXMESH_UTIL_H

#include <zjucad/matrix/matrix.h>
#include <jtflib/mesh/mesh.h>
#include "hexmesh.h"

namespace jtf{
  namespace hexmesh{
    typedef zjucad::matrix::matrix<size_t> matrixst;
    typedef zjucad::matrix::matrix<double> matrixd;

    class hex_singularity_extractor
    {
    public:
      hex_singularity_extractor(const jtf::hex_mesh &hm):hm_(hm){
        extractor();
      }

      const zjucad::matrix::matrix<size_t> get_all_singularity_edges()const{
        if(singularity_edges_.empty()) return zjucad::matrix::matrix<size_t>(0,0);
        zjucad::matrix::itr_matrix<const size_t*> sem(2, singularity_edges_.size()/2, &singularity_edges_[0]);
        return sem;
      }
      const zjucad::matrix::matrix<size_t> get_surface_singularity_edges()const{
        if(singularity_edges_.empty() || surface_edges_idx_.empty())
          return zjucad::matrix::matrix<size_t>(0,0);
        zjucad::matrix::itr_matrix<const size_t*> sem(2, singularity_edges_.size()/2, &singularity_edges_[0]);
        zjucad::matrix::itr_matrix<const size_t*> surface_sem(surface_edges_idx_.size(),1, &surface_edges_idx_[0]);
        return sem(zjucad::matrix::colon(), surface_sem);
      }
      const zjucad::matrix::matrix<size_t> get_inner_singularity_edges()const{
        if(singularity_edges_.empty() || inner_edges_idx_.empty())
          return zjucad::matrix::matrix<size_t>(0,0);
        zjucad::matrix::itr_matrix<const size_t*> sem(2, singularity_edges_.size()/2, &singularity_edges_[0]);
        zjucad::matrix::itr_matrix<const size_t*> inner_sem(inner_edges_idx_.size(),1, &inner_edges_idx_[0]);
        return sem(zjucad::matrix::colon(), inner_sem);
      }

      const zjucad::matrix::matrix<size_t> get_all_singularity_edges_degree()const{
        if(singularity_edges_.empty()) return zjucad::matrix::matrix<size_t>(0,0);
        zjucad::matrix::itr_matrix<const size_t*> dm(edge_degree_.size(),1, &edge_degree_[0]);
        return dm;
      }
      const zjucad::matrix::matrix<size_t> get_surface_singularity_edges_degree()const{
        if(singularity_edges_.empty() || surface_edges_idx_.empty())
          return zjucad::matrix::matrix<size_t>(0,0);
        zjucad::matrix::itr_matrix<const size_t*> dm(edge_degree_.size(),1,&edge_degree_[0]);
        zjucad::matrix::itr_matrix<const size_t*> surface_dm(surface_edges_idx_.size(),1, &surface_edges_idx_[0]);
        return dm(surface_dm, 0);
      }
      const zjucad::matrix::matrix<size_t> get_inner_singularity_edges_degree()const{
        if(singularity_edges_.empty() || inner_edges_idx_.empty())
          return zjucad::matrix::matrix<size_t>(0,0);
        zjucad::matrix::itr_matrix<const size_t*> dm(edge_degree_.size(),1,&edge_degree_[0]);
        zjucad::matrix::itr_matrix<const size_t*> inner_dm(inner_edges_idx_.size(),1, &inner_edges_idx_[0]);
        return dm(inner_dm, 0);
      }

    protected:
      void extractor(){
        clear();
        for(const auto & one_edge : hm_.orhae_.e2h_){
            const std::pair<size_t,size_t> & edge = one_edge.first;
            const std::vector<size_t> & around_hexs = one_edge.second;
            if(around_hexs.front() == -1 || around_hexs.back() == -1){//surface edge
                size_t degree = around_hexs.size();
                if(around_hexs.front() == -1) --degree;
                if(around_hexs.back() == -1) --degree;
                if(degree != 2){
                    singularity_edges_.push_back(edge.first);
                    singularity_edges_.push_back(edge.second);
                    surface_edges_idx_.push_back(singularity_edges_.size()/2-1);
                    edge_degree_.push_back(degree);
                  }
              }else{
                size_t degree = around_hexs.size()-1;
                if(degree != 4){
                    singularity_edges_.push_back(edge.first);
                    singularity_edges_.push_back(edge.second);
                    inner_edges_idx_.push_back(singularity_edges_.size()/2-1);
                    edge_degree_.push_back(degree);
                  }
              }
          }
      }
      void clear(){
        singularity_edges_.clear();
        inner_edges_idx_.clear();
        surface_edges_idx_.clear();
      }
    private:
      const jtf::hex_mesh & hm_;
      std::vector<size_t> singularity_edges_;
      std::vector<size_t> edge_degree_;
      std::vector<size_t> inner_edges_idx_;
      std::vector<size_t> surface_edges_idx_;
    };

    int subdivide_hexmesh(matrixst & hex, matrixd & node);

    ///
    /// @brief orient_hex
    /// @param hex
    /// @param node
    ///
    void orient_hex(zjucad::matrix::matrix<size_t> & hex,
                    const zjucad::matrix::matrix<double> & node);

    class hexmesh_quality_improver
    {
    public:
      hexmesh_quality_improver(const zjucad::matrix::matrix<size_t> & hex,
                               zjucad::matrix::matrix<double> & node)
        :hex_(hex), node_(node){
        init();
      }

      int improve(const std::string strategy = "shape");
    private:
      int init();
      int improve_suface_by_one_ring();

      std::unique_ptr<jtf::mesh::face2hex_adjacent> fa_;
      const zjucad::matrix::matrix<size_t> & hex_;
      zjucad::matrix::matrix<double> & node_;
      zjucad::matrix::matrix<size_t> outside_face_;
      zjucad::matrix::matrix<size_t> outside_face_idx_;
    };


		void extend_hexmesh_based_on_faces(const zjucad::matrix::matrix<size_t> & faces,
																		const zjucad::matrix::matrix<double> & normal,
																		const jtf::mesh::meshes & hexmesh,
																		jtf::mesh::meshes & new_hexmesh,
																		const double len_percent);
  }
}

#endif // HEXMESH_UTIL_H
