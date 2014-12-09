#ifndef HJ_TET_MESH_H_
#define HJ_TET_MESH_H_
#include <utility>
#include <vector>
#include <map>
#include <iostream>
#include <memory>
#include <boost/unordered_map.hpp>
#include <zjucad/matrix/itr_matrix.h>
#include <zjucad/matrix/io.h>
#include <jtflib/mesh/mesh.h>
#include <jtflib/mesh/io.h>
#include <jtflib/mesh/util.h>

#include "../common/def.h"

//! @param type 0: undirectional, 1:increase, -1:decrease
int tet2tet_adj_matrix(
    const matrixst &tet,
    matrixst &ptr,
    matrixst &idx,
    const jtf::mesh::face2tet_adjacent *fa = 0,
    int type = 0
    );

int face2face_adj_matrix(
    const matrixst &tet,
    matrixst &ptr,
    matrixst &idx,
    const jtf::mesh::face2tet_adjacent *fa = 0,
    int type = 0
    );

template<typename T1, typename T2>
int orient_tet_raw(const double * node_array, T2 node_num,
                   T1 * tet_array, T2 tet_num)
{
  zjucad::matrix::itr_matrix<const double *> node(3, node_num, node_array);
  zjucad::matrix::itr_matrix<T1 *> tet(4, tet_num, tet_array);

  for(T2 ti = 0; ti < tet.size(2); ++ti) {
      matrixd ele(3, 3);
      for(T2 ni = 0; ni < 3; ++ni)
        ele(zjucad::matrix::colon(), ni) =
            node(zjucad::matrix::colon(), tet(ni+1, ti))
            - node(zjucad::matrix::colon(), tet(0, ti));

      if(zjucad::matrix::dot(
           zjucad::matrix::cross(ele(zjucad::matrix::colon(), 0),
                                 ele(zjucad::matrix::colon(), 1)),
           ele(zjucad::matrix::colon(), 2)) < 0) {
          std::swap(tet(1, ti), tet(2, ti));
        }
    }
  return 0;
}

//! @brief <(p1-p0)x(p2-p0),(p3-p0)> >= 0
template <typename T>
int orient_tet(const matrixd &node,
               zjucad::matrix::matrix<T> &tet)
{
  return orient_tet_raw(&node[0], node.size(2), &tet[0], tet.size(2));
}

namespace jtf{
  class tet_mesh
  {
  public:
    tet_mesh(const char *file){
      if(init(file))
        throw std::invalid_argument("wrong tet mesh.");
    }

    tet_mesh(){}

    tet_mesh(const jtf::mesh::meshes &tm_){
      if(init(tm_))
        throw std::invalid_argument("wrong tet mesh.");
    }

    tet_mesh(const tet_mesh &tm){
      if(init(tm))
        throw std::invalid_argument("wrong tet mesh.");
    }
    tet_mesh (const zjucad::matrix::matrix<size_t> & mesh,
              const zjucad::matrix::matrix<double> & node){
      if(init(mesh, node)){
          throw std::invalid_argument("wrong tet mesh.");
        }
    }
    tet_mesh & operator = (const tet_mesh &tm)
    {
      if(&tm == this) return *this;
      if(init(tm))
        throw std::invalid_argument("wrong tet mesh.");
      return *this;
    }
    int load(const char * file){return init(file);}
    int load(const jtf::mesh::meshes &tm){return init(tm);}
    int load(const zjucad::matrix::matrix<size_t> & mesh,
             const zjucad::matrix::matrix<double> & node){
      return init(mesh, node);
    }

    //////////////////////////////////////////////////////////////////////////////
    jtf::mesh::meshes tetmesh_;
    zjucad::matrix::matrix<double> vol_;
    zjucad::matrix::matrix<size_t> outside_face_;
    zjucad::matrix::matrix<size_t> outside_face_idx_;
    zjucad::matrix::matrix<double> outside_face_normal_;
    zjucad::matrix::matrix<double> outside_face_area_;
    std::shared_ptr<jtf::mesh::face2tet_adjacent> fa_;
    std::shared_ptr<jtf::mesh::edge2cell_adjacent> ea_outside_;
    jtf::mesh::one_ring_tet_at_edge ortae_;

  private:
    int init(const char * file){
      using namespace std;
      if(jtf::mesh::tet_mesh_read_from_zjumat(file, &tetmesh_.node_, &tetmesh_.mesh_)){
          cerr << "# [error] can not load tetmesh." << endl;
          return 1;
        }
      return init(tetmesh_);
    }
    int init(const zjucad::matrix::matrix<size_t> & mesh,
             const zjucad::matrix::matrix<double> & node){
     tetmesh_.mesh_ = mesh;
     tetmesh_.node_ = node;
     return init(tetmesh_);
    }
    int init(const jtf::mesh::meshes &tm){
      using namespace zjucad::matrix;
      using namespace std;
      if(&tm != &tetmesh_)
        tetmesh_ = tm;
      vol_.resize(tetmesh_.mesh_.size(2),1);
      for(size_t ti = 0; ti < tetmesh_.mesh_.size(2); ++ti){
          vol_[ti] = jtf::mesh::cal_tet_vol(tetmesh_.node_(colon(), tetmesh_.mesh_(colon(),ti)));
        }
      fa_.reset(jtf::mesh::face2tet_adjacent::create(tetmesh_.mesh_));
      if(!fa_.get()){
          cerr << "# [error] invalide tetmesh." << endl;
          return __LINE__;
        }
      jtf::mesh::get_outside_face(*fa_, outside_face_,true, &tetmesh_.node_);
      jtf::mesh::get_outside_face_idx(*fa_, outside_face_idx_);
      jtf::mesh::cal_face_normal(outside_face_,tetmesh_.node_,outside_face_normal_, true);
      outside_face_area_.resize(outside_face_.size(2),1);
      for(size_t fi = 0; fi < outside_face_.size(2); ++fi){
          outside_face_area_[fi] = jtf::mesh::cal_face_area(outside_face_(colon(),fi),tetmesh_.node_);
        }
      ea_outside_.reset(jtf::mesh::edge2cell_adjacent::create(outside_face_));
      if(!ea_outside_.get()){
          cerr << "# [error] invalide outside_face." << endl;
          return __LINE__;
        }
      ortae_.e2t_.clear();
      ortae_.add_tets(tetmesh_.mesh_, *fa_);
      ortae_.sort_into_loop(tetmesh_.mesh_, tetmesh_.node_);
      return 0;
    }
    int init(const jtf::tet_mesh &tm){
      tetmesh_ = tm.tetmesh_;
      vol_ = tm.vol_;
      outside_face_ = tm.outside_face_;
      outside_face_idx_ = tm.outside_face_idx_;
      outside_face_normal_ = tm.outside_face_normal_;
      outside_face_area_ = tm.outside_face_area_;
      ortae_ = tm.ortae_;
      fa_.reset(jtf::mesh::face2tet_adjacent::create(tetmesh_.mesh_));
      if(!fa_.get()){
          std::cerr << "# [error] invalide tetmesh." << std::endl;
          return __LINE__;
        }
      ea_outside_.reset(jtf::mesh::edge2cell_adjacent::create(outside_face_));
      if(!ea_outside_.get()){
          std::cerr << "# [error] invalide outside_face." << std::endl;
          return __LINE__;
        }
      return 0;
    }
  };
}

#endif
