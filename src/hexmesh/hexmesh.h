#ifndef HEXMESH_H
#define HEXMESH_H

#include <jtflib/mesh/mesh.h>
#include <jtflib/mesh/util.h>
#include "io.h"
namespace jtf{
  class hex_mesh
  { public:
    hex_mesh(const char *file):hex_format_(1){
      if(init(file))
        throw std::invalid_argument("wrong hex mesh.");
    }

    hex_mesh():hex_format_(1){}

    hex_mesh(const jtf::mesh::meshes &hm_):hex_format_(1){
      if(init(hm_))
        throw std::invalid_argument("wrong hex mesh.");
    }

    hex_mesh(const hex_mesh &hm):hex_format_(1){
      if(init(hm))
        throw std::invalid_argument("wrong hex mesh.");
    }
    hex_mesh & operator = (const hex_mesh &hm)
    {
      if(&hm == this) return *this;
      if(init(hm))
        throw std::invalid_argument("wrong hex mesh.");
      return *this;
    }
    int load(const char * file){return init(file);}
    int load(const jtf::mesh::meshes &hm){return init(hm);}
    //////////////////////////////////////////////////////////////////////////////
    jtf::mesh::meshes hexmesh_;
    //zjucad::matrix::matrix<double> vol_;
    zjucad::matrix::matrix<size_t> outside_face_;
    zjucad::matrix::matrix<size_t> outside_face_idx_;
    zjucad::matrix::matrix<double> outside_face_normal_;
    zjucad::matrix::matrix<double> outside_face_area_;
    std::shared_ptr<jtf::mesh::face2hex_adjacent> fa_;
    std::shared_ptr<jtf::mesh::edge2cell_adjacent> ea_outside_;
    jtf::mesh::one_ring_hex_at_edge orhae_;
    size_t hex_format_;
  private:
    int init(const char * file){
      using namespace zjucad::matrix;
      using namespace std;
      if(jtf::hexmesh::hex_mesh_read_from_wyz(file, hexmesh_.mesh_, hexmesh_.node_, hex_format_)){
          cerr << "# [error] can not load hexmesh." << endl;
          return 1;
        }
      return init(hexmesh_);
    }
    int init(const jtf::mesh::meshes &hm){
      using namespace zjucad::matrix;
      using namespace std;
      if(&hm != & hexmesh_)
        hexmesh_ = hm;
      fa_.reset(jtf::mesh::face2hex_adjacent::create(hexmesh_.mesh_));
      if(!fa_.get()){
          cerr << "# [error] invalide hexmesh." << endl;
          return __LINE__;
        }
      jtf::mesh::get_outside_face(*fa_, outside_face_);
      jtf::mesh::get_outside_face_idx(*fa_, outside_face_idx_);
      jtf::mesh::cal_face_normal(outside_face_,hexmesh_.node_,outside_face_normal_, true);
      outside_face_area_.resize(outside_face_.size(2),1);
      for(size_t fi = 0; fi < outside_face_.size(2); ++fi){
          outside_face_area_[fi] = jtf::mesh::cal_face_area(outside_face_(colon(),fi),hexmesh_.node_);
        }
      ea_outside_.reset(jtf::mesh::edge2cell_adjacent::create(outside_face_));
      if(!ea_outside_.get()){
          cerr << "# [error] invalide outside_face." << endl;
          return __LINE__;
        }
      orhae_.e2h_.clear();
      orhae_.add_hex(hexmesh_.mesh_, hexmesh_.node_, *fa_);
      orhae_.sort_into_loop(hexmesh_.mesh_, hexmesh_.node_, *fa_);

      return 0;
    }
    int init(const jtf::hex_mesh &hm){
      hexmesh_ = hm.hexmesh_;
      outside_face_ = hm.outside_face_;
      outside_face_idx_ = hm.outside_face_idx_;
      outside_face_normal_ = hm.outside_face_normal_;
      outside_face_area_ = hm.outside_face_area_;
      fa_.reset(jtf::mesh::face2hex_adjacent::create(hexmesh_.mesh_));
      if(!fa_.get()){
          std::cerr << "# [error] invalide tetmesh." << std::endl;
          return __LINE__;
        }
      ea_outside_.reset(jtf::mesh::edge2cell_adjacent::create(outside_face_));
      if(!ea_outside_.get()){
          std::cerr << "# [error] invalide outside_face." << std::endl;
          return __LINE__;
        }
      orhae_.e2h_.clear();
      orhae_.add_hex(hexmesh_.mesh_, hexmesh_.node_, *fa_);
      orhae_.sort_into_loop(hexmesh_.mesh_, hexmesh_.node_, *fa_);
      return 0;
    }
  };
}
#endif // HEXMESH_H
