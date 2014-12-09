#ifndef CONNECTION_GENERATOR_H
#define CONNECTION_GENERATOR_H

#include <jtflib/mesh/mesh.h>
#include "tet_mesh_fxz.h"
#include "../tetmesh/tetmesh.h"

namespace fxz {
  typedef zjucad::matrix::matrix<double> matrixd;
  typedef zjucad::matrix::matrix<size_t> matrixst;

  class connection_generator
  {
 public:
    connection_generator(jtf::tet_mesh& tm,
                         const matrixd& size_field,
                         std::map<std::pair<size_t,size_t>,matrixd>& faces_conn)
        : tm_(tm), size_field_(size_field), faces_conn_(faces_conn) { }
    ~connection_generator() { }

    int run();
    
 protected:

    int cal_size_grad();
    int cal_tets_inv_grad();
    int cal_faces_connection();

 private:
    double cal_point_size(size_t tid, const matrixd& p);
    int cal_W(size_t tid, const matrixd& p0, const matrixd& p1, matrixd& conn);

 private:
    connection_generator(const connection_generator&) = delete;
    connection_generator& operator=(const connection_generator&) = delete;

 private:
    jtf::tet_mesh& tm_;
    const matrixd& size_field_;
    std::map<std::pair<size_t,size_t>, matrixd>& faces_conn_;
    
    std::shared_ptr<tet_mesh> mesh_ptr_;
    std::vector<matrixd> grad_inv_; // 3*3 matrixd vector
    std::vector<matrixd> size_grad_; // 1*3 matrixd vector
  };


} // namespace fxz

#endif
