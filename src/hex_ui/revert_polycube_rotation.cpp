#include <zjucad/ptree/ptree.h>
#include <jtflib/mesh/mesh.h>
#include <jtflib/mesh/io.h>
#include <hjlib/math/blas_lapack.h>
#include <zjucad/matrix/lapack.h>
#include <zjucad/matrix/io.h>
#include <hjlib/math/polar.h>
#include <functional>
#include "../common/zyz.h"
using namespace jtf::mesh;
using boost::property_tree::ptree;
using namespace std;
using namespace zjucad::matrix;

inline int calc_tet_def_grad_op(const double *tet, double *grad_op)
{
  using namespace zjucad::matrix;
  const static double node2edge[] = {
    1, 0, 0, -1,
    0, 1, 0, -1,
    0, 0, 1, -1
  };
  const static itr_matrix<const double *> E(4, 3, node2edge);

  itr_matrix<const double *> T(3, 4, tet);
  itr_matrix<double *> G(4, 3, grad_op);
  matrix<double> A = T*E;

  if(inv(A))
    return __LINE__;
  G = E*A;
  return 0;
}

int revert_polycube_rotation(ptree &pt)
{

  jtf::mesh::meshes orig_mesh, polycube_mesh;
  if(jtf::mesh::tet_mesh_read_from_zjumat(pt.get<string>("orig_tet.value").c_str(),
                                          &orig_mesh.node_, &orig_mesh.mesh_))
    return __LINE__;

  if(jtf::mesh::tet_mesh_read_from_zjumat(pt.get<string>("polycube_tet.value").c_str(),
                                          &polycube_mesh.node_, &polycube_mesh.mesh_))
    return __LINE__;

  assert(orig_mesh.mesh_.size(2) == polycube_mesh.mesh_.size(2));
  matrix<matrix<double> > frame(polycube_mesh.mesh_.size(2),1);
  matrix<double>  zyz(3, orig_mesh.mesh_.size(2));
  matrix<double> tet_node_orig(3,4), tet_node_poly(3,4), grad_op(4,3);
  for(size_t ti = 0; ti < orig_mesh.mesh_.size(2); ++ti){
      tet_node_orig = orig_mesh.node_(colon(), orig_mesh.mesh_(colon(),ti));
      calc_tet_def_grad_op(&tet_node_orig[0], &grad_op[0]);
      tet_node_poly = polycube_mesh.node_(colon(), polycube_mesh.mesh_(colon(),ti));
      frame[ti].resize(3,3);
      frame[ti] = tet_node_poly * grad_op;
      hj::polar3d p;
      p(frame[ti], 2);
      rotation_matrix_2_zyz_angle(&frame[ti][0], &zyz(0,ti),0);
    }

  const string zyz_str = pt.get<string>("zyz.value");
  if(write_matrix(zyz_str.c_str(), zyz))
    return __LINE__;

  return 0;
}
