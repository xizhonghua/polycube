#include "arap.h"
#include <zjucad/matrix/io.h>
#include <hjlib/math/blas_lapack.h>
#include <zjucad/matrix/lapack.h>
#include <jtflib/mesh/util.h>
#include <hjlib/math_func/func_aux.h>
#include <hjlib/math_func/operation.h>

using namespace std;
using namespace zjucad::matrix;

int calc_tet_def_grad_op(const double *tet, double *grad_op)
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

std::shared_ptr<hj::math_func::math_func_t<double,int32_t> >
build_arap_math_func(const zjucad::matrix::matrix<size_t> & tet,
                     const zjucad::matrix::matrix<double> & node,
                     const double weight)
{
  typedef hj::math_func::math_func_t<double, int32_t> math_func_type;
  typedef std::shared_ptr<math_func_type> math_func_ptr;
  shared_ptr<vector<math_func_ptr> > all_func(new vector<math_func_ptr>);

  vector<double> vol(tet.size(2));
  for(size_t ti = 0; ti < tet.size(2); ++ti){
      vol[ti] = fabs(jtf::mesh::cal_tet_vol(node(colon(), tet(colon(), ti))));
    }

  const double total_vol = std::accumulate(vol.begin(), vol.end(), 0.0);

  for(size_t ti = 0; ti < tet.size(2); ++ti)
    all_func->push_back(math_func_ptr(
                          new arap_math_func<double,int32_t>(tet, node, ti, sqrt(weight*vol[ti]/total_vol))));

  cerr << "# [info] arap function number " << all_func->size() << endl;
  math_func_ptr fun_cat(new hj::math_func::fcat<double,int32_t, vector<math_func_ptr> >(all_func));
  return fun_cat;
}
