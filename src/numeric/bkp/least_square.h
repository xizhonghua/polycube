#ifndef HJ_LEAST_SQUARE_H_
#define HJ_LEAST_SQUARE_H_

#include <hjlib/function/function.h>

#include "func_opt.h"

#include <zjucad/matrix/matrix.h>
#include <hjlib/sparse/sparse.h>

//! assume the dim and pattern of f is constant
class least_square_Gauss_Newton : public hj_func_opt::function
{
public:
  least_square_Gauss_Newton(const hj::function::function_t<double, int32_t> &f);
  virtual size_t dim(void) const;
  virtual int val(const double *x, double &v);

  //! @param nnz=0 means dense, g=0 or h=0 means query nnz
  virtual int gra(const double *x, double *g);
  virtual int hes(const double *x, size_t &nnz, size_t &format, double *h, int32_t *ptr, int32_t *idx, double alpha = 1);
  virtual int hes_block(const double *x, double *h, double alpha = 1) {
    return -1;
  }
protected:
  const hj::function::function_t<double, int32_t> &f_;
  zjucad::matrix::matrix<double> r_;
  hj::sparse::csc<double, int32_t> JT_, hes_;
  bool is_JT_sorted_;
};


#endif
