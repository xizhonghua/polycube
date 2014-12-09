#include "least_square.h"

#include <iostream>

#include <hjlib/sparse/cached_AAT.h>
#include <zjucad/matrix/itr_matrix.h>
#include <hjlib/sparse/fast_AAT.h>

using namespace std;
using namespace zjucad::matrix;
using namespace hj::sparse;

least_square_Gauss_Newton::least_square_Gauss_Newton(const hj::function::function_t<double, int32_t> &f)
  :f_(f), is_JT_sorted_(false)
{
  r_.resize(f_.dim_of_f(), 1);
	JT_.resize(f_.dim_of_x(), f_.dim_of_f(), f_.jac_nnz());
	JT_.ptr()[0] = 0;
}

size_t least_square_Gauss_Newton::dim(void) const
{
  return f_.dim_of_x();
}

int least_square_Gauss_Newton::val(const double *x, double &v)
{
  if(f_.val(x, &r_[0]))
    return __LINE__;
  v += dot(r_, r_)/2;
  return 0;
}

int least_square_Gauss_Newton::gra(const double *x, double *g)
{
  assert(g);
  if(f_.val(x, &r_[0]))
    return __LINE__;
  if(f_.jac(x, &JT_.val()[0], &JT_.ptr()[0], &JT_.idx()[0]))
    return __LINE__;
  itr_matrix<double *> g0(f_.dim_of_x(), 1, g);
  mv(false, JT_, r_, g);
}

int least_square_Gauss_Newton::hes(const double *x, size_t &nnz, size_t &format, double *h, int32_t *ptr, int32_t *idx, double alpha)
{
  format = 1;
	if(hj::sparse::nnz(hes_))
		fast_AAT(JT_, hes_, is_JT_sorted_);
	else {
		is_JT_sorted_ = is_sorted_csc(JT_);
		AAT<map_by_sorted_vector>(JT_, hes_);
	}
  if(h == 0 && ptr == 0 && idx == 0) { // query nnz
    nnz = hj::sparse::nnz(hes_);
    return 0;
  }
  if(h == 0 && ptr !=0 && idx != 0) { // query patten
    if(nnz < hj::sparse::nnz(hes_)) {
      cerr << "nnz err: " << nnz << " " << hj::sparse::nnz(hes_) << endl;
      return __LINE__;
    }
    copy(hes_.ptr().begin(), hes_.ptr().end(), ptr);
    copy(hes_.idx().begin(), hes_.idx().end(), idx);
    return 0;
  }
  if(h != 0 && ptr !=0 && idx != 0) { // accumulate
    for(size_t ci = 0; ci < hes_.size(2); ++ci) {
      for(size_t nzi = hes_.ptr()[ci]; nzi < hes_.ptr()[ci+1]; ++nzi) {
        size_t off;
        for(off = ptr[ci]; off < ptr[ci+1]; ++off) {
          if(idx[off] == hes_.idx()[nzi])
            break;
        }
        if(off == ptr[ci+1]) {
          cerr << "incorrect sparse patten in least_square_Gauss_Newton: "
               << ci << " " << hes_.idx()[nzi] << endl;
          return __LINE__;
        }
        h[off] += hes_.val()[nzi]*alpha;
      }
    }
    return 0;
  }
  return __LINE__;
}
