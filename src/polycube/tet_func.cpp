#include "tet_func.h"
#include "../numeric/util.h"
#include "util.h"
#include <numeric>
#include <iostream>
#include <vector>

#include <hjlib/math/polar.h>
#include <hjlib/math/blas_lapack.h>

#include <zjucad/matrix/itr_matrix.h>
#include <zjucad/matrix/io.h>
#include <zjucad/matrix/lapack.h>

#include <hjlib/function/func_aux.h>
#include <hjlib/sparse/sparse.h>
#include <jtflib/math/math.h>


using namespace std;
using namespace zjucad::matrix;
using namespace hj::function;

int calc_tet_def_grad_op(const double *tet, double *grad_op)
{
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

int calc_tet_arap_distortion(const double *orig_tet, const double * new_tet,
                             double *dis)
{
  matrix<double> grad_op(4,3);
  itr_matrix<const double *> orig_T(3, 4, orig_tet);
  itr_matrix<const double *> new_T(3, 4, new_tet);
  if(calc_tet_def_grad_op(&orig_T[0], &grad_op[0])) {
    cerr << "# degenerated tet" << __LINE__ << orig_T << endl;
  }
  itr_matrix<double*> f(3,3,dis);

  f = new_T * grad_op;
  matrix<double> R = f;
  hj::polar3d p;
  p(R,2);
  f -= R;
  return 0;
}

class tet_grad_func : public hj::function::function_t<double, int32_t>
{
public:
  tet_grad_func(const matrix<double> &node, const matrix<size_t> &tet, size_t idx)
    :grad_op_(4, 3), tet_(tet), node_num_(node.size(2)), idx_(idx) {
    matrix<double> tet_node = node(colon(), tet_(colon(), idx_));
    if(calc_tet_def_grad_op(&tet_node[0], &grad_op_[0]))
      cerr << "# degenerated tet" << __LINE__ << tet_node << endl;
  }
  virtual size_t dim_of_x(void) const {
    return node_num_*3;
  }
  virtual size_t dim_of_f(void) const {
    return 9;
  }
  virtual int val(const double *x, double *f, func_ctx *ctx = 0) const {
    itr_matrix<double *> f0(3, 3, f);
    const itr_matrix<const double *> T(3, node_num_, x);
    f0 = T(colon(), tet_(colon(), idx_))*grad_op_;
    for(size_t i = 0; i < f0.size(); ++i) jtf::math::erase_nan_inf(f0[i]);
    return 0;
  }
  virtual int jac(const double *x, double *val, int32_t *ptr = 0,
                  int32_t *idx = 0, func_ctx *ctx = 0) const {
    for(int c = 0; c < 3; ++c) {
      for(int r = 0; r < 3; ++r) {
        const int fi = c*3+r;
        ptr[fi+1] = ptr[fi] + 4;
        for(int k = 0; k < 4; ++k) {
          idx[ptr[fi]+k] = tet_(k, idx_)*3+r;
          val[ptr[fi]+k] = grad_op_(k, c);
          jtf::math::erase_nan_inf(val[ptr[fi]+k]);
        }
      }
    }
    return 0;
  }
  virtual size_t jac_nnz(void) const {
    return 4*dim_of_f();
  }
private:
  const size_t idx_, node_num_;
  const matrix<size_t> &tet_;
  matrix<double> grad_op_;
};

function_t<double, int32_t> *
build_tetmesh_grad_func(const zjucad::matrix::matrix<double> &node,
                        const zjucad::matrix::matrix<size_t> &tet)
{
  const size_t tn = tet.size(2);

  boost::shared_ptr<vector<boost::shared_ptr<function_t<double, int32_t> > > >
      funcs(new vector<boost::shared_ptr<function_t<double, int32_t> > >(tn));
  for(size_t ti = 0; ti < tn; ++ti) {
    (*funcs)[ti].reset(new tet_grad_func(node, tet, ti));
#if 0 // validate jac
    {
      const double err = jac_err(*(*funcs)[ti], &node[0]);
      if(err > 1e-7) {
        cerr << "large error in jac: " << err << endl;
      }
    }
#endif
  }
  return new_catenated_function<double, int32_t>(funcs);
}

class tetmesh_arap_with_weight_func :
    public hj::function::function_t<double, int32_t>
{
public:
  tetmesh_arap_with_weight_func(const matrix<double> &node,
                                const matrix<size_t> &tet,
                                const size_t &idx,
                                const double &weight)
    :grad_op_(4, 3), tet_(tet), node_num_(node.size(2)), idx_(idx) , weight_(weight){
    matrix<double> tet_node = node(colon(), tet_(colon(), idx_));
    if(calc_tet_def_grad_op(&tet_node[0], &grad_op_[0])) {
      cerr << "# degenerated tet" << __LINE__ << tet_node << endl;
    }
  }
  virtual size_t dim_of_x(void) const {
    return node_num_*3;
  }
  virtual size_t dim_of_f(void) const {
    return 9;
  }
  virtual int val(const double *x, double *f, func_ctx *ctx = 0) const {
    itr_matrix<double *> f0(3, 3, f);
    const itr_matrix<const double *> T(3, node_num_, x);
    f0 = T(colon(), tet_(colon(), idx_))*grad_op_;
	
    matrix<double> R = f0;
    hj::polar3d p;
    p(R, 2);
    const double alpha = 1;
    R = alpha*R+(1-alpha)*eye<double>(3);
    f0 -= R;
//	if(idx_ == 4640){
//	  cerr << "# x" << T(colon(), tet_(colon(), idx_)) << endl;
//	  cerr << "# f0 " << f0 << endl;
//	  cerr << "# R " << R << endl;
//}
    f0 *= weight_;

    for(size_t i = 0; i < f0.size(); ++i) jtf::math::erase_nan_inf(f0[i]);
    return 0;
  }
  virtual int jac(const double *x, double *val, int32_t *ptr = 0,
                  int32_t *idx = 0, func_ctx *ctx = 0) const {
    for(int c = 0; c < 3; ++c) {
      for(int r = 0; r < 3; ++r) {
        const int fi = c*3+r;
        ptr[fi+1] = ptr[fi] + 4;
        for(int k = 0; k < 4; ++k) {
          idx[ptr[fi]+k] = tet_(k, idx_)*3+r;
          val[ptr[fi]+k] = grad_op_(k, c) * weight_;
          jtf::math::erase_nan_inf(val[ptr[fi]+k]);
        }
      }
    }
    return 0;
  }
  virtual size_t jac_nnz(void) const {
    return 4*dim_of_f();
  }
private:
  const size_t idx_, node_num_;
  const matrix<size_t> &tet_;
  matrix<double> grad_op_;
  const double weight_;
};

int nnz(const zjucad::matrix::matrix<double> & m){
  int a = 0;
  for(size_t i = 0; i < m.size(); ++i){
    if(fabs(m[i]) > 1e-6)
      ++a;
  }
  return a;
}


class tetmesh_arap_func : public hj::function::function_t<double, int32_t>
{
public:
  tetmesh_arap_func(const matrix<double> &node, const matrix<size_t> &tet)
    :grad_(build_tetmesh_grad_func(node, tet)), tn_(tet.size(2)){
  }
  virtual size_t dim_of_x(void) const {
    return grad_->dim_of_x();
  }
  virtual size_t dim_of_f(void) const {
    return grad_->dim_of_f();
  }
  virtual int val(const double *x, double *f, func_ctx *ctx = 0) const {
    grad_->val(x, f, ctx);
    for(size_t ti = 0; ti < tn_; ++ti, f+=9) {
      itr_matrix<double *> f0(3, 3, f);
      matrix<double> R = f0;
      hj::polar3d p;
      p(R, 2);
      const double alpha = 1;
      R = alpha*R+(1-alpha)*eye<double>(3);
      f0 -= R;
      for(size_t i = 0; i < f0.size(); ++i) jtf::math::erase_nan_inf(f0[i]);
    }
    return 0;
  }
  virtual int jac(const double *x, double *val, int32_t *ptr = 0, int32_t *idx = 0, func_ctx *ctx = 0) const {
    return grad_->jac(x, val, ptr, idx, ctx);
  }
  virtual size_t jac_nnz(void) const {
    return grad_->jac_nnz();
  }
private:
  const size_t tn_;
  unique_ptr<function_t<double, int32_t> > grad_;
};


hj::function::function_t<double, int32_t> *
build_tetmesh_arap_func(const zjucad::matrix::matrix<double> &node,
                        const zjucad::matrix::matrix<size_t> &tet,
                        const double arap_w)
{
  matrix<double> v_weight_grad = ones<double>(tet.size(2),1);

  matrix<double> one_tet_node = zeros<double>(3,4);
  for(size_t ti = 0; ti < tet.size(2); ++ti){
    one_tet_node = node(colon(), tet(colon(), ti));
    v_weight_grad[ti] = jtf::mesh::cal_tet_vol(one_tet_node);
    if(v_weight_grad[ti] < 0){
      v_weight_grad[ti] *= -1; // strange, why should v * -1?
    }
  }

  const double total_v =
      std::accumulate(v_weight_grad.begin(), v_weight_grad.end(), 0.0);

  for(size_t ti = 0; ti < tet.size(2); ++ti){
    v_weight_grad[ti] = sqrt(arap_w * v_weight_grad[ti] / total_v);
  }

  const size_t tn = tet.size(2);
  boost::shared_ptr<vector<boost::shared_ptr<function_t<double, int32_t> > > >
      funcs(new vector<boost::shared_ptr<function_t<double, int32_t> > >(tn));
  for(size_t ti = 0; ti < tn; ++ti) {
    (*funcs)[ti].reset(
          new tetmesh_arap_with_weight_func( node, tet, ti, v_weight_grad[ti]));
  }
  return new_catenated_function<double, int32_t>(funcs);

  //return new tetmesh_arap_func(node, tet, weight);
}
