#ifndef FUNCTION_TERM_H
#define FUNCTION_TERM_H

#include <string.h>

#include <iostream>
#include <fstream>
#include <algorithm>
#include <numeric>

#include <omp.h>

#include "../common/zyz.h"
#include "../common/IO.h"
#include "../common/util.h"
#include "../common/vtk.h"

#include "../spherical_harmonics/rot_cubic_f_SH.h"

#include <hjlib/sparse/sparse.h>
//#include <hjlib/arg_opts/arg_opts.h>
#include <hjlib/function/function.h>
//#include <hjlib/optimizer/optimizer.h>
#include <zjucad/optimizer/optimizer.h>

#include <hjlib/math/blas_lapack.h>
#include <jtflib/function/function.h>
#include <jtflib/function/func_aux.h>
#include <zjucad/matrix/io.h>
#include <zjucad/matrix/itr_matrix.h>
#include <zjucad/matrix/lapack.h>

#include "../tetmesh/tetmesh.h"
#include "../tetmesh/hex_io.h"

#include "../hex_frame_opt/hex_frame_opt.h"
#include "../hex_param/hex_param.h"
#include "../hex_param/common.h"
#include "../numeric/util.h"
#include <boost/property_tree/ptree.hpp>
#include <zjucad/ptree/ptree.h>

#include <minpack.h>
#include "Ax.h"
using namespace std;
using namespace zjucad::matrix;
using namespace hj::sparse;
using namespace hj::function;

using boost::property_tree::ptree;

class s2i_function : public function_t<double, int32_t>
{
public:
  s2i_function(const jtf::mesh::face2tet_adjacent &fa)
    :num_x_(fa.faces_.size()*9) {
  }
  virtual ~s2i_function(){}

  virtual size_t dim_of_x(void) const {
    return num_x_;
  }
protected:
  const size_t num_x_;
};

class s2i_smooth : public s2i_function
{
public:
  s2i_smooth(const jtf::mesh::face2tet_adjacent &fa, const matrixst &tet, const matrixd &node, double weight)
    :s2i_function(fa) {
    matrixd P = ones<double>(4, 4);
    node_idx_.resize(4);
    for(int i = 0; i < 4; ++i) {
        node_idx_[i] = fa.get_face_idx(tet[i], tet[(i+1)%4], tet[(i+2)%4]);
        if(node_idx_[i] >= fa.faces_.size())
          cerr << "wrong face index in fa." << endl;
        matrixd face_center = zeros<double>(3, 1);
        for(int j = 0; j < 3; ++j)
          face_center += node(colon(), tet[(j+i)%4]);
        P(colon(0, 2), i) = face_center/3.0;
      }
    if(inv(P))
      cerr << "inverse fail." << endl;
    M_ = P(colon(), colon(0, 2))*weight;
  }

  virtual size_t dim_of_f(void) const {
    return 9*3;
  }
  virtual int val(const double *x, double *f, func_ctx *ctx = 0) const {
    const itr_matrix<const double *> x0(9, num_x_/9, x);
    itr_matrix<double *> f0(9, 3, f);
    const matrixd x1 = x0(colon(), node_idx_);
    f0 = x1*M_;
    return 0;
  }
  virtual int jac(const double *x, double *val, int32_t *ptr = 0, int32_t *idx = 0, func_ctx *ctx = 0) const {
    for(int32_t d = 0, i = 0; d < 3; ++d) {
        for(int32_t shi = 0; shi < 9; ++shi, ++i) { // for each function component
            ptr[i+1] = ptr[i] + 4;
            for(int ni = 0; ni < 4; ++ni) { // for each face node
                idx[ptr[i]+ni] = node_idx_[ni]*9+shi;
                val[ptr[i]+ni] = M_(ni, d);
              }
          }
      }
    return 0;
  }
  virtual size_t jac_nnz(void) const {
    return 4*dim_of_f();
  }
protected:
  matrixd M_;
  matrixst node_idx_;
};

class s2i_fix : public s2i_function
{
public:
  s2i_fix(const jtf::mesh::face2tet_adjacent &fa, size_t node_idx, const matrixd &sh, double weight)
    :s2i_function(fa), node_idx_(node_idx), sh_(sh), weight_(weight) {
  }

  virtual size_t dim_of_f(void) const {
    return 9;
  }
  virtual int val(const double *x, double *f, func_ctx *ctx = 0) const {
    const itr_matrix<const double *> x0(9, num_x_/9, x);
    itr_matrix<double *> f0(9, 1, f);
    const matrixd x1 = x0(colon(), node_idx_);
    f0 = (x1-sh_)*weight_;
    return 0;
  }
  virtual int jac(const double *x, double *val, int32_t *ptr = 0, int32_t *idx = 0, func_ctx *ctx = 0) const {
    for(size_t i = 0; i < 9; ++i) {
        ptr[i+1] = ptr[i]+1;
        idx[ptr[i]] = node_idx_*9+i;
        val[ptr[i]] = weight_;
      }
    return 0;
  }
  virtual size_t jac_nnz(void) const {
    return dim_of_f();
  }
protected:
  const matrixd sh_;
  const size_t node_idx_;
  const double weight_;
};

class inner_smooth : public function_t<double, int32_t>
{
public:
  inner_smooth(const size_t node_num,
               const size_t adjacent_tet_i,
               const size_t adjacent_tet_j,
               const double weight)
    :num_x(node_num),  weight_(weight)//,M_[0](1), M_[1](-1)
  {
    v[0] = adjacent_tet_i;
    v[1] = adjacent_tet_j;
    M_[0] = 1;
    M_[1] = -1;
  }
  virtual size_t dim_of_x() const
  {
    return num_x * 9;
  }
  virtual size_t dim_of_f() const
  {
    return 9 * 1;
  }
  virtual int val(const double *x, double *f, func_ctx *ctx = 0) const
  {
    itr_matrix<const double*> x_value(9,num_x,x);
    itr_matrix<double*> f0(9,1,f);
    f0 = x_value(colon(),v[0])  - x_value(colon(),v[1]) ;
    f0 *= weight_;
    return 0;
  }
  virtual int jac(const double *x, double *val, int32_t *ptr=0, int32_t *idx = 0, func_ctx *ctx = 0) const
  {
    itr_matrix<const double*> x0(9,num_x,x);
    for(int32_t shi = 0, i = 0; shi < 9; ++i, ++shi)
      {
        ptr[i+1] = ptr[i] + 2;
        for(int nd = 0; nd < 2; ++nd)
          {
            idx[ptr[i] + nd] = v[nd] * 9 + shi;
            val[ptr[i] + nd] = weight_ * M_[nd];
          }
      }
    return 0;
  }
  virtual size_t jac_nnz(void) const {
    return 2 * dim_of_f();
  }
protected:
  size_t num_x;
  size_t v[2];
  double M_[2];
  const double weight_;
};

class align_sh_inner : public function_t<double, int32_t>
{
public:
  align_sh_inner(const size_t tet_num,
                 size_t node_idx,
                 const double *normal,
                 double weight)
    :node_idx_(node_idx),  weight_(weight), sqrt_7_(sqrt(7.0)) {
    matrixd zyz(3), Rnz_sh(9, 9);
    rot_n_2_z_by_zyz(normal, &zyz[0]);
    calc_rot_cubic_f_sh_mat_(&Rnz_sh[0], &zyz[0]);
    sh_ = trans(Rnz_sh(4, colon()));
    num_x_ = tet_num * 9;
  }
  virtual size_t dim_of_x(void) const {
    return num_x_;
  }
  virtual size_t dim_of_f(void) const {
    return 1;
  }
  virtual int val(const double *x, double *f, func_ctx *ctx = 0) const {
    const itr_matrix<const double *> x0(9, num_x_/9, x);
    const matrixd x1 = x0(colon(), node_idx_);
    *f = (dot(x1, sh_)-sqrt_7_)*weight_;
    return 0;
  }
  virtual int jac(const double *x, double *val, int32_t *ptr = 0, int32_t *idx = 0, func_ctx *ctx = 0) const {

    ptr[1] = ptr[0] + 9;
    for(int32_t d = 0; d < 9; ++d) {
        idx[ptr[0]+d] = node_idx_*9+d;
        val[ptr[0]+d] = sh_[d]*weight_;
      }
    return 0;
  }
  virtual size_t jac_nnz(void) const {
    return 9;
  }
protected:
  size_t num_x_;
  matrixd sh_;
  const size_t node_idx_;
  const double weight_;
  const double sqrt_7_;
};

class surface_smooth_func : public function_t<double,int32_t>
{
public:
  surface_smooth_func(const matrixd & zyz,
                      const pair<size_t,size_t> & tet_pair,
                      const size_t tet_num,
                      const double w)
    :tet_pair_(tet_pair), tet_num_(tet_num), w_(w){
    assert(tet_pair.first < tet_num && tet_pair.second < tet_num);
    zyz_sh_.resize(9,9);
    calc_rot_cubic_f_sh_mat_(&zyz_sh_[0], &zyz[0]);
  }

  virtual size_t dim_of_f() const
  {
    return 9;
  }
  virtual size_t dim_of_x() const
  {
    return 9 * tet_num_;
  }

  virtual int val(const value_type *x, value_type *f, hj::function::func_ctx *ctx) const
  {
    itr_matrix<const value_type*> x0(9,tet_num_, x);
    itr_matrix<value_type*> f0(9,1,f);
    f0 = w_*(zyz_sh_ * x0(colon(),tet_pair_.first) - x0(colon(), tet_pair_.second));
    return 0;
  }
  virtual int jac(const value_type *x, value_type *val, int_type *ptr, int_type *idx,
                  hj::function::func_ctx *ctx) const
  {
    for(size_t i = 0; i < 9; ++i){
        ptr[i+1] = ptr[i] + 10;
        for(size_t j = 0;  j < 9; ++j){
            idx[ptr[i]+j] = 9*tet_pair_.first + j;
            val[ptr[i]+j] = w_ * zyz_sh_(i,j);
          }
        idx[ptr[i]+9] = 9*tet_pair_.second + i;
        val[ptr[i]+9] = -1*w_;
      }
    return 0;
  }
  virtual size_t jac_nnz()
  {
    return 90;
  }
private:
  matrixd zyz_sh_;
  const pair<size_t,size_t> tet_pair_;
  const size_t tet_num_;
  const double w_;
};

class surface_smooth_jtf_func : public jtf::function::functionN1_t<double,int32_t>
{
public:
  surface_smooth_jtf_func(const matrixd & zyz,
                          const pair<size_t,size_t> & tet_pair,
                          const size_t tet_num,
                          const double w)
    :tet_pair_(tet_pair), tet_num_(tet_num), w_(w){
    assert(tet_pair.first < tet_num && tet_pair.second < tet_num);
    zyz_sh_.resize(9,9);
    calc_rot_cubic_f_sh_mat_(&zyz_sh_[0], &zyz[0]);
  }

  virtual ~surface_smooth_jtf_func(){}

  virtual size_t dim() const
  {
    return 9 * tet_num_;
  }

  virtual int val(const double *x, double &v)
  {
    itr_matrix<const double*> x0(9, tet_num_, x);
    matrix<double> Ax1_x2 = zyz_sh_ * x0(colon(), tet_pair_.first) - x0(colon(), tet_pair_.second);
    v += 0.5*w_*dot(Ax1_x2, Ax1_x2);
    return 0;
  }

  virtual int gra(const double *x, double *g)
  {
    itr_matrix<const double*> x0(9,tet_num_,x);
    static matrix<double> x2(9,2);
    static matrix<double> jac(9,2);
    x2(colon(),0) = x0(colon(), tet_pair_.first);
    x2(colon(),1) = x0(colon(), tet_pair_.second);
    ax_jac_(&jac[0], &x2[0], &zyz_sh_[0]);
    for(size_t i = 0; i < 9; ++i){
        g[9*tet_pair_.first+i] += w_*jac(i,0);
        g[9*tet_pair_.second+i] += w_*jac(i,1);
      }
    return 0;
  }
  virtual int gra(const double *x, size_t &nnz, double *g, int_type *idx)
  {
    if(g == 0 && idx == 0){
        nnz = 18;
        return 0;
      }

    itr_matrix<const double*> x0(9,tet_num_,x);
    static matrix<double> x2(9,2);
    static matrix<double> jac(9,2);
    x2(colon(),0) = x0(colon(), tet_pair_.first);
    x2(colon(),1) = x0(colon(), tet_pair_.second);
    ax_jac_(&jac[0], &x2[0], &zyz_sh_[0]);

    for(size_t i = 0; i < 9; ++i){
        idx[i] = 9*tet_pair_.first + i;
        idx[i+9] = 9*tet_pair_.second+i;
        g[i] = jac(i,0)*w_;
        g[i+9] = jac(i,1)*w_;
      }
    return 0;
  }
  virtual int hes(const double *x, size_t &nnz, size_t &format,
                  double *h, int32_t *ptr, int32_t *idx, double alpha)
  {
    if(h == 0 && ptr == 0 && idx == 0){
        nnz = 18*18;
        format =1 ;
        return 0;
      }
    if(h == 0&& ptr != 0 && idx!= 0){
        size_t two_tet_idx[] = {tet_pair_.first, tet_pair_.second};
        if(two_tet_idx[0] > two_tet_idx[1])
          std::swap(two_tet_idx[0], two_tet_idx[1]);
        for(size_t i = 0; i < 9; ++i){
            ptr[9*two_tet_idx[0] + i+1] = ptr[9*two_tet_idx[0]+i] + 18;
            for(size_t j = 0; j < 18; ++j){
                idx[ptr[9*two_tet_idx[0]+i]+j] = 9*two_tet_idx[j/9]+j%9;
              }
          }
        if(two_tet_idx[1] != two_tet_idx[0]+1)
          ptr[9*two_tet_idx[1]] = ptr[9*two_tet_idx[0]+9];
        for(size_t i = 0; i < 9; ++i){
            ptr[9*two_tet_idx[1] + i+1] = ptr[9*two_tet_idx[1]+i] + 18;
            for(size_t j = 0; j < 18; ++j){
                idx[ptr[9*two_tet_idx[1]+i]+j] = 9*two_tet_idx[j/9]+j%9;
              }
          }
        return 0;
      }
    if(h !=0 && ptr != 0 && idx != 0){
        matrix<double> hessian(18,18);
        ax_hes_(&hessian[0], x, &zyz_sh_[0]);
        hessian *= w_;
        for(size_t i = 0 ; i < 18; ++i){
            for(size_t j = 0; j < 18; ++j){
                jtf::function::add_to_csc(h, ptr, idx,
                                          9*(i/9==0?tet_pair_.first:tet_pair_.second)+i%9,
                                          9*(j/9==0?tet_pair_.first:tet_pair_.second)+j%9,
                                          hessian(i,j));
              }
          }
      }
    return 0;
  }
  virtual int hes_block(const double *x, double *h, double alpha)
  {
    return 0;
  }
private:
  matrixd zyz_sh_;
  const pair<size_t,size_t> tet_pair_;
  const size_t tet_num_;
  const double w_;
};

class tet_frame_fix : public function_t<double, int32_t>
{
public:
  tet_frame_fix(const size_t tet_num, size_t node_idx, const matrixd &sh, double weight)
    :node_idx_(node_idx), sh_(sh), weight_(weight) { num_x_ = tet_num * 9;  }
  virtual size_t dim_of_x() const{ return num_x_;}
  virtual size_t dim_of_f(void) const {
    return 9;
  }
  virtual int val(const double *x, double *f, func_ctx *ctx = 0) const {
    const itr_matrix<const double *> x0(9, num_x_/9, x);
    itr_matrix<double *> f0(9, 1, f);
    const matrixd x1 = x0(colon(), node_idx_);
    f0 = (x1-sh_)*weight_;
    return 0;
  }
  virtual int jac(const double *x, double *val, int32_t *ptr = 0, int32_t *idx = 0, func_ctx *ctx = 0) const {
    for(size_t i = 0; i < 9; ++i) {
        ptr[i+1] = ptr[i]+1;
        idx[ptr[i]] = node_idx_*9+i;
        val[ptr[i]] = weight_;
      }
    return 0;
  }
  virtual size_t jac_nnz(void) const {
    return dim_of_f();
  }
protected:
  size_t num_x_;
  const matrixd sh_;
  const size_t node_idx_;
  const double weight_;
};

class align_sh : public s2i_function
{
public:
  align_sh(const jtf::mesh::face2tet_adjacent &fa, size_t node_idx, const double *normal, double weight)
    :s2i_function(fa), node_idx_(node_idx),  weight_(weight),sqrt_7_(sqrt(7.0)) {
    matrixd zyz(3), Rnz_sh(9, 9);
    rot_n_2_z_by_zyz(normal, &zyz[0]);
    calc_rot_cubic_f_sh_mat_(&Rnz_sh[0], &zyz[0]);
    sh_ = trans(Rnz_sh(4, colon()));
  }

  virtual size_t dim_of_f(void) const {
    return 1;
  }
  virtual int val(const double *x, double *f, func_ctx *ctx = 0) const {
    const itr_matrix<const double *> x0(9, num_x_/9, x);
    const matrixd x1 = x0(colon(), node_idx_);
    *f = (dot(x1, sh_)-sqrt_7_)*weight_;
    //cerr << "# node_idx_ = " << node_idx_ << ", f value = " << *f << endl;
    return 0;
  }
  virtual int jac(const double *x, double *val, int32_t *ptr = 0, int32_t *idx = 0, func_ctx *ctx = 0) const {
    ptr[1] = ptr[0] + 9;
    for(int32_t d = 0; d < 9; ++d) {
        idx[ptr[0]+d] = node_idx_*9+d;
        val[ptr[0]+d] = sh_[d]*weight_;
      }
    return 0;
  }
  virtual size_t jac_nnz(void) const {
    return 9;
  }
protected:
  matrixd sh_;
  const size_t node_idx_;
  const double weight_;
  const double sqrt_7_;
};
//const double align_sh::sqrt_7_ = sqrt(7.0);

class normalize_sh : public s2i_function
{
public:
  normalize_sh(const jtf::mesh::face2tet_adjacent &fa, size_t node_idx, double weight)
    :s2i_function(fa), node_idx_(node_idx),  weight_(weight) {
  }

  virtual size_t dim_of_f(void) const {
    return 1;
  }
  virtual int val(const double *x, double *f, func_ctx *ctx = 0) const {
    const itr_matrix<const double *> x0(9, num_x_/9, x);
    const matrixd x1 = x0(colon(), node_idx_);
    *f = (dot(x1, x1)-12.0)*weight_;
    return 0;
  }
  virtual int jac(const double *x, double *val, int32_t *ptr = 0, int32_t *idx = 0, func_ctx *ctx = 0) const {
    const itr_matrix<const double *> x0(9, num_x_/9, x);
    ptr[1] = ptr[0] + 9;
    for(int32_t d = 0; d < 9; ++d) {
        idx[ptr[0]+d] = node_idx_*9+d;
        val[ptr[0]+d] = 2*x0(d, node_idx_)*weight_;
      }
    return 0;
  }
  virtual size_t jac_nnz(void) const {
    return 9;
  }
protected:
  matrixd sh_;
  const size_t node_idx_;
  const double weight_;
};

#endif // FUNCTION_TERM_H
