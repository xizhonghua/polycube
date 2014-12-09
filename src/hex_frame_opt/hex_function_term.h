#ifndef HEX_FUNCTION_TERM_H
#define HEX_FUNCTION_TERM_H

#include <vector>
#include <iostream>
#include <fstream>
#include <numeric>

#include <omp.h>

#include <zjucad/matrix/itr_matrix.h>
#include <zjucad/matrix/io.h>

#include "../tetmesh/tetmesh.h"

#include <hjlib/math/blas_lapack.h>
#include <zjucad/matrix/lapack.h>

#include <hjlib/function/func_aux.h>
#include <hjlib/function/operation.h>
#include <hjlib/sparse/sparse.h>

#include "../common/transition_type.h"
#include "../common/zyz.h"
#include "../common/util.h"
#include "../common/transition.h"

extern "C" {
#include "../spherical_harmonics/rot_cubic_f_SH.h"
}

using namespace zjucad::matrix;
using namespace hj::sparse;
using namespace hj::function;
using namespace std;

template <typename VAL_TYPE, typename INT_TYPE, template <typename CON> class PTR, typename CON>
class catenated_function_omp : public catenated_function<VAL_TYPE, INT_TYPE, PTR, CON>
{
public:
  typedef VAL_TYPE value_type;
  typedef INT_TYPE int_type;

  catenated_function_omp(const CON &funcs)
    :catenated_function<VAL_TYPE, INT_TYPE, PTR, CON>(funcs) {
    setup_dim_of_f();
  }
  catenated_function_omp(PTR<const CON> funcs)
    :catenated_function<VAL_TYPE, INT_TYPE, PTR, CON>(funcs) {
    setup_dim_of_f();
  }

  typedef typename catenated_function<VAL_TYPE, INT_TYPE, PTR, CON>::catenated_func_ctx catenated_func_ctx;
  virtual int val(const value_type *x, value_type *f, func_ctx *ctx = 0) const {
    if(ctx) {
        catenated_func_ctx *ctxs = dynamic_cast<catenated_func_ctx *>(ctx);
        assert(ctxs && ctxs->size() == this->funcs_.size());
        func_ctx **ctx = &(*ctxs->begin());
        size_t i;
#pragma omp parallel for private(i)
        for(i = 0; i < this->funcs_.size(); ++i) {
            this->funcs_[i]->val(x, f+dim_of_f_[i], ctx[i]);
          }
      }
    else {
        size_t i;
#pragma omp parallel for private(i)
        for(i = 0; i < this->funcs_.size(); ++i) {
            this->funcs_[i]->val(x, f+dim_of_f_[i], 0);
          }
      }
    return 0;
  }
  virtual int jac(const value_type *x, value_type *val, int_type *ptr = 0, int_type *idx = 0, func_ctx *ctx = 0) const {
    if(ctx) {
        catenated_func_ctx *ctxs = dynamic_cast<catenated_func_ctx *>(ctx);
        func_ctx **ctx = &(*ctxs->begin());
        assert(ctxs && ctxs->size() == this->funcs_.size());
        size_t i;
#pragma omp parallel for private(i)
        for(i = 0; i < this->funcs_.size(); ++i) {
            ptr[dim_of_f_[i]] = nnz_[i];
            this->funcs_[i]->jac(x, val, ptr+dim_of_f_[i], idx, ctx[i]);
          }
      }
    else {
        size_t i;
#pragma omp parallel for private(i)
        for(i = 0; i < this->funcs_.size(); ++i) {
            ptr[dim_of_f_[i]] = nnz_[i];
            this->funcs_[i]->jac(x, val, ptr+dim_of_f_[i], idx, 0);
          }
      }
    return 0;
  }
protected:
  void setup_dim_of_f(void) {
    dim_of_f_.resize(this->funcs_.size());
    nnz_.resize(this->funcs_.size());
    dim_of_f_[0] = 0;
    nnz_[0] = 0;
    for(size_t i = 1; i < this->funcs_.size(); ++i) {
        dim_of_f_[i] = dim_of_f_[i-1] + this->funcs_[i-1]->dim_of_f();
        nnz_[i] = nnz_[i-1] + this->funcs_[i-1]->jac_nnz();
      }
  }
  vector<size_t> dim_of_f_;
  vector<size_t> nnz_;
};

template <typename VAL_TYPE, typename INT_TYPE, typename CON>
function_t<VAL_TYPE, INT_TYPE> *new_catenated_function_omp(const CON &con)
{
  return new catenated_function_omp<VAL_TYPE, INT_TYPE, nil_ptr, CON>(con);
}

template <typename VAL_TYPE, typename INT_TYPE, template <typename FUNC> class PTR, typename CON>
function_t<VAL_TYPE, INT_TYPE> *new_catenated_function_omp(PTR<CON> &con)
{
  return new catenated_function_omp<VAL_TYPE, INT_TYPE, PTR, CON>(PTR<const CON>(con));
}

class sym_frame_opt_sys
{
public:
  sym_frame_opt_sys(size_t node_num) {
    sh_.resize(9, node_num);
    jac_sh_.resize(node_num);
    for(size_t i = 0; i < node_num; ++i)
      jac_sh_[i].resize(9, 3);
    is_val_jac_dirty[0] = true;
    is_val_jac_dirty[1] = true;
  }
  virtual ~sym_frame_opt_sys(){}
  void val(const double *zyz) {
    if(is_val_jac_dirty[0] == false) return;
    size_t ni;
#pragma omp parallel for private(ni)
    for(ni = 0; ni < sh_.size(2); ++ni)
      calc_rot_cubic_f_sh_(&sh_(0, ni), zyz+ni*3);
    is_val_jac_dirty[0] = false;
  }
  void jac(const double *zyz) {
    if(is_val_jac_dirty[1] == false) return;
    size_t ni;
#pragma omp parallel for private(ni)
    for(ni = 0; ni < sh_.size(2); ++ni)
      calc_jac_rot_cubic_f_sh_(&jac_sh_[ni][0], zyz+ni*3);
    is_val_jac_dirty[1] = false;
  }
  matrixd sh_;
  zjucad::matrix::matrix<matrixd > jac_sh_;
  bool is_val_jac_dirty[2];
};

class sym_frame_opt_sys_func_ctx : public func_ctx
{
public:
  sym_frame_opt_sys_func_ctx(sym_frame_opt_sys &sys)
    :sys_(sys) {
  }
  sym_frame_opt_sys &sys_;
};

class sym_frame_opt_function : public function_t<double, int32_t>
{
public:
  sym_frame_opt_function(sym_frame_opt_sys &sys)
    :sys_(sys) {
  }
  // virtual void x_changed(void) {
  // 	sys_.is_val_jac_dirty[0] = true;
  // 	sys_.is_val_jac_dirty[1] = true;
  // }

  virtual func_ctx *new_ctx(const double *x) const {
    //cout << "call new_ctx" << endl;
    unique_ptr<sym_frame_opt_sys_func_ctx> ctx(new sym_frame_opt_sys_func_ctx(sys_));
    ctx->sys_.is_val_jac_dirty[0] = true;
    ctx->sys_.is_val_jac_dirty[1] = true;
    return ctx.release();
  }
protected:
  sym_frame_opt_sys &sys_;
};



class fix_function : public sym_frame_opt_function
{
public:
  fix_function(sym_frame_opt_sys &sys, int32_t idx, const double *sh)
    :sym_frame_opt_function(sys), idx_(idx), sh_(9, 1, sh) {
  }
  virtual size_t dim_of_x(void) const {
    return sys_.sh_.size(2)*3;
  }
  virtual size_t dim_of_f(void) const {
    return 9;
  }
  virtual int val(const double *x, double *f, func_ctx *ctx = 0) const {
    sys_.val(x);
    itr_matrix<double *> f0(9, 1, f);
    f0 = sys_.sh_(colon(), idx_)-sh_;
    return 0;
  }
  virtual int jac(const double *x, double *val, int32_t *ptr = 0, int32_t *idx = 0, func_ctx *ctx = 0) const {
    sys_.jac(x);
    for(int32_t i = 0; i < dim_of_f(); ++i) {
        ptr[i+1] = ptr[i] + 3;
        for(int32_t d = 0; d < 3; ++d) {
            idx[ptr[i]+d] = idx_*3+d;
            val[ptr[i]+d] = sys_.jac_sh_[idx_](i, d);
          }
      }
    return 0;
  }
  virtual size_t jac_nnz(void) const {
    return 3*dim_of_f();
  }
  const int32_t idx_;
  const itr_matrix<const double *> sh_;
};

class fix_function_inner : public sym_frame_opt_function
{
public:
  fix_function_inner(sym_frame_opt_sys &sys, int32_t idx, const double *sh)
    :sym_frame_opt_function(sys), idx_(idx), sh_(9, 1, sh) {
  }
  virtual size_t dim_of_x(void) const {
    return sys_.sh_.size(2)*3;
  }
  virtual size_t dim_of_f(void) const {
    return 9;
  }
  virtual int val(const double *x, double *f, func_ctx *ctx = 0) const {
    sys_.val(x);
    itr_matrix<double *> f0(9, 1, f);
    f0 = sys_.sh_(colon(), idx_)-sh_;
    return 0;
  }
  virtual int jac(const double *x, double *val, int32_t *ptr = 0, int32_t *idx = 0, func_ctx *ctx = 0) const {
    sys_.jac(x);
    for(int32_t i = 0; i < dim_of_f(); ++i) {
        ptr[i+1] = ptr[i] + 3;
        for(int32_t d = 0; d < 3; ++d) {
            idx[ptr[i]+d] = idx_*3+d;
            val[ptr[i]+d] = sys_.jac_sh_[idx_](i, d);
          }
      }
    return 0;
  }
  virtual size_t jac_nnz(void) const {
    return 3*dim_of_f();
  }
  const int32_t idx_;
  const itr_matrix<const double *> sh_;
};

//class align_function : public sym_frame_opt_function
//{
//public:
//    align_function(sym_frame_opt_sys &sys, int32_t idx, const double *n)
//        :sym_frame_opt_function(sys), idx_(idx), sqrt_7_(sqrt(7.0)) {
//        matrixd zyz(3), Rnz_sh(9, 9);
//        rot_n_2_z_by_zyz(n, &zyz[0]);
//        calc_rot_cubic_f_sh_mat_(&Rnz_sh[0], &zyz[0]);
//        sh_ = trans(Rnz_sh(4, colon()));
//    }
//    virtual size_t dim_of_x(void) const {
//        return sys_.sh_.size(2)*3;
//    }
//    virtual size_t dim_of_f(void) const {
//        return 1;
//    }
//    virtual int val(const double *x, double *f, func_ctx *ctx = 0) const {
//        sys_.val(x);
//        *f = (dot(sys_.sh_(colon(), idx_), sh_)-sqrt_7_);
//        return 0;
//    }
//    virtual int jac(const double *x, double *val, int32_t *ptr = 0, int32_t *idx = 0, func_ctx *ctx = 0) const {
//        sys_.jac(x);
//        ptr[1] = ptr[0] + 3;
//        for(int32_t d = 0; d < 3; ++d) {
//            idx[ptr[0]+d] = idx_*3+d;
//            val[ptr[0]+d] = dot(sys_.jac_sh_[idx_](colon(), d), sh_);
//        }
//        return 0;
//    }
//    virtual size_t jac_nnz(void) const {
//        return 3*dim_of_f();
//    }

//    int32_t idx_;
//    matrixd sh_;
//    const double sqrt_7_;
//};

class align_function : public sym_frame_opt_function//hj::function::function_t<double, int32_t>
{
public:
  align_function(const std::shared_ptr<sym_frame_opt_sys> sfptr,
                 size_t tet_idx,
                 const double *normal,
                 const double weight)
    :sym_frame_opt_function(*sfptr), idx_(tet_idx), sqrt_7_(sqrt(7.0)),  w_(weight), sfptr_(sfptr) {
    matrixd zyz(3), Rnz_sh(9, 9);
    rot_n_2_z_by_zyz(normal, &zyz[0]);
    calc_rot_cubic_f_sh_mat_(&Rnz_sh[0], &zyz[0]);
    sh_ = trans(Rnz_sh(4, colon()));
  }
  virtual size_t dim_of_x(void) const {
    return sys_.sh_.size(2)*3;
  }
  virtual size_t dim_of_f(void) const {
    return 1;
  }
  virtual int val(const double *x, double *f, func_ctx *ctx = 0) const {
    sys_.val(x);
    *f = (dot(sys_.sh_(colon(), idx_), sh_)-sqrt_7_) * w_;
    return 0;
  }
  virtual int jac(const double *x, double *val, int32_t *ptr = 0, int32_t *idx = 0, func_ctx *ctx = 0) const {
    sys_.jac(x);
    ptr[1] = ptr[0] + 3;
    for(int32_t d = 0; d < 3; ++d) {
        idx[ptr[0]+d] = idx_*3+d;
        val[ptr[0]+d] = dot(sys_.jac_sh_[idx_](colon(), d), sh_) * w_;
      }
    return 0;
  }
  virtual size_t jac_nnz(void) const {
    return 3*dim_of_f();
  }

private:
  size_t idx_;
  matrixd sh_;
  const double w_;
  const double sqrt_7_;
  const std::shared_ptr<sym_frame_opt_sys> sfptr_;
};

//const double align_function::sqrt_7_ = sqrt(7.0);


class align_function_LP : public sym_frame_opt_function
{
public:
  align_function_LP(sym_frame_opt_sys &sys, int32_t idx, const double *n,double LP_surface)
    :sym_frame_opt_function(sys), idx_(idx), LP(LP_surface), sqrt_7_(sqrt(7.0)){
    matrixd zyz(3), Rnz_sh(9, 9);
    rot_n_2_z_by_zyz(n, &zyz[0]);
    calc_rot_cubic_f_sh_mat_(&Rnz_sh[0], &zyz[0]);
    sh_ = trans(Rnz_sh(4, colon()));
  }
  virtual size_t dim_of_x(void) const {
    return sys_.sh_.size(2)*3;
  }
  virtual size_t dim_of_f(void) const {
    return 1;
  }
  virtual int val(const double *x, double *f, func_ctx *ctx = 0) const {
    sys_.val(x);
    *f = (dot(sys_.sh_(colon(), idx_), sh_)-sqrt_7_);
    *f = pow(*f,LP);
    return 0;
  }
  virtual int jac(const double *x, double *val, int32_t *ptr = 0, int32_t *idx = 0, func_ctx *ctx = 0) const {
    sys_.jac(x);
    ptr[1] = ptr[0] + 3;
    const double f = (dot(sys_.sh_(colon(), idx_), sh_)-sqrt_7_);
    for(int32_t d = 0; d < 3; ++d) {
        idx[ptr[0]+d] = idx_*3+d;
        val[ptr[0]+d] = dot(sys_.jac_sh_[idx_](colon(), d), sh_) * LP * pow(f,LP-1);
      }
    return 0;
  }
  virtual size_t jac_nnz(void) const {
    return 3*dim_of_f();
  }

  int32_t idx_;
  matrixd sh_;
  double LP;
  const double sqrt_7_;
};

//const double align_function_LP::sqrt_7_ = sqrt(7.0);

class gradient_function_f : public sym_frame_opt_function
{
public:
  gradient_function_f(sym_frame_opt_sys &sys,
                      const matrixst &tet,
                      const matrixd &node,
                      const jtf::mesh::face2tet_adjacent &fa,
                      double weight)
    :sym_frame_opt_function(sys), weight_(weight) {
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
    M_ = P(colon(), colon(0, 2))*weight_;
  }
  virtual size_t dim_of_x(void) const {
    return sys_.sh_.size(2)*3;
  }
  virtual size_t dim_of_f(void) const {
    return 9*3;
  }
  virtual int val(const double *x, double *f, func_ctx *ctx = 0) const {
    sys_.val(x);
    itr_matrix<double *> f0(9, 3, f);
    f0 = sys_.sh_(colon(), node_idx_)*M_;
#if defined L1
    // sqrt(|f0|)
    f0 = temp(sqrt(fabs(f0)));
#endif
    return 0;
  }
  virtual int jac(const double *x, double *val, int32_t *ptr = 0, int32_t *idx = 0, func_ctx *ctx = 0) const {
    sys_.jac(x);
    for(int32_t d = 0, i = 0; d < 3; ++d) {
        for(int32_t shi = 0; shi < 9; ++shi, ++i) { // for each function component
            ptr[i+1] = ptr[i] + 12;
            for(int ni = 0; ni < 4; ++ni) {
                for(int nd = 0; nd < 3; ++nd) { // for each variable
                    idx[ptr[i]+ni*3+nd] = node_idx_[ni]*3+nd;
                    val[ptr[i]+ni*3+nd] =
                        sys_.jac_sh_[node_idx_[ni]](shi, nd)*M_(ni, d);
                  }
              }
          }
      }
#if defined L1
    // diff(sqrt(|f0|))
    sys_.val(x);
    matrixd f = sys_.sh_(colon(), node_idx_)*M_;
    for(size_t ci = 0; ci < 3*9; ++ci) { // for each function component
        for(size_t nzi = ptr[ci]; nzi < ptr[ci+1]; ++nzi) {
            const double eps = 1e-16;
            if(fabs(f[ci]) < eps)
              return 0;
            val[nzi] *= 0.5*((f[ci] > 0)?1:-1)/sqrt(fabs(f[ci]));
          }
      }
#endif
    return 0;
  }
  virtual size_t jac_nnz(void) const {
    return 3*4*dim_of_f();
  }
protected:
  matrixd M_;
  matrixst node_idx_;
  const double weight_;
};

class inner_smooth_sh : public sym_frame_opt_function//public hj::function::function_t<double, int32_t>
{
public:
  inner_smooth_sh(const std::shared_ptr<sym_frame_opt_sys> sfptr,
                  const size_t ti, const size_t tj,
                  const double weight)
    :ti_(ti), tj_(tj), w_(weight), sym_frame_opt_function(*sfptr), sfptr_(sfptr){}

  virtual size_t dim_of_x() const {
    return 3 * sys_.sh_.size(2);
  }
  virtual size_t dim_of_f() const {
    return 9 * 1;
  }
  virtual int val(const double* x, double *f, func_ctx *ctx = 0) const
  {
    sys_.val(x);
    itr_matrix<double*> f0(9,1,f);
    f0 = (sys_.sh_(colon(), ti_) - sys_.sh_(colon(), tj_))* w_;

    return 0;
  }
  virtual int jac(const double *x, double *val, int32_t *ptr = 0, int32_t *idx = 0, func_ctx *ctx = 0) const
  {
    sys_.jac(x);
    // Jf: 9 * 3
    // Jf_k: (\partial(f_i-f_j)/(\alpha)  \partial(f_i-f_j)/(\beta)   \partial(f_i-f_j)/(\gamma))_k
    const size_t t_v[] = {ti_, tj_};
    const int M_[] = {1,-1};

    for(int32_t shi = 0, i = 0; shi < 9; ++i, ++shi) {
        ptr[i+1] = ptr[i] + 6;

        for(int ni = 0 ; ni < 2; ++ni)
          for(int nd = 0; nd < 3; ++nd) {
              idx[ptr[i] + ni * 3 + nd] = t_v[ni] * 3 + nd;
              val[ptr[i] + ni * 3 + nd] = sys_.jac_sh_[t_v[ni]](shi,nd) * M_[ni] * w_;
            }
      }
    return 0;
  }
  virtual size_t jac_nnz(void) const
  {
    return 2 * 3 * dim_of_f();
  }


private:
  const std::shared_ptr<sym_frame_opt_sys> sfptr_;
  const size_t ti_, tj_;
  const double w_;
};


class frame_fix_func : public hj::function::function_t<double, int32_t>
{
public:
  frame_fix_func(const size_t node_num,
                 const matrixd & aligned_frame,
                 const size_t & aligned_tet_idx,
                 const double weight)
    : node_num_(node_num), aligned_tet_idx_(aligned_tet_idx), weight_(weight){
    aligned_frame_sh_.resize(9,1);
    calc_rot_cubic_f_sh_(&aligned_frame_sh_[0], &aligned_frame[0]);
  }
  frame_fix_func(const size_t node_num, const size_t & aligned_tet_idx,
                 const matrixd & sh, const double weight)
    : node_num_(node_num), aligned_tet_idx_(aligned_tet_idx), aligned_frame_sh_(sh), weight_(weight){}
  virtual size_t dim_of_x() const {
    return 9 * node_num_;
  }
  virtual size_t dim_of_f() const {
    return 9 * 1;
  }

  virtual int val(const double * x, double * f, func_ctx * ctx = 0) const{
    zjucad::matrix::itr_matrix<const double *> x0(9, node_num_, x);
    zjucad::matrix::itr_matrix<double *> f0(9,1, f);
    f0 = x0(colon(), aligned_tet_idx_) - aligned_frame_sh_;
    f0 *= weight_;
    return 0;
  }

  virtual int jac(const double * x, double * val,  int32_t * ptr = 0, int32_t * idx = 0,
                  func_ctx * ctx = 0) const {
    for(size_t fi = 0; fi < 9; ++fi){
        ptr[fi + 1] = ptr[fi] + 1;
        idx[ptr[fi]] = 9 * aligned_tet_idx_ + fi;
        val[ptr[fi]] = weight_;
      }
    return 0;
  }
  virtual size_t jac_nnz(void) const {
    return 9;
  }
private:
  const size_t node_num_;
  matrixd aligned_frame_sh_;
  const size_t aligned_tet_idx_;
  const double weight_;
};

class gradient_function_f_inner_4face_adjacent : public sym_frame_opt_function
{
public:
  gradient_function_f_inner_4face_adjacent(sym_frame_opt_sys &sys,
                                           const vector<size_t> &tet_adjacent_idx, //TODO: assume the tet_adjacent tet has been found. so the searching process should be handled first.
                                           const matrixst &tet,
                                           const matrixd &node,
                                           const jtf::mesh::face2tet_adjacent &fa,
                                           double weight)
    :sym_frame_opt_function(sys), weight_(weight) {
    matrixd P = ones<double>(4, 4);
    node_idx_.resize(4);
    matrixd inner_center = zeros<double>(3,1);

    if(tet_adjacent_idx.size() != 4)
      {
        std::cerr << "tet adjacent num is not enought" << endl;
        return;
      }
    size_t i = 0;

    for(vector<size_t>::const_iterator it= tet_adjacent_idx.begin();
        it != tet_adjacent_idx.end(); ++it,++i)
      {
        int tet_idx = *it;

        for(int j = 0; j < 4; ++j)
          {
            int tet_vertex_idx = tet(j,tet_idx);
            inner_center += node(colon(), tet_vertex_idx);
          }
        P(colon(0, 2), i) = inner_center/4.0;
        node_idx_[i] = tet_idx;
      }
    if(inv(P))
      cerr << "inverse fail." << endl;
    M_ = P(colon(), colon(0, 2))*weight_;
  }
  virtual size_t dim_of_x(void) const {
    return sys_.sh_.size(2)*3;
  }
  virtual size_t dim_of_f(void) const {
    return 9*3;
  }
  virtual int val(const double *x, double *f, func_ctx *ctx = 0) const {
    sys_.val(x);
    itr_matrix<double *> f0(9, 3, f);
    f0 = sys_.sh_(colon(), node_idx_)*M_;
#if defined L1
    // sqrt(|f0|)
    f0 = temp(sqrt(fabs(f0)));
#endif
    return 0;
  }
  virtual int jac(const double *x, double *val, int32_t *ptr = 0, int32_t *idx = 0, func_ctx *ctx = 0) const {
    sys_.jac(x);
    for(int32_t d = 0, i = 0; d < 3; ++d) {
        for(int32_t shi = 0; shi < 9; ++shi, ++i) { // for each function component
            ptr[i+1] = ptr[i] + 12;
            for(int ni = 0; ni < 4; ++ni) {
                for(int nd = 0; nd < 3; ++nd) { // for each variable
                    idx[ptr[i]+ni*3+nd] = node_idx_[ni]*3+nd;
                    val[ptr[i]+ni*3+nd] =
                        sys_.jac_sh_[node_idx_[ni]](shi, nd)
                        *M_(ni, d);
                  }
              }
          }
      }
#if defined L1
    // diff(sqrt(|f0|))
    sys_.val(x);
    matrixd f = sys_.sh_(colon(), node_idx_)*M_;
    for(size_t ci = 0; ci < 3*9; ++ci) { // for each function component
        for(size_t nzi = ptr[ci]; nzi < ptr[ci+1]; ++nzi) {
            const double eps = 1e-16;
            if(fabs(f[ci]) < eps)
              return 0;
            val[nzi] *= 0.5*((f[ci] > 0)?1:-1)/sqrt(fabs(f[ci]));
          }
      }
#endif
    return 0;
  }
  virtual size_t jac_nnz(void) const {
    return 3*4*dim_of_f();
  }
protected:
  matrixd M_;
  matrixst node_idx_;
  const double weight_;
};

class gradient_function_f_inner_barycenter_adjacent_face: public sym_frame_opt_function
{
  // F =(f_i f_j)*(1 -1)T
public:
  gradient_function_f_inner_barycenter_adjacent_face(sym_frame_opt_sys &sys,
                                                     const size_t adjacent_tet_i,
                                                     const size_t adjacent_tet_j
                                                     ):sym_frame_opt_function(sys),vi(adjacent_tet_i),vj(adjacent_tet_j){}
  virtual size_t dim_of_x(void) const{
    return sys_.sh_.size(2)*3;
  }
  virtual size_t dim_of_f(void) const
  {
    return 9 * 1;
  }
  virtual int val(const double* x, double *f, func_ctx *ctx = 0) const
  {
    sys_.val(x);
    itr_matrix<double*> f0(9,1,f);
    f0 = sys_.sh_(colon(), vi) - sys_.sh_(colon(), vj);
    return 0;
  }
  virtual int jac(const double *x, double *val, int32_t *ptr = 0, int32_t *idx = 0, func_ctx *ctx = 0) const
  {
    sys_.jac(x);
    // jac is const (E -E) 9*18 : wrong!
    // Jf: 9 * 3
    // Jf_k: (\partial(f_i-f_j)/(\alpha)  \partial(f_i-f_j)/(\beta)   \partial(f_i-f_j)/(\gamma))_k
    vector<size_t> t_v;
    t_v.push_back(vi);
    t_v.push_back(vj);
    vector<int> M_;
    M_.push_back(1);
    M_.push_back(-1);

    for(int32_t shi = 0, i = 0; shi < 9; ++i, ++shi)
      {
        ptr[i+1] = ptr[i] + 6;

        for(int ni = 0 ; ni < 2; ++ni)
          for(int nd = 0; nd < 3; ++nd)
            {
              idx[ptr[i] + ni * 3 + nd] = t_v[ni] * 3 + nd;
              val[ptr[i] + ni * 3 + nd] = sys_.jac_sh_[t_v[ni]](shi,nd) * M_[ni];
            }
      }
    return 0;
  }
  virtual size_t jac_nnz(void) const
  {
    return 2 * 3 * dim_of_f();
  }
private:
  const size_t vi;
  const size_t vj;
};

class smooth_function_inner_barycenter_adjacent_face_for_polycube: public sym_frame_opt_function
{
  // F = (fi-fj)/d_ij
public:
  smooth_function_inner_barycenter_adjacent_face_for_polycube(
      sym_frame_opt_sys &sys,
      const size_t adjacent_tet_i,
      const size_t adjacent_tet_j
      ):sym_frame_opt_function(sys),vi(adjacent_tet_i),vj(adjacent_tet_j){
    //vector<size_t> t_v;
    t_v.push_back(vi);
    t_v.push_back(vj);
    //vector<int> M_;
    M_.push_back(1);
    M_.push_back(-1);
  }

  template<typename value_type, typename idx_type>
  value_type polycube_jac(const value_type* zyz, idx_type row, idx_type col)const
  {
    assert(row < 9 && row >= 0);
    assert(col < 3 && col >= 0);
    assert(zyz);

    const value_type &a = zyz[0];
    const value_type &b = zyz[1];
    const value_type &c = zyz[2];


    const value_type jac_matrix[9][3] =
    {
      {-cos(a)*sin(c)-sin(a)*cos(b)*cos(c),-cos(a)*sin(b)*cos(c),-cos(a)*cos(b)*sin(c)-sin(a)*cos(c)},
      {cos(a)*cos(c)-sin(a)*cos(b)*sin(c),-cos(a)*sin(b)*sin(c),cos(a)*cos(b)*cos(c)-sin(a)*sin(c)},
      {sin(a)*sin(b),-cos(a)*cos(b),0},
      {sin(a)*sin(c)-cos(a)*cos(b)*cos(c),sin(a)*sin(b)*cos(c),sin(a)*cos(b)*sin(c)-cos(a)*cos(c)},
      {-cos(a)*cos(b)*sin(c)-sin(a)*cos(c),sin(a)*sin(b)*sin(c),-cos(a)*sin(c)-sin(a)*cos(b)*cos(c)},
      {cos(a)*sin(b),sin(a)*cos(b),0},
      {0,cos(b)*cos(c),-sin(b)*sin(c)},
      {0,cos(b)*sin(c),sin(b)*cos(c)},
      {0,-sin(b),0}
    };
    return jac_matrix[row][col];
  }

  virtual size_t dim_of_x(void) const{
    return sys_.sh_.size(2)*3; // tets number * 3
  }
  virtual size_t dim_of_f(void) const
  {
    return 9 * 1;
  }
  virtual int val(const double* x, double *f, func_ctx *ctx = 0) const
  {
    itr_matrix<const double*> x0(3,sys_.sh_.size(2),x);
    itr_matrix<double*> f0(9,1,f);
    matrixd rot_vi(3,3),rot_vj(3,3);

    // TODO: need to speed up
    zyz_angle_2_rotation_matrix1(&x0(0,vi),&rot_vi[0]);
    zyz_angle_2_rotation_matrix1(&x0(0,vj),&rot_vj[0]);
    f0 = (rot_vi - rot_vj)(colon());
    //f0 = x0(colon(),vi) - x0(colon(),vj);
    return 0;
  }
  virtual int jac(const double *x, double *val, int32_t *ptr = 0, int32_t *idx = 0, func_ctx *ctx = 0) const
  {
    itr_matrix<const double*> x0(3,sys_.sh_.size(2),x);
    for(int32_t rot_i = 0, i = 0; rot_i < 9; ++i, ++rot_i)
      {
        ptr[i+1] = ptr[i] + 6;

        for(int ni = 0 ; ni < 2; ++ni)
          for(int32_t nd = 0; nd < 3; ++nd)
            {
              idx[ptr[i] + ni * 3 + nd] = t_v[ni] * 3 + nd;
              val[ptr[i] + ni * 3 + nd] = polycube_jac(&x0(0,t_v[ni]),rot_i,nd) * M_[ni];
            }
      }

    return 0;
  }
  virtual size_t jac_nnz(void) const
  {
    return 2  * 3 * dim_of_f();
  }


private:
  const size_t vi;
  const size_t vj;
  vector<size_t> t_v;
  vector<int> M_;
  //static double jac_rot[9][3];
};

// non sym frame optimzation

template<typename value_type, typename idx_type>
value_type polycube_jac(const value_type* zyz, idx_type row, idx_type col)
{
  assert(row < 9 && row >= 0);
  assert(col < 3 && col >= 0);
  assert(zyz);

  const value_type &a = zyz[0];
  const value_type &b = zyz[1];
  const value_type &c = zyz[2];

  const value_type jac_matrix[9][3] =
  {
    {-cos(a)*sin(c)-sin(a)*cos(b)*cos(c),-cos(a)*sin(b)*cos(c),-cos(a)*cos(b)*sin(c)-sin(a)*cos(c)},
    {cos(a)*cos(c)-sin(a)*cos(b)*sin(c),-cos(a)*sin(b)*sin(c),cos(a)*cos(b)*cos(c)-sin(a)*sin(c)},
    {sin(a)*sin(b),-cos(a)*cos(b),0},
    {sin(a)*sin(c)-cos(a)*cos(b)*cos(c),sin(a)*sin(b)*cos(c),sin(a)*cos(b)*sin(c)-cos(a)*cos(c)},
    {-cos(a)*cos(b)*sin(c)-sin(a)*cos(c),sin(a)*sin(b)*sin(c),-cos(a)*sin(c)-sin(a)*cos(b)*cos(c)},
    {cos(a)*sin(b),sin(a)*cos(b),0},
    {0,cos(b)*cos(c),-sin(b)*sin(c)},
    {0,cos(b)*sin(c),sin(b)*cos(c)},
    {0,-sin(b),0}
  };

  return jac_matrix[row][col];
}

class cut_inner_smooth_function: public function_t<double,int32_t>
{
  // F = (fi-fj)/d_ij
public:
  cut_inner_smooth_function(
      const size_t adjacent_tet_i,
      const size_t adjacent_tet_j,
      const size_t tet_num,
      const double LP_N
      ):vi(adjacent_tet_i),vj(adjacent_tet_j),tet_num_(tet_num),LP(LP_N){
    //vector<size_t> t_v;
    t_v.push_back(vi);
    t_v.push_back(vj);
    //vector<int> M_;
    M_.push_back(1);
    M_.push_back(-1);
  }

  virtual size_t dim_of_x(void) const{
    return tet_num_*3; // tets number * 3
  }
  virtual size_t dim_of_f(void) const
  {
    return 9 * 1;
  }
  virtual int val(const double* x, double *f, func_ctx *ctx = 0) const
  {
    itr_matrix<const double*> x0(3,tet_num_,x);
    itr_matrix<double*> f0(9,1,f);
    matrixd rot_vi(3,3),rot_vj(3,3);

    // TODO: need to speed up
    zyz_angle_2_rotation_matrix1(&x0(0,vi),&rot_vi[0]);
    zyz_angle_2_rotation_matrix1(&x0(0,vj),&rot_vj[0]);

    for(size_t t = 0; t < 9; ++t)
      f0[t] = pow((rot_vi[t] - rot_vj[t]),LP);

    //cerr << f0 << endl;
    return 0;
  }
  virtual int jac(const double *x, double *val, int32_t *ptr = 0, int32_t *idx = 0, func_ctx *ctx = 0) const
  {
    itr_matrix<const double*> x0(3,tet_num_,x);

    matrixd rot_vi(3,3),rot_vj(3,3);

    // TODO: need to speed up
    zyz_angle_2_rotation_matrix1(&x0(0,vi),&rot_vi[0]);
    zyz_angle_2_rotation_matrix1(&x0(0,vj),&rot_vj[0]);

    for(int32_t rot_i = 0, i = 0; rot_i < 9; ++i, ++rot_i)
      {
        ptr[i+1] = ptr[i] + 6;

        for(int ni = 0 ; ni < 2; ++ni)
          for(int32_t nd = 0; nd < 3; ++nd)
            {
              idx[ptr[i] + ni * 3 + nd] = t_v[ni] * 3 + nd;

              val[ptr[i] + ni * 3 + nd] = polycube_jac(&x0(0,t_v[ni]),rot_i,nd) * M_[ni] * LP * pow((rot_vi[rot_i] - rot_vj[rot_i]),LP-1);
            }
      }

    return 0;
  }
  virtual size_t jac_nnz(void) const
  {
    return 2  * 3 * dim_of_f();
  }

private:
  const size_t vi;
  const size_t vj;
  const double LP;
  const size_t tet_num_;
  vector<size_t> t_v;
  vector<int> M_;
  //static double jac_rot[9][3];
};

template <typename value_type>
double assign_type_sign(value_type t)
{
  if(t%2 == 0) return 1;
  else
    return -1;
}


class cut_surface_align_function: public function_t<double,int32_t>
{
public:
  cut_surface_align_function(
      const size_t tet_idx,
      const size_t surface_align_type,
      const size_t tet_num,
      const matrixd &normal,
      const double LP_N
      ):tet_idx_(tet_idx),type_(surface_align_type),tet_num_(tet_num),normal_(normal),LP(LP_N){}

  virtual size_t dim_of_x(void) const{
    return tet_num_*3; // tets number * 3
  }
  virtual size_t dim_of_f(void) const
  {
    return 3 * 1;
  }
  virtual int val(const double* x, double *f, func_ctx *ctx = 0) const
  {
    itr_matrix<const double*> x0(3,tet_num_,x);
    itr_matrix<double*> f0(3,1,f);
    matrixd rot_v(3,3);


    zyz_angle_2_rotation_matrix1(&x0(0,tet_idx_),&rot_v[0]);


    matrixd rotation(3,3);
    rotation = type_transition2(type_);
    //type_matrix_map(type_,rotation);
    matrixd frame_rot = rot_v * rotation;
    for(size_t t = 0; t < 3; ++t)
      f0[t] = pow((frame_rot[t] - normal_[t]),LP);
    return 0;
  }
  virtual int jac(const double *x, double *val, int32_t *ptr = 0, int32_t *idx = 0, func_ctx *ctx = 0) const
  {
    itr_matrix<const double*> x0(3,tet_num_,x);
    matrixd rot_v(3,3);
    zyz_angle_2_rotation_matrix1(&x0(0,tet_idx_),&rot_v[0]);

    matrixd f0(3,1);
    matrixd rotation(3,3);
    rotation = type_transition2(type_);
    //type_matrix_map(type_,rotation);
    matrixd frame_rot = rot_v * rotation;
    f0 = frame_rot(colon(),0) - normal_;

    for(int32_t rot_i = 0; rot_i < 3; ++rot_i)
      {
        ptr[rot_i+1] = ptr[rot_i] +  3;

        for(int32_t nd = 0; nd < 3; ++nd)
          {
            idx[ptr[rot_i]  + nd] = tet_idx_ * 3 + nd;
            val[ptr[rot_i]  + nd] = 0;
            for(int32_t nj = 0; nj < 3; ++nj)
              val[ptr[rot_i]  + nd] += rotation[rot_i] * polycube_jac(&x0(0,tet_idx_),rot_i + nj * 3 ,nd);
            val[ptr[rot_i]  + nd] *= LP * pow(f0[rot_i%3],LP - 1);
          }
      }

    return 0;
  }
  virtual size_t jac_nnz(void) const
  {
    return  3 * dim_of_f();
  }

private:
  const size_t tet_idx_;
  const size_t type_;
  const double LP;
  const size_t tet_num_;
  const matrixd normal_;
};

//#define LP 1
class cut_jump_smooth_function: public function_t<double,int32_t>
{
  // F = (fi-fj)/d_ij
public:
  cut_jump_smooth_function(
      const size_t jump_type,
      const size_t adjacent_tet_i,
      const size_t adjacent_tet_j,
      const size_t tet_num,
      const double LP_N
      ):type_(jump_type),vi(adjacent_tet_i),vj(adjacent_tet_j),tet_num_(tet_num),LP(LP_N){
    //vector<size_t> t_v;
    //assert(type_ != 9 && type_ != -1); // type can not not black(9) and unknow(-1)
    assert(adjacent_tet_i != -1 && adjacent_tet_j != -1);
    assert(adjacent_tet_j < tet_num && adjacent_tet_i < tet_num);

    t_v.push_back(vi);
    t_v.push_back(vj);
    //vector<int> M_;
    M_.push_back(1);
    M_.push_back(-1);

    type_m = type_transition2(type_);
    //type_matrix_map(type_,type_m);
    //cerr << type_m;
  }

  virtual size_t dim_of_x(void) const{
    return tet_num_*3; // tets number * 3
  }
  virtual size_t dim_of_f(void) const
  {
    return 9 * 1;
  }
  virtual int val(const double* x, double *f, func_ctx *ctx = 0) const
  {
    itr_matrix<const double*> x0(3,tet_num_,x);
    itr_matrix<double*> f0(9,1,f);
    matrixd rot_vi(3,3),rot_vj(3,3);

    zyz_angle_2_rotation_matrix1(&x0(0,vi),&rot_vi[0]);
    zyz_angle_2_rotation_matrix1(&x0(0,vj),&rot_vj[0]);
#if 0   // check
    {
      matrixd jump_rot(3,3);
      get_best_alignment(&rot_vi[0],&rot_vj[0],&jump_rot[0]);
      cerr << "# jump error " << norm(jump_rot - type_m) << endl;
      //          cerr << jump_rot << endl;
      //          cerr << type_m << endl;

    }
#endif
    f0 = (rot_vi * type_m- rot_vj)(colon());

    for(size_t t = 0; t < 9; ++t)
      f0[t] = pow(f0[t],LP);


    return 0;
  }

  virtual int jac(const double *x, double *val, int32_t *ptr = 0, int32_t *idx = 0, func_ctx *ctx = 0) const
  {
    itr_matrix<const double*> x0(3,tet_num_,x);
    //cerr << "x 22: " << x0(colon(),22) << endl;
    matrixd rot_vi(3,3),rot_vj(3,3);

    zyz_angle_2_rotation_matrix1(&x0(0,vi),&rot_vi[0]);
    zyz_angle_2_rotation_matrix1(&x0(0,vj),&rot_vj[0]);

    matrixd f0 = rot_vi * type_m - rot_vj;

    for(int32_t rot_i = 0, i = 0; rot_i < 9; ++i, ++rot_i)
      {
        ptr[i+1] = ptr[i] + 6;

        for(int ni = 0 ; ni < 2; ++ni){
            for(int32_t nd = 0; nd < 3; ++nd) // three angle of zyz
              {
                idx[ptr[i] + ni * 3 + nd] = t_v[ni] * 3 + nd;
                if(ni == 0){
                    //val[ptr[i] + ni * 3 + nd] = polycube_jac(&x0(0,t_v[ni]),rot_i,nd) * M_[ni] * type_m[rot_i];
                    val[ptr[i] + ni * 3 + nd] = 0;
                    for(int32_t j = 0; j < 3; ++j) {
                        val[ptr[i] + ni * 3 + nd] += polycube_jac(&x0(0,t_v[ni]),j * 3 + rot_i % 3,nd) * M_[ni] * type_m(j,rot_i / 3)
                            * LP * pow(f0[rot_i],LP-1);
                      }

                  }
                else{
                    val[ptr[i] + ni * 3 + nd] = polycube_jac(&x0(0,t_v[ni]),rot_i,nd) * M_[ni]
                        * LP * pow(f0[rot_i],LP-1);
                  }

              }
          }
      }

    return 0;
  }

  virtual size_t jac_nnz(void) const
  {
    return 2  * 3 * dim_of_f();
  }

private:
  const size_t vi;
  const size_t vj;
  const size_t tet_num_;
  const size_t type_;
  vector<size_t> t_v;
  vector<int> M_;
  matrixd type_m;
  const double LP;
  //static double jac_rot[9][3];
};


class cut_inner_smooth_linear_function: public function_t<double,int32_t>
{
  // F = (fi-fj)/d_ij
public:
  cut_inner_smooth_linear_function(
      const size_t adjacent_tet_i,
      const size_t adjacent_tet_j,
      const size_t tet_num,
      const double LP_smooth
      ):vi(adjacent_tet_i),vj(adjacent_tet_j),tet_num_(tet_num),LP(LP_smooth){
    //vector<size_t> t_v;
    t_v.push_back(vi);
    t_v.push_back(vj);
    //vector<int> M_;
    M_.push_back(1);
    M_.push_back(-1);
  }

  virtual size_t dim_of_x(void) const{
    return tet_num_ * 9; // tets number * 9
  }
  virtual size_t dim_of_f(void) const
  {
    return 9 * 1;
  }
  virtual int val(const double* x, double *f, func_ctx *ctx = 0) const
  {
    itr_matrix<const double*> x0(9,tet_num_,x);
    itr_matrix<double*> f0(9,1,f);


    for(size_t t = 0; t < 9; ++t)
      f0[t] = pow((x0(t,vi) - x0(t,vj)),LP);

    //        cerr << " vi - vj " << vi << " " << vj << endl;
    //        cerr << x0(colon(),vi) << x0(colon(),vj) << endl;
    //        cerr << f0 << endl;
    return 0;
  }
  virtual int jac(const double *x, double *val, int32_t *ptr = 0, int32_t *idx = 0, func_ctx *ctx = 0) const
  {
    itr_matrix<const double*> x0(9,tet_num_,x);

    for(int32_t rot_i = 0, i = 0; rot_i < 9; ++i, ++rot_i)
      {
        ptr[i+1] = ptr[i] + 2 ;

        for(int ni = 0 ; ni < 2; ++ni)
          {
            idx[ptr[i] + ni ] = t_v[ni] * 9 + rot_i;
            val[ptr[i] + ni ] =  M_[ni] * LP * pow((x0(rot_i,vi) - x0(rot_i,vj)),LP-1);
          }
      }

    return 0;
  }
  virtual size_t jac_nnz(void) const
  {
    return 2  * dim_of_f();
  }

private:
  const size_t vi;
  const size_t vj;
  const size_t tet_num_;
  const double LP;
  vector<size_t> t_v;
  vector<int> M_;
  //static double jac_rot[9][3];
};

class cut_surface_align_linear_function_new: public function_t<double,int32_t>
{
public:
  cut_surface_align_linear_function_new(
      const size_t tet_idx,
      const size_t surface_align_type,
      const size_t tet_num,
      const matrixd &normal,
      const double LP_surface):tet_idx_(tet_idx),type_(surface_align_type),tet_num_(tet_num),normal_(normal),LP(LP_surface)
  {}
  virtual size_t dim_of_x(void) const{
    return tet_num_ * 9;
  }
  virtual size_t dim_of_f(void) const{
    return 3;
  }
  virtual int val(const double *x, double *f, func_ctx * ctx = 0) const{
    itr_matrix<const double*> x0(9,tet_num_,x);
    itr_matrix<double *> f0(3,1,f);

    matrixd x_(3,3);
    x_(colon()) = x0(colon(),tet_idx_);
    x_ = temp(x_ * type_transition2(type_));
    f0 = x_(colon(),0) - normal_;

    for(size_t t = 0; t < 3; ++t)
      {
        f0[t] = pow(f0[t],LP);
      }
    return 0;
  }
  virtual int jac(const double *x, double *val, int32_t *ptr = 0, int32_t *idx = 0, func_ctx *ctx = 0) const
  {
    itr_matrix<const double*> x0(9,tet_num_,x);
    matrixd x_(3,3);
    x_(colon()) = x0(colon(),tet_idx_);

    matrixd Pi = type_transition2(type_);
    x_ = temp(x_ * Pi);
    matrixd f0 = x_(colon(),0) - normal_;

    for(int32_t t = 0; t < 3; ++t){
        ptr[t+1] = ptr[t] + 3;
        for(int32_t i = 0; i < 3; ++i){
            idx[ptr[t] + i] = tet_idx_ * 9 + t + i * 3;
            val[ptr[t] + i] = Pi(i,0) * LP * pow(f0[t],LP-1);
          }
      }
    return 0;
  }

  virtual size_t jac_nnz(void) const
  {
    return  3 * dim_of_f();
  }
private:
  const size_t tet_idx_;
  const size_t type_;
  const size_t tet_num_;
  const double LP;
  const matrixd normal_;

};


class cut_surface_align_linear_function: public function_t<double,int32_t>
{
public:
  cut_surface_align_linear_function(
      const size_t tet_idx,
      const size_t surface_align_type,
      const size_t tet_num,
      const matrixd &normal,
      const double LP_surface
      ):tet_idx_(tet_idx),type_(surface_align_type),tet_num_(tet_num),normal_(normal),LP(LP_surface){}

  virtual size_t dim_of_x(void) const{
    return tet_num_ * 9; // tets number * 9
  }
  virtual size_t dim_of_f(void) const
  {
    return 3 * 1;
  }
  virtual int val(const double* x, double *f, func_ctx *ctx = 0) const
  {
    itr_matrix<const double*> x0(9,tet_num_,x);
    itr_matrix<double*> f0(3,1,f);
    //matrixd rot(3,3);
    //rot(colon()) = x0(colon(),tet_idx_);

    for(size_t t = 0; t < 3; ++t)
      f0[t] = pow((assign_type_sign(type_) * x0(type_/2 * 3 + t,tet_idx_) - normal_[t]),LP);
    return 0;
  }
  virtual int jac(const double *x, double *val, int32_t *ptr = 0, int32_t *idx = 0, func_ctx *ctx = 0) const
  {
    itr_matrix<const double*> x0(9,tet_num_,x);

    for(int32_t rot_i = type_/2 * 3, i = 0; rot_i < (type_/2 + 1) * 3; ++i, ++rot_i)
      {
        ptr[i+1] = ptr[i] + 1;

        idx[ptr[i]] = tet_idx_ * 9 + rot_i;
        val[ptr[i]] = assign_type_sign(type_) * LP * pow((assign_type_sign(type_) * x0(rot_i,tet_idx_)- normal_[rot_i % 3]),LP-1);
      }

    return 0;
  }
  virtual size_t jac_nnz(void) const
  {
    return  2 * dim_of_f();
  }

private:
  const size_t tet_idx_;
  const size_t type_;
  const size_t tet_num_;
  const double LP;
  const matrixd normal_;
};


class cut_jump_smooth_linear_function_new: public function_t<double,int32_t>
{
  // F = (fi * \Pi -fj)
public:
  cut_jump_smooth_linear_function_new(
      const size_t jump_type,
      const size_t adjacent_tet_i,
      const size_t adjacent_tet_j,
      const size_t tet_num,
      const size_t LP_jump
      ):type_(jump_type),vi(adjacent_tet_i),vj(adjacent_tet_j),tet_num_(tet_num),LP(LP_jump){
    //vector<size_t> t_v;
    //assert(type_ != 9 && type_ != -1); // type can not not black(9) and unknow(-1)
    assert(adjacent_tet_i != -1 && adjacent_tet_j != -1);
    assert(adjacent_tet_j < tet_num && adjacent_tet_i < tet_num);

    t_v.push_back(vi);
    t_v.push_back(vj);
    //vector<int> M_;
    M_.push_back(1);
    M_.push_back(-1);
    type_m = type_transition2(type_);
    //type_matrix_map(type_,type_m);
    //cerr << type_m;
  }

  virtual size_t dim_of_x(void) const{
    return tet_num_ * 9; // tets number * 9
  }
  virtual size_t dim_of_f(void) const
  {
    return 9 * 1;
  }
  virtual int val(const double* x, double *f, func_ctx *ctx = 0) const
  {
    itr_matrix<const double*> x0(9,tet_num_,x);
    itr_matrix<double*> f0(9,1,f);

#if 0   // check
    {
      matrixd jump_rot(3,3);
      get_best_alignment(&rot_vi[0],&rot_vj[0],&jump_rot[0]);
      cerr << "# jump error " << norm(jump_rot - type_m) << endl;
      //          cerr << jump_rot << endl;
      //          cerr << type_m << endl;

    }
#endif

    matrixd xi(3,3),xj(3,3);
    xi(colon()) = x0(colon(),vi);
    xj(colon()) = x0(colon(),vj);

    f0 = (xi * type_m - xj)(colon());

    for(size_t t = 0; t < 9; ++t)
      f0[t] = pow(f0[t],LP);
    return 0;
  }

  virtual int jac(const double *x, double *val, int32_t *ptr = 0, int32_t *idx = 0, func_ctx *ctx = 0) const
  {
    itr_matrix<const double*> x0(9,tet_num_,x);

    matrixd xi(3,3),xj(3,3);
    xi(colon()) = x0(colon(),vi);
    xj(colon()) = x0(colon(),vj);

    matrixd f0 = xi * type_m - xj;

    for(int32_t rot_i = 0; rot_i < 9;++rot_i)
      {
        ptr[rot_i + 1] = ptr[rot_i] + 4;
        for(size_t t = 0; t < 3; ++t)
          {
            idx[ptr[rot_i] + t] = vi * 9 + (rot_i%3) + t * 3;
            val[ptr[rot_i] + t] = type_m[t + rot_i/3 * 3];
          }
        idx[ptr[rot_i] + 3] = vj* 9 + rot_i;
        val[ptr[rot_i] + 3] = -1;
      }

    return 0;
  }

  virtual size_t jac_nnz(void) const
  {
    return 4 * dim_of_f();
  }

private:
  const size_t vi;
  const size_t vj;
  const size_t tet_num_;
  const size_t type_;
  const double LP;
  vector<size_t> t_v;
  vector<int> M_;
  matrixd type_m;
  //static double jac_rot[9][3];
};

//#define LP 1
// no use
class cut_jump_smooth_linear_function: public function_t<double,int32_t>
{
  // F = (fi-fj)/d_ij
public:
  cut_jump_smooth_linear_function(
      const size_t jump_type,
      const size_t adjacent_tet_i,
      const size_t adjacent_tet_j,
      const size_t tet_num,
      const size_t LP_jump
      ):type_(jump_type),vi(adjacent_tet_i),vj(adjacent_tet_j),tet_num_(tet_num),LP(LP_jump){
    //vector<size_t> t_v;
    //assert(type_ != 9 && type_ != -1); // type can not not black(9) and unknow(-1)
    assert(adjacent_tet_i != -1 && adjacent_tet_j != -1);
    assert(adjacent_tet_j < tet_num && adjacent_tet_i < tet_num);

    t_v.push_back(vi);
    t_v.push_back(vj);
    //vector<int> M_;
    M_.push_back(1);
    M_.push_back(-1);
    type_m = type_transition2(type_);
    //type_matrix_map(type_,type_m);
    //cerr << type_m;
  }

  virtual size_t dim_of_x(void) const{
    return tet_num_ * 9; // tets number * 9
  }
  virtual size_t dim_of_f(void) const
  {
    return 9 * 1;
  }
  virtual int val(const double* x, double *f, func_ctx *ctx = 0) const
  {
    itr_matrix<const double*> x0(9,tet_num_,x);
    itr_matrix<double*> f0(9,1,f);

#if 0   // check
    {
      matrixd jump_rot(3,3);
      get_best_alignment(&rot_vi[0],&rot_vj[0],&jump_rot[0]);
      cerr << "# jump error " << norm(jump_rot - type_m) << endl;
      //          cerr << jump_rot << endl;
      //          cerr << type_m << endl;

    }
#endif

    matrixd xi(3,3),xj(3,3);
    xi(colon()) = x0(colon(),vi);
    xj(colon()) = x0(colon(),vj);

    f0 = (xi * type_m - xj)(colon());

    for(size_t t = 0; t < 9; ++t)
      f0[t] = pow(f0[t],LP);
    return 0;
  }

  virtual int jac(const double *x, double *val, int32_t *ptr = 0, int32_t *idx = 0, func_ctx *ctx = 0) const
  {
    itr_matrix<const double*> x0(9,tet_num_,x);

    matrixd xi(3,3),xj(3,3);
    xi(colon()) = x0(colon(),vi);
    xj(colon()) = x0(colon(),vj);

    matrixd f0 = xi * type_m - xj;
    for(int32_t rot_i = 0, i = 0; rot_i < 9; ++i, ++rot_i)
      {
        ptr[i+1] = ptr[i] + 2 ;

        for(int ni = 0 ; ni < 2; ++ni){
            // for(int32_t nd = 0; nd < 9; ++nd) // three angle of zyz
            //{
            idx[ptr[i] + ni ] = t_v[ni] * 9 + rot_i;
            if(ni == 0){
                //val[ptr[i] + ni * 3 + nd] = polycube_jac(&x0(0,t_v[ni]),rot_i,nd) * M_[ni] * type_m[rot_i];
                val[ptr[i] + ni ] = 0;
                for(int32_t j = 0; j < 3; ++j) {
                    val[ptr[i] + ni] +=  M_[ni] * type_m(j,rot_i / 3) * LP * pow(f0[rot_i],LP-1);
                  }

              }
            else{
                val[ptr[i] + ni ] = M_[ni] * LP * pow(f0[rot_i],LP-1);
              }
          }
      }

    return 0;
  }

  virtual size_t jac_nnz(void) const
  {
    return 2 * dim_of_f();
  }

private:
  const size_t vi;
  const size_t vj;
  const size_t tet_num_;
  const size_t type_;
  const double LP;
  vector<size_t> t_v;
  vector<int> M_;
  matrixd type_m;
  //static double jac_rot[9][3];
};


class cut_RTR_linear_function: public function_t<double,int32_t>
{
  // F = RTR-I
public:
  cut_RTR_linear_function(const size_t tet_num,
                          const size_t tet_idx,
                          const double LP_RTR):tet_num_(tet_num),tet_idx_(tet_idx),LP(LP_RTR){}

  virtual size_t dim_of_x(void) const{
    return tet_num_ * 9; // tets number * 9
  }
  virtual size_t dim_of_f(void) const
  {
    return 9 * 1;
  }
  virtual int val(const double* x, double *f, func_ctx *ctx = 0) const
  {
    itr_matrix<const double*> x0(9,tet_num_,x);
    itr_matrix<double*> f0(9,1,f);

    matrixd x_(3,3);
    x_(colon()) = x0(colon(),tet_idx_);

    f0(colon()) = (trans(x_) * x_ - eye<double>(3))(colon());

    for(size_t t = 0; t < 9; ++t)
      f0[t] = pow(f0[t],LP);
    return 0;
  }

  virtual int jac(const double *x, double *val, int32_t *ptr = 0, int32_t *idx = 0, func_ctx *ctx = 0) const
  {
    itr_matrix<const double*> x0(9,tet_num_,x);

    matrixd x_(3,3);
    x_(colon()) = x0(colon(),tet_idx_);

    matrixd f0 = trans(x_) * x_ - eye<double>(3);

    for(int32_t rot_i = 0, i = 0; rot_i < 9; ++i, ++rot_i)
      {
        // row 0,4 8 is 2(x0,x1,x2); 2(x3,x4,x5); 2(x6,x7,x8)
        if(rot_i % 4 == 0)
          ptr[i+1] = ptr[i] + 3;
        else
          ptr[i+1] = ptr[i] + 6 ;

        for(size_t ni = 0; ni < ptr[i+1] - ptr[i]; ++ni){
            if(rot_i % 4 == 0)
              {
                idx[ptr[i] + ni] = tet_idx_ * 9 + rot_i / 4 * 3 + ni;
                val[ptr[i] + ni] = 2 * x_(ni % 3,rot_i / 4) * LP * pow(f0[rot_i],LP-1);
              }else if(rot_i % 4 == 2)
              {
                idx[ptr[i] + ni] = tet_idx_ * 9 + ni / 3 * 3 * 2 + ni % 3;
                val[ptr[i] + ni] = x_(ni % 3, (ni / 3 == 0)?2:0) * LP * pow(f0[rot_i],LP-1);
              }else if(rot_i % 2 == 1)
              {
                if(rot_i == 1 || rot_i == 3){
                    idx[ptr[i] + ni] = tet_idx_ * 9 + ni;
                    val[ptr[i] + ni] = x_(ni % 3, 1 - ni / 3) * LP * pow(f0[rot_i],LP-1);
                  }else if(rot_i == 5 || rot_i == 7){
                    idx[ptr[i] + ni] = tet_idx_ * 9 + 3 + ni;
                    val[ptr[i] + ni] = x_(ni % 3, 2 - ni / 3) * LP * pow(f0[rot_i],LP-1);
                  }
              }
          }
      }

    return 0;
  }

  virtual size_t jac_nnz(void) const
  {
    return 45; // 3 * 3 + 6 * 6
  }

private:
  const size_t tet_num_;
  const size_t tet_idx_;
  const double LP;
};

// not correct
class cut_RTR_linear_function_new: public function_t<double,int32_t>
{
  // F = RTR-I
public:
  cut_RTR_linear_function_new(const size_t tet_num,
                              const size_t tet_idx,
                              const double LP_RTR):tet_num_(tet_num),tet_idx_(tet_idx),LP(LP_RTR){}

  virtual size_t dim_of_x(void) const{
    return tet_num_ * 9; // tets number * 9
  }
  virtual size_t dim_of_f(void) const
  {
    return 9 * 1;
  }
  virtual int val(const double* x, double *f, func_ctx *ctx = 0) const
  {
    itr_matrix<const double*> x0(9,tet_num_,x);
    itr_matrix<double*> f0(9,1,f);

    matrixd x_(3,3);
    x_(colon()) = x0(colon(),tet_idx_);

    f0(colon()) = (trans(x_) * x_ - eye<double>(3))(colon());

    for(size_t t = 0; t < 9; ++t)
      f0[t] = pow(f0[t],LP);
    return 0;
  }

  virtual int jac(const double *x, double *val, int32_t *ptr = 0, int32_t *idx = 0, func_ctx *ctx = 0) const
  {
    itr_matrix<const double*> x0(9,tet_num_,x);

    matrixd x_(3,3);
    x_(colon()) = x0(colon(),tet_idx_);

    matrixd f0 = trans(x_) * x_ - eye<double>(3);

    for(int32_t i = 0; i < 9; ++i)
      {
        ptr[i+1] = ptr[i] + 1;
        idx[ptr[i] + 1] = tet_idx_ * 9 + i;
        val[ptr[i] + 1] = 2 * x_[i] * LP * pow(f0[i],LP-1);
      }

    return 0;
  }

  virtual size_t jac_nnz(void) const
  {
    return 9;
  }

private:
  const size_t tet_num_;
  const size_t tet_idx_;
  const double LP;
};

#endif // HEX_FUNCTION_TERM_H
