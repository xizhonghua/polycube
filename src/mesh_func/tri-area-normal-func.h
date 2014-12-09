#ifndef TRIAREANORMALFUNC_H
#define TRIAREANORMALFUNC_H

#include <jtflib/function/function.h>
#include <jtflib/function/func_aux.h>

#include <hjlib/math/blas_lapack.h>
#include <zjucad/matrix/lapack.h>
#include <zjucad/matrix/matrix.h>
#include "tri-area-normal.h"
#include "tri-area.h"


static inline bool is_degenerate_tri(const zjucad::matrix::matrix<double> &tri)
{
  using namespace zjucad::matrix;
    matrix<double> e1 = tri(colon(), 0)-tri(colon(), 2),
            e2 = tri(colon(), 1)-tri(colon(), 2);
    if(norm(cross(e1, e2)) < 1e-8) {
        return true;
    }
    return false;
}

class area_sum : public jtf::function::functionN1_t<double,int32_t>
{
public:
  area_sum(size_t node_num, const zjucad::matrix::matrix<size_t> &faces, double ori_area)
    :node_num_(node_num), faces_(faces), ori_area_(ori_area) {
  }
  virtual size_t dim(void) const { return node_num_*3; }
  //NOTE: add to v
  virtual int val(const double *x, double &v) {
    using namespace zjucad::matrix;
    const zjucad::matrix::itr_matrix<const double *> T(3, node_num_, x);
    matrix<double> tri(3, 3);
    v = -ori_area_;
    for(size_t fi = 0; fi < faces_.size(2); ++fi) {
        tri = T(colon(), faces_(colon(), fi));
        double a = 0;
        calc_tri_area_(&a, &tri[0]);
        v += a;
      }
    return 0;
  }
  //NOTE: add to v
  virtual int gra(const double *x, double *g) {
    using namespace zjucad::matrix;
    const zjucad::matrix::itr_matrix<const double *> T(3, node_num_, x);
    matrix<double> tri(3, 3), g9(9);
    for(size_t fi = 0; fi < faces_.size(2); ++fi) {
        tri = T(colon(), faces_(colon(), fi));
        g9(colon()) = 0;
        calc_tri_area_jac_(&g9[0], &tri[0]);
        for(size_t i = 0; i < 9; ++i) {
            g[i%3+faces_(i/3, fi)*3] += g9[i];
          }
      }
    return 0;
  }
  virtual int gra(const double *x, size_t & nnz, double * val, int32_t * idx){
    using namespace zjucad::matrix;
    if(val == 0 && idx == 0){
        nnz = node_num_*3;
      }else{
        itr_matrix<double*> val_(node_num_*3,1,val);
        val_ *= 0;
        gra(x, val);
        for(size_t i = 0; i < node_num_*3; ++i) idx[i] = i;
      }
    return 0;
  }

  virtual int hes(const double *x, size_t &nnz, size_t &format, double *h, int32_t *ptr, int32_t *idx, double alpha = 1) {
    std::cerr << "should not be here." << std::endl;
    exit(1);
    return 1;
  }
  virtual int hes_block(const double *x, double *h, double alpha = 1) {return -1;}
protected:
  const zjucad::matrix::matrix<size_t> &faces_;
  const size_t node_num_;
  double ori_area_;
};

class face_area_normal : public jtf::function::functionN1_t<double,int32_t>
{
public:
  face_area_normal(const zjucad::matrix::matrix<double> &node,
                   const zjucad::matrix::matrix<size_t> &tri,
                   size_t xyz,
                   const double weight)
    :tri_(tri), node_num_(node.size(2)), weight_(weight), xyz_(xyz) {
  }
  virtual size_t dim(void) const {
    return node_num_*3;
  }
  virtual int val(const double *x, double &v) {
    const zjucad::matrix::itr_matrix<const double *> T(3, node_num_, x);
    zjucad::matrix::matrix<double> tri = T(zjucad::matrix::colon(), tri_);
    double f[3];
    calc_tri_area_normal_(f, &tri[0]);
    jtf::math::erase_nan_inf(f[xyz_]);
    v += f[xyz_] * weight_;
    return 0;
  }

  virtual int gra(const double *x, double *g) {
    double sp_g[9];
    int32_t idx[9];
    size_t nnz = 9;
    gra(x, nnz, sp_g, idx);
    for(size_t i = 0; i < 9; ++i)
      g[idx[i]] += sp_g[i];
    return 0;
  }

  virtual int gra(const double *x, size_t &nnz, double *g, int32_t *idx) {
    if(g == 0) {
        nnz = 9;
        return 0;
      }
    if(nnz != 9) {
        std::cerr << "strange." << std::endl;
      }
    const zjucad::matrix::itr_matrix<const double *> T(3, node_num_, x);
    zjucad::matrix::matrix<double> tri =
        T(zjucad::matrix::colon(), tri_), jac(3, 9);
    calc_tri_area_normal_jac_(&jac[0], &tri[0]);
    for(size_t i = 0; i < 9; ++i) {
        idx[i] = tri_[i/3]*3+i%3;
        g[i] = jac(xyz_, i)*weight_;
        jtf::math::erase_nan_inf(g[i]);
      }
    return 0;
  }
  virtual int hes(const double *x, size_t &nnz, size_t &format, double *h, int32_t *ptr, int32_t *idx, double alpha = 1) {
    using namespace zjucad::matrix;
    format = 1;
    if(h == 0 && ptr == 0 && idx == 0) { // query nnz
        nnz = 9*9;
        return 0;
      }
    if(h == 0 && ptr != 0 && idx != 0) { // query pattern
        for(size_t ci = 0; ci < 9; ++ci) {
            const size_t var_ci = tri_[ci/3]*3+ci%3;
            ptr[var_ci+1] = ptr[var_ci]+9;
            for(size_t ri = 0; ri < 9; ++ri) {
                const size_t var_ri = tri_[ri/3]*3+ri%3;
                idx[ptr[var_ci]+ri] = var_ri;
              }
          }
        return 0;
      }
    if(h != 0 && ptr != 0 && idx != 0) { // accumulate
        const zjucad::matrix::itr_matrix<const double *> T(3, node_num_, x);
        zjucad::matrix::matrix<double> tri =
            T(zjucad::matrix::colon(), tri_), hes(27, 9);
        if(is_degenerate_tri(tri)) return 0;
        calc_tri_area_normal_hes_(&hes[0], &tri[0]);
        zjucad::matrix::matrix<double> H(9, 9);
        for(size_t ci = 0; ci < 9; ++ci) {
            const size_t var_ci = tri_[ci/3]*3+ci%3;
            for(size_t ri = 0; ri < 9; ++ri) {
                const size_t var_ri = tri_[ri/3]*3+ri%3;
                H(ri, ci) = hes(xyz_+ri*3, ci)*weight_*alpha;
              }
          }
        const matrix<double> save_H = H;
        matrix<double> e(9), diag_e = zeros<double>(9, 9);
        zjucad::matrix::eig(H, e);
        for(size_t ei = 0; ei < 9; ++ei) {
            if(e[ei] > 0)
              diag_e(ei, ei) = e[ei];
            else
              diag_e(ei, ei) = 0;
          }
        H = temp(H*temp(diag_e*trans(H)));
        for(size_t ci = 0; ci < 9; ++ci) {
            const size_t var_ci = tri_[ci/3]*3+ci%3;
            for(size_t ri = 0; ri < 9; ++ri) {
                const size_t var_ri = tri_[ri/3]*3+ri%3;
                if(jtf::function::add_to_csc(h, ptr, idx, var_ri, var_ci, H(ri, ci))) {
                    return __LINE__;
                  }
              }
          }
        return 0;
      }
    return __LINE__;
  }
  // accumulate
  virtual int hes_block(const double *x, double *h, double alpha = 1) {
    using namespace zjucad::matrix;
    const zjucad::matrix::itr_matrix<const double *> T(3, node_num_, x);
    zjucad::matrix::matrix<double> tri =
        T(zjucad::matrix::colon(), tri_), hes(27, 9);
    if(is_degenerate_tri(tri)) return 0;
    calc_tri_area_normal_hes_(&hes[0], &tri[0]);
    itr_matrix<double *> H(9, 9, h);
    for(size_t ci = 0; ci < 9; ++ci) {
        const size_t var_ci = tri_[ci/3]*3+ci%3;
        for(size_t ri = 0; ri < 9; ++ri) {
            const size_t var_ri = tri_[ri/3]*3+ri%3;
            H(ri, ci) += hes(xyz_+ri*3, ci)*weight_*alpha;
          }
      }
    return 0;
  }
private:
  const size_t node_num_, xyz_;
  const zjucad::matrix::matrix<size_t> tri_;
  const double weight_;
};

#endif // TRIAREANORMALFUNC_H
