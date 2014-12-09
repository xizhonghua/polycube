#ifndef HJ_SMOOTH_L1_H_
#define HJ_SMOOTH_L1_H_

#include <vector>
#include <set>
#include <map>
#include <iostream>

#include <boost/shared_ptr.hpp>
#include <zjucad/matrix/matrix.h>
#include <zjucad/matrix/itr_matrix.h>
#include <zjucad/matrix/io.h>
#include <hjlib/math/blas_lapack.h>
#include <zjucad/matrix/lapack.h>
#include <jtflib/math/math.h>
#include <jtflib/function/function.h>
#include <jtflib/function/func_aux.h>

using namespace zjucad::matrix;
using namespace std;

#include "util.h"

#define USE_SQRT 1

class smooth_L1 : public jtf::function::functionN1_t<double,int32_t>
{
public:
  // NOTICE: assume src provide sparse gradient and hessian
  smooth_L1(std::shared_ptr<jtf::function::functionN1_t<double,int32_t> > &src, double w, double eps=1e-2)
    :src_(src), w_(w), eps2_(eps*eps), eps_(eps) {
  }
  virtual size_t dim(void) const {
    return src_->dim();
  }
  virtual int val(const double *x, double &v) {
    double t = 0;
    src_->val(x, t);

#if USE_SQRT
    v += sqrt(t*t+eps2_)*w_;
#else
    if(fabs(t) > eps_)
      v += (fabs(t)-eps_/2)*w_;
    else
      v += t*t/(2*eps_)*w_;
#endif
    //v += t;
    return 0;
  }

  virtual int gra(const double *x, size_t & nnz, double * val, int32_t * idx)
  {
    if(val == 0 && idx == 0){
        src_->gra(x, nnz, 0, 0);
      }else{
        eval_g(x);
        std:copy(g_.begin(), g_.end(), val);
        std::copy(g_idx_.begin(), g_idx_.end(), idx);
      }
    return 0;
  }

  //! @param nnz=0 means dense, g=0 or h=0 means query nnz
  virtual int gra(const double *x, double *g) {
    double t = 0;
    src_->val(x, t);
    //const double L1_x = t/sqrt(t*t+eps2_);
    double L1_x = t/sqrt(t*t+eps2_);
		jtf::math::erase_nan_inf(L1_x);
    eval_g(x);
    for(size_t xi = 0; xi < g_idx_.size(); ++xi) {
#if USE_SQRT
      g[g_idx_[xi]] += g_[xi]*L1_x*w_;
#else
      if(fabs(t) > eps_) {
        g[g_idx_[xi]] += g_[xi]*((t>0)?1:-1)*w_;
      }
      else {
        g[g_idx_[xi]] += g_[xi]*(t/eps_)*w_;
      }
#endif
    }
    //g[g_idx_[xi]] += g_[xi];
    return 0;
  }
  virtual int hes(const double *x, size_t &nnz, size_t &format, double *h, int32_t *ptr, int32_t *idx, double alpha) {
    format = 2;
    eval_g(x);
    if(h == 0 && ptr == 0 && idx == 0) { // query nnz
      std::vector<int32_t> ptr_, idx_;
      size_t format = -1;
      src_->hes(x, nnz, format, 0, 0, 0);
      ptr_.resize(dim()+1);
      ptr_[0] = 0;
      idx_.resize(nnz);
      h_.resize(nnz);
      src_->hes(x, nnz, format, 0, &ptr_[0], &idx_[0]);

      if(nnz > g_idx_.size()*g_idx_.size()) {
        // TODO: this acceleration need seriuos check, assuming ggt cover all the entries
        for(size_t xi = 0; xi < dim(); ++xi) {
          for(size_t nzi = ptr_[xi]; nzi < ptr_[xi+1]; ++nzi) {
            pattern_[xi].insert(idx_[nzi]);
          }
        }
      }
      for(size_t gci = 0; gci < g_idx_.size(); ++gci) {
        for(size_t gri = 0; gri < g_idx_.size(); ++gri) {
          pattern_[g_idx_[gci]].insert(g_idx_[gri]);
        }
      }
      nnz = 0;
      for(std::map<int32_t, std::set<int32_t> >::const_iterator xi = pattern_.begin();
          xi != pattern_.end(); ++xi) {
        nnz += xi->second.size();
      }
      return 0;
    }
    if(h == 0 && ptr != 0 && idx != 0) { // query pattern
      size_t nzi = 0;
      for(std::map<int32_t, std::set<int32_t> >::const_iterator xi = pattern_.begin();
          xi != pattern_.end(); ++xi) {
        for(std::set<int32_t>::const_iterator ri = xi->second.begin();
            ri != xi->second.end(); ++ri, ++nzi) {
          ptr[nzi] = xi->first;
          idx[nzi] = *ri;
        }
      }

      // size_t ci = 0;
      // for(std::map<int32_t, std::set<int32_t> >::const_iterator xi = pattern_.begin();
      //     xi != pattern_.end(); ++xi) {
      //   for(; ci < xi->first; ++ci) {
      //     ptr[ci+1] = ptr[ci];
      //   }
      //   ptr[ci+1] = ptr[ci]+xi->second.size();
      //   copy(xi->second.begin(), xi->second.end(), idx+ptr[ci]);
      //   ++ci;
      // }
      // for(; ci < dim(); ++ci)
      //   ptr[ci+1] = ptr[ci];

      std::map<int32_t, std::set<int32_t> > tmp;
      swap(pattern_, tmp);
      return 0;
    }
    if(h != 0 && ptr != 0 && idx != 0) { // accumulate
      double t = 0;
      src_->val(x, t);
#if USE_SQRT
      const double L1v = sqrt(t*t+eps2_);
      double L1_x = t/L1v, L1_xx = eps2_/(L1v*L1v*L1v);
			jtf::math::erase_nan_inf(L1_x);
			jtf::math::erase_nan_inf(L1_xx);
      matrix<double> H = (L1_xx*w_)*g_*trans(g_);
      if(src_->hes_block(x, &H[0], L1_x*w_)) {
        cerr << "hes block error." << endl;
        exit(0);
      }

      matrix<double> e(H.size(1)), diag_e = zeros<double>(H.size(1), H.size(1));
      eig(H, e);
      for(size_t ei = 0; ei < e.size(); ++ei) {
        if(e[ei] > 0)
          diag_e(ei, ei) = e[ei];
        else
          diag_e(ei, ei) = 0;
      }
      H = temp(H*temp(diag_e*trans(H)));
      for(size_t gci = 0; gci < g_idx_.size(); ++gci) {
        for(size_t gri = 0; gri < g_idx_.size(); ++gri) {
          if(jtf::function::add_to_csc(h, ptr, idx, g_idx_[gri], g_idx_[gci], H(gci, gri)))
            return __LINE__;
        }
      }
#else
      if(fabs(t) > eps_) {
        src_->hes(x, nnz, &h[0], ptr, idx, ((t>0)?1:-1)*w_);
      }
      else {
        src_->hes(x, nnz, &h[0], ptr, idx, t/eps_*w_);
        for(size_t gci = 0; gci < g_idx_.size(); ++gci) {
          for(size_t gri = 0; gri < g_idx_.size(); ++gri) {
            if(hj_func_opt::add_to_csc(h, ptr, idx, g_idx_[gri], g_idx_[gci], g_[gci]*g_[gri]*w_/eps_))
              return __LINE__;
          }
        }
      }
#endif
      return 0;
    }
    return __LINE__;
  }
  virtual int hes_block(const double *x, double *h, double alpha = 1) {return -1;}
protected:
  int eval_g(const double *x) {
    if(g_.size() == 0) {
      size_t nnz = 0;
      src_->gra(x, nnz, 0, 0);
      g_.resize(nnz, 1);
      g_idx_.resize(nnz);
    }
    fill(g_.begin(), g_.end(), 0);
    size_t nnz = g_.size();
    src_->gra(x, nnz, &g_[0], &g_idx_[0]);
    return 0;
  }
  std::shared_ptr<jtf::function::functionN1_t<double,int32_t> > src_;
  zjucad::matrix::matrix<double> g_, h_;
  std::vector<int32_t> g_idx_;
  std::map<int32_t, std::set<int32_t> > pattern_;
  const double eps2_, w_, eps_;
};

#endif
