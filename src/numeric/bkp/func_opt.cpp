#include "func_opt.h"
#include "SQP.h"

#include <limits>
#include <iostream>
#include <set>
#include <zjucad/matrix/itr_matrix.h>
#include <zjucad/matrix/matrix.h>
#include <zjucad/matrix/io.h>

#include <hjlib/sparse/cached_AAT.h>
#include <hjlib/sparse/fast_AAT.h>
#include <jtflib/math/math.h>

#include "util.h"

using namespace std;
using namespace zjucad::matrix;

namespace hj_func_opt {

sum_function::sum_function(const container children)
  :children_(children), nnz_(-1)
{
  const size_t fn = children_.size();
  assert(fn);
  for(size_t i = 1; i < fn; ++i) {
    if(children_[i]->dim() != children_[0]->dim()) {
      cerr << "incompatible functions." << children_[0]->dim()
           << " " << children_[i]->dim() << endl;
    }
  }
}

size_t sum_function::dim(void) const
{
  return children_[0]->dim();
}

int sum_function::val(const double *x, double &v)
{
  const size_t fn = children_.size();
  for(size_t i = 0; i < fn; ++i) {
    double tmp = 0;
    children_[i]->val(x, v);
  }
  return 0;
}

int sum_function::gra(const double *x, double *g)
{
  const size_t fn = children_.size();
  for(size_t i = 0; i < fn; ++i) {
    children_[i]->gra(x, g);
  }
  return 0;
}

int sum_function::hes(const double *x, size_t &nnz, size_t &format, double *h, int32_t *ptr, int32_t *idx, double alpha)
{
  format = 1;
  if(h == 0 && ptr == 0 && idx == 0) {// query nnz
    pattern_.resize(dim());
    pair<vector<int32_t>, vector<int32_t> > ptr_idx;
    for(size_t fi = 0; fi < children_.size(); ++fi) {
      size_t nnz0, format = -1;
      if(children_[fi]->hes(x, nnz0, format, 0, 0, 0))
        return __LINE__;
      if(format == 1) { // csc
        ptr_idx.first.resize(dim()+1);
        ptr_idx.first[0] = 0;
        ptr_idx.second.resize(nnz0);
        if(children_[fi]->hes(x, nnz0, format, 0, &ptr_idx.first[0], &ptr_idx.second[0]))
          return __LINE__;
        for(size_t ci = 0; ci < dim(); ++ci) {
          for(size_t nzi = ptr_idx.first[ci]; nzi < ptr_idx.first[ci+1]; ++nzi) {
            pattern_[ci].insert(ptr_idx.second[nzi]);
          }
        }
      }
      else if(format == 2) {// pair
        ptr_idx.first.resize(nnz0);
        ptr_idx.second.resize(nnz0);
        if(children_[fi]->hes(x, nnz0, format, 0, &ptr_idx.first[0], &ptr_idx.second[0]))
          return __LINE__;
        for(size_t nzi = 0; nzi < nnz0; ++nzi) {
          pattern_[ptr_idx.first[nzi]].insert(ptr_idx.second[nzi]);
        }
      }
    }
    nnz = 0;
    for(size_t xi = 0; xi < dim(); ++xi)
      nnz += pattern_[xi].size();
    nnz_ = nnz;
    return 0;
  }
  if(h == 0 && ptr != 0 && idx != 0) {// query patten
    if(nnz < nnz_) {
      cerr << "incorrect input at query pattern: " << nnz << " " << nnz_;
      return __LINE__;
    }
    for(size_t xi = 0; xi < dim(); ++xi) {
      ptr[xi+1] = ptr[xi] + pattern_[xi].size();
      size_t nzi = ptr[xi];
      for(set<int32_t>::const_iterator i = pattern_[xi].begin();
          i != pattern_[xi].end(); ++i, ++nzi) {
        idx[nzi] = *i;
      }
    }
    std::vector<std::set<int32_t> > tmp;
    swap(pattern_, tmp);
    return 0;
  }
  if(h != 0 && ptr != 0 && idx != 0) {// accumulate
    if(nnz < nnz_ && nnz_ != -1) { // when nnz_ == -1, client know the pattern already
      cerr << "incorrect input at accumulate: " << nnz << " " << nnz_;
      return __LINE__;
    }
    size_t format = -1;
    for(size_t fi = 0; fi < children_.size(); ++fi)
      if(children_[fi]->hes(x, nnz, format, h, ptr, idx, alpha))
        return __LINE__;
    return 0;
  }
  return __LINE__;
}

bool sum_function::is_valid(const double *x) const
{
  for(size_t fi = 0; fi < children_.size(); ++fi)
    if(!children_[fi]->is_valid(x))
      return false;
  return true;
}

neg_log_sum_func::neg_log_sum_func(src_func_t src, const zjucad::matrix::matrix<double> &weight)
  :src_(src), weight_(weight), is_JT_sorted_(false)
{
  assert(src_.get());
  src_val_.resize(src_->dim_of_f());
  src_jac_.resize(src_->dim_of_x(), src_->dim_of_f(), src_->jac_nnz());
}

size_t neg_log_sum_func::dim(void) const
{
  return src_->dim_of_x();
}

const double NEG_LOG_MIN = 1e-2;
int neg_log_sum_func::val(const double *x, double &v)
{
  src_->val(x, &src_val_[0]);
  for(size_t fi = 0; fi < src_val_.size(); ++fi) {
    if(!is_interior_point(src_val_[fi])) continue;
    //      src_val_[fi] = NEG_LOG_MIN;
    double t = -log(src_val_[fi])*weight_[fi];
    jtf::math::erase_nan_inf(t);
    v += t;
  }
  return 0;
}

int neg_log_sum_func::gra(const double *x, double *g)
{
  src_->val(x, &src_val_[0]);
  src_->jac(x, &src_jac_.val()[0], &src_jac_.ptr()[0], &src_jac_.idx()[0]);
  for(size_t fi = 0; fi < src_jac_.size(2); ++fi) {
    if(!is_interior_point(src_val_[fi])) continue;
    //      src_val_[fi] = NEG_LOG_MIN;
    for(size_t nzi = src_jac_.ptr()[fi]; nzi < src_jac_.ptr()[fi+1]; ++nzi) {
      double acc = -src_jac_.val()[nzi]/src_val_[fi]*weight_[fi];
      jtf::math::erase_nan_inf(acc);
      g[src_jac_.idx()[nzi]] += acc;
    }
  }
  return 0;
}

int neg_log_sum_func::hes(const double *x, size_t &nnz, size_t &format, double *h, int32_t *ptr, int32_t *idx, double alpha)
{
  format = 1;
  src_->jac(x, &src_jac_.val()[0], &src_jac_.ptr()[0], &src_jac_.idx()[0]);
	if(hj::sparse::nnz(src_hes_))
		fast_AAT(src_jac_, src_hes_, is_JT_sorted_);
	else {
		is_JT_sorted_ = is_sorted_csc(src_jac_);
    hj::sparse::AAT<hj::sparse::map_by_sorted_vector>(src_jac_, src_hes_);
	}
  if(h == 0 && ptr == 0 && idx == 0) { // query nnz
    nnz = hj::sparse::nnz(src_hes_);
    return 0;
  }
  if(h == 0 && ptr !=0 && idx != 0) { // query patten
    if(nnz < hj::sparse::nnz(src_hes_)) {
      cerr << "nnz err: " << nnz << " " << hj::sparse::nnz(src_hes_) << endl;
      return __LINE__;
    }
    copy(src_hes_.ptr().begin(), src_hes_.ptr().end(), ptr);
    copy(src_hes_.idx().begin(), src_hes_.idx().end(), idx);
    return 0;
  }
  if(h != 0 && ptr !=0 && idx != 0) { // accumulate
    src_->val(x, &src_val_[0]);
    matrix<double> src_g, src_H;
    matrix<int32_t> src_g_idx;
    for(size_t fi = 0; fi < src_jac_.size(2); ++fi) {
      if(!is_interior_point(src_val_[fi])) continue;
      //        src_val_[fi] = NEG_LOG_MIN;
      src_g_idx = src_jac_.idx()(colon(src_jac_.ptr()[fi], src_jac_.ptr()[fi+1]-1));
      src_g = src_jac_.val()(colon(src_jac_.ptr()[fi], src_jac_.ptr()[fi+1]-1));
      const double scale = alpha*weight_[fi]/(src_val_[fi]*src_val_[fi]);
      src_H = src_g*trans(src_g)*scale;
      for(int ci = 0; ci < src_g.size(); ++ci) {
        for(int ri = 0; ri < src_g.size(); ++ri) {
          double t = src_H(ri, ci);
           jtf::math::erase_nan_inf(t);
          add_to_csc(h, ptr, idx, src_g_idx[ri], src_g_idx[ci], t);
        }
      }
    }
    return 0;
  }
  return __LINE__;
}

bool neg_log_sum_func::is_interior_point(double v) const
{
  //  const double eps = std::numeric_limits<double>::min();
  const double eps = 1e-8;
  if(v < eps) {
    //    cerr << "# not interior of neg_log_sum_func: " << v << endl;
    return false;
  }
  return true;
}

bool neg_log_sum_func::is_valid(const double *x) const
{
  double *srv_val_non_cost = const_cast<double *>(&src_val_[0]);
  src_->val(x, srv_val_non_cost);
  for(size_t fi = 0; fi < src_val_.size(); ++fi) {
    //    if(fabs(weight_[fi]) < 1e-5) continue;
    if(!is_interior_point(src_val_[fi])) return false;
  }
  return true;
}

double gra_err(function &f, double *x)
{
  const size_t dim = f.dim();
  double val = 0;
  f.val(x, val);
  matrix<double> g = zeros<double>(dim, 1);
  f.gra(x, &g[0]);
  cerr << "max g: " << max(fabs(g)) << endl;
  const double eps = 1e-6;
  for(size_t xi = 0; xi < dim; ++xi) {
    const double save = x[xi];
    double v[2] = {0, 0};
    x[xi] = save-eps;
    f.val(x, v[0]);
    x[xi] = save+eps;
    f.val(x, v[1]);
    g[xi] -= (v[1]-v[0])/(2*eps);
    x[xi] = save;
  }
  return max(fabs(g));
}

// assume grad is accurate
double hes_err(function &f, double *x)
{
  const size_t dim = f.dim();
  matrix<double> hes(dim, dim);
  size_t nnz, format = -1;
  f.hes(x, nnz, format, 0, 0, 0);
  matrix<double> h = zeros<double>(nnz, 1);
  matrix<int32_t> ptr(dim+1), idx(nnz);
  f.hes(x, nnz, format, 0, &ptr[0], &idx[0]);
  f.hes(x, nnz, format, &h[0], &ptr[0], &idx[0]);
  for(size_t ci = 0; ci < dim; ++ci) {
    for(size_t nzi = ptr[ci]; nzi < ptr[ci+1]; ++nzi)
      hes(idx[nzi], ci) = h[nzi];
  }
  cout << "max hes: " << max(fabs(hes)) << endl;

  const double eps = 1e-6;
  matrix<double> g0 = zeros<double>(dim, 1), ga = zeros<double>(dim, 1), gb = zeros<double>(dim, 1);;
  f.gra(x, &g0[0]);
  for(size_t xi = 0; xi < dim; ++xi) {
    const double x0 = x[xi];

    x[xi] = x0+eps;
    ga(colon()) = 0;
    f.gra(x, &ga[0]);

    x[xi] = x0-eps;
    gb(colon()) = 0;
    f.gra(x, &gb[0]);

    hes(colon(), xi) -= (ga-gb)/(2*eps);

    x[xi] = x0;
  }
  return max(fabs(hes));
}

opt* opt::create(const char *name)
{
  if(!strcmp(name, "SQP")) {
    return new SQP;
  }
  return 0;
}

}
