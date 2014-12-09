#ifndef HJ_FUNC_OPT_NUMERIC_H_
#define HJ_FUNC_OPT_NUMERIC_H_

#include <iostream>
#include <vector>
#include <set>

#include <boost/property_tree/ptree.hpp>
#include <boost/shared_ptr.hpp>
#include <hjlib/function/function.h>
#include <hjlib/sparse/sparse.h>

namespace hj_func_opt {

inline int add_to_csc(double *val, const int32_t *ptr, const int32_t *idx,
                      int32_t r, int32_t c, double v) {
  size_t off;
  for(off = ptr[c]; off < ptr[c+1]; ++off) {
    if(idx[off] == r)
      break;
  }
  if(off == ptr[c+1]) {
    std::cerr << "# incorrect sparse pattern in add_to_csc: "
              << r << ":" << c;
    for(off = ptr[c]; off < ptr[c+1]; ++off) {
      std::cerr << " " << idx[off];
    }
    std::cerr << std::endl;
    return __LINE__;
  }
  val[off] += v;
  return 0;
}

//! @brief: R^n -> R
class function
{
public:
  virtual ~function(){}

  virtual size_t dim(void) const = 0;
  //NOTE: add to v
  virtual int val(const double *x, double &v) = 0;
  //NOTE: add to v
  virtual int gra(const double *x, double *g) = 0;
  virtual int gra(const double *x, size_t &nnz, double *g, int32_t *idx) {return __LINE__;}
  //NOTE: add to hes
  //! @param h,ptr,idx == 0 means query nnz
  //! @param h == 0, ptr,idx != 0 means query patten
  //! @param h,ptr,idx != 0 means accumulate
  //! @param format for query pattern: 1 csc, 2 pair of row, col
  virtual int hes(const double *x, size_t &nnz, size_t &format, double *h, int32_t *ptr, int32_t *idx, double alpha = 1) = 0;
  virtual int hes_block(const double *x, double *h, double alpha=1)  = 0 ;//{ assert(0); return -1; }
  virtual bool is_valid(const double *x) const {return true;}
};

double gra_err(function &f, double *x);
double hes_err(function &f, double *x);

class sum_function : public function
{
public:
  typedef std::vector<std::shared_ptr<function> > container;
  sum_function(const container children);
  virtual size_t dim(void) const;
  virtual int val(const double *x, double &v);
  virtual int gra(const double *x, double *g);
  virtual int hes(const double *x, size_t &nnz, size_t &format, double *h, int32_t *ptr, int32_t *idx, double alpha = 1);
  virtual int hes_block(const double *x, double *h, double alpha = 1) {return -1;}
  virtual bool is_valid(const double *x) const;
protected:
  container children_;
  std::vector<std::set<int32_t> > pattern_;
  size_t nnz_;
};

//! @brief: for inequality, -log(f)
class neg_log_sum_func : public function
{
public:
  typedef std::shared_ptr<hj::function::function_t<double, int32_t> > src_func_t;

  neg_log_sum_func(src_func_t src, const zjucad::matrix::matrix<double> &weight);

  virtual size_t dim(void) const;
  virtual int val(const double *x, double &v);
  virtual int gra(const double *x, double *g);
  virtual int hes(const double *x, size_t &nnz, size_t &format, double *h, int32_t *ptr, int32_t *idx, double alpha = 1);
  virtual int hes_block(const double *x, double *h, double alpha = 1) {return -1;}
  virtual bool is_valid(const double *x) const;
protected:
  bool is_interior_point(double v) const;
  src_func_t src_;
  std::vector<double> src_val_;
  hj::sparse::csc<double, int32_t> src_jac_, src_hes_;
  bool is_JT_sorted_;
  const zjucad::matrix::matrix<double> weight_;
};

class opt
{
public:
  static opt* create(const char *name);

  opt():cb_(0){}
  virtual ~opt(){}

  class callbacks
  {
  public:
    virtual ~callbacks(){}
    //! @return non-zero means break the optimization
    virtual int at_point(const double *x) = 0;
  };
  virtual void set_callbacks(callbacks *cb) { cb_ = cb; }

  virtual int set_f(function &f) = 0;
  virtual int set_c(function *c) {return 1;}
  virtual int solve(double *x, boost::property_tree::ptree &pt) = 0;
protected:
  callbacks *cb_;
};

}

#endif
