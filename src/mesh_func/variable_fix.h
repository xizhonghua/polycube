#ifndef VARIABLE_FIX_H
#define VARIABLE_FIX_H

#include <hjlib/math_func/math_func.h>
#include <zjucad/matrix/matrix.h>
#include <zjucad/matrix/itr_matrix.h>
#include <hjlib/math_func/operation.h>

template <typename val_type, typename int_type>
class variable_fix_math_func: public hj::math_func::math_func_t<val_type,int_type>
{
public:
  variable_fix_math_func(const int_type v_num, const int_type idx,
                         const val_type target, const val_type weight)
    :v_num_(v_num), target_(target), weight_(weight), idx_(idx){}

  virtual size_t nx(void) const{
    return v_num_;
  }
  virtual size_t nf(void) const{
    return 1;
  }
  virtual int eval(size_t k, const val_type *x,
                   const hj::math_func::coo2val_t<val_type, int_type> &cv,
                   hj::math_func::func_ctx *ctx = 0) const{
    using namespace zjucad::matrix;
    if(k == 0){
        int_type c[] = {0};
        cv[c] += weight_*(x[idx_] - target_);
      }
    if(k == 1){
        int_type c[] = {0, idx_};
        cv[c] += weight_;
      }
    return 0;
  }
  virtual int patt(size_t k, hj::math_func::coo_set<int_type> &cs,
                   const hj::math_func::coo_l2g &l2g,
                   hj::math_func::func_ctx *ctx = 0) const{
    if(k == 1){
        int_type c[] = {0, idx_};
        l2g.add(cs,c);
      }
    return 0;
  }
  virtual size_t nnz(size_t k) const{
    if(k == 0) return -1;
    if(k == 1) return 1;
  }

private:
  const int_type idx_;
  const int_type v_num_;
  const val_type target_;
  const val_type weight_;
};

template <typename it_type>
std::shared_ptr<const hj::math_func::math_func_t<double,int32_t> >
build_variable_fix_math_func(const size_t v_num, it_type begin, it_type end,
                             const double target, const double weight)
{
  typedef hj::math_func::math_func_t<double, int32_t> math_func_type;
  typedef std::shared_ptr<math_func_type> math_func_ptr;
  std::shared_ptr<std::vector<math_func_ptr> > all_func(new std::vector<math_func_ptr>);

  for(auto it = begin; it != end; ++it){
       all_func->push_back(math_func_ptr(new variable_fix_math_func<double,int32_t>(v_num, *it, target, sqrt(weight))));
    }

  std::cerr << "# [info] fix variable function number " << all_func->size() << std::endl;
  math_func_ptr fun_cat(new hj::math_func::fcat<double,int32_t, std::vector<math_func_ptr> >(all_func));
  return fun_cat;
}

#endif // VARIABLE_FIX_H
