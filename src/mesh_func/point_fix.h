#ifndef POINT_FIX_H
#define POINT_FIX_H

#include <hjlib/math_func/math_func.h>
#include <zjucad/matrix/matrix.h>
#include <zjucad/matrix/itr_matrix.h>
#include <hjlib/math_func/operation.h>

template <typename val_type, typename int_type>
class point_fix_math_func: public hj::math_func::math_func_t<val_type,int_type>
{
public:
  point_fix_math_func(const size_t point_num, const size_t point_idx,
                      const zjucad::matrix::matrix<val_type> &point_target,
                      const val_type weight):point_num_(point_num),
    point_idx_(point_idx), point_target_(point_target), weight_(weight){}

  virtual size_t nx(void) const{
    return 3*point_num_;
  }
  virtual size_t nf(void) const{
    return 3;
  }
  virtual int eval(size_t k, const val_type *x,
                   const hj::math_func::coo2val_t<val_type, int_type> &cv,
                   hj::math_func::func_ctx *ctx = 0) const{
    using namespace zjucad::matrix;
    const itr_matrix<const val_type*> T(3, point_num_, x);
    if(k == 0){
        for(int_type i = 0; i < 3; ++i){
            int_type c[] = {i};
            cv[c] += weight_*(T(i,point_idx_) - point_target_[i]);
          }
      }
    if(k == 1){
        for(int_type i = 0; i < 3; ++i){
            int_type c[] = {i, int_type(3*point_idx_+i)};
            cv[c] += weight_;
          }
      }
    return 0;
  }
  virtual int patt(size_t k, hj::math_func::coo_set<int_type> &cs,
                   const hj::math_func::coo_l2g &l2g,
                   hj::math_func::func_ctx *ctx = 0) const{
    if(k == 1){
        for(int_type i = 0; i < 3; ++i){
            int_type c[] = {i, int_type(3*point_idx_+i)};
            l2g.add(cs,c);
          }
      }
    return 0;
  }
  virtual size_t nnz(size_t k) const{
    if(k == 0) return -1;
    if(k == 1) return 3;
  }

private:
  const size_t point_num_;
  const size_t point_idx_;
  const zjucad::matrix::matrix<val_type> point_target_;
  const double weight_;
};

template <typename val_type>
std::shared_ptr<const hj::math_func::math_func_t<val_type,int32_t> >
build_point_fix_math_func(const zjucad::matrix::matrix<val_type> & target,
                          const zjucad::matrix::matrix<size_t> & point_idx,
                          const size_t point_num,
                          const val_type weight)
{
  using namespace zjucad::matrix;
  typedef hj::math_func::math_func_t<val_type, int32_t> math_func_type;
  typedef std::shared_ptr<math_func_type> math_func_ptr;
  std::shared_ptr<std::vector<math_func_ptr> > all_func(new std::vector<math_func_ptr>);

  for(size_t it = 0; it != point_idx.size(); ++it){
      all_func->push_back(
            math_func_ptr(
              new point_fix_math_func<val_type,int32_t>(
                point_num, point_idx[it], target(colon(),point_idx[it]), sqrt(weight))));
    }

  std::cerr << "# [info] fix point function number " << all_func->size() << std::endl;
  math_func_ptr fun_cat(new hj::math_func::fcat<val_type,int32_t, std::vector<math_func_ptr> >(all_func));
  return fun_cat;
}

template <typename val_type>
std::shared_ptr<const hj::math_func::math_func_t<val_type,int32_t> >
build_point_fix_math_func(const zjucad::matrix::matrix<val_type> & target,
                          const val_type weight)
{
  using namespace zjucad::matrix;
  typedef hj::math_func::math_func_t<val_type, int32_t> math_func_type;
  typedef std::shared_ptr<math_func_type> math_func_ptr;
  std::shared_ptr<std::vector<math_func_ptr> > all_func(new std::vector<math_func_ptr>);

  for(size_t it = 0; it != target.size(2); ++it){
      all_func->push_back(
            math_func_ptr(
              new point_fix_math_func<double,int32_t>(target.size(2), it, target(colon(),it), sqrt(weight))));
    }

  std::cerr << "# [info] fix point function number " << all_func->size() << std::endl;
  math_func_ptr fun_cat(new hj::math_func::fcat<double,int32_t, std::vector<math_func_ptr> >(all_func));
  return fun_cat;
}


#endif // POINT_FIX_H
