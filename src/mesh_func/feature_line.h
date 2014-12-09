#ifndef FEATURE_LINE_H
#define FEATURE_LINE_H

#include <hjlib/math_func/math_func.h>
#include <zjucad/matrix/matrix.h>
#include <zjucad/matrix/itr_matrix.h>
#include <zjucad/matrix/io.h>
#include <hjlib/math_func/operation.h>
#include <jtflib/mesh/mesh.h>
#include <jtflib/mesh/util.h>
#include <vector>
#include <deque>

template <typename val_type, typename int_type>
class feature_align_func: public hj::math_func::math_func_t<val_type,int_type>
{
public:
  feature_align_func(const size_t point_num, const size_t point_idx,
                     const zjucad::matrix::matrix<val_type> & point,
                     const zjucad::matrix::matrix<val_type> & dir,
                     const val_type weight):point_num_(point_num),
    point_idx_(point_idx), dir_(dir), weight_(weight),point_node_(point){}

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
        // cross(\delta x, dir)
        matrix<val_type> delta_x = T(colon(), point_idx_) - point_node_;
        matrix<val_type> cross_dir = cross(delta_x, dir_);

        for(int_type i = 0; i < 3; ++i){
            int_type c[] = {i};
            cv[c] += weight_ * cross_dir[i];
          }
      }
    if(k == 1){
        // cross_dir = (x1d2-x2d1)i+(x2d0-x0d2)i+(x0d1-x1d0)k
        for(int_type i = 0; i < 3; ++i){
            for(int_type j = 1; j < 3; ++j){
                int_type c[] = {i, int_type(3*point_idx_+(i+j)%3)};
                cv[c] += weight_* dir_[(i+3-j)%3]*(j==1?1.0:-1.0);
              }
          }
      }
    return 0;
  }
  virtual int patt(size_t k, hj::math_func::coo_set<int_type> &cs,
                   const hj::math_func::coo_l2g &l2g,
                   hj::math_func::func_ctx *ctx = 0) const{
    if(k == 1){
        // cross_dir = (x1d2-x2d1)i+(x2d0-x0d2)i+(x0d1-x1d0)k
        for(int_type i = 0; i < 3; ++i){
            for(int_type j = 1; j < 3; ++j){
                int_type c[] = {i, int_type(3*point_idx_+(i+j)%3)};
                l2g.add(cs,c);
              }
          }
      }
    return 0;
  }
  virtual size_t nnz(size_t k) const{
    if(k == 0) return -1;
    if(k == 1) return 3*2;
  }

private:
  const size_t point_num_;
  const size_t point_idx_;
  const zjucad::matrix::matrix<val_type> dir_;
  const zjucad::matrix::matrix<val_type> point_node_;
  const double weight_;
};


template <typename val_type>
std::shared_ptr<const hj::math_func::math_func_t<val_type, int32_t> >
build_feature_line_math_func(const std::vector<std::deque<std::pair<size_t,size_t> > > & fl,
                             const zjucad::matrix::matrix<val_type> &node,
                             const val_type weight)
{
  using namespace std;
  using namespace zjucad::matrix;
  typedef hj::math_func::math_func_t<val_type, int32_t> math_func_type;
  typedef std::shared_ptr<math_func_type> math_func_ptr;
  std::shared_ptr<std::vector<math_func_ptr> > all_func(new std::vector<math_func_ptr>);

  matrix<val_type> dir(3,1);
  for(size_t fi = 0; fi < fl.size(); ++fi){
      const deque<pair<size_t,size_t> > & one_line = fl[fi];
      for(const auto & one_edge : one_line){
          dir = node(colon(), one_edge.second) - node(colon(), one_edge.first);
          all_func->push_back(
                math_func_ptr(
                  new feature_align_func<val_type, int32_t>(
                    node.size(2), one_edge.first, node(colon(), one_edge.first), dir, sqrt(weight)
                    )));
        }
    }
  std::cerr << "# [info] smooth point function number " << all_func->size() << std::endl;
  math_func_ptr fun_cat(new hj::math_func::fcat<double,int32_t, std::vector<math_func_ptr> >(all_func));
  return fun_cat;
}


#endif // FEATURE_LINE_H
