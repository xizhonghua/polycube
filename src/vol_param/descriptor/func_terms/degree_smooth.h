#ifndef DEGREE_SMOOTH_H
#define DEGREE_SMOOTH_H

#include "../def.h"

#include <zjucad/matrix/matrix.h>
#include <zjucad/matrix/itr_matrix.h>
#include <hjlib/sparse/sparse.h>

#include <map>

class degree_smooth_func : public hj_func
{
public:
  degree_smooth_func(const size_t point_num,
                     const size_t point_idx,
                     const std::vector<size_t> & one_ring)
    :point_num_(point_num), point_idx_(point_idx), one_ring_(one_ring), dim_(3){}
  virtual ~degree_smooth_func(){}

public:
  virtual size_t dim_of_x(void) const {
    return point_num_ * dim_;
  }
  virtual size_t dim_of_f(void) const {
    return dim_;
  }
  virtual int val(const double *x, double *f,
                  hj::function::func_ctx *ctx = 0)const{
    assert(x);
    using namespace zjucad::matrix;
    zjucad::matrix::itr_matrix<double* > f_m(dim_,1,f);
    zjucad::matrix::matrix<double> one_ring_node;
    zjucad::matrix::itr_matrix<const double *> x_m(dim_, point_num_, x);
    itr_matrix<const size_t*> one_ring_mat(one_ring_.size(),1, &one_ring_[0]);
    one_ring_node = x_m(zjucad::matrix::colon(), one_ring_mat);

    f_m  *= 0;
    for(size_t i = 0; i < one_ring_.size(); ++i){
        f_m -= one_ring_node(colon(),i);
      }
    zjucad::matrix::matrix<double> point_(dim_,1);

    point_ = x_m(zjucad::matrix::colon(), point_idx_);

    f_m /= one_ring_.size();
    f_m += point_;
    return 0;
  }
  virtual int jac(const double *x, double *val, int32_t *ptr = 0,
                  int32_t *idx = 0, hj::function::func_ctx *ctx = 0) const {

    for(size_t di = 0; di < dim_; ++di){
        ptr[di+1] = ptr[di] + one_ring_.size() + 1;
        for(size_t i = 0; i < one_ring_.size(); ++i){
            idx[ptr[di] + i] = dim_*one_ring_[i]+di;
            val[ptr[di] + i] = -1.0/one_ring_.size();
          }
        idx[ptr[di] + one_ring_.size()] = dim_*point_idx_+di;
        val[ptr[di] + one_ring_.size()] = 1.0;
      }
    return 0;
  }
  virtual size_t jac_nnz(void) const {
    return (one_ring_.size()+1)*dim_of_f();
  }
private:
  const std::vector<size_t> one_ring_;
  const size_t point_idx_;
  const size_t point_num_;
  const size_t dim_;
};


template <typename T1>
hj_func * build_degree_smooth_func(
    const zjucad::matrix::matrix_expression<T1> & mesh)
{
  using namespace std;
  using namespace zjucad::matrix;

  std::map<size_t,set<size_t> > point2neighbour;
  for(size_t ci = 0; ci < mesh().size(2); ++ci){
      for(size_t pi = 0; pi < mesh().size(1); ++pi){
          for(size_t pj = 1; pj < mesh().size(1); ++pj)
            point2neighbour[mesh()(pi,ci)].insert(mesh()((pj+pi)%mesh().size(1),ci));
        }
    }

  std::shared_ptr<vector<std::shared_ptr<const hj::function::function_t<double,int32_t> > > >
      funcs(new vector<std::shared_ptr<const hj::function::function_t<double,int32_t> > >);

  vector<size_t> neighbour;
  for(const auto & pi2n : point2neighbour){
      neighbour.resize(pi2n.second.size());
      std::copy(pi2n.second.begin(), pi2n.second.end(), neighbour.begin());
      funcs->push_back(hj_func_cons_ptr(new degree_smooth_func(
                                          point2neighbour.size(), pi2n.first,neighbour)));
    }

  return hj::function::new_catenated_function<double, int32_t>(funcs);
}

#endif // DEGREE_SMOOTH_H
