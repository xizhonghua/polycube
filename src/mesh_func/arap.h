#ifndef MESH_FUNC_ARAP_H
#define MESH_FUNC_ARAP_H

#include <hjlib/math_func/math_func.h>
#include <zjucad/matrix/matrix.h>
#include <zjucad/matrix/itr_matrix.h>
#include <hjlib/math/polar.h>

int calc_tet_def_grad_op(const double *tet, double *grad_op);


template <typename val_type, typename int_type>
class arap_math_func: public hj::math_func::math_func_t<val_type,int_type>
{
public:
  arap_math_func(const zjucad::matrix::matrix<size_t> & tet,
                 const zjucad::matrix::matrix<double> & node,
                 const size_t idx, const double weight)
    :grad_op_(4,3), one_tet_points_(tet(zjucad::matrix::colon(), idx)),
      node_num_(node.size(2)), idx_(idx), weight_(weight){
    using namespace zjucad::matrix;
    matrix<val_type> tet_node = node(colon(), one_tet_points_);
    if(calc_tet_def_grad_op(&tet_node[0], &grad_op_[0])){
        std::cerr << "# [error] degenerated tet."  << std::endl;
      }
  }

  virtual size_t nx(void) const{
    return node_num_*3;
  }
  virtual size_t nf(void) const{
    return 3*3;
  }
  virtual int eval(size_t k, const val_type *x,
                   const hj::math_func::coo2val_t<val_type, int_type> &cv,
                   hj::math_func::func_ctx *ctx = 0) const{
    using namespace zjucad::matrix;
    if(k == 0){
        matrix<val_type> f0(3,3);
        const itr_matrix<const val_type*> T(3, node_num_, x);
        f0 = T(colon(), one_tet_points_) * grad_op_;
        matrix<val_type> R = f0;
        hj::polar3d p;
        p(R,2);
        f0 -= R;
        f0 *= weight_;
        for(int_type i = 0; i < 9; ++i){
            int_type c[] = {i};
            cv[c] += f0[i];
          }
      }
    if(k == 1){
        for(int_type c = 0; c < 3; ++c){
            for(int_type r = 0; r < 3; ++r){
                const int_type fi = c*3+r;
                for(int_type k = 0; k < 4; ++k){
                    int_type c1[2] = {fi, static_cast<int_type>(3*one_tet_points_[k]+r)};
                    cv[c1] += grad_op_(k,c) * weight_;
                  }
              }
          }
      }
    return 0;
  }
  virtual int patt(size_t k, hj::math_func::coo_set<int_type> &cs,
                   const hj::math_func::coo_l2g &l2g,
                   hj::math_func::func_ctx *ctx = 0) const{
    if(k == 1){
        for(int_type c = 0; c < 3; ++c){
            for(int_type r = 0; r < 3; ++r){
                const int_type fi = c*3+r;
                for(int_type k = 0; k < 4; ++k){
                    int_type c1[2] = {fi, static_cast<int_type>(3*one_tet_points_[k]+r)};
                    l2g.add(cs, c1);
                  }
              }
          }
      }
    return 0;
  }
  virtual size_t nnz(size_t k) const{
    if(k == 0) return -1;
    if(k == 1) return 9*4;
  }

private:
  zjucad::matrix::matrix<double> grad_op_;
  const zjucad::matrix::matrix<size_t> one_tet_points_;
  const int_type idx_;
  const int_type node_num_;
  const val_type weight_;
};

std::shared_ptr<hj::math_func::math_func_t<double,int32_t> >
build_arap_math_func(const zjucad::matrix::matrix<size_t> & tet,
                     const zjucad::matrix::matrix<double> & node,
                     const double weight);

#endif // MESH_FUNC_ARAP_H
