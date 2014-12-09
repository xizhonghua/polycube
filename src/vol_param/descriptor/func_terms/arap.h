#ifndef VOL_PARAM_FUNC_ARAP_H
#define VOL_PARAM_FUNC_ARAP_H

#include <vector>
#include <iostream>
#include <set>

#include <zjucad/matrix/matrix.h>
#include <zjucad/matrix/io.h>
#include <zjucad/matrix/itr_matrix.h>
#include <jtflib/mesh/util.h>
#include <jtflib/algorithm/equation.h>
#include <hjlib/math/blas_lapack.h>
#include <zjucad/matrix/lapack.h>

#include <hjlib/math/polar.h>

#include <hjlib/sparse/sparse.h>
#include <hjlib/function/function.h>

#include "../def.h"
#include "../../common/util.h"

///
/// @brief calc_tet_def_grad_op
/// @param tet input tet nodes
/// @param grad_op output gradient operator
/// @return  0 if succeed, or non-zeros
///
inline int calc_tet_def_grad_op(const double *tet, double *grad_op)
{
  using namespace zjucad::matrix;
  const static double node2edge[] = {
    1, 0, 0, -1,
    0, 1, 0, -1,
    0, 0, 1, -1
  };
  const static itr_matrix<const double *> E(4, 3, node2edge);

  itr_matrix<const double *> T(3, 4, tet);
  itr_matrix<double *> G(4, 3, grad_op);
  matrix<double> A = T*E;

  if(inv(A))
    return __LINE__;
  G = E*A;
  return 0;
}

///
/// @brief calc_tet_def_grad_op
/// @param tet input tet nodes
/// @param grad_op output gradient operator
/// @return  0 if succeed, or non-zeros
///
inline int calc_tet_def_grad_size_op(const double *tet, double *grad_op, const double *size_metric)
{
  using namespace zjucad::matrix;
  const static double node2edge[] = {
    1, 0, 0, -1,
    0, 1, 0, -1,
    0, 0, 1, -1
  };
  const static itr_matrix<const double *> E(4, 3, node2edge);
  const itr_matrix<const double *> S(3,3, size_metric);

  itr_matrix<const double *> T(3, 4, tet);
  itr_matrix<double *> G(4, 3, grad_op);
  matrix<double> A = T*E;

  if(inv(A))
    return __LINE__;

  G = E*S*A;
  return 0;
}

// eqn0 += eqn1 * v
inline void add_to_eqn(jtf::algorithm::equation<double> & eqn0,
                       const std::vector<std::pair<size_t,double> > & eqn1,
                       const double & v)
{
  for(size_t exp_i = 0; exp_i < eqn1.size(); ++exp_i){
      eqn0.add_expression(jtf::algorithm::make_expression(eqn1[exp_i].first,
                                                          eqn1[exp_i].second * v));
    }
  eqn0.sort_equation();
  eqn0.merge_similar();
}
template<typename T1, typename T2, SIMPLEX>
class arap_func : public hj::function::function_t<double,int32_t> {};

template<typename T1, typename T2, SIMPLEX>
class asap_func : public hj::function::function_t<double,int32_t> {};

template<typename T1, typename T2>
class arap_func <T1,T2,TET> : public hj::function::function_t<double,int32_t>
{
public:
  arap_func(const zjucad::matrix::matrix_expression<T1> & node,
            const zjucad::matrix::matrix_expression<T2> & tet,
            const size_t idx, const double w,
            const node_mapping * NM = 0)
    :grad_op_(4, 3), tet_(tet()), node_num_(node().size(2)), idx_(idx) , weight_(w),
      node_mapping_(NM), nnz_(-1){
    using namespace zjucad::matrix;
    one_tet_ = tet_(colon(),idx_);

    zjucad::matrix::matrix<double> tet_node = node()(colon(), one_tet_);
    if(calc_tet_def_grad_op(&tet_node[0], &grad_op_[0])) {
        std::cerr << "# degenerated tet" << idx << tet_node << std::endl;
      }
    assemble();
  }
  virtual ~arap_func(){}
public:

  virtual size_t dim_of_x(void) const {
    if(!node_mapping_)
      return node_num_*3;
    else
      return node_mapping_->ZT.size(1);
  }
  virtual size_t dim_of_f(void) const {
    return 9;
  }
  virtual int val(const double *x, double *f,
                  hj::function::func_ctx *ctx = 0) const {
    using namespace zjucad::matrix;
    itr_matrix<double *> f0(3, 3, f);

    if(!node_mapping_){
        const itr_matrix<const double *> T(3, node_num_, x);
        f0 = T(colon(), one_tet_)*grad_op_;
      }else{
        matrix<double> tet_node(3,4);
        get_cell_node(x, one_tet_, *node_mapping_, tet_node);
        f0 = tet_node * grad_op_;
      }

    matrix<double> R = f0;
    hj::polar3d p;
    p(R, 2);
    const double alpha = 1;
    R = alpha*R+(1-alpha)*eye<double>(3);
    f0 -= R;
    f0 *= weight_;

    for(size_t i = 0; i < f0.size(); ++i) jtf::math::erase_nan_inf(f0[i]);
    return 0;
  }
  virtual int jac(const double *x, double *val, int32_t *ptr = 0,
                  int32_t *idx = 0, hj::function::func_ctx *ctx = 0) const {
    if(!node_mapping_){
        for(int c = 0; c < 3; ++c) {
            for(int r = 0; r < 3; ++r) {
                const int fi = c*3+r;
                ptr[fi+1] = ptr[fi] + tet_.size(1);
                for(int k = 0; k < tet_.size(1); ++k) {
                    idx[ptr[fi]+k] = tet_(k, idx_)*3+r;
                    val[ptr[fi]+k] = grad_op_(k, c) * weight_;
                    jtf::math::erase_nan_inf(val[ptr[fi]+k]);
                  }
              }
          }
      }else{
        assert(eqns_of_each_f.size() == 9);
        for(size_t fi = 0; fi < 9; ++fi){
            ptr[fi+1] = ptr[fi] + eqns_of_each_f[fi].e_vec_.size();
            size_t k = ptr[fi];
            for(const auto & exp : eqns_of_each_f[fi].e_vec_){
                assert(k != ptr[fi+1]);
                idx[k] = exp.index;
                val[k] = exp.coefficient * weight_;
                jtf::math::erase_nan_inf(val[k]);
                ++k;
              }
          }
      }
    return 0;
  }
  virtual size_t jac_nnz(void) const {
    if(!node_mapping_)
      return 4*dim_of_f();
    else
      return nnz_;
  }
  void assemble(){
    using  namespace zjucad::matrix;
    if(node_mapping_){
        get_node_mapping(tet_(colon(),idx_), node_mapping_->ZT, node_mapping_of_each_variabe_);
        nnz_ = 0;
        eqns_of_each_f.resize(9); // dim_of_f
        //f(i,j) = \sum_k^4 tn_i,3-k * grad_3-k, j
        for(size_t item = 0; item < 9; ++item){
            const size_t i = item % 3;
            const size_t j = item / 3;
            jtf::algorithm::equation<double> one_eqn;
            for(size_t k = 0; k < 4; ++k){
                add_to_eqn(one_eqn, node_mapping_of_each_variabe_[i + 3*(3-k)],grad_op_(3-k,j));
              }
            nnz_ += one_eqn.e_vec_.size();
            eqns_of_each_f[item] = one_eqn;
          }
      }
  }

private:
  const size_t idx_, node_num_;
  const zjucad::matrix::matrix<size_t> &tet_;
  zjucad::matrix::matrix<double> grad_op_;
  const double weight_;
  const node_mapping *node_mapping_;
  zjucad::matrix::matrix<size_t> one_tet_;
  std::vector<std::vector<std::pair<size_t,double> > > node_mapping_of_each_variabe_;
  std::vector<jtf::algorithm::equation<double> > eqns_of_each_f;
  size_t nnz_ ;

};

//template<typename T1, typename T2>
//class arap_func <T1,T2,TRI> : public hj::function::function_t<double,int32_t>
//{
//public:
//  arap_func(const zjucad::matrix::matrix_expression<T1> & node,
//            const zjucad::matrix::matrix_expression<T2> & tri,
//            const size_t idx, const double w,
//            const hj::sparse::csc<double,size_t> * node_mapping = 0)
//    :node_num_(node().size(2)), tri_(tri()), w_(w), idx_(idx),
//      node_mapping_(node_mapping){
//    using namespace zjucad::matrix;
//    matrix<double> E = zeros<double>(3,2);
//    E(0,0) = E(1,1) = 1.0;
//    E(2,0) = E(2,1) = -1.0;

//    one_tri_ = tri_(colon(),idx);

//    matrix<double> M = node()(colon(), one_tri_) * E;

//    matrix<double> MTM = trans(M) * M;
//    if(inv(MTM)){
//        std::cerr << "# [error] degenerated triangle." << std::endl;
//      }
//    grad_ = E * MTM * trans(M);
//  }
//  virtual size_t dim_of_x(void) const {
//    if(!node_mapping_)
//      return node_num_*3;
//    else
//      return node_mapping_->size(2);
//  }
//  virtual size_t dim_of_f(void) const {
//    return 3*3;
//  }
//  virtual int val(const double *x, double *f, hj::function::func_ctx *ctx = 0) const {

//    using namespace zjucad::matrix;
//    itr_matrix<double*> f0(3,3,f);

//    if(!node_mapping_){
//        const itr_matrix<const double *> T(3, node_num_, x);
//        f0 = T(colon(),one_tri_) * grad_;
//      }else{
//        matrix<double> tri_node(3,3);
//        for(size_t pi = 0; pi < one_tri_.size(); ++pi){
//            for(size_t di = 0; di < 3; ++di){
//                const size_t variable_idx = 3 * one_tri_[pi] + di;
//                tri_node(di,pi) = //x[(*node_mapping_)()[3 * one_tri_[pi] + di]];
//                    csc_one_line_multy(
//                      &(node_mapping_->idx()[0]) + node_mapping_->ptr()[variable_idx],
//                    &(node_mapping_->idx()[0]) + node_mapping_->ptr()[variable_idx+1],
//                    &(node_mapping_->val()[0]) + node_mapping_->ptr()[variable_idx],
//                    &(node_mapping_->val()[0]) + node_mapping_->ptr()[variable_idx+1],
//                    x);
//              }
//          }
//        f0 = tri_node * grad_;
//      }

//    matrix<double> R = f0;
//    hj::polar3d p;
//    p(R, 2);

//    f0 -= R;
//    f0 *= w_;

//    return 0;
//  }
//  virtual int jac(const double *x, double *val, int32_t *ptr = 0,
//                  int32_t *idx = 0, hj::function::func_ctx *ctx = 0) const {
//    using namespace zjucad::matrix;
//    const itr_matrix<const double *> T(3, node_num_, x);
//    for(size_t c = 0; c < 3; ++c){
//        for(size_t r = 0; r < 3; ++r){
//            ptr[3*c+r+1] = ptr[3*c+r] + one_tri_.size();
//            for(size_t k = 0; k < one_tri_.size(); ++k){
//                idx[ptr[3*c+r]+k] =
//                    (node_mapping_?(*node_mapping_)()[one_tri_[k] * 3 + r]:one_tri_[k]*3+r);
//                val[ptr[3*c+r]+k] = grad_(k,c) * w_;
//                //jtf::util::jtf::math::erase_nan_inf(val[ptr[3*c+r]+k]);
//              }
//          }
//      }
//    return 0;
//  }
//  virtual size_t jac_nnz(void) const {
//    return 3 * dim_of_f();
//  }
//private:
//  const size_t node_num_;
//  const size_t idx_;
//  const zjucad::matrix::matrix<size_t> &tri_;
//  zjucad::matrix::matrix<double> grad_;
//  const double w_;
//  const  hj::sparse::csc<double,size_t> *node_mapping_;
//  zjucad::matrix::matrix<size_t> one_tri_;
//};

///////////////////////////////////////////////////////////////
//template<typename T1, typename T2>
//class asap_func <T1,T2,TET> : public hj::function::function_t<double,int32_t>
//{
//public:
//  asap_func(const zjucad::matrix::matrix_expression<T1> & node,
//            const zjucad::matrix::matrix_expression<T2> & tet,
//            const size_t idx, const double w,
//            const zjucad::matrix::matrix_expression<T2> * node_mapping = 0)
//    :grad_op_(4, 3), tet_(tet()), node_num_(node().size(2)), idx_(idx) , weight_(w),
//      node_mapping_(&((*node_mapping)())){
//    using namespace zjucad::matrix;
//    one_tet_ = tet_(colon(),idx_);

//    zjucad::matrix::matrix<double> tet_node = node()(colon(), one_tet_);
//    if(calc_tet_def_grad_op(&tet_node[0], &grad_op_[0])) {
//        std::cerr << "# degenerated tet" << __LINE__ << tet_node << std::endl;
//      }
//  }
//  virtual ~asap_func(){}
//public:

//  virtual size_t dim_of_x(void) const {
//    if(!node_mapping_)
//      return node_num_*3;
//    else
//      return zjucad::matrix::max((*node_mapping_)())+1;
//  }
//  virtual size_t dim_of_f(void) const {
//    return 9;
//  }
//  virtual int val(const double *x, double *f,
//                  hj::function::func_ctx *ctx = 0) const {
//    using namespace zjucad::matrix;
//    itr_matrix<double *> f0(3, 3, f);

//    if(!node_mapping_){
//        const itr_matrix<const double *> T(3, node_num_, x);
//        f0 = T(colon(), one_tet_)*grad_op_;
//      }else{
//        matrix<double> tet_node(3,4);
//        for(size_t pi = 0; pi < one_tet_.size(); ++pi){
//            for(size_t di = 0; di < 3; ++di){
//                tet_node(di,pi) = x[(*node_mapping_)()[3 * one_tet_[pi] + di]];
//              }
//          }
//        f0 = tet_node * grad_op_;
//      }
//    matrix<double> R = f0, A = f0;
//    hj::polar3d p;
//    p(A, R, 2);
//    const double alpha = 1;
//    R = alpha*R+(1-alpha)*eye<double>(3);
//    const double sacle_v = (fabs(A(0,0)) + fabs(A(1,1)) + fabs(A(2,2)))/A.size(1);
//    f0 -= R * sacle_v;
//    f0 *= weight_;

//    for(size_t i = 0; i < f0.size(); ++i) jtf::math::erase_nan_inf(f0[i]);
//    return 0;
//  }
//  virtual int jac(const double *x, double *val, int32_t *ptr = 0,
//                  int32_t *idx = 0, hj::function::func_ctx *ctx = 0) const {
//    for(int c = 0; c < 3; ++c) {
//        for(int r = 0; r < 3; ++r) {
//            const int fi = c*3+r;
//            ptr[fi+1] = ptr[fi] + one_tet_.size();
//            for(int k = 0; k < one_tet_.size(); ++k) {
//                idx[ptr[fi]+k] = (node_mapping_?(*node_mapping_)()[one_tet_[k]*3+r]:one_tet_[k]*3+r);
//                val[ptr[fi]+k] = grad_op_(k, c) * weight_;
//                jtf::math::erase_nan_inf(val[ptr[fi]+k]);
//              }
//          }
//      }
//    return 0;
//  }
//  virtual size_t jac_nnz(void) const {
//    return 4*dim_of_f();
//  }
//private:
//  const size_t idx_, node_num_;
//  const zjucad::matrix::matrix<size_t> &tet_;
//  zjucad::matrix::matrix<double> grad_op_;
//  const double weight_;
//  const zjucad::matrix::matrix_expression<T2> *node_mapping_;
//  zjucad::matrix::matrix<size_t> one_tet_;
//};

//template<typename T1, typename T2>
//class asap_func <T1,T2,TRI> : public hj::function::function_t<double,int32_t>
//{
//public:
//  asap_func(const zjucad::matrix::matrix_expression<T1> & node,
//            const zjucad::matrix::matrix_expression<T2> & tri,
//            const size_t idx, const double w,
//            const zjucad::matrix::matrix_expression<T2>* node_mapping = 0)
//    :node_num_(node().size(2)), tri_(tri()), w_(w), idx_(idx),
//      node_mapping_(&((*node_mapping)())){
//    using namespace zjucad::matrix;
//    matrix<double> E = zeros<double>(3,2);
//    E(0,0) = E(1,1) = 1.0;
//    E(2,0) = E(2,1) = -1.0;

//    one_tri_ = tri_(colon(),idx);

//    matrix<double> M = node()(colon(), one_tri_) * E;

//    matrix<double> MTM = trans(M) * M;
//    if(inv(MTM)){
//        std::cerr << "# [error] degenerated triangle." << std::endl;
//      }
//    grad_ = E * MTM * trans(M);
//  }
//  virtual size_t dim_of_x(void) const {
//    if(!node_mapping_)
//      return node_num_*3;
//    else
//      return zjucad::matrix::max((*node_mapping_)())+1;
//  }
//  virtual size_t dim_of_f(void) const {
//    return 3*3;
//  }
//  virtual int val(const double *x, double *f, hj::function::func_ctx *ctx = 0) const {
//    using namespace zjucad::matrix;

//    itr_matrix<double*> f0(3,3,f);

//    if(!node_mapping_){
//        const itr_matrix<const double *> T(3, node_num_, x);
//        matrix<double> tri_node = T(colon(), one_tri_);
//        f0 = T(colon(),one_tri_) * grad_;
//      }else{
//        matrix<double> tet_node(3,3);
//        for(size_t pi = 0; pi < one_tri_.size(); ++pi){
//            for(size_t di = 0; di < 3; ++di){
//                tet_node(di,pi) = x[(*node_mapping_)()[3 * one_tri_[pi] + di]];
//              }
//          }
//        f0 = tet_node * grad_;
//      }
//    matrix<double> R = f0;
//    hj::polar3d p;
//    p(R, 2);
//    const double scale_v = zjucad::matrix::dot(
//          f0(zjucad::matrix::colon()),R(zjucad::matrix::colon()))/
//        dot(R(zjucad::matrix::colon()), R(zjucad::matrix::colon()));
//    f0 -= R * scale_v;
//    f0 *= w_;

//    return 0;
//  }
//  virtual int jac(const double *x, double *val, int32_t *ptr = 0,
//                  int32_t *idx = 0, hj::function::func_ctx *ctx = 0) const {
//    using namespace zjucad::matrix;
//    const itr_matrix<const double *> T(3, node_num_, x);
//    for(size_t c = 0; c < 3; ++c){
//        for(size_t r = 0; r < 3; ++r){
//            ptr[3*c+r+1] = ptr[3*c+r] + one_tri_.size();
//            for(size_t k = 0; k < one_tri_.size(); ++k){
//                idx[ptr[3*c+r]+k] = (node_mapping_?(*node_mapping_)()[one_tri_[k] * 3 + r]:one_tri_[k]*3+r);
//                val[ptr[3*c+r]+k] = grad_(k,c) * w_;
//                //jtf::util::jtf::math::erase_nan_inf(val[ptr[3*c+r]+k]);
//              }
//          }
//      }
//    return 0;
//  }
//  virtual size_t jac_nnz(void) const {
//    return 3 * dim_of_f();
//  }
//  virtual size_t hes_nnz(void) const {
//    return 0;
//  }
//private:
//  const size_t node_num_;
//  const size_t idx_;
//  const zjucad::matrix::matrix<size_t> &tri_;
//  zjucad::matrix::matrix<double> grad_;
//  const double w_;
//  const zjucad::matrix::matrix_expression<T2> *node_mapping_;
//  zjucad::matrix::matrix<size_t> one_tri_;
//};

/////////////////////////////////////////////////////////////////
///
/// @brief build_arap_func build arap function on tet/tri mesh.
/// @param tet
/// @param node
/// @param bi
/// @param sim
/// @return
///
template <typename T1, typename T2>
hj_func * build_arap_func(
    const zjucad::matrix::matrix_expression<T1> & mesh,
    const zjucad::matrix::matrix_expression<T2> & node,
    const SIMPLEX sim,
    const node_mapping * node_mapping_ = 0)
{
  using namespace std;
  using namespace zjucad::matrix;

  if(node_mapping_){
      if((*node_mapping_).ZT.size(1) == 0 || (*node_mapping_).ZT.size(2) == 0) {
          cerr << "# [error] node_mapping should not be empty." << endl;
          return 0;
        }
    }

  matrix<double> weight = ones<double>(mesh().size(2),1);

  matrix<double> one_simplex_node = zeros<double>(3,mesh().size(1));
  for(size_t ti = 0; ti < mesh().size(2); ++ti){
      one_simplex_node = node()(colon(), mesh()(colon(),ti));
      weight[ti] = fabs((sim==TET?jtf::mesh::cal_tet_vol(one_simplex_node)
                                :jtf::mesh::cal_face_area(one_simplex_node)));
    }

  const double total_w = std::accumulate(weight.begin(), weight.end(), 0.0);

  std::shared_ptr<vector<std::shared_ptr<hj::function::function_t<double,int32_t> > > >
      funcs(new vector<std::shared_ptr<hj::function::function_t<double,int32_t> > >(mesh().size(2)));

  for(size_t ti = 0; ti < mesh().size(2); ++ti) {
      if(sim == TET)
        (*funcs)[ti].reset(
            new arap_func<T2,T1,TET>(
              node(), mesh(), ti, sqrt(weight[ti]/total_w), node_mapping_));
      else{
          std::cerr << "# [error] unsupported simplex type." << std::endl;
          return 0;
        }
    }

  return hj::function::new_catenated_function<double, int32_t>(funcs);
}

#endif // ARAP_H
