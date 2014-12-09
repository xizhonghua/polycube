#ifndef HEX_DEFORM_FUNCTION_H
#define HEX_DEFORM_FUNCTION_H
#include "../hexmesh/util.h"
#include "../common/util.h"
#include <numeric>

#include <zjucad/matrix/matrix.h>
#include <zjucad/matrix/itr_matrix.h>
#include <hjlib/function/function.h>
#include <hjlib/function/func_aux.h>
#include <hjlib/math/blas_lapack.h>
#include <hjlib/math/polar.h>
#include <zjucad/matrix/lapack.h>

//! @brief: rigid matching function should use this one, function is assigned on
//  each point
class rigid_matching_function_p: public hj::function::function_t<double, int32_t>
{
public:
  rigid_matching_function_p(const matrixst & one_hex,
                            const matrixd & node,
                            const matrixd & weight,
                            const size_t & hex_node_idx,
                            const double w = 1)
    :one_hex_(one_hex),node_num_(node.size(2)),weight_(weight),
      point_idx_(hex_node_idx), w_(w)
  {
    using namespace zjucad::matrix;
    assert(weight.size() == 8 && one_hex.size() == 8);

    const matrixd hex_node = node(colon(), one_hex);
    const double d = calc_bounding_sphere_size(hex_node);
    //average_edge /= edges.size(2);
    const double average_edge = 0.57735 * d; // sqrt(3)/3 * d
    rigid_node_ = zjucad::matrix::zeros<double>(3,8);
    rigid_node_(1,1) = -1 * average_edge; // node 1: (0,-a,0)
    rigid_node_(0,2) = -1 * average_edge; // node 2: (-a,0,0)
    rigid_node_(0,3) = rigid_node_(1,3) = -1 * average_edge; // node 3: (-a,-a,0)
    rigid_node_(2,4) = -1 * average_edge; // node 4: (0,0,-a)
    rigid_node_(1,5) = rigid_node_(2,5) = -1 * average_edge; // node 5: (0,-a,-a)
    rigid_node_(0,6) = rigid_node_(2,6) = -1 * average_edge; // node 6: (-a,0,-a)
    rigid_node_(zjucad::matrix::colon(),7) =
        zjucad::matrix::ones<double>(3,1) * (-1*average_edge);//node 7:(-a,-a,-a)
    //////////////////////////////

    total_weight_ = std::accumulate(weight_.begin(), weight_.end(),0.0);
    assert(fabs(total_weight_)> 1e-6);
    rigid_center_ = zeros<double>(3,1);
    for(size_t pi = 0; pi < 8; ++pi){
      rigid_center_ += weight_[pi] * rigid_node_(colon(), pi) ;
    }
    rigid_center_ /= total_weight_;

    Aqq_ = zjucad::matrix::zeros<double>(3,3);
    for(size_t pi = 0; pi < rigid_node_.size(2); ++pi){
      Aqq_ += weight_[pi] *
              (rigid_node_(colon(),pi) - rigid_center_)
              * trans((rigid_node_(colon(),pi)) - rigid_center_);
    }

    if(inv(Aqq_)){
      std::cerr << "# [error] strange degenerate rigid hex." << std::endl;
    }
  }

  virtual ~rigid_matching_function_p(){}
  virtual size_t dim_of_x(void) const {
    return 3 * node_num_;
  }
  virtual size_t dim_of_f(void) const {
    return 3 ;
  }
  virtual int val(const double *x, double *f, hj::function::func_ctx *ctx = 0) const
  {
    using namespace zjucad::matrix;
    itr_matrix<const double*> x_mat(3, node_num_, x);
    itr_matrix<double *> f_mat(3,1, f);
    matrixd center_new = zeros<double>(3,1);

    for(size_t pi = 0; pi < one_hex_.size(); ++pi){
      // TODO: this format is not support in current matrix version
      // center_new += x_mat(colon(), one_hex_[pi]) * weight_[pi];
      matrixd test = x_mat(colon(), one_hex_[pi]);
      test *= weight_[pi];
      center_new += test;
    }
    center_new /= total_weight_;

    matrixd Apq = zeros<double>(3,3);

    for(size_t pi = 0; pi < one_hex_.size(); ++pi){
      matrixd test = x_mat(colon(), one_hex_[pi]);
      Apq += weight_[pi]
             // TODO: this format is not support in current matrix version
             //* (x_mat(colon(), one_hex_[pi]) - center_new )
             * (test  - center_new)
             * trans((rigid_node_(colon(),pi) - rigid_center_)) ;
    }
    matrixd R = Apq ;//* Aqq_;
    hj::polar3d p;
    if(p(R, 2) < 0)
      std::cerr << "# [error] polar fail." << Apq << Aqq_ << std::endl;

    assert(R.size(1) == R.size(2));
    assert(R.size(2) == 3);
    //for(size_t pi = 0; pi < 8; ++pi)
    f_mat = R * (rigid_node_(colon(),point_idx_) - rigid_center_ )
            - (x_mat(colon(), one_hex_[point_idx_]) - center_new);
    f_mat *= w_;
    return 0;
  }

  // TODO: is it correct??
  virtual int jac(const double *x, double *val, int32_t *ptr = 0, int32_t *idx = 0,
                  hj::function::func_ctx *ctx = 0) const
  {
    using namespace zjucad::matrix;
    itr_matrix<const double*> x_mat(3, node_num_, x);
    size_t pi = 0;
    for(size_t fr = 0; fr < 3; ++fr, ++pi){
      ptr[pi+1] = ptr[pi]  + 8;
      for(size_t fx = 0; fx < 8; ++fx){
        idx[ptr[pi] + fx] =  3 * one_hex_[fx] + fr;
        val[ptr[pi] + fx] = weight_[fx]/total_weight_;
        if(fx == point_idx_)
          val[ptr[pi] + fx] -= 1;
        val[ptr[pi]+fx] *= w_;
      }
    }

    return 0;
  }
  virtual size_t jac_nnz(void) const {
    return 3 * 8;
  }

private:
  const matrixst one_hex_;
  const matrixd weight_;
  const size_t point_idx_;
  double total_weight_;
  matrixd rigid_node_;
  matrixd rigid_center_;
  matrixd Aqq_;
  const double w_;
  const size_t node_num_;

};
//// for each hex, design a rest cube whose edge length equals the average length
//// of input hex
//class rigid_matching_function: public hj::function::function_t<double, int32_t>
//{
//public:
//  rigid_matching_function(const matrixst & one_hex,
//                          const matrixd & node,
//                          const matrixd & weight)
//    :one_hex_(one_hex),node_num_(node.size(2)),weight_(weight)
//  {
//    using namespace zjucad::matrix;
//    assert(weight.size() == 8 && one_hex.size() == 8);
//    //// setting rigid hex  ////
//    //    matrixst edges;
//    //    jtf::hexmesh::get_edges_for_one_hex(one_hex, edges);
//    //double average_edge = 0.0;
//    //    for(size_t ei = 0; ei < edges.size(2); ++ei){
//    //      average_edge +=
//    //          norm(node(colon(), edges(0, ei)) - node(colon(), edges(1, ei)));
//    //    }
//    const double d = calc_bounding_sphere_size(node(colon(), one_hex));
//    //average_edge /= edges.size(2);
//    const double average_edge = 0.57735 * d; // sqrt(3)/3 * d
//    rigid_node_ = zjucad::matrix::zeros<double>(3,8);
//    rigid_node_(1,1) = -1 * average_edge; // node 1: (0,-a,0)
//    rigid_node_(0,2) = -1 * average_edge; // node 2: (-a,0,0)
//    rigid_node_(0,3) = rigid_node_(1,3) = -1 * average_edge; // node 3: (-a,-a,0)
//    rigid_node_(2,4) = -1 * average_edge; // node 4: (0,0,-a)
//    rigid_node_(1,5) = rigid_node_(2,5) = -1 * average_edge; // node 5: (0,-a,-a)
//    rigid_node_(0,6) = rigid_node_(2,6) = -1 * average_edge; // node 6: (-a,0,-a)
//    rigid_node_(zjucad::matrix::colon(),7) =
//        zjucad::matrix::ones<double>(3,1) * (-1*average_edge);//node 7:(-a,-a,-a)
//    //////////////////////////////

//    total_weight_ = std::accumulate(weight_.begin(), weight_.end(),0.0);
//    assert(fabs(total_weight_)> 1e-6);
//    rigid_center_ = zeros<double>(3,1);
//    for(size_t pi = 0; pi < 8; ++pi){
//      rigid_center_ += weight_[pi] * rigid_node_(colon(), pi) ;
//    }
//    rigid_center_ /= total_weight_;

//    Aqq_ = zjucad::matrix ::zeros<double>(3,3);
//    for(size_t pi = 0; pi < rigid_node_.size(2); ++pi){
//      Aqq_ += weight_[pi] *
//              (rigid_node_(colon(),pi) - rigid_center_)
//              * trans((rigid_node_(colon(),pi)) - rigid_center_);
//    }
//    if(inv(Aqq_) == false){
//      std::cerr << "# [error] strange degenerate rigid hex." << std::endl;
//    }
//  }

//  virtual ~rigid_matching_function(){}
//  virtual size_t dim_of_x(void) const {
//    return 3 * node_num_;
//  }
//  virtual size_t dim_of_f(void) const {
//    return 3 * 8;
//  }
//  virtual int val(const double *x, double *f, hj::function::func_ctx *ctx = 0) const
//  {
//    using namespace zjucad::matrix;
//    itr_matrix<const double*> x_mat(3, node_num_, x);
//    itr_matrix<double *> f_mat(3,8, f);
//    matrixd center_new = zeros<double>(3,1);

//    for(size_t pi = 0; pi < one_hex_.size(); ++pi){
//      // TODO: this format is not support in current matrix version
//      // center_new += x_mat(colon(), one_hex_[pi]) * weight_[pi];
//      matrixd test = x_mat(colon(), one_hex_[pi]);
//      test *= weight_[pi];
//      center_new += test;
//    }
//    center_new /= total_weight_;

//    matrixd Apq = zeros<double>(3,3);

//    for(size_t pi = 0; pi < one_hex_.size(); ++pi){
//      matrixd test = x_mat(colon(), one_hex_[pi]);
//      Apq += weight_[pi]
//             // TODO: this format is not support in current matrix version
//             //* (x_mat(colon(), one_hex_[pi]) - center_new )
//             * (test  - center_new)
//             * trans((rigid_node_(colon(),pi) - rigid_center_)) ;
//    }
//    matrixd R = Apq ;//* Aqq_;
//    hj::polar3d p;
//    if(p(R, 2) < 0)
//      std::cerr << "# [error] polar fail." << Apq << Aqq_ << std::endl;

//    assert(R.size(1) == R.size(2) == 3);
//    for(size_t pi = 0; pi < 8; ++pi)
//      f_mat(colon(), pi) =
//          R * (rigid_node_(colon(),pi) - rigid_center_ )
//          - (x_mat(colon(), one_hex_[pi]) - center_new);
//    return 0;
//  }

//  // TODO: is it correct??
//  virtual int jac(const double *x, double *val, int32_t *ptr = 0, int32_t *idx = 0,
//                  hj::function::func_ctx *ctx = 0) const
//  {
//    using namespace zjucad::matrix;
//    itr_matrix<const double*> x_mat(3, node_num_, x);
//    size_t pi = 0;
//    for(size_t fr = 0; fr < 3; ++fr){
//      for(size_t fc = 0; fc < 8; ++fc, ++pi){
//        ptr[pi+1] = ptr[pi]  + 8;
//        for(size_t fx = 0; fx < 8; ++fx){
//          idx[ptr[pi] + fx] =  3 * one_hex_[fc] + fr;
//          val[ptr[pi] + fx] = weight_[fc]/total_weight_;
//          if(fx == fc)
//            val[ptr[pi] + fx] -= 1;
//        }
//      }
//    }

//    return 0;
//  }
//  virtual size_t jac_nnz(void) const {
//    return 3 * 8 * 8;
//  }

//private:
//  const matrixst one_hex_;
//  const matrixd weight_;
//  double total_weight_;
//  matrixd rigid_node_;
//  matrixd rigid_center_;
//  matrixd Aqq_;
//  const size_t node_num_;
//};

class laplacian_smooth_function
    : public hj::function::function_t<double, int32_t>
{
public:
  laplacian_smooth_function(const size_t & point_idx,
                            const size_t & node_num,
                            const std::vector<size_t> & adj_points)
    :point_idx_(point_idx), node_num_(node_num), adj_points_(adj_points){}
  virtual ~laplacian_smooth_function(){}
  virtual size_t dim_of_x(void) const {
    return 3 * node_num_;
  }
  virtual size_t dim_of_f(void) const {
    return 3;
  }
  virtual int val(const double *x, double *f, hj::function::func_ctx *ctx = 0) const {
    using namespace zjucad::matrix;
    itr_matrix<const double*> x_mat(3, node_num_, x);
    matrixd center_new = zeros<double>(3,1);

    for(size_t pi = 0; pi < adj_points_.size(); ++pi){
      center_new += x_mat(colon(), adj_points_[pi]);
    }
    center_new /= adj_points_.size();

    itr_matrix<double*> f_mat(3,1, f);
    f_mat = x_mat(colon(), point_idx_) - center_new;
    return 0;
  }

  virtual int jac(const double *x, double *val, int32_t *ptr = 0, int32_t *idx = 0,
                  hj::function::func_ctx *ctx = 0) const {
    for(size_t fi = 0; fi < 3; ++fi){
      ptr[fi+1] = ptr[fi] + adj_points_.size() + 1;
      size_t pi = 0;
      idx[ptr[fi] + pi] = 3 * point_idx_ + fi;
      val[ptr[fi] + pi] = 1;

      for(size_t api = 0, pi = 1; api < adj_points_.size(); ++api, ++pi){
        idx[ptr[fi] + pi] = 3 * adj_points_[api] + fi;
        val[ptr[fi] + pi] = -1.0 / adj_points_.size();
      }
    }
    return 0;
  }
  virtual size_t jac_nnz(void) const {
    return 3 * (adj_points_.size() + 1);
  }

private:
  const size_t point_idx_;
  const size_t node_num_;
  const std::vector<size_t> adj_points_;

};

class tangent_moveing_function
    : public hj::function::function_t<double, int32_t>
{
public:
  tangent_moveing_function(
      const size_t & point_idx,
      const size_t & node_num,
      const matrixd & point_node,
      const matrixd & normal)
    : point_idx_(point_idx), node_num_(node_num),
      point_node_(point_node),normal_(normal){}
  virtual ~tangent_moveing_function(){}
  virtual size_t dim_of_x(void) const {
    return 3 * node_num_;
  }
  virtual size_t dim_of_f(void) const {
    return 1;
  }
  virtual int val(const double *x, double *f, hj::function::func_ctx *ctx = 0) const {
    using namespace zjucad::matrix;
    itr_matrix<const double*> x_mat(3, node_num_, x);

    *f = dot(normal_, point_node_ - x_mat(colon(), point_idx_));

    return 0;
  }

  virtual int jac(const double *x, double *val, int32_t *ptr = 0, int32_t *idx = 0,
                  hj::function::func_ctx *ctx = 0) const {
    ptr[1] = ptr[0] + 3;
    for(size_t t = 0; t < 3; ++t){
      idx[ptr[0] + t] = 3 * point_idx_ + t;
      val[ptr[0] + t] = -1 * normal_[t];
    }
    return 0;
  }
  virtual size_t jac_nnz(void) const {
    return 3 ;
  }
private:
  const size_t point_idx_;
  const size_t node_num_;
  const matrixd point_node_;
  const matrixd normal_;
};

#endif // HEX_DEFORM_FUNCTION_H
