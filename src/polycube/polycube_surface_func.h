#ifndef POLYCUBE_SURFACE_FUNC_H
#define POLYCUBE_SURFACE_FUNC_H

#include <numeric>
#include <set>
#include <boost/unordered_set.hpp>


#include <hjlib/function/func_aux.h>
#include <hjlib/function/function.h>
#include <zjucad/matrix/matrix.h>
#include <zjucad/matrix/io.h>
#include <hjlib/math/blas_lapack.h>
#include <zjucad/matrix/lapack.h>
#include <zjucad/matrix/itr_matrix.h>
#include <hjlib/math/polar.h>
#include <jtflib/mesh/mesh.h>
#include "../tetmesh/tetmesh.h"
#include "../mesh_func/tri-area.h"
#include "../mesh_func/tri-normal.h"
#include "../mesh_func/so3.h"
#include "../numeric/util.h"
#include "smooth_L1.h"

/////////// high level func ////////////////////////

hj::function::function_t<double, int32_t> *
build_fix_zero_node_func(
    const zjucad::matrix::matrix<double> &node,
    const zjucad::matrix::matrix<double> &zero_pos,
    const double weight);

hj::function::function_t<double, int32_t> *
build_surface_patch_align_func(
    const zjucad::matrix::matrix<double> &node,
    const std::vector<zjucad::matrix::matrix<size_t> > &surface_patch,
    const double patch_w);

hj::function::function_t<double, int32_t> *
build_feature_line_func(
    const zjucad::matrix::matrix<double> &node,
    const zjucad::matrix::matrix<size_t> &outside_face,
    const jtf::mesh::edge2cell_adjacent &ea,
    const boost::unordered_set<std::pair<size_t,size_t> > &feature_line,
    const double feature_line_w);

hj::function::function_t<double, int32_t> *
build_group_equation_func(
    const zjucad::matrix::matrix<double> &node,
    const std::vector<boost::unordered_set<size_t> > &fnode_group,
    const double group_weight);

jtf::function::functionN1_t<double,int32_t> *
build_polycube_rot_func2(const zjucad::matrix::matrix<double> &node,
                         const zjucad::matrix::matrix<size_t> &faces,
                         const zjucad::matrix::matrix<double> &areas);

hj::function::function_t<double, int32_t> *
build_patch_normal_fix_func(
    const zjucad::matrix::matrix<double> &node,
    const std::vector<zjucad::matrix::matrix<double> > & patch_normal,
    const std::vector<zjucad::matrix::matrix<size_t> > &surface_patches,
    const double normal_fix_w,
    const std::vector<std::vector<size_t> > * select_ptr = 0);


jtf::function::functionN1_t<double,int32_t> *
build_smooth_L1_area_normal(const matrix<double> &node, const matrix<size_t> &tri,
                            const matrix<double> &areas, double normal_align_w);


hj::function::function_t<double, int32_t> *
build_polycube_function(
    const matrix<double> &orig_node,
    const matrix<double> &deform_node,
    const matrix<double> &zero_pos,
    const matrix<size_t> &tet,
    const matrix<size_t> & faces,
    hj::function::function_t<double,int32_t> *&diagnose,
    const boost::unordered_set<pair<size_t,size_t> > *feature_line = nullptr,
    const vector<matrix<size_t> > *surface_patch = nullptr,
    const double adj_normal_w = 1, const double fl_w = 0, const double patch_w = 0);

///////////////////////////////////////////////////
/////////// low level func ////////////////////////

class polycube_edge_dir_func : public hj::function::function_t<double, int32_t>
{
public:
  polycube_edge_dir_func(const zjucad::matrix::matrix<double> &node,
                         const zjucad::matrix::matrix<size_t> &tri_a,
                         const zjucad::matrix::matrix<size_t> &tri_b,
                         const double &weight)
    :tri_a_(tri_a), tri_b_(tri_b), node_num_(node.size(2)), weight_(weight) {
  }
  virtual size_t dim_of_x(void) const {
    return node_num_*3;
  }
  virtual size_t dim_of_f(void) const {
    return 3;
  }
  virtual int val(const double *x, double *f, hj::function::func_ctx *ctx = 0) const {
    const zjucad::matrix::itr_matrix<const double *> T(3, node_num_, x);
    zjucad::matrix::matrix<double> tri_a = T(zjucad::matrix::colon(), tri_a_);
    zjucad::matrix::matrix<double> tri_b = T(zjucad::matrix::colon(), tri_b_);
    zjucad::matrix::matrix<double> fa(3,1) , fb(3,1);
    calc_tri_normal_(&fa[0], &tri_a[0]);
    calc_tri_normal_(&fb[0], &tri_b[0]);
    if(fabs(zjucad::matrix::norm(fa)-1) > 0.5){
        std::cerr << "# [info] degenerated tri." << std::endl;
        std::cerr << fa << std::endl;
        fa *= 0.0;
      }
    if(fabs(zjucad::matrix::norm(fb)-1) > 0.5){
        std::cerr << "# [info] degenerated tri." << std::endl;
        std::cerr << fb << std::endl;
        fb *= 0.0;
      }
    for(int i = 0; i < 3; ++i) {
        f[i] = (fa[i] - fb[i]) * weight_;
        jtf::math::erase_nan_inf(f[i]);
      }
    return 0;
  }
  virtual int jac(const double *x, double *val, int32_t *ptr = 0,
                  int32_t *idx = 0, hj::function::func_ctx *ctx = 0) const {
    const zjucad::matrix::itr_matrix<const double *> T(3, node_num_, x);
    zjucad::matrix::matrix<double> tri_a =
        T(zjucad::matrix::colon(), tri_a_), jac_a(3, 9);
    zjucad::matrix::matrix<double> tri_b =
        T(zjucad::matrix::colon(), tri_b_), jac_b(3, 9);

    calc_tri_normal_jac_(&jac_a[0], &tri_a[0]);
    calc_tri_normal_jac_(&jac_b[0], &tri_b[0]);

    for(int fi = 0; fi < 3; ++fi) {
        ptr[fi+1] = ptr[fi] + 18;
        for(int xi = 0; xi < 9; ++xi) {
            idx[ptr[fi]+xi] = tri_a_[xi/3]*3+xi%3;
            val[ptr[fi]+xi] = jac_a(fi, xi)*weight_;
            jtf::math::erase_nan_inf(val[ptr[fi]+xi]);
          }

        for(int xi = 0; xi < 9; ++xi) {
            idx[ptr[fi]+xi+9] = tri_b_[xi/3]*3+xi%3;
            val[ptr[fi]+xi+9] = -1* jac_b(fi, xi)*weight_;
            jtf::math::erase_nan_inf(val[ptr[fi]+xi+9]);
          }
      }
    return 0;
  }
  virtual size_t jac_nnz(void) const {
    return 18 * dim_of_f();
  }
private:
  const size_t node_num_;
  const zjucad::matrix::matrix<size_t> tri_a_, tri_b_;
  const double weight_;
};

//class face_normal_align_func : public hj::function::function_t<double, int32_t>
//{
//public:
//  face_normal_align_func(const zjucad::matrix::matrix<double> &node,
//                         const zjucad::matrix::matrix<size_t> &tri,
//                         const size_t & type, const double &weight)
//    :tri_(tri), type_(type), node_num_(node.size(2)), weight_(weight) {
//    dir_.resize(3,1);
//  }
//  virtual size_t dim_of_x(void) const {
//    return node_num_*3;
//  }
//  virtual size_t dim_of_f(void) const {
//    return 3;
//  }
//  virtual int val(const double *x, double *f, hj::function::func_ctx *ctx = 0) const {
//    zjucad::matrix::itr_matrix<const double*> node_(3, node_num_, x);
//    calc_tri_normal_(&dir_[0], const_cast<double*>(&node_(zjucad::matrix::colon(),
//                                      tri_(zjucad::matrix::colon()))[0]));
//    const zjucad::matrix::matrix<double>Identity = zjucad::matrix::eye<double>(3);
//    if(zjucad::matrix::dot(dir_,
//                           Identity(zjucad::matrix::colon(),type_)) >
//       zjucad::matrix::dot(-1*dir_, Identity(zjucad::matrix::colon(),type_))){
//      dir_ = Identity(zjucad::matrix::colon(), type_);
//    }
//    else{
//      dir_ = -1*Identity(zjucad::matrix::colon(), type_);
//    }

//    const zjucad::matrix::itr_matrix<const double *> T(3, node_num_, x);
//    zjucad::matrix::matrix<double> tri = T(zjucad::matrix::colon(), tri_);
//    calc_tri_normal_(f, &tri[0]);
//    for(int i = 0; i < 3; ++i){
//      f[i] = (f[i] - dir_[i])* weight_;
//			jtf::math::erase_nan_inf(f[i]);
//    }
//    return 0;
//  }
//  virtual int jac(const double *x, double *val, int32_t *ptr = 0,
//                  int32_t *idx = 0, hj::function::func_ctx *ctx = 0) const {
//    const zjucad::matrix::itr_matrix<const double *> T(3, node_num_, x);
//    zjucad::matrix::matrix<double> tri =
//        T(zjucad::matrix::colon(), tri_), jac(3, 9);
//    calc_tri_normal_jac_(&jac[0], &tri[0]);
//    for(int fi = 0; fi < 3; ++fi) {
//      ptr[fi+1] = ptr[fi] + 9;
//      for(int xi = 0; xi < 9; ++xi) {
//        idx[ptr[fi]+xi] = tri_[xi/3]*3+xi%3;
//        val[ptr[fi]+xi] = jac(fi, xi)*weight_;
//				jtf::math::erase_nan_inf(val[ptr[fi]+xi]);
//      }
//    }
//    return 0;
//  }
//  virtual size_t jac_nnz(void) const {
//    return 9*dim_of_f();
//  }
//private:
//  const size_t type_;
//  const size_t node_num_;
//  const zjucad::matrix::matrix<size_t> tri_;
//  mutable zjucad::matrix::matrix<double> dir_;
//  const double weight_;
//};

class surface_anti_flipped_func : public hj::function::function_t<double, int32_t>
{
public:
  surface_anti_flipped_func(const zjucad::matrix::matrix<double> &node,
                            const zjucad::matrix::matrix<size_t> &tri,
                            const double &weight)
    :tri_(tri), node_num_(node.size(2)), weight_(weight), delta_(0.1){
    zjucad::matrix::matrix<double> tri_nodes =
        node(zjucad::matrix::colon(), tri_);
    calc_tri_area_(&ori_area_, &tri_nodes[0]);
  }
  virtual size_t dim_of_x(void) const {
    return node_num_*3;
  }
  virtual size_t dim_of_f(void) const {
    return 1;
  }
  virtual int val(const double *x, double *f, hj::function::func_ctx *ctx = 0) const {
    const zjucad::matrix::itr_matrix<const double *> node(3, node_num_, x);
    zjucad::matrix::matrix<double> tri_nodes =
        node(zjucad::matrix::colon(), tri_);
    double area;

    calc_tri_area_(&area, &tri_nodes[0]);

    *f = area / ori_area_;
    *f += delta_;
    *f = weight_ * (1.0/fabs(*f));
    jtf::math::erase_nan_inf(*f);
    return 0;
  }

  virtual int jac(const double *x, double *val, int32_t *ptr = 0,
                  int32_t *idx = 0, hj::function::func_ctx *ctx = 0) const {
    const zjucad::matrix::itr_matrix<const double *> node(3, node_num_, x);
    zjucad::matrix::matrix<double> tri_nodes =
        node(zjucad::matrix::colon(), tri_);
    double area;

    calc_tri_area_jac_(val+ptr[0], &tri_nodes[0]);
    calc_tri_area_(&area, &tri_nodes[0]);

    area /= ori_area_;
    area += delta_;
    for(int i = 0; i < 9; ++i) {
        val[ptr[0]+i] /= ori_area_;
        val[ptr[0]+i] *= -1.0/(area*area);
        val[ptr[0]+i] *= area < 0?-1*weight_:1 * weight_;
        jtf::math::erase_nan_inf(val[ptr[0]+i]);
      }
    ptr[1] = ptr[0]+9;
    for(int ni = 0; ni < 3; ++ni) {
        for(int di = 0; di < 3; ++di) {
            idx[ptr[0]+ni*3+di] = tri_[ni]*3+di;
          }
      }
    return 0;
  }
  virtual size_t jac_nnz(void) const {
    return 9;
  }
private:
  const size_t node_num_;
  const zjucad::matrix::matrix<size_t> tri_;
  double ori_area_;
  double weight_;
  const double delta_;
};

class polycube_face_normal_func : public hj::function::function_t<double, int32_t>
{
public:
  polycube_face_normal_func(const zjucad::matrix::matrix<double> &node,
                            const zjucad::matrix::matrix<size_t> &tri,
                            const double &weight)
    :tri_(tri), node_num_(node.size(2)), weight_(weight){
  }
  virtual size_t dim_of_x(void) const {
    return node_num_*3;
  }
  virtual size_t dim_of_f(void) const {
    return 3;
  }
  virtual int val(const double *x, double *f,
                  hj::function::func_ctx *ctx = 0) const {
    const zjucad::matrix::itr_matrix<const double *> T(3, node_num_, x);
    zjucad::matrix::matrix<double> tri = T(zjucad::matrix::colon(), tri_);
    calc_tri_normal_(f, &tri[0]);
    zjucad::matrix::itr_matrix<double*> f0(3,1,f);
    if(fabs(zjucad::matrix::norm(f0)-1) > 0.5){
        std::cerr << "# [info] degenerated tri." << std::endl;
        std::cerr << f0 << std::endl;
        f0 *= 0.0;
      }
    for(int i = 0; i < 3; ++i){
        f[i] *= weight_;
        jtf::math::erase_nan_inf(f[i]);
      }
    return 0;
  }
  virtual int jac(const double *x, double *val, int32_t *ptr = 0,
                  int32_t *idx = 0, hj::function::func_ctx *ctx = 0) const {
    const zjucad::matrix::itr_matrix<const double *> T(3, node_num_, x);
    zjucad::matrix::matrix<double> tri =
        T(zjucad::matrix::colon(), tri_), jac(3, 9);
    calc_tri_normal_jac_(&jac[0], &tri[0]);
    for(int fi = 0; fi < 3; ++fi) {
        ptr[fi+1] = ptr[fi] + 9;
        for(int xi = 0; xi < 9; ++xi) {
            idx[ptr[fi]+xi] = tri_[xi/3]*3+xi%3;
            val[ptr[fi]+xi] = jac(fi, xi)*weight_;
            jtf::math::erase_nan_inf(val[ptr[fi]+xi]);
          }
      }
    return 0;
  }
  virtual size_t jac_nnz(void) const {
    return 9*dim_of_f();
  }
private:
  const size_t node_num_;
  const zjucad::matrix::matrix<size_t> tri_;
  const double weight_;
};

class anti_flip_edge_func : public hj::function::function_t<double, int32_t>
{
public:
  anti_flip_edge_func(const zjucad::matrix::matrix<double> & node,
                      const zjucad::matrix::matrix<size_t> & tri_a,
                      const zjucad::matrix::matrix<size_t> & tri_b,
                      const double weight)
    : node_num_(node.size(2)), w_(weight), delta_(-1){
    std::vector<size_t> sort_points;
    for(size_t i = 0; i < tri_a.size(); ++i)
      for(size_t j = 0; j < tri_b.size(); ++j){
          if(tri_a[i] == tri_b[j]) {
              sort_points.push_back(tri_a[i]);
              break;
            }
        }
    assert(sort_points.size() == 2);
    const size_t shared_points =
        std::accumulate(sort_points.begin(), sort_points.end(), static_cast<size_t>(0));
    sort_points.push_back(
          std::accumulate(tri_a.begin(), tri_a.end(), static_cast<size_t>(0)) - shared_points);
    sort_points.push_back(
          std::accumulate(tri_b.begin(), tri_b.end(), static_cast<size_t>(0)) - shared_points);

    tp_.resize(4,1);
    std::copy(sort_points.begin(), sort_points.end(), tp_.begin());
  }

  virtual size_t dim_of_x(void) const {
    return node_num_*3;
  }
  virtual size_t dim_of_f(void) const {
    return 1;
  }
  virtual int val(const double *x, double *f, hj::function::func_ctx *ctx = 0) const {
    const zjucad::matrix::itr_matrix<const double *> T(3, node_num_, x);
    zjucad::matrix::matrix<double> tp_node = T(zjucad::matrix::colon(), tp_);

    calc_tri_normal_dot_(f, &tp_node[0]);

    *f -= delta_;
    *f = 1.0/(*f);
    jtf::math::erase_nan_inf(*f);
    return 0;
  }
  virtual int jac(const double *x, double *val, int32_t *ptr = 0,
                  int32_t *idx = 0, hj::function::func_ctx *ctx = 0) const {
    const zjucad::matrix::itr_matrix<const double *> T(3, node_num_, x);
    zjucad::matrix::matrix<double> tp_node = T(zjucad::matrix::colon(), tp_);
    double f = 0;
    calc_tri_normal_dot_(&f, &tp_node[0]);
    double a = -1.0 / ((f-delta_) * (f-delta_));
    zjucad::matrix::matrix<double> jac_tp(1,12);

    calc_tri_normal_dot_jac_(&jac_tp[0], &tp_node[0]);

    ptr[1] = ptr[0] + 12;
    for(size_t pi = 0; pi < 4; ++pi)
      for(size_t di = 0; di < 3; ++di){
          idx[ptr[0] + pi * 3 + di] = tp_[pi] * 3 + di;
          val[ptr[0] + pi * 3 + di] = jac_tp[pi * 3 + di] * w_ * a;
          jtf::math::erase_nan_inf(val[ptr[0] + pi * 3 + di]);
        }
    return 0;
  }
  virtual size_t jac_nnz(void) const {
    return 12 * dim_of_f();
  }
private:
  const size_t node_num_;
  zjucad::matrix::matrix<size_t> tp_;
  const double w_;
  const double delta_;
};

class normal_dot_func : public hj::function::function_t<double,int32_t>
{
public:
  normal_dot_func(const zjucad::matrix::matrix<double> &node,
                  const zjucad::matrix::matrix<size_t> &tri_a,
                  const zjucad::matrix::matrix<size_t> &tri_b,
                  const double &weight,
                  const double epsilon = 0)
    :node_num_(node.size(2)), weight_(weight),  epsilon_(epsilon) {
    std::vector<size_t> sort_points;
    for(size_t i = 0; i < tri_a.size(); ++i)
      for(size_t j = 0; j < tri_b.size(); ++j){
          if(tri_a[i] == tri_b[j]) {
              sort_points.push_back(tri_a[i]);
              break;
            }
        }
    assert(sort_points.size() == 2);
    const size_t shared_points =
        std::accumulate(sort_points.begin(), sort_points.end(), static_cast<size_t>(0));
    sort_points.push_back(
          std::accumulate(tri_a.begin(), tri_a.end(), static_cast<size_t>(0)) - shared_points);
    sort_points.push_back(
          std::accumulate(tri_b.begin(), tri_b.end(), static_cast<size_t>(0)) - shared_points);

    tp_.resize(4,1);
    std::copy(sort_points.begin(), sort_points.end(), tp_.begin());

#ifdef debug
    {
      zjucad::matrix::matrix<double> tri_a_node = node(zjucad::matrix::colon(), tri_a);
      zjucad::matrix::matrix<double> tri_b_node = node(zjucad::matrix::colon(), tri_b);
      double area[] = {0,0};
      calc_tri_area_(&area[0], &tri_a_node[0]);
      calc_tri_area_(&area[1], &tri_b_node[0]);
      if(fabs(area[0]) <  1e-4 || fabs(area[1]) < 1e-4)
        std::cerr << "# [error] degenerated triangle." << std::endl;
    }
#endif
  }
  virtual size_t dim_of_x(void) const {
    return node_num_*3;
  }
  virtual size_t dim_of_f(void) const {
    return 1;
  }
  virtual int val(const double *x, double *f, hj::function::func_ctx *ctx = 0) const {
    const zjucad::matrix::itr_matrix<const double *> T(3, node_num_, x);
    zjucad::matrix::matrix<double> tp_node = T(zjucad::matrix::colon(), tp_);

    calc_tri_normal_dot_(f, &tp_node[0]);

    *f -= epsilon_;
    *f *= weight_;
    jtf::math::erase_nan_inf(*f);
    return 0;
  }
  virtual int jac(const double *x, double *val, int32_t *ptr = 0,
                  int32_t *idx = 0, hj::function::func_ctx *ctx = 0) const {
    const zjucad::matrix::itr_matrix<const double *> T(3, node_num_, x);
    zjucad::matrix::matrix<double> tp_node = T(zjucad::matrix::colon(), tp_);

    zjucad::matrix::matrix<double> jac_tp(1,12);

    calc_tri_normal_dot_jac_(&jac_tp[0], &tp_node[0]);

    ptr[1] = ptr[0] + 12;
    for(size_t pi = 0; pi < 4; ++pi)
      for(size_t di = 0; di < 3; ++di){
          idx[ptr[0] + pi * 3 + di] = tp_[pi] * 3 + di;
          val[ptr[0] + pi * 3 + di] = jac_tp[pi * 3 + di] * weight_;
          jtf::math::erase_nan_inf(val[ptr[0] + pi * 3 + di]);
        }
    return 0;
  }
  virtual size_t jac_nnz(void) const {
    return 12 * dim_of_f();
  }
private:
  const double epsilon_;
  const size_t node_num_;
  //const zjucad::matrix::matrix<size_t> tri_a_, tri_b_;
  zjucad::matrix::matrix<size_t> tp_;
  const double weight_;
};

class fix_node : public hj::function::function_t<double,int32_t>
{
public:
  fix_node(const zjucad::matrix::matrix<double> & node,
           const size_t & p,
           const double w)
    : position_(node(colon(),p)), node_num_(node.size(2)), p_(p), w_(w){}

  virtual size_t dim_of_x(void) const {
    return node_num_*3;
  }
  virtual size_t dim_of_f(void) const {
    return 3;
  }
  virtual int val(const double *x, double *f, hj::function::func_ctx *ctx = 0) const {
    const zjucad::matrix::itr_matrix<const double*> T(3, node_num_, x);
    zjucad::matrix::itr_matrix<double*> f0(3,1,f);
    f0  = T(colon(),p_) - position_;
    f0 *= w_;
    return 0;
  }
  virtual int jac(const double *x, double *val, int32_t *ptr = 0,
                  int32_t *idx = 0, hj::function::func_ctx *ctx = 0) const {
    for(int fi = 0; fi < 3; ++fi) {
        ptr[fi+1] = ptr[fi]+1;
        idx[ptr[fi]] = 3 * p_ + fi;
        val[ptr[fi]] = w_;
      }
    return 0;
  }
  virtual size_t jac_nnz(void) const {
    return 1*dim_of_f();
  }
private:
  const zjucad::matrix::matrix<double> position_;
  const size_t node_num_;
  const size_t p_;
  const double w_;
};

class fix_zero_node_func : public hj::function::function_t<double, int32_t>
{
public:
  fix_zero_node_func(const zjucad::matrix::matrix<double> &node,
                     const zjucad::matrix::matrix<double> &zero_pos,
                     const double weight)
    :node_num_(node.size(2)) , weight_(weight){
    zero_node_pos_ = zero_pos;
  }
  virtual size_t dim_of_x(void) const {
    return node_num_*3;
  }
  virtual size_t dim_of_f(void) const {
    return 3;
  }
  virtual int val(const double *x, double *f, hj::function::func_ctx *ctx = 0) const {
    std::copy(x, x+3, f);
    zjucad::matrix::itr_matrix<double*> f0(3,1,f);
    f0 -= zero_node_pos_;
    f0 *= weight_;
    return 0;
  }
  virtual int jac(const double *x, double *val, int32_t *ptr = 0,
                  int32_t *idx = 0, hj::function::func_ctx *ctx = 0) const {
    for(int fi = 0; fi < 3; ++fi) {
        ptr[fi+1] = ptr[fi]+1;
        idx[ptr[fi]] = fi;
        val[ptr[fi]] = weight_;
      }
    return 0;
  }
  virtual size_t jac_nnz(void) const {
    return 9*dim_of_f();
  }
protected:
  size_t node_num_;
  zjucad::matrix::matrix<double> zero_node_pos_;
  const double weight_;
};


class fix_node_func : public hj::function::function_t<double, int32_t>
{
public:
  fix_node_func(const zjucad::matrix::matrix<double> &node,
                const size_t p,
                const double weight)
    :fix_node_(node(zjucad::matrix::colon(),p)), p_(p), node_num_(node.size(2)) ,
      weight_(weight){
    //fix_node_ = zeros<double>(3,1);
  }
  virtual size_t dim_of_x(void) const {
    return node_num_*3;
  }
  virtual size_t dim_of_f(void) const {
    return 3;
  }
  virtual int val(const double *x, double *f, hj::function::func_ctx *ctx = 0) const {
    zjucad::matrix::itr_matrix<const double*> T(3, node_num_, x);
    zjucad::matrix::itr_matrix<double*> f0(3,1,f);
    f0 = (T(zjucad::matrix::colon(),p_) - fix_node_) * weight_;
    return 0;
  }
  virtual int jac(const double *x, double *val, int32_t *ptr = 0,
                  int32_t *idx = 0, hj::function::func_ctx *ctx = 0) const {
    for(int fi = 0; fi < 3; ++fi) {
        ptr[fi+1] = ptr[fi]+1;
        idx[ptr[fi]] = 3 * p_ + fi;
        val[ptr[fi]] = weight_;
      }
    return 0;
  }
  virtual size_t jac_nnz(void) const {
    return 1*dim_of_f();
  }
protected:
  zjucad::matrix::matrix<double> fix_node_;
  const size_t p_;
  size_t node_num_;
  const double weight_;
};

class equation_func : public hj::function::function_t<double, int32_t>
{
public:
  equation_func(const zjucad::matrix::matrix<double> & node,
                const size_t & fi,
                const size_t & fj,
                const double weight)
    : node_num_(node.size(2)), fi_(fi), fj_(fj), weight_(weight){}
  virtual size_t dim_of_x()const{
    return 3 * node_num_;
  }
  virtual size_t dim_of_f()const{
    return 1;
  }
  virtual int val(const double *x, double *f,
                  hj::function::func_ctx *ctx = 0) const {
    const zjucad::matrix::itr_matrix<const double *> T(3, node_num_, x);
    *f = (T[fi_] - T[fj_]) * weight_;
    return 0;
  }
  virtual int jac(const double *x, double *val, int32_t *ptr = 0,
                  int32_t *idx = 0, hj::function::func_ctx *ctx = 0) const {
    ptr[1] = ptr[0] + 2;

    idx[ptr[0] + 0] = fi_;
    val[ptr[0] + 0] = weight_;

    idx[ptr[0] + 1] = fj_;
    val[ptr[0] + 1] = -1*weight_;

    return 0;
  }
  virtual size_t jac_nnz() const{
    return 2 * dim_of_f();
  }

private:
  const size_t node_num_;
  const size_t fi_;
  const size_t fj_;
  const double weight_;
};

class global_rot : public hj::function::function_t<double, int32_t>
{
public:
  global_rot(const zjucad::matrix::matrix<double> &normal,
             const zjucad::matrix::matrix<double> &weight)
    :normal_(normal), weight_(weight) {
  }
  global_rot(const zjucad::matrix::matrix<double> &normal)
    :normal_(normal),weight_(zjucad::matrix::ones<double>(normal_.size(2),1)){
  }
  virtual size_t dim_of_x(void) const {
    return 9;
  }
  virtual size_t dim_of_f(void) const {
    return normal_.size();
  }
  virtual int val(const double *x, double *f, hj::function::func_ctx *ctx = 0) const {
    zjucad::matrix::itr_matrix<const double *> R(3, 3, x);
    zjucad::matrix::itr_matrix<double *> err(3, normal_.size(2), f);
    err = R*normal_;
    for(size_t ni = 0; ni < normal_.size(2); ++ni){
        err(colon(),ni) *= weight_[ni];
      }
    return 0;
  }
  virtual int jac(const double *x, double *val, int32_t *ptr = 0,
                  int32_t *idx = 0, hj::function::func_ctx *ctx = 0) const {
    zjucad::matrix::itr_matrix<const double *> R(3, 3, x);
    for(size_t ni = 0; ni < normal_.size(2); ++ni) {
        for(int fi = 0; fi < 3; ++fi, ptr+=1) {
            ptr[1] = ptr[0]+3;
            for(int i = 0; i < 3; ++i) {
                idx[ptr[0]+i] = fi+i*3;
                val[ptr[0]+i] = normal_(i, ni) * weight_[ni];
              }
          }
      }
    return 0;
  }
  virtual size_t jac_nnz(void) const {
    return dim_of_f()*3;
  }
protected:
  const zjucad::matrix::matrix<double> normal_;
  const zjucad::matrix::matrix<double> weight_;
};

class so3 : public hj::function::function_t<double, int32_t>
{
public:
  so3()
    :weight_(1e0){}
  so3(double w)
    :weight_(w){}
  virtual size_t dim_of_x(void) const {
    return 9;
  }
  virtual size_t dim_of_f(void) const {
    return 9;
  }
  virtual int val(const double *x, double *f, hj::function::func_ctx *ctx = 0) const {
    zjucad::matrix::itr_matrix<const double *> R(3, 3, x);
    zjucad::matrix::itr_matrix<double *> err(3, 3, f);
    calc_so3_(&err[0], const_cast<double*>(&R[0]));
    err *= weight_;
    return 0;
  }
  virtual int jac(const double *x, double *val, int32_t *ptr = 0,
                  int32_t *idx = 0, hj::function::func_ctx *ctx = 0) const {
    calc_so3_jac_(val, const_cast<double*>(x));
    zjucad::matrix::itr_matrix<double *> v(9, 9, val);
    v = zjucad::matrix::temp(zjucad::matrix::trans(v))*weight_;
    for(size_t fi = 0; fi < 9; ++fi) {
        ptr[fi+1] = ptr[fi]+9;
        for(size_t xi = 0; xi < 9; ++xi) {
            idx[ptr[fi]+xi] = xi;
          }
      }
    return 0;
  }
  virtual size_t jac_nnz(void) const {
    return 9*9;
  }
  double weight_;
};

class normal_fix_func:public hj::function::function_t<double,int32_t>
{
public:
  normal_fix_func(const zjucad::matrix::matrix<double> & node,
                  const zjucad::matrix::matrix<double> & normal,
                  const zjucad::matrix::matrix<size_t> & tri,
                  const double weight)
    :node_num_(node.size(2)), normal_(normal), tri_(tri), w_(weight){}
  virtual size_t dim_of_x(void) const {
    return node_num_*3;
  }
  virtual size_t dim_of_f(void) const {
    return 3;
  }
  virtual int val(const double *x, double *f, hj::function::func_ctx *ctx = 0) const {
    const zjucad::matrix::itr_matrix<const double *> T(3, node_num_, x);
    zjucad::matrix::matrix<double> tri_node = T(zjucad::matrix::colon(), tri_);
    zjucad::matrix::matrix<double> fa(3,1);
    calc_tri_normal_(&fa[0], &tri_node[0]);

    for(int i = 0; i < 3; ++i) {
        f[i] = (fa[i] - normal_[i]) * w_;
        jtf::math::erase_nan_inf(f[i]);
      }
    return 0;
  }
  virtual int jac(const double *x, double *val, int32_t *ptr = 0,
                  int32_t *idx = 0, hj::function::func_ctx *ctx = 0) const {
    const zjucad::matrix::itr_matrix<const double *> T(3, node_num_, x);
    zjucad::matrix::matrix<double> tri_node = T(zjucad::matrix::colon(), tri_),
        jac_(3, 9);

    calc_tri_normal_jac_(&jac_[0], &tri_node[0]);

    for(int fi = 0; fi < 3; ++fi) {
        ptr[fi+1] = ptr[fi] + 9;
        for(int xi = 0; xi < 9; ++xi) {
            idx[ptr[fi]+xi] = tri_[xi/3]*3+xi%3;
            val[ptr[fi]+xi] = jac_(fi, xi)*w_;
            jtf::math::erase_nan_inf(val[ptr[fi]+xi]);
          }

      }
    return 0;
  }
  virtual size_t jac_nnz(void) const {
    return 9 * dim_of_f();
  }
private:
  const size_t node_num_;
  const zjucad::matrix::matrix<double> normal_;
  const zjucad::matrix::matrix<size_t> tri_;
  const double w_;
};

class surface_arap_func: public hj::function::function_t<double,int32_t>
{
public:
  surface_arap_func(const zjucad::matrix::matrix<double> & node,
                    const zjucad::matrix::matrix<size_t> & one_face,
                    const double weight)
    :node_num_(node.size(2)), tri_(one_face), w_(weight){
    matrix<double> E = zeros<double>(3,2);
    E(0,0) = E(1,1) = 1.0;
    E(2,0) = E(2,1) = -1.0;

    matrix<double> M = node(colon(), one_face) * E;

    matrix<double> MTM = trans(M) * M;
    if(inv(MTM)){
        cerr << "# [error] degenerated triangle." << endl;
      }
    grad_ = E * MTM * trans(M);
  }
  virtual size_t dim_of_x(void) const {
    return node_num_*3;
  }
  virtual size_t dim_of_f(void) const {
    return 3*3;
  }
  virtual int val(const double *x, double *f, hj::function::func_ctx *ctx = 0) const {
    const zjucad::matrix::itr_matrix<const double *> T(3, node_num_, x);
    zjucad::matrix::matrix<double> tri_node = T(zjucad::matrix::colon(), tri_);
    itr_matrix<double*> f0(3,3,f);
    f0 = T(colon(),tri_) * grad_;
    matrix<double> R = f0;
    hj::polar3d p;
    p(R, 2);

    f0 -= R;

    f0 *= w_;
    for(size_t i = 0; i < f0.size(); ++i) {
        jtf::math::erase_nan_inf(f0[i]);
      }
    return 0;
  }
  virtual int jac(const double *x, double *val, int32_t *ptr = 0,
                  int32_t *idx = 0, hj::function::func_ctx *ctx = 0) const {
    const zjucad::matrix::itr_matrix<const double *> T(3, node_num_, x);
    for(size_t c = 0; c < 3; ++c){
        for(size_t r = 0; r < 3; ++r){
            ptr[3*c+r+1] = ptr[3*c+r] + 3;
            for(size_t k = 0; k < 3; ++k){
                idx[ptr[3*c+r]+k] = tri_(k) * 3 + r;
                val[ptr[3*c+r]+k] = grad_(k,c) * w_;
                jtf::math::erase_nan_inf(val[ptr[3*c+r]+k]);
              }
          }
      }
    return 0;
  }
  virtual size_t jac_nnz(void) const {
    return 3 * dim_of_f();
  }
private:
  const size_t node_num_;
  const zjucad::matrix::matrix<size_t> tri_;
  zjucad::matrix::matrix<double> grad_;
  const double w_;
};

hj::function::function_t<double, int32_t> *
build_adj_normal_func(
    const matrix<double> & variant,
    const vector<vector<pair<size_t,double> > > & var_map,
    const boost::unordered_map<pair<size_t,size_t>,size_t> & face_rot,
    const jtf::mesh::face2tet_adjacent & fa_cut,
    const matrix<size_t> & cut_tet2tet,
    matrix<double> &weight);

hj::function::function_t<double, int32_t> *
build_adj_normal_func(const zjucad::matrix::matrix<double> & node,
                      const zjucad::matrix::matrix<size_t> & faces,
                      zjucad::matrix::matrix<double> & weight);
//                     const double adj_normal_w,
//     const matrix<double> *point_weight_diffused = 0);

/* hj::function::function_t<double, int32_t> * */
/* build_adj_normal_dynamical_func(const zjucad::matrix::matrix<double> & orig_node, */
/*                                 const zjucad::matrix::matrix<double> & deform_node, */
/*                                 const zjucad::matrix::matrix<size_t> & faces, */
/*                                 const double adj_normal_w, */
/*                                 const matrix<double> *point_weight_diffused = 0); */

hj::function::function_t<double, int32_t> *
build_adj_normal_sub_func(
    const zjucad::matrix::matrix<double> & node,
    const zjucad::matrix::matrix<size_t> & faces,
    double weight);



hj::function::function_t<double, int32_t> *
build_surface_arap_func(const zjucad::matrix::matrix<double> &node,
                        const zjucad::matrix::matrix<size_t> &faces);

hj::function::function_t<double,int32_t> *
build_area_preserve_func(const size_t node_num,
                         const matrix<size_t> & faces,
                         const double total_area,
                         const double percent,
                         const double weight = 1.0);

#endif

