#ifndef SURFACE_SMOOTH_H
#define SURFACE_SMOOTH_H

#include <hjlib/math_func/math_func.h>
#include <zjucad/matrix/matrix.h>
#include <zjucad/matrix/itr_matrix.h>
#include <hjlib/math_func/operation.h>
#include "../numeric/csc_filler.h"
#include <jtflib/mesh/mesh.h>
#include <jtflib/mesh/util.h>

template <typename val_type, typename int_type>
class center_smooth_math_func: public hj::math_func::math_func_t<val_type,int_type>
{
public:
  center_smooth_math_func(const size_t num, const size_t idx,
                          const std::vector<size_t> & one_ring_points,
                          const double weight)
    :point_num_(num), point_idx_(idx), one_ring_points_(one_ring_points), weight_(weight){}

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
        matrix<val_type> f = T(colon(), point_idx_);
        for(size_t i  = 0; i < one_ring_points_.size(); ++i){
            f -= T(colon(), one_ring_points_[i]);
          }
        f *= weight_;
        for(int_type di = 0; di < 3; ++di){
            int_type c[] = {di};
            cv[c] += f[di];
          }
      }
    if(k == 1){
        for(int_type di = 0; di < 3; ++di){
            for(int_type pi = 0; pi < one_ring_points_.size(); ++pi){
                int_type c[] = {di, int_type(3*one_ring_points_[pi]+di)};
                cv[c] += -1*weight_;
              }
            int_type c[] = {di, int_type(3*point_idx_+di)};
            cv[c] += weight_;
          }
      }
    return 0;
  }
  virtual int patt(size_t k, hj::math_func::coo_set<int_type> &cs,
                   const hj::math_func::coo_l2g &l2g,
                   hj::math_func::func_ctx *ctx = 0) const{
    if(k == 1){
        for(int_type di = 0; di < 3; ++di){
            for(int_type pi = 0; pi < one_ring_points_.size(); ++pi){
                int_type c[] = {di, int_type(3*one_ring_points_[pi]+di)};
                l2g.add(cs,c);
              }
            int_type c[] = {di, int_type(3*point_idx_+di)};
            l2g.add(cs,c);
          }
      }
    return 0;
  }
  virtual size_t nnz(size_t k) const{
    if(k == 0) return -1;
    if(k == 1) return 3*(one_ring_points_.size() + 1);
  }

private:
  const size_t point_num_;
  const size_t point_idx_;
  const std::vector<size_t> one_ring_points_;
  const double weight_;
};


template <typename val_type>
std::shared_ptr<const hj::math_func::math_func_t<val_type, int32_t> >
build_surface_smooth_math_func(const zjucad::matrix::matrix<size_t> &mesh,
                               const zjucad::matrix::matrix<val_type> &node,
                               const val_type weight)
{
  using namespace zjucad::matrix;
  typedef hj::math_func::math_func_t<val_type, int32_t> math_func_type;
  typedef std::shared_ptr<math_func_type> math_func_ptr;
  std::shared_ptr<std::vector<math_func_ptr> > all_func(new std::vector<math_func_ptr>);

  matrix<val_type> face_area(mesh.size(2),1);
  for(size_t fi = 0; fi < mesh.size(2); ++fi)
    face_area[fi] = jtf::mesh::cal_face_area(mesh(zjucad::matrix::colon(),fi), node);

  const val_type total_area = std::accumulate(face_area.begin(), face_area.end(), 0.0);
  matrix<val_type> point_area_ = zjucad::matrix::zeros<val_type>(node.size(2),1);

  for(size_t fi = 0; fi < mesh.size(2); ++fi){
      for(size_t pi = 0; pi < mesh.size(1); ++pi){
          point_area_[mesh(pi,fi)] += fabs(face_area[fi])/3.0;
        }
    }

  std::shared_ptr<jtf::mesh::one_ring_point_at_point> orpap(
        jtf::mesh::one_ring_point_at_point::create(mesh));
  assert(orpap.get());
  for(size_t pi = 0; pi < node.size(2); ++pi){
      auto it = orpap->p2p_.find(pi);
      assert(it != orpap->p2p_.end());
      all_func->push_back(
            math_func_ptr(
              new center_smooth_math_func<double,int32_t>(node.size(2), pi, it->second, sqrt(weight*point_area_[pi]/total_area))));
    }
  std::cerr << "# [info] smooth point function number " << all_func->size() << std::endl;
  math_func_ptr fun_cat(new hj::math_func::fcat<double,int32_t, std::vector<math_func_ptr> >(all_func));
  return fun_cat;
}

#endif // SURFACE_SMOOTH_H
