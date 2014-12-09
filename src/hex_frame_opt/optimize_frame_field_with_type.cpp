#include <zjucad/ptree/ptree.h>
#include <jtflib/optimizer/optimizer.h>
#include <hjlib/math_func/math_func.h>
#include <hjlib/math_func/operation.h>
#include <hjlib/math_func/func_aux.h>
#include "optimize_frame_field_with_type.h"
#include "../common/zyz.h"
#include "../common/transition_type.h"

using namespace std;
using namespace zjucad::matrix;

template <typename val_type, typename int_type>
class orthog_func : public hj::math_func::math_func_t<val_type, int_type>
{
public:
  orthog_func(const size_t tet_num, const size_t ti, const double w)
    :tet_num_(tet_num), ti_(ti), w_(w){}
  virtual ~orthog_func(){}
  virtual size_t nx() const{
    return 9 * tet_num_;
  }
  virtual size_t nf() const{
    return 6;
  }
  virtual int eval(size_t k, const val_type *x, const hj::math_func::coo2val_t<val_type, int_type> &cv,
                   hj::math_func::func_ctx *ctx = 0) const
  {
    if(k == 0){
        itr_matrix<const val_type*> x0(9,tet_num_, x);
        itr_matrix<const val_type*> rot(3,3 ,&x0(0,ti_));
        matrix<val_type> diff = trans(rot) *rot - eye<val_type>(3);
        for(int_type i = 0; i < 3; ++i){
            for(int_type j = i; j < 3; ++j){
                int_type c[1] = {i*3+j-i*(i+1)/2};
                cv[c] += w_*diff(i,j);
              }
          }
      }
    if(k == 1){
        itr_matrix<const val_type*> x0(9,tet_num_, x);
        itr_matrix<const val_type*> rot(3,3 ,&x0(0,ti_));
        for(int_type i = 0; i < 3; ++i){
            for(int_type j = i; j < 3; ++j){
                if(i == j){
                    for(int_type k = 0; k < 3; ++k){
                        int_type c[2] = {i*3+j-i*(i+1)/2, 9*ti_+3*i+k};
                        cv[c] += w_*2*rot[3*i+k];
                      }
                  }else if(i == 0 && j == 1){
                    for(int_type k = 0; k < 5; ++k){
                        int_type c[2] = {i*3+j-i*(i+1)/2, 9*ti_+k};
                        cv[c] += w_*rot[(k+3)%6];
                      }
                  }else if(i == 0 && j == 2){
                    for(int_type k = 0; k < 5; ++k){
                        int_type c[2] = {i*3+j-i*(i+1)/2, 9*ti_+k};
                        cv[c] += w_*rot[(k+6)%12];
                      }
                  }else if(i == 1 && j == 2){
                    for(int_type k = 3; k < 9; ++k){
                        int_type c[2] = {i*3+j-i*(i+1)/2, 9*ti_+k};
                        if(k < 5) cv[c] += w_*rot[k+3];
                        else cv[c] += w_*rot[k-3];
                      }
                  }
              }
          }
      }
    return 0;
  }

  virtual int patt(size_t k, hj::math_func::coo_set<int_type> &cs,
                   const hj::math_func::coo_l2g &l2g,
                   hj::math_func::func_ctx *ctx = 0) const
  {
    if(k == 1){
        for(int_type i = 0; i < 3; ++i){
            for(int_type j = i; j < 3; ++j){
                if(i == j){
                    for(int_type k = 0; k < 3; ++k){
                        int_type c[2] = {i*3+j-i*(i+1)/2, 9*ti_+3*i+k};
                        l2g.add(cs,c);
                      }
                  }else if(i == 0 && j == 1){
                    for(int_type k = 0; k < 5; ++k){
                        int_type c[2] = {i*3+j-i*(i+1)/2, 9*ti_+k};
                        l2g.add(cs,c);
                      }
                  }else if(i == 0 && j == 2){
                    for(int_type k = 0; k < 5; ++k){
                        int_type c[2] = {i*3+j-i*(i+1)/2, 9*ti_+k};
                        l2g.add(cs,c);
                      }
                  }else if(i == 1 && j == 2){
                    for(int_type k = 3; k < 9; ++k){
                        int_type c[2] = {i*3+j-i*(i+1)/2, 9*ti_+k};
                        l2g.add(cs,c);
                      }
                  }
              }
          }
      }
    return 0;
  }
  virtual size_t nnz(size_t k) const
  {
    if(k == 0) return -1;
    if(k == 1) return 27;
  }
private:
  const size_t tet_num_;
  const int_type ti_;
  const double w_;
};

template <typename val_type, typename int_type>
class inner_smooth_func_with_type : public hj::math_func::math_func_t<val_type, int_type>
{
public:
  inner_smooth_func_with_type(const size_t tet_num, const size_t ti, const size_t tj,
                              const zjucad::matrix::matrix<size_t> &rot, const double w)
    :tet_num_(tet_num), ti_(ti), tj_(tj), rot_(rot), w_(w){}
  virtual ~inner_smooth_func_with_type(){}
  virtual size_t nx() const{
    return 9 * tet_num_;
  }
  virtual size_t nf() const{
    return 9;
  }
  virtual int eval(size_t k, const val_type *x, const hj::math_func::coo2val_t<val_type, int_type> &cv,
                   hj::math_func::func_ctx *ctx = 0) const
  {
    if(k == 0){
        itr_matrix<const val_type*> x0(9,tet_num_, x);
        itr_matrix<const val_type*> rot1(3,3 ,&x0(0,ti_));
        itr_matrix<const val_type*> rot2(3,3, &x0(0,tj_));
        matrix<val_type> diff = rot1 * rot_ - rot2;
        for(int_type i = 0; i < 9; ++i){
            int_type c[1] = {i};
            cv[c] += w_*diff[i];
          }
      }
    if(k == 1){
        itr_matrix<const val_type*> x0(9,tet_num_, x);
        for(int_type i = 0; i < 9; ++i){
            int_type col = i/3;
            int_type row = i%3;
            for(int_type j = 0; j < 3; ++j){
                int_type c[2] = {i,9*ti_+j*3+row};
                cv[c] += w_*rot_(j,col);
              }
            int_type c[2] = {i,9*tj_+i};
            cv[c] += -w_;
          }
      }
    return 0;
  }

  virtual int patt(size_t k, hj::math_func::coo_set<int_type> &cs,
                   const hj::math_func::coo_l2g &l2g,
                   hj::math_func::func_ctx *ctx = 0) const
  {
    if(k == 1){
        for(int_type i = 0; i < 9; ++i){
            int_type col = i/3;
            int_type row = i%3;
            for(int_type j = 0; j < 3; ++j){
                int_type c[2] = {i,9*ti_+j*3+row};
                l2g.add(cs,c);
              }
            int_type c[2] = {i,9*tj_+i};
            l2g.add(cs,c);
          }
      }
    return 0;
  }
  virtual size_t nnz(size_t k) const
  {
    if(k == 0) return -1;
    if(k == 1) return 36;
  }
private:
  const size_t tet_num_;
  const int_type ti_, tj_;
  const zjucad::matrix::matrix<double> rot_;
  const double w_;
};

shared_ptr<const hj::math_func::math_func_t<double, int32_t> >
build_inner_smooth_func_with_type(
    const zjucad::matrix::matrix<size_t> &tet_mesh,
    const jtf::mesh::face2tet_adjacent &fa,
    const zjucad::matrix::matrix<double> &vol,
    const boost::unordered_map<pair<size_t,size_t>,size_t> &inner_type,
    const double w)
{
  typedef hj::math_func::math_func_t<double,int32_t> math_func_type;
  typedef std::shared_ptr<const math_func_type> math_func_ptr;

  const double total_v = std::accumulate(zjucad::matrix::fabs(vol).begin(),
                                         zjucad::matrix::fabs(vol).end(), 0.0);

  shared_ptr<vector<math_func_ptr> > all_func(new vector<math_func_ptr>);
  matrix<double> rot = eye<double>(3);
  for(size_t fi = 0; fi < fa.face2tet_.size(); ++fi){
      const pair<size_t,size_t> & tet_pair = fa.face2tet_[fi];
      if(fa.is_outside_face(tet_pair)) continue;
      auto it = inner_type.find(tet_pair);
      if(it == inner_type.end()) rot = eye<double>(3);
      else rot = type_transition2(it->second);
      double weight = sqrt(w*(vol[tet_pair.first]+vol[tet_pair.second])/(4.0*total_v));
      all_func->push_back(math_func_ptr(
                            new inner_smooth_func_with_type<double,int32_t>(
                              tet_mesh.size(2), tet_pair.first, tet_pair.second,
                              rot, weight)));
    }
  math_func_ptr fun_cat(new hj::math_func::fcat<double,int32_t, vector<math_func_ptr> >(all_func));
  return fun_cat;
}

shared_ptr<const hj::math_func::math_func_t<double, int32_t> >
build_orthogonal_func(
    const zjucad::matrix::matrix<size_t> &tet_mesh,
    const zjucad::matrix::matrix<double> &vol,
    const double w)
{
  typedef hj::math_func::math_func_t<double,int32_t> math_func_type;
  typedef std::shared_ptr<const math_func_type> math_func_ptr;

  const double total_v = std::accumulate(zjucad::matrix::fabs(vol).begin(),
                                         zjucad::matrix::fabs(vol).end(), 0.0);

  shared_ptr<vector<math_func_ptr> > all_func(new vector<math_func_ptr>);

  for(size_t ti = 0; ti < tet_mesh.size(2); ++ti){
      const double weight = sqrt(w*vol[ti]/total_v);
      all_func->push_back(math_func_ptr(
                            new orthog_func<double,int32_t>(
                              tet_mesh.size(2), ti, weight)));
    }
  math_func_ptr fun_cat(new hj::math_func::fcat<double,int32_t, vector<math_func_ptr> >(all_func));
  return fun_cat;
}

void optimize_frame_field_with_type(const jtf::tet_mesh &tm,
                                    const boost::unordered_map<std::pair<size_t,size_t>,size_t> &inner_type,
                                    zjucad::matrix::matrix<double> &zyz)
{
  typedef hj::math_func::math_func_t<double,int32_t> math_func_type;
  typedef std::shared_ptr<math_func_type> math_func_ptr;

  if(zyz.size(2) != tm.tetmesh_.mesh_.size(2))
    throw std::logic_error("zyz size is not compatible with tetmesh.");

  matrix<double> rot_m(9,zyz.size(2));
  for(size_t ti = 0; ti < zyz.size(2); ++ti){
      zyz_angle_2_rotation_matrix1(&zyz(0,ti), &rot_m(0,ti));
    }

  rot_m += rand<double>(9,zyz.size(2));
  shared_ptr<vector<math_func_ptr> > func(new vector<math_func_ptr>);
  {
    func->push_back(math_func_ptr(
                      new hj::math_func::sumsqr<double,int32_t>(
                        build_inner_smooth_func_with_type(tm.tetmesh_.mesh_, *tm.fa_, tm.vol_, inner_type,1))));
    func->push_back(math_func_ptr(
                      new hj::math_func::sumsqr<double,int32_t>(
                        build_orthogonal_func(tm.tetmesh_.mesh_, tm.vol_, 10))));
  }

  math_func_ptr func_cat(new hj::math_func::fcat<double, int32_t, vector<math_func_ptr> >(func));
  math_func_ptr obj(new hj::math_func::sum<double,int32_t>(func_cat));

  boost::property_tree::ptree pt;
  pt.put("package.value","jtf");
  pt.put("alg.value", "SQP");
  pt.put("iter.value",10);
  jtf::optimize(*obj, rot_m, pt, nullptr, nullptr, nullptr);

  for(size_t ti = 0; ti < rot_m.size(2); ++ti){
      rotation_matrix_2_zyz_angle(&rot_m(0,ti), &zyz(0,ti),0);
    }
}

