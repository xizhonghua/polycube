/* detail control param
/usr/bin/time ../../bin/polycube prog=polycube3 linear_solver/type=direct linear_solver/name=cholmod iter_w=20 output="a.tet" iter=100 tet=../../dat/kitty-4.8k.mesh-split-surface-tet.mesh normal_align_w=5e-2 L1_sqrt_eps=5e-1 epsg=1e-3 adj_normal_w=2e1
 */
#include <boost/property_tree/ptree.hpp>

#include <string>
#include <fstream>
#include <numeric>
#include <hjlib/function/func_aux.h>
#include <hjlib/math/polar.h>
#include <hjlib/sparse/sparse.h>


#include <zjucad/optimizer/optimizer.h>
#include <zjucad/ptree/ptree.h>
#include <hjlib/math/blas_lapack.h>
#include <zjucad/matrix/lapack.h>

#include "../tetmesh/tetmesh.h"
#include "../common/vtk.h"
#include "../common/util.h"
#include "../common/IO.h"
#include <jtflib/mesh/mesh.h>

#include "io.h"
#include "polycube_surface_func.h"
#include "smooth_L1.h"

#include "polycube_surface_func.h"
#include "../mesh_func/tri-area.h"
#include "../mesh_func/tri-area-normal.h"
#include "tet_func.h"
#include "util.h"

#include <jtflib/function/function.h>
#include <jtflib/optimizer/optimizer.h>
#include <jtflib/function/func_aux.h>
#include <hjlib/math_func/math_func.h>
#include <hjlib/math_func/operation.h>

using namespace std;
using boost::property_tree::ptree;
using namespace hj::function;
using namespace zjucad::matrix;

extern double L1_sqrt_eps;
extern double arap_w_control;

template <typename val_type, typename int_type>
class normal_align_func : public hj::math_func::math_func_t<val_type, int_type>
{
public:
  normal_align_func(const zjucad::matrix::matrix<double> &node,
                    const zjucad::matrix::matrix<size_t> &tri,
                    const zjucad::matrix::matrix<double> & target_area_normal,
                    const double weight)
    :node_num_(node.size(2)), tri_(tri), target_area_normal_(target_area_normal),
      weight_(weight){
    for(size_t di = 0; di < 3; ++di)
      jtf::math::erase_nan_inf(target_area_normal_[di]);
  }
  virtual ~normal_align_func(){}
  virtual size_t nx() const{
    return 3*node_num_;
  }
  virtual size_t nf() const{
    return 1;
  }
  virtual int eval(size_t k, const val_type *x, const hj::math_func::coo2val_t<val_type, int_type> &cv,
                   hj::math_func::func_ctx *ctx = 0) const
  {
    if(k == 0){
        zjucad::matrix::itr_matrix<const val_type *> T(3, node_num_, x);
        zjucad::matrix::matrix<val_type> tri = T(zjucad::matrix::colon(), tri_);
        val_type diff = 0;

        {
          const val_type area = jtf::mesh::cal_face_area(tri);
          matrix<val_type> area_normal = area*target_area_normal_;
          tri_area_normal_diff_(&diff, &tri[0], &area_normal[0]);
        }
        jtf::math::erase_nan_inf(diff);
        int_type c[1] = {0};
        cv[c] +=  diff * weight_;
      }
    if(k == 1){
        zjucad::matrix::itr_matrix<const val_type *> T(3, node_num_, x);
        zjucad::matrix::matrix<val_type> tri = T(zjucad::matrix::colon(), tri_);
        zjucad::matrix::matrix<val_type> fjac(9,1);
        {
          const val_type area = jtf::mesh::cal_face_area(tri);
          matrix<val_type> area_normal = area*target_area_normal_;
          tri_area_normal_diff_jac_(&fjac[0], &tri[0], &area_normal[0]);
        }

        for(size_t i  = 0; i < fjac.size(); ++i) jtf::math::erase_nan_inf(fjac[i]);
        for(int_type i = 0; i < 9; ++i){
            int_type c1[2] = {0,tri_[i/3]*3+i%3};
            cv[c1] += weight_*fjac[i];
          }
      }
    if(k == 2){
        matrix<val_type> H(9, 9);
        zjucad::matrix::itr_matrix<const val_type *> T(3, node_num_, x);
        zjucad::matrix::matrix<val_type> tri = T(zjucad::matrix::colon(), tri_);
        {
          const val_type area = jtf::mesh::cal_face_area(tri);
          matrix<val_type> area_normal = area*target_area_normal_;
          tri_area_normal_diff_hes_(&H[0], &tri[0], &area_normal[0]);
        }

        for(size_t i  = 0; i < H.size(); ++i) jtf::math::erase_nan_inf(H[i]);
        H *= weight_;
        matrix<val_type> e(9), diag_e = zeros<val_type>(9, 9);
        eig(H, e);
        for(size_t ei = 0; ei < 9; ++ei) {
            if(e[ei] > 0)
              diag_e(ei, ei) = e[ei];
            else
              diag_e(ei, ei) = 0;
          }
        H = temp(H*temp(diag_e*trans(H)));
        for(int_type i = 0 ; i < 9; ++i){
            for(int_type j = 0; j < 9; ++j){
                int_type c[3] = {0, tri_[i/3]*3+i%3,tri_[j/3]*3+j%3};
                cv[c] += H(i,j);
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
        for(int_type i = 0; i < 9; ++i){
            int_type c1[2] = {0,tri_[i/3]*3+i%3};
            l2g.add(cs,c1);
          }
      }
    if(k == 2){
        for(int_type i = 0 ; i < 9; ++i){
            for(int_type j = 0; j < 9; ++j){
                int_type c[3] = {0, tri_[i/3]*3+i%3,tri_[j/3]*3+j%3};
                l2g.add(cs,c);
              }
          }
      }
    return 0;
  }
  virtual size_t nnz(size_t k) const
  {
    if(k == 0) return -1;
    if(k == 1) return 9;
    if(k == 2) return 9*9;
  }
private:
  const size_t node_num_;
  zjucad::matrix::matrix<val_type> target_area_normal_;
  const zjucad::matrix::matrix<size_t> tri_;
  const double weight_;
};


shared_ptr<hj::math_func::math_func_t<double,int32_t> >
build_normal_align_func(
    const zjucad::matrix::matrix<double> &node,
    const zjucad::matrix::matrix<size_t> &face,
    const zjucad::matrix::matrix<double> &area,
    const zjucad::matrix::matrix<double> &target_normal,
    const double normal_align)
{
  typedef hj::math_func::math_func_t<double,int32_t> math_func_type;
  typedef std::shared_ptr<math_func_type> math_func_ptr;

  const double total_area = std::accumulate(fabs(area).begin(), fabs(area).end(), 0.0);

  shared_ptr<vector<math_func_ptr> > all_func(new vector<math_func_ptr>);

  for(size_t fi = 0; fi < face.size(2); ++fi) {
      all_func->push_back(
            math_func_ptr(
              new normal_align_func<double,int32_t>(
                node, face(colon(),fi),target_normal(colon(),fi),
                sqrt(normal_align/total_area))));
    }
  math_func_ptr func_cat(new hj::math_func::fcat<double,int32_t, vector<math_func_ptr> >(all_func));
  math_func_ptr obj(new hj::math_func::sum<double,int32_t>(func_cat));

  return obj;
}



class face_area_normal_diff : public jtf::function::functionN1_t<double,int32_t>
{
public:
  face_area_normal_diff(const zjucad::matrix::matrix<double> &node,
                        const zjucad::matrix::matrix<size_t> &tri,
                        const zjucad::matrix::matrix<double> &target_normal,
                        size_t xyz,
                        const double weight)
    :tri_(tri), node_num_(node.size(2)), weight_(weight), xyz_(xyz), target_normal_(target_normal){
  }
  virtual size_t dim(void) const {
    return node_num_*3;
  }
  virtual int val(const double *x, double &v) {
    const zjucad::matrix::itr_matrix<const double *> T(3, node_num_, x);
    zjucad::matrix::matrix<double> tri = T(zjucad::matrix::colon(), tri_);
    double f[3], area;
    area = jtf::mesh::cal_face_area(tri);
    calc_tri_area_normal_(f, &tri[0]);
    jtf::math::erase_nan_inf(f[xyz_]);
    v += (f[xyz_] - fabs(area)*target_normal_[xyz_]) * weight_;
    return 0;
  }

  virtual int gra(const double *x, double *g) {
    double sp_g[9];
    int32_t idx[9];
    size_t nnz = 9;
    gra(x, nnz, sp_g, idx);
    for(size_t i = 0; i < 9; ++i)
      g[idx[i]] += sp_g[i];
    return 0;
  }

  virtual int gra(const double *x, size_t &nnz, double *g, int32_t *idx) {
    if(g == 0) {
        nnz = 9;
        return 0;
      }
    if(nnz != 9) {
        cerr << "strange." << endl;
      }
    const zjucad::matrix::itr_matrix<const double *> T(3, node_num_, x);
    zjucad::matrix::matrix<double> tri =
        T(zjucad::matrix::colon(), tri_), jac(3, 9);
    calc_tri_area_normal_jac_(&jac[0], &tri[0]);
    for(size_t i = 0; i < 9; ++i) {
        idx[i] = tri_[i/3]*3+i%3;
        g[i] = jac(xyz_, i)*weight_;
        jtf::math::erase_nan_inf(g[i]);
      }
    return 0;
  }
  virtual int hes(const double *x, size_t &nnz, size_t &format, double *h, int32_t *ptr, int32_t *idx, double alpha = 1) {
    format = 1;
    if(h == 0 && ptr == 0 && idx == 0) { // query nnz
        nnz = 9*9;
        return 0;
      }
    if(h == 0 && ptr != 0 && idx != 0) { // query pattern
        for(size_t ci = 0; ci < 9; ++ci) {
            const size_t var_ci = tri_[ci/3]*3+ci%3;
            ptr[var_ci+1] = ptr[var_ci]+9;
            for(size_t ri = 0; ri < 9; ++ri) {
                const size_t var_ri = tri_[ri/3]*3+ri%3;
                idx[ptr[var_ci]+ri] = var_ri;
              }
          }
        return 0;
      }
    if(h != 0 && ptr != 0 && idx != 0) { // accumulate
        const zjucad::matrix::itr_matrix<const double *> T(3, node_num_, x);
        zjucad::matrix::matrix<double> tri =
            T(zjucad::matrix::colon(), tri_), hes(27, 9);
        if(is_degenerate(tri)) return 0;
        calc_tri_area_normal_hes_(&hes[0], &tri[0]);
        matrix<double> H(9, 9);
        for(size_t ci = 0; ci < 9; ++ci) {
            const size_t var_ci = tri_[ci/3]*3+ci%3;
            for(size_t ri = 0; ri < 9; ++ri) {
                const size_t var_ri = tri_[ri/3]*3+ri%3;
                H(ri, ci) = hes(xyz_+ri*3, ci)*weight_*alpha;
              }
          }
        const matrix<double> save_H = H;
        matrix<double> e(9), diag_e = zeros<double>(9, 9);
        eig(H, e);
        for(size_t ei = 0; ei < 9; ++ei) {
            if(e[ei] > 0)
              diag_e(ei, ei) = e[ei];
            else
              diag_e(ei, ei) = 0;
          }
        H = temp(H*temp(diag_e*trans(H)));
        for(size_t ci = 0; ci < 9; ++ci) {
            const size_t var_ci = tri_[ci/3]*3+ci%3;
            for(size_t ri = 0; ri < 9; ++ri) {
                const size_t var_ri = tri_[ri/3]*3+ri%3;
                if(jtf::function::add_to_csc(h, ptr, idx, var_ri, var_ci, H(ri, ci))) {
                    return __LINE__;
                  }
              }
          }
        return 0;
      }
    return __LINE__;
  }
  // accumulate
  virtual int hes_block(const double *x, double *h, double alpha = 1) {
    const zjucad::matrix::itr_matrix<const double *> T(3, node_num_, x);
    zjucad::matrix::matrix<double> tri =
        T(zjucad::matrix::colon(), tri_), hes(27, 9);
    if(is_degenerate(tri)) return 0;
    calc_tri_area_normal_hes_(&hes[0], &tri[0]);
    itr_matrix<double *> H(9, 9, h);
    for(size_t ci = 0; ci < 9; ++ci) {
        const size_t var_ci = tri_[ci/3]*3+ci%3;
        for(size_t ri = 0; ri < 9; ++ri) {
            const size_t var_ri = tri_[ri/3]*3+ri%3;
            H(ri, ci) += hes(xyz_+ri*3, ci)*weight_*alpha;
          }
      }
    return 0;
  }
private:
  const size_t node_num_, xyz_;
  const zjucad::matrix::matrix<size_t> tri_;
  const double weight_;
  const zjucad::matrix::matrix<double> target_normal_;
};

jtf::function::functionN1_t<double,int32_t> *
build_smooth_L1_area_normal_diff(const matrix<double> &node,
                                 const matrix<size_t> &tri,
                                 const matrix<double> &areas,
                                 const matrix<double> &target_normal,
                                 double normal_align_w)
{
  const matrixd &surface_area = areas;

  const double total_area =
      std::accumulate(surface_area.begin(), surface_area.end(), 0.0);
  cerr << "total_area in build_smooth_L1_normal: " << total_area << endl;

  const size_t fn = tri.size(2);
  vector<shared_ptr<jtf::function::functionN1_t<double,int32_t> > > sum(fn*3);
  const double weight = normal_align_w / total_area;
  for(size_t fi = 0; fi < fn; ++fi) {
      for(size_t i = 0; i < 3; ++i) {
          shared_ptr<jtf::function::functionN1_t<double,int32_t> > normal_func(new face_area_normal_diff(node, tri(colon(), fi), target_normal(colon(),fi), i, 1));
          sum[fi*3+i].reset(new smooth_L1(normal_func, weight, L1_sqrt_eps*surface_area[fi]));
        }
    }

  unique_ptr<jtf::function::functionN1_t<double,int32_t> > sum_func(new jtf::function::sum_function<double,int32_t, jtf::function::SMART_STD>(sum));

  return sum_func.release();
}

hj::function::function_t<double, int32_t> *
build_fix_normal_polycube_function_control(const jtf::tet_mesh &tm)
{
  //////////////////////////////////////////////////////////////////
  ////////////////  add function term  /////////////////////////////

  boost::shared_ptr<function_t<double, int32_t> > tetmesh_arap(
        build_tetmesh_arap_func(tm.tetmesh_.node_, tm.tetmesh_.mesh_, arap_w_control));

  matrix<double> zero_pos = zeros<double>(3,1);
  boost::shared_ptr<function_t<double, int32_t> > fix_pos(
        new fix_zero_node_func(tm.tetmesh_.node_, zero_pos, 1));

  boost::shared_ptr<vector<boost::shared_ptr<function_t<double, int32_t> > > >
      funcs(new vector<boost::shared_ptr<function_t<double, int32_t> > >);

  if(tetmesh_arap.get() && non_zero(arap_w_control)){
      funcs->push_back(tetmesh_arap);
      cerr << "# [info] add jtf::mesh::mesh arap function, weight " << arap_w_control << endl;
    }

  //  if(fix_pos.get()){
  //      funcs->push_back(fix_pos);
  //      cerr << "# [info] add fix first node zero function" << endl;
  //    }

  return new_catenated_function<double, int32_t>(funcs);
}

int polycube_fix_normal(ptree &pt)
{
  jtf::tet_mesh tm(pt.get<string>("input/tet.value").c_str());

  pt.put("weight/normal_align_w.desc", "surface L1 weight, default is 1e2.");
  pt.put("wegiht/fix_zero_w.desc", "fix first node to zero, default is 1");

  double normal_align_w = pt.get<double>("weight/normal_align_w.value", 0.1);

  matrix<double> target_normal(3, tm.outside_face_.size(2));
  if(jtf::mesh::read_matrix(pt.get<string>("input/normal.value").c_str(), target_normal)){
      cerr << "# [error] can not read normal matrix." << endl;
      return __LINE__;
    }

  L1_sqrt_eps = 0.5;
  arap_w_control = 1.0;


  const size_t iter_w = pt.get<size_t>("weight/iter_num.value");
  for(size_t i = 0; i < iter_w; ++i){
      shared_ptr<function_t<double, int32_t> > func(build_fix_normal_polycube_function_control(tm));
      vector<shared_ptr<jtf::function::functionN1_t<double,int32_t> > > sum;
      sum.push_back(shared_ptr<jtf::function::functionN1_t<double,int32_t> >(jtf::function::least_square_warpper(func)));
      sum.push_back(shared_ptr<jtf::function::functionN1_t<double,int32_t> >(
                      build_smooth_L1_area_normal_diff(tm.tetmesh_.node_, tm.outside_face_,
                                                       tm.outside_face_area_, target_normal, normal_align_w)));

      shared_ptr<jtf::function::functionN1_t<double,int32_t> > target(new jtf::function::sum_function<double,int32_t, jtf::function::SMART_STD>(sum));
      jtf::optimize(*target,tm.tetmesh_.node_,pt,nullptr, nullptr, nullptr);
      normal_align_w *= 2;
      {
        stringstream ss;
        ss << "deform_" << i <<".vtk";
        ofstream ofs(ss.str().c_str());
        tet2vtk(ofs, &tm.tetmesh_.node_[0], tm.tetmesh_.node_.size(2), &tm.tetmesh_.mesh_[0], tm.tetmesh_.mesh_.size(2));
      }
    }
  jtf::mesh::tet_mesh_write_to_zjumat(pt.get<string>("output/tet.value").c_str(),
                                      &tm.tetmesh_.node_, &tm.tetmesh_.mesh_);

  cerr << "success." << endl;
  return 0;
}
