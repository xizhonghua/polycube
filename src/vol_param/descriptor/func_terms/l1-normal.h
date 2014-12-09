#ifndef L1NORMAL_H
#define L1NORMAL_H

#include <map>
#include <zjucad/matrix/matrix.h>
#include <zjucad/matrix/itr_matrix.h>
#include <hjlib/math/blas_lapack.h>
#include <zjucad/matrix/lapack.h>
#include <hjlib/sparse/sparse.h>
#include <jtflib/function/func_aux.h>
#include "../def.h"
#include "../../common/util.h"
#include "../../../mesh_func/tri-area-normal.h"

#define USE_SQRT 1

class smooth_L1 : public jtf_func
{
public:
  // NOTICE: assume src provide sparse gradient and hessian
  smooth_L1(jtf_func_ptr &src, double w, double eps=1e-2)
    :src_(src), w_(w), eps_(eps) {
  }
  virtual size_t dim(void) const {
    return src_->dim();
  }
  virtual int val(const double *x, double &v) {
    double t = 0;
    src_->val(x, t);

#if USE_SQRT
    const double eps2_ = pow(eps_ * L1_sqrt_eps_,2);
    v += sqrt(t*t+eps2_)*w_ * scalar_;
#else
    if(fabs(t) > eps_)
      v += (fabs(t)-eps_/2)*w_;
    else
      v += t*t/(2*eps_)*w_;
#endif
    return 0;
  }

  //! @param nnz=0 means dense, g=0 or h=0 means query nnz
  virtual int gra(const double *x, double *g) {
    double t = 0;
    src_->val(x, t);
    //const double L1_x = t/sqrt(t*t+eps2_);
    const double eps2_ = pow(eps_ * L1_sqrt_eps_,2);
    double L1_x = t/sqrt(t*t+eps2_);
    jtf::math::erase_nan_inf(L1_x);
    eval_g(x);
    for(size_t xi = 0; xi < g_idx_.size(); ++xi) {
#if USE_SQRT
        g[g_idx_[xi]] += g_[xi]*L1_x*w_ * scalar_;
#else
        if(fabs(t) > eps_) {
            g[g_idx_[xi]] += g_[xi]*((t>0)?1:-1)*w_*scalar_;
          }
        else {
            g[g_idx_[xi]] += g_[xi]*(t/eps_)*w_*scalar_;
          }
#endif
      }
    return 0;
  }

  //! @param nnz=0 means dense, g=0 or h=0 means query nnz
  virtual int gra(const double *x, size_t &nnz , double *g, int32_t *idx) {

    eval_g(x);

    if(g == 0 && idx == 0){
        nnz = g_idx_.size();
        return 0;
      }else{
        assert(g != 0 && idx != 0);
        double t = 0;
        src_->val(x, t);
        //const double L1_x = t/sqrt(t*t+eps2_);
        const double eps2_ = pow(eps_ * L1_sqrt_eps_,2);
        double L1_x = t/sqrt(t*t+eps2_);
        jtf::math::erase_nan_inf(L1_x);

        for(size_t xi = 0; xi < g_idx_.size(); ++xi) {
            idx[xi] = g_idx_[xi];

#if USE_SQRT
            g[xi] = g_[xi]*L1_x*w_ * scalar_;
#else
            if(fabs(t) > eps_) {
                g[xi] = g_[xi]*((t>0)?1:-1)*w_*scalar_;
              }
            else {
                g[xi] += g_[xi]*(t/eps_)*w_*scalar_;
              }
#endif
          }
      }
    return 0;
  }

  virtual int hes(const double *x, size_t &nnz, size_t &format, double *h,
                  int32_t *ptr, int32_t *idx, double alpha) {
    format = 2;
    eval_g(x);
    if(h == 0 && ptr == 0 && idx == 0) { // query nnz
        std::vector<int32_t> ptr_, idx_;
        size_t format = -1;
        src_->hes(x, nnz, format, 0, 0, 0);
        ptr_.resize(dim()+1);
        ptr_[0] = 0;
        idx_.resize(nnz);
        h_.resize(nnz);
        src_->hes(x, nnz, format, 0, &ptr_[0], &idx_[0]);

        if(nnz > g_idx_.size()*g_idx_.size()) {
            // TODO: this acceleration need seriuos check, assuming ggt cover all the entries
            for(size_t xi = 0; xi < dim(); ++xi) {
                for(size_t nzi = ptr_[xi]; nzi < ptr_[xi+1]; ++nzi) {
                    pattern_[xi].insert(idx_[nzi]);
                  }
              }
          }
        for(size_t gci = 0; gci < g_idx_.size(); ++gci) {
            for(size_t gri = 0; gri < g_idx_.size(); ++gri) {
                pattern_[g_idx_[gci]].insert(g_idx_[gri]);
              }
          }
        nnz = 0;
        for(std::map<int32_t, std::set<int32_t> >::const_iterator xi = pattern_.begin();
            xi != pattern_.end(); ++xi) {
            nnz += xi->second.size();
          }
        return 0;
      }
    if(h == 0 && ptr != 0 && idx != 0) { // query pattern
        size_t nzi = 0;
        for(std::map<int32_t, std::set<int32_t> >::const_iterator xi = pattern_.begin();
            xi != pattern_.end(); ++xi) {
            for(std::set<int32_t>::const_iterator ri = xi->second.begin();
                ri != xi->second.end(); ++ri, ++nzi) {
                ptr[nzi] = xi->first;
                idx[nzi] = *ri;
              }
          }

        std::map<int32_t, std::set<int32_t> > tmp;
        swap(pattern_, tmp);
        return 0;
      }
    if(h != 0 && ptr != 0 && idx != 0) { // accumulate
        double t = 0;
        src_->val(x, t);
#if USE_SQRT
        const double eps2_ = pow(eps_ * L1_sqrt_eps_,2);
        const double L1v = sqrt(t*t+eps2_);
        double L1_x = t/L1v, L1_xx = eps2_/(L1v*L1v*L1v);
        jtf::math::erase_nan_inf(L1_x);
        jtf::math::erase_nan_inf(L1_xx);
        zjucad::matrix::matrix<double> H = (L1_xx*w_*scalar_)*g_*zjucad::matrix::trans(g_);
        if(src_->hes_block(x, &H[0], L1_x*w_*scalar_)) {
            std::cerr << "hes block error." << std::endl;
            exit(0);
          }

        zjucad::matrix::matrix<double> e(H.size(1)), diag_e =
            zjucad::matrix::zeros<double>(H.size(1), H.size(1));
        zjucad::matrix::eig(H, e);
        for(size_t ei = 0; ei < e.size(); ++ei) {
            if(e[ei] > 0)
              diag_e(ei, ei) = e[ei];
            else
              diag_e(ei, ei) = 0;
          }
        H = zjucad::matrix::temp(H*zjucad::matrix::temp(diag_e*zjucad::matrix::trans(H)));
        for(size_t gci = 0; gci < g_idx_.size(); ++gci) {
            for(size_t gri = 0; gri < g_idx_.size(); ++gri) {
                if(jtf::function::add_to_csc(h, ptr, idx, g_idx_[gri], g_idx_[gci], H(gci, gri)))
                  return __LINE__;
              }
          }
#else
        if(fabs(t) > eps_) {
            src_->hes(x, nnz, &h[0], ptr, idx, ((t>0)?1:-1)*w_);
          }
        else {
            src_->hes(x, nnz, &h[0], ptr, idx, t/eps_*w_);
            for(size_t gci = 0; gci < g_idx_.size(); ++gci) {
                for(size_t gri = 0; gri < g_idx_.size(); ++gri) {
                    if(hj_func_opt::add_to_csc(h, ptr, idx, g_idx_[gri], g_idx_[gci], g_[gci]*g_[gri]*w_/eps_))
                      return __LINE__;
                  }
              }
          }
#endif
        return 0;
      }
    return __LINE__;
  }
  virtual int hes_block(const double *x, double *h, double alpha = 1) {return -1;}
protected:
  int eval_g(const double *x) {
    if(g_.size() == 0) {
        size_t nnz = 0;
        src_->gra(x, nnz, 0, 0);
        g_.resize(nnz, 1);
        g_idx_.resize(nnz);
      }
    std::fill(g_.begin(), g_.end(), 0);
    size_t nnz = g_.size();
    src_->gra(x, nnz, &g_[0], &g_idx_[0]);
    return 0;
  }
  jtf_func_ptr src_;
  zjucad::matrix::matrix<double> g_, h_;
  std::vector<int32_t> g_idx_;
  std::map<int32_t, std::set<int32_t> > pattern_;
  const double  w_, eps_;
public:
  static double scalar_; // usr tune;
  static double L1_sqrt_eps_;
};

inline void set_l1_normal_align_weight(double &v)
{
  smooth_L1::scalar_ = v;
}

inline void set_l1_eps_weight(double &v)
{
  smooth_L1::L1_sqrt_eps_ = v;
}

inline double get_l1_normal_align_weight()
{
  return smooth_L1::scalar_;
}
inline double get_l1_eps_weight()
{
  return smooth_L1::L1_sqrt_eps_;
}

class face_area_normal : public jtf_func
{
public:
  face_area_normal(const zjucad::matrix::matrix<double> &node,
                   const zjucad::matrix::matrix<size_t> &tri,
                   size_t xyz,
                   const double weight,
                   const node_mapping * NM = 0)
    :tri_(tri), node_num_(node.size(2)), weight_(weight), xyz_(xyz),
      node_mapping_(NM){
    if(node_mapping_){
        get_node_mapping(tri_, node_mapping_->ZT, node_mapping_of_each_variabe_);

        std::set<size_t> real_variables;
        for(const auto & one_node: node_mapping_of_each_variabe_){
            for(const auto & one_exp : one_node){
                real_variables.insert(one_exp.first);
              }
          }
        nnz_of_node_mapping_ = real_variables.size();

        for(size_t ni = 0; ni < node_mapping_of_each_variabe_.size(); ++ni){
            const std::vector<std::pair<size_t,double> > & one_eqn =
                node_mapping_of_each_variabe_[ni];
            for(const auto & one_exp : one_eqn){
                real_v_to_orig_v_coeff_[one_exp.first].push_back(
                      std::make_pair(ni, one_exp.second));
              }
          }
      }
  }
  virtual size_t dim(void) const {
    if(!node_mapping_)
      return node_num_*3;
    else
      return node_mapping_->ZT.size(1);
  }
  virtual int val(const double *x, double &v) {
    zjucad::matrix::matrix<double> tri;
    if(!node_mapping_){
        const zjucad::matrix::itr_matrix<const double *> T(3, node_num_, x);
        tri = T(zjucad::matrix::colon(), tri_);
      }else{
        tri.resize(3, tri_.size());
        get_cell_node(x, tri_, *node_mapping_, tri);
      }
    zjucad::matrix::matrix<double> f(3,1);
    calc_tri_area_normal_(&f[0], &tri[0]);
    jtf::math::erase_nan_inf(f[xyz_]);
    v += f[xyz_] * weight_;
    return 0;
  }

  virtual int gra(const double *x, double *g) {
    if(!node_mapping_){
        double sp_g[9];
        int32_t idx[9];
        size_t nnz = 9;
        gra(x, nnz, sp_g, idx);
        for(size_t i = 0; i < 9; ++i)
          g[idx[i]] += sp_g[i];
      }else{
        double sp_g[nnz_of_node_mapping_];
        int32_t idx[nnz_of_node_mapping_];
        size_t nnz = nnz_of_node_mapping_;
        gra(x, nnz, sp_g, idx);
        for(size_t i = 0; i < nnz_of_node_mapping_; ++i)
          g[idx[i]] += sp_g[i];
      }
    return 0;
  }

  virtual int gra(const double *x, size_t &nnz, double *g, int32_t *idx) {
    if(g == 0) {
        if(!node_mapping_)
          nnz = 9;
        else
          nnz = nnz_of_node_mapping_;
        return 0;
      }

    zjucad::matrix::matrix<double> tri, jac(3,9) ;
    if(!node_mapping_){
        const zjucad::matrix::itr_matrix<const double *> T(3, node_num_, x);
        tri = T(zjucad::matrix::colon(), tri_);
      }else{
        tri.resize(3, tri_.size());
        get_cell_node(x, tri_, *node_mapping_, tri);
      }
    calc_tri_area_normal_jac_(&jac[0], &tri[0]);
    if(!node_mapping_){
        for(size_t i = 0; i < 9; ++i) {
            idx[i] = tri_[i/3]*3+i%3;
            g[i] = jac(xyz_, i)*weight_;
            jtf::math::erase_nan_inf(g[i]);
          }
      }else{
        assert(nnz_of_node_mapping_ == real_v_to_orig_v_coeff_.size());
        size_t i = 0;
        for(const auto & one_real_v : real_v_to_orig_v_coeff_){
            idx[i] = one_real_v.first;
            g[i] = 0;
            const std::vector<std::pair<size_t,double> > & orig_v_coeff = one_real_v.second;
            for(const auto & one_real_v_coeff : orig_v_coeff){
                g[i] += jac(xyz_, one_real_v_coeff.first) * one_real_v_coeff.second * weight_;
              }
            ++i;
          }
      }
    return 0;
  }
  virtual int hes(const double *x, size_t &nnz, size_t &format, double *h, int32_t *ptr, int32_t *idx, double alpha = 1) {
    format = 1;
    if(h == 0 && ptr == 0 && idx == 0) { // query nnz
        if(!node_mapping_)
          nnz = 9*9;
        else
          nnz = nnz_of_node_mapping_*nnz_of_node_mapping_;
        return 0;
      }
    if(h == 0 && ptr != 0 && idx != 0) { // query pattern
        if(!node_mapping_){
            std::vector<size_t> ptr_vec;
            ptr_vec.push_back(0);

            for(size_t ci = 0; ci < 9; ++ci) {

                if(tri_[ci/3]*3+ci%3 != ptr_vec.back()+1){
                    ptr[tri_[ci/3]*3+ci%3] = ptr[ptr_vec.back()+1];
                  }
                ptr_vec.push_back(tri_[ci/3]*3+ci%3);

                ptr[tri_[ci/3]*3+ci%3+1] = ptr[tri_[ci/3]*3+ci%3]+9;
                for(size_t ri = 0; ri < 9; ++ri) {
                    idx[ptr[tri_[ci/3]*3+ci%3]+ri] = tri_[ri/3]*3+ri%3;
                  }
              }
          }else{// with node_mapping
            assert(nnz_of_node_mapping_ == real_v_to_orig_v_coeff_.size());

            std::vector<size_t> ptr_vec;
            ptr_vec.push_back(0);

            for(const auto & one_real_v : real_v_to_orig_v_coeff_){
                ptr[one_real_v.first+1] = ptr[one_real_v.first] + nnz_of_node_mapping_;
                if(one_real_v.first != ptr_vec.back()+1){
                    ptr[one_real_v.first] = ptr[ptr_vec.back()+1];
                  }
                ptr_vec.push_back(one_real_v.first);

                size_t ri = 0;
                for(const auto & one_real_v_ri : real_v_to_orig_v_coeff_){
                    idx[ptr[one_real_v.first] + ri] = one_real_v_ri.first;
                    ++ri;
                  }
              }
          }
        return 0;
      }
    if(h != 0 && ptr != 0 && idx != 0) { // accumulate

        zjucad::matrix::matrix<double> tri, hes(27, 9);// =
        if(!node_mapping_){
            const zjucad::matrix::itr_matrix<const double *> T(3, node_num_, x);
            tri = T(zjucad::matrix::colon(), tri_);
          }else{
            tri.resize(3, tri_.size());
            get_cell_node(x, tri_, *node_mapping_, tri);
          }
        if(is_degenerate(tri)) return 0;
        calc_tri_area_normal_hes_(&hes[0], &tri[0]);
        zjucad::matrix::matrix<double> H(9, 9);
        for(size_t ci = 0; ci < 9; ++ci) {
            for(size_t ri = 0; ri < 9; ++ri) {
                H(ri, ci) = hes(xyz_+ri*3, ci)*weight_*alpha;
              }
          }
        const zjucad::matrix::matrix<double> save_H = H;
        zjucad::matrix::matrix<double> e(9), diag_e = zjucad::matrix::zeros<double>(9, 9);
        zjucad::matrix::eig(H, e);
        for(size_t ei = 0; ei < 9; ++ei) {
            if(e[ei] > 0)
              diag_e(ei, ei) = e[ei];
            else
              diag_e(ei, ei) = 0;
          }
        H = zjucad::matrix::temp(H*zjucad::matrix::temp(diag_e*zjucad::matrix::trans(H)));
        if(!node_mapping_){
            for(size_t ci = 0; ci < 9; ++ci) {
                const size_t var_ci = tri_[ci/3]*3+ci%3;
                for(size_t ri = 0; ri < 9; ++ri) {
                    const size_t var_ri =tri_[ri/3]*3+ri%3;
                    if(jtf::function::add_to_csc(h, ptr, idx, var_ri, var_ci, H(ri, ci))) {
                        return __LINE__;
                      }
                  }
              }
          }else{
            zjucad::matrix::matrix<double> ZTHZ =
                zjucad::matrix::matrix<double>(nnz_of_node_mapping_,
                                               nnz_of_node_mapping_);
            get_ZTHZ(x, H, ZTHZ);

            size_t ri = 0, ci = 0;
            for(const auto & ri_variable : real_v_to_orig_v_coeff_){
                const size_t variable_ri = ri_variable.first;
                for(const auto & ci_variable : real_v_to_orig_v_coeff_){
                    const size_t variable_ci = ci_variable.first;
                    if(jtf::function::add_to_csc(h,ptr,idx, variable_ri, variable_ci, ZTHZ(ri,ci)))
                      return __LINE__;
                    ++ci;
                  }
                ++ri;
              }
          }
        return 0;
      }
    return __LINE__;
  }
  // accumulate
  virtual int hes_block(const double *x, double *h, double alpha = 1) {
    zjucad::matrix::matrix<double> tri, hes(27, 9);
    if(!node_mapping_){
        const zjucad::matrix::itr_matrix<const double *> T(3, node_num_, x);
        tri = T(zjucad::matrix::colon(), tri_);
      }else{
        tri.resize(3, tri_.size());
        get_cell_node(x, tri_, *node_mapping_, tri);
      }
    if(is_degenerate(tri)) return 0;
    calc_tri_area_normal_hes_(&hes[0], &tri[0]);
    if(!node_mapping_){
        zjucad::matrix::itr_matrix<double *> H(9, 9, h);
        for(size_t ci = 0; ci < 9; ++ci) {
            //        const size_t var_ci = tri_[ci/3]*3+ci%3;
            for(size_t ri = 0; ri < 9; ++ri) {
                //            const size_t var_ri = tri_[ri/3]*3+ri%3;
                H(ri, ci) += hes(xyz_+ri*3, ci)*weight_*alpha;
              }
          }
      }else{
        zjucad::matrix::itr_matrix<double *> H(nnz_of_node_mapping_, nnz_of_node_mapping_, h);
        zjucad::matrix::matrix<double> H_orig(9,9);
        for(size_t ci = 0; ci < 9; ++ci) {
            for(size_t ri = 0; ri < 9; ++ri) {
                H_orig(ri, ci) += hes(xyz_+ri*3, ci)*weight_*alpha;
              }
          }
        get_ZTHZ(x, H_orig, H);
      }
    return 0;
  }

  template <typename T1, typename T2>
  void get_ZTHZ(const double *x,
                const zjucad::matrix::matrix_expression<T1> & H,
                zjucad::matrix::matrix_expression<T2> &ZTHZ)const{
    if(ZTHZ().size(1) != nnz_of_node_mapping_ || ZTHZ().size(2) != nnz_of_node_mapping_)
      ZTHZ() = zjucad::matrix::zeros<double>(nnz_of_node_mapping_,nnz_of_node_mapping_);

    zjucad::matrix::matrix<double> ZTH = zjucad::matrix::zeros<double>(nnz_of_node_mapping_, 9);
    assert(nnz_of_node_mapping_ == real_v_to_orig_v_coeff_.size());
    size_t ri = 0, ci = 0;
    for(const auto & one_real_v : real_v_to_orig_v_coeff_){
        const std::vector<std::pair<size_t,double> > & to_orig_v = one_real_v.second;
        for(size_t ci = 0; ci < H().size(2); ++ci){
            for(const auto & one_orig_v : to_orig_v){
                ZTH(ri,ci) +=
                    one_orig_v.second *
                    H()(one_orig_v.first, ci);
              }
          }
        ++ri;
      }

    // get ZTHZ from ZTH and Z
    ri = 0;
    for(const auto & one_real_v : real_v_to_orig_v_coeff_){
        const std::vector<std::pair<size_t,double> > &to_orig_v = one_real_v.second;
        for(size_t ri = 0; ri < ZTH.size(1); ++ri){
            for(const auto & one_exp : to_orig_v){
                ZTHZ()(ri,ci) += ZTH(ri, one_exp.first)* one_exp.second;
              }
          }
        ++ci;
      }
  }
private:
  const size_t node_num_, xyz_;
  const zjucad::matrix::matrix<size_t> tri_;
  const double weight_;
  const node_mapping * node_mapping_;
  std::vector<std::vector<std::pair<size_t,double> > > node_mapping_of_each_variabe_;
  size_t nnz_of_node_mapping_;
  std::map<size_t, std::vector<std::pair<size_t,double> > > real_v_to_orig_v_coeff_;
};


template <typename T1, typename T2>
jtf_func_ptr
build_smooth_L1_area_normal(
    const zjucad::matrix::matrix_expression<T1> & node,
    const zjucad::matrix::matrix_expression<T2> &tri,
    const zjucad::matrix::matrix_expression<T1> &areas,
    const node_mapping * node_mapping = 0)
{
  using namespace zjucad::matrix;
  using namespace std;

  const zjucad::matrix::matrix<double> &surface_area = areas;

  const double total_area =
      std::accumulate(surface_area.begin(), surface_area.end(), 0.0);
  cerr << "total_area in build_smooth_L1_normal: " << total_area << endl;

  const size_t fn = tri().size(2);
  vector<jtf_func_ptr>  sum(fn*3);
  const double weight = 1 / total_area;
  for(size_t fi = 0; fi < fn; ++fi) {
      for(size_t i = 0; i < 3; ++i) {
          jtf_func_ptr normal_func(
                new face_area_normal(node(), tri()(colon(), fi), i, 1, node_mapping));
          sum[fi*3+i].reset(new smooth_L1(normal_func, weight, surface_area[fi]));
        }
    }
  jtf_func_ptr rtn(new jtf::function::sum_function<double,int32_t, jtf::function::SMART_STD>(sum));

  return rtn;
}

#endif // L1NORMAL_H
