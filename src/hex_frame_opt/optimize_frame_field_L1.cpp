#include <jtflib/function/function.h>
#include <jtflib/optimizer/optimizer.h>
#include "../tetmesh/tetmesh.h"
#include <zjucad/ptree/ptree.h>

#include <hjlib/math/polar.h>
#include "permutation.h"
#include "../vol_param/descriptor/func_terms/l1-normal.h"

using namespace std;
using boost::property_tree::ptree;
using namespace zjucad::matrix;

typedef jtf::function::functionN1_t<double,int32_t> jtf_func;
typedef std::shared_ptr<const jtf_func> jtf_func_cons_ptr;

typedef hj::function::function_t<double,int32_t> hj_func;
typedef std::shared_ptr<const hj_func> hj_func_cons_ptr;

double smooth_L1::scalar_ = 1.0;
double smooth_L1::L1_sqrt_eps_ = 0.5;

class normal_align_func : public jtf_func
{
public:
  normal_align_func(const size_t & tet_num,
                    const size_t idx,
                    const zjucad::matrix::matrix<double> & normal,
                    const size_t di,
                    const double weight)
    :tet_num_(tet_num), idx_(idx), normal_(normal), di_(di), weight_(weight){}
  virtual ~normal_align_func(){}

  virtual size_t dim(void) const{
    return tet_num_ * 9;
  }
  virtual int val(const double *x, double &v){
    using namespace zjucad::matrix;
    itr_matrix<const double*> x0(9, tet_num_, x);
    itr_matrix<const double*> F(3,3, &x0(0,idx_));
//    matrix<double> dn = F(colon(),di_);
//    double f = dot(dn, normal_);
    matrix<double> dn = F*normal_;
    double f = dn[di_];
    jtf::math::erase_nan_inf(f);
    v += f * weight_;
    return 0;
  }
  virtual int gra(const double *x, double *g){
    //    size_t nnz;
    //    gra(x,nnz,0,0);
    //    vector<double> gra_(nnz);
    //    vector<int32_t> idx_(nnz);
    //    gra(x,nnz, &gra_[0], &idx_[0]);
    //    for(size_t i = 0; i < nnz; ++i){
    //        g[idx_[i]] += gra_[i];
    //      }

    //////////////////
    ///  speed up ////
    static double gra_[3];
    static int32_t idx_[3];
    static size_t nnz = 3;
    gra(x,nnz, &gra_[0], &idx_[0]);
    for(size_t i = 0; i < nnz; ++i){
        g[idx_[i]] += gra_[i];
      }
    return 0;
  }
  virtual int gra(const double *x, size_t &nnz, double *g, int32_t *idx)
  {
    if(g == 0 && idx == 0){
        nnz = 3;
        return 0;
      }
    assert(g != 0 && idx != 0);
    for(size_t i = 0; i< 3; ++i){
        idx[i] = 9 * idx_ + i + 3*di_;
        g[i] = normal_[i] * weight_;
        jtf::math::erase_nan_inf(g[i]);
      }
    return 0;
  }
  virtual int hes(const double *x, size_t & nnz, size_t &format,
                  double *h, int32_t*ptr, int32_t *idx, double alpha = 1){
    nnz = 0;
    return 0;
  }
  virtual int hes_block(const double *x, double *h, double alpha = 1) {return 0;}
private:
  const size_t tet_num_;
  const size_t idx_;
  const size_t di_;
  const double weight_;
  const zjucad::matrix::matrix<double> normal_;
};


class direction_align : public jtf_func
{
public:
  direction_align(const size_t tet_num,  const size_t idx0,
                  const size_t di, // in idx0 frame
                  const size_t idx1, const size_t dj,
                  const double w)
    :tet_num_(tet_num), idx0_(idx0), idx1_(idx1), di_(di), dj_(dj), w_(w){
    if(idx0 > idx1){
        swap(idx0_,idx1_);
        swap(di_,dj_);
      }
  }
  virtual ~direction_align(){}

  virtual size_t dim(void) const{
    return 9 * tet_num_;
  }
  virtual int val(const double *x, double &v)
  {
    using namespace zjucad::matrix;
    itr_matrix<const double*> x0(9, tet_num_, x);
    itr_matrix<const double*> F0(3,3, &x0(0,idx0_));
    itr_matrix<const double*> F1(3,3, &x0(0,idx1_));
    matrix<double> di = F0(colon(),di_);
    matrix<double> dj = F1(colon(),dj_);
    double dot_v = dot(di, dj);
    jtf::math::erase_nan_inf(dot_v);
    v += dot_v * w_;
    return 0;
  }
  virtual int gra(const double *x, double *g){
    //    size_t nnz;
    //    gra(x,nnz, 0,0);
    //    vector<double> gra_(nnz);
    //    vector<int32_t> idx_(nnz);
    //    gra(x, nnz, &gra_[0], &idx_[0]);
    //    for(size_t i = 0; i < nnz; ++i)
    //      g[idx_[i]] += gra_[i];

    size_t nnz = 6;
    static vector<double> gra_(6);
    static vector<int32_t> idx_(6);
    gra(x, nnz, &gra_[0], &idx_[0]);
    for(size_t i = 0; i < nnz; ++i)
      g[idx_[i]] += gra_[i];
    return 0;
  }
  virtual int gra(const double *x, size_t & nnz, double *g, int32_t *idx){
    using namespace zjucad::matrix;
    if(g == 0 && idx == 0){
        nnz = 6;
        return 0;
      }
    assert(g != 0 && idx != 0);

    for(size_t i = 0; i < 6; ++i){
        idx[i] = (i<3?9*idx0_+3*di_+i%3:9*idx1_+3*dj_+i%3);
        g[i] = (i<3?x[9*idx1_+3*dj_+i%3]:x[9*idx0_+3*di_+i%3]);
        g[i] *= w_;
      }

    return 0;
  }
  virtual int hes(const double *x, size_t &nnz, size_t &format, double *h,
                  int32_t *ptr, int32_t *idx, double alpha = 1)
  {
    if(h == 0 && ptr == 0 && idx == 0){
        nnz = 6; format = 1;
        return 0;
      }
    if(h == 0 && idx != 0 && ptr != 0){
        ptr_vec_.clear();
        ptr_vec_.push_back(0);
        for(size_t i = 0; i < 6; ++i){
            const size_t pidx = i<3?9*idx0_+3*di_+i%3:9*idx1_+3*dj_+i%3;
            if(pidx != ptr_vec_.back()+1){
                ptr[pidx] = ptr[ptr_vec_.back()+1];
              }
            ptr_vec_.push_back(pidx);
            ptr[pidx+1]= ptr[pidx] + 1;
            idx[ptr[pidx]]= i<3?9*idx1_+3*dj_+i%3:9*idx0_+3*di_+i%3;
          }
        return 0;
      }
    if(h != 0 && ptr != 0 && idx != 0){
        for(size_t i = 0; i < 6; ++i){
            const size_t ri = i<3?9*idx0_+3*di_+i%3:9*idx1_+3*dj_+i%3;
            const size_t ci = i<3?9*idx1_+3*dj_+i%3:9*idx0_+3*di_+i%3;
            if(jtf::function::add_to_csc(h,ptr, idx, ri,ci,w_))
              return __LINE__;
          }
      }
    return 0;
  }
  virtual int hes_block(const double *x, double *h, double alpha = 1) {
    using namespace zjucad::matrix;
    itr_matrix<double*> H(6,6,h);
    for(size_t i = 0; i < 6; ++i){
        H(i,i) += w_;
      }
    return 0;
  }
private:
  const size_t tet_num_;
  size_t idx0_,idx1_;
  size_t di_, dj_;
  const double w_;
  vector<size_t> ptr_vec_;
};

class permutation_matrix_euler_func : public jtf_func
{
public:
  permutation_matrix_euler_func(const size_t tet_num,
                                const size_t idx0, const size_t idx1,
                                size_t order, const double w)
    :tet_num_(tet_num), idx0_(idx0), idx1_(idx1), order_(order), w_(w){}
  virtual ~permutation_matrix_euler_func(){}
  virtual size_t dim(void) const{}
  virtual int val(const double *x, double &v) {}
  virtual int gra(const double *x, double *g){}
  virtual int gra(const double *x, size_t & nnz, double *g, int32_t *idx){}
  virtual int hes(const double *x, size_t &nnz, size_t &format, double *h,
                  int32_t *ptr, int32_t *idx, double alpha=1){}
  virtual int hes_block(const double *x, double *h, double alpha =1){}
private:
  const size_t tet_num_;
  const size_t idx0_,idx1_, order_;
  const double w_;
};

class permutation_matrix_func : public jtf_func
{
public:
  permutation_matrix_func(const size_t tet_num, const size_t idx0,
                          const size_t idx1, size_t order,const double w)
    :tet_num_(tet_num),idx0_(idx0),idx1_(idx1), order_(order), w_(w){
    if(idx0_ > idx1_){
        swap(idx0_, idx1_);
      }
    assemble();}
  virtual ~permutation_matrix_func(){}
  virtual size_t dim(void)const{
    return 9 * tet_num_;
  }
  virtual int val(const double *x, double &v) {
    itr_matrix<const double*> x0(9,tet_num_, x);
    itr_matrix<const double*> F0(3,3, &x0(0,idx0_));
    itr_matrix<const double*> F1(3,3, &x0(0,idx1_));

    static matrix<double> two_frame(3,6);
    two_frame(colon(),colon(0,2)) = F0;
    two_frame(colon(),colon(3,5)) = F1;

    double f = 0;
    calc_permutation(&f, &two_frame[0], order_);

    v += f*w_;

    return 0;
  }
  virtual int gra(const double *x, double *g)
  {
    //    size_t nnz = 0;
    //    gra(x, nnz, 0,0);
    //    vector<double> g_(nnz);
    //    vector<int32_t> idx_(nnz);
    //    return gra(x,nnz, &g_[0], &idx_[0]);

    size_t nnz = 12;
    static vector<double> g_(nnz);
    static vector<int32_t> idx_(nnz);
    return gra(x,nnz, &g_[0], &idx_[0]);
  }
  virtual int gra(const double *x, size_t &nnz, double *g, int32_t *idx){
    if(g == 0 && idx == 0){
        nnz = 12;
        return 0;
      }
    assert(g != 0 && idx != 0);

    itr_matrix<const double*> x0(9,tet_num_, x);
    itr_matrix<const double*> F0(3,3, &x0(0,idx0_));
    itr_matrix<const double*> F1(3,3, &x0(0,idx1_));

    static matrix<double> two_frame(3,6);
    two_frame(colon(),colon(0,2)) = F0;
    two_frame(colon(),colon(3,5)) = F1;
    static matrix<double> gra_(18,1);
    calc_permutation_jac(&gra_[0], &two_frame[0], order_);

    itr_matrix<double*> g_m(12,1,g);
    g_m = gra_(idx_in_18_);
    itr_matrix<int32_t*> idx_m(12,1,idx);
    idx_m = idx_used_;
    return 0;
  }
  virtual int hes(const double *x, size_t &nnz, size_t &format, double *h,
                  int32_t *ptr, int32_t *idx, double alpha = 1)
  {
    if(h == 0 && ptr == 0 && idx == 0){
        nnz = 144; format = 1;
        return 0;
      }
    if(h == 0 && ptr != 0 && idx != 0){
        ptr_vec_.clear();
        ptr_vec_.push_back(0);

        for(size_t i = 0; i < 12; ++i){

            if(idx_used_[i] != ptr_vec_.back()+1){
                assert(i < idx_used_.size() && idx_used_[i] < dim()+1);
                assert(ptr_vec_.back()+1 < dim()+1);
                ptr[idx_used_[i]] = ptr[ptr_vec_.back()+1];
              }
            ptr_vec_.push_back(idx_used_[i]);

            assert(idx_used_[i] < dim()+1);
            ptr[idx_used_[i]+1] = ptr[idx_used_[i]] + 12;

            for(size_t j = 0; j < 12; ++j){
                assert(j < 12);
                assert(idx_used_[i] < dim()+1);
                assert(ptr[idx_used_[i]]+j < nnz) ;
                idx[ptr[idx_used_[i]]+j] = idx_used_[j];
              }
          }
        return 0;
      }
    if(h != 0 && ptr !=0 && idx != 0){
        itr_matrix<const double*> x0(9,tet_num_, x);
        itr_matrix<const double*> F0(3,3, &x0(0,idx0_));
        itr_matrix<const double*> F1(3,3, &x0(0,idx1_));

        static matrix<double> two_frame(3,6);
        two_frame(colon(),colon(0,2)) = F0;
        two_frame(colon(),colon(3,5)) = F1;

        static matrix<double> hm(18,18);
        calc_permutation_hes(&hm[0], &two_frame[0], order_);

        for(size_t i = 0; i < 12; ++i){
            for(size_t j = 0; j < 12; ++j){
                if(jtf::function::add_to_csc(h,ptr,idx,idx_used_[i],idx_used_[j],
                                             w_*hm(idx_in_18_[i],idx_in_18_[j])))
                  return __LINE__;
              }
          }
      }
    return 0;
  }
  virtual int hes_block(const double *x, double *h, double alpha = 1){
    return __LINE__;
  }
  void assemble(){
    idx_in_18_.resize(12,1);
    idx_used_.resize(12,1);
    if(order_ < 3){ // cols
        size_t i = 0;
        for(; i < 9; ++i) {
            idx_in_18_[i] = i;
            idx_used_[i] = 9*idx0_ + i;
          }
        for(; i < 12; ++i) {
            idx_in_18_[i] = i + 3 * order_ ;
            idx_used_[i] = 9*idx1_ + 3*order_+i%3;
          }
      }else{
        size_t i = 3;
        for(; i < 12; ++i) {
            idx_in_18_[i] = i+6;
            idx_used_[i] = 9*idx1_+i-3;
          }
        for(i = 0; i < 3; ++i) {
            idx_in_18_[i] = 3*i+order_;
            idx_used_[i] = 9*idx0_+3*i+order_;
          }
      }
  }
private:
  const size_t tet_num_;
  size_t idx0_,idx1_;
  const size_t order_; // order 0,1,2 (colomns); 3,4,5 (rows)
  const double w_;
  matrix<int32_t> idx_in_18_;
  matrix<int32_t> idx_used_;
  std::vector<int32_t> ptr_vec_;
};


class rotation_matrix_func : public hj_func
{
public:
  rotation_matrix_func(const size_t tet_num,
                       const size_t idx)
    :tet_num_(tet_num), idx_(idx){}

  virtual ~rotation_matrix_func(){}

  virtual size_t dim_of_x(void) const {
    return 9 * tet_num_;
  }
  virtual size_t dim_of_f(void) const {
    return 9;
  }
  virtual int val(const double *x, double *f,
                  hj::function::func_ctx *ctx = 0) const {
    using namespace zjucad::matrix;
    itr_matrix<const double*> x0(9,tet_num_, x);
    itr_matrix<const double*> F(3,3, &x0(0,idx_));
    itr_matrix<double*> f0(3,3,f);

    matrix<double> R = F;
    hj::polar3d p;
    p(R);
    matrix<double> eigv(3,1);
    matrix<double> RR=R;
    eig(RR,eigv);
    if(eigv[0]*eigv[1]*eigv[2]<0){//det < 0
        for(size_t i = 0; i < 3; ++i)
          if(eigv[i] < 0) R(colon(),i) *= -1;
      }

    f0 = F-R;
    f0 *= w_;
    for(size_t i = 0; i < f0.size(); ++i) jtf::math::erase_nan_inf(f0[i]);
    return 0;
  }
  virtual int jac(const double *x, double *val, int32_t *ptr = 0,
                  int32_t *idx = 0, hj::function::func_ctx *ctx = 0) const {

    for(int c = 0; c < 3; ++c) {
        for(int r = 0; r < 3; ++r) {
            const int fi = c*3+r;
            ptr[fi+1] = ptr[fi] + 1;
            idx[ptr[fi]] = 9*idx_ + fi;
            val[ptr[fi]] = w_;
          }
      }

    return 0;
  }
  virtual size_t jac_nnz(void) const {
    return 9;
  }
  static double w_;
private:
  const size_t tet_num_;
  const size_t idx_;
};

double rotation_matrix_func::w_ = 1.0;

jtf_func_ptr
build_surface_normal_align_func(const jtf::tet_mesh & tm, ptree & pt)
{
  const double align_w = pt.get<double>("weight/surface_align.value");
  if(fabs(align_w) < 1e-6) return nullptr;

  vector<jtf_func_ptr> all_funcs;
  const double total_area = std::accumulate(tm.outside_face_area_.begin(), tm.outside_face_area_.end(),0.0);
  const double w = 1.0/total_area;
  for(size_t fi = 0; fi < tm.outside_face_idx_.size(); ++fi){
      const pair<size_t,size_t> & tet_pair = tm.fa_->face2tet_[tm.outside_face_idx_[fi]];
      assert(tet_pair.first == -1 || tet_pair.second == -1);
      const size_t tet_idx = (tet_pair.first==-1?tet_pair.second:tet_pair.first);

      for(size_t di = 0; di < 3; ++di){
          jtf_func_ptr fptr(
                new normal_align_func(
                  tm.tetmesh_.mesh_.size(2),tet_idx,tm.outside_face_normal_(colon(),fi), di, align_w*tm.outside_face_area_[fi]));
          all_funcs.push_back(jtf_func_ptr(new smooth_L1(fptr, w,align_w*tm.outside_face_area_[fi])));
        }
    }
  if(all_funcs.empty()) return nullptr;
  jtf_func_ptr rtn(new jtf::function::sum_function<double,int32_t, jtf::function::SMART_STD>(all_funcs));
  return rtn;
}

jtf_func_ptr
build_inner_frame_smooth_func_euler(const jtf::tet_mesh & tm, ptree &pt)
{
  const double total_vol = std::accumulate(tm.vol_.begin(), tm.vol_.end(), 0.0);

  const double w = 1.0/total_vol;

  vector<jtf_func_ptr> all_funcs;
  for(const auto & one_edge: tm.ortae_.e2t_){
      const vector<size_t> & tet_loop = one_edge.second;
      if(!tm.ortae_.is_inner_edge(tet_loop)) continue;

      assert(tet_loop.front() == tet_loop.back());
      double v = 0;
      for(size_t ti = 0; ti < one_edge.second.size()-1; ++ti) v += tm.vol_[one_edge.second[ti]];
      for(size_t ti = 0; ti < one_edge.second.size()-1; ++ti){
          for(size_t tj = ti+1; tj < one_edge.second.size()-1; ++tj){
              for(size_t di = 0; di < 6; ++di){
                  jtf_func_ptr fptr(
                        new permutation_matrix_euler_func(
                          tm.tetmesh_.mesh_.size(2),
                          tet_loop[ti], tet_loop[tj],di, v*w/6.0));
                  all_funcs.push_back(fptr);
                }
            }
        }
    }

  cerr << "# [info] use Yufei Li strategy of inner smooth." << endl;

  jtf_func_ptr rtn(new jtf::function::sum_function<double,int32_t, jtf::function::SMART_STD>(all_funcs));

  return rtn;
}

jtf_func_ptr
build_inner_frame_smooth_func(const jtf::tet_mesh & tm, ptree &pt, int strategy = 0) // strategy = 0: L1, 1: Yufei Li use edge, 2: Yufei Li use face
{
  vector<jtf_func_ptr> all_funcs;
  const double total_vol = std::accumulate(tm.vol_.begin(), tm.vol_.end(), 0.0);

  const double w = 1.0/total_vol;
  if(strategy == 0){// L1 strategy
      for(size_t fi = 0; fi < tm.fa_->face2tet_.size(); ++fi){
          const pair<size_t,size_t> & tet_pairs = tm.fa_->face2tet_[fi];
          if(tm.fa_->is_outside_face(tet_pairs)) continue;

          const double vw = (fabs(tm.vol_[tet_pairs.first])+fabs(tm.vol_[tet_pairs.second]))/4.0;
          for(size_t di = 0; di < 3; ++di){
              for(size_t dj = 0; dj < 3; ++dj){
                  jtf_func_ptr fptr(
                        new direction_align(tm.tetmesh_.mesh_.size(2),
                                            tet_pairs.first, di,
                                            tet_pairs.second, dj,vw));
                  all_funcs.push_back(jtf_func_ptr(new smooth_L1(fptr,w, vw)));
                }
            }
        }
      cerr << "# [info] use L1 strategy of inner smooth ." << endl;
    }else if(strategy == 1){ // use Yufei Li's strategy for each edge
      for(const auto & one_edge: tm.ortae_.e2t_){
          const vector<size_t> & tet_loop = one_edge.second;
          if(!tm.ortae_.is_inner_edge(tet_loop)) continue;

          assert(tet_loop.front() == tet_loop.back());
          double v = 0;
          for(size_t ti = 0; ti < one_edge.second.size()-1; ++ti) v += tm.vol_[one_edge.second[ti]];
          for(size_t ti = 0; ti < one_edge.second.size()-1; ++ti){
              for(size_t tj = ti+1; tj < one_edge.second.size()-1; ++tj){
                  for(size_t di = 0; di < 6; ++di){
                      jtf_func_ptr fptr(
                            new permutation_matrix_func(
                              tm.tetmesh_.mesh_.size(2),
                              tet_loop[ti], tet_loop[tj],di, v*w/6.0));
                      all_funcs.push_back(fptr);
                    }
                }
            }
        }
      cerr << "# [info] use Yufei Li strategy of inner smooth." << endl;
    }else if(strategy == 2){//use Yufei Li's strategy for each face
      for(size_t fi = 0; fi < tm.fa_->face2tet_.size(); ++fi){
          const pair<size_t,size_t> & tet_pairs = tm.fa_->face2tet_[fi];
          if(tm.fa_->is_outside_face(tet_pairs)) continue;

          for(size_t di = 0; di < 6; ++di){
              jtf_func_ptr fptr(
                    new permutation_matrix_func(
                      tm.tetmesh_.mesh_.size(2),
                      tet_pairs.first, tet_pairs.second,di,
                      (fabs(tm.vol_[tet_pairs.first])+fabs(tm.vol_[tet_pairs.second]))*w/4.0));
              all_funcs.push_back(fptr);
            }
        }
      cerr << "# [info] use Yufei Li strategy of inner smooth." << endl;
    }

  if(all_funcs.empty()) return nullptr;
  jtf_func_ptr rtn(new jtf::function::sum_function<double,int32_t, jtf::function::SMART_STD>(all_funcs));

  return rtn;
}

jtf_func_ptr
build_rotation_matrix_restriction(const jtf::tet_mesh & tm, ptree &pt)
{
  shared_ptr<vector<hj_func_ptr> > all_funcs(new vector<hj_func_ptr>);
  for(size_t ti = 0; ti < tm.tetmesh_.mesh_.size(2); ++ti){
      all_funcs->push_back(
           hj_func_ptr(
              new rotation_matrix_func(tm.tetmesh_.mesh_.size(2), ti)));
    }
  hj_func_ptr funcs(hj::function::new_catenated_function<double,int32_t>(all_funcs));

  jtf_func_ptr jtf_funcs(jtf::function::least_square_warpper(funcs));
  return jtf_funcs;
}

void optimize_frame_field_L1(const jtf::tet_mesh &tm, matrix<double> & frame, ptree &pt)
{
  vector<jtf_func_ptr> all_funcs;
  {// surface normal align
    jtf_func_ptr fptr = build_surface_normal_align_func(tm, pt);
    if(fptr != nullptr){
        all_funcs.push_back(fptr);
        cerr << "# [info] add surface normal align func." << endl;
      }
  }
  {// inner smooth
    jtf_func_ptr fptr = build_inner_frame_smooth_func(tm, pt,0);
    if(fptr != nullptr){
        all_funcs.push_back(fptr);
        cerr << "# [info] add inner frame smooth func." << endl;
      }
  }
  {// restrict frame to be rotation matrix
    all_funcs.push_back(build_rotation_matrix_restriction(tm, pt));
    cerr << "# [info] add rotation matrix restriction func." << endl;
  }

  jtf_func_ptr funs(new jtf::function::sum_function<double,int32_t, jtf::function::SMART_STD>(all_funcs));


  double w = pt.get<double>("weight/rotation.value");
  size_t iter = pt.get<size_t>("input/iter.value");

  for(size_t i = 0; i < iter; ++i){
      w *= 2.0;
      rotation_matrix_func::w_ = w;
      cerr << "# [info] rotation matrix fix weight " << w << endl;
      jtf::optimize(*funs, frame, pt, nullptr, nullptr, nullptr);

      itr_matrix<const double*> frame0(3,3, &frame[0]);
      cerr << frame0 << endl;
      cerr << trans(frame0) * frame0 << endl;
    }
}


void optimize_frame_field_L1_euler(const jtf::tet_mesh & tm,
                                   zjucad::matrix::matrix<double> & zyz,
                                   boost::property_tree::ptree &pt)
{
  // this function exactly follow Yufei Li's strategy
  vector<jtf_func_ptr> all_funcs;
  {// inner smooth
    all_funcs.push_back(build_inner_frame_smooth_func_euler(tm, pt));
    cerr << "# [info] add inner frame smooth func." << endl;
  }
  {// restrict frame to be rotation matrix
    //all_funcs.push_back(build_rotation_matrix_restriction_euler(tm, pt));
    cerr << "# [info] add rotation matrix restriction func." << endl;
  }

  jtf_func_ptr funs(new jtf::function::sum_function<double,int32_t, jtf::function::SMART_STD>(all_funcs));

  jtf::optimize(*funs, zyz, pt, nullptr, nullptr, nullptr);
}
