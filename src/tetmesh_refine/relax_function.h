#ifndef RELAX_FUNCTION_H
#define RELAX_FUNCTION_H
#include <algorithm>
#include <map>
#include <hjlib/function/function.h>
#include <hjlib/function/func_aux.h>
#include <zjucad/matrix/matrix.h>
#include <zjucad/matrix/itr_matrix.h>
#include <hjlib/math/blas_lapack.h>
#include <hjlib/math/polar.h>
#include <zjucad/matrix/lapack.h>

class point_fix_function : public hj::function::function_t<double,int32_t>
{
public:
  point_fix_function(const size_t node_num,
                     const size_t node_idx,
                     const matrixd &node)
    :node_num_(node_num),node_idx_(node_idx),node_(node){}
  virtual size_t dim_of_x()const {return 3 * node_num_;}
  virtual size_t dim_of_f()const {return 3 * 1;}
  virtual int val(const double *x, double *f, hj::function::func_ctx *ctx = 0) const
  {
    zjucad::matrix::itr_matrix<const double*> x0(3,node_num_,x);
    zjucad::matrix::itr_matrix<double*> f0(3,1,f);
    f0 = x0(zjucad::matrix::colon(),node_idx_) - node_;
    return 0;
  }
  virtual int jac(const double *x, double *val, int32_t *ptr, int32_t *idx, hj::function::func_ctx *ctx=0) const
  {
    for(size_t t = 0; t < 3; ++t){
      ptr[t+1] = ptr[t] + 1;
      idx[ptr[t] + 0] = node_idx_ * 3 + t;
      val[ptr[t] + 0] = 1;
    }
    return 0;
  }
  virtual size_t jac_nnz(void)const{ return 3;}
private:
  const size_t node_num_;
  const size_t node_idx_;
  const matrixd node_;
};

class arap_func : public hj::function::function_t<double,int32_t>
{
public:
  arap_func(const matrixst & one_tet,
            const matrixd & node,
            const size_t &idx)
    :idx_(idx), node_num(node.size(2)),one_tet_(one_tet) {
    if(E_.size() == 0) {
      E_.resize(4, 3);
      double val[] = {
        1, 0, 0, -1,
        0, 1, 0, -1,
        0, 0, 1, -1
      };
      std::copy(val, val+E_.size(), E_.begin());
    }

    matrixd A =
        node(zjucad::matrix::colon(),one_tet)*E_;
    if(inv(A))
      std::cerr << "degenerated tet: " << idx_ << std::endl;
    G_ = E_*A;
  }
  virtual ~arap_func(){}
  virtual size_t dim_of_x(void) const {
    return 3 * node_num;
  }
  virtual size_t dim_of_f(void) const {
    return 9;
  }
  virtual int val(const double *x, double *f, hj::function::func_ctx *ctx = 0) const {
    //std::cerr << "# [info] node size " << node_.size(2) << std::endl;
    const zjucad::matrix::itr_matrix<const double *> px(3, node_num, x);
    V_ = px(zjucad::matrix::colon(), one_tet_);
    M_ = V_*G_;
    R_ = M_;
    hj::polar3d p;
    if(p(R_, 2) < 0)
      std::cerr << "# polar fail." << M_ << V_ << G_ << std::endl;
    zjucad::matrix::itr_matrix<double *> pf(3, 3, f);
    pf = M_-R_;
    return 0;
  }
  // NOTICE: only a proximation of the real jac, which ignore R_
  virtual int jac(const double *x, double *val, int32_t *ptr = 0, int32_t *idx = 0,
                  hj::function::func_ctx *ctx = 0) const {
    for(int c = 0, fi = 0; c < 3; ++c) {
      for(int r = 0; r < 3; ++r, ++fi) { // for each element in M_-R_, i.e. M_=V_*G_
        ptr[fi+1] = ptr[fi]+4;
        for(int i = 0; i < 4; ++i) {
          idx[ptr[fi]+i] = one_tet_[i]*3+r;
          val[ptr[fi]+i] = G_(i, c);
        }
      }
    }
    return 0;
  }
  virtual size_t jac_nnz(void) const {
    return 4*dim_of_f();
  }
protected:
  size_t idx_;
  matrixd G_;
  const matrixst  one_tet_;
  //const matrixd & node_;
  const size_t node_num;
  mutable matrixd M_, R_, V_;
  /*  static */matrixd E_;
};

class point_to_sphere_function : public hj::function::function_t<double, int32_t>
{
public:
  point_to_sphere_function(const size_t &node_num, const size_t &node_idx,
                           const double &radius,
                           const matrixd & center_node)
    :node_num_(node_num),node_idx_(node_idx),radius_(radius),
      center_node_(center_node){}
  virtual size_t dim_of_x()const {return 3 * node_num_;}
  virtual size_t dim_of_f()const {return 1;}

  virtual int val(const double *x, double *f, hj::function::func_ctx *ctx = 0) const
  {
    zjucad::matrix::itr_matrix<const double*> x0(3,node_num_,x);
    const double  X = x0(0,node_idx_) - center_node_[0];
    const double  Y = x0(1,node_idx_) - center_node_[1];
    const double  Z = x0(2,node_idx_) - center_node_[2];

    *f = X * X + Y * Y + Z * Z - radius_*radius_;
    return 0;
  }

  virtual int jac(const double *x, double *val, int32_t *ptr, int32_t *idx,
                  hj::function::func_ctx *ctx = 0) const
  {
    zjucad::matrix::itr_matrix<const double*> x0(3,node_num_,x);

    ptr[1] = ptr[0] + 3;
    for(size_t t = 0; t < 3; ++t){
      idx[ptr[0] + t] = 3 * node_idx_ + t;
      val[ptr[0] + t] = 2 * (x0(t, node_idx_) - center_node_[t]);
    }
    return 0;
  }
  virtual size_t jac_nnz(void)const{ return 3;}
private:
  const size_t node_num_;
  const size_t node_idx_;
  const double radius_;
  const matrixd center_node_;
};

//! @brief WARNING: not corret!!!!
class deformation_dest_func : public hj::function::function_t<double,int32_t>
{
  typedef zjucad::matrix::matrix<size_t> matrixi;
  typedef zjucad::matrix::matrix<double> matrixd;
public:
  deformation_dest_func(const matrixi & one_tet,
                        const matrixd &node,
                        const size_t tet_idx):one_tet_(one_tet),tet_idx_(tet_idx)
  {
    E_p.resize(3,3);
    assert(one_tet.size() == 4);
    for(size_t t = 1; t < one_tet.size(); ++t){
      E_p(zjucad::matrix::colon(),t-1) =
          node(zjucad::matrix::colon(),one_tet[t])
          - node(zjucad::matrix::colon(),one_tet[0]);
    }
    if(zjucad::matrix::inv(E_p)){
      std::cerr << "# [error] inverse fail." << std::endl;
    }
    node_num_ = node.size(2);
  }
  virtual size_t dim_of_x() const { return 3 * node_num_; }
  virtual size_t dim_of_f() const { return 3 * 3;}
  virtual int val(const double *x, double *f, hj::function::func_ctx *ctx = 0) const
  {
    zjucad::matrix::itr_matrix<const double *> x0(3,node_num_,x);
    zjucad::matrix::itr_matrix<double *> f0(3,3,f);
    matrixd E = zjucad::matrix::zeros<double>(3,3);
    for(size_t t = 1; t < one_tet_.size(); ++t){
      E(zjucad::matrix::colon(),t-1) = x0(zjucad::matrix::colon(),one_tet_[t])
                                       - x0(zjucad::matrix::colon(),one_tet_[0]);
    }

    matrixd G_; // deformation gradient
    G_ = E_p * E;

    matrixd U,S,VT,E_bkp;
    // static std::map<size_t,matrixd> map_E_bkp;
    //typedef std::map<size_t,matrixd>::iterator msmit;

    //msmit mm = map_E_bkp.find(tet_idx_);
    matrixd A_ = G_;

    //    if(mm == map_E_bkp.end()){
    //      map_E_bkp[tet_idx_] = E;
    //      E_bkp = E;
    //      A_ = zjucad::matrix::eye<double>(3);
    //    }else{
    //      E_bkp = mm->second;
    //      mm->second = E;
    //      A_ = E_p * E_bkp;
    //    }
    //    std::cerr << "# [info] show me static bool address " << &is_ebkp_init << std::endl;
    zjucad::matrix::svd(A_,U,S,VT);
#if 1
    matrixd U_ = U, VT_ = VT;
    const double det_U = zjucad::matrix::det(U_);
    const double det_VT = zjucad::matrix::det(VT_);

    assert(fabs(det_U - 1) < 1e-8 || fabs(det_U + 1) < 1e-8);
    assert(fabs(det_VT - 1) < 1e-8 || fabs(det_VT + 1) < 1e-8);

    std::vector<std::pair<double,size_t> > S_value_idx;
    for(size_t t = 0;t < S.size(2); ++t)
    {S_value_idx.push_back(std::make_pair(S(t,t),t));}
    std::sort(S_value_idx.begin(),S_value_idx.end());

    if(fabs(det_U + 1) < 1e-8)
    {
      const std::pair<double,size_t> &minimal_s = S_value_idx.front();
      U(minimal_s.second,zjucad::matrix::colon()) *= -1;
      S(minimal_s.second,minimal_s.second) *= -1;
    }
    if(fabs(det_VT + 1) < 1e-8)
    {
      const std::pair<double,size_t> &minimal_s = S_value_idx.front();
      VT(zjucad::matrix::colon(),minimal_s.second) *= -1;
      S(minimal_s.second,minimal_s.second) *= -1;
    }
#endif
    //    if(det_G  < 0){
    //      const std::pair<double,size_t> &minimal_s = S_value_idx.front();
    //      S(minimal_s.second,minimal_s.second) *= -1;
    //    }
    //    matrixd newG  = U * S * VT;

    // TODO: I do not know how to format G - R(U*VT),
    // is it necessary to multipl -1 in G
    f0 = G_ - U * VT;
    //    std::cerr << "# [info] det U/VT " << det_U << " " << det_VT << std::endl;
    //std::cerr << "# [info] S " << S << std::endl;
    //std::cerr << "# [info] G - UVT" << f0 << std::endl;
    return 0;
  }

  virtual int jac(const double *x, double *val, int32_t *ptr, int32_t *idx,
                  hj::function::func_ctx *ctx = 0) const
  {
    zjucad::matrix::itr_matrix<const double *> x0(3, node_num_,x);
    //    matrixd e2x(3,4); // store the edge to coordinate representation
    //    e2x(zjucad::matrix::colon(),0) = zjucad::matrix::ones<double>(3,1) * -1;
    //    e2x(zjucad::matrix::colon(),zjucad::matrix::colon(1,3))
    //        = zjucad::matrix::eye<double>(3);
    //    const matrixd coeff = E_p * e2x;
    //    assert(coeff.size(1) == 3 && coeff.size(2) == 4);

    //    for(size_t t = 0; t < 9; ++t){
    //      ptr[t + 1] = ptr[t] + 4;
    //      for(size_t i = 0; i < 4; ++i)
    //      {
    //        idx[ptr[t] + i] = one_tet_[i] * 3 + t%3; // t/3 to determine we use x/y/z
    //        val[ptr[t] + i] = coeff(t%3,i);
    //      }
    //    }

    for(size_t t = 0; t < 9; ++t)
    {
      ptr[t + 1] = ptr[t] + 6;
      for(size_t i = 0; i < 3; ++i){
        idx[ptr[t] + i] = one_tet_[t/3 + 1] * 3 + i;
        val[ptr[t] + i] = E_p(i,t%3);
      }

      for(size_t i = 3; i < 6; ++i){
        idx[ptr[t] + i] = one_tet_[0] * 3 + i - 3;
        val[ptr[t] + i] = -1 * E_p(i-3,t%3);
      }
    }
    return 0;
  }

  virtual size_t jac_nnz()const{
    return 9 * 6;
  }
private:
  matrixd E_p; // G = E_p * E ,
  size_t node_num_;
  const matrixi one_tet_;
  const size_t tet_idx_;
};

class equal_function: public hj::function::function_t<double,int32_t>
{
public:
  equal_function(const size_t & node_num,
                 const size_t & point_a,
                 const size_t & point_b)
    :point_a_(point_a),point_b_(point_b),node_num_(node_num){}

  virtual ~equal_function(){}

  virtual size_t dim_of_x()const {return 3 * node_num_;}
  virtual size_t dim_of_f()const {return 3 * 1;}

  virtual int val(const double *x, double *f, hj::function::func_ctx *ctx = 0) const
  {
    zjucad::matrix::itr_matrix<const double*> x0(3,node_num_,x);
    zjucad::matrix::itr_matrix<double*> f0(3,1,f);
    f0 = x0(zjucad::matrix::colon(),point_a_)
         - x0(zjucad::matrix::colon(),point_b_);
    return 0;
  }
  virtual int jac(const double *x, double *val, int32_t *ptr, int32_t *idx, hj::function::func_ctx *ctx=0) const
  {
    for(size_t t = 0; t < 3; ++t){
      ptr[t+1] = ptr[t] + 2;
      idx[ptr[t] + 0] = point_a_ * 3 + t;
      val[ptr[t] + 0] = 1;

      idx[ptr[t] + 1] = point_b_ * 3 + t;
      val[ptr[t] + 1] = -1;
    }
    return 0;
  }
  virtual size_t jac_nnz(void)const{ return 2*3;}

  const size_t point_a_,point_b_;
  const size_t node_num_;
};



#endif // RELAX_FUNCTION_H
