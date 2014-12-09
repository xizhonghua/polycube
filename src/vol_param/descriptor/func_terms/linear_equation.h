#ifndef LINEAR_EQUATION_H
#define LINEAR_EQUATION_H

#include "../def.h"
#include <jtflib/function/func_aux.h>
#include <hjlib/sparse/sparse.h>
#include <zjucad/matrix/itr_matrix.h>
#include <hjlib/sparse/fast_AAT.h>
#include <stdexcept>
#include "../../common/util.h"
class linear_equation_hj_func : public hj::function::function_t<double,int32_t>
{
public:
    linear_equation_hj_func(
            const size_t & node_number,
            const std::vector<size_t> & variable_idx,
            const std::vector<double> & variable_coeff,
            const double eqn_val)
        :idx_(variable_idx), coeff_(variable_coeff), val_(eqn_val),
          node_num_(node_number){
        assert(variable_coeff.size() == variable_idx.size());
        for(const auto & idx : variable_idx)
            if(idx > node_num_)
                throw std::invalid_argument("input equation idx overrange.");
    }
    virtual ~linear_equation_hj_func(){}

public:
    virtual size_t dim_of_x(void) const {
        return node_num_;
    }
    virtual size_t dim_of_f(void) const{
        return 1;
    }
    virtual int val(const double *x, double *f,
                    hj::function::func_ctx *ctx = 0) const
    {
        *f = -1*val_;
        for(size_t i = 0; i < idx_.size(); ++i)
            *f += x[idx_[i]] * coeff_[i];
        *f *= weight_;
        return 0;
    }
    virtual int jac(const double *x, double *val, int32_t *ptr = 0,
                    int32_t *idx = 0, hj::function::func_ctx *ctx = 0 )const
    {
        ptr[1] = ptr[0] + coeff_.size();
        for(size_t i = 0; i < coeff_.size(); ++i){
            idx[ptr[0]+i] = idx_[i];
            val[ptr[0]+i] = coeff_[i] * weight_;
        }
        return 0;
    }
    virtual size_t jac_nnz(void) const
    {
        return coeff_.size();
    }
public:
    static double weight_;
private:
    const std::vector<size_t> idx_;
    const std::vector<double> coeff_;
    const double val_;
    const size_t node_num_;
};

template <typename T>
void set_linear_equation_hj_func_weight(const T &v){
    linear_equation_hj_func::weight_ = v;
}

//class linear_equation_jtf_func : public jtf_func
//{
//public:
//  linear_equation_jtf_func(
//      const size_t & node_number,
//      const std::vector<size_t> & variable_idx,
//      const std::vector<double> & variable_coeff,
//      const double eqn_val,
//      const hj::sparse::csc<double,int32_t> * node_mapping = 0)
//    :idx_(variable_idx), coeff_(variable_coeff), val_(eqn_val),
//      node_num_(node_number), node_mapping_(node_mapping) {
//    assert(variable_coeff.size() == variable_idx.size());
//  }
//  virtual ~linear_equation_jtf_func(){}

//public:
//  virtual size_t dim(void) const {
//    if(!node_mapping_)
//      return node_num_;
//    else
//      return zjucad::matrix::max(*node_mapping_)+1;
//  }
//  virtual int val(const double *x, double &f) const
//  {
//    f += -1*val_;
//    for(size_t i = 0; i < idx_.size(); ++i){
//        f += (node_mapping_?x[(*node_mapping_)[idx_[i]]]*coeff_[i]: x[idx_[i]] * coeff_[i]);
//      }
//    return 0;
//  }
//  virtual int gra(const double *x, double *g)
//  {
//    for(size_t i = 0; i < coeff_.size(); ++i){
//        g[ node_mapping_?(*node_mapping_)[idx_[i]]:idx_[i]] += coeff_[i];
//      }
//    return 0;
//  }
//  virtual int gra(const double *x, size_t &nnz, double * g, int32_t *idx)
//  {

//    if(g == 0 || idx == 0) {
//        nnz = coeff_.size();
//        return 0;
//      }

//    if(g != 0 && idx != 0){
//        for(size_t i = 0; i < idx_.size(); ++i){
//            idx[i] = node_mapping_?(*node_mapping_)[idx_[i]]:idx_[i];
//            g[i] = coeff_[i];
//          }
//      }
//    return 0;
//  }
//  virtual int hes(const double *x, size_t &nnz, size_t & format, double * h,
//                  int32_t *ptr, int32_t *idx, double alpha = 1)
//  {
//    if(h == 0 && ptr == 0 && idx == 0) nnz = 0;
//    return 0;
//  }
//  virtual int hes_block(const double *x, double *h, double alpha=1)
//  {
//    return __LINE__;
//  }
//private:
//  const std::vector<size_t> idx_;
//  const std::vector<double> coeff_;
//  const double val_;
//  const size_t node_num_;
//  const hj::sparse::csc<double,int32_t> * node_mapping_;
//};

class point_fix_hj : public hj_func
{
public:
    point_fix_hj(const size_t node_num,
                 const size_t point_idx,
                 const zjucad::matrix::matrix<double> & target,
                 const double weight = 1.0)
        :nm_(node_num), pidx_(point_idx), target_(target), w_(weight){}
    virtual ~point_fix_hj(){}
public:
    virtual size_t dim_of_x(void)const{
        return 3 * nm_;
    }
    virtual size_t dim_of_f()const{
        return 3;
    }
    virtual int val(const double *x, double *f, hj::function::func_ctx *ctx) const{
        using namespace zjucad::matrix;
        itr_matrix<const double*> x0(3, nm_,x);
        itr_matrix<double*> f0(3,1,f);

        f0 = x0(colon(),pidx_) - target_;
        f0 *= w_;
        return 0;
    }
    virtual int jac(const double * x, double *val, int32_t * ptr, int32_t * idx,
                    hj::function::func_ctx *ctx = 0)const{
        for(size_t di = 0; di < 3; ++di){
            ptr[di+1] = ptr[di] + 1;
            idx[ptr[di]] = 3 * pidx_ + di;
            val[ptr[di]] = w_;
        }
        return 0;
    }
    virtual size_t jac_nnz() const{
        return dim_of_f() * 1;
    }
private:
    const size_t nm_;
    const size_t pidx_;
    const zjucad::matrix::matrix<double> target_;
    const double w_;
};

class variant_fix_hj : public hj_func
{
public:
    variant_fix_hj(const size_t v_num,
                   const size_t v_idx,
                   const double value,
                   const hj::sparse::csc<double,int32_t> * node_mapping = 0,
                   const double weight = 1.0)
        :v_num_(v_num), v_idx_(v_idx), value_(value), weight_(weight),
          node_mapping_(node_mapping){}
    virtual ~variant_fix_hj(){}
public:
    virtual size_t dim_of_x(void)const{
        if(!node_mapping_){
            return v_num_;
        }else{
            return node_mapping_->size(1);
        }
    }
    virtual size_t dim_of_f()const{
        return 1;
    }
    virtual int val(const double *x, double *f, hj::function::func_ctx *ctx) const{
        if(!node_mapping_)
            *f = weight_ * (x[v_idx_] - value_);
        else{
            const size_t begin =  node_mapping_->ptr()[v_idx_];
            const size_t end = node_mapping_->ptr()[v_idx_+1];
            assert(begin != end);
            *f = -1 * weight_ * value_;
            for(size_t i = begin; i != end; ++i){
                *f += weight_ * x[node_mapping_->idx()[i]] * node_mapping_->val()[i];
            }
        }
        return 0;
    }
    virtual int jac(const double * x, double *val, int32_t * ptr, int32_t * idx,
                    hj::function::func_ctx *ctx = 0)const{
        if(!node_mapping_){
            ptr[1] = ptr[0] + 1;
            idx[ptr[0]] = v_idx_;
            val[ptr[0]] = weight_;
        }else{
            const size_t begin =  node_mapping_->ptr()[v_idx_];
            const size_t end = node_mapping_->ptr()[v_idx_+1];
            assert(begin != end);
            ptr[1] = ptr[0] + end - begin;
            for(size_t i = begin; i != end; ++i){
                idx[ptr[0] + i - begin] = node_mapping_->idx()[i];
                val[ptr[0] + i - begin] = weight_ * node_mapping_->val()[i];
            }
        }
        return 0;
    }
    virtual size_t jac_nnz() const{
        if(!node_mapping_) return 1;

        const size_t begin =  node_mapping_->ptr()[v_idx_];
        const size_t end = node_mapping_->ptr()[v_idx_+1];
        return end - begin;
    }
private:
    const size_t v_num_;
    const size_t v_idx_;
    const double value_;
    const double weight_;
    const hj::sparse::csc<double,int32_t> * node_mapping_;
};

class variant_fix_hj2 : public hj::function::function_t<double,int32_t>
{
public:
    variant_fix_hj2(const size_t v_num,
                    const size_t v_idx,
                    const double value,
                    const node_mapping * nm = 0,
                    const double weight = 1.0)
        :v_num_(v_num), v_idx_(v_idx), value_(value), weight_(weight),
          node_mapping_(nm){}
    virtual ~variant_fix_hj2(){}

public:
    virtual size_t dim_of_x(void) const {
        if(!node_mapping_)
            return v_num_;
        else
            return node_mapping_->ZT.size(1);
    }
    virtual size_t dim_of_f(void) const{
        return 1;
    }

    virtual int val(const double *x, double *f, hj::function::func_ctx *ctx) const
    {
        if(!node_mapping_)
            *f = weight_ * (x[v_idx_] - value_);
        else{
            const size_t begin =  node_mapping_->ZT.ptr()[v_idx_];
            const size_t end = node_mapping_->ZT.ptr()[v_idx_+1];
            assert(begin != end);
            *f = 0;
            for(size_t i = begin; i != end; ++i){
                *f += weight_ * x[node_mapping_->ZT.idx()[i]] * node_mapping_->ZT.val()[i];
            }
            *f += weight_ * node_mapping_->q[v_idx_];
            *f -= weight_ * value_;
        }
        return 0;
    }
    virtual int jac(const double * x, double *val, int32_t * ptr, int32_t * idx,
                    hj::function::func_ctx *ctx = 0)const
    {
        if(!node_mapping_){
            ptr[1] = ptr[0] + 1;
            idx[ptr[0]] = v_idx_;
            val[ptr[0]] = weight_;
        }else{
            ptr[1] = ptr[0] + node_mapping_->ZT.ptr()[v_idx_+1]-node_mapping_->ZT.ptr()[v_idx_];
            const size_t begin =  node_mapping_->ZT.ptr()[v_idx_];
            const size_t end = node_mapping_->ZT.ptr()[v_idx_+1];
            assert(begin != end);
            for(size_t i = begin; i != end; ++i){
                idx[ptr[0]+i-begin] = node_mapping_->ZT.idx()[i];
                val[ptr[0]+i-begin] = weight_;
            }
        }
        return 0;
    }

    virtual size_t jac_nnz()const{
        if(!node_mapping_) return 1;
        const size_t begin =  node_mapping_->ZT.ptr()[v_idx_];
        const size_t end = node_mapping_->ZT.ptr()[v_idx_+1];
        return end - begin;
    }
private:
    const size_t v_num_;
    const size_t v_idx_;
    const double value_;
    const double weight_;
    const node_mapping * node_mapping_;
};

class variant_fix : public jtf_func
{
public:
    variant_fix(const size_t v_num,
                const size_t v_idx,
                const double value,
                const node_mapping * nm = 0,
                const double weight = 1.0)
        :v_num_(v_num), v_idx_(v_idx), value_(value), weight_(weight),
          node_mapping_(nm){}
    virtual ~variant_fix(){}

public:
    virtual size_t dim(void) const {
        if(!node_mapping_)
            return v_num_;
        else
            return node_mapping_->ZT.size(1);
    }
    virtual int val(const double *x, double &f)
    {
        if(!node_mapping_)
            f += weight_ * (x[v_idx_] - value_);
        else{
            const size_t begin =  node_mapping_->ZT.ptr()[v_idx_];
            const size_t end = node_mapping_->ZT.ptr()[v_idx_+1];
            assert(begin != end);
            for(size_t i = begin; i != end; ++i){
                f += weight_ * x[node_mapping_->ZT.idx()[i]] * node_mapping_->ZT.val()[i];
            }
            f += weight_ * node_mapping_->q[v_idx_];
            f -= weight_ * value_;
        }
        return 0;
    }
    virtual int gra(const double *x, double *g)
    {
        if(!node_mapping_)
            g[v_idx_] += weight_;
        else{
            const size_t begin =  node_mapping_->ZT.ptr()[v_idx_];
            const size_t end = node_mapping_->ZT.ptr()[v_idx_+1];
            assert(begin != end);
            for(size_t i = begin; i != end; ++i){
                g[node_mapping_->ZT.idx()[i]] += weight_;
            }
        }

        return 0;
    }
    virtual int gra(const double *x, size_t &nnz, double * g, int32_t *idx)
    {
        if(!node_mapping_){
            if(g == 0 || idx == 0) {
                nnz = 1;
                return 0;
            }

            if(g != 0 && idx != 0){
                g[0] = weight_;
                idx[0] = v_idx_;
            }
        }else{
            const size_t begin =  node_mapping_->ZT.ptr()[v_idx_];
            const size_t end = node_mapping_->ZT.ptr()[v_idx_+1];
            assert(begin != end);

            if(g == 0 || idx == 0){
                nnz = end - begin;
                return 0;
            }
            if(g != 0 && idx != 0){
                for(size_t i = begin; i != end; ++i){
                    g[i - begin] = node_mapping_->ZT.val()[i];
                    idx[i - begin] = node_mapping_->ZT.idx()[i];
                }
                return 0;
            }
        }
        return 0;
    }
    virtual int hes(const double *x, size_t &nnz, size_t & format, double * h,
                    int32_t *ptr, int32_t *idx, double alpha = 1)
    {
        if(h == 0 && ptr == 0 && idx == 0){
            nnz = 0;
            return 0;
        }

    }
    virtual int hes_block(const double *x, double *h, double alpha=1)
    {
        return __LINE__;
    }
private:
    const size_t v_num_;
    const size_t v_idx_;
    const double value_;
    const double weight_;
    const node_mapping * node_mapping_;
};

class X_c_func : public hj_func
{
public:
    X_c_func(const zjucad::matrix::matrix<double> & b,
             const double w = 1)
        :b_(b), nm_(b.size()), w_(w){}
    virtual size_t dim_of_x() const{return nm_;}
    virtual size_t dim_of_f() const{ return nm_;}
    virtual int val(const double *x, double *f, hj::function::func_ctx *ctx) const{
        using namespace zjucad::matrix;

        itr_matrix<const double *> x0(nm_,1,x);
        itr_matrix<double*> f0(nm_,1,f);
        f0 = x0 - b_;
        f0 *= w_;
        return 0;
    }
    virtual int jac(const double *x, double *val, int32_t *ptr, int32_t *idx,
                    hj::function::func_ctx *ctx) const
    {
        for(size_t i = 0; i < nm_; ++i){
            ptr[i+1] = ptr[i] + 1;
            idx[ptr[i]] = i;
            val[ptr[i]] = w_;
        }
        return 0;
    }
    virtual size_t jac_nnz() const{
        return nm_;
    }
private:
    const zjucad::matrix::matrix<double> &b_;
    const size_t nm_;
    const double w_;
};

///
/// \brief The AX_b_func class, ||Z * u + b - x||^2
class AX_b_func : public hj_func
{
public:
    AX_b_func(const node_mapping &ZT_q,
              const zjucad::matrix::matrix<double> &b)
        :ZT_q_(ZT_q), b_(b){
        assert(ZT_q_.ZT.size(2) == b.size());
        for(size_t i = 0; i < ZT_q_.ZT.ptr().size()-1; ++i){
            if(ZT_q_.ZT.ptr()[i+1] != ZT_q_.ZT.ptr()[i]) f_idx.push_back(i);
        }
        f_temp = zjucad::matrix::zeros<double>(b.size(),1);
    }
    virtual size_t dim_of_x() const{
        return ZT_q_.ZT.size(1);
    }
    virtual size_t dim_of_f() const{
        return f_idx.size();
    }
    virtual int val(const double *x, double *f, hj::function::func_ctx *ctx) const {
        using namespace zjucad::matrix;
        itr_matrix<const double*> x0(dim_of_x(),1, x);
        itr_matrix<double *> f0(dim_of_f(),1, f);
        itr_matrix<const size_t*> f_idx_(f_idx.size(),1, &f_idx[0]);
        f0 = -1*b_(f_idx_) + ZT_q_.q(f_idx_);

        f_temp *= 0;
        hj::sparse::mv(true, ZT_q_.ZT, x0, f_temp);
        f0 += f_temp(f_idx_);
        return 0;
    }
    virtual int jac(const double * x, double *val, int32_t * ptr, int32_t * idx,
                    hj::function::func_ctx *ctx = 0) const {
        using namespace zjucad::matrix;
        for(size_t fi = 0; fi < f_idx.size(); ++fi){
            ptr[fi+1] = ptr[fi] + ZT_q_.ZT.ptr()[f_idx[fi]+1] - ZT_q_.ZT.ptr()[f_idx[fi]];
            for(size_t i = ZT_q_.ZT.ptr()[f_idx[fi]]; i != ZT_q_.ZT.ptr()[f_idx[fi]+1]; ++i){
                idx[ptr[fi] + i - ZT_q_.ZT.ptr()[f_idx[fi]]] = ZT_q_.ZT.idx()[i];
                val[ptr[fi] + i - ZT_q_.ZT.ptr()[f_idx[fi]]] = ZT_q_.ZT.val()[i];
            }
        }

        return 0;
    }
    virtual size_t jac_nnz() const{
        return ZT_q_.ZT.idx().size(); // nnz
    }
private:
    const node_mapping &ZT_q_;
    const zjucad::matrix::matrix<double> &b_;
    std::vector<size_t> f_idx;
    mutable zjucad::matrix::matrix<double> f_temp;
};

//template <typename T1, typename T2>
//void build_linear_equation_funcs(
//    const size_t node_num,
//    const std::vector<std::vector<T1> > & variable_idx,
//    const std::vector<std::vector<T2> > & variable_coeff,
//    std::vector<jtf_func_cons_ptr> & cons,
//    bool use_least_square = false)
//{
//  assert(variable_coeff.size() == variable_idx.size());

//  for(size_t i = 0; i < variable_coeff.size(); ++i){
//      if(!use_least_square)
//        cons.push_back(
//              jtf_func_cons_ptr(new linear_equation_jtf_func(node_num, variable_idx[i],
//                                                             variable_coeff[i],0.0)));
//      else{
//          cons.push_back(jtf_func_cons_ptr(
//                           jtf::function::least_square_warpper(
//                             hj_func_cons_ptr(new linear_equation_hj_func(
//                                                node_num, variable_idx[i], variable_coeff[i], 0.0))
//                             )));
//        }
//    }
//}


//template <typename T1>
//void build_group_equation_funcs(
//    const size_t node_num,
//    const std::vector<std::vector<T1> > & groups,
//    const std::vector<jtf_func_cons_ptr> & cons)
//{
//  std::vector<T1> variable_idx(2);
//  std::vector<double> variable_coeff={1.0,-1.0};
//  for(const auto & g : groups) {
//      for(size_t i = 0; i < g.size()-1; ++i){
//          variable_idx[0] = g[i];
//          variable_idx[1] = g[i+1];
//          cons.push_back(
//                jtf_func_cons_ptr(
//                  new linear_equation_jtf_func(node_num, variable_idx,
//                                               variable_coeff,0.0)));
//        }
//    }
//}

template <typename ITERATION>
void add_zero_index_equation(std::vector<jtf_func_ptr> & obj_vec,
                             ITERATION begin, ITERATION end,
                             const size_t variable_number,
                             const node_mapping * node_mapping_)
{
    assert(node_mapping_);
    for(ITERATION it = begin; it != end; ++it){
        if(node_mapping_){
            if(node_mapping_->is_irrelevant_variable(*it))
                continue;
        }
        obj_vec.push_back(
                    jtf_func_ptr(
                        jtf::function::least_square_warpper(
                            hj_func_ptr(new variant_fix_hj(variable_number, *it, 0, &(node_mapping_->ZT), 1.0)))));
    }
}

#endif // LINEAR_EQUATION_H
