#ifndef POLYCUBE_PARAM_FUNCTION_TERM_H
#define POLYCUBE_PARAM_FUNCTION_TERM_H

#include "hex_param.h"
#include "../hex_frame_opt/hex_function_term.h"
#include <memory>
#include <numeric>
#include <hjlib/function/function.h>

#include <hjlib/math/blas_lapack.h>
#include <zjucad/optimizer/optimizer.h>

#include <zjucad/matrix/lapack.h>
#include <iostream>
#include <boost/static_assert.hpp>

using namespace std;
using namespace zjucad::matrix;
using namespace hj::function;

class frame_align_normal_function : public sym_frame_opt_function
{
public:
    frame_align_normal_function(
        sym_frame_opt_sys &sys,int32_t idx, const double *n,
        size_t vertex_num,
        size_t frame_num)
        :sym_frame_opt_function(sys),idx_(idx), vertex_num_(vertex_num), frame_num_(frame_num) {
        matrixd zyz(3), Rnz_sh(9, 9);
        rot_n_2_z_by_zyz(n, &zyz[0]);
        calc_rot_cubic_f_sh_mat_(&Rnz_sh[0], &zyz[0]);
        sh_ = trans(Rnz_sh(4, colon()));
    }
    virtual size_t dim_of_x(void) const {
        //return (frame_num_)*3;
        return (vertex_num_ + frame_num_)*3;
    }
    virtual size_t dim_of_f(void) const {
        return 1;
    }
    virtual int val(const double *x, double *f, func_ctx *ctx = 0) const {
        sys_.val(x); // the sys_.val only hanle the sh_.size() of x
        *f = (dot(sys_.sh_(colon(), idx_), sh_)-sqrt_7_);
        return 0;
    }
    virtual int jac(const double *x, double *val, int32_t *ptr = 0, int32_t *idx = 0, func_ctx *ctx = 0) const {
        sys_.jac(x);
        ptr[1] = ptr[0] + 3;
        for(int32_t d = 0; d < 3; ++d) {
            idx[ptr[0]+d] = idx_*3+d;
            val[ptr[0]+d] = dot(sys_.jac_sh_[idx_](colon(), d), sh_);
        }
        return 0;
    }
    virtual size_t jac_nnz(void) const {
        return 3*dim_of_f();
    }

    int32_t idx_;
    matrixd sh_;
    static const double sqrt_7_;
protected:
    const size_t vertex_num_;
    const size_t frame_num_;
};

class gradient_aligend: public function_t<double,int32_t>
{
public:
    gradient_aligend(
        const matrixst &tet,
        const matrixd &node,
        size_t tet_idx,
        const jtf::mesh::face2tet_adjacent & fa,
        size_t tets_num,
        size_t vertex_num):frame_num_(tets_num),tet_idx_(tet_idx), vertex_num_(vertex_num)
    {
        matrixd P = ones<double>(4,4);
        node_idx_.resize(4);
        for(int i = 0; i < 4; ++i) {
            node_idx_[i] = tet[i];
            P(colon(0,2),i) = node(colon(), tet[i]);
        }
        if(inv(P))
            cerr << "inverse fail." << endl;
        M_ = P(colon(), colon(0, 2));
    }

    virtual size_t dim_of_x() const
    {
        return (frame_num_ + vertex_num_) * 3;
    }
    virtual size_t dim_of_f() const
    {
        return 3 * 3;
    }
    virtual int val(const double *x, double *f, func_ctx *ctx=0) const
    {
        itr_matrix<const double*> x0(3,frame_num_ + vertex_num_,x) ;
        itr_matrix<double*> f0(3,3,f);
        matrixd zyz_rot(3,3);
        // TODO: it's strange that can not use x0(colon(),node_idx_) * M_
        zyz_angle_2_rotation_matrix1(&x0(0,tet_idx_),&zyz_rot[0]);
        matrixd x0tmp(3,4);
        for(size_t t = 0; t < 4; ++t)
            x0tmp(colon(),t) = x0(colon(), frame_num_ + node_idx_[t]);
        f0 = x0tmp * M_ - zyz_rot;
        return 0;
    }
    virtual int jac(const double *x, double *val, int_type *ptr = 0, int_type *idx = 0, func_ctx *ctx = 0) const
    {
        itr_matrix<const double*> x0(3,frame_num_ + vertex_num_,x) ;
        const matrixd I = eye<double>(3);
        for(size_t t = 0, i = 0; t < 9; ++t, ++i)
        {
            ptr[i+1] = ptr[i] + 15;
            for(size_t ni = 0; ni < 4; ++ni)
            {
                for(size_t nd = 0; nd < 3; ++nd)
                {
                    idx[ptr[i] + 3 * ni + nd] = frame_num_ * 3 + node_idx_[ni] * 3 + nd;
                    val[ptr[i] + 3 * ni + nd] = I(t% 3,nd) * M_(ni,nd);
                }
            }
            for(size_t nd = 0; nd < 3; ++nd)
            {
                idx[ptr[i] + 4 * 3 + nd] = tet_idx_ * 3 + nd;
                val[ptr[i] + 4 * 3 + nd] = -polycube_jac(&x0(0,tet_idx_),t,nd) ;
            }

        }
//        for(size_t t = 0, i=0 ; t < 3; ++i,++t)
//        {
//            ptr[i+1] = ptr[i] + 4;
//            for(size_t ni = 0; ni < 4; ++ni)
//            {
//                idx[ptr[i] + ni] = node_idx_[ni];
//            }
//        }
        return 0;
    }
    virtual size_t jac_nnz(void) const
    {
        return 15 * dim_of_f();
    }
protected:
    matrixd M_;
    matrixst node_idx_;
    const size_t frame_num_;
    const size_t tet_idx_;
    const size_t vertex_num_;
};

template<size_t a, size_t b, size_t c>
class planar_aligned : public function_t<double, int32_t>
{
    BOOST_STATIC_ASSERT(a < 3 && b < 3 && c < 3 && a != b && a != c && b != c);
public:
    planar_aligned(
        const matrixst &tet,
        const matrixd &node,
        const jtf::mesh::face2tet_adjacent & fa,
        const double *n,
        size_t tet_idx,
        size_t frame_num,
        size_t vertex_num):frame_num_(frame_num),vertex_num_(vertex_num),
        normal_(n),tet_idx_(tet_idx)
    {
        matrixd P = ones<double>(4,4);
        node_idx_.resize(4);
        for(int i = 0; i < 4; ++i) {
            node_idx_[i] = tet[i];
            P(colon(0,2),i) = node(colon(), tet[i]);
        }
        if(inv(P))
            cerr << "inverse fail." << endl;
        M_ = P(colon(), colon(0, 2));
    }

    virtual size_t dim_of_x() const
    {
        return (frame_num_ + vertex_num_) * 3;
    }

    virtual size_t dim_of_f() const
    {
        return 3 * 1;
    }

    virtual int val(const double *x, double *f, func_ctx *ctx=0) const
    {
        itr_matrix<const double*> n0(3,1,normal_);
        itr_matrix<const double*> x0(3,frame_num_ + vertex_num_,x);
        itr_matrix<double *> f0(3,1,f);
        matrixd rot_v(3,3);

        zyz_angle_2_rotation_matrix1(&x0(0,tet_idx_),&rot_v[0]);
        const matrixd &r0 = rot_v(colon(),a);
        const matrixd &r1 = rot_v(colon(),b);
        const matrixd &r2 = rot_v(colon(),c);

        //TODO: there is a bug about dot function , it can not handle
        //dot(const matrix<const double>, const matrixd)
        // the type of r will be const double, and can not be changed.
        const double nrIk = dot(r0, n0);
        const double nrIi = dot(r1, n0);
        const double nrIj = dot(r2, n0);
        const double weight = exp(0.01*nrIk * nrIk/(nrIi * nrIi + nrIj * nrIj)) - 1;
        cerr << "a,b,c" << a << "," << b << "," << c
             << " weight = " << weight << endl;
        matrixd x0tmp(3,4);
        for(size_t t = 0; t < 4; ++t)
            x0tmp(colon(),t) = x0(colon(), frame_num_ + node_idx_[t]);
        f0 = x0tmp * M_;
        f0 *= weight;
        return 0;
    }

    virtual int jac(const value_type *x, value_type *val, int_type *ptr = 0, int_type *idx = 0, func_ctx *ctx = 0) const
    {
        itr_matrix<const double*> x0(3,frame_num_ + vertex_num_,x);
        itr_matrix<const double*> n0(3,1,normal_);
        matrixd rot_v(3,3);

        zyz_angle_2_rotation_matrix1(&x0(0,tet_idx_),&rot_v[0]);
        const matrixd &r0 = rot_v(colon(),a);
        const matrixd &r1 = rot_v(colon(),b);
        const matrixd &r2 = rot_v(colon(),c);

        const double nrIk = dot(r0, n0);
        const double nrIi = dot(r1, n0);
        const double nrIj = dot(r2, n0);
        const double weight = exp(0.01*nrIk * nrIk/(nrIi * nrIi + nrIj * nrIj)) - 1;

        const matrixd I = eye<double>(3);

        matrixd x0tmp(3,4);
        for(size_t t = 0; t < 4; ++t)
            x0tmp(colon(),t) = x0(colon(), frame_num_ + node_idx_[t]);
        const matrixd Lambda = x0tmp * M_(colon(),a);
        const double E = weight + 1;
        const double D = nrIi * nrIi + nrIj * nrIj;
        const double Dk = nrIk * nrIk;
        matrixd N = eye<double>(3);
        N(0,0) = n0[0] * n0[0];
        N(1,1) = n0[1] * n0[1];
        N(2,2) = n0[2] * n0[2];

        size_t i = 0;
        for(size_t t=0; t < 3; ++t, ++i)
        {
            ptr[i+1] = ptr[i] + 15;
            // assemble the vertex jacobi
            for(size_t k = 0; k < 4; ++k)
            {
                for(size_t p = 0; p < 3; ++p)
                {
                    idx[ptr[i] + 3 * k + p] = frame_num_ * 3 + node_idx_[k] * 3 + p;
                    val[ptr[i] + 3 * k + p] = I(t,p) * M_(k,a) * weight;
                }
            }
            for(size_t p = 0; p < 3; ++p)
            {
                idx[ptr[i] + 3 * 4 + p] = tet_idx_ * 3 + p;
                val[ptr[i] + 3 * 4 + p] = ((2 * E * dot(N * r0, Lambda)/D) * polycube_jac(&x0(0,tet_idx_),t + 3 * a,p)
                                           -2 * E * dot(N * r1, Lambda) * Dk/(D*D) * polycube_jac(&x0(0,tet_idx_),t+ 3 * b,p)
                                           -2 * E * dot(N * r2, Lambda) * Dk/(D*D) * polycube_jac(&x0(0,tet_idx_),t+ 3 * c,p)
                                           );
            }
        }

        return 0;
    }

    virtual size_t jac_nnz(void) const
    {
        return 15 * dim_of_f();
    }

protected:
    matrixd M_;
    matrixst node_idx_;
    const size_t tet_idx_;
    const size_t frame_num_;
    const size_t vertex_num_;
    const double * normal_;

};

//const double align_function::sqrt_7_ = sqrt(7.0);
#endif // POLYCUBE_PARAM_FUNCTION_TERM_H
