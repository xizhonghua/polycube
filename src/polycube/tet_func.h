#ifndef HJ_POLYCUBE_TET_FUNC_H_
#define HJ_POLYCUBE_TET_FUNC_H_

#include <hjlib/function/function.h>
#include <zjucad/matrix/matrix.h>

//! @param: tet 3x4, grad_op 4x3, tet*grad_op = grad
int calc_tet_def_grad_op(const double *tet, double *grad_op);

//! @param : tet 3 * 4, dis 3 * 3
int calc_tet_arap_distortion(const double *orig_tet,
                             const double *new_tet,
                             double *dis);

//! @brief: a function of |node.size()| -> |tet.size(2)|*9
hj::function::function_t<double, int32_t> *
build_tetmesh_grad_func(const zjucad::matrix::matrix<double> &node,
                        const zjucad::matrix::matrix<size_t> &tet);

//! @brief: constant jac approximation
hj::function::function_t<double, int32_t> *
build_tetmesh_arap_func(const zjucad::matrix::matrix<double> &node,
                        const zjucad::matrix::matrix<size_t> &tet,
                        const double arap_w = 1);

//! @brief: constant jac approximation
hj::function::function_t<double, int32_t> *
build_tetmesh_arap_func(const zjucad::matrix::matrix<double> &node,
                        const zjucad::matrix::matrix<double> &compress_node,
                        const std::vector<std::vector<std::pair<size_t,double> > > & var_map,
                        const zjucad::matrix::matrix<size_t> &tet,
                        const double arap_w = 1);

hj::function::function_t<double, int32_t> *
build_tetmesh_arap_func(const zjucad::matrix::matrix<double> &node,
                        const zjucad::matrix::matrix<double> &compress_node,
                        const std::vector<std::vector<std::pair<size_t,double> > > & var_map,
                        const zjucad::matrix::matrix<size_t> &tet,
                        std::vector<boost::shared_ptr<hj::function::function_t<double,int32_t> > > &func_temp,
                        const double arap_w = 1);


#endif
