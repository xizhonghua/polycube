#ifndef MESH_FUNC_SO3_H_
#define MESH_FUNC_SO3_H_

//! @param M: 3x3

//! @param err: 3x3
extern "C" void calc_so3_(double *err, double *M);

//! @param err: 9*9
extern "C" void calc_so3_jac_(double *err_jac, double *M);

#endif
