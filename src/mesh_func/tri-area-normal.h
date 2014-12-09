#ifndef MESH_FUNC_TRI_AREA_NORMAL_H_
#define MESH_FUNC_TRI_AREA_NORMAL_H_

//! @param tri: 3x3

//! @param tri_area_normal: 3x1
extern "C" void calc_tri_area_normal_(double *tri_area_normal, double *tri);

//! @param tri_area_normal: 3x9
extern "C" void calc_tri_area_normal_jac_(double *tri_area_normal_jac, double *tri);

//! @param tri_area_normal: 27x9
extern "C" void calc_tri_area_normal_hes_(double *tri_area_normal_hes, double *tri);

//! @param tri_area_normal: 1*1
extern "C" void tri_area_normal_diff_(double *diff, const double *tri, const double *target_normal);

//! @param tri_area_normal: 1*9
extern "C" void tri_area_normal_diff_jac_(double *jac, const double *tri, const double *target_normal);

//! @param tri_area_normal: 9*9
extern "C" void tri_area_normal_diff_hes_(double *hes, const double *tri, const double *target_normal);
#endif
