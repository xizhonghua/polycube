#ifndef MESH_FUNC_TRI_AREA_H_
#define MESH_FUNC_TRI_AREA_H_

//! @param tri: 3x3

//! @param area: 1x1
extern "C" void calc_tri_area_(double *tri_area, double *tri);

//! @param tri_area_jac: 1x9
extern "C" void calc_tri_area_jac_(double *tri_area_jac, double *tri);

#endif
