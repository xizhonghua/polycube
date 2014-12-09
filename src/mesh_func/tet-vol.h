#ifndef MESH_FUNC_TETVOL_H
#define MESH_FUNC_TETVOL_H

//! @param tet_node: 3x4
extern "C" void calc_tet_vol_(double * tet_vol, double *tet_node);

extern "C" void calc_tet_vol_jac_(double *jac, double *tet_node);

#endif // TETVOL_H
