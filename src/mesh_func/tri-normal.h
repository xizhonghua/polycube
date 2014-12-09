#ifndef MESH_FUNC_TRI_NORMAL_H
#define MESH_FUNC_TRI_NORMAL_H

//! @param edge_dir: 3x1
extern "C" void calc_edge_dir_(double *edge_dir, double *edge);

//! @param edge_dir_jac: 3 * 6
extern "C" void calc_edge_dir_jac_(double * edge_dir_jac, double *edge);

//! @param edge_dir_hes: 18 * 6
extern "C" void calc_edge_dir_hes_(double *edge_dir_hes, double *edge);

//! @param edge_dir: 3x1
extern "C" void calc_edge_len_dir_(double *edge_len_dir, double *edge);

//! @param edge_dir_jac: 3 * 6
extern "C" void calc_edge_len_dir_jac_(double * edge_len_dir_jac, double *edge);

//! @param edge_dir_hes: 18 * 6
extern "C" void calc_edge_len_dir_hes_(double *edge_len_dir_hes, double *edge);

//! @param edge_dir: 3x1
extern "C" void calc_edge_len_(double *edge_len_, double *edge);

//! @param edge_dir_jac: 3 * 6
extern "C" void calc_edge_len_jac_(double * edge_len_jac, double *edge);

//! @param edge_dir_hes: 18 * 6
extern "C" void calc_edge_len_hes_(double *edge_len_hes, double *edge);


//! @param tri_normal: 3x1
extern "C" void calc_tri_normal_(double *tri_normal, double *tri);

//! @param tri_normal: 3x9
extern "C" void calc_tri_normal_jac_(double *tri_normal_jac, double *tri);

//! @param tri_normal: 27x9
extern "C" void calc_tri_normal_hes_(double *tri_normal_hes, double *tri);

extern "C" void calc_tri_normal_dot_(double *tri_normal, double *tp);

//! @param tri_normal: 3x9
extern "C" void calc_tri_normal_dot_jac_(double *tri_normal_jac, double *tp);


//! @ sub 3x1
extern "C" void calc_adj_norm_sub_(double *sub, double *tp);

//! @param jac 3x12
extern "C" void calc_adj_norm_sub_jac_(double *sub_jac, double *tp);


#endif // MESH_FUNC_TRI_NORMAL_H
