#ifndef HEX_FRAME_OPT_ROT_H
#define HEX_FRAME_OPT_ROT_H

// @param permu 1*1, two_orth_frame 3*3
extern "C" void rot_matrix_(double *diff, const double *rot);

// @param permu_jac 6*1, two_orth_frame 3*2
extern "C" void rot_matrix_jac_(double *jac, const double *rot);

// @param permu_hes 6*6, two_orth_frame 3*2
extern "C" void rot_matrix_hes_(double *hes, const double *rot);

// @param permu 1*1, two_orth_frame 3*3
extern "C" void rot_diff_(double *diff, const double *rot, const double  *p0, const double  *p1);

// @param permu_jac 6*1, two_orth_frame 3*2
extern "C" void rot_diff_jac_(double *jac, const double  *rot, const double  *p0, const double  *p1);

// @param permu_hes 6*6, two_orth_frame 3*2
extern "C" void rot_diff_hes_(double *hes, const double  *rot, const double  *p0, const double  *p1);

#endif // ROT_H
