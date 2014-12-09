#ifndef HEX_FRAME_SMOOTH_L1_H
#define HEX_FRAME_SMOOTH_L1_H

// @param permu 1*1, two_orth_frame 3*3
extern "C" void smooth_l1_(double *diff, const double *pp, const double *a);

extern "C" void smooth_l1_jac_(double *jac, const double *pp, const double *a);

extern "C" void smooth_l1_hes_(double *hes, const double *pp, const double *a);

#endif // SMOOTH_L1_H
