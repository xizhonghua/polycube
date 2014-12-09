#ifndef YUFEI_PERMU_H
#define YUFEI_PERMU_H

// @param permu 1*1, two_orth_frame 3*2
extern "C" void calc_yufei_permu_(double *permu, double *two_orth_frame);

// @param permu_jac 6*1, two_orth_frame 3*2
extern "C" void calc_yufei_permu_jac_(double *permu_jac, double *two_orth_frame);

// @param permu_hes 6*6, two_orth_frame 3*2
extern "C" void calc_yufei_permu_hes_(double *permu_hes, double *two_orth_frame);
#endif // YUFEI_PERMU_H
