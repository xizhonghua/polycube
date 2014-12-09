#ifndef PERMUTATION_MATRIX_H
#define PERMUTATION_MATRIX_H

/// @param two_frame 3 *6; n 1*1
extern "C" void calc_permutation_0_(double *n, double * two_frame);

/// @param jac 18*1, two_frame 3 * 6
extern "C" void calc_permutation_0_jac_(double *jac, double * two_frame);

/// @param hes 18*18, two_frame 3*6
extern "C" void calc_permutation_0_hes_(double *hes, double * two_frame);





/// @param two_frame 3 *6; n 1*1
extern "C" void calc_permutation_1_(double *n, double * two_frame);

/// @param jac 18*1, two_frame 3 * 6
extern "C" void calc_permutation_1_jac_(double *jac, double * two_frame);

/// @param hes 18*18, two_frame 3*6
extern "C" void calc_permutation_1_hes_(double *hes, double * two_frame);





/// @param two_frame 3 *6; n 1*1
extern "C" void calc_permutation_2_(double *n, double * two_frame);

/// @param jac 18*1, two_frame 3 * 6
extern "C" void calc_permutation_2_jac_(double *jac, double * two_frame);

/// @param hes 18*18, two_frame 3*6
extern "C" void calc_permutation_2_hes_(double *hes, double * two_frame);





/// @param two_frame 3 *6; n 1*1
extern "C" void calc_permutation_3_(double *n, double * two_frame);

/// @param jac 18*1, two_frame 3 * 6
extern "C" void calc_permutation_3_jac_(double *jac, double * two_frame);

/// @param hes 18*18, two_frame 3*6
extern "C" void calc_permutation_3_hes_(double *hes, double * two_frame);





/// @param two_frame 3 *6; n 1*1
extern "C" void calc_permutation_4_(double *n, double * two_frame);

/// @param jac 18*1, two_frame 3 * 6
extern "C" void calc_permutation_4_jac_(double *jac, double * two_frame);

/// @param hes 18*18, two_frame 3*6
extern "C" void calc_permutation_4_hes_(double *hes, double * two_frame);




/// @param two_frame 3 *6; n 1*1
extern "C" void calc_permutation_5_(double *n, double * two_frame);

/// @param jac 18*1, two_frame 3 * 6
extern "C" void calc_permutation_5_jac_(double *jac, double * two_frame);

/// @param hes 18*18, two_frame 3*6
extern "C" void calc_permutation_5_hes_(double *hes, double * two_frame);


inline void calc_permutation(double *n, double *two_frame, size_t order)
{
  switch(order){
    case 0: calc_permutation_0_(n,two_frame); break;
    case 1: calc_permutation_1_(n,two_frame); break;
    case 2: calc_permutation_2_(n,two_frame); break;
    case 3: calc_permutation_3_(n,two_frame); break;
    case 4: calc_permutation_4_(n,two_frame); break;
    case 5: calc_permutation_5_(n,two_frame); break;
    default:break;
    }
}

inline void calc_permutation_jac(double *jac, double *two_frame, size_t order)
{
  switch(order){
    case 0: calc_permutation_0_jac_(jac,two_frame); break;
    case 1: calc_permutation_1_jac_(jac,two_frame); break;
    case 2: calc_permutation_2_jac_(jac,two_frame); break;
    case 3: calc_permutation_3_jac_(jac,two_frame); break;
    case 4: calc_permutation_4_jac_(jac,two_frame); break;
    case 5: calc_permutation_5_jac_(jac,two_frame); break;
    default:break;
    }
}

inline void calc_permutation_hes(double *hes, double *two_frame, size_t order)
{
  switch(order){
    case 0: calc_permutation_0_hes_(hes,two_frame); break;
    case 1: calc_permutation_1_hes_(hes,two_frame); break;
    case 2: calc_permutation_2_hes_(hes,two_frame); break;
    case 3: calc_permutation_3_hes_(hes,two_frame); break;
    case 4: calc_permutation_4_hes_(hes,two_frame); break;
    case 5: calc_permutation_5_hes_(hes,two_frame); break;
    default:break;
    }
}
#endif // PERMUTATION_MATRIX_H
