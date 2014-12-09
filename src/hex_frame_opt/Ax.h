#ifndef JTF_AX_B_H
#define JTF_AX_B_H

//! @param x:9*2, A:9*9 v=0.5*A*norm(x(colon,1)-x(colon(),2))
extern "C" void ax_(double *v, const double *x, const double *A);

//! @param x:9*2, jac: 9*2, A:9*9
extern "C" void ax_jac_(double *jac, const double *x, const double *A);

//! @param x:9*2, hes:18*18, A:9*9
extern "C" void ax_hes_(double * hes, const double *x, const double *A);
#endif // AX_B_H
