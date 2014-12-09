#ifndef JTF_RELATION_H
#define JTF_RELATION_H


//! @param v: 1x1, x:11x1, rot:9x9
extern "C" void relate_diff_energy_(double *v, double *x, double *rot);

//! @param jac: 11x1, x:11x1, rot:9x9
extern "C" void relate_diff_energy_jac_(double *jac, double *x, double *rot);

//! @param hes:11x11, x:11x1, rot:9x9
extern "C" void relate_diff_energy_hes_(double * hes, double *x, double *rot);

#endif // RELATION_H
