#ifndef HJ_ROT_CUBIC_F_SH_H_
#define HJ_ROT_CUBIC_F_SH_H_

//! @param a, b, c denotes the Z(c)Y(b)Z(a)*vec angle respectively.

extern "C" {
void calc_rot_cubic_f_sh_mat_(double *rot_cubic_f_SH_mat,
							  const double *abc);

//! @param rot_cubic_f_SH 9x1
void calc_rot_cubic_f_sh_(double *rot_cubic_f_SH,
						  const double *abc);

//! @param Jac_rot_cubic_f_SH 9x3
void calc_jac_rot_cubic_f_sh_(double *Jac_rot_cubic_f_SH,
							  const double *abc);

//! @param SH is assume to be normalized, return lmdif info
int SH2zyz(double *zyz, const double *SH, unsigned int iter_num);
int SH2zyz_lmder1(double *zyz, const double *SH, unsigned int iter_num);
}

#endif
