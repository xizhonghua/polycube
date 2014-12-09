#ifndef HJ_HEX_ZYZ_H_
#define HJ_HEX_ZYZ_H_

#include "config.h"

#ifdef __cplusplus
extern "C" {
#endif

/**
   rot_matrix is a column-major 3x3
   Z(gamma)*Y(beta)*Z(alpha);
 */

HEXGEN_COMMON_API
void zyz_angle_2_rotation_matrix1(
	const double *alpha_beta_gamma,
	double *rot_matrix);

HEXGEN_COMMON_API
void zyz_angle_2_rotation_matrix(
	double alpha_, double beta_, double gamma_,
	double *rot_matrix);

HEXGEN_COMMON_API
void rotation_matrix_2_zyz_angle(const double *rot_matrix, double *abc, double *err);

HEXGEN_COMMON_API
void rot_n_2_z_by_zyz(const double *n, double *abc);

#ifdef __cplusplus
}
#endif

#include <zjucad/matrix/matrix.h>
void zyz2frame(const zjucad::matrix::matrix<double> & zyz,
	       zjucad::matrix::matrix<zjucad::matrix::matrix<double> > & frame);

void frame2zyz(const zjucad::matrix::matrix<zjucad::matrix::matrix<double> > & frame,
	       zjucad::matrix::matrix<double> & zyz);

#endif
