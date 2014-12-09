#ifndef HJ_MESH_LAPLACIAN_H_
#define HJ_MESH_LAPLACIAN_H_

#include "conf.h"

#include <zjucad/matrix/matrix.h>

namespace hj { namespace mesh {

/// Laplacian
// 0 --> divide area
// 1 --> divide sqrt(area)
// 2 --> not divide
bool HJ_MESH_API cot_Laplacian_val(const zjucad::matrix::matrix<double> &v, const zjucad::matrix::matrix<int> &t,
					const zjucad::matrix::matrix<int> &ptr, const zjucad::matrix::matrix<int> &idx,
					zjucad::matrix::matrix<double> &val, int weighting = 0);

bool HJ_MESH_API cot_Laplacian(const zjucad::matrix::matrix<double> &v, const zjucad::matrix::matrix<int> &t,
				   zjucad::matrix::matrix<int> &ptr, zjucad::matrix::matrix<int> &idx,
				   zjucad::matrix::matrix<double> &val, int weighting = 0);

}}

#endif
