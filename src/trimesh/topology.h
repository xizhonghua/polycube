#ifndef _HJ_MESH_TOPOLOGY_H_
#define _HJ_MESH_TOPOLOGY_H_

#include "conf.h"

#include <zjucad/matrix/matrix.h>
#include <string>
#include <map>
#include <vector>

namespace hj { namespace mesh {

/// cell2graph
void HJ_MESH_API cell2graph(const zjucad::matrix::matrix<int> &cell,
				zjucad::matrix::matrix<zjucad::matrix::matrix<int> > &graph,
				bool self);

void HJ_MESH_API cell2csc(const zjucad::matrix::matrix<int> &cell,
				zjucad::matrix::matrix<int> &ptr,
				zjucad::matrix::matrix<int> &idx,
				bool self);

}}

#endif
