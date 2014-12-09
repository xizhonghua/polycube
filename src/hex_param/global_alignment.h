#ifndef GLOBAL_ALIGNMENT_H
#define GLOBAL_ALIGNMENT_H

#include <zjucad/matrix/io.h>
#include <zjucad/matrix/itr_matrix.h>
#include <zjucad/matrix/matrix.h>
#include "../tetmesh/tetmesh.h"
#include <vector>
#include <deque>
#include <map>

void hj_frame_alignemt(const matrixst& tet,
                       const jtf::mesh::face2tet_adjacent &fa,
                       const matrixd &tet_zyz,
                       matrixd &frame_inner);
#endif // GLOBAL_ALIGNMENT_H
