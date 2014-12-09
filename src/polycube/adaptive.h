#ifndef HJ_POLYCUBE_ADAPTIVE_H_
#define HJ_POLYCUBE_ADAPTIVE_H_


#include <algorithm>
#include <zjucad/matrix/matrix.h>

//! NOTICE: assume that no tet has two boundary faces!
//! @output: child_node = [parent_node, new_node]
int boundary_uniform_subdivide(const zjucad::matrix::matrix<double> &parent_node,
                               const zjucad::matrix::matrix<size_t> &parent_tet,
                               zjucad::matrix::matrix<double> &child_node,
                               zjucad::matrix::matrix<size_t> &child_tet,
                               zjucad::matrix::matrix<size_t> &new_node_parent_id//nx2 to ori node
                               );

#endif
