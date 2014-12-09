#ifndef CELL_QUALITY_H
#define CELL_QUALITY_H

#include <verdict/verdict.h>
#include <zjucad/matrix/matrix.h>

//! @param tet_node 3x4
template <typename T>
double tet_scaled_jacobian(const zjucad::matrix::matrix_expression<T> & tet_node)
{
  // tet_node is 3 * 4
  VERDICT_REAL tet_node_coor[4][3];
  for(size_t pi = 0; pi < 4; ++pi){
    for(size_t di = 0; di < 3; ++di){
      tet_node_coor[pi][di] = tet_node()(di, pi);
    }
  }
  return v_tet_scaled_jacobian(tet_node().size(2), tet_node_coor);
}

//! @param hex_node 3x8
template <typename T>
double hex_scaled_jacobian(const zjucad::matrix::matrix_expression<T> & hex_node)
{
  // hex_node is 3 * 8
  VERDICT_REAL hex_node_coor[8][3];
  const size_t p_order[] = {7,5,4,6,3,1,0,2};
  for(size_t pi = 0; pi < 8; ++pi){
    for(size_t di = 0; di < 3; ++di){
      hex_node_coor[pi][di] = hex_node()(di, p_order[pi]);
    }
  }
  return v_hex_scaled_jacobian(hex_node().size(2), hex_node_coor);
}

//! @param hex_node 3x8
template <typename T>
double hex_volume(const zjucad::matrix::matrix_expression<T> & hex_node)
{
  VERDICT_REAL hex_node_coor[8][3];
  const size_t p_order[] = {7,5,4,6,3,1,0,2};
  for(size_t pi = 0; pi < 8; ++pi){
    for(size_t di = 0; di < 3; ++di){
      hex_node_coor[pi][di] = hex_node()(di, p_order[pi]);
    }
  }
  return v_hex_volume(hex_node().size(2), hex_node_coor);
}

#endif // CELL_QUALITY_H
