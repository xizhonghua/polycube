#ifndef HJ_HEX_IO_2_H_
#define HJ_HEX_IO_2_H_

#include <zjucad/matrix/matrix.h>
#include "../tetmesh/tetmesh.h"

int surface_frame2volme_frame_inner(const char *surf_node_uv_file, const char *surf_node2tet_node_file,
                              const char *obj_file, const jtf::mesh::face2tet_adjacent &fa,
                              matrixd &surface_zyz_frame,
                              matrixst &idx);
int load_from_tet_inner(
    const matrixd &node, const matrixst &tet,
    matrixd &fixed_frame, matrixst &fixed_frame_idx,
    matrixd &aligned, matrixst &aligned_idx,
    const char *surf_node_uv_file, const char *surf_node2tet_node_file, const char *obj_file
    );
#endif
