#ifndef IO_H
#define IO_H
#include "../common/def.h"

namespace jtf{
namespace hexmesh{
int hex_mesh_read_from_wyz(const char *path,
                           matrixst &hex,
                           matrixd &node,
                           const size_t hex_format,
                           bool is_remove_invalid_hex = false);

int hex_mesh_write_to_wyz(const char *path,
                           matrixst &hex,
                           matrixd &node);
///
/// @brief hex_mesh_read_from_vtk
/// @param path
/// @param node
/// @param hex
/// @return
///
int hex_mesh_read_from_vtk(
    const char *path,
    matrixd *node = 0,
    matrixst *hex = 0);

}
}

#endif // IO_H
