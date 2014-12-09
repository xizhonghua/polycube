#include <jtflib/mesh/io.h>

#include "read_tet.h"
#include "../../tetmesh/tetmesh.h"

namespace lq {


void read_tet::load_data(const char *path, matrixd *node,
                         matrixst *tet, matrixst *tri)
{
  jtf::mesh::tet_mesh_read_from_zjumat(path, node, tet, tri);
}

}
