#ifndef QUAD_UTIL_H
#define QUAD_UTIL_H

#include "../common/def.h"
#include <jtflib/mesh/mesh.h>
namespace jtf
{
namespace quadmesh
{
int extract_quadmesh_singularity_line(const char *path1,
                                      const char *path2);

void find_vertex_lines(size_t index,
                       std::vector<std::vector<size_t> > &lines,
                       const jtf::mesh::edge2cell_adjacent &edge2quad,
                       const std::vector<std::vector<size_t> > &adj_ver,
                       const matrixst &quad, std::vector<std::vector<bool> > &is_visited);

void choose_next_node(size_t &cur_node, size_t &next_node,
                      const jtf::mesh::edge2cell_adjacent &edge2quad,
                      const std::vector<std::vector<size_t> > &adj_ver,
                      const matrixst &quad);

int write_singularity_line(const char* path2,
                           const std::vector<std::vector<size_t> > &singularity_lines);


inline bool is_in_vec(size_t index, const std::vector<size_t> &vec)
{
  for(size_t i = 0; i < vec.size(); ++i)
    if(vec[i] == index)
      return true;
  return false;
}
}

}

#endif // UTIL_H
