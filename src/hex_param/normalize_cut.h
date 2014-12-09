#ifndef NORMALIZE_CUT_H
#define NORMALIZE_CUT_H

#include <vector>
#include <algorithm>
#include "../tetmesh/tetmesh.h"

/// @brief: is the index in c cover all elements from [0, max(c)-1]
template <typename Con>
bool is_compact_index(const Con &c)
{
  if(c.size() == 0) return true;

  std::vector<typename Con::value_type> v(c.size());
  std::copy(c.begin(), c.end(), v.begin());
  std::sort(v.begin(), v.end());
  if(v[0] != 0) return false;
  return v[v.size()-1] == (std::unique(v.begin(), v.end()) - v.begin() - 1);
}

/// @brief: assume tet is compact
bool is_normalized_cut_tet(const matrixst &tet,
                           const matrixst &cut);

int normalize_cut_tet_node_idx(const matrixst &tet,
                                matrixst &cut);
#endif // NORMALIZE_CUT_H
