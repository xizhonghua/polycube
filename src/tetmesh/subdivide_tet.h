#ifndef HJ_SUBDIVDE_TET_H_
#define HJ_SUBDIVDE_TET_H_

#include <vector>

#include <zjucad/matrix/matrix.h>

inline void make_id(size_t *id, size_t n) {
  std::sort(id, id+n);
}

//! NOTICE: assume id is sorted
inline size_t find(const zjucad::matrix::matrix<size_t> &ids, const size_t *id) {
  using namespace std;
  size_t beg = 0, end = ids.size(1);
  for(size_t i = 0; i < ids.size(2); ++i) {
    const pair<const size_t *, const size_t *> r = equal_range(&ids(0, i)+beg, &ids(0, i)+end, id[i]);
    beg = r.first - &ids(0, i);
    end = r.second - &ids(0, i);
  }
  if(end - beg != 1)
    return -1;
  return beg;
}

bool is_sorted(const zjucad::matrix::matrix<size_t> &ids);

void add_edge_node4subdivide(const zjucad::matrix::matrix<size_t> &edge_node,
                             zjucad::matrix::matrix<size_t> &new_node);

void validate_edge_nodes(const zjucad::matrix::matrix<size_t> &tetmesh,
                         const zjucad::matrix::matrix<size_t> &edge_node,
                         zjucad::matrix::matrix<size_t> &new_nodes);

//! @NOTICE: The following routine require valide edge_node (through
//! the above) as input.

//! @input: edge_node is 2xn matrix
//! @output: children is a 2xn matrix point to local node_idx
int subdivide_tet(const zjucad::matrix::matrix<size_t> &edge_node,
                  zjucad::matrix::matrix<size_t> &children);

//! @input: edge_node is 2xn matrix
//! @output: children is a 2xn matrix point to node_idx
int subdivide_tet(const zjucad::matrix::matrix<size_t> &tet,
                  const zjucad::matrix::matrix<size_t> &edge_node,
                  zjucad::matrix::matrix<size_t> &children);

//! @input: NOTICE that edge_node must be sorted nx2
int subdivide_tetmesh(const zjucad::matrix::matrix<size_t> &tetmesh,
                      size_t parent_node_num,
                      const zjucad::matrix::matrix<size_t> &edge_node,
                      zjucad::matrix::matrix<size_t> &children);

class ordered_table_builder
{
public:
  ordered_table_builder(size_t dim);
  ordered_table_builder &operator << (const size_t *id);
  void operator >> (zjucad::matrix::matrix<size_t> &tab); // num x dim
private:
  std::vector<std::vector<size_t> > buf_;
  const size_t dim_;
};


//! @output child_node = [parent_node, new_node]
void subdivide_top2geo(const zjucad::matrix::matrix<double> &parent_node,
                       const zjucad::matrix::matrix<size_t> &new_node_parent_id, //nx2 to ori node
                       zjucad::matrix::matrix<double> &child_node
                       );

//! @ param: eles is in num x dim, but cells is in dim x num, TODO: fix it
void create_ordered_table(const zjucad::matrix::matrix<size_t> &cells,
                          const zjucad::matrix::matrix<size_t> &pattern,
                          zjucad::matrix::matrix<size_t> &eles);

#endif
