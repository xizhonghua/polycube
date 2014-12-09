#include "normalize_cut.h"
using namespace std;
using namespace zjucad::matrix;

/// @brief: assume tet is compact
bool is_normalized_cut_tet(const matrixst &tet,
                           const matrixst &cut)
{
  is_compact_index(cut);
  const size_t ori_node_num = max(tet)+1;
  matrixst cut2tet(max(cut)+1);
  cut2tet(cut) = tet(colon());
  for(size_t i = 0; i < ori_node_num; ++i)
    if(cut2tet[i] != i)
      return false;
  return true;
}

/// @brief: assume that the node index in tet is compact and in range
/// of [0, n-1], make the cut tet node index is also compact, and the
/// original node index of ith node (i < n) in cut is also i.
int normalize_cut_tet_node_idx(const matrixst &tet,
                                matrixst &cut)
{
  if(!is_compact_index(tet)) {
    cerr << "index in tet is not compact." << endl;
    return __LINE__;
  }
  cout << "min max index: " << min(cut) << ", " << max(cut) << endl;
  const size_t max_cut_node_num = max(cut)+1;
  matrixst cut_tet2tet = zeros<size_t>(max_cut_node_num, 1)-1,
    cut_tet2new_cut_tet = zeros<size_t>(max_cut_node_num, 1)-1;
  cut_tet2tet(cut) = tet(colon());

  const size_t ori_node_num = max(tet)+1;
  matrix<char> is_index_occupied = zeros<char>(ori_node_num, 1);
  size_t available_index = ori_node_num;
  for(size_t i = 0; i < cut.size(); ++i) {
    if(cut_tet2new_cut_tet[cut[i]] == -1) { // still not assign a new index
      const size_t ori_node_idx = cut_tet2tet[cut[i]];
      if(is_index_occupied[ori_node_idx])
        cut_tet2new_cut_tet[cut[i]] = available_index++;
      else {
        cut_tet2new_cut_tet[cut[i]] = ori_node_idx;
        is_index_occupied[ori_node_idx] = 1;
      }
    }
    cut[i] = cut_tet2new_cut_tet[cut[i]];
  }
  cout << "min max index: " << min(cut) << ", " << max(cut) << endl;
#if 1 // check
  if(!is_normalized_cut_tet(tet, cut))
    cerr << "output is not normalized." << endl;
  size_t uniq_node = 0;
  for(size_t i; i < cut_tet2tet.size(); ++i) {
    if(cut_tet2tet[i] != -1)
      ++uniq_node;
  }
  if(uniq_node != available_index) {
      cerr << "wrong num of uniq node in cut." << endl;
  }
  cerr << "compact check over." << endl;
  matrixst new_cut_tet2tet(max(cut)+1);
  new_cut_tet2tet(cut) = tet(colon());
  for(size_t i = 0; i < ori_node_num; ++i) {
    if(new_cut_tet2tet[i] != i)
      cerr << __LINE__ << endl;
  }
  cerr << "range check over." << endl;
  for(size_t i = 0; i < cut_tet2new_cut_tet.size(); ++i) {
    const size_t new_idx = cut_tet2new_cut_tet[i];
    if(new_idx == -1) continue;
    if(cut_tet2tet[i] != new_cut_tet2tet[new_idx])
      cerr << __LINE__ << endl;
  }
  cerr << "correct mapping check over." << endl;
  cerr << "check over." << endl;
#endif
  return 0;
}
