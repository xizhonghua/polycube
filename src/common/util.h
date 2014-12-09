#ifndef HJ_HEX_UTIL_H_
#define HJ_HEX_UTIL_H_

#include <vector>
#include <deque>
#include <numeric>
#include <map>
#include <set>
#include <boost/unordered_set.hpp>
#include <boost/static_assert.hpp>
#include <zjucad/matrix/matrix.h>
#include <zjucad/matrix/itr_matrix.h>

#include "../common/def.h"

///
/// @brief calc_bounding_box
/// @param input nodes
/// @param output bounding box with min,max point
///
inline void calc_bounding_box(const matrixd &node, double *bb)
{
  using namespace zjucad::matrix;
  const size_t dim = node.size(1);
  for(int i = 0; i < node.size(1); ++i) {
      bb[i] = min(node(i, colon()));
      bb[i+dim] = max(node(i, colon()));
    }
}

///
/// @brief calc_bounding_box_sie
/// @param input nodes
/// @param output size, xyz axes length
///
inline void calc_bounding_box_size(const matrixd &node, double *size)
{
  std::vector<double> bb(node.size(1)*2);
  calc_bounding_box(node, &bb[0]);
  for(int i = 0; i < node.size(1); ++i)
    size[i] = bb[i+node.size(1)]-bb[i];
}

///
/// @brief calc_bounding_sphere_size: 2 * radius
/// @param input nodes
/// @param output size: 2 * radius
///
inline double calc_bounding_sphere_size(const matrixd &node)
{
  double rtn = 0;
  std::vector<double> size(node.size(1));
  calc_bounding_box_size(node, &size[0]);
  for(int i = 0; i < node.size(1); ++i)
    rtn += size[i]*size[i];
  return sqrt(rtn);
}

///
/// @brief calc max bounding box edge length
/// @param input nodes
/// @param output max size
///
inline double calc_max_bounding_box_edge_length(const matrixd &node)
{
  std::vector<double> size(node.size(1));
  calc_bounding_box_size(node, &size[0]);
  return *std::max_element(size.begin(), size.end());
}

///
/// @brief calc average node
/// @param node: input nodes
/// @param average_node: output average node
///
inline void cal_average_node(const matrixd & node, matrixd & average_node_)
{
  average_node_ = node * zjucad::matrix::ones<double>(node.size(2),1) / node.size(2);
}

inline void remove_extra_node(zjucad::matrix::matrix<size_t> & cell,
                              zjucad::matrix::matrix<double> & node,
                              zjucad::matrix::matrix<size_t> * orig2new_mapping = 0)
{
  using namespace std;
  using namespace zjucad::matrix;
  set<size_t> used_node_idx(cell.begin(), cell.end());
  if(used_node_idx.size() == node.size(2)) return;
  matrixst used_node_mat(used_node_idx.size(),1);
  copy(used_node_idx.begin(), used_node_idx.end(), used_node_mat.begin());

  map<size_t,size_t> p2p;

  matrixd new_node(3, used_node_mat.size());
  for(size_t pi = 0; pi < used_node_mat.size(); ++pi){
      new_node(colon(),pi) = node(colon(), used_node_mat[pi]);
      p2p[used_node_mat[pi]] = pi;
    }
  for(size_t pi = 0; pi < cell.size(); ++pi)
    cell[pi] = p2p[cell[pi]];

  if(orig2new_mapping != 0){
      *(orig2new_mapping) = zjucad::matrix::ones<size_t>(node.size(2),1) * -1;
      for(map<size_t,size_t>::const_iterator cit = p2p.begin(); cit != p2p.end();
          ++cit){
          (*orig2new_mapping)[cit->first] = cit->second;
        }
    }

  node = new_node;
}

///
/// \brief remove_extra_node
/// \param mesh     input mesh/output mesh
/// \param mapping  mapping original point to new one
///
inline void remove_extra_node(zjucad::matrix::matrix<size_t> & mesh,
                              std::map<size_t,size_t> &orig2new_mapping)
{
  const size_t max_size = *std::max_element(mesh.begin(), mesh.end());
  std::set<size_t> compress_v(mesh.begin(), mesh.end());
  std::vector<size_t> compress_vec(compress_v.size());
  std::copy(compress_v.begin(), compress_v.end(), compress_vec.begin());
  if(max_size == compress_vec.size()-1) return ;
  std::cerr << "orig mesh max point: " << max_size
            << " compressed mesh max point: " << compress_vec.size()-1 << std::endl;
  orig2new_mapping.clear();
  for(size_t i = 0; i < compress_vec.size(); ++i){
      orig2new_mapping[compress_vec[i]] = i;
    }
  for(size_t j = 0; j < mesh.size(); ++j)
    mesh[j] = orig2new_mapping[mesh[j]];
}

inline void remove_duplicated_node(zjucad::matrix::matrix<size_t> & mesh,
                                   zjucad::matrix::matrix<double> & node,
                                   const double epsilon = 1e-3)
{
  using namespace std;
  using namespace zjucad::matrix;

  const double d = calc_bounding_sphere_size(node);
  map<size_t,size_t> p2p;
  std::vector<bool> remain(node.size(2), true);
  for(size_t i = 0; i < remain.size(); ++i){
      if(remain[i] == false) continue;
      for(size_t j = i + 1; j < remain.size(); ++j){
          if(remain[j] == false) continue;
          if(norm(node(colon(), j) - node(colon(), i)) < d * epsilon) {
              remain[j]=false;
              p2p[j] = i;
            }
        }
    }

  for(size_t mi = 0; mi < mesh.size(); ++mi){
      const auto it = p2p.find(mesh[mi]);
      if(it ==  p2p.end()) continue;
      mesh[mi] = it->second;
    }
  remove_extra_node(mesh, node);
}
#endif
