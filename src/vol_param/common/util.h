#ifndef VOL_PARAM_COMMON_UTIL_H
#define VOL_PARAM_COMMON_UTIL_H

#include <zjucad/matrix/matrix.h>
#include <fstream>
#include <jtflib/mesh/io.h>
#include <jtflib/math/math.h>
#include <hjlib/sparse/sparse.h>
///
/// @brief str2int  Fast simple string hash (Bernstein?)
// "http://stackoverflow.com/questions/2111667/compile-time-string-hashing"
/// @param s
/// @param off
/// @return
///
///
constexpr unsigned int str2int(const char *s, int off = 0) {
  return !s[off]?5381:(str2int(s, off+1)*33)^s[off];
}

//inline int read_zyz(const char * zyz_file,
//                    zjucad::matrix::matrix<double> &tet_zyz)
//{
//  std::ifstream ifs(zyz_file, std::ifstream::binary);
//  if(ifs.fail()){
//      std::cerr << "# [error] can not read zyz file." << std::endl;
//      return __LINE__;
//    }
//  if(jtf::mesh::read_matrix(ifs, tet_zyz)) {
//      std::cerr << "# not a zyz file." << std::endl;
//      return __LINE__;
//    }
//  return 0;
//}


inline bool is_degenerate(const zjucad::matrix::matrix<double> &tri)
{
  using namespace zjucad::matrix;

  matrix<double> e1 = tri(colon(), 0)-tri(colon(), 2),
      e2 = tri(colon(), 1)-tri(colon(), 2);
  if(norm(cross(e1, e2)) < 1e-8) {
      return true;
    }
  return false;
}

template <typename ITERATION1, typename ITERATION2>
double csc_one_line_multy(ITERATION1 idx_begin,
                          ITERATION1 idx_end,
                          ITERATION2 val_begin,
                          ITERATION2 val_end,
                          const double * x)
{
  double v = 0;
  ITERATION2 v_it = val_begin;
  for(ITERATION1 it = idx_begin; it != idx_end; ++it){
      assert(v_it != val_end);
      v += x[*it] * (*v_it);
      ++v_it;
    }
  return v;
}

///
/// \brief get_node_mapping, x=Zu,
/// \param one_cell
/// \param NMT
/// \param node_mapping_of_each_variable
///
inline void get_node_mapping(
    const zjucad::matrix::matrix<size_t>& one_cell,
    const hj::sparse::csc<double,int32_t> & NMT,
    std::vector<std::vector<std::pair<size_t,double> > > & node_mapping_of_each_variable)
{
  const size_t dim = 3;
  node_mapping_of_each_variable.resize(dim * one_cell.size());
  for(size_t pi = 0; pi < one_cell.size(); ++pi){
      for(size_t di = 0; di < dim; ++di){
          size_t nnz_pi = NMT.ptr()[dim * one_cell[pi] + di + 1] -
              NMT.ptr()[dim * one_cell[pi] + di];

          node_mapping_of_each_variable[dim * pi + di].reserve(nnz_pi);
          size_t begin = NMT.ptr()[dim * one_cell[pi] + di];

          for(size_t pj = 0; pj < nnz_pi; ++pj){
              node_mapping_of_each_variable[dim * pi + di].push_back(
                    std::make_pair(NMT.idx()[begin+pj],
                    NMT.val()[begin+pj]));
            }
        }
    }
}



class node_mapping
{
public:
  node_mapping(){}
  ~node_mapping(){}
  // x = Zu + q;
  hj::sparse::csc<double,int32_t> ZT;
  zjucad::matrix::matrix<double> q;
  zjucad::matrix::matrix<int> is_fixed; // 0: no, 1: yes
  ///
  /// \brief is_irrelevant_variable, to check whether given variable is relevant with real variables
  /// \param idx
  /// \return
  ///
  bool is_irrelevant_variable(const size_t idx)const{
    if(idx > ZT.size(2)) return true;
    const size_t begin = ZT.ptr()[idx];
    const size_t end = ZT.ptr()[idx+1];
    if(begin == end) return true;
    return false;
  }
};



///
/// \brief get_cell_node, get cell node under node mapping
/// \param x
/// \param one_cell
/// \param node_mapping
/// \param cell_node
///
inline void get_cell_node(const double *x,
                          const zjucad::matrix::matrix<size_t> & one_cell,
                          const node_mapping & NM,
                          zjucad::matrix::matrix<double> & cell_node)
{
  const size_t dim = 3;
  cell_node.resize(dim, one_cell.size());
  assert(NM.q.size() == NM.ZT.size(2));
  for(size_t pi = 0; pi < one_cell.size(); ++pi){
      for(size_t di = 0; di < dim; ++di){
          const size_t variable_idx = dim * one_cell[pi] + di;
          cell_node(di,pi) = csc_one_line_multy(
                &(NM.ZT.idx()[0]) + NM.ZT.ptr()[variable_idx],
              &(NM.ZT.idx()[0]) + NM.ZT.ptr()[variable_idx+1],
              &(NM.ZT.val()[0]) + NM.ZT.ptr()[variable_idx],
              &(NM.ZT.val()[0]) + NM.ZT.ptr()[variable_idx+1],
              x) + NM.q[variable_idx];
        }
    }
}

// smooth triangle surface using degree laplacian.
void smooth_mesh_surface(const zjucad::matrix::matrix<size_t> & tri_faces,
                         zjucad::matrix::matrix<double> & node,
                         const size_t iter);

#endif
