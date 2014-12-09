#ifndef HJ_HEX_IO_H_
#define HJ_HEX_IO_H_

#include "config.h"

#include <stdint.h>

#include <fstream>
#include <iostream>

#include "graph.h"

#include "../common/def.h"
#include <hjlib/sparse/sparse.h>
#include <jtflib/mesh/io.h>
/**
   \param node_pos 3*node_num
 */
HEXGEN_COMMON_API
int load_graph(const char *path, fix_graph &fg,
               matrixd &node_pos);

/**
   \param rot 9*node_num matrix, each column is a packed 3*3 column
   major rotation matrix.
 */
HEXGEN_COMMON_API
int load_vec_field(const char *path, matrixd &rot);

HEXGEN_COMMON_API
int save_param(const char *path,
	       const matrixd &uvw,
	       const matrixd &trans_rot,
	       const matrixd &trans_offset);

template <typename VAL_TYPE, typename INT_TYPE>
int read_csc(std::istream &is, hj::sparse::csc<VAL_TYPE, INT_TYPE> &m)
{
  int32_t nrow, ncol, nnz;
  is.read((char *)&nrow, sizeof(int32_t));
  is.read((char *)&ncol, sizeof(int32_t));
  is.read((char *)&nnz, sizeof(int32_t));
  m.resize(nrow, ncol, nnz);
  jtf::mesh::read_matrix(is, m.ptr_);
  jtf::mesh::read_matrix(is, m.idx_);
  jtf::mesh::read_matrix(is, m.val_);
  return is.fail();
}

template <typename VAL_TYPE, typename INT_TYPE>
int write_csc(std::ostream &os, const hj::sparse::csc<VAL_TYPE, INT_TYPE> &m)
{
  int32_t nrow = m.size(1), ncol = m.size(2), nnz = hj::sparse::nnz(m);
  os.write((const char *)&nrow, sizeof(int32_t));
  os.write((const char *)&ncol, sizeof(int32_t));
  os.write((const char *)&nnz, sizeof(int32_t));
  jtf::mesh::write_matrix(os, m.ptr_);
  jtf::mesh::write_matrix(os, m.idx_);
  jtf::mesh::write_matrix(os, m.val_);
  return os.fail();
}

inline int read_zyz(const char * zyz_file,
                    zjucad::matrix::matrix<double> &tet_zyz)
{
  std::ifstream ifs(zyz_file, std::ifstream::binary);
  if(ifs.fail()){
      std::cerr << "# [error] can not read zyz file." << std::endl;
      return __LINE__;
    }
  if(jtf::mesh::read_matrix(ifs, tet_zyz)) {
      std::cerr << "# not a zyz file." << std::endl;
      return __LINE__;
    }
  return 0;
}

inline int write_zyz(const char * zyz_file,
                     zjucad::matrix::matrix<double> & zyz)
{
  std::ofstream ofs_aligned(zyz_file, std::ofstream::binary);
  if(ofs_aligned.fail()){
      std::cerr << "# [error] can not open zyz file." << std::endl;
      return __LINE__;
    }
  jtf::mesh::write_matrix(ofs_aligned, zyz);
  return 0;
}

#endif
