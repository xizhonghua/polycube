#ifndef WYZ_FORMAT_H
#define WYZ_FORMAT_H

#include <fstream>

#include <zjucad/matrix/matrix.h>
#include "../tetmesh/tetmesh.h"


/**
 * @brief load wyz's param file, there are 4 * tet_num vertex in
 *        parameterization field
 *
 * @param filename input wyz parameterization field file
 * @param tet_num  total tet number
 * @param param    output param matrix 3 * (4 * tet_num)
 * @return int     return 0 if load correctly, or non-zeros.
 */
template <typename T>
int load_wyz_param_file(const char * filename,
                        const size_t &tet_num,
                        zjucad::matrix::matrix<T> & param)
{
  std::ifstream ifs(filename,std::ios::binary);
  if(ifs.fail()){
    std::cerr << "# [error] can not read wyz's param file." << std::endl;
    return __LINE__;
  }
  param.resize(3, 4 * tet_num); //
  ifs.read((char*)&param[0],  sizeof(T)*param.size());
  return 0;
}

#endif // WYZ_FORMAT_H
