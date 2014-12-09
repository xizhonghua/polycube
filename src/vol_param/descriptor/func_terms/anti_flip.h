#ifndef ANTI_FLIP_H
#define ANTI_FLIP_H

#include <zjucad/matrix/matrix.h>
#include "def.h"

template<SIMPLEX>
class anti_flip_func : jtf_func{};

template <>
class anti_flip_func<TET> : jtf_func
{
public:
  anti_flip_func(
      const zjucad::matrix::matrix<double> & node,
      const zjucad::matrix::matrix<size_t> & tet,
      const size_t & idx, const double w)
    :tet_(tet), idx_(idx), weight_(w){}
  virtual ~anti_flip_func(){}
public:
  virtual size_t dim_of_x()const{
  }
private:
  const size_t idx_;
  const double weight_;
  const zjucad::matrix::matrix<size_t> &tet_;
};

#endif // ANTI_FLIP_H
