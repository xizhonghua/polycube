#ifndef TRANSITION_TYPE_H
#define TRANSITION_TYPE_H

#include <zjucad/matrix/matrix.h>

static const size_t TRIVIAL_TYPE = 9;

zjucad::matrix::matrix<double> type_transition2(const size_t type);

size_t type_transition1(const zjucad::matrix::matrix<double>& transition);

inline bool is_trivial_type(const size_t type){
  if(type == TRIVIAL_TYPE) return true;
  return false;
}

inline bool is_black_line_new(const size_t type){
  if(type > TRIVIAL_TYPE) return true;
  else
    return false;
}

inline bool is_regular_type(const size_t type){
  if(type < TRIVIAL_TYPE) return true;
  else
    return false;
}

inline size_t get_trans_type(const size_t  & type){
  assert(type < 24);
  return type_transition1(trans(type_transition2(type)));
}

//! @brief return axis which is around with,
// if rotation type > 8, then this rotation is not a legal rotation
// which should be one axis rotation, return -1
// else return 0,1,2 // which means u,v,w
inline size_t axis_to_around(const size_t rotation_type){
  if(rotation_type > 8) // 9,10,...,23 is error
    return -1;
  return rotation_type/3;
}

inline std::pair<size_t,size_t> get_compound_axis(const size_t rotation_type){
  assert(rotation_type < 24);
  if(rotation_type <= TRIVIAL_TYPE)
    return std::make_pair(-1,-1);
  if(rotation_type < 22 && rotation_type > 9)
    return std::make_pair(0,1); // u and v
  else
    return std::make_pair(0,2);
}
#endif // TRANSITION_TYPE_H
