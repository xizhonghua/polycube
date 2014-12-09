#include "face_orient.h"
#include <vector>
namespace dzw {

int align_axis(const matrixd& face_normal)
{
  auto max_it = std::max_element(zjucad::matrix::fabs(face_normal).begin(), zjucad::matrix::fabs(face_normal).end());
  int idx = 0;
  for(auto it = zjucad::matrix::fabs(face_normal).begin(); it != max_it; ++it, ++idx);
  idx *= *max_it>0?1.0:-1.0;
  return idx;
}

bool is_same_orient_two_faces(const matrixd& face_normal_1, const matrixd& face_normal_2)
{
  if(align_axis(face_normal_1) == align_axis(face_normal_2))
    return true;
  return false;
}

}
