#ifndef NUMERIC_UTIL_H
#define NUMERIC_UTIL_H
#include <cmath>
#include <deque>
#include <iostream>
#include <limits>

#include <boost/dynamic_bitset.hpp>
#include <boost/type_traits/is_same.hpp>
#include <boost/math/constants/constants.hpp>

#include <hjlib/sparse/sparse.h>
#include <jtflib/math/math.h>
#include "../common/def.h"

constexpr double My_PI()
{
  return 3.1415926535897932384626;
}

inline double float_mod(const double a,
                        const double b)
{
  assert(fabs(b) > 1e-6);
  int N = a/b;
  return a - N*b;
}

inline double calculate_dihedral_angle_degree(
    const zjucad::matrix::matrix<double> & n1,
    const zjucad::matrix::matrix<double> & n2)
{
  using namespace zjucad::matrix;
  const double len1 = norm(n1);
  const double len2 = norm(n2);
  matrix<double> n1_norm = n1, n2_norm = n2;
  if(len1 > 1e-6) n1_norm /= len1;
  if(len2 > 1e-6) n2_norm /= len2;

  return jtf::math::cos2deg(dot(n1_norm, n2_norm));
}
// U,V,W construct a triangle of tet,
// and u,v,w is the opposite edge of U,V,W
inline double tet_vol(const double U, const double V, const double W,
                      const double u, const double v, const double w)
{
  const double X = (w-U+v)*(U+v+w);
  const double x = (U-v+w)*(v-w+U);
  const double Y = (u-V+w)*(V+w+u);
  const double y = (V-w+u)*(w-u+V);
  const double Z = (v-W+u)*(W+u+v);
  const double z = (W-u+v)*(u-v+W);
  const double a = sqrt(x*Y*Z);
  const double b = sqrt(y*Z*X);
  const double c = sqrt(z*X*Y);
  const double d = sqrt(x*y*z);
  return sqrt((-a+b+c+d)*(a-b+c+d)*(a+b-c+d)*(a+b+c-d))/(192*u*v*w);
}

inline double tet_vol(const double *e)
{
  assert(e);
  return tet_vol(e[0],e[1],e[2],e[3],e[4],e[5]);
}

/**
 * @brief convert a rotation matrix into euler zyx representation
 * there are solutions when cos(y) not equal 0,
 * and infinity solutions when cos(y) = 0,
 * so as normal, we output two solutions
 * input should be formate 3*3, colon major
 * output will be 3 * 2, colon major
 * WARNING: it works for 24 basic rotation, but failed on rand rotation
 * TODO: fix this bug!
 *
 * @param rot input rotation matrix
 * @param zyx output zyx angle
 * @return int
 */
int convert_rot_to_euler_zyx(const double *rot,double *zyx);

/**
 * @brief convert euler angle to rotation matrix
 *
 * @param zyx input zyx angle
 * @param rot output rot matrix
 * @return int
 */
int convert_euler_angle_to_rot(const double *zyx,double *rot);


/**
 * @brief quat: [w, x, y, z]
 *         zyx: [z, y, x]
 *
 * @param quat
 * @param zyx
 * @return int
 */
int convert_quat_to_euler_angle(const double *quat, double *zyx);

/**
 * @brief quat: [w, x, y, z]
 *         zyx: [z, y, x]
 *
 * @param zyx euler angle
 * @param quat
 * @return int
 */
int convert_euler_angle_to_quat(const double *zyx, double *quat);

/**
 * @brief    calculate a rotation matrix which describe a rotation
 *           around axis_in, and the rotation angle is given by angle
 *
 * @param angle     input rotation angle around axis
 * @param axis_in   input rotation axis
 * @param rot       output rotation matrix
 * @return int
 */
int from_angle_to_rotation_matrix(
    const double &angle,
    const matrixd &axis_in,
    matrixd &rot);

//! @brief convert quaternion to rotation matrix
//! @param quat input quaternion [qx,qy,qz,qw]
//! @param rot  output rotation matrix
void convert_quat_to_rotation_matrix(const double * quat, matrixd & rot);


//! @brief convert a rotation matrix to quaternion, details in http://planning.cs.uiuc.edu/node153.html
//! @param rot input rotation matrix
//! @param quat output quaternion [qx,qy,qz,qw]
void convert_rotation_matrix_to_quat(const matrixd &rot, double *quat);


//! @param quat : [x,y,z,w]
template <typename T>
void convert_quat_to_axis_angle(T const *quat, T *axis, T &angle, bool anti_order = false)
{
  assert(quat && axis);
  T w,x,y,z;
  if(anti_order){
      w = quat[0];
      x = quat[1];
      y = quat[2];
      z = quat[3];
    }else{
      w = quat[3];
      x = quat[0];
      y = quat[1];
      z = quat[2];
    }
  T len = sqrt(w*w+x*x+y*y+z*z);
  w/=len; x/=len; y/=len; z/=len;

  angle = 2*std::acos(std::max(-1.0, std::min(w,1.0)));
  if(fabs(angle) < 1e-6) {
      axis[0] = 1.0;
      axis[1] = 0;
      axis[2] = 0;
    }else{
      const T sqrt_1_qw_qw = 1.0/std::sqrt(1-w*w);
      axis[0] = x * sqrt_1_qw_qw;
      axis[1] = y * sqrt_1_qw_qw;
      axis[2] = z * sqrt_1_qw_qw;
    }
}


template <typename T1, typename T2>
void convert_rotation_matrix_to_axis_angle(
    const zjucad::matrix::matrix_expression<T1> & rot,
    zjucad::matrix::matrix_expression<T2> & axis,
    double &angle)
{
  axis().resize(3,1);
  double quad[4];
  convert_rotation_matrix_to_quat(rot(),&quad[0]);
  convert_quat_to_axis_angle(&quad[0],&axis()[0],angle);
}


template <typename T1, typename T2>
void get_2d_rotation(T1 angle, zjucad::matrix::matrix_expression<T2> &rot)
{
  if(rot().size(1) != 2 || rot().size(2) != 2)
    rot().resize(2,2);
  rot()(0,0) = cos(angle);
  rot()(1,1) = cos(angle);
  rot()(0,1) = -sin(angle);
  rot()(1,0) = sin(angle);
}

template <typename T1, typename T2>
void get_3d_z_rotation(T1 angle, zjucad::matrix::matrix_expression<T2> &rot)
{
  if(rot().size(1) != 3 || rot().size(2) != 3)
    rot() = zjucad::matrix::eye<double>(3);
  rot()(0,0) = cos(angle);
  rot()(1,1) = cos(angle);
  rot()(0,1) = -sin(angle);
  rot()(1,0) = sin(angle);
}

template <typename T>
inline int number_rounding(const T & value)
{
  BOOST_STATIC_ASSERT((boost::is_same<T,double>::value)||
                      (boost::is_same<T,float>::value));
  return std::floor(value+0.5);
}

int long_number_add(std::deque<size_t> & data, size_t  a);

template <typename OS>
void print_long_number(OS &os, const std::deque<size_t> & data)
{
  for(size_t i = data.size() -1; i != -1; --i){
      os << data[i];
    }
  os << std::endl;
}


//double eigen_val_sm(const hj::sparse::csc<double, int32_t> & csc_);

template <typename VAL_TYPE, typename INT_TYPE>
inline void diag(const hj::sparse::csc<VAL_TYPE, INT_TYPE> &A, zjucad::matrix::matrix<VAL_TYPE> &d)
{
  d = zjucad::matrix::zeros<VAL_TYPE>(A.size(1), 1);
  for(INT_TYPE ci = 0; ci < A.size(2); ++ci) {
      for(INT_TYPE nzi = A.ptr()[ci]; nzi < A.ptr()[ci+1]; ++nzi) {
          if(A.idx()[nzi] == ci)
            d[ci] = A.val()[nzi];
        }
    }
}

template <typename VAL>
int round_integer(VAL &v)
{
  int low = std::floor(v);
  int ceil = std::ceil(v);
  if(std::fabs(v-low) < std::fabs(v-ceil))
    return low;
  else
    return ceil;
}

//! @brief find the eigenvector with smallest eigenvaluet
void solve_smallest_eigenvalue_problem(const hj::sparse::csc<double, int32_t> &A,
                                       zjucad::matrix::matrix<double> &eigen_vector);

void solve_smallest_eigenvalue_problem(const hj::sparse::csc<std::complex<double>, int32_t> &A,
                                       zjucad::matrix::matrix<std::complex<double> > &eigen_vector);

void cal_deformation_gradient(const zjucad::matrix::matrix<size_t> & mesh,
                              const zjucad::matrix::matrix<double> & node0,
                              const zjucad::matrix::matrix<double> & node1,
                              zjucad::matrix::matrix<zjucad::matrix::matrix<double> > & deform_gradient);
#endif // UTIL_H
