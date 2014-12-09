#include "zyz.h"

#include <algorithm>
#include <iostream>

#include "../common/def.h"
#include <zjucad/matrix/io.h>
#include <zjucad/matrix/itr_matrix.h>

using namespace std;
using namespace zjucad::matrix;

void zyz_angle_2_rotation_matrix(
    double alpha_, double beta_, double gamma_,
    double *rot_matrix)
{
  matrixd z1 = eye<double>(3);
  z1(0, 0) = cos(alpha_);
  z1(0, 1) = -sin(alpha_);
  z1(1, 0) = sin(alpha_);
  z1(1, 1) = cos(alpha_);

  matrixd y = eye<double>(3);
  y(0, 0) = cos(beta_);
  y(0, 2) = sin(beta_);
  y(2, 0) = -sin(beta_);
  y(2, 2) = cos(beta_);

  matrixd z2 = eye<double>(3);
  z2(0, 0) = cos(gamma_);
  z2(0, 1) = -sin(gamma_);
  z2(1, 0) = sin(gamma_);
  z2(1, 1) = cos(gamma_);

  matrixd mrot_matrix = temp(z2*y)*z1;
  std::copy(mrot_matrix.begin(), mrot_matrix.end(), rot_matrix);
}

void zyz_angle_2_rotation_matrix1(
    const double *alpha_beta_gamma,
    double *rot_matrix)
{
  zyz_angle_2_rotation_matrix(
        alpha_beta_gamma[0],
      alpha_beta_gamma[1],
      alpha_beta_gamma[2],
      rot_matrix);
}

void rotation_matrix_2_zyz_angle(const double *rot_matrix, double *abc, double *err)
{
  double &alpha = abc[0];
  double &beta = abc[1];
  double &gamma = abc[2];
  if(fabs(fabs(rot_matrix[8])-1) < 1e-9) {
      beta = acos((std::min)((std::max)(rot_matrix[8],-1.0),1.0));
      alpha = 0;
      gamma = atan2(rot_matrix[1], rot_matrix[0]);
    }
  else {
      beta = acos((std::min)((std::max)(rot_matrix[8],-1.0),1.0));
      alpha = atan2(rot_matrix[5], -rot_matrix[2]);
      gamma = atan2(rot_matrix[7], rot_matrix[6]);
    }
  if(!err)
    return;

  *err = 0;
  const double c[3] = {cos(alpha), cos(beta), cos(gamma)};
  const double s[3] = {sin(alpha), sin(beta), sin(gamma)};

  *err += fabs(c[0]*c[1]*c[2]-s[0]*s[2]-rot_matrix[0]);
  *err += fabs(c[0]*s[2]+c[2]*c[1]*s[0]-rot_matrix[1]);
  *err += fabs(-c[2]*s[1]-rot_matrix[2]);


  *err += fabs(-c[1]*s[2]*c[0]-c[2]*s[0]-rot_matrix[3]);
  *err += fabs(c[0]*c[2]-c[1]*s[0]*s[2]-rot_matrix[4]);
  *err += fabs(s[2]*s[1]-rot_matrix[5]);

  *err += fabs(c[0]*s[1]-rot_matrix[6]);
  *err += fabs(s[1]*s[0]-rot_matrix[7]);
  *err += fabs(c[1]-rot_matrix[8]);
}

void rot_n_2_z_by_zyz(const double *n, double *abc)
{
  abc[0] = -atan2(n[1], n[0]);
  abc[1] = -acos((std::min)((std::max)(n[2],-1.0),1.0));
  abc[2] = 0;
#if 0 // check
  matrixd Rnz(3, 3);
  zyz_angle_2_rotation_matrix(
        abc[0], abc[1], abc[2], &Rnz[0]);
  itr_matrix<const double *> n1(3, 1, n);
  matrixd z = Rnz*n1;
  z[2] -= 1;
  if(norm(z) > 1e-5)
    cout << "Rnz: " << Rnz*n1 << endl;
#endif

}


void zyz2frame(const zjucad::matrix::matrix<double> & zyz,
               zjucad::matrix::matrix<zjucad::matrix::matrix<double> > & frame)
{
  if(frame.size() != zyz.size(2)) frame.resize(zyz.size(2),1);
  for(size_t ti = 0; ti < zyz.size(2); ++ti){
      if(frame[ti].size(1) != 3 || frame[ti].size(2) != 3)
        frame[ti].resize(3,3);
      zyz_angle_2_rotation_matrix1(&zyz(0,ti), &frame[ti][0]);
    }
}

void frame2zyz(const zjucad::matrix::matrix<zjucad::matrix::matrix<double> > & frame,
               zjucad::matrix::matrix<double> & zyz)
{
  if(zyz.size(1) != 3 || zyz.size(2) != frame.size())
    zyz.resize(3,frame.size());
  for(size_t ti = 0; ti < frame.size(); ++ti){
      rotation_matrix_2_zyz_angle(&frame[ti][0], &zyz(0,ti), 0);
    }
}
