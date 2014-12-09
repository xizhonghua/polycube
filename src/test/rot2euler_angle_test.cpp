#include "rot2euler_angle_test.h"
#include <zjucad/matrix/matrix.h>
#include <zjucad/matrix/io.h>
#include <hjlib/math/quaternion.h>
#include "../common/transition.h"
#include "../common/transition_type.h"

using namespace std;
using namespace zjucad::matrix;
using namespace hj::quaternion;

void rot2euler_angle_test::convert_rot_to_euler_zyx_test()
{
  //matrixd A = rand<double>(3,3) + 3.1415926535897/4;
  for(size_t t = 0; t < 0; ++t){
    cerr << "# [info] test " << t << endl;
    matrixd A = type_transition2(t) /*+ rand<double>(3,3) * 0.1*/;
    matrixd angle_zyx = zeros<double>(3,2);

    //matrixd quartenion = zeros<double>(4,1);
    //m332quat<double>(A,quartenion);
    //  matrixd angle = zeros<double>(3,1);
    //    convert_quat_to_euler_angle(&quartenion[0],&angle[0]);

    convert_rot_to_euler_zyx(&A[0],&angle_zyx[0]);

    //cerr << "# [info] angle " << angle_zyx << endl;
    matrixd B0 = zeros<double>(3,3);
    convert_euler_angle_to_rot(&angle_zyx(0,0),&B0[0]);
    cerr << "# [info] matrix A " << A << endl;
    cerr << "# [info] matrix B0 " << B0 << endl;

    CPPUNIT_ASSERT_MESSAGE("Two matrix do not match", fabs(norm(A - B0)) < 1e-8);
  }

  matrixd rot = zeros<double>(3,3);
  rot(1,0) = -1;
  rot(2,1) = 1;
  rot(0,2) = -1;
  cerr << "# [info] rot type " << type_transition1(rot) << endl;
  matrixd zyx = zeros<double>(3,2);
  convert_rot_to_euler_zyx(&rot[0],&zyx[0]);
  cerr << "# [info] zyx " << zyx << endl;
  matrixd zyx_m = zeros<double>(3,1);
  zyx_m[2] = zyx(2,0) * -1;
  matrixd new_ = zeros<double>(3,3);
  convert_euler_angle_to_rot(&zyx_m[0],&new_[0]);
  cerr << "# [info] rot type " << type_transition1(new_) << endl;
  cerr << "# [info] new_ " << new_ << endl;
}

void rot2euler_angle_test::delete_axis_test()
{
  matrixd rot = eye<double>(3);
  rot = temp(rot * type_transition2(0));
  rot = temp(rot * type_transition2(3));
  matrixd angle_zyx = zeros<double>(3,2);
  convert_rot_to_euler_zyx(&rot[0],&angle_zyx[0]);
  matrixd angle_z = zeros<double>(3,1);
  angle_z[0] = angle_zyx(0,0) * -1;
  matrixd angle_x = zeros<double>(3,1);
  angle_x[2] = angle_zyx(2,0) * -1;
  matrixd delete_z = eye<double>(3), delete_x = eye<double>(3);
  convert_euler_angle_to_rot(&angle_z[0],&delete_z[0]);
  cerr << "# delete angle_z type " << type_transition1(delete_z) ;
  convert_euler_angle_to_rot(&angle_x[0],&delete_x[0]);
  cerr << "# [info] rot = = " << delete_z * rot << endl;
  cerr << "# [info] type: " << type_transition1(delete_z * rot) << endl;
  cerr << "# [info] type: " << type_transition1( rot * delete_x) << endl;
}

