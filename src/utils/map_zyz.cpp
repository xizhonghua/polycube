#include <fstream>
#include <iostream>
#include <hjlib/math/polar.h>
#include "../tetmesh/tetmesh.h"
#include "../common/zyz.h"
#include "../numeric/util.h"


using namespace std;
using namespace zjucad::matrix;

int map_zyz(int argc, char * argv[])
{
  if(argc != 5){
      cerr << "# [usage] map_zyz tet_0 zyz_0 tet_1 zyz1" << endl;
      return __LINE__;
    }

  jtf::tet_mesh t0(argv[1]);
  jtf::tet_mesh t1(argv[3]);

  if(norm(t0.tetmesh_.mesh_ - t1.tetmesh_.mesh_) > 1e-6){
      cerr << "# [error] mesh is not compatible." << endl;
      return __LINE__;
    }

  matrix<double> zyz;
  if(jtf::mesh::read_matrix(argv[2], zyz)){
      cerr << "# [error] can not open zyz file." << endl;
      return __LINE__;
    }
  if(zyz.size(2) != t0.tetmesh_.mesh_.size(2)){
      cerr << "# [error] zyz file is not compatible with tet." << endl;
      return __LINE__;
    }

  matrix<matrix<double> > rot0(zyz.size(2),1);
  matrix<matrix<double> > rot1(zyz.size(2),1);
  matrix<matrix<double> > deform_gradient;
  cal_deformation_gradient(t0.tetmesh_.mesh_, t1.tetmesh_.node_, t0.tetmesh_.node_, deform_gradient);

  matrix<double> R;
  for(size_t ti = 0; ti < t0.tetmesh_.mesh_.size(2); ++ti){
      rot0[ti].resize(3,3);
      rot1[ti].resize(3,3);
      zyz_angle_2_rotation_matrix1(&zyz(0,ti), &rot0[ti][0]);
      R = deform_gradient[ti];
      hj::polar3d p;
      p(R,2);
      rot1[ti] = R * rot0[ti];
    }

  matrix<double> zyz1(3,t1.tetmesh_.mesh_.size(2));
  for(size_t ti = 0; ti < t1.tetmesh_.mesh_.size(2); ++ti){
      rotation_matrix_2_zyz_angle(&rot1[ti][0], &zyz1(0,ti),0);
    }

  jtf::mesh::write_matrix(argv[4], zyz1);

  return 0;
}
