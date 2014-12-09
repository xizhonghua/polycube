#include <iostream>
#include <fstream>
#include <jtflib/mesh/io.h>
#include <jtflib/mesh/mesh.h>
#include <jtflib/mesh/util.h>
#include "../tetmesh/util.h"

using namespace std;
using namespace zjucad::matrix;

int fv2nsys(int argc, char *argv[]){
  if(argc != 2){
    cerr << "#  [usage] fv2nsys fv" << endl;
    return __LINE__;
  }
  ifstream ifs(argv[1]);
  if(ifs.fail()){
    cerr << "# [error] can not open fv file." << endl;
    return __LINE__;
  }
  size_t face_num = 0;
  ifs >> face_num;
  cerr << face_num << " " << 4 << endl;
  vector<double> dir(6);
  for(size_t fi = 0; fi < face_num; ++fi){
    for(size_t i = 0; i < 6; ++i)
      ifs >> dir[i];
    cerr << dir[0] << " " << dir[1] << " " << dir[2] << endl;
  }
  return 0;
}

int nsys2fv(int argc, char * argv[])
{
  if(argc != 3){
    cerr << "# [usage] nsys2fv obj nsys" << endl;
    return __LINE__;
  }

  jtf::mesh::meshes trim;

  if(jtf::mesh::load_obj(argv[1], trim.mesh_, trim.node_)){
    cerr << "# [error] can not open obj file." << endl;
    return __LINE__;
  }

  ifstream ifs(argv[2]);
  if(ifs.fail()){
    cerr << "# [error] can not open nsys file." << endl;
    return __LINE__;
  }

  size_t face_num = 0;
  ifs >> face_num;
  if(face_num != trim.mesh_.size(2)){
    cerr << "# [error] nsys face num does not match obj" << endl;
    return __LINE__;
  }
  size_t nsys_type = -1;
  ifs >> nsys_type;
  if(nsys_type != 4){
    cerr << "# [error] only 4-sys frame field is supported." << endl;
    return __LINE__;
  }

  matrixd frame_field(3, face_num);
  for(size_t fi = 0; fi < face_num; ++fi){
    ifs >> frame_field(0,fi) >> frame_field(1,fi) >> frame_field(2, fi);
  }

  matrixd face_normal;
  jtf::mesh::cal_face_normal(trim.mesh_,trim.node_, face_normal);

  cout << face_num << endl;
  matrixd other_dir = zeros<double>(3,1);
  for(size_t fi = 0; fi < face_num; ++fi){
    cout << frame_field(0, fi) << " " << frame_field(1, fi)
         << " " << frame_field(2, fi) << endl;
    other_dir = cross(frame_field(colon(),fi), face_normal(colon(),fi));
    const double len = norm(other_dir);
    if(len > 1e-6)
      other_dir /= len;
    cout << other_dir[0] << " " << other_dir[1] << other_dir[2] << endl;
    cout << endl;
  }

  return 0;
}
