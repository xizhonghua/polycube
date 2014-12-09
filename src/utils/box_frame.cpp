#include <fstream>
#include <iostream>
#include <vector>
#include <algorithm>

#include "../common/vtk.h"
#include "../common/IO.h"
#include "../common/zyz.h"

#include "../tetmesh/tetmesh.h"

#include <zjucad/matrix/io.h>

using namespace std;
using namespace zjucad::matrix;

int box_frame(const matrixd &pts, const matrixd &F_9xn, double size)
{
  matrixd box_pts(3, 8);
  for(int z = 0; z < 2; ++z) {
      for(int y = 0; y < 2; ++y) {
          for(int x = 0; x < 2; ++x) {
              box_pts(0, z*4+y*2+x) = (x-0.5)*size;
              box_pts(1, z*4+y*2+x) = (y-0.5)*size;
              box_pts(2, z*4+y*2+x) = (z-0.5)*size;
            }
        }
    }

  // face_in_hex is ordered by normal x, -x, y, -y, z, -z.
  static const int face_in_hex[][4] = {
    7,5,1,3,
    0,4,6,2,
    7,3,2,6,
    0,1,5,4,
    7,6,4,5,
    0,2,3,1,
  };
  static const char *mat_name[] = {
    "+x", "-x", "+y", "-y", "+z", "-z"
  };

  cout << "mtllib box-frame-mat.mtl\n";
  matrixd frame(3, 3), corners(3, 8);
  char box_name[256];
  for(size_t pi = 0; pi < pts.size(2); ++pi) {
      sprintf(box_name, "box-%04ld", pi);

      frame(colon()) = F_9xn(colon(), pi);
      corners = frame*box_pts+pts(colon(), pi)*ones<double>(1, 8);
      for(int ci = 0; ci < 8; ++ci) {
          cout << "v " << corners(0, ci)
               << ' ' << corners(1, ci)
               << ' ' << corners(2, ci) << '\n';
        }
      cout << "g " << box_name << '\n';
      for(int fi = 0; fi < 6; ++fi) {
          cout << "usemtl box-frame-mat" << mat_name[fi] << '\n';
          cout << 'f';
          for(int ni = 0; ni < 4; ++ni)
            cout << ' ' << face_in_hex[fi][ni]+pi*8+1;
          cout << '\n';
        }

    }
  return 0;
}

int box_frame(int argc, char *argv[])
{
  if(argc < 4) {
      cerr << "box_frame point frame size" << endl;
      return __LINE__;
    }

  matrixd pts;
  matrixst mesh;
  if(jtf::mesh::load_obj(argv[1], mesh, pts)){
      cerr << "load_points_from_obj fail." << endl;
      return __LINE__;
    }

  matrixd frame_9xn(9, pts.size(2)); {
    ifstream ifs(argv[2], ifstream::binary);
    if(ifs.fail()) {
        cerr << "open " << argv[2] << " fail." << endl;
        return __LINE__;
      }
    matrixd zyz;
    jtf::mesh::read_matrix(ifs, zyz);
    if(frame_9xn.size(2) != zyz.size(2)) {
        cerr << "# incompatible frame and pts." << endl;
        return __LINE__;
      }

    for(size_t pi = 0; pi < zyz.size(2); ++pi)
      zyz_angle_2_rotation_matrix1(&zyz(0, pi), &frame_9xn(0, pi));
  }

  return box_frame(pts, frame_9xn, atof(argv[3]));
}
