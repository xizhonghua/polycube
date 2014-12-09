#include <fstream>
#include <iostream>
#include <vector>
#include <algorithm>

using namespace std;

#include "../common/zyz.h"
#include "../common/IO.h"
#include "../spherical_harmonics/rot_cubic_f_SH.h"

#include <zjucad/matrix/io.h>
using namespace zjucad::matrix;

template <typename T>
inline void HSV_to_RGB(T h, T s, T v, T *r, T *g, T *b)
{
  // H is given on [0, 6] or UNDEFINED. S and V are given on [0, 1].
  // RGB are each returned on [0, 1].
  T m, n, f;
  int i;

  if (h < 0 || h > 6) {
      *r = *g = *b = v;
      return;
      //RETURN_RGB(v, v, v);
    }
  i = floor(h);
  f = h - i;
  if ( !(i&1) ) f = 1 - f; // if i is even
  m = v * (1 - s);
  n = v * (1 - s * f);
  switch (i) {
    case 6:
    case 0:
      *r = v; *g = n; *b = m; return; //RETURN_RGB(v, n, m);
    case 1:
      *r = n; *g = v; *b = m; return; //RETURN_RGB(n, v, m);
    case 2:
      *r = m; *g = v; *b = n; return; //RETURN_RGB(m, v, n)
    case 3:
      *r = m; *g = n; *b = v; return; //RETURN_RGB(m, n, v);
    case 4:
      *r = n; *g = m; *b = v; return; //RETURN_RGB(n, m, v);
    case 5:
      *r = v; *g = m; *b = n; return; //RETURN_RGB(v, m, n);
    }
} 

// err must be in [0, 1]
void scalar2rgb(const matrixd &err, matrixd &rgb)
{
  assert(err.size() == rgb.size(2));
  for(size_t i = 0; i < err.size(); ++i) {
      HSV_to_RGB((1-err[i])*4, 1.0, 1.0,
                 &rgb(0, i), &rgb(1, i), &rgb(2, i));
    }
}

int zyz2color(int argc, char *argv[])
{
  if(argc < 2) {
      cerr << "zyz2color zyz." << endl;
      return __LINE__;
    }

  string zyz_file = argv[1], color_file = zyz_file+".color.mat";

  matrixd zyz;
  ifstream ifs(zyz_file.c_str(), ifstream::binary);
  jtf::mesh::read_matrix(ifs, zyz);

  const size_t n = zyz.size()/3;
  {
    matrixd t = zyz;
    zyz.resize(3, n);
    zyz(colon()) = t(colon());
  }

  matrixd err(zyz.size(2));
  matrixd sh(9), sh0(9), I=zeros<double>(3, 1);
  calc_rot_cubic_f_sh_(&sh0[0], &I[0]);
  for(size_t i = 0; i < zyz.size(2); ++i) {
      calc_rot_cubic_f_sh_(&sh[0], &zyz(0, i));
      err[i] = norm(sh-sh0);
    }
  err /= norm(sh0)*2;

  matrixd rgb(3, zyz.size(2));

  scalar2rgb(err, rgb);

  ofstream ofs(color_file.c_str(), ofstream::binary);
  jtf::mesh::write_matrix(ofs, rgb);

  return 0;
}


int zyz2frame(int argc, char * argv[])
{
  if(argc != 3){
      cerr << "# [usage] zyz2frame zyz output_frame" << endl;
      return __LINE__;
    }
  matrixd zyz;
  ifstream ifs(argv[1], ifstream::binary);
  jtf::mesh::read_matrix(ifs, zyz);

  zjucad::matrix::matrix<zjucad::matrix::matrix<double> > frame(zyz.size(2),1);
  for(size_t fi = 0; fi < zyz.size(2); ++fi){
      frame[fi].resize(3,3);
      zyz_angle_2_rotation_matrix1(&zyz(0,fi), &frame[fi][0]);
    }

  ofstream ofs(argv[2]);
  ofs << "frame " << zyz.size(2) << endl;
  for(size_t fi = 0; fi < frame.size(); ++fi){
      for(size_t i = 0; i < frame[fi].size(); ++i)
        ofs << frame[fi][i] << " ";
      ofs << endl;
    }
  return 0;
}
