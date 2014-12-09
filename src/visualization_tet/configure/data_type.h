/**
 *@file type_data.h
 *@brief this is file support some base function and struct
 *
 *@author Li Quan
 *@data 2012-11-20
 *@version1.0
 */

#ifndef DATA_TYPE_H
#define DATA_TYPE_H
#include <zjucad/matrix/matrix.h>
#include <zjucad/matrix/io.h>
#include <cmath>
#include <cstring>
#include <string.h>

namespace lq {
 
   
  const unsigned char red[3] = {255, 0, 0};
  const unsigned char blue[3] = {0, 0, 255};
  const unsigned char white[3] = {255, 255,255};
  const unsigned char COLOR[3][3] = { {0 , 0, 255},
                                      {255, 255, 255},
                                      {255, 0, 0}};
  const size_t RED = 2;
  const size_t BLUE = 0;
  const size_t WHITE = 1;
  const double ACCURACY = 1e-5;
  enum WIN{LEFT = 1, RIGHT};
  using namespace std;
  typedef zjucad::matrix::matrix<double> matrixd;
  typedef zjucad::matrix::matrix<float> matrixf;
  typedef zjucad::matrix::matrix<size_t> matrixst;
  struct vertex_normal {
    double nx, ny, nz;
    vertex_normal() {}
   vertex_normal(double nx0, double ny0, double nz0):
    nx(nx0), ny(ny0), nz(nz0) {}
  };
  
  struct point {
    double x, y, z;
    point() {}
   point(double x0, double y0, double z0):
    x(x0), y(y0), z(z0) {}
   point(const point &p):x(p.x), y(p.y), z(p.z){}
   point(const double p[3]):x(p[0]), y(p[1]), z(p[2]){}
   point(const vertex_normal &normal):x(normal.nx), y(normal.ny), z(normal.nz){}
    
  };

 

  struct triangle {
    std::vector<int> vertex;
    std::vector<int> texture;
    std::vector<int> point_normal;
    vertex_normal triangle_normal;
  };
  template <typename T>
  inline bool equal(const T &x, const T &y)
  {
    return (fabs(x - y) < ACCURACY);
  }

  template <typename T>
  inline bool is_zero(const T &x)
  {
    return (fabs(x) < ACCURACY);
  }    
  
  template <typename T>
  inline bool is_less_equal(const T &x, const T &y)
  {
    return ((equal(x, y)) || (x < y));
  }

  template <typename T>
  inline bool is_larger_equal(const T &x, const T &y)
  {
    return ((equal(x, y)) || (x > y));
  }

  //get area of a triangle 
  template <typename T>
  inline double tri_area(const T &p1, const T &p2, const T &p3)
  {
    return norm(cross((p1 - p2), (p3 - p2))) / 2;
  }

 
}


#endif
