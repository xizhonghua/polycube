#include <fstream>

#include <jtflib/mesh/io.h>
#include <jtflib/mesh/mesh.h>
#include <jtflib/mesh/util.h>

#include <zjucad/matrix/matrix.h>

#include "../common/vtk.h"

using namespace std;
using namespace zjucad::matrix;

int normal_map(int argc,char * argv[])
{
  if(argc != 4){
      cerr << "# [usage] normal_map obj_a obj_a output_vtk" << endl;
      return __LINE__;
    }
  using namespace jtf::mesh;
  meshes a,b;

  if(load_obj(argv[1], a.mesh_, a.node_))
    return __LINE__;
  if(load_obj(argv[2], b.mesh_, b.node_))
    return __LINE__;

  if(a.mesh_.size(2) != b.mesh_.size(2)){
      cerr << "# [error] incompatiable input obj mesh." << endl;
      return __LINE__;
    }

  matrix<double> face_normal(3, a.mesh_.size(2));

  cal_face_normal(a.mesh_, a.node_, face_normal);

  {
    ofstream ofs(argv[3]);
    if(ofs.fail()){
        cerr << "# [error] can not open file." << endl;
        return __LINE__;
      }

    matrix<double> color = ones<double>(4, a.mesh_.size(2));
    color(colon(0,2),colon()) = face_normal;

    tri2vtk(ofs, &b.node_[0], b.node_.size(2), &b.mesh_[0], b.mesh_.size(2));
    cell_data_rgba(ofs, &color[0], color.size(2), "normal_color");
    matrix<double> x = face_normal(0,colon()),
        y = face_normal(1,colon()),
        z = face_normal(2,colon());
    vtk_data(ofs, &x[0], x.size(), "normal_x", "normal_x");
    vtk_data(ofs, &y[0], y.size(), "normal_y", "normal_y");
    vtk_data(ofs, &z[0], z.size(), "normal_z", "normal_z");
  }

  {
    string point_normal_str = argv[3];
    point_normal_str += ".point.vtk";

    ofstream ofs(point_normal_str.c_str());
    if(ofs.fail()){
        cerr << "# [error] can not open file." << endl;
        return __LINE__;
      }

    matrix<double> point_normal;
    jtf::mesh::cal_point_normal(a.mesh_, a.node_, point_normal);

    matrix<double> color = ones<double>(4,point_normal.size(2));
    color(colon(0,2),colon()) = point_normal;

    tri2vtk(ofs, &b.node_[0], b.node_.size(2), &b.mesh_[0], b.mesh_.size(2));
    matrix<double> x = point_normal(0,colon()),
        y = point_normal(1,colon()),
        z = point_normal(2,colon());
    matrix<double> xyz = point_normal(0,colon());
    for(size_t i = 0; i < xyz.size(); ++i)
        xyz[i] = x[i] * y[i] *z[i];
    point_data(ofs, &x[0], x.size(), "normal_x", "normal_x");
    vtk_data(ofs, &y[0], y.size(), "normal_y", "normal_y");
    vtk_data(ofs, &z[0], z.size(), "normal_z", "normal_z");
    vtk_data(ofs, &xyz[0], xyz.size(), "normal_xyz", "normal_xyz");
    vtk_data_rgba(ofs, &color[0], color.size(2), "normal_color", "normal_color");
  }

  return 0;
}
