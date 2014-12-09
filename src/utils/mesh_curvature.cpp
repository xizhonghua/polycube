#include <jtflib/mesh/mesh.h>
#include <jtflib/math/math.h>
#include <jtflib/mesh/util.h>
#include <jtflib/mesh/io.h>
#include <fstream>
#include "../common/vtk.h"

using namespace std;
using namespace zjucad::matrix;

int mesh_curvature(int argc, char * argv[])
{
  if(argc != 2){
      cerr << "# [usage] mesh_curvature obj" << endl;
      return __LINE__;
    }

  jtf::mesh::meshes tri;
  if(jtf::mesh::load_obj(argv[1], tri.mesh_, tri.node_)){
      return __LINE__;
    }

  jtf::mesh::curvature cur(tri.mesh_, tri.node_);
  cur.generate();

  matrix<double> gauss_curvature = cur.get_gauss_curvature();
  matrix<double> mean_curvature  = cur.get_mean_curvature();
  matrix<double> k1k2 = cur.get_principle_curvature();
  matrix<double> kmax = k1k2(0,colon());
  matrix<double> kmin = k1k2(1,colon());
  matrix<double> d1d2 = cur.get_principle_curvature_direction();

  //  {
  //    ofstream ofs("gauss_curvature.vtk");
  //    tri2vtk(ofs, &tri.node_[0], tri.node_.size(2), &tri.mesh_[0], tri.mesh_.size(2));
  //    point_data(ofs, &gauss_curvature[0], gauss_curvature.size(), "gauss");
  //  }
  //  {
  //    ofstream ofs("mean_curvature.vtk");
  //    tri2vtk(ofs, &tri.node_[0], tri.node_.size(2), &tri.mesh_[0], tri.mesh_.size(2));
  //    point_data(ofs, &mean_curvature[0], mean_curvature.size(), "mean");
  //  }
    {
      ofstream ofs("kmax_curvature.vtk");
      tri2vtk(ofs, &tri.node_[0], tri.node_.size(2), &tri.mesh_[0], tri.mesh_.size(2));
      point_data(ofs, &kmax[0], kmax.size(), "kmax");
    }
    {
      ofstream ofs("kmin_curvature.vtk");
      tri2vtk(ofs, &tri.node_[0], tri.node_.size(2), &tri.mesh_[0], tri.mesh_.size(2));
      point_data(ofs, &kmin[0], kmin.size(), "kmin");
    }
  const double delta = 0.02;
  {
    ofstream ofs("d1.vtk");
    matrix<size_t> edges(2, 2*tri.node_.size(2));
    matrix<double> new_node(3, tri.node_.size(2)*3);
    for(size_t i = 0; i < tri.node_.size(2); ++i){
        new_node(colon(),i) = tri.node_(colon(),i);
        new_node(colon(),i+tri.node_.size(2)) = tri.node_(colon(),i) + delta * d1d2(colon(0,2),i);
        new_node(colon(),i+2*tri.node_.size(2)) = tri.node_(colon(),i) - delta * d1d2(colon(0,2),i);
        edges(0,2*i) = i;
        edges(1,2*i) = i+tri.node_.size(2);
        edges(0,2*i+1) = i;
        edges(1,2*i+1) = i+2*tri.node_.size(2);
      }
    line2vtk(ofs, &new_node[0], new_node.size(2), &edges[0], edges.size(2));
  }
  {
    ofstream ofs("d2.vtk");
    matrix<size_t> edges(2, 2*tri.node_.size(2));
    matrix<double> new_node(3, tri.node_.size(2)*3);
    for(size_t i = 0; i < tri.node_.size(2); ++i){
        new_node(colon(),i) = tri.node_(colon(),i);
        new_node(colon(),i+tri.node_.size(2)) = tri.node_(colon(),i) + delta * d1d2(colon(3,5),i);
        new_node(colon(),i+2*tri.node_.size(2)) = tri.node_(colon(),i) - delta * d1d2(colon(3,5),i);
        edges(0,2*i) = i;
        edges(1,2*i) = i+tri.node_.size(2);
        edges(0,2*i+1) = i;
        edges(1,2*i+1) = i+2*tri.node_.size(2);
      }
    line2vtk(ofs, &new_node[0], new_node.size(2), &edges[0], edges.size(2));
  }
  return 0;
}
