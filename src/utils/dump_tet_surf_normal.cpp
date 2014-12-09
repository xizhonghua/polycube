#include <iostream>
#include "../tetmesh/util.h"
#include "../tetmesh/hex_io.h"
#include "../common/IO.h"
#include "../common/zyz.h"
#include "../common/def.h"
#include "../common/vtk.h"
#include <jtflib/mesh/util.h>

using namespace std;
using namespace zjucad::matrix;

int dump_tet_surf_normal(int argc, char * argv[])
{
  if(argc != 3){
    cerr << "# [usage] dump_tet_surf_normal tet normal_vtk" << endl;
    return __LINE__;
  }

  matrixst tet;
  matrixd node;
  if(jtf::mesh::tet_mesh_read_from_zjumat(argv[1], &node, & tet)){
    cerr << "# [error] can not read tet file." << endl;
    return __LINE__;
  }

  unique_ptr<jtf::mesh::face2tet_adjacent> fa(jtf::mesh::face2tet_adjacent::create(tet));
  if(!fa.get()){
    cerr << "# [error] can not buildjtf::mesh::face2tet_adjacent" << endl;
    return __LINE__;
  }

  matrixst outside_face_idx, outside_face;
  get_outside_face_idx(*fa, outside_face_idx);
  get_outside_face(*fa, outside_face);
  matrixd face_normal;

  jtf::mesh::cal_face_normal(outside_face, node, face_normal);
  jtf::tetmesh::orient_face_normal_outside_tetmesh(
        tet, node, outside_face, outside_face_idx, *fa, face_normal);

  const double av_len = jtf::mesh::cal_average_edge(outside_face, node);

  matrixd new_node(3, outside_face.size(2) * 2);
  matrixd center;
  vector<size_t> edges;
  edges.reserve(2 * outside_face.size(2));
  for(size_t fi = 0; fi < outside_face.size(2); ++fi){
    center = zeros<double>(3,1);
    for(size_t pi = 0; pi < outside_face.size(1); ++pi){
      center += node(colon(), outside_face(pi, fi));
    }
    center /= outside_face.size(1);
    new_node(colon(), 2 * fi) = center;
    new_node(colon(), 2 * fi + 1) = center + face_normal(colon(), fi) * av_len;
    edges.push_back(2 * fi);
    edges.push_back(2 * fi+1);
  }

  ofstream ofs(argv[2]);
  if(ofs.fail()){
    cerr << "# [error] can not open normal_vtk file" << endl;
    return __LINE__;
  }
  line2vtk(ofs, &new_node[0], new_node.size(2), &edges[0], edges.size()/2);
  return 0;
}
