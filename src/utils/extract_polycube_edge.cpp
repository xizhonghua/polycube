#include <fstream>

#include <jtflib/mesh/mesh.h>
#include <jtflib/mesh/util.h>
#include <jtflib/mesh/io.h>

#include "../common/vtk.h"
#include "../tetmesh/tetmesh.h"
#include "../tetmesh/util.h"

using namespace std;
using namespace zjucad::matrix;
using namespace jtf::mesh;

int extract_polycube_edge(int argc, char * argv[])
{
  if(argc != 2){
    cerr << "# [usage] extract_polycube_edge polycube_tet." << endl;
    return __LINE__;
  }

  jtf::mesh::meshes tm;
  if(jtf::mesh::tet_mesh_read_from_zjumat(argv[1], &tm.node_, &tm.mesh_))
    return __LINE__;

  unique_ptr<jtf::mesh::face2tet_adjacent> fa(jtf::mesh::face2tet_adjacent::create(tm.mesh_));
  if(!fa.get()){
    cerr << "# [error] can not buildjtf::mesh::face2tet_adjacent." << endl;
    return __LINE__;
  }

  matrix<size_t> outside_face;
  get_outside_face(*fa, outside_face);

  matrix<double> face_normal;
  jtf::mesh::cal_face_normal(outside_face, tm.node_, face_normal);


  unique_ptr<edge2cell_adjacent> ea(edge2cell_adjacent::create(outside_face));
  if(!ea.get()){
    cerr << "# [error] can not build edge2cell_adjacent." << endl;
    return __LINE__;
  }

  vector<size_t> edges;
  for(size_t ei = 0; ei < ea->edges_.size(); ++ei){
    const pair<size_t,size_t> & tri_pair = ea->edge2cell_[ei];
    if(ea->is_boundary_edge(tri_pair)) continue;
    const double dot_val = dot(face_normal(colon(),tri_pair.first),
                               face_normal(colon(),tri_pair.second));
    if(fabs(dot_val) < 0.3) {
      edges.push_back(ea->edges_[ei].first);
      edges.push_back(ea->edges_[ei].second);
    }
  }

  ofstream ofs("polycube_edge.vtk");
  line2vtk(ofs, &tm.node_[0], tm.node_.size(2), &edges[0], edges.size()/2);
  return 0;
}
