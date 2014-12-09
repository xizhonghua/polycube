#include <jtflib/mesh/io.h>
#include <fstream>
#include <jtflib/util/util.h>
#include <jtflib/mesh/mesh.h>

#include "../hexmesh/util.h"
#include "../common/vtk.h"
#include "../common/util.h"
using namespace std;
using namespace zjucad::matrix;

int draw_quad_type(int argc, char * argv[])
{
  if(argc != 3)  {
    cerr << "# [usage] draw_quad_type quad output_vtk" << endl;
    return __LINE__;
  }

  jtf::mesh::meshes qm;
  //  if(jtf::mesh::load_obj(argv[1], qm.quad, qm.node_,4))
  //    return __LINE__;

  ifstream fake_ofs(argv[1]);
  size_t node_num, face_num;
  fake_ofs >> node_num;
  qm.node_.resize(3, node_num);
  for(size_t pi = 0; pi < node_num; ++pi)
    for(size_t di = 0; di < 3; ++di)
      fake_ofs >> qm.node_(di,pi);

  fake_ofs >> face_num;
  qm.mesh_.resize(4, face_num);
  size_t trash;
  for(size_t fi = 0; fi < face_num; ++fi){
    fake_ofs >> trash;
    for(size_t pi = 0; pi < 4; ++pi)
      fake_ofs >> qm.mesh_(pi, fi);
  }

  unique_ptr<jtf::mesh::edge2cell_adjacent> eq(
        jtf::mesh::edge2cell_adjacent::create(qm.mesh_));
  if(!eq.get()){
    cerr << "# [error] can not build edge2cell_adjacent." << endl;
    return __LINE__;
  }

  vector<pair<size_t,size_t> > edges;
  for(size_t ei = 0; ei < eq->edges_.size();++ei){
    const pair<size_t,size_t> &one_edge = eq->edges_[ei];
    const pair<size_t,size_t> &quad_pair = eq->edge2cell_[ei];
    if(eq->is_boundary_edge(quad_pair)) edges.push_back(one_edge);
  }

  vector<deque<pair<size_t,size_t> > > chains;
  jtf::util::extract_chain_from_edges(edges, chains);
  for(size_t ci = 0; ci < chains.size(); ++ci){
    cerr << "# chain " << ci << endl;
    const deque<pair<size_t,size_t> > & one_chain = chains[ci];
    for(size_t ei = 0; ei < one_chain.size(); ++ei)
      cerr << one_chain[ei].first << " ";
    cerr << one_chain.back().second << endl;
  }
  //  jtf::mesh::save_to_obj("fake_quad_from_vtk.obj", qm.node_, qm.quad);
  //  matrix<double> face_normal;
  //  jtf::hexmesh::cal_face_normal(qm.quad,qm.node_,face_normal);

  //  matrix<double> eye_ = eye<double>(3);
  //  matrix<size_t> face_type(qm.quad.size(2),1);
  //  vector<pair<double,size_t> > error(6);
  //  for(size_t fi = 0; fi < qm.quad.size(2); ++fi){
  //    for(size_t ai = 0; ai < 6; ++ai){
  //      error[ai] = make_pair(dot(face_normal(colon(),fi),
  //                                 (ai%2==0?1.0:-1.0)*eye_(colon(),ai/2)),
  //                            ai/2);
  //    }
  //    sort(error.begin(), error.end());
  //    face_type[fi] = error.back().second;
  //  }

  //  ofstream ofs(argv[2]);
  //  quad2vtk(ofs, &qm.node[0], qm.node_.size(2), &qm.quad[0], qm.quad.size(2));
  //  cell_data(ofs, &face_type[0], face_type.size(), "normal_dir");

  return 0;
}
