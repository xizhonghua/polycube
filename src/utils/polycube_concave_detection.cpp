#include <numeric>
#include <fstream>
#include "../common/vtk.h"
#include "../tetmesh/hex_io.h"
#include <jtflib/mesh/mesh.h>
#include <jtflib/mesh/io.h>

using namespace std;
using namespace zjucad::matrix;

int polycube_concave_detection(int argc, char * argv[])
{
  if(argc != 3){
    cerr << "# [usage] polycube_concave_detection orig_tet polycube_tet" << endl;
    return __LINE__;
  }

  jtf::mesh::meshes tm_orig, tm_polycube;
  if(jtf::mesh::tet_mesh_read_from_zjumat(argv[1], &tm_orig.node_, &tm_orig.mesh_))
    return __LINE__;

  if(jtf::mesh::tet_mesh_read_from_zjumat(argv[2], &tm_polycube.node_, &tm_polycube.mesh_))
    return __LINE__;

  if(norm(tm_orig.mesh_ - tm_polycube.mesh_) > 1e-5) {
    cerr << "# [error] input orig_tet is not the same as polycube in topolgy." << endl;
    return __LINE__;
  }

  std::unique_ptr<jtf::mesh::face2tet_adjacent> fa(jtf::mesh::face2tet_adjacent::create(tm_orig.mesh_));
  if(!fa.get()){
    cerr << "# [error] can not buildjtf::mesh::face2tet_adjacent." << endl;
    return __LINE__;
  }

  matrix<size_t> outside_face;
  get_outside_face(*fa, outside_face);

  unique_ptr<jtf::mesh::edge2cell_adjacent> ea(
        jtf::mesh::edge2cell_adjacent::create(outside_face));
  if(!ea.get()){
    cerr << "# [error] can not build edge2cell_adjacent." << endl;
    return __LINE__;
  }

  // for each edge, check wether it has the same concave state between original tet
  // and polycube tet
  vector<size_t> concave_strange_edgs;
  matrix<double> e_orig(3,3), e_polycube(3,3);
  for(size_t ei = 0; ei < ea->edges_.size(); ++ei){
    const pair<size_t,size_t> & one_edge = ea->edges_[ei];
    if((one_edge.first == 123 && one_edge.second == 1214) ||
       (one_edge.first == 1214 && one_edge.second == 123))
      cerr << endl;
    const pair<size_t,size_t> & tri_pair = ea->edge2cell_[ei];
    if(ea->is_boundary_edge(tri_pair)) continue;
    const size_t left_point =
        std::accumulate(outside_face(colon(), tri_pair.first).begin(),
                        outside_face(colon(), tri_pair.first).end(),
                        static_cast<size_t>(0))
        - one_edge.first - one_edge.second;

    const size_t right_point =
        std::accumulate(outside_face(colon(), tri_pair.second).begin(),
                        outside_face(colon(), tri_pair.second).end(),
                        static_cast<size_t>(0))
        - one_edge.first - one_edge.second;


    e_orig(colon(),0) = tm_orig.node_(colon(),one_edge.second)
                        - tm_orig.node_(colon(),one_edge.first);
    e_orig(colon(),1) = tm_orig.node_(colon(), left_point)
                        - tm_orig.node_(colon(), one_edge.first);
    e_orig(colon(),2) = tm_orig.node_(colon(), right_point)
                        -tm_orig.node_(colon(), one_edge.first);

    e_polycube(colon(),0) = tm_polycube.node_(colon(),one_edge.second)
                            - tm_polycube.node_(colon(),one_edge.first);
    e_polycube(colon(),1) = tm_polycube.node_(colon(), left_point)
                            - tm_polycube.node_(colon(), one_edge.first);
    e_polycube(colon(),2) = tm_polycube.node_(colon(), right_point)
                            -tm_polycube.node_(colon(), one_edge.first);

    matrix<double> N_left_orig = cross(e_orig(colon(),0), e_orig(colon(),1));

    double N_left_len = norm(N_left_orig);
    if(N_left_len < 1e-6) N_left_len = 1.0;
    N_left_orig /= N_left_len;

    matrix<double> N_right_orig = cross(e_orig(colon(),2), e_orig(colon(),0));

    double N_right_len = norm(N_right_orig);
    if(N_right_len < 1e-6) N_right_len = 1.0;
    N_right_orig /= N_right_len;

    double e_orig_len = norm(e_orig(colon(),0));
    if(e_orig_len <1e-6) e_orig_len = 1.0;
    e_orig(colon(),0) /= e_orig_len;

    const double orig_val = dot(cross(N_left_orig, N_right_orig), e_orig(colon(),0));

    matrix<double> N_left_polycube = cross(e_polycube(colon(),0), e_polycube(colon(),1));

    double N_left_p = norm(N_left_polycube);
    if(N_left_p < 1e-6) N_left_p = 1.0;
    N_left_polycube /= N_left_p;

    matrix<double> N_right_polycube = cross(e_polycube(colon(),2), e_polycube(colon(),0));

    double N_right_p = norm(N_right_polycube);
    if(N_right_p < 1e-6) N_right_p = 1.0;
    N_right_polycube /= N_right_p;

    double e_polycube_len = norm(e_polycube(colon(),0));
    if(e_polycube_len < 1e-6) e_polycube_len = 1.0;
    e_polycube(colon(),0) /= e_polycube_len;

    const double polycube_val = dot(cross(N_left_polycube, N_right_polycube), e_polycube(colon(),0));
    if(orig_val * polycube_val < 0 && fabs(orig_val * polycube_val) > 1e-3){ // concave_changed
      concave_strange_edgs.push_back(one_edge.first);
      concave_strange_edgs.push_back(one_edge.second);
    }
  }


  cerr << "# [info] find concave_strange_edge number " << concave_strange_edgs.size()/2 << endl;

  {
    ofstream ofs_orig("concave_strange_edge_in_orig.vtk");
    ofstream ofs_polycube("concave_strange_edge_in_polycube.vtk");

    line2vtk(ofs_orig, &tm_orig.node_[0], tm_orig.node_.size(2),
             &concave_strange_edgs[0], concave_strange_edgs.size()/2);

    line2vtk(ofs_polycube, &tm_polycube.node_[0], tm_polycube.node_.size(2),
             &concave_strange_edgs[0], concave_strange_edgs.size()/2);
  }
  return 0;

}
