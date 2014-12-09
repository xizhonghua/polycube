#include "../tetmesh/tetmesh.h"
#include "../tetmesh/util.h"
#include "../tetmesh/hex_io.h"
#include <jtflib/mesh/mesh.h>

#include "../common/util.h"
#include "../common/vtk.h"
#include "../common/visualize_tool.h"
#include <fstream>
using namespace std;
using namespace zjucad::matrix;

int load_surface_restricted_type_static(
    const char * filename,
    boost::unordered_map<size_t,size_t> & result_surface_type);

int dump_surface_singularity_point(
    const vector<vector<size_t> > & singularity_point,
    const matrix<double> & node,
    const double radius,
    const char * sphere_file,
    const string &singularity_point_str);

int find_polycube_edge(int argc, char * argv[])
{
  if(argc != 4 && argc != 5){
    cerr << "# [usage] find_polycube_edge polycube_tet surface_type input_unit_obj [radius]" << endl;
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

  boost::unordered_map<size_t,size_t> surface_type;
  if(load_surface_restricted_type_static(argv[2], surface_type))
    return __LINE__;

  matrix<size_t> outside_face, outside_face_idx;
  get_outside_face(*fa, outside_face);
  get_outside_face_idx(*fa, outside_face_idx);

  unique_ptr<jtf::mesh::edge2cell_adjacent> ea(
        jtf::mesh::edge2cell_adjacent::create(outside_face));
  if(!ea.get()){
    cerr << "# [error] can not build edge2cell_adjacent." << endl;
    return __LINE__;
  }

  vector<pair<size_t,size_t> > polycube_edges;
  for(size_t ei = 0; ei < ea->edges_.size(); ++ei){
    const pair<size_t,size_t> & one_edge = ea->edges_[ei];
    const pair<size_t,size_t> & tri_pair = ea->edge2cell_[ei];
    if(ea->is_boundary_edge(tri_pair)) {
      polycube_edges.push_back(one_edge);
      continue;
    }
    boost::unordered_map<size_t,size_t>::const_iterator left =
        surface_type.find(outside_face_idx[tri_pair.first]);
    boost::unordered_map<size_t,size_t>::const_iterator right =
        surface_type.find(outside_face_idx[tri_pair.second]);
    if(left == surface_type.end() || right == surface_type.end()){
      cerr << "# [error] can not find surface_type." << endl;
      return __LINE__;
    }
    if(left->second != right->second)
      polycube_edges.push_back(one_edge);
  }

  {
    ofstream ofs("line2vtk.vtk") ;
    vector<size_t> lines(polycube_edges.size() * 2);
    for(size_t ei = 0; ei < polycube_edges.size(); ++ei){
      lines.push_back(polycube_edges[ei].first);
      lines.push_back(polycube_edges[ei].second);
    }
    line2vtk(ofs, &tm.node_[0], tm.node_.size(2), &lines[0], lines.size()/2);
  }

  vector<deque<pair<size_t,size_t> > > chains;
  jtf::util::extract_chain_from_edges(polycube_edges, chains);

  dump_singularity_to_vtk("chains.vtk", tm.node_, chains);

  boost::unordered_map<size_t, set<size_t> > corner_to_chain_idx;
  for(size_t ci = 0; ci < chains.size(); ++ci){
    const deque<pair<size_t,size_t> > & one_chain = chains[ci];
    corner_to_chain_idx[one_chain.front().first].insert(ci);
    corner_to_chain_idx[one_chain.back().second].insert(ci);
  }


  double radius = 0.002;
  if(argc == 5)
    radius = atof(argv[4]);

  vector<vector<size_t> > singularity_points(2);
  for(boost::unordered_map<size_t,set<size_t> >::const_iterator cit =
      corner_to_chain_idx.begin(); cit != corner_to_chain_idx.end(); ++cit){
    if(cit->second.size() < 4)
      singularity_points[0].push_back(cit->first);
    if(cit->second.size() > 4)
      singularity_points[1].push_back(cit->first);
  }

  cerr << "# [info] positive singularity " << singularity_points[0].size() << endl;
  cerr << "# [info] negative singularity " << singularity_points[1].size() << endl;

//  vector<deque<pair<size_t,size_t> > > temp_chains ;//= chains[78];
//  temp_chains.push_back(chains[78]);
  const string tet_edge_str = "shar_edge.obj";
  dump_singularity_to_cylinder(
        tet_edge_str.c_str(), tm.node_, chains,radius);
  const string singularity_obj = "surface_singularity_point";
  dump_surface_singularity_point(singularity_points, tm.node_,
                                 radius, argv[3], singularity_obj );

  return 0;
}

