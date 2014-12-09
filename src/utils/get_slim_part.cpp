#include <jtflib/util/vertex_connection.h>
#include <jtflib/mesh/io.h>
#include <jtflib/mesh/util.h>

#include "../tetmesh/tetmesh.h"
#include "../tetmesh/hex_io.h"
#include "../tetmesh/util.h"
#include "../hex_param/remove_surface_wedge.h"
#include "../common/vtk.h"
#include "../common/util.h"

/////////////////////////////////////////////////////////////////////
/// WARNING !!!  This function takes two assumption:
///               1. input jtf::mesh::meshes is almost aligned with x,y,z axes
///               2. input jtf::mesh::meshes is not a cut-open mesh
////////////////////////////////////////////////////////////////////

using namespace std;
using namespace zjucad::matrix;

///////////////////  interaction  /////////////////////////////////
int get_slim_part(const matrix<double> & node, const matrix<size_t> & tet,
                  std::vector<vector<size_t> > & slim_faces);
///////////////////////////////////////////////////////////////////

int extract_polycube_surface_type(
    const matrixst & tet,
    const matrixd & node,
    const matrixst & outside_face,
    const matrixst & outside_face_idx,
    boost::unordered_map<size_t,size_t> & surface_type)
{
  matrixd face_normal_in_polycube;

  jtf::mesh::cal_face_normal(outside_face, node, face_normal_in_polycube);

  const matrixd axes = eye<double>(3);
  vector<pair<double,int> > axis_choice(6);

  matrixst outside_face_type(outside_face.size(2));

  for(size_t fi = 0; fi < outside_face.size(2); ++fi){
    for(size_t ai = 0; ai < 3; ++ai){
      axis_choice[2*ai] =
          make_pair( dot(face_normal_in_polycube(colon(),fi), axes(colon(),ai)),
                     2*ai);
      axis_choice[2*ai+1] =
          make_pair(-1*dot(face_normal_in_polycube(colon(),fi), axes(colon(),ai)),
                    2*ai+1);
    }
    sort(axis_choice.begin(), axis_choice.end());
    outside_face_type[fi] = axis_choice.back().second/2;
  }

  surface_type.clear();
  for(size_t fi = 0; fi < outside_face.size(2); ++fi){
    surface_type[outside_face_idx[fi]] = outside_face_type[fi];
  }

  return 0;
}

int get_patch_arounded_by_boundary_and_corner(
    const zjucad::matrix::matrix<double> &node,
    const jtf::mesh::edge2cell_adjacent & ea,
    const matrix<size_t> & patch_face,
    const matrix<size_t> & face_idx,
    const deque<pair<size_t,size_t> > &one_chain,
    const size_t p0, const size_t p1,
    const pair<size_t,size_t> & extreme_edge,
    vector<vector<size_t> > &slim_faces)
{
  map<pair<size_t,size_t>, double > edge_weight;
  for(size_t ei = 0; ei < ea.edges_.size(); ++ei){
    const pair<size_t,size_t> & one_edge = ea.edges_[ei];
    const pair<size_t,size_t> & tri_pair = ea.edge2cell_[ei];
    edge_weight[one_edge] = norm(node(colon(), one_edge.first)
                                 - node(colon(), one_edge.second));
  }

  unique_ptr<vertex_connection<UNDIRECT> > vc(
        vertex_connection<UNDIRECT>::create(edge_weight));
  if(!vc.get()){
    cerr << "# [error] can not build vertex_connection." << endl;
    return __LINE__;
  }

  vector<size_t> path;
  vc->get_shortest_path(p0,p1, path);
  if(path.empty()){
    cerr << "# [error] strange can not find path from " << p0 << " " << p1
         << endl;
    return __LINE__;
  }

  // check, if the shortes edges are on boundary, then this patch is discardes.
  size_t boundary_short_edge = 0;
  for(size_t p = 0; p < path.size()-1; ++p){
    const size_t edge_idx = ea.get_edge_idx(path[p], path[p+1]);
    if(edge_idx == -1){
      cerr << "# [error] strange, can not find edge " << path[p] << " "
           << path[p+1] << endl;
      return __LINE__;
    }
    if(ea.is_boundary_edge(ea.edge2cell_[edge_idx])) ++boundary_short_edge;
  }

  if(boundary_short_edge == path.size()-1) // this patch should be addressed by others
    return 0;

  // extract patch with given boundary
  boost::unordered_set<pair<size_t,size_t> > boundary_edges;
  for(size_t ei = 0; ei < one_chain.size(); ++ei){
    const pair<size_t,size_t> & one_edge = one_chain[ei];
    if(one_edge.first > one_edge.second)
      boundary_edges.insert(make_pair(one_edge.second, one_edge.first));
    else
      boundary_edges.insert(one_edge);
  }
  for(size_t p = 0; p < path.size()-1; ++p){
    if(path[p] > path[p+1])
      boundary_edges.insert(make_pair(path[p+1], path[p]));
    else
      boundary_edges.insert(make_pair(path[p], path[p+1]));
  }

  stack<size_t> face_stak;
  // assign the face seed
  const pair<size_t,size_t> tri_pair =
      ea.query(extreme_edge.first, extreme_edge.second);
  assert(ea.is_boundary_edge(tri_pair));
  const size_t seed_face_idx = (tri_pair.first == -1?tri_pair.second: tri_pair.first);

  face_stak.push(seed_face_idx);
  vector<size_t> is_face_visited(face_idx.size(),false);
  while(!face_stak.empty()){
    const size_t face_idx_ = face_stak.top();
    face_stak.pop();
    is_face_visited[face_idx_] = true;
    for(size_t p = 0; p < patch_face.size(1); ++p){
      pair<size_t,size_t> edge(
            patch_face(p, face_idx_),
            patch_face((p+1)%patch_face.size(1), face_idx_));
      if(edge.first > edge.second)
        swap(edge.first, edge.second);
      if(boundary_edges.find(edge) != boundary_edges.end()) continue;
      const pair<size_t,size_t> tri_pair = ea.query(edge.first, edge.second);
      assert(tri_pair.first == face_idx_ || tri_pair.second == face_idx_);
      const size_t other_face_idx = tri_pair.first + tri_pair.second - face_idx_;
      if(!is_face_visited[other_face_idx])
        face_stak.push(other_face_idx);
    }
  }

  vector<size_t> slim_face_patch;
  jtf::mesh::euler_number eun;
  for(size_t fi = 0; fi < face_idx.size(); ++fi){
    if(is_face_visited[fi]) {
      slim_face_patch.push_back(face_idx[fi]);
      eun.add_face(&patch_face(0,fi), 3);
    }
  }

  // if the slim face is a simple patch, and size is not large,
  // or this wedge should be detected in another patch
  if(eun() == 1 && slim_face_patch.size() * 2 < face_idx.size())
    slim_faces.push_back(slim_face_patch);

  return 0;
}

int get_slim_face_by_checking_patch_boundary(
    const boost::unordered_set<size_t> & one_patch,
    const matrix<double> &node,
    const jtf::mesh::face2tet_adjacent & fa,
    vector<vector<size_t> > &slim_faces)
{
  matrix<size_t> one_patch_face(3, one_patch.size());
  matrix<size_t> face_idx(one_patch.size(), 1);
  size_t i = 0;
  for(boost::unordered_set<size_t>::const_iterator sit = one_patch.begin();
      sit != one_patch.end(); ++sit, ++i){
    face_idx[i] = *sit;
    const vector<size_t> & one_face = fa.faces_[*sit];
    copy(one_face.begin(), one_face.end(), one_patch_face(colon(),i).begin());
  }

  vector<deque<pair<size_t,size_t> > > chains;
  unique_ptr<jtf::mesh::edge2cell_adjacent> ea(
        jtf::mesh::edge2cell_adjacent::create(one_patch_face));
  {
    if(!ea.get()){
      cerr << "# [error] can not build edge2cell_adjacent." << endl;
      return __LINE__;
    }

    vector<pair<size_t,size_t> > boundary_edges;
    for(size_t ei = 0; ei < ea->edges_.size(); ++ei){
      if(ea->is_boundary_edge(ea->edge2cell_[ei]))
        boundary_edges.push_back(ea->edges_[ei]);
    }
    jtf::util::extract_chain_from_edges(boundary_edges, chains);
  }

  cerr << "# [info] this patch contains " << chains.size() << " chains." << endl;
  vector<int> corner;
  for(size_t ci = 0; ci < chains.size(); ++ci){
    const deque<pair<size_t,size_t> > & one_chain = chains[ci];
    //assert(one_chain.front().first == one_chain.back().second);
    corner.resize(one_chain.size());

    for(size_t ai = 0; ai < one_chain.size(); ++ai){
      const pair<size_t,size_t> & current_edge = one_chain[ai];
      const pair<size_t,size_t> & next_edge = one_chain[(ai+1)%one_chain.size()];
      double dot_val = dot(node(colon(), current_edge.second)
                           - node(colon(), current_edge.first),
                           node(colon(), next_edge.second)
                           - node(colon(), next_edge.first));
      double len_0 = norm(node(colon(), current_edge.second)
                          - node(colon(), current_edge.first));
      if(len_0 < 1e-6) len_0 = 1.0;
      double len_1 = norm(node(colon(), next_edge.second)
                          - node(colon(), next_edge.first));

      if(len_1 < 1e-6) len_1 = 1.0;
      dot_val /= len_0 * len_1;

      if(fabs(dot_val) < sqrt(2.0)/2.0){
        corner[ai] = 1; // orthagnoal angle
      }else if(dot_val < 0){
        corner[ai] = -1; // extrema angle
      }else{
        corner[ai] = 0; // flat angle
      }
    }

    vector<int>::const_iterator extrema_cit = find(corner.begin(), corner.end(), -1);
    if(extrema_cit == corner.end()) continue;

    deque<pair<size_t,size_t> > wedge;
    while(extrema_cit != corner.end()){
      wedge.clear();
      const size_t p = extrema_cit - corner.begin();
      wedge.push_back(one_chain[p]);
      int next = -1, prev = -1;
      for(int i = p+1; i % corner.size() != p; ++i){
        wedge.push_back(one_chain[i%corner.size()]);
        if(corner[i] != 0){
          next = i;
          break;
        }
      }

      for(int i = p-1; i % corner.size() != p; --i){
        if(i < 0) i += corner.size();
        wedge.push_front(one_chain[i%corner.size()]);
        if(corner[i] !=0){
          prev = i;
          break;
        }
      }

      if(next == -1 || prev == -1){
        // can not find an orthognal corner for such extrema point, may be with noise.
        break;
      }

      wedge.pop_front();
      get_patch_arounded_by_boundary_and_corner(
            node, *ea, one_patch_face, face_idx, wedge, wedge.front().first,
            wedge.back().second, one_chain[p] ,slim_faces);

      extrema_cit = find(corner.begin() + p + 1, corner.end(), -1);
    }
  }

  return 0;
}

int get_slim_part(const matrix<double> & node, const matrix<size_t> & tet,
                  std::vector<vector<size_t> > & slim_faces)
{
  unique_ptr<jtf::mesh::face2tet_adjacent> fa(jtf::mesh::face2tet_adjacent::create(tet));
  if(!fa.get()){
    cerr << "# [error] can not build face2tet_adjacnet." << endl;
    return __LINE__;
  }

  matrixst outside_face, outside_face_idx;
  get_outside_face(*fa, outside_face);
  get_outside_face_idx(*fa, outside_face_idx);

  boost::unordered_map<size_t,size_t> surface_type;
  extract_polycube_surface_type(tet, node, outside_face, outside_face_idx, surface_type);

  unique_ptr<jtf::mesh::edge2cell_adjacent> ea(
        jtf::mesh::edge2cell_adjacent::create(outside_face));
  if(!ea.get()){
    cerr << "# [error] can not build edge2cell_adjacent." << endl;
    return __LINE__;
  }

  // to fit the function
  std::vector<std::pair<size_t,size_t> >  g_unknown_face_pair;
  matrixst cut_tet2tet(max(tet)+1);
  cut_tet2tet(tet) = tet(colon());


  std::vector<boost::unordered_set<size_t> > patches;
  boost::unordered_map<size_t,boost::unordered_set<size_t> > group_linking_info;
  extract_surface_patch_graph(
        tet, tet, node, cut_tet2tet, outside_face, outside_face_idx, *fa, *fa,
        *ea, g_unknown_face_pair, surface_type, patches, group_linking_info);


  // patch degree less than 4 is a slim part, and parts of the other patch may be
  // slim too. For each other patch, we travel along their boundary,
  // and check the edge direction, if dot(e1,e2) < 0, it's an extrema,
  // we will find the nearest two edge corner, and link them to make a slim patch
  slim_faces.clear();
  vector<size_t> one_patch;
  for(boost::unordered_map<size_t, boost::unordered_set<size_t> >::const_iterator
      cit = group_linking_info.begin(); cit != group_linking_info.end(); ++cit){
    one_patch.clear();
    if(cit->second.size() < 4){ // pure topology check
      one_patch.resize(patches[cit->first].size());
      copy(patches[cit->first].begin(), patches[cit->first].end(), one_patch.begin());
      slim_faces.push_back(one_patch);
    }else{// use geometry to check other patches
      get_slim_face_by_checking_patch_boundary(patches[cit->first], node, *fa, slim_faces);
    }
  }

  return 0;
}


/////////////////////  example  ///////////////////////////////////
int vis_slim_part(int argc, char * argv[])
{
  if(argc != 2){
    cerr << "# [usage] get_slim_part polycube_tet" << endl;
    return __LINE__;
  }

  jtf::mesh::meshes tm;
  if(jtf::mesh::tet_mesh_read_from_zjumat(argv[1], &tm.node_, &tm.mesh_))
    return __LINE__;

  vector<vector<size_t> > slim_faces;
  get_slim_part(tm.node_, tm.mesh_, slim_faces);

  {// visual
    ofstream ofs("slim_face.vtk");
    unique_ptr<jtf::mesh::face2tet_adjacent> fa(jtf::mesh::face2tet_adjacent::create(tm.mesh_));
    vector<size_t> faces;
    vector<size_t> face_patch_idx;
    for(size_t pi = 0; pi < slim_faces.size(); ++pi){
      const vector<size_t> & one_patch = slim_faces[pi];
      for(size_t fi = 0; fi < one_patch.size(); ++fi){
        const size_t & face_idx = one_patch[fi];
        const vector<size_t> & one_face = fa->faces_[face_idx];
        faces.insert(faces.end(), one_face.begin(), one_face.end());
        face_patch_idx.push_back(pi);
      }
    }

    tri2vtk(ofs, &tm.node_[0], tm.node_.size(2), &faces[0], faces.size()/3);
    cell_data(ofs, &face_patch_idx[0], face_patch_idx.size(), "patch_idx");
  }
  return 0;
}
