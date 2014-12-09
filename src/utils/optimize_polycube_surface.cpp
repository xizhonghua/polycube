#include <map>
#include <stack>
#include <fstream>

#include <jtflib/util/vertex_connection.h>
#include <jtflib/mesh/mesh.h>
#include <jtflib/mesh/io.h>
#include <jtflib/util/util.h>

#include "../tetmesh/hex_io.h"
#include "../tetmesh/tetmesh.h"
#include "../common/vtk.h"
#include "../hex_param/io.h"
using namespace std;
using namespace zjucad::matrix;

//static void load_surface_type(
//    const char * filename,
//    map<size_t,size_t> & surface_type)
//{
//  ifstream ifs(filename);
//  if(ifs.fail())
//    throw std::logic_error("# [error] can not open surface restricted type.");

//  size_t face_idx, face_type;
//  while(!ifs.eof()){
//      ifs >> face_idx >> face_type;
//      surface_type[face_idx] = face_type;
//    }
//  return;
//}

void build_face_patch(const matrix<size_t> & outside_face,
                      const matrix<size_t> &outside_face_type,
                      const jtf::mesh::edge2cell_adjacent &ea,
                      vector<set<size_t> > &surface_patches)
{
  using namespace jtf::mesh;

  vector<bool> is_visited_face(outside_face_type.size(), false);

  stack<size_t> face_stack;
  surface_patches.clear();
  vector<bool>::const_iterator it =
      find(is_visited_face.begin(), is_visited_face.end(), false);

  while(it != is_visited_face.end()){
      if(face_stack.empty()){
          face_stack.push(static_cast<size_t>(it-is_visited_face.begin()));
          is_visited_face.at(it-is_visited_face.begin()) = true;
        }

      const size_t group_type = outside_face_type[face_stack.top()];

      set<size_t> one_group;

      one_group.insert(face_stack.top());

      while(!face_stack.empty()){
          const size_t current_face_idx = face_stack.top();
          face_stack.pop();

          for(size_t i = 0; i < outside_face.size(1); ++i){
              pair<size_t,size_t> edge(outside_face(i, current_face_idx),
                                       outside_face((i+1)%outside_face.size(1), current_face_idx));

              if(edge.first > edge.second)
                swap(edge.first, edge.second);

              const pair<size_t,size_t> face_pair = ea.query(edge.first,edge.second);

              const size_t other_face_idx =
                  face_pair.first + face_pair.second -current_face_idx;

              if(outside_face_type[other_face_idx] != group_type
                 || is_visited_face.at(other_face_idx) == true)
                continue;

              is_visited_face.at(other_face_idx) = true;
              face_stack.push(other_face_idx);
              one_group.insert(other_face_idx);
            }
        }
      surface_patches.push_back(one_group);
      it = find(is_visited_face.begin(), is_visited_face.end(), false);
    }
}

int update_patch_boundary(const size_t patch_idx,
                          const jtf::mesh::edge2cell_adjacent &ea,
                          const set<size_t> &corners,
                          const matrix<double> & node,
                          const matrix<size_t> & outside_face,
                          vector<set<size_t> > &surface_patches)
{
  set<size_t> &current_patch = surface_patches[patch_idx];
  map<pair<size_t,size_t>, double> edge_len;
  set<pair<size_t,size_t> > boundary_edges;
  for(set<size_t>::const_iterator cit = current_patch.begin();
      cit != current_patch.end(); ++cit){
      for(size_t pi = 0; pi < outside_face.size(1); ++pi){
          pair<size_t,size_t> one_edge(outside_face(0,*cit),
                                       outside_face(1,*cit));
          if(one_edge.first > one_edge.second) swap(one_edge.first, one_edge.second);
          edge_len[one_edge] = norm(node(colon(), one_edge.first) -
                                    node(colon(), one_edge.second));
          const pair<size_t,size_t> face_pair = ea.query(one_edge.first,
                                                         one_edge.second);
          if(find(current_patch.begin(), current_patch.end(), face_pair.first)
             != current_patch.end() &&
             find(current_patch.begin(), current_patch.end(), face_pair.second)
                          != current_patch.end() )
          continue;
          boundary_edges.insert(one_edge);
        }
    }

  vector<pair<size_t,size_t> > boundary_edge_vec(boundary_edges.begin(), boundary_edges.end());
  vector<deque<pair<size_t,size_t> > > chains;
  jtf::util::extract_chain_from_edges(boundary_edge_vec, chains);

  std::shared_ptr<vertex_connection<UNDIRECT> > vc(vertex_connection<UNDIRECT>::create(edge_len));
  if(!vc.get()){
      cerr << "# [error] can not build vertex_connection." << endl;
      return __LINE__;
    }

  for(size_t ci = 0; ci < chains.size(); ++ci){
      deque<pair<size_t,size_t> > & one_chain = chains[ci];
      while(find(corners.begin(), corners.end(), one_chain.front().first)
            == corners.end()) {
          one_chain.push_back(one_chain.front());
          one_chain.pop_front();
        }
    }
  return 0;
}

int optimize_patch_boundary(
    const matrix<size_t> & outside_face,
    matrix<size_t> &outside_face_type,
    const vertex_connection<UNDIRECT> & vc,
    const jtf::mesh::edge2cell_adjacent & ea,
    vector<set<size_t> > &surface_patches)
{
  // this function fix all patch corner and strength path between them in each patch
  // step 0: gather all corners
  set<size_t> corners;
  {
    map<size_t, set<size_t> > p2types;
    for(size_t fi = 0; fi < outside_face.size(2); ++fi){
        for(size_t pi = 0; pi < outside_face.size(1); ++pi){
            p2types[outside_face(pi,fi)].insert(outside_face_type[fi]);
          }
      }
    for(map<size_t,set<size_t> >::const_iterator cit = p2types.begin();
        cit != p2types.end(); ++cit){
        if(cit->second.size() > 2)
          corners.insert(cit->first);
      }
  }

  // step 1: for each patch, extract its boundarty and calculate each neighbour corner on it
  for(size_t pi = 0; pi < surface_patches.size(); ++pi){
//      const set<size_t> & one_patch = surface_patches[pi];
//      matrix<size_t> patch_faces(outside_face.size(1), one_patch.size());
//      {
//        size_t i = 0;
//        for(set<size_t>::const_iterator scit = one_patch.begin();
//            scit != one_patch.end(); ++scit, ++i){
//            patch_faces(colon(),i) = outside_face(colon(), *scit);
//          }
//        unique_ptr<jtf::mesh::edge2cell_adjacent> ea(jtf::mesh::edge2cell_adjacent::create(patch_faces));
//        matrix<size_t> boundary_edges;
//        jtf::mesh::get_boundary_edge(*ea, boundary_edges);

//        // find all corner associated to this boundary
//        set<size_t> associated_corners;
//        for(size_t ppi = 0; ppi < boundary_edges.size(); ++ppi){
//            if(find(corners.begin(), corners.end(), boundary_edges[ppi])
//               != corners.end()){
//                associated_corners.insert(boundary_edges[ppi]);
//              }
//          }

//        vector<pair<size_t,size_t> > boundary_edges_pair;
//        for(size_t ei = 0; ei < boundary_edges.size(2); ++ei)
//          boundary_edges_pair.insert(make_pair(boundary_edges(0,ei), boundary_edges)(1,ei));

//        vector<deque<pair<size_t,size_t> > > chains;
//        jtf::util::extract_chain_from_edges(boundary_edges_pair, chains);

        //update_patch_boundary(pi, ea, corners, surface_patches);
      }

  return 0;
}

int optimize_polycube_surface(const jtf::mesh::meshes & tm,
                              boost::unordered_map<size_t,size_t> & surface_type)
{
  std::shared_ptr<jtf::mesh::face2tet_adjacent> fa(jtf::mesh::face2tet_adjacent::create(tm.mesh_));
  if(!fa.get()){
      cerr << "# [error] can not buildjtf::mesh::face2tet_adjacent." << endl;
      return __LINE__;
    }
  matrix<size_t> outside_face, outside_face_idx;
  get_outside_face(*fa, outside_face);
  get_outside_face_idx(*fa, outside_face_idx);

  matrix<size_t> outside_face_type(outside_face_idx.size(),1);
  for(size_t fi = 0; fi < outside_face_idx.size(); ++fi)
    outside_face_type[fi] = surface_type[outside_face_idx[fi]];

  std::shared_ptr<jtf::mesh::edge2cell_adjacent> ea(
        jtf::mesh::edge2cell_adjacent::create(outside_face));
  if(!ea.get()){
      cerr << "# [error] can not build edge2cell_adjacent." << endl;
      return __LINE__;
    }

  map<pair<size_t,size_t>, double> edge_weight;
  for(size_t ei = 0; ei < ea->edges_.size(); ++ei){
      const pair<size_t,size_t> & one_edge = ea->edges_[ei];
      const double len = norm(tm.node_(colon(),one_edge.first)
                              - tm.node_(colon(), one_edge.second));
      edge_weight[one_edge] = len;
    }
  std::shared_ptr<vertex_connection<UNDIRECT> > vc(
        vertex_connection<UNDIRECT>::create(edge_weight));
  if(!vc.get()){
      cerr << "# [error] can not build vertex connection." << endl;
      return __LINE__;
    }

  vector<set<size_t> > surface_patches;
  build_face_patch(outside_face, outside_face_type, *ea, surface_patches);

  //optimize_patch_boundary(outside_face_type, *ea, surface_patches);

  return 0;
}


int optimize_polycube_surface(int argc, char * argv[])
{
  cerr << "# [error] this function is not available." << endl;
  return __LINE__;

  if(argc != 4){
      cerr << "# [usage] optimize_polycube_surface orig_tet surface_type new_surface_type" << endl;
      return __LINE__;
    }

  jtf::mesh::meshes tm;
  if(jtf::mesh::tet_mesh_read_from_zjumat(argv[1], &tm.node_, & tm.mesh_))
    return __LINE__;

  boost::unordered_map<size_t,size_t> surface_type;
  load_surface_type(argv[2], surface_type);

  optimize_polycube_surface(tm, surface_type);
  return 0;
}
