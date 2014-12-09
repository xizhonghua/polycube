#include <stack>
#include <fstream>
#include <sstream>
#include <numeric>
#include "../mesh_func/tri-normal.h"
#include "util.h"
#include "../common/vtk.h"
#include "../tetmesh/util.h"
#include "../common/util.h"
#include <jtflib/mesh/mesh.h>
#include <jtflib/mesh/util.h>
#include <jtflib/util/util.h>
#include "../tetmesh/tetmesh.h"

#include <hjlib/math/polar.h>
#include <hjlib/math/blas_lapack.h>

#include <verdict/verdict.h>


using namespace std;
using namespace zjucad::matrix;

bool is_degenerate(const zjucad::matrix::matrix<double> &tri)
{
  matrix<double> e1 = tri(colon(), 0)-tri(colon(), 2),
      e2 = tri(colon(), 1)-tri(colon(), 2);
  if(norm(cross(e1, e2)) < 1e-8) {
    return true;
  }
  return false;
}

int extra_surface_patch_according_to_type(
    const boost::unordered_map<size_t,size_t> &surface_type,
    const matrix<size_t> & outside_face_idx,
    const jtf::mesh::edge2cell_adjacent & ea,
    const jtf::mesh::face2tet_adjacent & fa,
    vector<boost::unordered_set<size_t> >  &surface_patches,
    boost::unordered_set<pair<size_t,size_t> > &patch_boundary)
{
  vector<bool> is_visited_face(outside_face_idx.size(), false);

  boost::unordered_map<size_t,size_t> face_idx2vec_idx;
  for(size_t fi = 0; fi < outside_face_idx.size(); ++fi){
    face_idx2vec_idx[outside_face_idx[fi]] = fi;
  }

  stack<size_t> face_stack;
  surface_patches.clear();
  vector<bool>::const_iterator it =
      find(is_visited_face.begin(), is_visited_face.end(), false);

  while(it != is_visited_face.end()){
    if(face_stack.empty()){
      face_stack.push(outside_face_idx[it-is_visited_face.begin()]);
      is_visited_face.at(it-is_visited_face.begin()) = true;
    }

    boost::unordered_map<size_t,size_t>::const_iterator bumcit =
        surface_type.find(face_stack.top());
    if(bumcit == surface_type.end()){
      cerr << "# [error] can not find surface in surface type." << endl;
      return __LINE__;
    }

    const size_t group_type = bumcit->second;

    boost::unordered_set<size_t> one_group;
    one_group.insert(face_stack.top());
    while(!face_stack.empty()){
      const size_t current_face = face_stack.top();
      face_stack.pop();
      const vector<size_t> & face_vec = fa.faces_[current_face];
      for(size_t i = 0; i < face_vec.size(); ++i){
        pair<size_t,size_t> edge(face_vec[i], face_vec[(i+1)%face_vec.size()]);
        if(edge.first > edge.second)
          swap(edge.first, edge.second);

        const pair<size_t,size_t> face_pair = ea.query(edge.first,edge.second);
        assert(outside_face_idx[face_pair.first] == current_face ||
               outside_face_idx[face_pair.second] == current_face);

        const size_t other_face_idx =
            (outside_face_idx[face_pair.first] + outside_face_idx[face_pair.second])
            - current_face;

        boost::unordered_map<size_t,size_t>::const_iterator cit =
            face_idx2vec_idx.find(other_face_idx);
        if(cit == face_idx2vec_idx.end()){
          cerr << "# [error] strange can not find face_idx " << other_face_idx
               << " to idx_vec." << endl;
          return __LINE__;
        }
        if(is_visited_face.at(cit->second) == true) {
          continue; // boundary edge
        }

        boost::unordered_map<size_t,size_t>::const_iterator type_cit =
            surface_type.find(other_face_idx);
        if(type_cit == surface_type.end()){
          cerr << "# [error] strange can not find surface type of "
               << other_face_idx << endl;
          return __LINE__;
        }
        if(type_cit->second != group_type){
          patch_boundary.insert(edge);
          continue; // boundary edge
        }else{
          is_visited_face.at(cit->second) = true;
          face_stack.push(other_face_idx);
          one_group.insert(other_face_idx);
        }
      }
    }

    surface_patches.push_back(one_group);
    it = find(is_visited_face.begin(), is_visited_face.end(), false);
  }
  return 0;
}

int sort_edge(pair<size_t,size_t> & edge)
{
  if(edge.first > edge.second)
    swap(edge.first, edge.second);
  return 0;
}



// return dot(p1-p0,p2-p0)/(norm(p1-p0, p2-p0))
double get_cos_angle(const matrix<double> & node,
                     const size_t p0, const size_t p1, const size_t p2)
{
  matrix<double> e1 = node(colon(), p1) - node(colon(),p0);
  matrix<double> e2 = node(colon(), p2) - node(colon(),p0);

  double len1 = norm(e1);
  if(len1 < 1e-6) len1 = 1.0;

  double len2 = norm(e2);
  if(len2 < 1e-6) len2 = 1.0;

  return dot(e1,e2)/(len1 * len2);
}

int convert_surface_type_to_surface_patches(
    const matrix<size_t> &tet,
    const boost::unordered_map<size_t,size_t> &surface_type,
    vector<matrix<size_t> > &surface_patches,
    boost::unordered_set<pair<size_t,size_t> > &patch_boundary)
{
  unique_ptr<jtf::mesh::face2tet_adjacent> fa(jtf::mesh::face2tet_adjacent::create(tet));
  if(!fa.get()){
    cerr << "# [error] can not buildjtf::mesh::face2tet_adjacent." << endl;
    return __LINE__;
  }

  matrix<size_t> outside_face, outside_face_idx;
  get_outside_face(*fa, outside_face,true);
  get_outside_face_idx(*fa, outside_face_idx);

  unique_ptr<jtf::mesh::edge2cell_adjacent> ea(
        jtf::mesh::edge2cell_adjacent::create(outside_face));
  if(!ea.get()){
    cerr << "# [error] can not build edge2cell_adjacent." << endl;
    return __LINE__;
  }

  vector<boost::unordered_set<size_t> > surface_patches_set;

  extra_surface_patch_according_to_type(surface_type, outside_face_idx, *ea, *fa,
                                        surface_patches_set, patch_boundary);

  surface_patches.resize(surface_patches_set.size());
  for(size_t pi = 0; pi < surface_patches.size(); ++pi){
    surface_patches[pi].resize(3,surface_patches_set[pi].size());
    const boost::unordered_set<size_t> & one_patch = surface_patches_set[pi];
    size_t fi = 0;
    for(boost::unordered_set<size_t>::const_iterator cit = one_patch.begin();
        cit != one_patch.end(); ++cit, ++fi){
      const vector<size_t> & one_face = fa->faces_[*cit];
      copy(one_face.begin(), one_face.end(), surface_patches[pi](colon(),fi).begin());
    }
  }

  return 0;
}

int smooth_boundary(
    const boost::unordered_set<pair<size_t,size_t> > & boundary,
    const matrix<double> & orig_node,
    matrix<double> & polycube_node,
    const size_t iter)
{
  if(iter == 0) return 0;
  vector<deque<pair<size_t,size_t> > > chains;
  vector<pair<size_t,size_t> > boundary_edges(boundary.begin(), boundary.end());
  jtf::util::extract_chain_from_edges(boundary_edges, chains);

  // smooth each boundary
  //for(size_t it = 0; it < iter; ++it){
  for(size_t li = 0; li < chains.size(); ++li){
    const deque<pair<size_t,size_t> > & one_chain = chains[li];
    vector<double> length_seg(one_chain.size()+1,0);
    for(size_t ei = 0; ei < one_chain.size(); ++ei) {
      length_seg[ei+1] = norm(orig_node(colon(), one_chain[ei].first) -
                              orig_node(colon(), one_chain[ei].second)) +
                         length_seg[ei];
    }

    matrix<double> dis = polycube_node(colon(), one_chain.back().second)
                         - polycube_node(colon(), one_chain.front().first);
    for(size_t ei = 0; ei < one_chain.size(); ++ei){
      const pair<size_t,size_t> & edge = one_chain[ei];
      //const pair<size_t,size_t> & edge_next = one_chain[ei+1];

      polycube_node(colon(), edge.second) =
          polycube_node(colon(), one_chain.front().first) +
          dis*length_seg[ei+1]/length_seg.back();
    }
  }
  //}
  return 0;
}

int smooth_patch(const vector<matrixst> & surface_patch,
                 const boost::unordered_set<pair<size_t,size_t> > &boundary,
                 matrixd & node,
                 const size_t iter)
{
  if(iter == 0) return 0;
  boost::unordered_set<size_t> fix_node;
  for(boost::unordered_set<pair<size_t,size_t> >::const_iterator cit =
      boundary.begin();  cit != boundary.end(); ++cit){
    fix_node.insert(cit->first);
    fix_node.insert(cit->second);
  }

  typedef boost::unordered_map<size_t, vector<size_t> > p2p_type;
  vector<p2p_type> p2p_group(surface_patch.size());

  for(size_t pi = 0; pi < surface_patch.size(); ++pi){
    unique_ptr<jtf::mesh::one_ring_point_at_point> orpap(
          jtf::mesh::one_ring_point_at_point::create(surface_patch[pi]));
    p2p_group[pi] = orpap->p2p_;
  }

  // remove fix_node
  for(size_t pi = 0; pi < p2p_group.size(); ++pi){
    p2p_type & p2p_ = p2p_group[pi];
    for(boost::unordered_set<size_t>::const_iterator cit = fix_node.begin();
        cit != fix_node.end(); ++cit){
      const size_t & node_idx = *cit;
      p2p_type::iterator pit = p2p_.find(node_idx);
      if(pit != p2p_.end()) p2p_.erase(pit);
    }
  }


  matrixd center_node = zeros<double>(3,1);
  for(size_t it = 0; it < iter; ++it){
    for(size_t pi = 0; pi < p2p_group.size(); ++pi){
      const p2p_type & p2p_ = p2p_group[pi];
      for(p2p_type::const_iterator pcit = p2p_.begin();
          pcit != p2p_.end(); ++pcit) {
        center_node = zeros<double>(3,1);
        const vector<size_t> & linked_node = pcit->second;
        for(size_t scidx = 0; scidx < linked_node.size(); ++scidx){
            center_node += node(colon(), linked_node[scidx]);
        }
        center_node /= linked_node.size();
        node(colon(), pcit->first) = center_node;
      }
    }
  }
  return 0;
}

int smooth_volume(
    const matrix<size_t> & tet,
    const vector<matrix<size_t> > &surface_patches,
    matrix<double> &node,
    const size_t iter_v )
{
  boost::unordered_set<size_t> fix_node;
  for(size_t pi = 0; pi < surface_patches.size(); ++pi){
    fix_node.insert(surface_patches[pi].begin(), surface_patches[pi].end());
  }

  jtf::tetmesh::one_ring_point_at_point orpap;
  orpap.build(tet);

  boost::unordered_map<size_t, boost::unordered_set<size_t> > p2p ;

  //p2p = orpap.p2p_;
  // only smooth flipped tet
  {
    boost::unordered_set<size_t> associated_nodes;
    double vol;
    for(size_t ti = 0; ti < tet.size(2); ++ti){
      vol = jtf::mesh::cal_tet_vol(node(colon(), tet(colon(),ti)));
      if(vol < 0){
        // flipped_tet.push_back(ti);
        associated_nodes.insert(tet(colon(),ti).begin(), tet(colon(),ti).end());
      }
    }
    cerr << "# [info] left " << associated_nodes.size() << " flipped nodes need"
         << " to be smoothed." << endl;
    for(boost::unordered_set<size_t>::const_iterator cit = associated_nodes.begin();
        cit != associated_nodes.end(); ++cit){
      boost::unordered_map<size_t, boost::unordered_set<size_t> >::const_iterator
          mscit = orpap.p2p_.find(*cit);
      if(mscit == orpap.p2p_.end()){
        cerr << "# [error] can not find p2p mapping." << endl;
        return __LINE__;
      }
      p2p[*cit] = mscit->second;
    }

  }


  for(boost::unordered_set<size_t>::const_iterator cit = fix_node.begin();
      cit != fix_node.end(); ++cit){
    boost::unordered_map<size_t, boost::unordered_set<size_t> >::iterator it
        = p2p.find(*cit);
    if(it  != p2p.end()) p2p.erase(it);
  }

  matrix<double> center_node = zeros<double>(3,1);
  for(size_t it = 0; it < iter_v; ++it){
    for(boost::unordered_map<size_t,boost::unordered_set<size_t> >::const_iterator
        cit = p2p.begin(); cit != p2p.end(); ++cit){
      center_node = zeros<double>(3,1);
      const boost::unordered_set<size_t> & linked_nodes = cit->second;
      for(boost::unordered_set<size_t>::const_iterator lit = linked_nodes.begin();
          lit != linked_nodes.end(); ++lit){
        center_node += node(colon(), *lit);
      }
      center_node /= linked_nodes.size();
      node(colon(), cit->first) = center_node;
    }
  }
  return 0;
}

int dump_out_normal_flipped_face(
    const zjucad::matrix::matrix<double> &node,
    const std::vector<zjucad::matrix::matrix<size_t> > &surface_patches,
    const std::vector<zjucad::matrix::matrix<double> > &patch_normal,
    const size_t i,
    const zjucad::matrix::matrix<double> &vtk_node)
{
  ostringstream os;
  os << i;
  string name = "flipped_face_iter_" + os.str() + ".vtk";

  vector<size_t> flipped_face;
  matrix<double> normal = zeros<double>(3,1), node_;
  for(size_t pi = 0; pi < surface_patches.size(); ++pi){
    const matrix<size_t> & one_patch = surface_patches[pi];
    for(size_t fi = 0; fi < one_patch.size(2); ++fi){
      node_ = node(colon(), one_patch(colon(),fi));
      calc_tri_normal_(&normal[0], &node_[0]);
      if(dot(normal, patch_normal[pi]) < 0.1)
        flipped_face.insert(flipped_face.end(),
                            one_patch(colon(),fi).begin(),
                            one_patch(colon(),fi).end());
    }
  }

  ofstream ofs(name.c_str());
  tri2vtk(ofs, &vtk_node[0], vtk_node.size(2), &flipped_face[0], flipped_face.size()/3);
  return 0;
}

int get_face_patches_according_to_boundary(
    const zjucad::matrix::matrix<size_t> &outside_face,
    const jtf::mesh::edge2cell_adjacent & ea,
    const boost::unordered_set<std::pair<size_t,size_t> > &boundary_edges,
    std::vector<std::vector<size_t> > &patches)
{
  // this function assume each edge in boundary_edges is increasing order
  vector<bool> is_face_visited(outside_face.size(2), false);

  stack<size_t> face_stack;

  patches.clear();

  vector<bool>::const_iterator cit = find(is_face_visited.begin(),
                                          is_face_visited.end(), false);
  while(cit != is_face_visited.end()){
    vector<size_t> one_patch;
    face_stack.push(cit-is_face_visited.begin());
    while(!face_stack.empty()){
      const size_t f_idx = face_stack.top();
      face_stack.pop();
      is_face_visited[f_idx] = true;
      one_patch.push_back(f_idx);

      for(size_t pi = 0; pi < outside_face.size(1); ++pi){
        pair<size_t,size_t> one_edge(
              outside_face(pi,f_idx),
              outside_face((pi+1)%outside_face.size(1),f_idx));
        if(one_edge.first > one_edge.second)
          swap(one_edge.first, one_edge.second);
        if(boundary_edges.find(one_edge) != boundary_edges.end()) continue;

        const size_t edge_idx = ea.get_edge_idx(one_edge.first, one_edge.second);
        if(edge_idx == -1){
          cerr << "# [error] strange can not find edge idx of "
               << one_edge.first << " " << one_edge.second << endl;
          return __LINE__;
        }
        const pair<size_t,size_t> & tri_pair = ea.edge2cell_[edge_idx];
        if(ea.is_boundary_edge(tri_pair)) continue;
        if(f_idx != tri_pair.first && f_idx != tri_pair.second){
          cerr << "# [error] strange face pair of edge "
               << one_edge.first << " " << one_edge.second << " is "
               << tri_pair.first << " " << tri_pair.second << " without "
               << f_idx << endl;
          return __LINE__;
        }
        const size_t other_face_idx = tri_pair.first + tri_pair.second - f_idx;
        if(is_face_visited[other_face_idx]) continue;
        face_stack.push(other_face_idx);
      }
    }
    patches.push_back(one_patch);
    cit = find(is_face_visited.begin(), is_face_visited.end(), false);
  }

  return 0;
}
