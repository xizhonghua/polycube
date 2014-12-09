#include <iostream>
#include <numeric>
#include <stack>
#include <jtflib/mesh/io.h>
#include <jtflib/mesh/mesh.h>
#include <jtflib/util/util.h>

#include "../common/vtk.h"
#include "../common/util.h"
#include "../tetmesh/tetmesh.h"
#include "../tetmesh/hex_io.h"
#include "../tetmesh/util.h"
#include "../common/IO.h"

using namespace std;
using namespace zjucad::matrix;

int load_surface_restricted_type_static(
    const char * filename,
    boost::unordered_map<size_t,size_t> & result_surface_type);

int convert_surface_type_to_surface_patches(
    const zjucad::matrix::matrix<size_t> &tet,
    const boost::unordered_map<size_t,size_t> &surface_type,
    std::vector<zjucad::matrix::matrix<size_t> > &surface_patches,
    boost::unordered_set<std::pair<size_t,size_t> > &patch_boundary);

class polygon_patch
{
public:
  polygon_patch(const vector<matrix<size_t> >& patches,
                const vector<size_t> & patch_type,
                const matrix<double> & orig_node,
                const matrix<double> & polycube_node,
                const vector<deque<pair<size_t,size_t> > > & chains)
    : patches_(patches), patch_type_(patch_type),orig_node_(orig_node),
      polycube_node_(polycube_node), chains_(chains){}
  void analysis();
  void save_out(const char * file);
private:
  int unit_one_order_for_all_patches();

  int reverse_patch_corner_chain(const size_t & patch_idx);

  int check_patch_corner_edge_order(
      const size_t &other_patch_idx,
      const pair<size_t,size_t> & one_edge);

  int get_patch_boundary_chains(
      const matrix<double> & node,
      const matrix<size_t> &one_patch,
      const boost::unordered_set<size_t> &corners,
      vector<vector<pair<size_t,int> > > & chain_idx_rand_order,
      vector<deque<pair<size_t,size_t> > > & corner_chain,
      const size_t patch_idx);

  int reorder_chains_in_one_patch(
      vector<deque<pair<size_t,size_t> > > &corner_chain_of_patch,
      const matrix<size_t> &one_patch,
      const matrix<double> &node,
      const jtf::mesh::edge2cell_adjacent &ea);

  int get_normal_based_on_edge(
      const pair<size_t,size_t> &first_edge,
      const jtf::mesh::edge2cell_adjacent &ea,
      const matrix<size_t> &one_patch,
      const matrix<double> &node,
      matrix<double> &normal);

public:
  const matrix<double> &orig_node_;
  const matrix<double> &polycube_node_;
  vector<vector<deque<pair<size_t,size_t> > > > corner_chains_; // for each patch
  vector<size_t> corner_;
  vector<vector<vector<pair<size_t,int> > > > chain_idx_; // for each patch

  boost::unordered_map<pair<size_t,size_t>,vector<size_t> > corner_edge_to_patch_idx_;
  // pair0: chain_ends; pair1: chain_idx and order, 1: original, -1: reverse
  boost::unordered_map<pair<size_t,size_t>,pair<size_t,int> > chain_end_to_idx_;
  const vector<matrix<size_t> >& patches_;
  const vector<size_t> & patch_type_;
  const vector<deque<pair<size_t,size_t> > > & chains_;
};

void sort_edge(pair<size_t,size_t> & edge){
  if(edge.first > edge.second)
    swap(edge.first, edge.second);
}

int reverse_chain(deque<pair<size_t,size_t> > & chain)
{
  reverse(chain.begin(),chain.end());
  for(size_t t = 0; t < chain.size(); ++t){
    pair<size_t,size_t> & edge = chain[t];
    swap(edge.first,edge.second);
  }
  return 0;
}

int polygon_patch::get_normal_based_on_edge(
    const pair<size_t,size_t> &first_edge,
    const jtf::mesh::edge2cell_adjacent &ea,
    const matrix<size_t> &one_patch,
    const matrix<double> &node,
    matrix<double> &normal)
{
  const size_t edge_idx = ea.get_edge_idx(first_edge.first, first_edge.second);
  if(edge_idx == -1){
    cerr << "# [error] can not find edge of " << first_edge.first << " "
         << first_edge.second << endl;
    return __LINE__;
  }

  const pair<size_t,size_t> & tri_pair = ea.edge2cell_[edge_idx];
  if(!ea.is_boundary_edge(tri_pair)){
    cerr << "# [error] strange, this edge is not boundary edge "
         << first_edge.first << "," << first_edge.second << endl;
    return __LINE__;
  }

  const size_t adj_face_idx = (tri_pair.first == -1?tri_pair.second:tri_pair.first);
  const size_t other_point_idx =
      std::accumulate(one_patch(colon(), adj_face_idx).begin(),
                      one_patch(colon(), adj_face_idx).end(),0) -
      first_edge.first - first_edge.second;

  matrix<double> e0,e1;
  e0 = node(colon(),first_edge.second) - node(colon(), first_edge.first);
  e1 = node(colon(),other_point_idx) - node(colon(), first_edge.second);

  normal = cross(e0,e1);

  return 0;
}

int polygon_patch::reorder_chains_in_one_patch(
    vector<deque<pair<size_t,size_t> > > &corner_chain_of_patch,
    const matrix<size_t> &one_patch,
    const matrix<double> &node,
    const jtf::mesh::edge2cell_adjacent &ea)
{
  vector<size_t> chain_idx;
  pair<size_t,size_t> first_edge(-1,-1);
  matrix<double> normal_of_first_chain = zeros<double>(3,1);
  int is_reversed = 1; // 1: not reversed, -1: rversed
  for(size_t ci = 0; ci < corner_chain_of_patch.size(); ++ci){
    is_reversed = 1;
    deque<pair<size_t,size_t> > & one_corner_chain = corner_chain_of_patch[ci];
    pair<size_t,size_t> first_corner_edge = one_corner_chain.front();
    sort_edge(first_corner_edge);
    boost::unordered_map<pair<size_t,size_t>,pair<size_t,int> >::const_iterator cit =
        chain_end_to_idx_.find(first_corner_edge);
    if(cit == chain_end_to_idx_.end()){
      cerr << "# [error] can not find chain_idx of corner edge " << cit->first.first
           << " " << cit->first.second << endl;
      return __LINE__;
    }

    if(cit->second.second == 1) // in regular order
      first_edge = chains_[cit->second.first].front();
    else{
      // in reverse order
      first_edge = chains_[cit->second.first].back();
      swap(first_edge.first, first_edge.second);
      is_reversed = -1;
    }

    matrix<double> normal = zeros<double>(3,1);
    get_normal_based_on_edge(first_edge, ea, one_patch, node, normal);
    if(ci == 0){
      normal_of_first_chain = normal;
    }else{
      if(dot(normal_of_first_chain, normal) < 0){ // this definition satisify the boundary order in geneous != 0
        is_reversed *= -1;
      }
      if(is_reversed == -1){
        reverse_chain(one_corner_chain);
      }
    }
  }

  return 0;
}


int polygon_patch::get_patch_boundary_chains(
    const matrix<double> & node,
    const matrix<size_t> &one_patch,
    const boost::unordered_set<size_t> &corners,
    vector<vector<pair<size_t,int> > > & chain_idx_rand_order,
    vector<deque<pair<size_t,size_t> > > & corner_chain,
    const size_t patch_idx) // chain_idx, order
{
  chain_idx_rand_order.clear();

  unique_ptr<jtf::mesh::edge2cell_adjacent> ea(
        jtf::mesh::edge2cell_adjacent::create(one_patch));
  if(!ea.get()){
    cerr << "# [error] can not build edge2cell_adjacent." << endl;
    return __LINE__;
  }

  boost::unordered_set<size_t> patch_boundary_corner;
  for(size_t ei =0; ei < ea->edges_.size(); ++ei){
    const pair<size_t,size_t> & one_edge = ea->edges_[ei];
    const pair<size_t,size_t> & tri_pair = ea->edge2cell_[ei];
    if(ea->is_boundary_edge(tri_pair)) {
      if(corners.find(one_edge.first) != corners.end()){
        patch_boundary_corner.insert(one_edge.first);
      }
      if(corners.find(one_edge.second) != corners.end()){
        patch_boundary_corner.insert(one_edge.second);
      }
    }
  }

  vector<pair<size_t,size_t> > possible_edges;
  for(boost::unordered_set<size_t>::const_iterator cit =
      patch_boundary_corner.begin(); cit != patch_boundary_corner.end(); ++cit){
    boost::unordered_set<size_t>::const_iterator next = cit;
    ++next;
    for(; next != patch_boundary_corner.end(); ++next){
      possible_edges.push_back(make_pair(*cit, *next));
    }
  }

  vector<pair<size_t,size_t> > real_edges;
  for(size_t ei = 0; ei < possible_edges.size(); ++ei){
    pair<size_t,size_t> & edge = possible_edges[ei];
    sort_edge(edge);
    if(chain_end_to_idx_.find(edge) != chain_end_to_idx_.end())
      real_edges.push_back(edge);
  }

  //vector<deque<pair<size_t,size_t> > > chains;
  jtf::util::extract_chain_from_edges(real_edges, corner_chain);

  if(corner_chain.size() > 1)
    reorder_chains_in_one_patch(corner_chain, one_patch, node, *ea);

  for(size_t ci = 0; ci < corner_chain.size(); ++ci){
    const deque<pair<size_t,size_t> > & one_chain = corner_chain[ci];
    for(size_t ei = 0; ei < one_chain.size(); ++ei){
      pair<size_t,size_t> one_edge = one_chain[ei];
      sort_edge(one_edge);
      corner_edge_to_patch_idx_[one_edge].push_back(patch_idx);
    }
  }
  return 0;
}

void polygon_patch::save_out(const char * file)
{
  matrix<size_t> v2s = zeros<size_t>(orig_node_.size(2),1) * -1;

  {
    set<size_t> used_node;
    for(size_t pi = 0; pi < patches_.size(); ++pi){
      used_node.insert(patches_[pi].begin(), patches_[pi].end());
    }
    matrix<int> s2v(used_node.size(),1);
    copy(used_node.begin(), used_node.end(), s2v.begin());

    const string s2v_str = "s2v.mat";
    jtf::mesh::write_matrix(s2v_str.c_str(), s2v);

    matrix<double> orig_compress_node = zeros<double>(3, used_node.size());
    matrix<double> polycube_compress_node = zeros<double>(3, used_node.size());
    size_t p = 0;
    for(set<size_t>::const_iterator cit = used_node.begin();
        cit != used_node.end(); ++cit, ++p){
      orig_compress_node(colon(), p) = orig_node_(colon(),*cit);
      polycube_compress_node(colon(), p) = polycube_node_(colon(),*cit);
      v2s[*cit] = p;
    }
    {
      // dump to obj
      string obj_name = "new_tet_patch.obj";//file;
      //obj_name += ".new.obj";
      ofstream ofs(obj_name.c_str());
      for(size_t v = 0; v < orig_compress_node.size(2); ++v){
        ofs << "v " ;
        for(size_t di = 0; di < orig_compress_node.size(1); ++di){
          ofs << orig_compress_node(di,v) << " ";
        }
        ofs << endl;
      }

      // dump to obj
      string poly_obj_name = "init_tet_patch.obj";//file;
      //poly_obj_name += ".polycube.obj";
      ofstream ofs_polycube(poly_obj_name.c_str());
      for(size_t v = 0; v < polycube_compress_node.size(2); ++v){
        ofs_polycube << "v " ;
        for(size_t di = 0; di < polycube_compress_node.size(1); ++di){
          ofs_polycube << polycube_compress_node(di,v) << " ";
        }
        ofs_polycube << endl;
      }

      for(size_t p = 0; p < patches_.size(); ++p){
        const matrix<size_t>  & one_patch = patches_[p];

        for(size_t fi = 0; fi < one_patch.size(2); ++fi){
          ofs << "f ";
          ofs_polycube << "f ";
          for(size_t di = 0; di < one_patch.size(1); ++di){
            ofs << v2s[one_patch(di, fi)] + 1 << " ";
            ofs_polycube << v2s[one_patch(di, fi)] + 1 << " ";
          }
          ofs << endl;
          ofs_polycube << endl;
        }
      }
    }
  }


  { // dump to quad
    string quad_name = "init_tet_patch.quad";// file;
    //quad_name += ".quad";
    ofstream ofs(quad_name.c_str());
    // phase 0; dump out chains
    ofs << chains_.size() << endl;
    for(size_t ci = 0; ci < chains_.size(); ++ci){
      const deque<pair<size_t,size_t> > & one_chain = chains_[ci];
      ofs << one_chain.size() + 1 << endl; // points
      for(size_t ei = 0; ei < one_chain.size(); ++ei){
        ofs << v2s[one_chain[ei].first] << " ";
      }
      ofs << v2s[one_chain.back().second] << endl;
    }

    size_t passed_face_num = 0;
    // phase 1; dump out patch
    ofs << patches_.size() << endl;
    for(size_t pi = 0; pi < patches_.size(); ++pi){
      size_t corner_point_num = 0;
      const vector<deque<pair<size_t,size_t> > > & corner_chain = corner_chains_[pi];
      for(size_t di = 0; di < corner_chain.size(); ++di){
        corner_point_num += corner_chain[di].size();
      }
      ofs << corner_point_num << endl;
      for(size_t di = 0; di < corner_chain.size(); ++di){
        for(size_t ei = 0; ei < corner_chain[di].size(); ++ei)
          ofs << v2s[corner_chain[di][ei].first] << " " << corner_chain[di][ei].first
              << " " << patch_type_[pi] << endl;
      }
      ofs << corner_point_num << endl; // chain number must equal the corner number
      for(size_t di = 0; di < corner_chain.size(); ++di){
        for(size_t ei = 0; ei < corner_chain[di].size(); ++ei){
          pair<size_t,size_t> one_edge = corner_chain[di][ei];
          sort_edge(one_edge);
          boost::unordered_map<pair<size_t,size_t>, pair<size_t,int> >::const_iterator
              ppcit = chain_end_to_idx_.find(one_edge);
          if(ppcit == chain_end_to_idx_.end()){
            cerr << "# [error] strange can not find chain_end2idx "
                 << one_edge.first << " " << one_edge.second << endl;
            return ;
          }
          ofs << ppcit->second.first << " ";
        }
      }
      ofs << endl;

      // dump out all faces which belong to this patch
      ofs << patches_[pi].size(2) << endl;
      for(size_t fi = 0; fi < patches_[pi].size(2); ++fi){
        ofs << passed_face_num + fi << " ";
      }
      ofs << endl;
      passed_face_num += patches_[pi].size(2);
    }
  }

}

void polygon_patch::analysis()
{
  boost::unordered_set<size_t> corners;
  for(size_t ci = 0; ci < chains_.size(); ++ci){
    const deque<pair<size_t,size_t> > & one_chain = chains_[ci] ;
    pair<size_t,size_t> chain_ends(one_chain.front().first, one_chain.back().second);

    corners.insert(chain_ends.first);
    corners.insert(chain_ends.second);

    sort_edge(chain_ends);
    if(chain_ends.first == one_chain.front().first)
      chain_end_to_idx_[chain_ends] = make_pair(ci, 1); // in regular order
    else
      chain_end_to_idx_[chain_ends] = make_pair(ci, -1); // in reverse order
  }

  chain_idx_.resize(patches_.size());
  corner_chains_.resize(patches_.size());
  for(size_t pi = 0; pi < patches_.size(); ++pi){
    if(pi == 6)
      cerr << endl;
    get_patch_boundary_chains(polycube_node_, patches_[pi], corners, chain_idx_[pi],
                              corner_chains_[pi], pi);
  }

  unit_one_order_for_all_patches();

  // check whether the corner chain is correct for each patch
#if 1
  {// debug
    vector<size_t> lines;
    for(size_t pi = 0; pi < patches_.size(); ++pi){
      cerr << "# patch " << pi << endl;
      const vector<deque<pair<size_t,size_t> > > & corner_chains = corner_chains_[pi];
      for(size_t ci = 0; ci < corner_chains.size(); ++ci){
        cerr << "# -- chain " << ci << endl;
        const deque<pair<size_t,size_t> > & one_chain = corner_chains[ci];
        cerr << "# ---- edges ";
        for(size_t ei = 0; ei < one_chain.size(); ++ei){
          lines.push_back(one_chain[ei].first);
          lines.push_back(one_chain[ei].second);
          cerr << one_chain[ei].first << " " << one_chain[ei].second << " ";
        }
        cerr << endl;
      }
    }

    ofstream ofs("corner_edges.vtk");
    line2vtk(ofs, &polycube_node_[0], polycube_node_.size(2), &lines[0], lines.size()/2);
  }
#endif

}

int polygon_patch::reverse_patch_corner_chain(const size_t & patch_idx)
{
  vector<deque<pair<size_t,size_t> > > &patch_corners = corner_chains_[patch_idx];
  for(size_t ci = 0; ci < patch_corners.size(); ++ci){
    reverse_chain(patch_corners[ci]);
  }

  return 0;
}

// to check whether the corner edges of other_patch is the same as one_edge
int polygon_patch::check_patch_corner_edge_order(
    const size_t &other_patch_idx,
    const pair<size_t,size_t> & one_edge)
{
  for(size_t ci = 0; ci < corner_chains_[other_patch_idx].size(); ++ci){
    const deque<pair<size_t,size_t> > & one_chain = corner_chains_[other_patch_idx][ci];
    deque<pair<size_t,size_t> >::const_iterator eit =
        find(one_chain.begin(), one_chain.end(), one_edge);
    deque<pair<size_t,size_t> >::const_iterator eit_rev =
        find(one_chain.begin(), one_chain.end(), make_pair(one_edge.second, one_edge.first));

    if(eit == one_chain.end() && eit_rev == one_chain.end()) continue;
    if(eit != one_chain.end()) return 0; // the same order as one_edge
    if(eit_rev != one_chain.end()) return 1; // reverse order of one_edge
  }
  return -1; // can not find the edge
}

int polygon_patch::unit_one_order_for_all_patches()
{
  vector<bool> is_patch_visited(patches_.size(), false);
  stack<size_t> ps;
  ps.push(0);

  while(!ps.empty()){
    const size_t patch_idx = ps.top();
    ps.pop();
    is_patch_visited[patch_idx] = true;
    for(size_t ci = 0; ci < corner_chains_[patch_idx].size(); ++ci){
      const deque<pair<size_t,size_t> > &one_chain = corner_chains_[patch_idx][ci];
      for(size_t ei = 0; ei < one_chain.size(); ++ei){
        pair<size_t,size_t> one_edge = one_chain[ei];
        sort_edge(one_edge);
        boost::unordered_map<pair<size_t,size_t>, vector<size_t> >::const_iterator
            cit = corner_edge_to_patch_idx_.find(one_edge);
        if(cit == corner_edge_to_patch_idx_.end() ||
           cit->second.size() != 2){
          cerr << "# [error] strange corner edge must links two patches " << endl;
          return __LINE__;
        }

        size_t other_patch_idx = -1;
        for(size_t pi = 0; pi < cit->second.size(); ++pi)
          if(!is_patch_visited[cit->second[pi]]){
            other_patch_idx = cit->second[pi];
            break;
          }
        if(other_patch_idx == -1) continue;

        int rtn = check_patch_corner_edge_order(other_patch_idx, one_chain[ei]);
        if(rtn == 0) // the same as one_chain, need to reverse
          reverse_patch_corner_chain(other_patch_idx);
        if(rtn == -1){
          cerr << "# [error] strange can not find corner edge "
               << one_chain[ei].first << " " << one_chain[ei].second
               << " in patch " << other_patch_idx << endl;
          return __LINE__;
        }
        ps.push(other_patch_idx);
      }
    }
  }

  return 0;
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

int dump_surface_patch_for_szy(int argc, char * argv[])
{
  if(argc != 4){
    cerr << "# [usage] dump_surface_patch_for_szy tet  polycube_tet surface_type" << endl;
    return __LINE__;
  }

  jtf::mesh::meshes tm_polycube, tm;
  if(jtf::mesh::tet_mesh_read_from_zjumat(argv[1], &tm.node_, &tm.mesh_))
    return __LINE__;

  if(jtf::mesh::tet_mesh_read_from_zjumat(argv[2], &tm_polycube.node_, &tm_polycube.mesh_))
    return __LINE__;

  boost::unordered_map<size_t,size_t> surface_type;
  if(load_surface_restricted_type_static(argv[3], surface_type))
    return __LINE__;

  vector<matrix<size_t> > surface_patches;
  boost::unordered_set<pair<size_t,size_t>  > patch_boundary;
  convert_surface_type_to_surface_patches(tm.mesh_, surface_type,
                                          surface_patches, patch_boundary);

  {
    ofstream ofs("patch_from_surface_type.vtk");
    vector<size_t> surface_patch_vec;
    vector<size_t> surface_patch_type_vec;
    for(size_t pi = 0; pi < surface_patches.size(); ++pi){
      const matrix<size_t> & one_patch = surface_patches[pi];
      surface_patch_vec.insert(surface_patch_vec.end(),
                               one_patch.begin(), one_patch.end());
      for(size_t fi = 0; fi < one_patch.size(2); ++fi)
        surface_patch_type_vec.push_back(pi);
    }
    tri2vtk(ofs, &tm_polycube.node_[0], tm_polycube.node_.size(2), &surface_patch_vec[0], surface_patch_vec.size()/3);
    cell_data(ofs, &surface_patch_type_vec[0], surface_patch_type_vec.size(), "patch");
  }

  vector<deque<pair<size_t,size_t> > > chains;
  {
    vector<pair<size_t,size_t> > edges(patch_boundary.begin(),
                                       patch_boundary.end());
    jtf::util::extract_chain_from_edges(edges,  chains);
  }

  unique_ptr<jtf::mesh::face2tet_adjacent> fa(jtf::mesh::face2tet_adjacent::create(tm.mesh_));
  if(!fa.get()){
    cerr << "# [error] can not buildjtf::mesh::face2tet_adjacent." << endl;
    return __LINE__;
  }


  vector<size_t> patch_type(surface_patches.size());
  for(size_t pi = 0; pi < surface_patches.size(); ++pi){
    const size_t face_idx = fa->get_face_idx(&surface_patches[pi](0,0));
    if(face_idx == -1){
      cerr << "# [error] can not find face " << surface_patches[pi](colon(),0) << endl;
      return __LINE__;
    }
    boost::unordered_map<size_t,size_t>::const_iterator cit = surface_type.find(face_idx);
    if(cit == surface_type.end()){
      cerr << "# [error] can not find surface_type of " << face_idx << endl;
      return __LINE__;
    }
    patch_type[pi] = cit->second;
  }

  polygon_patch pp(surface_patches, patch_type, tm.node_, tm_polycube.node_, chains);
  pp.analysis();
  string name = argv[1];
  name += "_patch";
  pp.save_out(name.c_str());

  return 0;
}
