#include "hex_param.h"
#include "topology_operation.h"
#include "common.h"
#include <stack>
#include <iostream>
#include <fstream>
#include <jtflib/mesh/mesh.h>
#include <jtflib/util/vertex_connection.h>

#include <jtflib/util/container_operation.h>
#include "../common/transition_type.h"
#include "../common/vtk.h"


using namespace std;
using namespace zjucad::matrix;
using namespace jtf::mesh;

bool is_edges_in_one_tet(const pair<size_t,size_t> &edge0,
                         const pair<size_t,size_t> &edge1,
                         const matrixst &new_tet)
{
  assert(edge0.second == edge1.first); // to ensure the edge is linked
  size_t vertex[3] = {edge0.first,edge0.second,edge1.second};
  for(size_t t = 0; t < new_tet.size(2); ++t){
    const matrixst &tet  = new_tet(colon(),t);
    if(find(tet.begin(),tet.end(),vertex[0]) != tet.end()
       && find(tet.begin(),tet.end(),vertex[1]) != tet.end()
       && find(tet.begin(),tet.end(),vertex[2]) != tet.end()){
      return true;
    }
  }
  return false;
}

int find_adjacent_tet(const pair<size_t,size_t> &edge,
                      const matrixst &new_tet,
                      vector<size_t> &adjacent_tet)
{
  for(size_t t = 0; t < new_tet.size(2); ++t){
    const matrixst &tet = new_tet(colon(),t);
    if(find(tet.begin(),tet.end(),edge.first) != tet.end()
       && find(tet.begin(),tet.end(),edge.second) != tet.end())
      adjacent_tet.push_back(t);
  }
  return 0;
}

int split_tet_and_copy_zyz_node(const pair<size_t,size_t> &edge,
                                const vector<size_t> &adjacent_tet,
                                const matrixst &tet,
                                const matrixd &node,
                                const matrixd &zyz,
                                matrixst &new_tet,
                                matrixd &new_node,
                                matrixd &new_zyz)
{
  // add new_vertex
  new_node.resize(3,node.size(2) + 1);
  new_node(colon(),colon(0,node.size(2) - 1)) = node;
  matrixd new_vertex = (node(colon(),edge.first) + node(colon(),edge.second)) / 2.0;
  new_node(colon(), node.size(2)) = new_vertex; // add a new vertex

  // split the tet mesh
  new_tet.resize(4,tet.size(2) + adjacent_tet.size());
  new_tet(colon(),colon(0,tet.size(2) - 1)) = tet;
  for(size_t t = 0; t < adjacent_tet.size(); ++t){
    new_tet(colon(),tet.size(2) + t) = tet(colon(),adjacent_tet[t]); // fill the new tet
  }
  const size_t new_vertex_idx = node.size(2);
  for(size_t t = 0; t < adjacent_tet.size(); ++t){
    //        matrixst &old_tet_  = new_tet(colon(),adjacent_tet[t]);
    //        matrixst &new_tet_  = new_tet(colon(),tet.size(2) + t);

    // replace edge.first in old tet with new_vertex_idx and replace edge.second in new_tet with new_vertex_idx
    for(size_t i = 0; i < 4; ++i){
      if(new_tet(i,adjacent_tet[t]) == edge.first) new_tet(i,adjacent_tet[t]) = new_vertex_idx;
      if(new_tet(i,tet.size(2) + t) == edge.second) new_tet(i,tet.size(2) + t)= new_vertex_idx;
    }
  }

  // copy the zyz
  new_zyz.resize(3,zyz.size(2) + adjacent_tet.size());
  new_zyz(colon(),colon(0,zyz.size(2) - 1)) = zyz;
  for(size_t t = 0 ; t < adjacent_tet.size(); ++t){
    new_zyz(colon(),zyz.size(2) + t ) = zyz(colon(),adjacent_tet[t]);
  }
  return 0;

}

int split_tet_at_zigzag(const matrixst &tet,
                        const matrixd &node,
                        const matrixd &zyz,
                        const vector<deque<pair<size_t,size_t> > > &singularity_chain,
                        matrixst &new_tet,
                        matrixd &new_node,
                        matrixd &new_zyz)
{

  matrixst new_tet_tmp ;
  matrixd new_node_tmp;
  matrixd new_zyz_tmp ;

  new_tet = tet;
  new_node = node;
  new_zyz = zyz;

  for(size_t t = 0; t < singularity_chain.size(); ++t){ // for each singularity chain
    for(size_t i = 0; i < singularity_chain[t].size() - 1; ++i){
      if(!is_edges_in_one_tet(singularity_chain[t][i],singularity_chain[t][i+1],new_tet)) continue; // not in one tet
      vector<size_t> adjacent_tet;
      assert(singularity_chain[t][i].second == singularity_chain[t][i+1].first);
      pair<size_t,size_t> edge(singularity_chain[t][i].first,singularity_chain[t][i+1].second);
      find_adjacent_tet(edge,new_tet,adjacent_tet);

      split_tet_and_copy_zyz_node(edge,adjacent_tet,new_tet,new_node,new_zyz,new_tet_tmp,new_node_tmp,new_zyz_tmp);
      swap(new_tet_tmp,new_tet);
      swap(new_node_tmp,new_node);
      swap(new_zyz_tmp,new_zyz);
    }
  }
  return 0;
}

int modidy_edge_type_by_change_face_type(
    const pair<size_t,size_t> & edge,
    const pair<size_t,size_t> & tet_pair,
    const matrixst & tet,
    const jtf::mesh::one_ring_tet_at_edge & ortae,
    boost::unordered_map<std::pair<size_t,size_t>,size_t> & inner_face_jump_type)
{
  assert(tet_pair.first < tet.size(2) && tet_pair.second < tet.size(2));
  typedef jtf::mesh::one_ring_tet_at_edge::e2tet_type::const_iterator oecit;
  oecit it = ortae.e2t_.find(edge);
  if(it == ortae.e2t_.end()) it = ortae.e2t_.find(make_pair(edge.second,
                                                            edge.first));
  if(it == ortae.e2t_.end()) {
    cerr << "# [error] can not find edge in one_ring_tet_at_edge <"
         <<  edge.first << "," << edge.second << ">." << endl;
    return __LINE__;
  }

  const vector<size_t> &loop = it->second;
  if(find(loop.begin(),loop.end(),-1) != loop.end()){
    cerr << "# [error] given edge is surface edge." << endl;
    return __LINE__;
  }

  assert(loop.front() == loop.back());
  deque<size_t> orient_loop(loop.size() - 1);
  copy(loop.begin(),loop.end()-1,orient_loop.begin());


  if(find(orient_loop.begin(),orient_loop.end(),tet_pair.first) == orient_loop.end() ||
     find(orient_loop.begin(),orient_loop.end(),tet_pair.second) == orient_loop.end()){
    cerr << "# [error] can not find tet pair in around tet loop." << endl;
    return __LINE__;
  }

  //while(loop.back())
  while((orient_loop.front() != tet_pair.first &&
         orient_loop.front() != tet_pair.second) ||
        (orient_loop.back() != tet_pair.first &&
         orient_loop.back() != tet_pair.second))
  {
    orient_loop.push_front(orient_loop.back());
    orient_loop.pop_back();
  }

  assert((orient_loop.back() == tet_pair.first ||
          orient_loop.back() == tet_pair.second) &&
         (orient_loop.front() == tet_pair.first ||
          orient_loop.front() == tet_pair.second));

  typedef boost::unordered_map<pair<size_t,size_t>,size_t>::const_iterator mpscit;
  typedef boost::unordered_map<pair<size_t,size_t>,size_t>::iterator mpsit;
  matrixd rot = eye<double>(3);
  // cerr << "# edge <" << edge.first << " " << edge.second << endl;
  for(size_t t = 0; t < orient_loop.size(); ++t){
    mpscit mit = inner_face_jump_type.find(
                   make_pair(orient_loop[t],
                             orient_loop[(t+1)%orient_loop.size()]));
    if(mit == inner_face_jump_type.end()) {
      //  cerr << "9 ";
      continue;
    }
    // cerr << mit->second << " ";
    rot = temp(rot * type_transition2(mit->second));
  }
  //cerr << endl;

  size_t compensate_face_type = type_transition1(trans(rot));
  assert(is_regular_type(compensate_face_type));
  if(orient_loop.back() == tet_pair.second){ // needs to trans face type
    compensate_face_type = type_transition1(rot);
  }
  mpsit mit = inner_face_jump_type.find(tet_pair);
  if(mit == inner_face_jump_type.end()){
    inner_face_jump_type[tet_pair] = compensate_face_type;
    inner_face_jump_type[make_pair(tet_pair.second,
                                   tet_pair.first)]
        = type_transition1(trans(type_transition2(compensate_face_type)));
  }else{
    mit->second = type_transition1(type_transition2(mit->second)*
                                   type_transition2(compensate_face_type));
    inner_face_jump_type[make_pair(tet_pair.second,
                                   tet_pair.first)]
        = type_transition1(trans(type_transition2(mit->second)));
  }

#if 1 // check
  {
    matrixd rot = eye<double>(3);
    //cerr << "# after modify edge <" << edge.first << " " << edge.second << endl;
    for(size_t t = 0; t < orient_loop.size(); ++t){
      mpscit mit = inner_face_jump_type.find(
                     make_pair(orient_loop[t],
                               orient_loop[(t+1)%orient_loop.size()]));
      if(mit == inner_face_jump_type.end()) {
        //  cerr << "9 ";
        continue;
      }
      // cerr << mit->second << " ";
      rot = temp(rot * type_transition2(mit->second));
    }
    //cerr << endl;
    if(!is_trivial_type(type_transition1(rot))){
      cerr << "# [error] modification error." << endl;
      return __LINE__;
    }
  }
#endif

  return 0;
}

int extract_singulairty_edges(
    const jtf::mesh::one_ring_tet_at_edge &ortae,
    const map<pair<size_t,size_t>,size_t> &inner_face_jump_type,
    vector<size_t> &singularity_edges)
{
  typedef jtf::mesh::one_ring_tet_at_edge::e2tet_type::const_iterator oecit;
  typedef map<pair<size_t,size_t>,size_t>::const_iterator mcit;
  matrixd rot = eye<double>(3);
  for(oecit it = ortae.e2t_.begin(); it != ortae.e2t_.end(); ++it){
    const vector<size_t> &chain = it->second;
    if(chain.front() == -1) continue;
    if(chain.front() != chain.back()) continue;
    rot = eye<double>(3);
    for(size_t t = 0; t < chain.size() - 1; ++t){
      mcit mit = inner_face_jump_type.find(make_pair(chain[t],chain[t+1]));
      if(mit == inner_face_jump_type.end()) continue;
      rot = temp(rot * type_transition2(mit->second));
    }
    if(fabs(norm(rot - eye<double>(3))) > 1e-8){
      singularity_edges.push_back(it->first.first);
      singularity_edges.push_back(it->first.second);
    }
  }
  return 0;
}

int modify_face_jump_type_at_given_tets_edge(
    const size_t tet_idx_a,
    const size_t tet_idx_b,
    boost::unordered_map<pair<size_t,size_t>,size_t> &inner_face_jump_type,
    const vector<size_t> &edge0_tet)
{
  for(size_t t = 0; t < edge0_tet.size() - 1; ++t)
  {
    if((edge0_tet[t] == tet_idx_a && (edge0_tet[t+1] == tet_idx_b))
       || (edge0_tet[t] == tet_idx_b && (edge0_tet[t+1] == tet_idx_a))){ // match the two tets
      vector<size_t> reorder_tets;
      reorder_tets.reserve(edge0_tet.size());
      for(size_t i = t; i < edge0_tet.size() - 1; ++i) reorder_tets.push_back(edge0_tet[i]);
      for(size_t i = 0; i < t + 1; ++i) reorder_tets.push_back(edge0_tet[i]);
      modify_face_jump_type_to_remove_singularity_edge(edge0_tet[t],edge0_tet[t+1],inner_face_jump_type,reorder_tets);
      return 0;
    }
  }
  return __LINE__;
}


int modify_face_jump_type_to_remove_singularity_edge(
    const size_t tet_idx_a,
    const size_t tet_idx_b,
    boost::unordered_map<pair<size_t,size_t>,size_t> &inner_face_jump_type,
    const vector<size_t> &edge0_tets)
{
  typedef boost::unordered_map<pair<size_t,size_t>,size_t>::const_iterator mcit;
  matrixd rot = eye<double>(3);
  matrixd jump_ab = eye<double>(3);

  for(size_t t = 0; t < edge0_tets.size() - 1; ++t){
    mcit it = inner_face_jump_type.find(make_pair(edge0_tets[t],edge0_tets[t+1]));
    if(t == 0){
      if(it != inner_face_jump_type.end())
        jump_ab = type_transition2(it->second);
    }
    if(it == inner_face_jump_type.end()) continue; // stands for Identity jump
    rot = temp(rot * type_transition2(it->second));
  }

  matrixd test_ab = jump_ab;
  jump_ab = temp(trans(rot) * jump_ab);
  inner_face_jump_type[make_pair(edge0_tets[0],edge0_tets[1])] = type_transition1(jump_ab);
  inner_face_jump_type[make_pair(edge0_tets[1],edge0_tets[0])] = type_transition1(trans(jump_ab));

#if 1
  matrixd rot_ = eye<double>(3);
  for(size_t t = 0; t < edge0_tets.size() - 1; ++t){
    mcit it = inner_face_jump_type.find(make_pair(edge0_tets[t],edge0_tets[t+1]));
    if(it == inner_face_jump_type.end()) continue; // means this face type equals identity
    rot_ = temp(rot_ * type_transition2(it->second));
  }
  if(fabs(norm(rot_ - eye<double>(3))) > 1e-8){
    matrixd original_rot = test_ab;
    cerr << "# original type list ";
    cerr << type_transition1(test_ab) << " ";
    for(size_t i = 1 ; i < edge0_tets.size() - 1; ++i)
    {
      original_rot = temp(original_rot * type_transition2(inner_face_jump_type[make_pair(edge0_tets[i],edge0_tets[i+1])]));
      //cerr << type_transition2(inner_face_jump_type[make_pair(edge0_tets[i],edge0_tets[i+1])]) << " ";
    }
    cerr << endl;
    cerr << "# original_rot " << original_rot << endl;

    cerr << "# new type list " ;
    matrixd new_rot = eye<double>(3);
    for(size_t i = 0 ; i < edge0_tets.size() - 1; ++i)
    {
      new_rot = temp(new_rot * type_transition2(inner_face_jump_type[make_pair(edge0_tets[i],edge0_tets[i+1])]));
      //cerr << type_transition1(inner_face_jump_type[make_pair(edge0_tets[i],edge0_tets[i+1])]) << " ";
    }
    cerr << endl;
    cerr << "# new_rot " << new_rot << endl;
    cerr << "# error, this modification is not correct on the modified edge ." << endl;
  }
#endif
  return 0;
}

int find_triangle_fan_path(
    const size_t begin_point_of_path,
    const matrixst &tet,
    const matrixd &node,
    const std::deque<size_t> &black_chain,
    const jtf::mesh::one_ring_tet_at_edge & ortae,
    const jtf::mesh::face2tet_adjacent & fa,
    std::deque<size_t> &path)
{
  path.clear();
  path.push_back(begin_point_of_path);
  vector<size_t> one_ring_point;
  vector<pair<double,size_t> > dis_points;

  typedef  deque<size_t>::const_iterator dcit;
  dcit iter = black_chain.begin();
  while(iter != black_chain.end() - 1){
    dis_points.clear();
    pair<size_t,size_t> edge(*iter, path.back());
    const size_t &next_point = *(iter+1);
    ortae.get_one_ring_points_of_edge(tet, edge, one_ring_point);
    for(size_t t = 0; t < one_ring_point.size(); ++t){
      // have not met such point
      if(find(path.begin(), path.end(), one_ring_point[t]) == path.end() &&
         find(black_chain.begin(), iter, one_ring_point[t])== iter){
        dis_points.push_back(
              make_pair(norm(node(colon(), one_ring_point[t]) -
                             node(colon(), next_point)),
                        t));
      }
    }
    sort(dis_points.begin(), dis_points.end());
    size_t pi = 0;
    //if(iter + 1 != black_chain.end()){
    for(; pi < dis_points.size(); ++pi){
      size_t face_idx = fa.get_face_idx(next_point,*iter, path.back());
      if(face_idx != -1) {
        cerr << "# [info] face idx " << face_idx << endl;
        break;
      }
    }
    if(pi < dis_points.size()){
      ++iter;
    }else
      path.push_back(one_ring_point[dis_points.front().second]);
  }
  path.push_back(black_chain.back());
  return 0;
}

int find_shortest_path_new(const size_t begin_,
                           const size_t end_,
                           const matrixst &tet,
                           const matrixd &node,
                           const deque<size_t> &black_chain,
                           deque<size_t> &path)
{
  using namespace boost;
  typedef adjacency_list < vecS, vecS, undirectedS,
      no_property, property < edge_weight_t, double > > Graph;
  typedef graph_traits<Graph>::vertex_descriptor Vertex;
  typedef graph_traits < Graph >::vertex_descriptor vertex_descriptor;

  set<size_t> around_tets;
  { // collect all the tets around the black_chain
    //for(size)
    set<size_t> black_chain_vertex;
    for(size_t t = 0 ; t < black_chain.size(); ++t){
      black_chain_vertex.insert(black_chain[t]);
    }

    // find all around tet
    for(set<size_t>::const_iterator sci = black_chain_vertex.begin();
        sci != black_chain_vertex.end(); ++sci)
      for(size_t t = 0; t < tet.size(2); ++t){
        if(find(tet(colon(),t).begin(),
                tet(colon(),t).end(),
                *sci) != tet(colon(),t).end())
          around_tets.insert(t);
      }
    // find all vertex in graph
    set<size_t> vertex_set;
    vertex_set.insert(end_);
    for(set<size_t>::const_iterator sci = around_tets.begin();
        sci != around_tets.end(); ++sci){
      for(size_t t = 0; t < 4; ++t){
        if(find(black_chain.begin(),black_chain.end(),tet(t,*sci)) == black_chain.end())
          vertex_set.insert(tet(t,*sci));
      }
    }

    vector<size_t> vertex_name(vertex_set.size());
    copy(vertex_set.begin(),vertex_set.end(),vertex_name.begin());

    map<size_t,size_t> smp;
    for(size_t t = 0; t < vertex_name.size(); ++t) smp[vertex_name[t]] = t;

    // add edge into graph
    set<pair<size_t,size_t> > edges;

    for(set<size_t>::const_iterator sci = around_tets.begin();
        sci != around_tets.end(); ++sci){
      for(size_t t = 0; t < 4; ++t) {
        size_t v0 = tet(t,*sci);
        size_t v1 = tet((t + 1) %4,*sci);
        if((v0 == end_ && (find(black_chain.begin(),black_chain.end(),v1) == black_chain.end()))
           || (v1 == end_ && (find(black_chain.begin(),black_chain.end(),v0) == black_chain.end()))
           || (find(black_chain.begin(),black_chain.end(),v0) == black_chain.end()
               &&  find(black_chain.begin(),black_chain.end(),v1) == black_chain.end())) // both v0 and v1 are not on black_chain
        {
          if(v0 > v1) swap(v0,v1);
          edges.insert(make_pair(smp[v0],smp[v1]));
        }
      }
    }

    matrixd weights = ones<double>(edges.size(),1);

    Graph g(edges.begin(),edges.end(), &weights[0],vertex_set.size());
    std::vector<size_t> d(num_vertices(g));
    vertex_descriptor s = vertex(smp[begin_], g);
    std::vector<vertex_descriptor> p(num_vertices(g));
    dijkstra_shortest_paths(g, s, predecessor_map(&p[0]).distance_map(&d[0]));

    //std::cout << "distances and parents:" << std::endl;
    //graph_traits < Graph >::vertex_iterator vi;//, vend;
    //tie(vi,vend) = vertices(g);
    //graph_traits < Graph >::vertex_descriptor s_end = vertex(end_,g);

#if 0 // modified by jtf 2012-05-03 can not be complied in boost-4.6

#endif
    vertex_descriptor vi = vertex(smp[end_],g);
    while(vi != s)
    {
      path.push_front(vertex_name[vi]);
      vi = p[vi];
      //cerr << vertex_name[*vi] << endl;
    }
    path.push_front(begin_);

    //    *vi = vertex(smp[end_],g);

    //    //tie(vi, vend) = vertices(g);
    //    while(vertex_name[*vi] != begin_)
    //    {
    //      path.push_front(vertex_name[*vi]);
    //      vi = p[*vi];
    //      //cerr << vertex_name[*vi] << endl;
    //    }
    //    path.push_front(begin_);



  }
  return 0;
}

#if 0
int find_shortest_path(const size_t &begin_,
                       const size_t &end_,
                       const matrixst & tet,
                       const matrixd &node,
                       const deque<size_t>  &black_chain,
                       map<pair<size_t,size_t>, bool > &map_changed_edge)
{
  using namespace boost;
  typedef adjacency_list < vecS, vecS, undirectedS,
      no_property, property < edge_weight_t, double > > Graph;
  typedef graph_traits<Graph>::vertex_descriptor Vertex;
  typedef graph_traits < Graph >::vertex_descriptor vertex_descriptor;

  size_t black_begin;
  if(end_ == black_chain.front())
    black_begin = black_chain.back();
  else
    black_begin = black_chain.front();
  map_changed_edge.insert(make_pair(make_pair(begin_,black_begin),false));

  set<size_t> around_tets;
  { // collect all the tets around the black_chain
    //for(size)
    set<size_t> black_chain_vertex;
    for(size_t t = 0 ; t < black_chain.size(); ++t){
      black_chain_vertex.insert(black_chain[t]);
    }

    // find all around tet
    for(set<size_t>::const_iterator sci = black_chain_vertex.begin();
        sci != black_chain_vertex.end(); ++sci)
      for(size_t t = 0; t < tet.size(2); ++t){
        if(find(tet(colon(),t).begin(),
                tet(colon(),t).end(),
                *sci) != tet(colon(),t).end())
          around_tets.insert(t);
      }
    // find all vertex in graph
    set<size_t> vertex_set;
    vertex_set.insert(end_);
    for(set<size_t>::const_iterator sci = around_tets.begin();
        sci != around_tets.end(); ++sci){
      for(size_t t = 0; t < 4; ++t){
        if(find(black_chain.begin(),black_chain.end(),tet(t,*sci)) == black_chain.end())
          vertex_set.insert(tet(t,*sci));
      }
    }

    vector<size_t> vertex_name(vertex_set.size());
    copy(vertex_set.begin(),vertex_set.end(),vertex_name.begin());

    map<size_t,size_t> smp;
    for(size_t t = 0; t < vertex_name.size(); ++t) smp[vertex_name[t]] = t;


    // add edge into graph
    set<pair<size_t,size_t> > edges;

    for(set<size_t>::const_iterator sci = around_tets.begin();
        sci != around_tets.end(); ++sci){
      for(size_t t = 0; t < 4; ++t) {
        size_t v0 = tet(t,*sci);
        size_t v1 = tet((t + 1) %4,*sci);
        if((v0 == end_ && (find(black_chain.begin(),black_chain.end(),v1) == black_chain.end()))
           || (v1 == end_ && (find(black_chain.begin(),black_chain.end(),v0) == black_chain.end()))
           || (find(black_chain.begin(),black_chain.end(),v0) == black_chain.end()
               &&  find(black_chain.begin(),black_chain.end(),v1) == black_chain.end())) // both v0 and v1 are not on black_chain
        {
          if(v0 > v1) swap(v0,v1);
          edges.insert(make_pair(smp[v0],smp[v1]));
        }
      }
    }

    matrixd weights = ones<double>(edges.size(),1);
    //matrixd weights = zeros<double>(edge)
    //        vector<pair<size_t,size_t> > edges_v(edges.size());
    //        copy(edges.begin(),edges.end(),edges_v.begin());
    //        for(size_t t = 0; t < edges_v.size(); ++t)
    //            weights[t] = norm(node(colon(),vertex_name[edges_v[t].first])
    //                            - node(colon(),vertex_name[edges_v[t].second]));
    Graph g(edges.begin(),edges.end(), &weights[0],vertex_set.size());
    std::vector<size_t> d(num_vertices(g));
    vertex_descriptor s = vertex(smp[begin_], g);
    std::vector<vertex_descriptor> p(num_vertices(g));
    dijkstra_shortest_paths(g, s, predecessor_map(&p[0]).distance_map(&d[0]));

    //std::cout << "distances and parents:" << std::endl;
    //graph_traits < Graph >::vertex_iterator vi;//, vend;
    //tie(vi,vend) = vertices(g);
    //graph_traits < Graph >::vertex_descriptor s_end = vertex(end_,g);

#if 0
    for (tie(vi, vend) = vertices(g); vi != vend; ++vi)
    {
      std::cout << "distance(" << vertex_name[*vi] << ") = " << d[*vi] << ", ";
      std::cout << "parent(" << vertex_name[*vi] << ") = " << vertex_name[p[*vi]] << std::endl;
    }
#endif

    vertex_descriptor vi = vertex(smp[end_],g);
    //tie(vi, vend) = vertices(g);
    while(vertex_name[*vi] != begin_)
    {
      map_changed_edge.insert(make_pair(make_pair(vertex_name[p[*vi]],vertex_name[*vi]),true));
      vi = p[*vi];
      //cerr << vertex_name[*vi] << endl;
    }
  }
  return 0;
}
#endif

size_t is_zigzag(const pair<size_t,size_t> &edge0,
                 const pair<size_t,size_t> &edge1,
                 const jtf::mesh::face2tet_adjacent &fa)
{
  set<size_t> vertex;
  vertex.insert(edge0.first);
  vertex.insert(edge0.second);
  vertex.insert(edge1.first);
  vertex.insert(edge1.second);

  if(vertex.size() != 3) return -1;

  vector<size_t> face(3,0);
  copy(vertex.begin(),vertex.end(),face.begin());
  const size_t face_idx = fa.get_face_idx(&face[0]);

  if( face_idx > fa.faces_.size())
    return -1;
  else
    return face_idx;
}


bool is_zigzag_edge(const std::pair<size_t,size_t> &edge0,
                    const std::pair<size_t,size_t> &edge1,
                    const size_t * tet_array,
                    size_t tet_num)
{
  itr_matrix<const size_t *> tet(4, tet_num, tet_array);
  set<size_t> points;
  points.insert(edge0.first);
  points.insert(edge0.second);
  points.insert(edge1.first);
  points.insert(edge1.second);
  if(points.size() != 3) return false;

  size_t count_extra_point_num = 0;
  for(size_t t = 0; t < tet.size(2); ++t){
    count_extra_point_num = 0;
    for(set<size_t>::const_iterator scit = points.begin();
        scit != points.end(); ++scit){
      if(    tet(0,t) != *scit
             && tet(1,t) != *scit
             && tet(2,t) != *scit
             && tet(3,t) != *scit)
        ++count_extra_point_num;

      if(count_extra_point_num == 1) // three point of two edges must be in the tet, or two edges do not belong to this tet
        break;
    }
    if(count_extra_point_num == 0)
      return true;
  }
  return false;
}

std::pair<size_t,size_t> get_adj_tets(const size_t point_a,
                                      const size_t point_b,
                                      const size_t point_c,
                                      const size_t *tet_array,
                                      const size_t tet_num)
{
  assert(tet_array);
  itr_matrix<const size_t *>  tet(4, tet_num, tet_array);
  size_t face_point_num = 0;
  vector<size_t> tets_pair;
  for(size_t t = 0; t < tet.size(2); ++t){
    face_point_num = 0;
    for(size_t i = 0; i < tet.size(1); ++i){
      if(    point_a == tet(i,t)
             || point_b == tet(i,t)
             || point_c == tet(i,t))
        ++face_point_num;

      if(face_point_num == 3){
        tets_pair.push_back(t);

        if(tets_pair.size() == 2)
          return make_pair(tets_pair.front(),tets_pair.back());

        break;
      }
    }
  }
  while(tets_pair.size() < 2) tets_pair.push_back(-1);
  return make_pair(tets_pair.front(), tets_pair.back());
}

// this function only define the chain which has all vertex one ring near surface as near surface chain.
bool is_one_ring_near_surface_chain(const deque<pair<size_t,size_t> > &chain,
                                    const matrixst &tet,
                                    const matrixst &outside_face)
{
  const size_t &begin_vertex = chain.front().first;
  const size_t &end_vertex = chain.back().second;
  if(!is_outside_vertex(begin_vertex,outside_face)
     || !is_outside_vertex(end_vertex,outside_face) )
    return false;

  if(chain.size() < 4) return true; // s0 --> a -- > b --> s1

  for(size_t t = 1; t < chain.size() - 2; ++t)
  {
    if(!is_one_ring_near_surface(chain[t].second,tet,outside_face))
      return false;
  }
  return true;
}

bool is_one_ring_near_surface(const size_t vertex,
                              const matrixst &tet,
                              const matrixst &outside_face)
{
  set<size_t> one_ring_vertex;
  for(size_t t = 0; t < tet.size(2); ++t){
    if(find(tet(colon(),t).begin(),tet(colon(),t).end(),t) != tet(colon(),t).end()){
      for(size_t i = 0; i < 4; ++i){
        if(tet(i,t) != vertex)
          one_ring_vertex.insert(tet(i,t));
      }
    }
  }

  for(set<size_t>::const_iterator sci = one_ring_vertex.begin(); sci != one_ring_vertex.end(); ++sci){
    if(find(outside_face.begin(),outside_face.end(),*sci) != outside_face.end())
      return true;
  }
  return false;
}

int get_surface_normal_align_type(
    vector<size_t> & surface_normal_align_type,
    const matrixst &tet,
    const matrixd &node,
    const matrixst &outside_face,
    const matrixst &outside_face_idx,
    const jtf::mesh::face2tet_adjacent &fa,
    const vector<matrixd > &frame_array)
{
  matrixd normal(3,1);
  surface_normal_align_type.resize(outside_face_idx.size(), -1);
  //surface_normal_align_type = ones<size_t>(1,outside_face_idx.size()) * -1;
  //  vector<pair<double,size_t> > surface_normal_type(6);

  for(size_t i = 0; i < outside_face_idx.size(); ++i){
    const pair<size_t,size_t> & tets = fa.query(&outside_face(0,i));
    const size_t &tet_idx = (tets.first == -1?tets.second:tets.first);
    //    if(tet_idx == 2571)
    //      cerr << "pause";
    jtf::tetmesh::calculate_face_normal(tet,node,tet_idx,outside_face(colon(),i),normal);
    double diff = 0;
    size_t idx = -1;
    for(size_t j = 0; j < 6; ++j){
      const double current_diff = dot(normal,
                                      frame_array[tet_idx](colon(),j/2)*(j%2==0?1:-1));
      if(current_diff > diff){
        diff = current_diff;
        idx = j;
      }
    }
    if(idx == -1)
      cerr << "# [error] calculate surface normal align error." << endl;
    surface_normal_align_type[i] = idx;
  }
  return 0;
}

int get_surface_normal_align_type(
    boost::unordered_map<size_t,size_t> &surface_normal_align_type,
    const matrixst &tet,
    const matrixd &node,
    const matrixst &outside_face,
    const matrixst &outside_face_idx,
    const jtf::mesh::face2tet_adjacent &fa,
    const zjucad::matrix::matrix<matrixd>& frame)
{
  matrixd normal(3,1);

  for(size_t i = 0; i < outside_face_idx.size(); ++i){
    const pair<size_t,size_t> & tets = fa.query(&outside_face(0,i));
    const size_t &tet_idx = (tets.first == -1?tets.second:tets.first);
    jtf::tetmesh::calculate_face_normal(tet,node,tet_idx,outside_face(colon(),i),normal);
    double diff = 0;
    size_t idx = -1;
    for(size_t j = 0; j < 6; ++j){
      const double current_diff =
          dot(normal, frame[tet_idx](colon(),j/2)*(j%2==0?1:-1));
      if(current_diff > diff){
        diff = current_diff;
        idx = j;
      }
    }
    if(idx == -1)
      cerr << "# [error] calculate surface normal align error." << endl;
    surface_normal_align_type[outside_face_idx[i]] = idx;
  }
  return 0;
}

int get_surface_patch_type_mat(matrixst &surface_patch_idx,
                               const matrixst &surface_normal_align_type,
                               const matrixst & tet,
                               const matrixd & node,
                               const matrixst & outside_face,
                               const matrixst & outside_face_idx,
                               const jtf::mesh::face2tet_adjacent &fa_new)
{
  surface_patch_idx.resize(1,outside_face_idx.size());
  //vector<size_t> surface_normal_align_type;
  //  surface_normal_align_type.clear();
  //  get_surface_normal_align_type(surface_normal_align_type,tet,node,
  //                                outside_face,outside_face_idx,fa_new,frame_array);

  // u,v,w axis is not unique to diff
  unique_ptr<edge2cell_adjacent> ea(edge2cell_adjacent::create(outside_face));
  vector<bool> is_visited(outside_face_idx.size(),false);

  vector<vector<size_t> > patch_list;
  stack<size_t> possible_face;
  possible_face.push(0);
  // TODO: need to speed up
  while(find(is_visited.begin(),is_visited.end(),false) != is_visited.end())
  {
    vector<size_t> one_patch;
    while(!possible_face.empty()){
      size_t one_face_idx = possible_face.top();
      possible_face.pop();
      one_patch.push_back(one_face_idx);
      is_visited[one_face_idx] = true;
      for(size_t t = 0; t <  3; ++t){ // for each adjacent face
        size_t other_face = -1;
        const pair<size_t,size_t> face_pair =
            ea->query(outside_face(t,one_face_idx),outside_face((t+1)%3,one_face_idx));
        if(face_pair.first != one_face_idx) other_face = face_pair.first;
        else other_face = face_pair.second;
        if(other_face == -1) continue;
        if(is_visited[other_face]) continue;
        if(surface_normal_align_type[other_face] == surface_normal_align_type[one_face_idx])
          possible_face.push(other_face);
      }
    }
    patch_list.push_back(one_patch);
    vector<bool>::const_iterator vcit = find(is_visited.begin(),is_visited.end(),false);
    if(vcit == is_visited.end()) break;
    possible_face.push(static_cast<size_t>(vcit-is_visited.begin()));
  }

  for(size_t t = 0; t < patch_list.size(); ++t){
    for(size_t i = 0; i < patch_list[t].size(); ++i){
      surface_patch_idx[patch_list[t][i]] = t;
    }
  }

  {// for visualization
    ofstream ofs("surface_patch_type.vtk");
    tri2vtk(ofs,&node[0],node.size(2),&outside_face[0],outside_face.size(2));
    cell_data(ofs,&surface_patch_idx[0],surface_patch_idx.size(),"patch_idx");

    ofstream ofs_n("surface_normal_align_type.vtk");
    tri2vtk(ofs_n,&node[0],node.size(2),&outside_face[0],outside_face.size(2));
    cell_data(ofs_n,&surface_normal_align_type[0],surface_normal_align_type.size(),"normal_align_type");
  }

  return 0;
}


int get_surface_patch_type(vector<size_t> &surface_patch_idx,
                           const vector<size_t> &surface_normal_align_type,
                           const matrixst & tet,
                           const matrixd & node,
                           const matrixst & outside_face,
                           const matrixst & outside_face_idx,
                           const jtf::mesh::face2tet_adjacent &fa_new)
{
  surface_patch_idx.resize(outside_face_idx.size());
  //vector<size_t> surface_normal_align_type;
  //  surface_normal_align_type.clear();
  //  get_surface_normal_align_type(surface_normal_align_type,tet,node,
  //                                outside_face,outside_face_idx,fa_new,frame_array);

  // u,v,w axis is not unique to diff
  unique_ptr<edge2cell_adjacent> ea(edge2cell_adjacent::create(outside_face));
  vector<bool> is_visited(outside_face_idx.size(),false);

  vector<vector<size_t> > patch_list;
  stack<size_t> possible_face;
  possible_face.push(0);
  // TODO: need to speed up
  while(find(is_visited.begin(),is_visited.end(),false) != is_visited.end())
  {
    vector<size_t> one_patch;
    while(!possible_face.empty()){
      size_t one_face_idx = possible_face.top();
      possible_face.pop();
      one_patch.push_back(one_face_idx);
      is_visited[one_face_idx] = true;
      for(size_t t = 0; t <  3; ++t){ // for each adjacent face
        size_t other_face = -1;
        const pair<size_t,size_t> face_pair =
            ea->query(outside_face(t,one_face_idx),outside_face((t+1)%3,one_face_idx));
        if(face_pair.first != one_face_idx) other_face = face_pair.first;
        else other_face = face_pair.second;
        if(other_face == -1) continue;
        if(is_visited[other_face]) continue;
        if(surface_normal_align_type[other_face] == surface_normal_align_type[one_face_idx])
          possible_face.push(other_face);
      }
    }
    patch_list.push_back(one_patch);
    vector<bool>::const_iterator vcit = find(is_visited.begin(),is_visited.end(),false);
    if(vcit == is_visited.end()) break;
    possible_face.push(static_cast<size_t>(vcit-is_visited.begin()));
  }

  for(size_t t = 0; t < patch_list.size(); ++t){
    for(size_t i = 0; i < patch_list[t].size(); ++i){
      surface_patch_idx[patch_list[t][i]] = t;
    }
  }

  {// for visualization
    ofstream ofs("surface_patch_type.vtk");
    tri2vtk(ofs,&node[0],node.size(2),&outside_face[0],outside_face.size(2));
    cell_data(ofs,&surface_patch_idx[0],surface_patch_idx.size(),"patch_idx");

    ofstream ofs_n("surface_normal_align_type.vtk");
    tri2vtk(ofs_n,&node[0],node.size(2),&outside_face[0],outside_face.size(2));
    cell_data(ofs_n,&surface_normal_align_type[0],surface_normal_align_type.size(),"normal_align_type");
  }

  return 0;
}

int find_adjacent_face_of_point(const matrixst &outside_face,
                                const edge2cell_adjacent &ea,
                                map<size_t,set<size_t> > & point_adjacent_tri)
{
  point_adjacent_tri.clear();

  for(size_t t = 0; t < ea.edge2cell_.size(); ++t){
    const pair<size_t,size_t> & edge = ea.edges_[t];
    const pair<size_t,size_t> & faces = ea.edge2cell_[t];
    if(faces.first != -1) {
      point_adjacent_tri[edge.first].insert(faces.first);
      point_adjacent_tri[edge.second].insert(faces.first);
    }
    if(faces.second != -1){
      point_adjacent_tri[edge.first].insert(faces.second);
      point_adjacent_tri[edge.second].insert(faces.second);
    }
  }
  return 0;
}

bool is_compatible_situation(const size_t v0,const size_t v1,
                             const std::set<size_t> &adj_tri_0,
                             const std::set<size_t> & adj_tri_1,
                             const std::vector<size_t> &surface_path_idx,
                             const std::vector<size_t> &surface_normal_align_uvw_idx)
{
  set<size_t> patch0,patch1;

  for(set<size_t>::const_iterator scit = adj_tri_0.begin();
      scit != adj_tri_0.end(); ++scit)
    patch0.insert(surface_path_idx[*scit]);

  for(set<size_t>::const_iterator scit = adj_tri_1.begin();
      scit != adj_tri_1.end(); ++scit)
    patch1.insert(surface_path_idx[*scit]);

  vector<size_t> intersection(patch0.size() + patch1.size());
  vector<size_t>::iterator it =
      find_intersection_set_with_sorted_input(patch0.begin(),patch0.end(),
                                              patch1.begin(),patch1.end(),
                                              intersection.begin());
  if(it == intersection.begin()) {
    //    if(surface_normal_align_uvw_idx.size() == surface_path_idx.size())
    //    {// test if the both ends of this singularity touch three kinds of patches,
    //      // which must be removed, or will result in degeneration
    //      size_t around_uvw[3] = {0,0,0};
    //      for(set<size_t>::const_iterator scit = adj_tri_0.begin();
    //          scit != adj_tri_0.end(); ++scit){
    //        around_uvw[surface_normal_align_uvw_idx[*scit]/2] = 1;
    //      }
    //      for(set<size_t>::const_iterator scit = adj_tri_1.begin();
    //          scit != adj_tri_1.end(); ++scit){
    //        around_uvw[surface_normal_align_uvw_idx[*scit]/2] = 1;
    //      }
    //      if(around_uvw[0] + around_uvw[1] + around_uvw[2] == 3) return true;
    //    }
    return false;
  }
  return true;
}

int find_near_surface_singularity(const matrixst &tet,
                                  const matrixd &node_array,
                                  const matrixst & outside_face,
                                  const vector<deque<pair<size_t,size_t> > > &chain_list,
                                  const vector<size_t> &surface_path_idx,
                                  const vector<size_t> & surface_path_align_uvw_idx, // this record which axis is this face most like
                                  vector<size_t> &near_surface_singularity_list)
{
  unique_ptr<edge2cell_adjacent> ea(edge2cell_adjacent::create(outside_face));

  map<size_t,set<size_t > > point_adjacent_tri;
  find_adjacent_face_of_point(outside_face,*ea,point_adjacent_tri);

#if 1 // check
  set<size_t> outside_vertex(outside_face.begin(),outside_face.end());
  if(outside_vertex.size() != point_adjacent_tri.size()) {
    cerr << "# [error] outside_point num is wrong." << endl;
    return __LINE__;
  }
#endif

  typedef map<size_t,set<size_t> >::const_iterator mcit;
  for(size_t t = 0; t < chain_list.size(); ++t){
    const deque<pair<size_t,size_t> > & one_chain = chain_list[t];
    if(one_chain.front().first == one_chain.back().second) continue; // inner loop

    mcit it0 = point_adjacent_tri.find(one_chain.front().first);
    mcit it1 = point_adjacent_tri.find(one_chain.back().second);

    if(it0 == point_adjacent_tri.end() || it1 == point_adjacent_tri.end()) // one or two end is inside
      continue;

    if(is_compatible_situation(one_chain.front().first,one_chain.back().second,
                               it0->second,it1->second,surface_path_idx,surface_path_align_uvw_idx))
      near_surface_singularity_list.push_back(t);
  }
  return 0;
}

///////////////////////////////////////////////////////////////////////////////////////////////
//template<DIRECTION DIRECTION_TYPE>
//vertex_connection<DIRECTION_TYPE>* vertex_connection<DIRECTION_TYPE>::create(const std::map<std::pair<size_t,size_t>,double> & edge_weight)
//{
//    unique_ptr<vertex_connection> vc(new vertex_connection);

//    if(vc->init(edge_weight)) // fail
//        return 0;
//    return vc.release();
//}

//template<DIRECTION DIRECTION_TYPE>
//int vertex_connection<DIRECTION_TYPE>::init(const std::map<std::pair<size_t,size_t>,double> & edge_weight)
//{
//    typedef map<pair<size_t,size_t>,double>::const_iterator mcit;
//    set<size_t> vertex_;
//    for(mcit it = edge_weight.begin(); it != edge_weight.end(); ++it){
//        vertex_.insert(it->first.first);
//        vertex_.insert(it->first.second);
//    }

//    vertex_name_.resize(vertex_.size());
//    copy(vertex_.begin(),vertex_.end(),vertex_name_.begin());

//    edges_.reserve(edge_weight.size());//(edge_weight.size());
//    weights_.reserve(edge_weight.size());

//    for(size_t t = 0; t < vertex_name_.size(); ++t) smp_[vertex_name_[t]] = t;

//    for(mcit it = edge_weight.begin(); it != edge_weight.end(); ++it){
//        edges_.push_back(make_pair(smp_[it->first.first],smp_[it->first.second]));
//        weights_.push_back(it->second);
//    }


//    graph_.reset(new Graph(edges_.begin(),edges_.end(),weights_.begin(),vertex_name_.size()));
//    return 0;
//}

//template<DIRECTION DIRECTION_TYPE>
//int vertex_connection<DIRECTION_TYPE>::get_shortest_path(const size_t begin_vertex,
//                                         const size_t end_vertex,
//                                         std::vector<size_t> &path) const
//{
//    typedef std::map<size_t,size_t>::const_iterator mcit;
//    const Graph & g = *graph_;

//    vector<double> distance(num_vertices(g));
//    mcit begin_vertex_name_iter = smp_.find(begin_vertex);
//    mcit end_vertex_name_iter = smp_.find(end_vertex);
//    if(begin_vertex_name_iter == smp_.end() || end_vertex_name_iter == smp_.end()) return __LINE__;

//    vertex_descriptor s = vertex(begin_vertex_name_iter->second, g);
//    std::vector<vertex_descriptor> p(num_vertices(g));
//    dijkstra_shortest_paths(g, s, predecessor_map(&p[0]).distance_map(&distance[0]));

//    typename graph_traits < Graph >::vertex_iterator vi;
//    vi = vertex(end_vertex_name_iter->second,g);


//    path.clear();

//    while(*vi != s){
//        if(p[*vi] == *vi){
//            return __LINE__;
//        }
//        path.push_back(vertex_name_[*vi]);
//        vi = p[*vi];
//    }
//    path.push_back(vertex_name_[*vi]);

//    reverse(path.begin(),path.end());
//    return 0;
//}

int find_one_ring_vertex_around_edge(
    const matrixst &tet,
    const std::pair<size_t,size_t>& edge,
    const jtf::mesh::one_ring_tet_at_edge & ortae,
    std::vector<size_t> &one_ring_vertex)
{
  typedef jtf::mesh::one_ring_tet_at_edge::e2tet_type::const_iterator oecit;
  oecit it = ortae.e2t_.find(edge);
  if(it == ortae.e2t_.end()) it = ortae.e2t_.find(make_pair(edge.second,edge.first));
  if(it == ortae.e2t_.end()) return __LINE__;
  const vector<size_t> & tets = it->second;
  set<size_t> around_vertex;
  for(size_t t = 0; t < tets.size() - 1; ++t){ // the tet front and back are the same
    if(tets[t] == -1) continue;
    for(size_t i = 0; i < 4; ++i){
      if(tet(i,tets[t]) != edge.first
         && tet(i,tets[t]) != edge.second )
      {
        around_vertex.insert(tet(i,tets[t]));
      }
    }
  }
  one_ring_vertex.resize(around_vertex.size());
  copy(around_vertex.begin(),around_vertex.end(),one_ring_vertex.begin());
  return 0;
}

int find_one_ring_edges_around_edge(
    const matrixst &tet,
    const matrixd &node,
    const std::pair<size_t,size_t>& around_edge,
    const jtf::mesh::one_ring_tet_at_edge & ortae,
    std::vector<std::pair<size_t,size_t> > &one_ring_edges)
{
  typedef jtf::mesh::one_ring_tet_at_edge::e2tet_type::const_iterator oecit;
  oecit it = ortae.e2t_.find(around_edge);
  bool is_reversed = false;
  if(it == ortae.e2t_.end()){
    it = ortae.e2t_.find(make_pair(around_edge.second,around_edge.first));
    is_reversed = true;
  }
  if(it == ortae.e2t_.end()) {
    cerr << "# [error] can not find this edge: " << around_edge.first << "-->" << around_edge.second << endl;
    return __LINE__;
  }
  vector<size_t>  tets_loop = it->second;
  if(is_reversed) reverse(tets_loop.begin(),tets_loop.end());

  vector<size_t> tmp_edge(2);
  if(tets_loop.size() == 2) { // it's a sharp edge
    const size_t tet_idx = (tets_loop.front() == -1)?tets_loop.back():tets_loop.front();
    adjust_cross_edges_right_hand_order_of_tet(tet(colon(),tet_idx),node,around_edge,tmp_edge);
    one_ring_edges.push_back(make_pair(tmp_edge[1],tmp_edge[0]));
  }else{
    for(size_t t = 0; t < tets_loop.size() - 1; ++t){
      if(tets_loop[t] == -1) continue;
      adjust_cross_edges_right_hand_order_of_tet(tet(colon(),tets_loop[t]),node,
                                                 around_edge,tmp_edge);
      one_ring_edges.push_back(make_pair(tmp_edge[1],tmp_edge[0]));
    }
#if 0 //check the order of around edges, must be linked right, this check is not correct, didn't consider open loop
    for(size_t t = 0; t < one_ring_edges.size(); ++t){
      if(one_ring_edges[t].second != one_ring_edges[(t + 1)%one_ring_edges.size()].first)
      {
        cerr << "# [error] arounding edges of edge <" << around_edge.first << "," << around_edge.second << "> is not correct:" << endl;
        for(size_t i = 0; i < one_ring_edges.size(); ++i){
          cerr << one_ring_edges[i].first << " --> " << one_ring_edges[i].second << endl;
        }
        return __LINE__;
      }
    }
#endif
  }
  return 0;
}

int adjust_cross_edges_right_hand_order_of_tet(const matrixst &tet,
                                               const matrixd &node,
                                               const std::pair<size_t,size_t> &around_edge,
                                               std::vector<size_t> &other_edge)
{
  other_edge.clear();
  for(size_t t = 0; t < 4; ++t)
    if(tet[t] != around_edge.first
       && tet[t] != around_edge.second){
      other_edge.push_back(tet[t]);
    }
  assert(other_edge.size() == 2);
  // to get right hand order if(dot(<D,C>,cross(<B,A>,<B,D>)) < 0) <D,C> reverse to <C,D>
  //       A
  //     / | \
  //    C--|--D
  //     \ | /
  //       B
  matrixd BA = node(colon(),around_edge.second) - node(colon(),around_edge.first);
  matrixd BD = node(colon(),other_edge[1]) - node(colon(),around_edge.first);
  matrixd DC = node(colon(),other_edge[0]) - node(colon(),other_edge[1]);
  if(dot(DC,cross(BA,BD)) < 0) swap(other_edge[0],other_edge[1]);
  other_edge.resize(2);
  return 0;
}

int remove_edge_in_tet_raw(
    const std::pair<size_t,size_t> & edge,
    std::vector<size_t> & tet_vec,
    const std::vector<double> & node_vec,
    std::vector<matrixd > & zyz_vec)
{
  itr_matrix<size_t*> tet_mat(4,tet_vec.size()/4,&tet_vec[0]);
  for(size_t t = 0; t < tet_vec.size(); ++t) {
    if(tet_vec[t] == edge.first)
      tet_vec[t] = edge.second;
  }
  vector<bool> is_tet_remained(tet_mat.size(2),true);
  set<size_t> one_tet;

  map<set<size_t>,vector<size_t> > duplicated_tets;
  for(size_t t = 0; t < tet_mat.size(2); ++t){
    one_tet.clear();
    one_tet.insert(tet_mat(colon(),t).begin(),
                   tet_mat(colon(),t).end());
    if(one_tet.size() != 4) {
      is_tet_remained[t] = false;
      continue;
    }
    duplicated_tets[one_tet].push_back(t);
  }

  // remove duplicated tets
  for(map<set<size_t>,vector<size_t> >::const_iterator msvcit =
      duplicated_tets.begin(); msvcit != duplicated_tets.end(); ++msvcit){
    const vector<size_t> & tets_vec = msvcit->second;
    if(tets_vec.size() != 1){
      for(size_t i = 0; i < tets_vec.size(); ++i)
        is_tet_remained[tets_vec[i]] = false;
    }
  }

  vector<size_t> new_tet_vec;
  vector<matrixd > new_zyz_vec;
  new_tet_vec.reserve(tet_vec.size());
  new_zyz_vec.reserve(zyz_vec.size());

  for(size_t t = 0; t < is_tet_remained.size(); ++t){
    if(is_tet_remained[t]){
      new_tet_vec.insert(new_tet_vec.end(),
                         tet_mat(colon(),t).begin(),
                         tet_mat(colon(),t).end());
      new_zyz_vec.push_back(zyz_vec[t]);
    }
  }

  tet_vec = new_tet_vec;
  zyz_vec = new_zyz_vec;
  return 0;
}

int remove_edge_and_update_info(
    const pair<size_t,size_t> & edge,
    std::vector<size_t> & tet_vec,
    const std::vector<double> & node_vec,
    std::vector<matrixd > & zyz_vec,
    jtf::mesh::one_ring_tet_at_edge & ortae)
{
  itr_matrix<size_t*> ori_tet_mat(4,tet_vec.size()/4,&tet_vec[0]);
  unique_ptr<jtf::mesh::face2tet_adjacent> fa_orign(
       jtf::mesh::face2tet_adjacent::create(ori_tet_mat,"topology"));

  typedef jtf::mesh::one_ring_tet_at_edge::e2tet_type::const_iterator oecit;
  typedef jtf::mesh::one_ring_tet_at_edge::e2tet_type::iterator oeit;

  oecit it = ortae.e2t_.find(edge);
  if(it == ortae.e2t_.end())
    it = ortae.e2t_.find(make_pair(edge.second,edge.first));
  if(it == ortae.e2t_.end()){
    cerr << "# [error] can not find edge < " << edge.first << ","
         << edge.second << "> in ortae." << endl;
    return __LINE__;
  }

  const vector<size_t> & around_tets = it->second;
  set<size_t> deleted_tets;
  for(size_t t = 0; t < around_tets.size(); ++t){
    if(around_tets[t] != -1)
      deleted_tets.insert(around_tets[t]);
  }

#if 1 // check each around tet should contain this edge
  {
    for(set<size_t>::const_iterator scit = deleted_tets.begin();
        scit != deleted_tets.end(); ++scit){
      if(find(ori_tet_mat(colon(),*scit).begin(),
              ori_tet_mat(colon(),*scit).end(),
              edge.first) == ori_tet_mat(colon(),*scit).end()
         || find(ori_tet_mat(colon(),*scit).begin(),
                 ori_tet_mat(colon(),*scit).end(),
                 edge.second) == ori_tet_mat(colon(),*scit).end())
      {
        cerr << "# [error] invalid around tet " << *scit << " ";
        copy(ori_tet_mat(colon(),*scit).begin(),
             ori_tet_mat(colon(),*scit).end(),
             ostream_iterator<size_t>(cerr,","));
        cerr << endl;
      }
    }
  }
#endif


  vector<bool> tet_remain_flag(ori_tet_mat.size(2),true);
  for(set<size_t>::const_iterator scit = deleted_tets.begin();
      scit != deleted_tets.end(); ++scit){
    tet_remain_flag[*scit] = false;
  }
  {
#if 0
    // update one_ring_tet_at_tet

    // remove all degenerated edges:
    // edge and all <edge.first, P> P is one ring point around edge
    // and update arounding tets information of all <P, edge.second>

    // this map store these tets which should be degenerated mapping to its upper
    // tet, here upper means the opposite tet throungh a face which linked edge.first
    //    map<size_t,size_t> degenerated_tet_to_up_tet;
    //    {
    //      vector<size_t> face(3);
    //      for(set<size_t>::const_iterator scit = deleted_tets.begin();
    //          scit != deleted_tets.end(); ++scit){
    //        size_t point_idx = 0;
    //        for(size_t i = 0; i < ori_tet_mat.size(1); ++i){
    //          if(ori_tet_mat(i,*scit) != edge.second)
    //            face[point_idx++] = ori_tet_mat(i,*scit);
    //        }
    //        assert(point_idx == 3);
    //        const pair<size_t,size_t> & tet_pair = fa_orign->query(&face[0]);
    //        if(tet_pair.first == -1 && tet_pair.second == -1){
    //          cerr << "# [error] face error <" ;
    //          copy(face.begin(),face.end(),ostream_iterator<size_t>(cerr,","));
    //          cerr << ">" << endl;
    //          return __LINE__;
    //        }
    //        degenerated_tet_to_up_tet[*scit] =
    //            (tet_pair.first == *scit? tet_pair.second: tet_pair.first);
    //      }
    //    } // end build degenerated tet to up tet map

    // suppose the edge is <first,second>, need to delete this edge,
    // and all surrounding edges <P,first>, where P is surrounding points of edge
    // <first,second>
    vector<size_t> one_ring_vertex;
    find_one_ring_vertex_around_edge(ori_tet_mat,edge,ortae,one_ring_vertex);

    {// remove degenerated edge
      oeit it = ortae.e2t_.find(edge);
      if(it == ortae.e2t_.end()) it = ortae.e2t_.find(make_pair(edge.second,edge.first));
      if(it == ortae.e2t_.end()) {
        cerr << "# [error] can not find edge < " << edge.first << ","
             << edge.second << "> in ortae." << endl;
        return __LINE__;
      }
      ortae.e2t_.erase(it);
    }

    { // merge around tets of <P,edge.first> and <P,edge.second>, and then remove
      // edge <P,edge.first>
      for(size_t t = 0; t < one_ring_vertex.size(); ++t){
        // get tets around <P,edge.first>
        oeit it_first = ortae.e2t_.find(make_pair(one_ring_vertex[t],edge.first));
        if(it_first == ortae.e2t_.end()){
          it_first = ortae.e2t_.find(make_pair(edge.first,one_ring_vertex[t]));
          if(it_first == ortae.e2t_.end()){
            cerr << "# [error] can not find edge " << one_ring_vertex[t]
                 << "," << edge.first << " in one ring tet at edge." << endl;
            return __LINE__;
          }
          vector<size_t> & around_tets = it_first->second;
          std::reverse(around_tets.begin(),around_tets.end());
        }
        vector<size_t> & around_tets_first = it_first->second;

        if(around_tets_first.front() != around_tets_first.back()) {
          if(around_tets_first.front() != -1 && around_tets_first.back() != -1){
            cerr << "# [error] invalid tets around edge " << one_ring_vertex[t]
                 << "," << edge.first;
            copy(around_tets_first.begin(),around_tets_first.end(),
                 ostream_iterator<size_t>(cerr,","));
            cerr << endl;
            return __LINE__;
          }
          if(around_tets_first.front() != -1)
            around_tets_first.insert(around_tets_first.begin(),-1);
          if(around_tets_first.back() != -1)
            around_tets_first.push_back(-1);
        }

        assert(around_tets_first.front() == around_tets_first.back());
        around_tets_first.pop_back(); // no need the duplicated tet idx

        // get tets around <P,edge.second>
        oeit it_second = ortae.e2t_.find(make_pair(one_ring_vertex[t],edge.second));
        if(it_second == ortae.e2t_.end()){
          it_second = ortae.e2t_.find(make_pair(edge.second,one_ring_vertex[t]));
          if(it_second == ortae.e2t_.end()){
            cerr << "# [error] can not find edge " << one_ring_vertex[t]
                 << "," << edge.second << " in one ring tet at edge." << endl;
            return __LINE__;
          }
          vector<size_t> & around_tets = it_second->second;
          std::reverse(around_tets.begin(),around_tets.end());
        }
        vector<size_t> & around_tets_second = it_second->second;

        if(around_tets_second.front() != around_tets_second.back()) {
          if(around_tets_second.front() != -1 && around_tets_second.back() != -1){
            cerr << "# [error] invalid tets around edge " << one_ring_vertex[t]
                 << "," << edge.first ;
            copy(around_tets_second.begin(),around_tets_second.end(),
                 ostream_iterator<size_t>(cerr,","));
            cerr << endl;
            return __LINE__;
          }
          if(around_tets_second.front() != -1)
            around_tets_second.insert(around_tets_second.begin(),-1);
          if(around_tets_second.back() != -1)
            around_tets_second.push_back(-1);
        }

        assert(around_tets_second.front() == around_tets_second.back());
        around_tets_second.pop_back();

        // find the tet pair of face (P, edge.first ,edge.second)
        const size_t face_idx  = fa_orign->get_face_idx(one_ring_vertex[t],
                                                        edge.first,edge.second);
        if(face_idx == -1){
          cerr << "# [error] invalid face " << one_ring_vertex[t]
               << "," << edge.first << "," << edge.second << endl;
          return __LINE__;
        }
        const pair<size_t,size_t> & tet_pair = fa_orign->face2tet_[face_idx];

        // merge around_tets_first and around_tets_second
        {
          deque<size_t> de_around_tets_first(around_tets_first.size());
          deque<size_t> de_around_tets_second(around_tets_second.size());
          copy(around_tets_first.begin(),around_tets_first.end(),
               de_around_tets_first.begin());
          copy(around_tets_second.begin(),around_tets_second.end(),
               de_around_tets_second.begin());
          while((de_around_tets_first.back() != tet_pair.first &&
                 de_around_tets_first.back() != tet_pair.second) ||
                (de_around_tets_first[de_around_tets_first.size()-2] != tet_pair.first &&
                 de_around_tets_first[de_around_tets_first.size()-2] != tet_pair.second))
          {
            de_around_tets_first.push_front(de_around_tets_first.back());
            de_around_tets_first.pop_back();
          }
          while((de_around_tets_second.front() != tet_pair.first &&
                 de_around_tets_second.front() != tet_pair.second) ||
                (de_around_tets_second[1] != tet_pair.first &&
                 de_around_tets_second[1] != tet_pair.second))
          {
            de_around_tets_second.push_front(de_around_tets_second.back());
            de_around_tets_second.pop_back();
          }
          assert(de_around_tets_first.back() == de_around_tets_second.front());
          assert(de_around_tets_first[de_around_tets_first.size()-2]
                 == de_around_tets_second[1]);

          // first may be A B C t0 t1 (A), second may be t1 t0 E F G (t1)
          // t0 and t1 are tets should be ne deleted
          // the merged loop should be A B C E F G A
          de_around_tets_second.pop_front();
          de_around_tets_second.pop_front();

          de_around_tets_first.pop_back();
          de_around_tets_first.pop_back();
          de_around_tets_second.insert(de_around_tets_second.begin(),
                                       de_around_tets_first.begin(),
                                       de_around_tets_first.end());
          de_around_tets_second.push_back(de_around_tets_second.front());

          around_tets_second.resize(de_around_tets_second.size());
          copy(de_around_tets_second.begin(),de_around_tets_second.end(),
               around_tets_second.begin());
          ortae.e2t_.erase(it_first);
        }
      }
    }
#endif
#if 0
    {// remove all <one ring vertex,edge.first> edges
      for(size_t t = 0; t < one_ring_vertex.size(); ++t){
        oeit it = ortae.e2t_.find(make_pair(one_ring_vertex[t],edge.first));
        if(it == ortae.e2t_.end())
          it = ortae.e2t_.find(make_pair(edge.first,one_ring_vertex[t]));
        if(it == ortae.e2t_.end()) {
          cerr << "# [error] can not find edge < " << edge.first << ","
               << one_ring_vertex[t] << "> in ortae." << endl;
          return __LINE__;
        }
        ortae.e2t_.erase(it);
      }
    }
    {// modify the other edges <one_ring_vertex, edge.second>, change the
      for(size_t t = 0; t  < one_ring_vertex.size(); ++t){
        oeit it = ortae.e2t_.find(make_pair(one_ring_vertex[t],edge.second));
        if(it == ortae.e2t_.end())
          it = ortae.e2t_.find(make_pair(edge.second,one_ring_vertex[t]));
        if(it == ortae.e2t_.end()) {
          cerr << "# [error] can not find edge < " << edge.second << ","
               << one_ring_vertex[t] << "> in ortae." << endl;
          return __LINE__;
        }
        vector<size_t> & around_tets = it->second;

        for(size_t i = 0; i < around_tets.size(); ++i){
          map<size_t,size_t>::const_iterator msscit =
              degenerated_tet_to_up_tet.find(around_tets[i]);
          if(msscit == degenerated_tet_to_up_tet.end()) continue;
          around_tets[i] = msscit->second;
        }
        // check whether there is -1 in  around_tets
        size_t outside_tet_num = count(around_tets.begin(),around_tets.end(),-1);
        if(outside_tet_num == 0)  continue;

        if(outside_tet_num > 0) {
          // if outside tets are more than 2, need to remove the extra ones,
          // if there are some outside tets beside each other, remove them and only
          // leave one
          vector<size_t> new_around_tets;
          for(size_t t = 0; t < around_tets.size(); ++t){
            if(new_around_tets.empty()) new_around_tets.push_back(around_tets[t]);
            if(new_around_tets.back() == -1 && around_tets[t] == -1) continue;
            new_around_tets.push_back(around_tets[t]);
          }

          if(new_around_tets.front() == -1 && new_around_tets.back() != -1)
            new_around_tets.push_back(-1);
          if(new_around_tets.front() != -1 && new_around_tets.back() == -1)
            new_around_tets.insert(new_around_tets.begin(),-1);

          assert(new_around_tets.front() == new_around_tets.back());

          around_tets = new_around_tets;
          outside_tet_num = count(new_around_tets.begin(),new_around_tets.end(),-1);
          assert(outside_tet_num <= 2);
        }

        if(outside_tet_num == 1) { // if there are one outside tet inside the sequence, need to rotate it
          assert(around_tets.front() == around_tets.back());
          assert(around_tets.back() != -1);
          around_tets.pop_back();
          vector<size_t>::iterator out_tet_ptr =
              find(around_tets.begin(), around_tets.end(), -1);
          assert(out_tet_ptr != around_tets.end());

          std::rotate(around_tets.begin(),out_tet_ptr,around_tets.end());
          around_tets.push_back(-1);
          continue;
        }

        if(outside_tet_num == 2) { // if there are two outside tets, they must be on both ends
          if(around_tets.front() != -1 && around_tets.back() != -1){
            cerr << "# [error] strange: the around tet loop contains two outside tets,"
                 << " but inside the sequence loop" << endl;
            return __LINE__;
          }
          continue;
        }
      }
    }
#endif
#if 0
    for(one_ring_tet_at_edge::e2tet_type::iterator oeit = ortae.e2t_.begin();
        oeit != ortae.e2t_.end();){
      const pair<size_t,size_t> & edge_ = oeit->first;
      vector<size_t> & around_tets = oeit->second;
      if((edge_.first == edge.first && edge_.second == edge.second) ||
         (edge_.first == edge.second && edge_.second == edge.first)){
        ortae.e2t_.erase(oeit++);
        continue;
      }

      // check each tet in around_tets, if there are degenerated tets, need to
      // replace this with the opposite one
      for(size_t i = 0; i < around_tets.size(); ++i){
        map<size_t,size_t>::const_iterator msscit =
            degenerated_tet_to_up_tet.find(around_tets[i]);
        if(msscit == degenerated_tet_to_up_tet.end()) continue;
        around_tets[i] = msscit->second;
      }
      // check whether there is -1 in  around_tets
      size_t outside_tet_num = count(around_tets.begin(),around_tets.end(),-1);
      if(outside_tet_num == 0) {
        ++oeit;
        continue;
      }

      if(outside_tet_num > 0) {
        // if outside tets are more than 2, need to remove the extra ones,
        // if there are some outside tets beside each other, remove them and only
        // leave one
        vector<size_t> new_around_tets;
        for(size_t t = 0; t < around_tets.size(); ++t){
          if(new_around_tets.empty()) new_around_tets.push_back(around_tets[t]);
          if(new_around_tets.back() == -1 && around_tets[t] == -1) continue;
          new_around_tets.push_back(around_tets[t]);
        }

        if(new_around_tets.front() == -1 && new_around_tets.back() != -1)
          new_around_tets.push_back(-1);
        if(new_around_tets.front() != -1 && new_around_tets.back() == -1)
          new_around_tets.insert(new_around_tets.begin(),-1);

        assert(new_around_tets.front() == new_around_tets.back());

        around_tets = new_around_tets;
        outside_tet_num = count(new_around_tets.begin(),new_around_tets.end(),-1);
        assert(outside_tet_num <= 2);
      }

      if(outside_tet_num == 1) { // if there are one outside tet inside the sequence, need to rotate it
        assert(around_tets.front() == around_tets.back());
        around_tets.pop_back();
        vector<size_t>::iterator out_tet_ptr =
            find(around_tets.begin(), around_tets.end(), -1);
        assert(out_tet_ptr != around_tets.end());

        std::rotate(around_tets.begin(),out_tet_ptr,around_tets.end());
        around_tets.push_back(-1);
        ++oeit;
        continue;
      }

      if(outside_tet_num == 2) { // if there are two outside tets, they must be on both ends
        if(around_tets.front() != -1 &&
           around_tets.back() != -1){
          cerr << "# [error] strange: the around tet loop contains two outside tets,"
               << " but inside the sequence loop" << endl;
          return __LINE__;
        }
        ++oeit;
      }
    }
#endif
  }


  {// delete all arounding tets and modify the tets which point to edge.first
    // to edge.second
    vector<size_t> new_tet_vec(tet_vec.size() - 4 * deleted_tets.size());
    vector<matrixd > new_zyz_vec(tet_vec.size()/4 - deleted_tets.size());

    itr_matrix<size_t*> new_tet_mat(4,new_tet_vec.size()/4,&new_tet_vec[0]);
    itr_matrix<size_t*> tet_mat(4,tet_vec.size()/4,&tet_vec[0]);
    size_t new_tet_idx = 0;
    for(size_t t = 0; t < tet_mat.size(2); ++t){
      if(tet_remain_flag[t]){
        new_tet_mat(colon(),new_tet_idx) = tet_mat(colon(),t);
        new_zyz_vec[new_tet_idx] = zyz_vec[t];
        ++new_tet_idx;
      }
    }
    assert(new_tet_idx == new_tet_mat.size(2));

    for(size_t t = 0; t < new_tet_vec.size(); ++t){
      if(new_tet_vec[t] == edge.first) new_tet_vec[t] = edge.second;
    }
    tet_vec = new_tet_vec;
    zyz_vec = new_zyz_vec;

#if 1
    assert(find(tet_vec.begin(),tet_vec.end(),edge.first) == tet_vec.end());
#endif

    //orient_tet_raw(&node_vec[0],node_vec.size()/3, &tet_vec[0],tet_vec.size()/4);
  }

  {// update the one ring tet info
    ortae.e2t_.clear();
    itr_matrix<size_t*> tet_mat(4,tet_vec.size()/4,&tet_vec[0]);
    itr_matrix<const double*> node_mat(3,node_vec.size()/3,&node_vec[0]);

    //      size_t negative_tet_idx = find_first_negative_tet(node_mat,tet_mat);
    //      while(negative_tet_idx != -1){
    //        cerr << "# [info] find negative tet idx " << negative_tet_idx << endl;
    //        swap(tet_mat(1,negative_tet_idx), tet_mat(2,negative_tet_idx));
    //      }
    //orient_tet_exp(node_mat,tet_mat);
#if 0 // debug
    {
      ofstream ofs("debug_tet.vtk");
      tet2vtk(ofs,&node_mat[0],node_mat.size(2),&tet_mat[0],tet_mat.size(2));
    }
#endif
    unique_ptr<jtf::mesh::face2tet_adjacent> fa(
         jtf::mesh::face2tet_adjacent::create(tet_mat,"topology"));
    ortae.add_tets(tet_mat,*fa);
    ortae.sort_into_loop(tet_mat,node_mat);
  }
  return 0;
}

int remove_one_chain_and_update_info(
    const std::deque<std::pair<size_t,size_t> > & one_chain,
    std::vector<size_t> & tet_vec,
    const std::vector<double> & node_vec,
    std::vector<matrixd > & zyz_vec,
    jtf::mesh::one_ring_tet_at_edge & ortae)
{

  cerr << "# [info] begin to remove one chain: edge num " << one_chain.size() << endl;
  for(size_t t = 0; t < one_chain.size(); ++t){
    const pair<size_t,size_t > & edge = one_chain[t];
    remove_edge_in_tet_raw(edge,tet_vec,node_vec,zyz_vec);
    cerr << "# [info] ------ edge " << t << " finished." << endl;
  }
  return 0;
}

int get_face_patches_according_to_type(
    const zjucad::matrix::matrix<size_t> & faces,
    const zjucad::matrix::matrix<size_t> & face_type,
    const jtf::mesh::one_ring_face_at_point &orfap,
    std::vector<std::vector<size_t> > &patches)
{
  vector<bool> is_face_visited(faces.size(2), false);

  stack<size_t> face_stack;

  patches.clear();

  vector<bool>::const_iterator cit = find(is_face_visited.begin(),
                                          is_face_visited.end(), false);
  while(cit != is_face_visited.end()){
    set<size_t> one_patch;
    face_stack.push(cit-is_face_visited.begin());
    while(!face_stack.empty()){
      const size_t f_idx = face_stack.top();
      face_stack.pop();
      if(is_face_visited[f_idx]) continue;
      is_face_visited[f_idx] = true;
      one_patch.insert(f_idx);

      for(size_t pi = 0; pi < faces.size(1); ++pi){
          const auto it = orfap.p2f_.find(faces(pi,f_idx));
          if(it == orfap.p2f_.end()){
              throw std::logic_error("can not find point in one_ring_face_at_point");
            }
          const vector<size_t> & one_ring_face = it->second;
          for(size_t fi = 0; fi < one_ring_face.size(); ++fi){
              if(one_ring_face[fi] == -1) continue;
              if(is_face_visited[one_ring_face[fi]]) continue;
              if(face_type[one_ring_face[fi]] == face_type[f_idx])
                face_stack.push(one_ring_face[fi]);
            }
        }
    }
    vector<size_t> one_patch_vec(one_patch.size());
    std::copy(one_patch.begin(), one_patch.end(), one_patch_vec.begin());
    patches.push_back(one_patch_vec);
    cit = find(is_face_visited.begin(), is_face_visited.end(), false);
  }

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
    set<size_t> one_patch;
    face_stack.push(cit-is_face_visited.begin());
    while(!face_stack.empty()){
      const size_t f_idx = face_stack.top();
      face_stack.pop();
      is_face_visited[f_idx] = true;
      one_patch.insert(f_idx);

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
    vector<size_t> one_patch_vec(one_patch.size());
    std::copy(one_patch.begin(), one_patch.end(), one_patch_vec.begin());
    patches.push_back(one_patch_vec);
    cit = find(is_face_visited.begin(), is_face_visited.end(), false);
  }

  return 0;
}

