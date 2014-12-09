#include "extract_loop_from_undirected_edges.h"

#include <stack>

#include <boost/unordered_map.hpp>

#include <jtflib/util/util.h>
#include <jtflib/util/vertex_connection.h>
#include "../common/util.h"

using namespace std;

int extract_chains_from_undirected_edges(
    const boost::unordered_set<pair<size_t,size_t> > &edges,
    vector<deque<pair<size_t,size_t> > >  &chains,
    vector<deque<pair<size_t,size_t> > >  &loops,
    const size_t mode)
{
  typedef boost::unordered_set<size_t> linked_point_type;
  typedef boost::unordered_map<size_t, linked_point_type > point_adj_type;
  typedef linked_point_type::const_iterator lpt_cit;
  typedef linked_point_type::iterator lpt_it;
  typedef point_adj_type::const_iterator pat_cit;
  typedef point_adj_type::iterator pat_it;

  point_adj_type point_adj;
  for(boost::unordered_set<pair<size_t,size_t> >::const_iterator buscit
      = edges.begin(); buscit != edges.end(); ++buscit){
    const pair<size_t,size_t> & one_edge = *buscit;
    point_adj[one_edge.first].insert(one_edge.second);
    point_adj[one_edge.second].insert(one_edge.first);
  }

  vector<std::shared_ptr<vertex_connection<UNDIRECT> > > graph_trees;
  boost::unordered_set<size_t> passed_points;

  boost::unordered_set<pair<size_t,size_t> > edges_to_extract_chains;

  // remove degree 1 points
  while(1){
    bool find_degree_one = false;
    for(pat_it pit = point_adj.begin(); pit != point_adj.end();){
      if(pit->second.size() == 1){
        find_degree_one = true;
        const size_t linked_point = *(pit->second.begin());

        pair<size_t,size_t> edge_left(pit->first, linked_point);
        if(edge_left.first > edge_left.second)
          swap(edge_left.first, edge_left.second);
        edges_to_extract_chains.insert(edge_left);

        pat_it p_it = point_adj.find(linked_point);
        if(p_it == point_adj.end()){
          point_adj.erase(pit++);
        }else{
          linked_point_type & lpt = p_it->second;
          lpt_it lit = lpt.find(pit->first);
          if(lit != lpt.end()){
            lpt.erase(lit);
          }
          point_adj.erase(pit++);
        }
      }else{
        ++pit;
      }
    }
    if(!find_degree_one) break;
  }
  for(pat_it pit = point_adj.begin(); pit != point_adj.end();){
    if(pit->second.empty())
      point_adj.erase(pit++);
    else
      ++pit;
  }
  //  {
  //    for(pat_it pit = point_adj.begin(); pit != point_adj.end(); ++pit){
  //      if(pit->second.size() < 2)
  //        cerr << "# [error] point " << pit->first << " size = "
  //             << pit->second.size() << endl;
  //    }

  //  }

  if(mode == 1){ // fast return if there are loops
    if(!point_adj.empty()){
      return __LINE__;
    }
  }

  {
    vector<pair<size_t,size_t> > edges_vec(edges_to_extract_chains.size());
    copy(edges_to_extract_chains.begin(), edges_to_extract_chains.end(),
         edges_vec.begin());
    jtf::util::extract_chain_from_edges(edges_vec, chains);
  }
  while(!point_adj.empty()){
    stack<size_t> tree;
    const size_t passed_points_num = passed_points.size();
    for(pat_cit cit = point_adj.begin(); cit != point_adj.end(); ++cit){
      passed_points.insert(cit->first);
      if(passed_points.size() != passed_points_num){
        tree.push(cit->first);
        break;
      }
    }
    if(tree.empty()) break;
    map<pair<size_t,size_t>,double> graph_edges;
    while(!tree.empty()){
      pat_it it = point_adj.find(tree.top());
      if(it == point_adj.end()){
        tree.pop();
        continue;
      }
      linked_point_type & lp = it->second;
      const size_t tree_size = tree.size();
      for(lpt_it lit = lp.begin(); lit != lp.end(); ){
        const size_t passed_points_num = passed_points.size();
        passed_points.insert(*lit);
        if(passed_points.size() != passed_points_num){ // have not passed yet
          graph_edges.insert(make_pair(make_pair(it->first,*lit),1));
          tree.push(*lit);
          pat_it it_adj = point_adj.find(*lit);
          if(it_adj != point_adj.end()){
            linked_point_type & lp_adj = it_adj->second;
            lpt_it lp_adj_it = lp_adj.find(it->first);
            assert(lp_adj_it != lp_adj.end());
            lp_adj.erase(lp_adj_it);
          }
          lp.erase(lit++);
        }else
          ++lit;
      }
      if(tree.size() == tree_size){ // can not add more edges
        tree.pop();
      }
    }
    std::shared_ptr<vertex_connection<UNDIRECT> > graph_(
          vertex_connection<UNDIRECT>::create(graph_edges));
    graph_trees.push_back(graph_);
  }


  // gather the left edges
  boost::unordered_set<pair<size_t,size_t> > left_edges;
  for(pat_cit pcit = point_adj.begin(); pcit != point_adj.end(); ++pcit){
    const linked_point_type & lp = pcit->second;
    for(lpt_cit lcit = lp.begin(); lcit != lp.end(); ++lcit){
      pair<size_t,size_t> left_edge(pcit->first, *lcit);
      if(left_edge.first > left_edge.second)
        swap(left_edge.first, left_edge.second);
      left_edges.insert(left_edge);
    }
  }

  loops.clear();
  // check each left edge for every sub-graph
  for(boost::unordered_set<pair<size_t,size_t> >::const_iterator pcit
      = left_edges.begin(); pcit != left_edges.end(); ++pcit){
    const pair<size_t,size_t> & one_edge = *pcit;

    for(size_t gi = 0; gi < graph_trees.size(); ++gi){
      vector<size_t> path;
      std::shared_ptr<vertex_connection<UNDIRECT> > pt = graph_trees[gi];
      pt->get_shortest_path(one_edge.first, one_edge.second, path);
      if(!path.empty()){
        path.push_back(one_edge.first);
        deque<pair<size_t,size_t> > loop;
        for(size_t pi = 0; pi < path.size()-1; ++pi){
          loop.push_back(make_pair(path[pi],path[pi+1]));
        }
        assert(loop.front().first == loop.back().second);
        loops.push_back(loop);
        break;
      }
    }
  }


  return 0;
}
