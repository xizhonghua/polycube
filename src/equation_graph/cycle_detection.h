#ifndef CYCLE_DETECTION_H
#define CYCLE_DETECTION_H

#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/kruskal_min_spanning_tree.hpp>
#include <boost/graph/dijkstra_shortest_paths.hpp>
#include <boost/graph/graph_traits.hpp>
#include <boost/static_assert.hpp>
#include <iostream>

#include "util.h"

enum DIRECTION2{UNDIRECT2,DIRECT2};


class cycle_detection
{
public:

  //! @brief check whether the graph has a cycle
  //! @param edges input graph
  //! @return return true if has cycles or false
  template <DIRECTION2 DIRECTION2_TYPE>
  static bool has_cycle(const std::vector<std::pair<size_t,size_t> > &edges)
  {
    std::cerr << "# [error] should not be here." << std::endl;
    return false;
  }

};

//! @brief check cycles in an UNDIRECT2 graph: remove all degree one points
//!       iteratively, if all points are removed, no cycle or not
//! @param edges input edges which consist a graph
//! @return return true if find a graph or false
template<>
bool cycle_detection::has_cycle<UNDIRECT2>(const std::vector<std::pair<size_t,size_t> > &edges)
{
  using namespace std;
  map<size_t,set<size_t> > point_degree;
  for(size_t ei = 0; ei < edges.size(); ++ei){
      point_degree[edges[ei].first].insert(edges[ei].second);
      point_degree[edges[ei].second].insert(edges[ei].first);
    }

  while(1){
      bool find_degree_one_point = false;
      for(map<size_t, set<size_t> >::iterator it = point_degree.begin();
          it != point_degree.end(); ++it){
          // a--> b
          // b--> (a,c,...)
          if(it->second.size() == 1) { // need to remove degree 1 points
              find_degree_one_point = true;
              map<size_t,set<size_t> >::iterator itb =
                  point_degree.find(*(it->second.begin()));
              assert_iterator(point_degree, itb,
                              "# [error] strange can not find point in map");
              assert(itb->second.find(it->first) != itb->second.end());
              set<size_t>::iterator itbs = itb->second.find(it->first);
              assert_iterator(itb->second, itbs, "# [error] strange can not find point in set");
              itb->second.erase(itbs);
              if(itb->second.empty())
                point_degree.erase(itb);
              point_degree.erase(it);
              break;
            }
        }
      if(!find_degree_one_point)
        break;
    }
  if(point_degree.empty())
    return false;
  else
    return true;
}

//! @brief check cycles in an DIRECT graph: remove all outgoing degree one points
//!       iteratively, if all points are removed, no cycle or not
//! @param edges input edges which consist a graph
//! @return return true if find a graph or false
template<>
bool cycle_detection::has_cycle<DIRECT2>(const std::vector<std::pair<size_t,size_t> > &edges)
{
  using namespace std;
  typedef map<size_t, set<pair<size_t,int> > > point_degree_type;
  point_degree_type point_degree;
  // 1: outgoing point; -1 incoming point
  for(size_t ei = 0; ei < edges.size(); ++ei){
      point_degree[edges[ei].first].insert(make_pair(edges[ei].second, 1));
      point_degree[edges[ei].second].insert(make_pair(edges[ei].first, -1));
    }

  while(1){
      bool find_outgoing_degree_one_point = false;

      for(point_degree_type::iterator it = point_degree.begin();
          it != point_degree.end(); ++it){
          const set<pair<size_t,int> > & connected_points = it->second;

          bool is_outgoing_point = true;
          for(set<pair<size_t,int> >::const_iterator spsit =
              connected_points.begin();
              spsit != connected_points.end(); ++spsit){
              if(spsit->second == -1) {
                  is_outgoing_point = false;
                  break;
                }
            }
          if(is_outgoing_point){
              find_outgoing_degree_one_point = true;
              for(set<pair<size_t,int> >::const_iterator
                  spsit = connected_points.begin();
                  spsit != connected_points.end(); ++spsit){
                  point_degree_type::iterator itb = point_degree.find(spsit->first);
                  assert_iterator(point_degree, itb,
                                  "# [error] strange can not find point");
                  assert(itb->second.find(make_pair(it->first, -1)) != itb->second.end());
                  set<pair<size_t,int> >::iterator spit =
                      itb->second.find(make_pair(it->first,-1));
                  assert_iterator(itb->second, spit,
                                  "# [error] can not find point in set.");
                  itb->second.erase(spit);
                  if(itb->second.empty())
                    point_degree.erase(itb);
                }
              point_degree.erase(it);
              break;
            }
        }

      if(!find_outgoing_degree_one_point)
        break;
    }
  if(point_degree.empty())
    return false;
  else
    return true;
}

#endif // CYCLE_DETECTION
