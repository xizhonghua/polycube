#ifndef EXTRACT_LOOP_FROM_UNDIRECTED_EDGES_H
#define EXTRACT_LOOP_FROM_UNDIRECTED_EDGES_H

#include <boost/unordered_set.hpp>
#include <vector>
#include <deque>

/**
 * @brief This function is used to extract loop from given edges
 *        First, we build minimal spanning trees
 *        Then, check left edges, insert each edge will make a loop
 *
 * @param edges undirected edges
 * @param loops output loops
 * @param mode  mode = 0: extract all loops
 *              mode = 1: fast return when meet loops
 * @return int  return 0 if works fine or return non-zeros
 */
int extract_chains_from_undirected_edges(
    const boost::unordered_set<std::pair<size_t,size_t> > & edges,
    std::vector<std::deque<std::pair<size_t,size_t> > > & chains,
    std::vector<std::deque<std::pair<size_t,size_t> > > & loops,
    const size_t mode = 0);

#endif // EXTRACT_LOOP_FROM_UNDIRECTED_EDGES_H
