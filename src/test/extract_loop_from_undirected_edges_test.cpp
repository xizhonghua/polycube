#include "extrac_loop_from_undirected_edges_test.h"
#include <map>

using namespace std;

void extract_loop_from_undirected_edges_test::extract_loops()
{
  boost::unordered_set<pair<size_t,size_t> > edges;
  edges.insert(make_pair(3676,4201));
  edges.insert(make_pair(3676,3712));
  edges.insert(make_pair(3676,11488));
  edges.insert(make_pair(4201,4204));
  edges.insert(make_pair(12064,7285));
  edges.insert(make_pair(12064,6703));
  edges.insert(make_pair(4204,11488));
  edges.insert(make_pair(6703,7492));
  edges.insert(make_pair(6703,5188));
  edges.insert(make_pair(6703,12064));
  edges.insert(make_pair(7492,8218));
  edges.insert(make_pair(6061,6331));
  edges.insert(make_pair(6061,12892));
  edges.insert(make_pair(6061,6742));
  edges.insert(make_pair(6331,6334));
  edges.insert(make_pair(3712,604));
  edges.insert(make_pair(6334,12892));
  edges.insert(make_pair(6334,7285));
  edges.insert(make_pair(6334,6742));
  edges.insert(make_pair(5188,8218));
  edges.insert(make_pair(604,11488));

  vector<deque<pair<size_t,size_t> > > loops;
  vector<deque<pair<size_t,size_t> > > chains;
  extract_chains_from_undirected_edges(edges,chains, loops);

  for(size_t li = 0; li < loops.size(); ++li){
    cerr << "# loop " << li << endl;
    const deque<pair<size_t,size_t> > & one_loop = loops[li];
    for(size_t ei = 0; ei < one_loop.size(); ++ei){
      cerr << one_loop[ei].first << " --> " << one_loop[ei].second << endl;
    }
  }

}
