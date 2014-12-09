#include "topology.h"

using namespace zjucad::matrix;

#include <vector>
#include <set>
using namespace std;

namespace hj { namespace mesh {

/// cell2graph
static void cell2graph(const matrix<int> &cell, vector<set<int> > &graph, bool self)
{
	const int vn = max(cell)+1, cn = cell.size(2), vpc = cell.size(1);
  int i, j, k;
	graph.resize(vn);
  for(i = 0; i < cn; ++i)	// for each cell
    for(j = 0; j < vpc; ++j)	// for each cell vertex
      for(k = 0; k < vpc; ++k)	// add others
        if(j != k) graph[cell(j, i)].insert(cell(k, i));
	if(self) {
		for(i = 0; i < vn; ++i)
			graph[i].insert(i);
	}
}

static void convert(vector<set<int> > &vsgraph, matrix<matrix<int> > &mmgraph)
{
    const int vn = static_cast<int>(vsgraph.size());
    vsgraph.resize(vn);
    for(int i = 0; i < vn; ++i) {
        mmgraph[i].resize(static_cast<int>(vsgraph[i].size()));
        copy(vsgraph[i].begin(), vsgraph[i].end(), mmgraph[i].begin());
		vsgraph[i].clear();	// release the memory, different to vector
    }
}

void cell2graph(const matrix<int> &cell, matrix<matrix<int> > &graph, bool self)
{
    vector<set<int> > vsgraph;	// vsgraph edge
	cell2graph(cell, vsgraph, self);
	graph.resize(vsgraph.size());
	convert(vsgraph, graph);
}

static void convert(vector<set<int> > &vsgraph, matrix<int> &ptr, matrix<int> &idx)
{
    const int vn = static_cast<int>(vsgraph.size());
	int i;
    ptr.resize(vn+1);
	ptr[0] = 0;
    for(i = 0; i < vn; ++i)
		ptr[i+1] = ptr[i]+static_cast<int>(vsgraph[i].size());

	idx.resize(ptr[vn]);
	for(i = 0; i < vn; ++i) {
        copy(vsgraph[i].begin(), vsgraph[i].end(), &idx[ptr[i]]);
		vsgraph[i].clear();
    }
}

void cell2csc(const matrix<int> &cell, matrix<int> &ptr, matrix<int> &idx, bool self)
{
  vector<set<int> > vsgraph;
	cell2graph(cell, vsgraph, self);
	convert(vsgraph, ptr, idx);
}

}}
