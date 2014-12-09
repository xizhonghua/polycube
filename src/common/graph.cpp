#include "graph.h"

using namespace std;

void convert(const dyn_graph &dg, fix_graph &fg)
{
	const size_t nn = dg.size();
	fg.ptr_.resize(nn+1);
	fg.ptr_[0] = 0;
	for(size_t ni = 0; ni < nn; ++ni) {
		fg.ptr_[ni+1] = fg.ptr_[ni] + dg[ni].size();
	}
	fg.idx_.resize(fg.ptr_[nn]);
	for(size_t ni = 0; ni < nn; ++ni) {
		set<node_idx>::const_iterator iter = dg[ni].begin();
		for(size_t ei = fg.ptr_[ni]; ei < fg.ptr_[ni+1]; ++ei, ++iter) {
			fg.idx_[ei] = *iter;
		}
	}
}

void find_undirect_edges(
	const fix_graph &fg,
	const std::vector<std::pair<node_idx, node_idx> > &edges,
	std::vector<edge_idx> &idx)
{
	idx.reserve(edges.size());
	idx.clear();

	size_t i;
	for(size_t ei = 0; ei < edges.size(); ++ei) {
		i = fg.find_edge(edges[ei].first,
						 edges[ei].second);
		assert(fg.is_valid_idx(i));
		idx.push_back(i);

		i = fg.find_edge(edges[ei].second,
						 edges[ei].first);
		assert(fg.is_valid_idx(i));
		idx.push_back(i);
	}
}

