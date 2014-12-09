#ifndef _HJ_HEX_GRAPH_H_
#define _HJ_HEX_GRAPH_H_

#include "config.h"

#include <vector>
#include <set>
#include <cassert>
#include <algorithm>

typedef size_t node_idx;
typedef size_t edge_idx;
typedef std::vector<std::set<node_idx> > dyn_graph;

/**
   \brief directional graph, store both i->j and j->i
 */
class HEXGEN_COMMON_API fix_graph {
public:
	inline size_t node_num (void) const {return ptr_.size()-1; }
	inline size_t edge_num (void) const {return idx_.size(); }

	typedef std::vector<node_idx>::const_iterator edge_const_iterator;
	typedef std::vector<node_idx>::iterator edge_iterator;

	inline edge_const_iterator begin(node_idx i) const {
		return idx_.begin()+ptr_[i];
	}
	inline edge_const_iterator end(node_idx i) const {
		return idx_.begin()+ptr_[i+1];
	}
	inline edge_iterator begin(node_idx i) {
		return idx_.begin()+ptr_[i];
	}
	inline edge_iterator end(node_idx i) {
		return idx_.begin()+ptr_[i+1];
	}
	inline edge_idx find_edge(node_idx i, node_idx j) const {
		assert(i >= 0 && i < node_num()
			   && j >= 0 && j < node_num());
		for(size_t ei = ptr_[i]; ei < ptr_[i+1]; ++ei) {
			if(idx_[ei] == j)
				return ei;
		}
		return -1;
	}
	inline std::pair<node_idx, node_idx> get_edge(edge_idx ei) const {
		assert(is_valid_idx(ei) && ei < edge_num());
		const node_idx beg = std::upper_bound(ptr_.begin(), ptr_.end(), ei)-ptr_.begin()-1;
		assert(find_edge(beg, idx_[ei]) == ei);
		return std::make_pair(beg, idx_[ei]);
	}
	template <typename T> static inline bool is_valid_idx(T idx) {
		return idx != static_cast<T>(-1);
	}
	std::vector<size_t> ptr_;
	std::vector<node_idx> idx_;
};

HEXGEN_COMMON_API void convert(const dyn_graph &dg, fix_graph &fg);

/**
   find the index of ij and ji in edge
 */
HEXGEN_COMMON_API void find_undirect_edges(
	const fix_graph &fg,
	const std::vector<std::pair<node_idx, node_idx> > &edges,
	std::vector<edge_idx> &idx);

#endif
