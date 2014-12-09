#ifndef HJ_DIJKSTRA_H_
#define HJ_DIJKSTRA_H_

#include <limits>
#include <vector>

#include <hjlib/algorithm/heap.h>

template <typename T>
class length_one
{
public:
	template <typename Idx>
	T operator()(Idx a, Idx b) const { // assume (a, b) is an edge
		return T(1);
	}
};


template <typename DstType,
		  typename RandomIterator1, typename RandomIterator2, typename Iterator,
		  typename RandomIterator3, typename RandomIterator4,
		  typename EdgeLen>
void dijkstra(RandomIterator1 ptr_first, RandomIterator1 ptr_last,
			  RandomIterator2 idx_first, RandomIterator2 idx_last,
			  Iterator zero_first, Iterator zero_last,
			  RandomIterator3 dst_first, RandomIterator3 dst_last,
			  RandomIterator4 prv_first, RandomIterator4 prv_last,
			  const EdgeLen &edge_len = length_one<DstType>())
{
	using namespace std;

	for(RandomIterator3 i = dst_first; i != dst_last; ++i)
		*i = -10;//-std::numeric_limits<DstType>::infinity();
	for(Iterator i = zero_first; i != zero_last; ++i) {
		*(dst_first+*i) = 0;
		*(prv_first+*i) = *i;
	}

	std::vector<char> status(dst_last-dst_first, 0); // 0: unvisit, 1:visited
	hj::algorithm::heap_index<RandomIterator3, DstType> heap(dst_first, dst_last);
	heap.make();
	while(!heap.empty()) {
		const size_t active_node = heap.top();
		status[active_node] = 1;
		heap.pop();
		for(RandomIterator2 i = idx_first + *(ptr_first+active_node);
			i != idx_first + *(ptr_first+active_node+1); ++i) { // for each neighbor node
			if(status[*i] == 1) continue;
			const DstType new_dst = *(dst_first + active_node) - edge_len(active_node, *i);
			if(new_dst > *(dst_first + *i)) {
				*(dst_first + *i) = new_dst;
				if(prv_first != prv_last)
					*(prv_first + *i) = active_node;
			}
			heap.update(*i);
		}
	}
}

#endif
