#ifndef HJ_FAST_MARCHING_H_
#define HJ_FAST_MARCHING_H_

#include "../common/def.h"

// assume dist has been preallocated into node.size(2)
// alg = 0, = "dijkstra"
int fast_marching(const matrixd &node,
				  const matrixst &cell,
				  matrixd &dist,
				  matrixst *prev_node,
				  const char *alg = 0);

#endif
