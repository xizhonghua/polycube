#include "fast_marching.h"

#include <limits>
#include <memory>
#include <iostream>

#include "../common/graph.h"
#include "../tetmesh/tetmesh.h"

#include "dijkstra.h"

using namespace std;
using namespace zjucad::matrix;

void cell2fix_graph(const matrixd &node, const matrixst &cell, fix_graph &fg)
{
	dyn_graph dg(node.size(2));
	for(size_t ci = 0; ci < cell.size(2); ++ci) {
		for(size_t i = 0; i < cell.size(1); ++i) {
			for(size_t j = 0; j < cell.size(1); ++j) {
				if(i == j) continue;
				dg[cell(i, ci)].insert(cell(j, ci));
			}
		}
	}
	convert(dg, fg);
}

class edge_length
{
public:
	edge_length(const matrixd &node)
		:_node(node) {
	}
	double operator()(size_t a, size_t b) const {
		return norm(_node(colon(), a)-_node(colon(), b));
	}
private:
	matrixd _node;
};

int fast_marching(const matrixd &node,
				  const matrixst &cell,
				  matrixd &dist,
				  matrixst *prev_node,
				  const char *alg)
{
	fix_graph fg;
	cell2fix_graph(node, cell, fg);

	unique_ptr<jtf::mesh::face2tet_adjacent> fa(jtf::mesh::face2tet_adjacent::create(cell));
	if(!fa.get())
		return __LINE__;
//	if(validate_face2tet_adjacent(*fa, cell, &node)) {
//		cerr << "# fa error" << endl;
//		return __LINE__;
//	}

	matrixst face_idx;
	get_outside_face_idx(*fa, face_idx);
	set<size_t> surface_nodes;
	for(size_t fi = 0; fi < face_idx.size(); ++fi)
		for(size_t i = 0; i < 3; ++i)
			surface_nodes.insert(fa->faces_[face_idx[fi]][i]);
	cout << "# " << face_idx.size() << " surface triangles." << endl;
	cout << "# " << surface_nodes.size() << " surface nodes in " << node.size(2) << endl;

	dijkstra<double>(fg.ptr_.begin(), fg.ptr_.end(),
					 fg.idx_.begin(), fg.idx_.end(),
					 surface_nodes.begin(), surface_nodes.end(),
					 dist.begin(), dist.end(),
					 prev_node?prev_node->begin():0, prev_node?prev_node->end():0,
					 edge_length(node));

	return 0;
}
