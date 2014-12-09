#include <fstream>
#include <iostream>
#include <vector>
#include <algorithm>

#include <zjucad/matrix/io.h>

#include "../tetmesh/tetmesh.h"
#include "../common/util.h"
#include "../common/IO.h"
#include "../common/vtk.h"

#include "../fast_marching/fast_marching.h"

using namespace zjucad::matrix;
using namespace std;

int fast_marching(int argc, char *argv[])
{
	if(argc < 3) {
		cerr << "fast_marching tet dst [vtk] [sphere (for validation)]" << endl;
		return __LINE__;
	}

	jtf::mesh::meshes tm;
	if(jtf::mesh::tet_mesh_read_from_zjumat(
		   argv[1],
		   &tm.node_, &tm.mesh_, 0)) {
		cerr << "# read tet fail." << endl;
		return __LINE__;
	}

	matrixd dist(tm.node_.size(2));
	matrixst prev_node(tm.node_.size(2));
	int rtn = fast_marching(tm.node_, tm.mesh_, dist, &prev_node);
	if(rtn) return rtn;

	{
		ofstream ofs(argv[2], ofstream::binary);
		if(ofs.fail()) {
			cerr << "# write dist fail." << endl;
			return __LINE__;
		}
		jtf::mesh::write_matrix(ofs, dist);
		jtf::mesh::write_matrix(ofs, prev_node);
	}

	if(argc > 3) { // vtk
		ofstream ofs(argv[3]);
		point_data(ofs, dist.begin(), dist.size(), "depth");
	}

	if(argc < 5)
		return 0;

	// validate by sphere
	matrixd bb(2, 3);
	calc_bounding_box(tm.node_, &bb[0]);
	bb = temp(trans(bb));
	const matrixd ct = bb*ones<double>(2, 1)/2.0;
	const matrixd r0 = (bb(colon(), 1) - bb(colon(), 0))/2.0;
	const double r = sum(r0)/3.0;
	cout << "# r0: " << trans(r0) << trans(ct) << " " << r << endl;

	cout << min(dist) << " " << max(dist) << endl;
	for(size_t ni = 0; ni < dist.size(); ++ni)
		dist[ni] -= norm(tm.node_(colon(), ni)-ct)-r;

	cout << min(dist) << " " << max(dist) << endl;

	return 0;
}
