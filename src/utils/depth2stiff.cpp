#include <fstream>
#include <iostream>
#include <vector>
#include <algorithm>
#include <string.h>

#include <zjucad/matrix/io.h>

#include "../tetmesh/tetmesh.h"
#include "../common/util.h"
#include "../common/IO.h"
#include "../common/vtk.h"

using namespace zjucad::matrix;
using namespace std;

inline double depth2stiff(double depth, const char *type, double coef) {
	if(!strcmp(type, "linear")) {
		return (1-coef)*depth+1;
	}
	if(!strcmp(type, "exp")) {
		const double min = 0.1;
		return exp(-depth*depth*coef)*(1-min)+min;
	}
	return 1;
}

int depth2stiff(int argc, char *argv[])
{
	if(argc < 6) {
		cerr << "depth2stiff tet depth stiff type coef [vtk]" << endl;
		return __LINE__;
	}

	jtf::mesh::meshes tm;
	if(jtf::mesh::tet_mesh_read_from_zjumat(
		   argv[1],
		   &tm.node_, &tm.mesh_, 0)) {
		cerr << "# read tet fail." << endl;
		return __LINE__;
	}

	matrixd depth;
	ifstream ifs(argv[2], ifstream::binary);
	if(ifs.fail()) {
		cerr << "# open depth fail." << endl;
		return __LINE__;
	}
	jtf::mesh::read_matrix(ifs, depth);

	// normalize depth to [-1, 0]
	double range = min(depth);
	cerr << "# range: " << range << " " << max(depth) << endl;
	depth /= fabs(range);

	matrixd stiff(tm.mesh_.size(2), 1);
	const double coef = atof(argv[5]);
	for(size_t ti = 0; ti < stiff.size(); ++ti) {
		const double tet_d = sum(depth(tm.mesh_(colon(), ti)))/4.0;
		stiff[ti] = depth2stiff(tet_d, argv[4], coef);
	}

	{
		ofstream ofs(argv[3], ofstream::binary);
		if(ofs.fail()) {
			cerr << "# write stiff fail." << endl;
			return __LINE__;
		}
		jtf::mesh::write_matrix(ofs, stiff);
	}

	if(argc > 6) {
		cell_data(cout, stiff.begin(), stiff.size(), "stiff");
	}
	return 0;
}
