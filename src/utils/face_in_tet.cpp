#include <fstream>
#include <iostream>
#include <vector>
#include <algorithm>

#include "../common/IO.h"
#include "../tetmesh/tetmesh.h"

using namespace zjucad::matrix;
using namespace std;

int face_in_tet(int argc, char *argv[])
{
	if(argc < 2) {
		cerr << "face_in_tet tet" << endl;
		return __LINE__;
	}

	jtf::mesh::meshes tm;
	if(jtf::mesh::tet_mesh_read_from_zjumat(
		   argv[1],
		   &tm.node_, &tm.mesh_))
		return __LINE__;

	unique_ptr<jtf::mesh::face2tet_adjacent> fa(jtf::mesh::face2tet_adjacent::create(tm.mesh_));

	const size_t fn = fa->faces_.size();
	cout << fn << "\n";
	for(size_t fi = 0; fi < fn; ++fi) {
		cout << fa->faces_[fi][0] << ' '
			 << fa->faces_[fi][1] << ' '
			 << fa->faces_[fi][2] << ' '
			 << fa->face2tet_[fi].first << ' '
			 << fa->face2tet_[fi].second << '\n';
	}
	return 0;
}
