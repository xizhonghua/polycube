#include <fstream>
#include <iostream>
#include <vector>
#include <algorithm>
using namespace std;

#include "../common/IO.h"
extern "C" {
#include <hjlib/ANN_c.h>
}

#include "../tetmesh/tetmesh.h"
using namespace zjucad::matrix;

int build_nearest(const matrixd &src_pts, const matrixd &dst_pts,
				  matrix<uint32_t> &nearest)
{
	const size_t DIM = src_pts.size(1);
	assert(dst_pts.size(1) == DIM);
	vector<const double *> ppts(src_pts.size(2));
	for(size_t pi = 0; pi < ppts.size(); ++pi)
		ppts[pi] = &src_pts(0, pi);
	void *ANNkd_tree_handle = ANNkd_tree_new(const_cast<double **>(&ppts[0]), ppts.size(), DIM);
	for(size_t dpi = 0; dpi < dst_pts.size(2); ++dpi) {
		int idx;
		double dist2;
		ANNkd_tree_search(ANNkd_tree_handle,
						  const_cast<double *>(&dst_pts(0, dpi)), 1, &idx, &dist2);
		nearest[dpi] = idx;
	}
	ANNkd_tree_delete(ANNkd_tree_handle);
	return 0;
}

int build_pts_from_tet(const char *path, matrixd &pts)
{
	jtf::mesh::meshes tm;
	if(jtf::mesh::tet_mesh_read_from_zjumat(path, &tm.node_, &tm.mesh_))
		return __LINE__;
	unique_ptr<jtf::mesh::face2tet_adjacent> fa(jtf::mesh::face2tet_adjacent::create(tm.mesh_));
	pts = zeros<double>(3, fa->faces_.size());
	for(size_t fi = 0; fi < fa->faces_.size(); ++fi) {
		const std::vector<size_t> &facei = fa->faces_[fi];
		for(size_t ni = 0; ni < 3; ++ni)
			pts(colon(), fi) += tm.node_(colon(), facei[ni]);
	}
	pts /= 3;
	return 0;
}

int interp(int argc, char *argv[], size_t dim, const matrix<uint32_t> &nearest)
{
	for(int ai = 0; ai+1 < argc; ai += 2) {
		matrixd src_val;
		if(jtf::mesh::read_matrix(argv[ai], src_val))
			return __LINE__;
		if(src_val.size(1) != dim) {
			matrixd tmp = src_val;
			src_val.resize(dim, tmp.size()/dim);
			if(src_val.size() != tmp.size()) {
				cerr << "incompatible value file: " << dim << " " << tmp.size(2) << endl;
				return __LINE__;
			}
			copy(tmp.begin(), tmp.end(), src_val.begin());
		}
		matrixd dst_val = src_val(colon(), nearest);
		if(jtf::mesh::write_matrix(argv[ai+1], dst_val))
			return __LINE__;
	}
	return 0;
}

int interp(int argc, char *argv[])
{
	if(argc < 3) {
		cerr << "interp src_tet dst_tet nearest [dim value output]*" << endl;
		return __LINE__;
	}

	matrix<uint32_t> nearest;
	int has_nearest = (jtf::mesh::read_matrix(argv[3], nearest) == 0);

	if(!has_nearest) {
		matrixd src_pts;
		if(build_pts_from_tet(argv[1], src_pts))
			return __LINE__;
		matrixd dst_pts;
		if(build_pts_from_tet(argv[2], dst_pts))
			return __LINE__;
		nearest.resize(dst_pts.size(2));
		if(build_nearest(src_pts, dst_pts, nearest))
			return __LINE__;
		if(jtf::mesh::write_matrix(argv[3], nearest))
			return __LINE__;
	}
	if(argc > 6)
		return interp(argc-5, argv+5, atoi(argv[4]), nearest);
	return 0;
}
