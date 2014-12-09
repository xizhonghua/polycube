#include <fstream>
#include <iostream>
#include <vector>
#include <algorithm>

#include "../common/obj_like_format_mat.h"
#include "../common/zyz.h"
#include "../common/transition.h"
#include "../common/graph.h"
#include "../common/IO.h"

#include <zjucad/matrix/io.h>
using namespace std;
using namespace zjucad::matrix;

void orient_uvw(
	const matrixd &pos,
	const fix_graph &fg, const matrix<int> &tet, const matrixd &frame_at_node,
	const matrixd &uvw_at_node, const matrixd &trans_rot, const matrixd &trans_offset,
	matrixd &uvw_in_tet, matrixd &frame_in_tet)
{
	const int tn = tet.size(2);
	uvw_in_tet.resize(3, tn*4);
	frame_in_tet.resize(9, tn);
	const int REF_NODE = 0;
	matrixd rot(3, 3), I = eye<double>(3);
	for(int ti = 0; ti < tn; ++ti) { // for each tet
		frame_in_tet(colon(), ti) = frame_at_node(colon(), tet(REF_NODE, ti));
		for(int ni = 0; ni < 4; ++ni) {
			if(ni == REF_NODE) {
				rot = eye<double>(3);
				uvw_in_tet(colon(), ti*4+ni) = uvw_at_node(colon(), tet(ni, ti));
			}
			else {
				const size_t edge_idx = fg.find_edge(tet(REF_NODE, ti), tet(ni, ti));
				rot(colon()) = trans_rot(colon(), edge_idx);

				uvw_in_tet(colon(), ti*4+ni) =
					rot*uvw_at_node(colon(), tet(ni, ti))+trans_offset(colon(), edge_idx);
			}
		}
#if 1 // check
		const size_t CHECK_TET = 0;
		if(1 || ti == CHECK_TET) {
			matrixd proj, edge, uvw_offset;
			rot(colon()) = frame_in_tet(colon(), ti);
			for(int ni = 0; ni < 4; ++ni) {
				if(ni == 0) continue;
				edge = pos(colon(), tet(ni, ti))-pos(colon(), tet(REF_NODE, ti));
				proj = trans(rot)*edge/0.17706;
				uvw_offset = (uvw_in_tet(colon(), ti*4+ni)-uvw_in_tet(colon(), ti*4+REF_NODE));

				const size_t edge_idx = fg.find_edge(tet(REF_NODE, ti), tet(ni, ti));
				const double wrong = norm(trans(uvw_offset) - trans(proj));
				if(wrong > 1) {
					cerr << "check edge: " << tet(REF_NODE, ti) << "<-" << tet(ni, ti) << endl;
					cerr << "diff at edge " << edge_idx << " is: " << wrong << "\n"
						 << trans(uvw_offset) << trans(proj)
						 << trans(uvw_offset) - trans(proj);

					cerr << "node: " << trans(uvw_at_node(colon(), tet(REF_NODE, ti)))
						 << trans(uvw_at_node(colon(), tet(ni, ti)));
					cerr << "tet:  " << trans(uvw_in_tet(colon(), ti*4+REF_NODE))
						 << trans(uvw_in_tet(colon(), ti*4+ni));

					cerr << "transition is: " << trans(trans_rot(colon(), edge_idx))
						 << trans(trans_offset(colon(), edge_idx));

					rot(colon()) = trans_rot(colon(), edge_idx);
					cerr << "recompute: " << rot << trans_offset(colon(), edge_idx)
						 << rot*uvw_at_node(colon(), tet(ni, ti))+trans_offset(colon(), edge_idx);
				}
			}
		}
#endif
	}
}

int split_tet(int argc, char *argv[])
{
	if(argc < 2)
		cerr << "Usage: plit_tet tet [vg vf param vf_in_tet] " << endl;

	ifstream ifs(argv[1]);
	matrixd pos;
	matrix<int> tet;
	if(read_section(ifs, 3, pos))
		return __LINE__;
	if(read_section(ifs, 3, tet))
		return __LINE__;
	if(read_section(ifs, 4, tet))
		return __LINE__;

	static const int tet2tri[4][3] = {
		0, 1, 2,
		1, 3, 2,
		2, 3, 0,
		3, 1, 0
	};

	// output split into obj
	// for(int vi = 0; vi < pos.size(2); ++vi)
	// 	cout << "v " << pos(0, vi)
	// 		 << ' ' << pos(1, vi)
	// 		 << ' ' << pos(2, vi) << '\n';
	// for(int ti = 0; ti < tet.size(2); ++ti) {
	// 	for(int fi = 0; fi < 4; ++fi) {
	// 		cout << "f " << tet(tet2tri[fi][0], ti)
	// 			 << ' ' << tet(tet2tri[fi][1], ti)
	// 			 << ' ' << tet(tet2tri[fi][2], ti) << '\n';
	// 	}
	// }

	// split param into tet
	tet -= 1;
	if(argc >= 5) {
		fix_graph fg;
		{
			matrixd node_pos;
			if(load_graph(argv[2], fg, node_pos))
				return __LINE__;
			cerr << "num of edges in the vg: " << fg.edge_num() << endl;
		}
		matrixd frame_at_node(9, pos.size(2));
		{
			ifstream ifs(argv[3]);
			matrixd zyz, I = eye<double>(3);
			if(read_section(ifs, 3, zyz))
				return __LINE__;
			if(zyz.size(2) != frame_at_node.size(2))
				return __LINE__;
			for(int ni = 0; ni < frame_at_node.size(2); ++ni) {
				zyz_angle_2_rotation_matrix(
					zyz(0, ni), zyz(1, ni), zyz(2, ni),
					&frame_at_node(0, ni));
//				frame_at_node(colon(), ni) = I(colon());
			}
		}
		matrixd uvw_at_node, trans_rot, trans_offset;
		{
			ifstream ifs(argv[4]); // param
			if(read_section(ifs, 3, uvw_at_node))
				return __LINE__;
			if(uvw_at_node.size(2) != pos.size(2))
				return __LINE__;

			if(read_section(ifs, 9, trans_rot))
				return __LINE__;
			if(trans_rot.size(2) != fg.edge_num())
				return __LINE__;

			if(read_section(ifs, 3, trans_offset))
				return __LINE__;
			if(uvw_at_node.size(2) != pos.size(2))
				return __LINE__;
		}

		matrixd uvw_in_tet, frame_in_tet;
		orient_uvw(pos, fg, tet, frame_at_node,
				   uvw_at_node, trans_rot, trans_offset,
				   uvw_in_tet, frame_in_tet);

		cerr << "begin output to: " << argv[5] << endl;
		ofstream ofs(argv[5]);
		write_section(ofs, uvw_in_tet, "uvw", "uvw");
		write_section(ofs, frame_in_tet, "frame", "rot");
	}

	return 0;
}
