#include "IO.h"

#include "../common/zyz.h"
#include "../common/obj_like_format_mat.h"

#include <fstream>
#include <iostream>

using namespace std;
using namespace zjucad::matrix;

int load_graph(const char *path, fix_graph &fg, matrixd &node_pos)
{
	ifstream ifs(path);
	if(ifs.fail()) {
		cerr << "read file error." << endl;
		return __LINE__;
	}
	if(read_section(ifs, 3, node_pos, "node")) {
		return __LINE__;
	}

    matrixst edge;
	if(read_section(ifs, 2, edge, "edge")) {
		return __LINE__;
	}
	edge -= 1;

	dyn_graph dg(node_pos.size(2));
	
	int beg, end;
	for(int ei = 0; ei < edge.size(2); ++ei) {
		beg = edge(0, ei);
		end = edge(1, ei);
		dg[beg].insert(end);
		dg[end].insert(beg);
	}
	convert(dg, fg);

	return 0;
}

int load_vec_field(const char *path, matrixd &rot)
{
	ifstream ifs(path);
	if(ifs.fail()) {
		cerr << "read file error." << endl;
		return __LINE__;
	}
    matrixd zyz;
        if(read_section(ifs, 3, zyz, "frame_as_zyz"))
                return __LINE__;
        if(zyz.size(2) != rot.size(2)) {
                cerr << "incompatible zyz file." << endl;
                return __LINE__;
        }
        for(int ni = 0; ni < zyz.size(2); ++ni) {
                zyz_angle_2_rotation_matrix(
                        zyz(0, ni), zyz(1, ni), zyz(2, ni),
                        &rot(0, ni));
        }
        return 0;
}

int save_param(const char *path,
               const matrixd &uvw,
               const matrixd &trans_rot,
               const matrixd &trans_offset)
{
	ofstream ofs(path);
	if(ofs.fail()) {
		cerr << "open " << path << " fail." << endl;
		return __LINE__;
	}
	write_section(ofs, uvw, "uvw", "uvw");
	write_section(ofs, trans_rot, "trans_rot", "tr");
	write_section(ofs, trans_offset, "trans_offset", "to");
	return 0;
}
