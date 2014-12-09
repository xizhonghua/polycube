#include "../common/obj_like_format.h"

#include <fstream>
#include <iostream>
#include <string>
#include <sstream>

#include <zjucad/matrix/matrix.h>
#include "../common/def.h"
using namespace std;
using namespace zjucad::matrix;

template <typename T>
static int read_tetgen_output(
	const char *path, size_t ele_per_line, matrix<T> &dat)
{
	string line;

	ifstream ifs(path);

	size_t num;
	{
		getline(ifs, line);
		istringstream iss(line.c_str());
		iss >> num;
	}

	dat.resize(ele_per_line, num);
	size_t idx;
	for(size_t i = 0; i < num; ++i) {
		getline(ifs, line);
		istringstream iss(line.c_str());
		iss >> idx;
		for(int ei = 0; ei < ele_per_line; ++ei)
			iss >> dat(ei, i);
	}
	return ifs.fail();
}

void orient_tet(const matrixd &node, matrix<int> &tet)
{
	matrixd tmp(3);
	for(int ti = 0; ti < tet.size(2); ++ti) {
		tmp = cross(node(colon(), tet(1, ti)) - node(colon(), tet(0, ti)),
					node(colon(), tet(2, ti)) - node(colon(), tet(0, ti)));
		const double vol = dot(tmp, node(colon(), tet(3, ti)) - node(colon(), tet(0, ti)));
		if(vol < 0) {
			cerr << "find inversed tet " << ti << endl;
			swap(tet(0, ti), tet(1, ti));
		}
	}
}

int tetgen2tet(int argc, char *argv[])
{
	if(argc < 2) {
		cerr << "Usage: tetgen2tet prefix" << endl;
		return __LINE__;
	}

	string prefix = argv[1];
	string path[3] = {
		prefix + ".1.node",
		prefix + ".1.face",
		prefix + ".1.ele"
	};
	matrixd node;
	read_tetgen_output(path[0].c_str(), 3, node);

	matrix<int> face, tet;
	read_tetgen_output(path[1].c_str(), 3, face);
	read_tetgen_output(path[2].c_str(), 4, tet);
	tet -= 1;

	orient_tet(node, tet);
	tet += 1;

	write_section(cout, node.size(2), 3, &node[0], "node", "v");
	write_section(cout, face.size(2), 3, &face[0], "face", "f");
	write_section(cout, tet.size(2), 4, &tet[0], "tet", "tet");
	return 0;
}
