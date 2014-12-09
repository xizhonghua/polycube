#include <string>
#include <fstream>
#include <iostream>
#include <vector>

#include "binvox.h"
#include "../common/IO_util.h"

using namespace std;

int flatten_to_ascii(const binvox &dat, const char *path)
{
	ofstream *out = new ofstream(path);
	if(!out->good()) {
		cout << "Error opening [voxels.txt]" << endl << endl;
		return (1);
	}

	cout << "Writing voxel data to ASCII file..." << endl;
  
	*out << "#binvox ASCII data" << endl;
	*out << "dim " << dat.depth << " " << dat.height << " " << dat.width << endl;
	*out << "translate " << dat.tx << " " << dat.ty << " " << dat.tz << endl;
	*out << "scale " << dat.scale << endl;
	*out << "data" << endl;

	int count = 0;
	for(int i = 0; i < dat.entries.size(); i++) {
		for(binvox::byte v = 0; v < dat.entries[i].second; ++v, ++count) {
			*out << (char) (dat.entries[i].first + '0') << " ";
			if (((count + 1) % dat.width) == 0) *out << endl;
		}
	}

	out->close();

	cout << "done" << endl << endl;
	return 0;
}

void binvox2pts(const binvox &dat, vector<double> &pts)
{
	grid_count gc(dat.width, dat.height, dat.depth);
	double pos[3];
	for(int i = 0; i < dat.entries.size(); i++) {
		if(!dat.entries[i].first) {
			gc += dat.entries[i].second;
			continue;
		}
		for(binvox::byte v = 0; v < dat.entries[i].second;
			++v, gc += 1) {
			counter2position(gc.yzx(), dat, pos);
			for(int d = 0; d < 3; ++d)
				pts.push_back(pos[d]);
		}
	}
}

int binvox2pts(int argc, char *argv[])
{
	if (argc < 2) {
		cout << "Usage: binvox2pts binvox [pts]" << endl << endl;
		return (1);
	}

	// load
	binvox dat;
	if (read_binvox(argv[1], dat)) {
		cerr << "Error reading [" << argv[1] << "]" << endl << endl;
		return (1);
	}

	// compute
	vector<double> pts;
	binvox2pts(dat, pts);

	//output
	if(argc < 3) {
		cout << "# node " << pts.size()/3 << '\n';
		for(size_t pi = 0; pi < pts.size(); pi+=3)
			cout << "v " << pts[pi] << ' ' << pts[pi+1] << ' ' << pts[pi+2] << '\n';
	}
	else {
		ofstream ofs(argv[2], ofstream::binary);
		if(ofs.fail()) {
			cerr << "open " << argv[2] << " fail." << endl;
			return __LINE__;
		}
		write_val(ofs, pts.size()/3);
		write_arr(ofs, &pts[0], pts.size());
	}
	return 0;
}
