#include "binvox.h"

#include <fstream>
#include <iostream>
#include <string>
using namespace std;

int read_binvox(const char *path, binvox &dat)
{
	ifstream *input = new ifstream(path, ios::in | ios::binary);

	//
	// read header
	//
	string line;
	*input >> line;  // #binvox
	if (line.compare("#binvox") != 0) {
		cout << "Error: first line reads [" << line << "] instead of [#binvox]" << endl;
		delete input;
		return 0;
	}
	*input >> dat.version;
//	cout << "reading binvox version " << dat.version << endl;

	dat.depth = -1;
	int done = 0;
	while(input->good() && !done) {
		*input >> line;
		if (line.compare("data") == 0) done = 1;
		else if (line.compare("dim") == 0) {
			*input >> dat.depth >> dat.height >> dat.width;
		}
		else if (line.compare("translate") == 0) {
			*input >> dat.tx >> dat.ty >> dat.tz;
		}
		else if (line.compare("scale") == 0) {
			*input >> dat.scale;
		}
		else {
			cout << "  unrecognized keyword [" << line << "], skipping" << endl;
			char c;
			do {  // skip until end of line
				c = input->get();
			} while(input->good() && (c != '\n'));

		}
	}
	if (!done) {
		cout << "  error reading header" << endl;
		return 0;
	}
	if (dat.depth == -1) {
		cout << "  missing dimensions in header" << endl;
		return 0;
	}

	//
	// read voxel data
	//
	binvox::byte value, count;
	input->unsetf(ios::skipws);  // need to read every byte now (!)
	*input >> value;  // read the linefeed char

	while(input->good()) {
		*input >> value >> count;
		if (input->good())
			dat.entries.push_back(make_pair(value, count));
	}  // while

	input->close();

	return 0;
}
