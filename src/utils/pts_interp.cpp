#include <fstream>
#include <iostream>
#include <vector>
#include <iterator>
#include <stdlib.h>
#include <string>

extern "C" {
#include <hjlib/ANN_c.h>
}

using namespace std;

/**
   compute weight for interpolating value on pts dst from pts src
 */

int pts_interp(int argc, char *argv[])
{
	if(argc < 4) {
		cerr << "Usage: pts_interp pts-src pts-dst num-neighbor [output]" << endl;
		cerr << "if specified output, use binary format" << endl;
		return -1;
	}

	string str;
	size_t num;
	vector<double> pts;
	vector<double *> ppts;
	{
		ifstream ifs(argv[1]);
		ifs >> str >> str >> num;
		pts.resize(num*3);
		ppts.resize(num);
		for(size_t pi = 0; pi < num; ++pi) {
			ifs >> str;
			for(int d = 0; d < 3; ++d)
				ifs >> pts[pi*3+d];
			ppts[pi] = &pts[pi*3];
		}
	}
	void *ANNkd_tree_handle = ANNkd_tree_new(&ppts[0], num, 3);

	ifstream ifs(argv[2]);
	size_t dist_num;
	ifs >> str >> str >> dist_num;
	int num_neighbor = atoi(argv[3]);
	vector<int> idx(num_neighbor);
	vector<double> dist(num_neighbor);
	double xyz[3];

	bool is_bin_format = (argc == 5);
	ofstream ofs;
	if(is_bin_format) {
		ofs.open(argv[4], ofstream::binary);
		ofs.write((const char *)&dist_num, sizeof(size_t));
		ofs.write((const char *)&num_neighbor, sizeof(int));
	}
	else {
		cout << dist_num << " " << num_neighbor << endl;
	}
	while(ifs.good()) {
		ifs >> str >> xyz[0] >> xyz[1] >> xyz[2];
		if(ifs.fail())
			break;
		ANNkd_tree_search(ANNkd_tree_handle,
						  xyz, num_neighbor, &idx[0], &dist[0]);
		if(is_bin_format) {
			ofs.write((const char *)&idx[0], sizeof(int)*idx.size());
			ofs.write((const char *)&dist[0], sizeof(double)*idx.size());
		}
		else {
			copy(idx.begin(), idx.end(), ostream_iterator<int>(cout, " "));
			cout << '\n';
			copy(dist.begin(), dist.end(), ostream_iterator<double>(cout, " "));
			cout << '\n';
		}
	}

	ANNkd_tree_delete(ANNkd_tree_handle);
	return 0;
}
