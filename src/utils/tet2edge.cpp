#include <fstream>
#include <iostream>
#include <vector>
#include <algorithm>
using namespace std;

int tet2edge(int argc, char *argv[])
{
	vector<pair<size_t, size_t> > edges;

	size_t idx[4];

	while(cin) {
		for(size_t i = 0; i < 4; ++i)
			cin >> idx[i];
		if(cin.fail())
			break;
		sort(idx, idx+4);
		for(size_t i = 0; i < 4; ++i) {
			for(size_t j = i+1; j < 4; ++j) {
				edges.push_back(make_pair(idx[i], idx[j]));
			}
		}
	}

	sort(edges.begin(), edges.end());
	vector<pair<size_t, size_t> >::const_iterator end = unique(edges.begin(), edges.end());
	const size_t edge_num = end - edges.begin();

	cout << edge_num << '\n';
	for(size_t ei = 0; ei < edge_num; ++ei)
		cout << edges[ei].first << ' ' << edges[ei].second << '\n';

	return 0;
}
