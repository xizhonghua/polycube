/**
   \brief surface frame to boundary condition
 */

#include <iostream>
#include <string>
#include <sstream>
using namespace std;

int sf2bc(int argc, char *argv[])
{
	string line, word;
	size_t num;
	getline(cin, line);
	istringstream iss(line.c_str());
	iss >> word >> word >> num;
	cout << "# boundary_condition " << num << '\n';
	for(size_t ni = 0; ni < num; ++ni) {
		if(cin.fail()) {
			cerr << "bad sf file." << endl;
			return __LINE__;
		}
		cout << ni << ' ' << 1 << ' ' << 3 << '\n';
		for(int d = 0; d < 3; ++d) {
			getline(cin, line);
			cout << line << '\n';
		}
	}

	return 0;
}
