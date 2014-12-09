#include "transition.h"

using namespace zjucad::matrix;
using namespace std;

double calc_rot_diff(const matrixd &R0, const matrixd &R1)
{
	return norm(R0-R1);
}

void shuffle(const matrixd &F, int type,
			 matrixd &sF)
{
	int shuffle[6];
	type2shuffle(type, shuffle);
	for(int i = 0; i < 3; ++i)
		sF(colon(), i) = F(colon(), shuffle[i*2])*shuffle[i*2+1];
#if 0 //check
	double e[2] = {
		norm(sF*trans(sF)-eye<double>(3)),
		norm(cross(sF(colon(), 0), sF(colon(), 1))-sF(colon(), 2))
	};
	if(e[0] > 1e-10 || e[1] > 1e-10)
		cout << "error: " << e[0] << ", " << e[1] << endl;
#endif
}

void calc_rot_diff_shuffle(
	const double *R0, const double *R1,
	vector<pair<double, int> > &error)
{
	matrixd R[2] = {zeros<double>(3, 3), zeros<double>(3, 3)},
		diff(3, 3), sR1(3, 3);
	copy(R0, R0+9, R[0].begin());
	copy(R1, R1+9, R[1].begin());
	error.resize(24);
	for(int type = 0; type < 24; ++type) {
		shuffle(R[1], type, sR1);
		error[type].first = calc_rot_diff(R[0], sR1);
		error[type].second = type;
	}
	sort(error.begin(), error.end());
}

void get_best_alignment(
	const double* from, const double* to,
	double* rot)
{
	std::vector<std::pair<double, int> > error;
    calc_rot_diff_shuffle(to, from, error);
    matrixd R(3, 3);
    shuffle(eye<double>(3), error[0].second, R);
	copy(R.begin(), R.end(), rot);
}
