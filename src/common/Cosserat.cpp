#include <zjucad/matrix/matrix.h>
#include <zjucad/matrix/io.h>
using namespace zjucad::matrix;

#include <algorithm>
#include <iostream>
using namespace std;

#include <hjlib/math/quaternion.h>
using namespace hj::quaternion;

#include "../common/def.h"

// TODO: maybe using quaternion is better

// t0 and t1 must be normalized vector
int tangent_to_Bishop_rotation(const matrixd &t0,
							   const matrixd &t1,
							   matrixd &R)
{
	matrixd e, Fe[2] = {zeros<double>(3, 3), zeros<double>(3, 3)};
	const double eps = 1e-8;

	e = t0-t1;
	if(norm(e) < eps) {
		R = eye<double>(3);
		return 0;
	}

	e = t0+t1;
	if(norm(e) < eps) {
		return 1;
	}

	Fe[0](colon(), 0) = t0;
	Fe[1](colon(), 0) = t1;

	matrixd u = cross(t0, t1);
	if(norm(u) < eps) {
		cerr << "# nearly parallel." << endl;
	}
	u /= norm(u);
	Fe[0](colon(), 1) = u;
	Fe[1](colon(), 1) = u;

	Fe[0](colon(), 2) = cross(t0, u);
	Fe[1](colon(), 2) = cross(t1, u);

	R = Fe[1]*trans(Fe[0]);

	return 0;
}

// pos : 3*n
// R : n-2
// R[i] : 3*3
int position_to_Bishop_rotation(const matrixd &pos,
								matrix<matrixd > &R)
{
	const int n = pos.size(2);
	int last_err_node = 0;
	matrixd t[2];
	for(int pi = 1; pi < n-1; ++pi) {
		for(int i = 0; i < 2; ++i) {
			t[i] = pos(colon(), pi+i)-pos(colon(), pi+i-1);
			const double len = norm(t[i]);
			const double eps = 1e-8;
			if(len < eps) {
				return -pi;
			}
			t[i] /= len;
		}
		if(tangent_to_Bishop_rotation(t[0], t[1], R[pi-1]))
			last_err_node = pi;
	}
	return last_err_node;
}

// R : n
// R[i] : 3*3
// u : n+1
// u[i]: 3*n
// u[0] must be initialized
void parallel_transport(const matrix<matrixd > &R,
					   matrix<matrixd > &u)
{
	const int n = R.size();
	for(int i = 1; i < n+1; ++i)
		u[i] = R[i-1]*u[i-1];
}

// pos : 3*n
// Fe : n-1;
// Fe[i] : 3*3
int position_to_Bishop_frame_at_edge(const matrixd &pos,
									  matrix<matrixd > &Fe)
{
	const int n = pos.size(2);
	matrix<matrixd > R(n-2);

	int err = position_to_Bishop_rotation(pos, R);
	if(err)	return err;

	parallel_transport(R, Fe);

	return 0;
}

void find_orthogonal_vector(const matrixd &t,
							matrixd &u)
{
	matrixd fabs_t = fabs(t);
	const int axis = min_element(fabs_t.begin(), fabs_t.end())-fabs_t.begin();
	matrixd &v = fabs_t;
	v(colon()) = 0;
	v[axis] = 1;
	u = cross(t, v);
}

// pos : 3*n, n >= 2
// Fe : 3*3
int init_first_frame(const matrixd &pos,
					 matrixd &Fe)
{
	if(pos.size(2) < 2)
		return -1;

	matrixd t = pos(colon(), 1)-pos(colon(), 0);
	double len = norm(t);
	const double eps = 1e-8;
	if(len < eps)
		return 1;
	t /= len;
	Fe(colon(), 0) = t;

	matrixd u(3);
	find_orthogonal_vector(t, u);
	len = norm(u);
	if(len < eps)
		return 2;
	u /= len;
	Fe(colon(), 1) = u;

	Fe(colon(), 2) = cross(Fe(colon(), 0), Fe(colon(), 1));
	if(norm(Fe*trans(Fe)-eye<double>(3)) > 1e-8) {
		cout << "# init first frame error." << Fe << Fe*trans(Fe)-eye<double>(3) << endl;
	}
	return 0;
}

void interp_frame(double a, const matrixd &R0,
				  double b, const matrixd &R1,
				  matrixd &R)
{
	matrixd q[2] = {zeros<double>(4, 1), zeros<double>(4, 1)}, Rq[2];
	m332quat<double>(R0, q[0]);
	m332quat<double>(R1, q[1]);
	Rq[0] = a*q[0]+b*q[1];
	Rq[1] = -a*q[0]+b*q[1];
	const double len[2] = {norm(Rq[0]), norm(Rq[1])};
	const int best = len[0] < len[1];
	if(len[best] < 1e-8) {
		cerr << "# degenerated quat." << endl;
	}
	Rq[best] /= norm(Rq[best]);
	quat2m33<double>(Rq[best], R);
}

// interpolation using piecewise linear basis
int edge_frame_to_node_frame(const matrix<matrixd > &Fe,
							 matrix<matrixd > &Fn)
{
// 	for(int j = 0; j < 6; ++j)
// 		cout << j << " " << Fe[j] << endl;
	const int n = Fe.size()+1;
	assert(n == Fn.size());
	if(n == 2) {
		Fn[0] = Fe[0];
		Fn[1] = Fe[0];
		return 0;
	}
	interp_frame(1.5, Fe[0], -0.5, Fe[1], Fn[0]);
	for(int i = 1; i < n-1; ++i) // for each node
		interp_frame(0.5, Fe[i-1], 0.5, Fe[i], Fn[i]);
	interp_frame(-0.5, Fe[n-3], 1.5, Fe[n-2], Fn[n-1]);

// 	for(int j = 0; j < 6; ++j)
// 		cout << j << " " << Fe[j] << endl;
	return 0;
}

// pos : 3*n
// Fn : n
// Fn[i] : 3*3
int position_to_Bishop_frame_at_node(const matrixd &pos,
									 matrix<matrixd > &Fn)
{
	const int n = pos.size(2);
	matrix<matrixd > Fe(n-1);
	fill(Fe.begin(), Fe.end(), zeros<double>(3, 3));
	init_first_frame(pos, Fe[0]);
	if(position_to_Bishop_frame_at_edge(pos, Fe))
		return 1;
	if(edge_frame_to_node_frame(Fe, Fn))
		return 2;
	return 0;
}
