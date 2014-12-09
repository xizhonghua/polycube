#include "rot_cubic_f_SH.h"

#include <cassert>
#include <vector>

#include <minpack.h>

using namespace std;

struct lmdif_prb {
	lmdif_prb(const double *current)
		:m(9), n(3), tol(1e-8), iwa(n), current_(current) {
		lwa = m*n+5*n+m;
		wa.resize(lwa);
	}
	int info, lwa, m, n;
	const double *current_;
	double tol;
	vector<int> iwa;
	vector<double> wa;
};

void lmdif_zyz2SH(int *m, int *n, double *x, double *fvec, int *iflag)
{
	assert(*m == 9 && *n == 3);
	const double *current = *reinterpret_cast<const double **>(n+1);
	calc_rot_cubic_f_sh_(fvec, x);
	for(int i = 0; i < 9; ++i)
		fvec[i] -= current[i];
}

int SH2zyz(double *zyz, const double *SH, unsigned int iter_num)
{
	vector<double> tmp(9);
	int info;
	for(size_t i = 0; i < iter_num; ++i) {
		lmdif_prb prb(SH);
		lmdif1_(lmdif_zyz2SH, &prb.m, &prb.n, zyz, &tmp[0], &prb.tol,
				&prb.info, &prb.iwa[0], &prb.wa[0], &prb.lwa);
		if(prb.info < 4) {
			info = prb.info;
			break;
		}
	}
	return info;
}

void lmder_zyz2SH(int *m, int*n, double *x, double *fvec, double *fjac, int *ldfjac, int* iflag)
{
	if(*iflag == 1)
		lmdif_zyz2SH(m, n, x, fvec, iflag);
	else if(*iflag == 2)
		calc_jac_rot_cubic_f_sh_(fjac, x);
}

int SH2zyz_lmder1(double *zyz, const double *SH, unsigned int iter_num)
{
	vector<double> fvec(9), fjac(9*3);
	vector<int> ipvt(3);
	int info;
	for(size_t i = 0; i < iter_num; ++i) {
		lmdif_prb prb(SH);
		lmder1_(lmder_zyz2SH, &prb.m, &prb.n, zyz, &fvec[0], &fjac[0], &prb.m,
				&prb.tol, &prb.info, &ipvt[0], &prb.wa[0], &prb.lwa);
		if(prb.info < 4) {
			info = prb.info;
			break;
		}
	}
	return info;
}
