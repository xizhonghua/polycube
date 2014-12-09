#include "Laplacian.h"

#include <iostream>
using namespace std;
using namespace zjucad::matrix;

#include "topology.h"

template <typename C>
static int find(const C &con, const typename C::value_type &v) {
	typename C::const_iterator i = std::find(con.begin(), con.end(), v);
	if(i == con.end()) return -1;
	return int(i-con.begin());
}

namespace hj { namespace mesh {

/// Laplacian
// 0 --> divide area
// 1 --> divide sqrt(area)
// 2 --> not divide
bool cot_Laplacian_val(const matrix<double> &v, const matrix<int> &t,
					const matrix<int> &ptr, const matrix<int> &idx, matrix<double> &val,
					int weighting)
{
    const int vn = ptr.size()-1;
    int i, j, k;
	val = zeros<double>(idx.size(), 1);

	matrix<double> edge[2], areas = zeros<double>(vn, 1);
    for(i = 0; i < t.size(2); ++i) {	// for each triangle
        matrix<int> tri = t(colon(), i);
        for(j = 0; j < 3; ++j) {
            for(k = 0; k < 2; ++k)
                edge[k] = v(colon(), tri[k])-v(colon(), tri[2]);
            areas[tri[0]] += norm(cross(edge[0], edge[1]));
            for(k = 0; k < 2; ++k) {
                double len = norm(edge[k]);
                if(len < 1e-5) {	// zero length edge
                    return false;
                }
                edge[k] /= len;
            }
            const double c = dot(edge[0], edge[1]);
            const double s = sqrt(1-c*c);
			if(s < 1e-10) {
				cout << "baquality mesh for coincide edges." << endl;
			}
            const double cot = c/s;
            for(k = 0; k < 2; ++k) {
                int pos = find(idx(colon(ptr[tri[k]], ptr[tri[k]+1]-1)), tri[1-k]);
                assert(pos >= 0);
                val[pos+ptr[tri[k]]] -= cot;
            }
            rotate(&tri[0], &tri[0]+1, &tri[0]+3);
        }
    }
	for(i = 0; i < vn; ++i) {
		int pos = find(idx(colon(ptr[i], ptr[i+1]-1)), i);
		assert(pos >= 0);
		val[pos+ptr[i]] = -sum(val(colon(ptr[i], ptr[i+1]-1)));
	}
	if(weighting == 0 || weighting == 1) {
        for(i = 0; i < vn; ++i)
			if(fabs(areas[i]) < 1e-10) {
				cout << "bad quality mesh for local area." << endl;
			}
	}

    if(weighting == 0) {
        for(i = 0; i < vn; ++i)
            val(colon(ptr[i], ptr[i+1]-1)) /= areas[i];
    }
    else if(weighting == 1) {
        for(i = 0; i < vn; ++i)
            val(colon(ptr[i], ptr[i+1]-1)) /= sqrt(areas[i]);
    }
	return true;
}

bool cot_Laplacian(const matrix<double> &v, const matrix<int> &t,
				   matrix<int> &ptr, matrix<int> &idx, matrix<double> &val,
				   int weighting)
{
    cell2csc(t, ptr, idx, true);
    return cot_Laplacian_val(v, t, ptr, idx, val, weighting);
}

}}
