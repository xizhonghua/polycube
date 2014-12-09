#ifndef HJ_CACHED_AAT_H_
#define HJ_CACHED_AAT_H_

#include <hjlib/sparse/sparse.h>

#ifndef cached_AAT_API
#define cached_AAT_API
#endif

namespace hj { namespace sparse {

typedef void (*resize_matrix_i_func)(zjucad::matrix::matrix<ptrdiff_t> &m, size_t size);
typedef void (*resize_matrix_d_func)(zjucad::matrix::matrix<double> &m, size_t size);

inline void resize_matrix_i(zjucad::matrix::matrix<ptrdiff_t> &m, size_t size) {
	m.resize(size);
}
inline void resize_matrix_d(zjucad::matrix::matrix<double> &m, size_t size) {
	m.resize(size);
}

class cached_AAT_API cached_AAT
{
public:
	cached_AAT();
	~cached_AAT();
	void operator()(
		const csc<double, ptrdiff_t> &A, csc<double, ptrdiff_t> &AAT,
		resize_matrix_i_func rmi = resize_matrix_i,
		resize_matrix_d_func rmd = resize_matrix_d);
private:
	void *ctx_;
};

}}

#endif
