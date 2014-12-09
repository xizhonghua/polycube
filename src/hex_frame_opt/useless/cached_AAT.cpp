#include "cached_AAT.h"

#include <algorithm>

using namespace std;
using namespace hj::sparse;
using namespace zjucad::matrix;

template <template <typename K, typename T> class MAP_TYPE>
class cached_patten_AAT
{
public:
	cached_patten_AAT():pass_(0), patten_size_(0) {}
	template<typename T1, typename INT_TYPE1, typename T2, typename INT_TYPE2>
	int operator()(const csc<T1, INT_TYPE1> &A, csc<T2, INT_TYPE2> &AAT,
				   resize_matrix_i_func rmi,
				   resize_matrix_d_func rmd) {
		if(pass_ == 0) {
			csc_by_vm<T2, INT_TYPE2, MAP_TYPE> vm(A.size(1), A.size(1), 0); // NOTICE nz per col
			assert(A.size(1) == vm.size(1) && A.size(1) == vm.size(2));
			const INT_TYPE2 n = A.size(2);
			INT_TYPE2 i, col_i = 0, row_i = 0;
			for(i = 0; i < n; ++i) {	// vector outer product: vvT(i, j) = vi*vj, for each v
				for(col_i = A.ptr_[i]; col_i < A.ptr_[i+1]; ++col_i) {	// for each nz column
					if(A.val_[col_i] == 0) continue;
					map_vec<T2, INT_TYPE2, MAP_TYPE> &col = vm[A.idx_[col_i]];
					for(row_i = A.ptr_[i]; row_i < A.ptr_[i+1]; ++row_i) {	// scalar * sparse, for each nz row
						col[A.idx_[row_i]] += A.val_[col_i]*A.val_[row_i];
						++patten_size_;
					}
				}
			}
			//convert(vm, AAT);
			AAT.rows_ = vm.size(1);
			const INT_TYPE2 cols = static_cast<INT_TYPE2>(vm.size(2));
			rmi(AAT.ptr_, cols+1);
			AAT.ptr_[0] = 0;
			for(i = 0; i < cols; ++i) // for each column
				AAT.ptr_[i+1] = AAT.ptr_[i]+static_cast<INT_TYPE2>(nnz(vm[i]));
			rmi(AAT.idx_, AAT.ptr_[cols]);
			rmd(AAT.val_,AAT.ptr_[cols]);

			typename map_vec<T1, INT_TYPE1, MAP_TYPE>::nz_const_iterator mi;
			for(i = 0; i < cols; ++i) {	// for each column
				const map_vec<T1, INT_TYPE1, MAP_TYPE> &col = vm[i];
				for(row_i = AAT.ptr_[i], mi = col.begin_nz(); row_i < AAT.ptr_[i+1]; ++row_i, ++mi) {
					AAT.idx_[row_i] = mi->first;
					AAT.val_[row_i] = mi->second;
				}
			}

			++pass_;
		}
		else if(pass_ == 1) {
			AAT.val_(colon()) = 0;
			assert(A.size(1) == AAT.size(1) && A.size(1) == AAT.size(2));
			const INT_TYPE2 n = A.size(2);
			INT_TYPE2 i, col_i = 0, row_i = 0;
			patten_.reserve(patten_size_);
			for(i = 0; i < n; ++i) {	// vector outer product: vvT(i, j) = vi*vj, for each v
				for(col_i = A.ptr_[i]; col_i < A.ptr_[i+1]; ++col_i) {	// for each nz column
					const INT_TYPE2 col_idx_of_AAT = A.idx_[col_i];
					const INT_TYPE2 nz_beg = AAT.ptr_[col_idx_of_AAT], nz_end = AAT.ptr_[col_idx_of_AAT+1];
					for(row_i = A.ptr_[i]; row_i < A.ptr_[i+1]; ++row_i) {	// scalar * sparse, for each nz row
						const INT_TYPE2 nz_of_AAT =
							lower_bound(&AAT.idx_[nz_beg], &AAT.idx_[nz_end], A.idx_[row_i])
							-&AAT.idx_[0];
						AAT.val_[nz_of_AAT] += A.val_[col_i]*A.val_[row_i];
						patten_.push_back(nz_of_AAT);
					}
				}
			}
			++pass_;
		}
		else {
			AAT.val_(colon()) = 0;
			assert(A.size(1) == AAT.size(1) && A.size(1) == AAT.size(2));
			const INT_TYPE2 n = A.size(2);
			INT_TYPE2 i, col_i = 0, row_i = 0;
			size_t j = 0;
			for(i = 0; i < n; ++i) {	// vector outer product: vvT(i, j) = vi*vj, for each v
				for(col_i = A.ptr_[i]; col_i < A.ptr_[i+1]; ++col_i) {	// for each nz column
					for(row_i = A.ptr_[i]; row_i < A.ptr_[i+1]; ++row_i, ++j) {	// scalar * sparse, for each nz row
						AAT.val_[patten_[j]] += A.val_[col_i]*A.val_[row_i];
					}
				}
			}
		}
	}
private:
	int pass_;
	size_t patten_size_;
	std::vector<size_t> patten_;
};

namespace hj { namespace sparse {

cached_AAT::cached_AAT()
{
	ctx_ = new cached_patten_AAT<map_by_sorted_vector>();
}

cached_AAT::~cached_AAT()
{
	cached_patten_AAT<map_by_sorted_vector> *this_
		= reinterpret_cast<cached_patten_AAT<map_by_sorted_vector> *>(ctx_);
	delete this_;
}

void cached_AAT::operator()(
		const csc<double, ptrdiff_t> &A, csc<double, ptrdiff_t> &AAT,
		resize_matrix_i_func rmi,
		resize_matrix_d_func rmd)
{
	cached_patten_AAT<map_by_sorted_vector> *this_
		= reinterpret_cast<cached_patten_AAT<map_by_sorted_vector> *>(ctx_);
	(*this_)(A, AAT, rmi, rmd);
}

}}
