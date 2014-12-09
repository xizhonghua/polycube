class smooth_function : public sym_frame_opt_function
{
public:
	smooth_function(sym_frame_opt_sys &sys, size_t beg, size_t end)
		:sym_frame_opt_function(sys), beg_(beg), end_(end) {
		assert(beg_ < end_); // ensure the csc is sorted
	}
	virtual size_t dim_of_x(void) const {
		return sys_.sh_.size(2)*3;
	}
	virtual size_t dim_of_f(void) const {
		return 9;
	}
	virtual int val(const double *x, double *f) const {
		sys_.val(x);
		itr_matrix<double *> f0(9, 1, f);
		f0 = sys_.sh_(colon(), beg_)-sys_.sh_(colon(), end_);
		return 0;
	}
	virtual int jac(const double *x, double *val, function::idx_type *ptr = 0, function::idx_type *idx = 0) const {
		sys_.jac(x);
		for(function::idx_type i = 0; i < dim_of_f(); ++i) {
			ptr[i+1] = ptr[i] + 6;
			for(function::idx_type d = 0; d < 3; ++d) {
				idx[ptr[i]+d] = beg_*3+d;
				val[ptr[i]+d] = sys_.jac_sh_[beg_](i, d);
				idx[ptr[i]+d+3] = end_*3+d;
				val[ptr[i]+d+3] = -sys_.jac_sh_[end_](i, d);
			}
		}
		return 0;
	}
	virtual size_t jac_nnz(void) const {
		return 3*2*dim_of_f();
	}
protected:
	const function::idx_type beg_, end_;
};
