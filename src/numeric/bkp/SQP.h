#ifndef HJ_SQP_H_
#define HJ_SQP_H_

#include "func_opt.h"
#include <zjucad/linear_solver/linear_solver.h>

namespace hj_func_opt {

class SQP : public opt
{
public:
  SQP();
  virtual int set_f(function &f);
  virtual int set_c(function *c);
  virtual int solve(double *x, boost::property_tree::ptree &pt);
protected:
  int solve_by_xx_decomposition(double *x, boost::property_tree::ptree &pt);
  function *f_, *c_;
	std::auto_ptr<linear_solver> slv_;
	hj::sparse::csc<double, int32_t> H_, corrected_H_;
	zjucad::matrix::matrix<double> g_, s_, D_, tmp_;
};

}

#endif
