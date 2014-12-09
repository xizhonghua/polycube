#ifndef SOLVER_NLOPT_H
#define SOLVER_NLOPT_H

#include "../descriptor/descriptor_base.h"
#include "solver_base.h"

class solver_nlopt : public solver_base
{
public:
  solver_nlopt(){}
  virtual ~solver_nlopt(){}
public:
  virtual void solve(zjucad::matrix::matrix<double> &init_node,
                     const std::shared_ptr<descriptor_base> & desc,
                     boost::property_tree::ptree & pt) const;
};

#endif // NLOPT_H
