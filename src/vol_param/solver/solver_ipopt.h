#ifndef SOLVER_IPOPT_H
#define SOLVER_IPOPT_H

#include "../descriptor/descriptor_base.h"
#include "solver_base.h"

class solver_ipopt : public solver_base
{
public:
  solver_ipopt(){}
  virtual ~solver_ipopt(){}

  virtual std::string get_solver_type()const {return "ipopt";}
public:
  ///
  /// @brief optimization, solve with zjucad solver
  /// @param init_node init_node
  /// @param desc descriptor of problem
  /// @param pt   input configuration
  ///
  virtual int solve(zjucad::matrix::matrix<double> & init_node,
                     const std::shared_ptr<descriptor_base> & desc,
                     boost::property_tree::ptree & pt) const;
};

#endif // SOLVER_IPOPT_H
