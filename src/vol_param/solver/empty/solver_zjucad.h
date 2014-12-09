#ifndef SOLVER_ZJUCAD_H
#define SOLVER_ZJUCAD_H

#include "solver_base.h"

class solver_zjucad : public solver_base
{
public:
  virtual ~solver_zjucad(){}

public:

  ///
  /// @brief optimization, solve with zjucad solver
  /// @param init_node init_node
  /// @param desc descriptor of problem
  /// @param pt   input configuration
  ///
  virtual void solve(zjucad::matrix::matrix<double> & init_node,
                     const std::shared_ptr<descriptor_base> & desc,
                     boost::property_tree::ptree & pt) const;
};

#endif
