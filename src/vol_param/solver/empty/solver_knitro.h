#ifndef SOLVER_KNITRO_H
#define SOLVER_KNITRO_H

#include "../descriptor/descriptor_base.h"
#include "solver_base.h"

class solver_knitro : public solver_base
{
public:
  solver_knitro(){}
  virtual ~solver_knitro(){}

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


#endif // SOLVER_KNITRO_H
