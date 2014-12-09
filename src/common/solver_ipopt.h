#ifndef COMMON_SOLVER_IPOPT_H
#define COMMON_SOLVER_IPOPT_H

#include <zjucad/matrix/matrix.h>
#include <jtflib/function/function.h>
#include <zjucad/ptree/ptree.h>

int ipopt_solve(
    zjucad::matrix::matrix<double> &init_node,
    jtf::function::functionN1_t<double,int32_t> &func,
    std::shared_ptr<std::vector<std::shared_ptr<jtf::function::functionN1_t<double,int32_t> > > > eqn_cons,
    boost::property_tree::ptree &pt);

#endif // SOLVER_IPOPT_H
