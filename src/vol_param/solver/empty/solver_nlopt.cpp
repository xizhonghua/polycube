#include "solver_nlopt.h"

#include <nlopt.hpp>
#include <jtflib/function/func_aux.h>
#include "../common/util.h"

using namespace std;
using namespace zjucad::matrix;


double set_f(const vector<double> &x, vector<double> &gradient, void *func_data)
{
  const jtf_func * f = reinterpret_cast<const jtf_func*>(func_data);

  if(gradient.size()){
      f->gra(&x[0], &gradient[0]);
    }
  double v = 0;
  f->val(&x[0], v);
  return v;
}

void solver_nlopt::solve(matrix<double> &init_node,
			 const shared_ptr<descriptor_base> & desc,
			 boost::property_tree::ptree & pt)const
{
  jtf_func_cons_ptr obj_f = desc->get_objective();
  const vector<jtf_func_cons_ptr> & eqn_cons = desc->get_eqn_constraint();
  const vector<jtf_func_cons_ptr> &inequal_cons = desc->get_ineqn_constraint();

  const string algorithm = pt.get<string>("alg.value");
  pt.put("alg.desc", "auglag_eq");
  nlopt::algorithm am;
  switch(str2int(algorithm.c_str())){
    case str2int("auglag_eq"): am = nlopt::AUGLAG_EQ;break;
    default:
      throw std::invalid_argument("invalid alg type.");
    }

  nlopt::opt opt(nlopt::LD_MMA, init_node.size());
  opt.set_min_objective(set_f, obj_f.get());
  for(size_t i = 0; i < eqn_cons.size(); ++i)
    opt.add_equality_constraint(set_f, eqn_cons[i].get(), 1e-8);

  vector<double> init_node_vec(init_node.size());
  std::copy(init_node.begin(), init_node.end(), init_node_vec.begin());
  double minf;
  nlopt::result result = opt.optimize(init_node_vec, minf);
  std::copy(init_node_vec.begin(), init_node_vec.end(), init_node.begin());
}
