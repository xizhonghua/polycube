#include "solver_jtf.h"
#include "util.h"
#include "../descriptor/descriptor_vol.h"
#include <jtflib/optimizer/optimizer.h>
#include <hjlib/sparse/sparse.h>
#include <hjlib/sparse/format.h>
#include <jtflib/function/func_aux.h>
#include <zjucad/matrix/io.h>
#include <fstream>

using namespace std;
using namespace zjucad::matrix;

int solver_jtf::solve(matrix<double> &init_node,
		      const shared_ptr<descriptor_base> & desc,
		      boost::property_tree::ptree & pt)
{
  jtf_func_ptr obj_f = desc->get_objective();
  const vector<jtf_func_ptr> & eqn_cons = desc->get_eqn_constraint();
  const vector<jtf_func_ptr> &inequal_cons = desc->get_ineqn_constraint();

  shared_ptr<descriptor_vol> dv = dynamic_pointer_cast<descriptor_vol>(desc);

  matrix<double> init_node_handled;
  if(desc->has_node_mapping()){
      const node_mapping & NM = desc->get_node_mapping();
      init_node_handled.resize(NM.ZT.size(1),1);

      matrix<double> init_node_extend(init_node.size(1), NM.ZT.size(2)/3);
      copy(init_node.begin(), init_node.end(), init_node_extend.begin());
      get_init_node(init_node_handled, dv, NM, init_node_extend);
    }else
    init_node_handled = init_node;

  const string solver = pt.get<string>("alg.value");
  int rtn = 0;
  if(solver == "SQP"){
      auto callback_ptr = get_callback();
      rtn = jtf::optimize(*obj_f, init_node_handled, pt, eqn_cons.empty()?nullptr:&eqn_cons,nullptr,callback_ptr);
      if(rtn == -1){
          cerr << "# [info] sunk in infeasible domain, need remesh." << endl;
        }else if(rtn != 0){
          cerr << "# [error] solve fail." << endl;
        }
    }else{
      cerr << "# [error] in this app, I only support SQP solver." << endl;
    }

  if(desc->has_node_mapping()){
      const hj::sparse::csc<double,int32_t> & node_mapping = desc->get_node_mapping().ZT;
      matrix<double> init_node_extend = zeros<double>(init_node.size(1), node_mapping.size(2)/3);
      hj::sparse::mv(true, node_mapping, init_node_handled, init_node_extend);
      init_node_extend(colon()) += desc->get_node_mapping().q;
      itr_matrix<const double *> init_node_m(init_node.size(1), init_node.size(2),&init_node_extend[0]);
      init_node = init_node_m;
    }else
    init_node = init_node_handled;
  return rtn;
}
