#include <boost/algorithm/string.hpp>
#include <sstream>
#include <iostream>

#include <jtflib/function/operation.h>
#include <zjucad/matrix/matrix.h>
#include <jtflib/optimizer/optimizer.h>

#include "../descriptor/descriptor_base.h"
#include "solver_zjucad.h"
#include "util.h"

using namespace std;
using namespace zjucad::matrix;
using boost::property_tree::ptree;

void solver_zjucad::solve(
    zjucad::matrix::matrix<double> & init_node,
    const std::shared_ptr<descriptor_base> & desc,
    boost::property_tree::ptree & pt) const
{
//  std::shared_ptr<vector<hj_func_cons_ptr> >
//      all_funcs(new vector<hj_func_cons_ptr> );

//  {  // assemble functions
//    all_funcs->push_back(desc->get_objective());

//    if(desc->get_constraint().size() != 0){

//        matrix<double> cons_weight = ones<double>(desc->get_constraint().size(),1);

//        load_cons_weight(pt, desc->get_constraint_type(), cons_weight);

//        double w = 0;
//        for(size_t i = 0; i < cons_weight.size(); ++i){
//            if(cons_weight[i] < 0) {
//                cerr << "# [error] sorry, do not support hard constraints." << endl;
//                continue; // Here I record weight < 0 as hard constraints
//              }
//            all_funcs->push_back(
//                  std::shared_ptr<const jtf_func>(desc->get_constraint()[i] * w));
//          }
//      }
//  }

//  hj_func_cons_ptr func(
//        hj::function::new_catenated_function<double,int32_t>(all_funcs));

//  matrix<double> resid(func->dim_of_f(),1);

// zjucad::optimize(*func, init_node, resid, pt);
}
