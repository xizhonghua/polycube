#ifndef SOLVER_UTIL_H
#define SOLVER_UTIL_H

#include <zjucad/matrix/matrix.h>
#include <zjucad/ptree/ptree.h>
#include <hjlib/sparse/sparse.h>
#include <zjucad/optimizer/optimizer.h>
#include <iostream>

#include <jtflib/optimizer/optimizer.h>
#include "../descriptor/def.h"
#include "../descriptor/descriptor_vol.h"
#include "../descriptor/func_terms/linear_equation.h"
#include "../descriptor/func_terms/vol-anti-flip.h"
//#include "../descriptor/func_terms/arap.h"

///
/// @brief load_cons_weight, -1: hard constraint, others: soft with weight
///  assume cons_weight is allocated
/// @param pt
/// @param cons_type
/// @param cons_weight
///
void load_cons_weight(
    boost::property_tree::ptree &pt,
    const std::vector<std::tuple<std::string,size_t,size_t> > & cons_type,
    zjucad::matrix::matrix<double> &cons_weight);


///
/// \brief get_init_node: min |x - NMT * u|^2
/// \param u
/// \param NMT
/// \param x
///
template <typename T1>
void get_init_node(
    zjucad::matrix::matrix<T1> & u,
    std::shared_ptr<descriptor_vol> dv,
    const node_mapping & NMT,
    const zjucad::matrix::matrix<T1> & x)
{
  using namespace std;
  using namespace zjucad::matrix;
  boost::property_tree::ptree pt;

  pt.put("package.value","hj");
  pt.put("iter.value",10);

  matrix<double> x_temp = x;
  if(x.size() != NMT.ZT.size(2)){
      x_temp.resize(x.size(1),NMT.ZT.size(2)/x.size(1));
      std::copy(x.begin(), x.end(), x_temp.begin());
    }

  dv->recover_gap_node(x_temp);

  cerr << "fit u to x" << endl;
  hj_func_cons_ptr fit_u2x(new AX_b_func(NMT, x_temp));
  matrix<double> residuial(fit_u2x->dim_of_f(),1);
  zjucad::optimize(*fit_u2x, u, residuial, pt);


//  SIMPLEX sim_ = TET;
//  matrix<double> b(x_temp.size(),1);
//  vector<hj_func_cons_ptr> funcs;
//  double w = 0.1;
//  hj_func_cons_ptr all_funcs;
//  for(size_t i = 0; i < 5; ++i){
//      w *= 2;
//      dv->recover_gap_node(x_temp);

//      cerr << "fit u to x" << endl;
//      hj_func_cons_ptr fit_u2x(new AX_b_func(NMT, x_temp));
//      matrix<double> residuial(fit_u2x->dim_of_f(),1);
//      zjucad::optimize(*fit_u2x, u, residuial, pt);

//      // fit x to u
//      b *= 0;
//      b = dv->bi_.NM.q;
//      funcs.clear();
//      hj::sparse::mv(true, dv->bi_.NM.ZT, u, b);
//      funcs.push_back(
//            hj_func_cons_ptr(
//              build_arap_func(dv->bi_.cut_tet, x_temp, sim_, 0)));

//      funcs.push_back(hj_func_cons_ptr(new X_c_func(b,w)));
//      all_funcs.reset(hj::function::new_catenated_function<double,int32_t>(funcs));
//      matrix<double> residual_x(all_funcs->dim_of_f(),1);
//      zjucad::optimize(*all_funcs,x_temp, residual_x, pt);
//    }

  cerr << "# after initialize." << endl;
}

#endif // SOLVER_UTIL_H
