#ifndef LSCM_IPOPT_FOLDFREE_H
#define LSCM_IPOPT_FOLDFREE_H

#include <vector>
#include <map>
#include <memory>
#include <Ipopt/coin/IpTNLP.hpp>

#include "tri_mesh.h"

namespace jy{
typedef std::vector<std::map<size_t,double> >      rowmat;
typedef std::map<size_t, double>::const_iterator   row_citerd;
typedef std::vector<size_t>      vectori;
typedef std::vector<double>      vectord;

using Ipopt::Index;
using Ipopt::Number;
using Ipopt::SolverReturn;

class lscm_ipopt_foldfree:public Ipopt::TNLP
{
public:
  lscm_ipopt_foldfree(const char *file_mesh,
                      const char *file_fix_vert,
                      const char *output_uv);
  virtual ~lscm_ipopt_foldfree();
public:
  /**@name Overloaded from TNLP */
  //@{
  /* Method to return some info about the nlp */
  virtual bool get_nlp_info(Index& n,Index& m, Index& nnz_jac_g,
                            Index& nnz_h_lag, IndexStyleEnum& index_style);

  /* Method to return the bounds for my problem */
  virtual bool get_bounds_info(Index n, Number* x_l, Number* x_u,
                               Index m, Number* g_l, Number* g_u);

  /* Method to return the starting point for the algorithm */
  virtual bool get_starting_point(Index n, bool init_x, Number* x,
                                  bool init_z, Number* z_L, Number* z_U,
                                  Index m, bool init_lambda,
                                  Number* lambda);

  /* Method to return the objective value */
  virtual bool eval_f(Index n, const Number* x, bool new_x, Number& obj_value);

  /* Method to return the gradient of the objective */
  virtual bool eval_grad_f(Index n, const Number* x, bool new_x, Number* grad_f);

  /* Method to return the constraint residuals */
  virtual bool eval_g(Index n, const Number* x, bool new_x, Index m, Number* g);

  /* Method to return:
      *   1) The structure of the jacobian (if "values" is NULL)
      *   2) The values of the jacobian (if "values" is not NULL)
      */
  virtual bool eval_jac_g(Index n, const Number* x, bool new_x,
                          Index m, Index nele_jac, Index* iRow, Index *jCol,
                          Number* values);

  /* Method to return:
      *   1) The structure of the hessian of the lagrangian (if "values" is NULL)
      *   2) The values of the hessian of the lagrangian (if "values" is not NULL)
      */
  virtual bool eval_h(Index n, const Number* x, bool new_x,
                      Number obj_factor, Index m, const Number* lambda,
                      bool new_lambda, Index nele_hess, Index* iRow,
                      Index* jCol, Number* values);

  //@}

  /* @name Solution Methods */
  //@{
  /* This method is called when the algorithm is complete so the TNLP can store/write the solution */
  virtual void finalize_solution(SolverReturn status,
                                 Index n, const Number* x, const Number* z_L, const Number* z_U,
                                 Index m, const Number* g, const Number* lambda,
                                 Number obj_value, const Ipopt::IpoptData* ip_data,Ipopt::IpoptCalculatedQuantities* ip_cq);
  //@}
public:
  bool init(const char *file_mesh,
            const char *file_fix_vert);
  //virtual bool init_guess(matrixd &x);

  void update_hessian(const double *lambda,bool add);
  size_t check_fold_over(const double *texture);

  bool write_uv(const char *file_name, const double *uv);
  void write_texture_mesh(const char *filename,const double *texture_uv);

public:
protected:
  bool       reverse_;
  std::unique_ptr<tri_mesh> mesh_;
  rowmat     hess_;
  vectord    int_weight_;
  vectord    lscm_weight_;

  vectori    fix_ids_;
  vectord    fix_coords_;

  vectori    front_map_;
  vectori    back_map_;

  vectord    init_x_;
  const std::string output_uv_;
  double avg_area_;
public:
  //vectord    texture_;
private:
  /**@name Methods to block default compiler methods.
      * The compiler automatically generates the following three methods.
      *  Since the default compiler implementation is generally not what
      *  you want (for all but the most simple classes), we usually
      *  put the declarations of these methods in the private section
      *  and never implement them. This prevents the compiler from
      *  implementing an incorrect "default" behavior without us
      *  knowing. (See Scott Meyers book, "Effective C++")
      *
      */
  //@{
  lscm_ipopt_foldfree(const lscm_ipopt_foldfree&);
  lscm_ipopt_foldfree& operator=(const lscm_ipopt_foldfree&);


};

bool ipopt_foldfree_parameterization(
    const char *file_mesh,
    const char *file_fix_vert,
    const char * output_uv);


}

#endif // LSCM_IPOPT_FOLDFREE_H
