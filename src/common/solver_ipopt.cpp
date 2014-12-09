#include <Ipopt/coin/IpIpoptApplication.hpp>
#include <Ipopt/coin/IpTNLP.hpp>

#include <zjucad/matrix/io.h>
#include <zjucad/matrix/itr_matrix.h>
#include <hjlib/sparse/sparse.h>
#include <jtflib/function/func_aux.h>
#include <jtflib/function/operation.h>

#include <boost/property_tree/ptree.hpp>

using namespace std;
using namespace Ipopt;
using namespace zjucad::matrix;
using boost::property_tree::ptree;

class solver_use_ipopt : public Ipopt::TNLP
{
public:
  typedef jtf::function::functionN1_t<Ipopt::Number,int32_t> jtf_func_ipopt;
  typedef hj::function::function_t<Ipopt::Number,int32_t> hj_func_ipopt;

  explicit solver_use_ipopt(
      jtf_func_ipopt * func,
      const std::vector<std::shared_ptr<jtf_func_ipopt> > &eqn_func,
      zjucad::matrix::matrix<Ipopt::Number> &init_val)
    :func_(func), eqn_func_(eqn_func), init_val_(init_val){}
  virtual ~solver_use_ipopt(){}

public:
  virtual bool get_nlp_info(
      Ipopt::Index& n, Ipopt::Index& m, Ipopt::Index& nnz_jac_g,
      Ipopt::Index& nnz_h_lag, IndexStyleEnum& index_style);

  //! @brief Method to return the bounds for my problem
  virtual bool get_bounds_info(
      Ipopt::Index n, Ipopt::Number* x_l, Ipopt::Number* x_u,
      Ipopt::Index m, Ipopt::Number* g_l, Ipopt::Number* g_u);

  //! @brief Method to return the starting point for the algorithm
  virtual bool get_starting_point(
      Ipopt::Index n, bool init_x, Ipopt::Number* x,
      bool init_z, Ipopt::Number* z_L, Ipopt::Number* z_U,
      Ipopt::Index m, bool init_lambda,
      Ipopt::Number* lambda);

  //! @brief Method to return the objective value
  virtual bool eval_f(Ipopt::Index n, const Ipopt::Number* x, bool new_x,
                      Ipopt::Number& obj_value);

  //! @brief Method to return the gradient of the objective
  virtual bool eval_grad_f(Ipopt::Index n, const Ipopt::Number* x, bool new_x,
                           Ipopt::Number* grad_f);

  //! @brief Method to return the constraint residuals
  virtual bool eval_g(Ipopt::Index n, const Ipopt::Number* x, bool new_x,
                      Ipopt::Index m, Ipopt::Number* g);

  //! @brief The structure of the jacobian (if "values" is NULL)
  //!        The values of the jacobian (if "values" is not NULL)
  virtual bool eval_jac_g(
      Ipopt::Index n, const Ipopt::Number* x, bool new_x,
      Ipopt::Index m, Ipopt::Index nele_jac, Ipopt::Index* iRow,
      Ipopt::Index *jCol, Ipopt::Number* values);

  //! @brief The structure of the hessian of the lagrangian (if "values" is NULL)
  //!        The values of the hessian of the lagrangian (if "values" is not NULL)
  virtual bool eval_h(Ipopt::Index n, const Ipopt::Number* x, bool new_x,
                      Ipopt::Number obj_factor, Ipopt::Index m,
                      const Ipopt::Number* lambda, bool new_lambda,
                      Ipopt::Index nele_hess, Ipopt::Index* iRow,
                      Ipopt::Index* jCol, Ipopt::Number* values);

  //! @brief  This method is called when the algorithm is complete
  //!         so the TNLP can store/write the solution */
  virtual void finalize_solution(
      Ipopt::SolverReturn status,
      Ipopt::Index n, const Ipopt::Number* x,
      const Ipopt::Number* z_L, const Ipopt::Number* z_U,
      Ipopt::Index m, const Ipopt::Number* g,
      const Ipopt::Number* lambda,
      Ipopt::Number obj_value,
      const Ipopt::IpoptData* ip_data,
      Ipopt::IpoptCalculatedQuantities* ip_cq);

private:
  jtf::function::functionN1_t<Ipopt::Number,int32_t> * func_;
  std::vector<std::shared_ptr<jtf_func_ipopt> > eqn_func_;
  zjucad::matrix::matrix<Ipopt::Number> &init_val_;

  hj::sparse::csc<Number, int32_t> HT_;
  vector<hj::sparse::csc<Number, int32_t> > eqn_HT_;
  map<pair<size_t,size_t>, vector<pair<size_t,double*> > > ij2double_ptr_hes_;

  map<pair<size_t,size_t>, vector<double*> >  ij2double_ptr_g_jac_;
  vector<size_t> g_jac_nnz_;
private:
  solver_use_ipopt(const solver_use_ipopt&);

  /** Overloaded Equals Operator */
  void operator=(const solver_use_ipopt&);
};


bool solver_use_ipopt::get_nlp_info(
    Index& n, Index& m, Index& nnz_jac_g,
    Index& nnz_h_lag, IndexStyleEnum& index_style)
{
  // lack nnz_h_lag
  n = func_->dim();
  // contains equality constraint and inequality constraint
  m = eqn_func_.size();
  vector<double> x(n);

  nnz_jac_g = 0;
  size_t format, nnz;
  g_jac_nnz_.resize(eqn_func_.size());
  for(size_t eqi = 0; eqi < eqn_func_.size(); ++eqi){
      const_cast<jtf_func_ipopt*>(eqn_func_[eqi].get())->gra(&x[0], nnz, 0, 0);
      g_jac_nnz_[eqi] = nnz;
      nnz_jac_g += nnz;
    }

  {
    if(HT_.size(1) == 0 || HT_.size(2) == 0){
        const_cast<jtf_func_ipopt*>(func_)->hes(&x[0], nnz,format,0,0,0);
        if(nnz > 0){
            HT_.resize(func_->dim(), func_->dim(), nnz);
            const_cast<jtf_func_ipopt*>(func_)->hes(&x[0], nnz, format, 0, &HT_.ptr()[0], &HT_.idx()[0]);
          }
      }

    if(eqn_HT_.empty()){
        eqn_HT_.resize(eqn_func_.size());
        for(size_t eqi = 0; eqi < eqn_func_.size(); ++eqi){
            const_cast<jtf_func_ipopt*>(eqn_func_[eqi].get())->hes(&x[0], nnz, format, 0, 0, 0);
            if(nnz > 0){
                eqn_HT_[eqi].resize(n,n, nnz);
                const_cast<jtf_func_ipopt*>(eqn_func_[eqi].get())->hes(&x[0],nnz,format, 0, &eqn_HT_[eqi].ptr()[0],
                    &eqn_HT_[eqi].idx()[0]);
              }
          }
      }

    if(ij2double_ptr_hes_.empty()){
        if(HT_.size(1) != 0 && HT_.size(2) != 0){
            for(size_t i = 0; i < HT_.ptr().size() - 1; ++i){
                for(size_t pi = HT_.ptr()[i]; pi < HT_.ptr()[i+1]; ++pi){
                    ij2double_ptr_hes_[make_pair(i,HT_.idx()[pi])].push_back(
                          make_pair(-1,&HT_.val()[pi]));
                  }
              }
          }
        for(size_t eqi = 0; eqi < eqn_HT_.size(); ++eqi){
            if(eqn_HT_[eqi].size(1) == 0 || eqn_HT_[eqi].size(2) == 0) continue;
            for(size_t i = 0; i < eqn_HT_[eqi].ptr().size() - 1; ++i){
                for(size_t pi = eqn_HT_[eqi].ptr()[i]; pi < eqn_HT_[eqi].ptr()[i+1]; ++pi){
                    ij2double_ptr_hes_[make_pair(i,eqn_HT_[eqi].idx()[pi])].push_back(
                          make_pair(eqi,&eqn_HT_[eqi].val()[pi]));
                  }
              }
          }
      }
  }

  nnz_h_lag = ij2double_ptr_hes_.size();

  index_style = TNLP::C_STYLE;

  return true;
}

bool solver_use_ipopt::get_bounds_info(
    Index n, Number* x_l, Number* x_u,
    Index m, Number* g_l, Number* g_u)
{
  //check the value of n and m
  assert(n == func_->dim());
  assert(m == eqn_func_.size());

  // The default value of nlp_upper_bound_inf and nlp_lower_bound_inf
  // is 1e19 and can be changed through ipopt options.

  // the variables have none lower or upper bounds
  for (Index i = 0; i < n; ++i) {
      x_l[i] = -1e19;
      x_u[i] = 1e19;
    }

  if (!eqn_func_.empty()) {
      // equality constraint
      for (Index i = 0; i < eqn_func_.size(); ++i)
        g_l[i] = g_u[i] = 0;
    }

  return true;
}

//! @brief Method to return the starting point for the algorithm
bool solver_use_ipopt::get_starting_point(
    Index n, bool init_x, Number* x,
    bool init_z, Number* z_L, Number* z_U,
    Index m, bool init_lambda,
    Number* lambda)
{
  assert(init_x == true);
  assert(init_z == false);
  assert(init_lambda == false);

  std::copy(init_val_.begin(), init_val_.end(), x);

  return true;
}

//! @brief Method to return the objective value
bool solver_use_ipopt::eval_f(
    Index n, const Number* x, bool new_x,
    Number& obj_value)
{
  assert(n == func_->dim());
  obj_value = 0;
  func_->val(&x[0], obj_value);

  return true;
}

//! @brief Method to return the gradient of the objective
bool solver_use_ipopt::eval_grad_f(
    Index n, const Number* x, bool new_x, Number* grad_f)
{
  assert(n == func_->dim());
  itr_matrix<Number*> grad_f_m(func_->dim(),1,grad_f);
  grad_f_m *= 0;
  const_cast<jtf_func_ipopt*>(func_)->gra(x, grad_f);

  return true;
}

//! @brief Method to return the constraint residuals
bool solver_use_ipopt::eval_g(
    Index n, const Number* x, bool new_x,
    Index m, Number* g)
{
  if(eqn_func_.empty())
    return true;

  assert(n == func_->dim());
  assert(m == eqn_func_.size());

  itr_matrix<Number*> gm(m,1,g);
  gm *= 0;
  for(size_t eqi = 0; eqi < eqn_func_.size(); ++eqi){
      eqn_func_[eqi]->val(x, g[eqi]);
    }

  return true;
}

//! @brief The structure of the jacobian (if "values" is NULL)
//!        The values of the jacobian (if "values" is not NULL)
bool solver_use_ipopt::eval_jac_g(
    Index n, const Number* x, bool new_x,
    Index m, Index nele_jac, Index* iRow,
    Index *jCol, Number* values)
{
  static vector<Number> fake_x(n);
  vector<Number> g;
  vector<Index> idx;
  size_t count = 0;

  for(size_t eqi = 0; eqi < eqn_func_.size(); ++eqi){
      g.resize(g_jac_nnz_[eqi]);
      idx.resize(g_jac_nnz_[eqi]);

      if(values == NULL)
        const_cast<jtf_func_ipopt*>(eqn_func_[eqi].get())->gra(&fake_x[0], g_jac_nnz_[eqi],&g[0], &idx[0]);
      else
        const_cast<jtf_func_ipopt*>(eqn_func_[eqi].get())->gra(x, g_jac_nnz_[eqi], &g[0], &idx[0]);

      if(values == NULL) {
          for(size_t i = 0; i < idx.size(); ++i, ++count){
              iRow[count] = eqi;
              jCol[count] = idx[i];
            }
        }else{
          for(size_t i = 0; i < idx.size(); ++i, ++count){
              values[count] = g[i];
            }
        }
    }

  return true;
}

//! @brief The structure of the hessian of the lagrangian (if "values" is NULL)
//!        The values of the hessian of the lagrangian (if "values" is not NULL)
bool solver_use_ipopt::eval_h(
    Index n, const Number* x, bool new_x,
    Number obj_factor, Index m,
    const Number* lambda, bool new_lambda,
    Index nele_hess, Index* iRow,
    Index* jCol, Number* values)
{  
  size_t format;

  if(values != NULL){
      if(hj::sparse::nnz(HT_)){
          HT_.val() *= 0;
          size_t nnz = hj::sparse::nnz(HT_);
          const_cast<jtf_func_ipopt*>(func_)->hes(x, nnz, format, &HT_.val()[0], &HT_.ptr()[0], &HT_.idx()[0]);
          HT_.val() *= obj_factor;
        }

      for(size_t eqi = 0; eqi < eqn_func_.size(); ++eqi){
          if(eqn_HT_[eqi].size(1) != 0 && eqn_HT_[eqi].size(2) != 0){
              eqn_HT_[eqi].val() *= 0;
              size_t nnz = hj::sparse::nnz(eqn_HT_[eqi]);
              const_cast<jtf_func_ipopt*>(eqn_func_[eqi].get())->hes(
                    x, nnz , format,
                    &eqn_HT_[eqi].val()[0], &eqn_HT_[eqi].ptr()[0],
                  &eqn_HT_[eqi].idx()[0]);
            }
        }
    }


  if(!ij2double_ptr_hes_.empty()){
      size_t term_idx = 0;
      for(const auto & term : ij2double_ptr_hes_){
          if(values == NULL){
              iRow[term_idx] = term.first.first;
              jCol[term_idx] = term.first.second;
            }else{
              values[term_idx] = 0;
              for(const auto & ptr : term.second){
                  double temp = *ptr.second;
                  if(ptr.first != -1)// if is constraint
                    temp *= lambda[ptr.first];
                  values[term_idx] += temp;
                }
            }
          ++term_idx;
        }
    }

  return true;
}

//! @brief  This method is called when the algorithm is complete
//!         so the TNLP can store/write the solution */
void solver_use_ipopt::finalize_solution(
    SolverReturn status,
    Index n, const Number* x,
    const Number* z_L, const Number* z_U,
    Index m, const Number* g,
    const Number* lambda,
    Number obj_value,
    const IpoptData* ip_data,
    IpoptCalculatedQuantities* ip_cq)
{
  assert(n == init_val_.size());
  std::copy(x, x+ n, init_val_.begin());
}

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

int ipopt_solve(
    matrix<double> &init_node,
    jtf::function::functionN1_t<double,int32_t> &func,
    shared_ptr<vector<shared_ptr<jtf::function::functionN1_t<double,int32_t> > > > eqn_cons,
    boost::property_tree::ptree &pt)
{
  const size_t max_iter_num = pt.get<size_t>("iter.value");
  Ipopt::SmartPtr<Ipopt::IpoptApplication> app = IpoptApplicationFactory();
  app->Options()->SetNumericValue("tol", 1e-9);
  //app->Options()->SetStringValue("derivative_test", "second-order");
  app->Options()->SetStringValue("output_file", "ipopt.out");
  //app->Options()->SetStringValue("mu_strategy", "adaptive");
  app->Options()->SetIntegerValue("max_iter", max_iter_num);
  app->Options()->SetStringValue("linear_solver", "ma97");
  app->Options()->SetNumericValue("acceptable_obj_change_tol", 1e-4);
  //app->Options()->SetStringValue("hessian_approximation","limited-memory");
  //app->Options()->SetStringValue("dependency_detector", "mumps");

  Ipopt::ApplicationReturnStatus status;
  status = app->Initialize();
  if (status != Ipopt::Solve_Succeeded) {
      std::cout << "# [error] Error during initialization!" << std::endl;
      return __LINE__;
    }

  vector<std::shared_ptr<solver_use_ipopt::jtf_func_ipopt> > eqn_vec;
  {
    if(eqn_cons != nullptr){
        for(size_t i = 0; i < eqn_cons->size(); ++i)
          eqn_vec.push_back((*eqn_cons)[i]);
      }
  }

  Ipopt::SmartPtr<Ipopt::TNLP> mynlp(
        new solver_use_ipopt(&func,eqn_vec, init_node));

  status = app->OptimizeTNLP(mynlp);

  if (status == Ipopt::Solve_Succeeded)
    std::cout << "# [info]  The problem solved!" << std::endl;
  else
    std::cout << "# [error] The problem failed!" << std::endl;

  return 0;
}
