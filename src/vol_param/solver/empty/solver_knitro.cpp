#include "solver_knitro.h"
#include "nlpProblemDef.h"
#include <knitro/knitro.h>
#include <zjucad/matrix/itr_matrix.h>

#include <hjlib/sparse/sparse.h>
using namespace std;
using namespace zjucad::matrix;
using boost::property_tree::ptree;

class knitro_wrapper : public NlpProblemDef
{
private:
  knitro_wrapper (const knitro_wrapper &);
  knitro_wrapper & operator =(const knitro_wrapper &);
protected:
  knitro_wrapper (void);
public:
  static knitro_wrapper *  getTheInstance (void);

  static void init_instance(knitro_wrapper * isk,
                            jtf_func_cons_ptr & obj,
                            const std::vector<jtf_func_cons_ptr> &eqn,
                            const zjucad::matrix::matrix<double> & init_x);
  ~knitro_wrapper (void){}

  //++ Declare virtual base class methods that are implemented here.
  //++ See NlpProblemDef.h for descriptions.

  virtual int   getN (void) {return var_number_;}

  virtual int   getM (void) {return eqn_.size();}

  virtual void  getInitialX (double * const  dax) {
    if(!_daXInit){
        cerr << "# [error] getInitialX Must be called after 'loadProblemIntoKnitro' " << endl;
        exit( EXIT_FAILURE );
      }
    std::copy(_daXInit, _daXInit + var_number_, dax);
  }

  virtual bool  loadProblemIntoKnitro (KTR_context_ptr  kc);
  virtual bool  areDerivativesImplemented(
      const DerivativesImplementedType  nWhichDers);

  virtual int evalFC (const double * const  daX,
                      double * const  dObj,
                      double * const  daC,
                      void   *        userParams);

  virtual int evalGA (const double * const  daX,
                      double * const  daG,
                      double * const  daJ,
                      void   * userParams);

  virtual int  evalH (const double * const  daX,
                      const double * const  daLambda,
                      const double          dSigma,
                      double * const  daH,
                      void   *        userParams);
  virtual int  evalHV (const double * const  daX,
                       const double * const  daLambda,
                       const double          dSigma,
                       double * const  daHV,
                       void   *        userParams);

private:
  static knitro_wrapper *  _pTheNlpProblemDefInstance;

  //  double get_equal_equation_val(const int n, const double *x,
  //                                const equation<double> &eq) const;
public:
  zjucad::matrix::matrix<double> init_x_;
  jtf_func_cons_ptr func_;
  size_t var_number_;
  std::vector<jtf_func_cons_ptr> eqn_;
  size_t jac_nnz_;
  size_t h_nnz_;
  hj::sparse::csc<double,int32_t> Hessian_;
};

knitro_wrapper * knitro_wrapper::_pTheNlpProblemDefInstance = NULL;
static NlpProblemDef *  g_pOptProblem = NULL;

//--------------------------------------------------------------------
//  Internal Function:  wrapperEvalFC
//--------------------------------------------------------------------
/** By necessity this wrapper signature matches the function KTR_callback.
 *  It calls the current optimization problem's eval method.
 */
static int  wrapperEvalFC (
    const int evalRequestCode, const int n, const int m, const int nnzJ,
    const int nnzH, const double * const  daX, const double * const  daLambda,
    double * const  dObj, double * const  daC, double * const  daG,
    double * const  daJ, double * const  daH, double * const  daHV,
    void   * userParams)
{
  if (g_pOptProblem == NULL) {
      cout << "*** Problem not defined  <wrapperEvalFC>\n";
      return( -1 );
    }
  if (evalRequestCode != KTR_RC_EVALFC) {
      cout << "*** Bad request code " << evalRequestCode
           << "  <wrapperEvalFC>\n";
      return( -1 );
    }
  return( g_pOptProblem->evalFC (daX, dObj, daC, userParams) );
}


//--------------------------------------------------------------------
//  Internal Function:  wrapperEvalGA
//--------------------------------------------------------------------
/** By necessity this wrapper signature matches the function KTR_callback.
 *  It calls the current optimization problem's eval method.
 */
static int  wrapperEvalGA (
    const int evalRequestCode, const int n, const int m, const int nnzJ,
    const int nnzH, const double * const  daX, const double * const  daLambda,
    double * const  dObj, double * const  daC, double * const  daG,
    double * const  daJ,  double * const  daH, double * const  daHV,
    void   * userParams)
{
  if (g_pOptProblem == NULL){
      cout << "*** Problem not defined  <wrapperEvalGA>\n";
      return( -1 );
    }
  if (evalRequestCode != KTR_RC_EVALGA) {
      cout << "*** Bad request code " << evalRequestCode
           << "  <wrapperEvalGA>\n";
      return( -1 );
    }
  return( g_pOptProblem->evalGA (daX, daG, daJ, userParams) );
}


//--------------------------------------------------------------------
//  Internal Function:  wrapperEvalHorHV
//--------------------------------------------------------------------
/** By necessity this wrapper signature matches the function KTR_callback.
 *  It calls the current optimization problem's eval method.
 */
static int  wrapperEvalHorHV (
    const int evalRequestCode, const int n,  const int m, const int nnzJ,
    const int nnzH,const double * const  daX, const double * const  daLambda,
    double * const  dObj, double * const  daC, double * const  daG,
    double * const  daJ,  double * const  daH, double * const  daHV,
    void   * userParams)
{
  if (g_pOptProblem == NULL) {
      cout << "*** Problem not defined  <wrapperEvalHorHV>\n";
      return( -1 );
    }
  if (evalRequestCode == KTR_RC_EVALH) {
      if (g_pOptProblem->areDerivativesImplemented (nCAN_COMPUTE_H) == false)
        {
          cout << "*** This problem not evaluate H  <wrapperEvalHorHV>\n";
          return( -1 );
        }
      return( g_pOptProblem->evalH (daX, daLambda, 1.0, daH, userParams) );
    }
  else if (evalRequestCode == KTR_RC_EVALH_NO_F) {
      if (g_pOptProblem->areDerivativesImplemented (nCAN_COMPUTE_H) == false) {
          cout << "*** This problem not evaluate H  <wrapperEvalHorHV>\n";
          return( -1 );
        }
      return( g_pOptProblem->evalH (daX, daLambda, 0.0, daH, userParams) );
    } else if (evalRequestCode == KTR_RC_EVALHV) {
      if (g_pOptProblem->areDerivativesImplemented (nCAN_COMPUTE_HV) == false){
          cout << "*** This problem not evaluate H*v  <wrapperEvalHorHV>\n";
          return( -1 );
        }
      return( g_pOptProblem->evalHV (daX, daLambda, 1.0, daHV, userParams) );
    } else if (evalRequestCode == KTR_RC_EVALHV_NO_F) {
      if (g_pOptProblem->areDerivativesImplemented (nCAN_COMPUTE_HV) == false){
          cout << "*** This problem not evaluate H*v  <wrapperEvalHorHV>\n";
          return( -1 );
        }
      return( g_pOptProblem->evalHV (daX, daLambda, 0.0, daHV, userParams) );
    } else {
      cout << "*** Bad request code " << evalRequestCode
           << "  <wrapperEvalHorHV>\n";
      return( -1 );
    }
}

knitro_wrapper::knitro_wrapper(void)
{
  _daXInit = NULL;
}

void knitro_wrapper::init_instance(
    knitro_wrapper * isk,
    jtf_func_cons_ptr & obj,
    const std::vector<jtf_func_cons_ptr> &eqn,
    const zjucad::matrix::matrix<double> & init_x)
{
  isk->var_number_ = obj->dim();
  isk->eqn_ = eqn;
  isk->func_ = obj;
  isk->init_x_ = init_x;
}

knitro_wrapper * knitro_wrapper::getTheInstance (void)
{
  if(_pTheNlpProblemDefInstance == NULL)
    _pTheNlpProblemDefInstance = new knitro_wrapper();
  return _pTheNlpProblemDefInstance;
}

bool knitro_wrapper::loadProblemIntoKnitro(KTR_context_ptr kc)
{
  _nN = var_number_;
  _nM = eqn_.size();

  static vector<double> x_temp(_nN);
  _nNnzJ = 0;
  size_t nnz_i = 0;
  for (size_t e = 0; e < eqn_.size(); ++e) {
      jtf_func_cons_ptr eqi = eqn_[e];
      eqi->gra(&x_temp[0], nnz_i, 0,0);
      _nNnzJ += nnz_i;
    }

  size_t format;
  size_t nNnzH;
  func_->hes(&x_temp[0], nNnzH, format, 0,0,0);
  _nNnzH = nNnzH;

  //---- VARIABLES ARE BOUNDED FROM BELOW.
  _daXLo  = new double[_nN];
  _daXUp  = new double[_nN];

  for (int i = 0; i < _nN; i++){
      _daXLo[i] = -1e20;
      _daXUp[i] = 1e20;
    }

  //---- THE CONSTRAINTS IS A LINEAR INEQUALITY.
  //---- PUT THE CONSTANT TERM IN THE RIGHT-HAND SIDE.
  _naCType  = new int[_nM];
  _daCLo    = new double[_nM];
  _daCUp    = new double[_nM];

  // equality constraint
  const size_t num_equation = eqn_.size();
  for (size_t i = 0; i < num_equation; ++i) {
      _daCLo[i] = _daCUp[i] = 0;
      _naCType[i] = KTR_CONTYPE_LINEAR; //WARNING!!! knitro solver will use null space project to speed up linear constraints
    }


  //---- SPECIFY THE CONSTRAINT JACOBIAN SPARSITY STRUCTURE.
  _naJacIndexVars = new int[_nNnzJ];
  _naJacIndexCons = new int[_nNnzJ];

  size_t count = 0;
  vector<int32_t> eqi_idx;
  vector<double> eqi_val;
  for (size_t e = 0; e < eqn_.size(); ++e) {
      jtf_func_cons_ptr eqi = eqn_[e];

      eqi->gra(&x_temp[0], nnz_i,0,0);
      eqi_idx.resize(nnz_i);
      eqi_val.resize(nnz_i);
      eqi->gra(&x_temp[0], nnz_i, &eqi_val[0], &eqi_idx[0]);

      for (size_t ii = 0; ii < nnz_i; ++ii) {
          _naJacIndexCons[count] = e;
          _naJacIndexVars[count++] = eqi_idx[ii];
        }
    }

  //---- SPECIFY THE HESSIAN OF THE LAGRANGIAN SPARSITY STRUCTURE.
  _naHessRows = new int[_nNnzH];
  _naHessCols = new int[_nNnzH];

  Hessian_.resize(_nN,_nN, _nNnzH);

  func_->hes(&x_temp[0], nNnzH, format, 0, &Hessian_.ptr_[0], &Hessian_.idx_[0]);

  for(size_t nni = 0, idx = 0; nni < var_number_; ++nni){
      for(size_t ptri = Hessian_.ptr_[nni]; ptri < Hessian_.ptr_[nni+1]; ++ptri){
          _naHessRows[idx] = nni;
          _naHessCols[idx++] = Hessian_.idx_[ptri];
        }
    }

  //---- INITIAL GUESS FOR x AND lambda.
  _daXInit = new double[_nN];
  double *  daLambdaInit = new double[_nM + _nN];
  for (int i = 0; i < _nN; i++)
    _daXInit[i] = init_x_[i];
  for (int i = 0; i < _nM + _nN; i++)
    daLambdaInit[i] = 100.0;

  int rtn = KTR_init_problem (kc, _nN,
                              KTR_OBJGOAL_MINIMIZE, KTR_OBJTYPE_GENERAL,
                              _daXLo, _daXUp,
                              _nM, _naCType, _daCLo, _daCUp,
                              _nNnzJ, _naJacIndexVars, _naJacIndexCons,
                              _nNnzH, _naHessRows, _naHessCols,
                              _daXInit, daLambdaInit);

  delete [] _daXLo;
  delete [] _daXUp;
  delete [] _naCType;
  delete [] _daCLo;
  delete [] _daCUp;
  delete [] _naJacIndexVars;
  delete [] _naJacIndexCons;
  delete [] _naHessRows;
  delete [] _naHessCols;
  delete [] daLambdaInit;

  if (rtn != 0){
      cout << "*** KTR_init_problem() returned " << rtn << "\n";
      return( false );
    }

  return( true );
}


int knitro_wrapper::evalFC(const double *const daX,
                           double *const dObj,
                           double *const daC,
                           void *userParams)
{
  itr_matrix<const double *> daX0(var_number_,1, daX);

  func_->val(daX, *dObj);

  for(size_t ci = 0; ci < eqn_.size(); ++ci){
      eqn_[ci]->val(daX, daC[ci]);
    }

  return 0;
}


int knitro_wrapper::evalGA (const double * const  daX,
                            double * const  daG,
                            double * const  daJ,
                            void   * userParams)
{
  func_->gra(daX, daG);

  size_t count = 0;
  size_t nnz = 0;
  vector<double> g;
  vector<int32_t> idx;
  for(size_t ei = 0; ei < eqn_.size(); ++ei){
      jtf_func_cons_ptr eqi = eqn_[ei];
      eqi->gra(daX, nnz, 0,0);
      g.resize(nnz);
      idx.resize(nnz);
      eqi->gra(daX, nnz, &g[0], &idx[0]);
      for(size_t i = 0; i < nnz; ++i){
          daJ[count++] = idx[i];
        }
    }

  return 0;
}

int knitro_wrapper::evalH (const double * const  daX,
                           const double * const  daLambda,
                           const double          dSigma,
                           double * const  daH,
                           void   *        userParams)
{
  return -1;

  size_t nnz_H = _nNnzH, foramt;
  func_->hes(daX, nnz_H, foramt, &Hessian_.val_[0], &Hessian_.ptr_[0], &Hessian_.idx_[0]);

  Hessian_.val_ *= dSigma;
  return 0;
}

int knitro_wrapper::evalHV (const double * const  daX,
                     const double * const  daLambda,
                     const double          dSigma,
                     double * const  daHV,
                     void   *        userParams)
{
  return -1;
}

bool knitro_wrapper::areDerivativesImplemented(
    const DerivativesImplementedType nWhichDers)
{
  if (nWhichDers == nCAN_COMPUTE_GA)
    return( true );
  if (nWhichDers == nCAN_COMPUTE_H)
    return( true );
  if (nWhichDers == nCAN_COMPUTE_HV)
    return( true );
  return( false );
}

///////////////////////////////////////////////////////////////////////////////

void solver_knitro::solve(
    matrix<double> &init_node,
    const std::shared_ptr<descriptor_base> &desc,
    ptree &pt) const
{
  //---- LOAD THE TEST PROBLEM FROM THE COMMAND LINE ARGUMENTS.
  NlpProblemDef *  pOptProb = knitro_wrapper::getTheInstance();
  bool             bWantToSolve = true;

  jtf_func_cons_ptr obj = desc->get_objective();

  const vector<jtf_func_cons_ptr> & eqn_cons = desc->get_eqn_constraint();

  knitro_wrapper::init_instance(
        dynamic_cast<knitro_wrapper*>(pOptProb),obj, eqn_cons, init_node);


  //---- OPEN A NEW INSTANCE OF KNITRO.
  KTR_context_ptr  kc;
  kc = KTR_new();
  if (kc == NULL) {
      cout << "*** KTR_new failed, maybe a license issue?\n";
      exit( EXIT_FAILURE );
    }

  //---- APPLY ANY USER OPTIONS (PROCEED EVEN IF THERE IS AN ERROR).
  KTR_load_param_file (kc, "knitro.opt");

  //---- LOAD THE PROBLEM INTO KNITRO.
  if (pOptProb->loadProblemIntoKnitro (kc) == false){
      cout << "*** loadProblemIntoKnitro failed\n";
      exit( EXIT_FAILURE );
    }

  //---- SET CALLBACK POINTERS FOR EVALUATION OF PROBLEM INFORMATION.
  //---- IF THE TEST CODE DOES NOT SUPPLY DERIVATIVES, THEN THE
  //---- USER OPTIONS IN "knitro.opt" SHOULD REQUEST AN ALTERNATIVE,
  //---- SUCH AS FINITE DIFFERENCES.
  g_pOptProblem = pOptProb;
  if(KTR_set_func_callback(kc, (KTR_callback * const) wrapperEvalFC) != 0){
      cout << "*** KTR_set_func_callback failed\n";
      exit( EXIT_FAILURE );
    }

  if(pOptProb->areDerivativesImplemented (nCAN_COMPUTE_GA) == true){
      if (KTR_set_grad_callback(kc,(KTR_callback * const) wrapperEvalGA) != 0){
          cout << "*** KTR_set_grad_callback failed\n";
          exit( EXIT_FAILURE );
        }
    }

  if((pOptProb->areDerivativesImplemented (nCAN_COMPUTE_H) == true)
     || (pOptProb->areDerivativesImplemented (nCAN_COMPUTE_HV) == true)){
      //---- SPECIFY THAT THE USER IS ABLE TO PROVIDE EVALUATIONS
      //---- OF THE HESSIAN MATRIX WITHOUT THE OBJECTIVE COMPONENT.
      //---- TURNED OFF BY DEFAULT BUT SHOULD BE ENABLED IF POSSIBLE.
      if(KTR_set_int_param(kc, KTR_PARAM_HESSIAN_NO_F, KTR_HESSIAN_NO_F_ALLOW) != 0)
        exit(EXIT_FAILURE);

      if(KTR_set_hess_callback(kc,(KTR_callback *const) wrapperEvalHorHV) != 0) {
          cout << "*** KTR_set_hess_callback failed\n";
          exit( EXIT_FAILURE );
        }
    }

  //---- ALLOCATE ARRAYS
  double *  daX      = new double[pOptProb->getN()];
  double *  daLambda = new double[pOptProb->getM() + pOptProb->getN()];

  if(bWantToSolve == true){
      //---- CALL KNITRO AND SOLVE.
      double  dFinalObj;
      int  nStatus = KTR_solve(kc, daX, daLambda, 0, &dFinalObj,
                               NULL, NULL, NULL, NULL, NULL, NULL);
      if(nStatus != 0)
        cerr << "# [error] Knitro failed to solve the problem, final status = "
             << nStatus << endl;
      else
        cerr << "# [info]  Knitro sucessful, finial obj = " << dFinalObj << endl;
    }else{
      //---- USE KNITRO TO CHECK THE DERIVATIVES CODED IN THE TEST PROBLEM.
      pOptProb->getInitialX (daX);
      KTR_check_first_ders(kc, daX, 2, 1.0e-14, 1.0e-14,
                           0, 0.0, NULL, NULL, NULL, NULL);
    }

  itr_matrix<double*> dax0(pOptProb->getN(),1, daX);
  init_node = dax0;

  delete [] daX;
  delete [] daLambda;

  KTR_free (&kc);

  // delete pOptProb;

  return( EXIT_SUCCESS );
}
