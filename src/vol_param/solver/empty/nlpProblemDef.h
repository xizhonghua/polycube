/*******************************************************/
/* Copyright (c) 2006-2011 by Ziena Optimization LLC   */
/* All Rights Reserved                                 */
/*******************************************************/

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//++  Nonlinear optimization test problems (NLPs) must implement
//++  all the methods of this abstract interface:
//++    getN()
//++    getM()
//++    loadProblemIntoKnitro()
//++    areDerivativesImplemented()
//++    evalFC()
//++    evalGA()
//++    evalH()
//++    evalHV()
//++
//++  If the implementation does not include derivatives, then the
//++  corresponding "evaluate" methods can return NULL.  For example,
//++  if there are no second derivatives, then evalH() should return
//++  NULL for daH and evalHV() should return NULL for daHV.
//++
//++  Each test problems must also export a C function named
//++  "getNlpProblemDef()" that returns a pointer to the single instance
//++  of the constructed object.
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#ifndef NLPPROBLEMDEF_H
#define NLPPROBLEMDEF_H


//--------------------------------------------------------------------
//  Includes
//--------------------------------------------------------------------

#ifndef KNITRO_H__
#include  <knitro/knitro.h>
#endif


//--------------------------------------------------------------------
//  Defines and Typedefs
//--------------------------------------------------------------------

//---- CONSTANTS FOR THE METHOD areDerivativesImplemented.
typedef enum
{
    nCAN_COMPUTE_GA,
    nCAN_COMPUTE_H,
    nCAN_COMPUTE_HV
} DerivativesImplementedType;


//--------------------------------------------------------------------
//  Class Declaration:  NlpProblemDef
//--------------------------------------------------------------------
//++
//++ This abstract class defines a nonlinear optimization problem
//++ to be solved by KNITRO.  Several methods are declared pure virtual.
//++
class NlpProblemDef
{
  public:

    //++ Abstract classes require a virtual destructor.
    virtual ~NlpProblemDef(){}

    //++ Returns the number of unknowns.
    virtual int  getN (void) = 0;

    //++ Returns the number of constraints, counting equalities and
    //++ inequalities, but not variable bounds.
    virtual int  getM (void) = 0;

    //++ Fills the n-vector "daX" with the current value of the variables.
    //++ Do not call until after calling "loadProblemIntoKnitro".
    virtual void  getInitialX (double * const  daX) = 0;

    //++ Pass the fixed problem definition information to KNITRO by
    //++ calling KTR_init_problem.  This cannot be called more than once.
    //++ Returns true if successful.
    //++   I   kc - KNITRO context, already allocated
    virtual bool  loadProblemIntoKnitro (KTR_context_ptr  kc) = 0;

    //++ Returns true if the requested derivatives are implemented
    //++ as analytic operations.
    //++   I   nWhichDers - derivative being queried
    virtual bool  areDerivativesImplemented
                      (const DerivativesImplementedType  nWhichDers) = 0;

    //++ Evaluate "dObj" and "daC" at "daX".  KNITRO allows first derivatives
    //++ to also be computed at this time, but this implementation does not
    //++ provide that generality.
    //++ Returns 0 if successful.
    //++   I   daX        - n-vector at which to make the evaluation
    //++    O  dObj       - objective function at x
    //++    O  daC        - m-vector of constraint functions at x
    //++   IO  userParams - place to preserve state information between calls
    virtual int  evalFC (const double * const  daX,
                               double * const  dObj,
                               double * const  daC,
                               void   *        userParams) = 0;

    //++ Evaluate "daG" and "daJ" at "daX".  The objective gradient is a
    //++ dense vector, while the Jacobian gradients must match the sparsity
    //++ pattern defined by "loadProblemIntoKnitro".
    //++ Returns  0 if successful.
    //++ Returns -1 if first derivatives are not implemented.
    //++   I   daX        - n-vector at which to make the evaluation
    //++    O  daG        - objective gradient at x
    //++    O  daJ        - nnzJ-vector of constraint gradients at x
    //++   IO  userParams - place to preserve state information between calls
    virtual int  evalGA (const double * const  daX,
                               double * const  daG,
                               double * const  daJ,
                               void   *        userParams) = 0;

    //++ Evaluate "daH" at "daX, daLambda", where H is the Hessian of the
    //++ Lagrangian sigma*f(x) + lambda^T c(x).
    //++ Returns  0 if successful.
    //++ Returns -1 if second derivatives are not implemented.
    //++   I   daX        - n-vector at which to make the evaluation
    //++   I   daLambda   - m-vector at which to make the evaluation
    //++   I   dSigma     - scales f(x) in the Lagrangian (either 0 or 1)
    //++    O  daH        - nnzH-vector of Hessian values at (x,lambda)
    //++   IO  userParams - place to preserve state information between calls
    //++ Note that daLambda actually has length m+n, but only the first m
    //++ components should be used.
    virtual int  evalH (const double * const  daX,
                        const double * const  daLambda,
                        const double          dSigma,
                              double * const  daH,
                              void   *        userParams) = 0;

    //++ Evaluate "daHV" at "daX, daLambda", where H is the Hessian of the
    //++ Lagrangian sigma*f(x) + lambda^T c(x).
    //++ Returns  0 if successful.
    //++ Returns -1 if second derivatives are not implemented.
    //++   I   daX        - n-vector at which to make the evaluation
    //++   I   daLambda   - m-vector at which to make the evaluation
    //++   I   dSigma     - scales f(x) in the Lagrangian (either 0 or 1)
    //++   IO  daHV       - n-vector
    //++   IO  userParams - place to preserve state information between calls
    //++ Note that daLambda actually has length m+n, but only the first m
    //++ components should be used.
    virtual int  evalHV (const double * const  daX,
                         const double * const  daLambda,
                         const double          dSigma,
                               double * const  daHV,
                               void   *        userParams) = 0;

  protected:

    int       _nN;
    double *  _daXInit;
    double *  _daXLo;
    double *  _daXUp;
    int       _nM;
    int    *  _naCType;
    double *  _daCLo;
    double *  _daCUp;
    int       _nNnzJ;
    int    *  _naJacIndexVars;
    int    *  _naJacIndexCons;
    int       _nNnzH;
    int    *  _naHessCols;
    int    *  _naHessRows;
};


//--------------------------------------------------------------------
//  Class Definition:  NlpProblemDef
//--------------------------------------------------------------------
//++
//++ The class is abstract, but the virtual destructor requires an
//++ implementation.  This is the complete definition for the class.
//++


//---- DECLARE A FUNCTION THAT RETURNS A POINTER TO THE CONSTRUCTED OBJECT.
//---- SUBCLASSES MUST EXPORT A FUNCTION OF THIS TYPE THAT IS COMPATIBLE
//---- WITH THE C LANGUAGE, AND IS NAMED "getNlpProblemDef()".
typedef  NlpProblemDef * (*FnType_getNlpProblemDef) (void);


#endif     //-- NLPPROBLEMDEF_H

//++++++++++++++++ End of source code ++++++++++++++++++++++++++++++++
