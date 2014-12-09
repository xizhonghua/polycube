// =============================================================================
// SpinXForm -- EigenSolver.h
// Keenan Crane
// August 16, 2011
//
// EigenSolver is used to solve the sparse eigenvalue problem
//
//    Ax = cx
//
// for the eigenvector x with smallest eigenvalue c, where A is a
// positive-definite matrix.  It implements the most basic inverse
// iteration scheme
//
//    x_{n+1} = A^-1 * x_n,
//
// which is equivalent to iteratively solving the linear system
//
//    A x_{n+1} = x_n.
//
// Currently the solver applies a fixed number of iterations and no
// other convergence criteria are considered.  Better results might be
// obtained by using a more sophisticated eigensolver such as SLEPc,
// PRIMME, or ARPACK (available as "eigs" in MATLAB).  Additionally, if
// a good initial guess x0 for the eigenvector is available, faster
// convergence might be achieved by computing the eigenvalue estimate
//
//    c0 = (x'*A*x)/(x'*x)
//
// and applying the same inverse iteration scheme to the shifted
// matrix B = A-c0*I.
//

#ifndef SPINXFORM_EIGENSOLVER_H
#define SPINXFORM_EIGENSOLVER_H

#include "QuaternionMatrix.h"
#include <vector>

using namespace std;

class EigenSolver
{
   public:
      static void solve( QuaternionMatrix& A,
                         vector<Quaternion>& x );
      // solves the eigenvalue problem Ax = cx for the
      // eigenvector x with the smallest eigenvalue c

   protected:
      static void normalize( vector<Quaternion>& x );
      // rescales x to have unit length
};

#endif
