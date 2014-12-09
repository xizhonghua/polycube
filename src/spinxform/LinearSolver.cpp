// =============================================================================
// SpinXForm -- LinearSolver.cpp
// Keenan Crane
// August 16, 2011
//

#include "LinearSolver.h"
#include "pcg/pcg_solver.h"

void LinearSolver :: solve( QuaternionMatrix& A,
                     vector<Quaternion>& x,
                     const vector<Quaternion>& b,
                     bool precondition )
// solves the linear system Ax = b where A is positive-semidefinite
{
   if( precondition == false )
   {
      cerr << "WARNING: using basic CG solver with diagonal preconditioner -- may be (very) slow!" << endl;
   }

   int n = x.size();
   vector<double> result( n*4 );
   vector<double> rhs( n*4 );

   // convert right-hand side to real values
   toReal( b, rhs );

   // setup solver
   PCGSolver<double> pcg;
   const int max_iterations = 2000;
   const double tolerance_factor = 1e-7;
   const double mic_parameter = .97;
   const double min_diagonal_ratio = .25;
   double residual;
   int iterations;

   pcg.set_solver_parameters( tolerance_factor, max_iterations, mic_parameter, min_diagonal_ratio, precondition );
   
   // solve real linear system
   pcg.solve( A.toReal(), rhs, result, residual, iterations );

   cout << "Linear solver achieved a residual of " << residual;
   cout << " after " << iterations << " iterations." << endl;

   // convert solution back to quaternions
   toQuat( result, x );
}

void LinearSolver :: toReal( const vector<Quaternion>& uQuat,
                             vector<double>& uReal )
// converts vector from quaternion- to real-valued entries
{
   for( size_t i = 0; i < uQuat.size(); i++ )
   {
      uReal[i*4+0] = uQuat[i].re();   // real
      uReal[i*4+1] = uQuat[i].im().x; // i
      uReal[i*4+2] = uQuat[i].im().y; // j
      uReal[i*4+3] = uQuat[i].im().z; // k
   }
}

void LinearSolver :: toQuat( const vector<double>& uReal,
                             vector<Quaternion>& uQuat )
// converts vector from real- to quaternion-valued entries
{
   for( size_t i = 0; i < uQuat.size(); i++ )
   {
      uQuat[i] = Quaternion( uReal[i*4+0],   // real
                             uReal[i*4+1],   // i
                             uReal[i*4+2],   // j
                             uReal[i*4+3] ); // k
   }
}

