// =============================================================================
// SpinXForm -- QuaternionMatrix.h
// Keenan Crane
// August 16, 2011
//
// QuaternionMatrix is a convenience class for representing sparse matrices with
// quaternion-valued entries.  Entries can be accessed via expressions like
//
//    A( i, j ) = Quaternion( 1., 2., 3., 4. );
//
// A QuaternionMatrix can be converted to a sparse matrix with real-valued
// entries by calling toReal().  Better performance could probably be achieved
// by building matrices directly in real, compressed-column format, but
// currently matrix construction is not a bottleneck.
//

#ifndef SPINXFORM_QUATERNIONMATRIX_H
#define SPINXFORM_QUATERNIONMATRIX_H

#include <map>
#include <vector>
#include <iostream>
#include "Quaternion.h"
#include "pcg/sparse_matrix.h"

class QuaternionMatrix
{
   public:
      void resize( int m, int n );
      // allocates an mxn matrix of zeros
      
      int size( int dim ) const;
      // returns the size of the dimension specified by scalar dim

            Quaternion& operator()( int row, int col );
      const Quaternion& operator()( int row, int col ) const;
      // access element (row,col)
      // note: uses 0-based indexing

      const SparseMatrixd& toReal( void );
      // returns real matrix where each quaternion becomes a 4x4 block

   protected:
      typedef std::pair<int,int> EntryIndex; // NOTE: column THEN row! (makes it easier to build compressed format)
      typedef std::map<EntryIndex,Quaternion> EntryMap;

      EntryMap data;
      // non-zero entries

      int m, n;
      // rows, columns

      static Quaternion zero;
      // dummy value for const access of zeros

      SparseMatrixd A;
      // real representation
};

#endif
