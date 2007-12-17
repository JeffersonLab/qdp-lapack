// -*- C++ -*-
// $Id: qdp-lapack.h,v 1.2 2007-12-17 19:33:08 kostas Exp $
/*! \file
 *  \brief QDP interface to Lapack lib
 */

#ifndef __qdp_lapack_h__
#define __qdp_lapack_h__

#include "qdp.h"

namespace QDPLapack
{
  // QDP Interface to Lapack functions
  /*!
   * @{
   */
  int zheev(char& jobz,
	    char& uplo,
	    //const int& N, // These are the dimensions of the array A
	    multi2d<DComplex>& A,
	    const int& lda,  // These are the dimensions of the array A
	    multi1d<Double>& w,
	    multi1d<DComplex>& Work, // Should be length LWork >= max(1,2*N-1)
	    //const int& LWork,
	    multi1d<Double>& RWork // Should be length max(1,3*N-2)
    ) ; 

  int zheev(char& jobz,
	    char& uplo,
	    //const int& N, // These are the dimensions of the array A
	    multi2d<DComplex>& A,
	    //const int& lda,  // These are the dimensions of the array A
	    multi1d<Double>& w
    ) ;

  int zheev(char& jobz,
	    char& uplo,
	    const int& N, // These are the dimensions of the array A
	    multi2d<DComplex>& A,
	    //const int& lda,  // These are the dimensions of the array A
	    multi1d<Double>& w
    ) ;

  int zgeqrf(const int M, // The vector length
	     const int N, // The number of vectors
	     multi2d<DComplex>& A, // the array containing the vectors
	     multi1d<DComplex>& TAU // some strange LAPACK beast
    );

  //! Applies the Q part of a QR factorization 
  /*!
   *  Purpose
   *  =======
   *
   *  ZUNMQR overwrites the general complex M-by-N matrix C with
   *
   *                  SIDE = 'L'     SIDE = 'R'
   *  TRANS = 'N':      Q * C          C * Q
   *  TRANS = 'C':      Q**H * C       C * Q**H
   *
   *  where Q is a complex unitary matrix defined as the product of k
   *  elementary reflectors
   *
   *        Q = H(1) H(2) . . . H(k)
   *
   *  as returned by ZGEQRF. Q is of order M if SIDE = 'L' and of order N
   *  if SIDE = 'R'.
   *
   */
  int zunmqr(char& side,//'L' or 'R': apply Q or Q^\dagger from the Left
	     char& trans,//'N' or 'C':  apply Q or Q^\dagger
	     const int M, // The vector length
	     const int N,  // The number of vectors
	     multi2d<DComplex>& A, //input MxN matrix
	     multi1d<DComplex>& TAU, // some strange LAPACK beast
	     multi2d<DComplex>& C // input,output
    );    


  /*--------------------------------------------------------------------
   *  ZHETRF( UPLO, N, A, IPIV)
   *  Double precision Cholesky factorization from Lapack
   *    	A = L*D*L**H
   *
   *  UPLO  'U' stores upper triangular / 'L' stores lower triangular
   *  N     order of A
   *  A     (input)  The matrix 
   *  	    (output) The upper/lower triangular factor stored in A
   *  IPIV  (output) INTEGER array, dimension (N) the pivot array
   *
   *  Work is allocated locally
   *--------------------------------------------------------------------*/
  int zhetrf( char &uplo, const int& n,
	      multi2d<DComplex>& A, 
	      multi1d<int>& ipiv
	      ) ;

/*--------------------------------------------------------------------
   * ZHETRS( UPLO, N, A, IPIV, B)
   *  
   *  Double precision solution of *only NRHS=1* linear system 
   *  using Cholesky factors A = L*D*L**H from ZHETRF
   *
   *  UPLO  'U' stores upper triangular / 'L' stores lower triangular
   *  N     order of A
   *  A     (input) The factor in U or L
   *  IPIV  (input) INTEGER array, dimension (N) the pivot array
   *  B	    one right hand size (NRHS = 1)
   *
   *--------------------------------------------------------------------*/
  int zhetrs( char &uplo, const int& n, 
	      multi2d<DComplex>& A, 
	      multi1d<int>& ipiv,
	      multi1d<DComplex>& B
	      ) ;

  // This is a BLAS routine 
  //     C := alpha*op( A )*op( B ) + beta*C,
  //
  // TRANSA TRANSB 'N' 'T' or 'C' (normal transpose dagger )
  int LatFermMat_x_Mat_cgemm(char& TRANSA,  
			     char& TRANSB,
			     const int& M,  // rows of A
			     const int& N,  // columns of B
			     const int& K,  // columns of A and rows of B
			     const Complex& ALPHA, 
			     const multi1d<LatticeFermion>& A, 
			     // The LDA is known: int *LDA,
			     const multi2d<Complex>& B, const int& LDB,
			     const Complex& BETA,
			     multi1d<LatticeFermion>& C
			     // The LDB is known: int *LDC
			     ) ;

  //Mixed precision
  int LatFermMat_x_Mat_cgemm(char& TRANSA, 
                             char& TRANSB,
                             const int& M, const int& N, const int& K, 
                             const DComplex& ALPHA,
                             const multi1d<LatticeFermion>& A,
                             // The LDA is known: int *LDA,
                             const multi2d<DComplex>& B, const int& LDB, 
                             const Complex& BETA,
                             multi1d<LatticeFermion>& C
                             // The LDB is known: int *LDC
                             ) ;

  /*--------------------------------------------------------------------
   * cgemm     	** this is a BLAS not a Lapack function **
   * 		COMPLEX 
   *            Full interface
   *
   *            C = alpha*A*B + beta*C
   *
   * 		A mxk,	B kxn, 	C mxn 
   * 		(or more accurately op(A),op(B),op(C) have these dimensions)
   * 		
   *--------------------------------------------------------------------*/
  int cgemm(char &transa,
	    char &transb,
	    int m,
	    int n,
	    int k,
	    Complex alpha,
	    multi2d<Complex>& A,  // input
	    int lda,
	    multi2d<Complex>& B,  // input
	    int ldb,
	    Complex beta,
	    multi2d<Complex>& C,  //input/output
	    int ldc   ) ;


/*--------------------------------------------------------------------
   * zgemm     	** this is a BLAS not a Lapack function **
   * 		DOUBLE COMPLEX 
   *            Full interface
   *
   *            C = alpha*A*B + beta*C
   *
   * 		A mxk,	B kxn, 	C mxn
   * 		(or more accurately op(A),op(B),op(C) have these dimensions)
   * 		
   *--------------------------------------------------------------------*/
  int zgemm(char &transa,
	    char &transb,
	    int m,
	    int n,
	    int k,
	    DComplex alpha,
	    multi2d<DComplex>& A,  // input
	    int lda,
	    multi2d<DComplex>& B,  // input
	    int ldb,
	    DComplex beta,
	    multi2d<DComplex>& C,  //input/output
	    int ldc   	    ) ;

  /*--------------------------------------------------------------------
   * CGEMV  	** this is a BLAS not a Lapack function **
   * 		single precision interface for 
   * 			y=A*x
   *  
   *  TRANS   'N' y := A*x,   'T' y := A'*x,   'C' y := conjg( A' )*x 
   *  M     rows of A
   *  N     columns of A
   *  A     The m x n matrix
   *  X     the vector of dim (n)
   *  Y     (output) the result
   *
   *--------------------------------------------------------------------*/
  int cgemv(char &trans, 
	    const int& m, const int& n,
	    multi2d<Complex>& A, 
	    multi1d<Complex>& x,
	    multi1d<Complex>& y
	    );


} // namespace QDPLapack

#endif 
