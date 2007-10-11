// -*- C++ -*-
// $Id: qdp-lapack.h,v 1.1 2007-10-11 20:22:00 edwards Exp $
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

} // namespace QDPLapack

#endif 
