// -*- C++ -*-
// $Id: lapack_stub.cc,v 1.1 2007-10-11 20:21:11 edwards Exp $
/*! \file
 *  \brief Stubs for QDP interface to Lapack lib
 */

#include "qdp-lapack.h"

using namespace QDP;

namespace Lapack
{
  int zheev(char& jobz,
	    char& uplo,
	    const int& N,
	    multi2d<DComplex>& A,
	    //const int& lda,
	    multi1d<Double>& w)
  {
    QDP_error_exit("QDPLapack: %s not implemented", __func__);
    return 0;
  }

  int zheev(char& jobz,
	    char& uplo,
	    //const int& N,
	    multi2d<DComplex>& A,
	    //const int& lda,
	    multi1d<Double>& w,
	    multi1d<DComplex>& Work,
	    //const int& LWork,
	    multi1d<Double>& RWork)
  {
    QDP_error_exit("QDPLapack: %s not implemented", __func__);
    return 0;
  }

  int zheev(char& jobz,
	    char& uplo,
	    //const int& N, // These are the dimensions of the array A
	    multi2d<DComplex>& A,
	    //const int& lda,
	    multi1d<Double>& w)
  {
    QDP_error_exit("QDPLapack: %s not implemented", __func__);
    return 0;
  }

  int zgeqrf(const int M,
	     const int N,
	     multi2d<DComplex>& A,
	     multi1d<DComplex>& TAU)
  {
    QDP_error_exit("QDPLapack: %s not implemented", __func__);
    return 0;
  }

  int zunmqr(char& side,
	     char& trans,
	     const int M,
	     const int N,
	     multi2d<DComplex>& A,
	     multi1d<DComplex>& TAU,
	     multi2d<DComplex>& C)
  {
    QDP_error_exit("QDPLapack: %s not implemented", __func__);
    return 0;
  }
  
} // namespace Lapack


