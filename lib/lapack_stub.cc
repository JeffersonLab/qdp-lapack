// -*- C++ -*-
// $Id: lapack_stub.cc,v 1.3 2009-02-04 21:26:59 kostas Exp $
/*! \file
 *  \brief Stubs for QDP interface to Lapack lib
 */

#include "qdp-lapack.h"

using namespace QDP;

namespace QDPLapack
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


  
  
  int zpotrf(char &uplo, int N,  multi2d<DComplex>& A, int LDA, int& info){
    
    QDP_error_exit("QDPLapack: %s not implemented", __func__);
    return 0;

  }
  int zpotrf(char &uplo, multi2d<DComplex>& A, int& info){
    int LDA = A.size1() ; // Assumes square matrix
    
    return zpotrf(uplo,LDA,A,LDA,info);
    
  }


  int cpotrf(char &uplo, int N,  multi2d<Complex>& A, int LDA, int& info){
    
    QDP_error_exit("QDPLapack: %s not implemented", __func__);
    return 0;

  }
  int cpotrf(char &uplo, multi2d<Complex>& A, int& info){
    int LDA = A.size1() ; // Assumes square matrix
    
    return cpotrf(uplo,LDA,A,LDA,info);
    
  }
  

  
  int zpotrs(char &uplo, int N, int nrhs,  
	     multi2d<DComplex>& A, int LDA, 
	     multi2d<DComplex>& B, int LDB, int& info){

    QDP_error_exit("QDPLapack: %s not implemented", __func__);
    return 0;
  }
  
  int zpotrs(char &uplo, multi2d<DComplex>& A,  multi2d<DComplex>& B,  
	     int& info){
    
    int LDA = A.size1() ; // Assumes square matrix
    int LDB = B.size1() ; // second index is fast and its size is size1
    int nrhs = B.size2() ;
    
    zpotrs(uplo,LDA,nrhs,A,LDA,B,LDB,info);
    
  }


  int cpotrs(char &uplo, int N, int nrhs,  
	     multi2d<Complex>& A, int LDA, 
	     multi2d<Complex>& B, int LDB, int& info){

    QDP_error_exit("QDPLapack: %s not implemented", __func__);
    return 0;
  
  }
  
  int cpotrs(char &uplo, multi2d<Complex>& A,  multi2d<Complex>& B,  
	     int& info){
    
    int LDA = A.size1() ; // Assumes square matrix
    int LDB = B.size1() ; // second index is fast and its size is size1
    int nrhs = B.size2() ;
    
    cpotrs(uplo,LDA,nrhs,A,LDA,B,LDB,info);
    
  }
  
  
  
} // namespace Lapack


