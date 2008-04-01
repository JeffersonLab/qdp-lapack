//MAC OS X
//#include <vecLib/vBLAS.h>
//#include <veclib/clapack.h>

// May need other includes for other OS 
#include "qdp.h"

extern "C" {

  //LAPACK routines
  int zheev_(char *,char *,
	     int *, DComplex *, int *, Double *, DComplex *, 
	     int *, Double *, int *);
  
  int  zgeqrf_(int *, int *,
	       DComplex *, int *,
	       DComplex *, DComplex *,
	       int *, int *) ;
  
  int zunmqr_(char *, char *,
	      int *, int *,
	      int *, DComplex *,
	      int *, DComplex *,  DComplex *,
	      int *, DComplex *,
	      int *, int *);

  int zhetrf_(char *uplo, int *n, DComplex *a, 
	      int *lda, int *ipiv, DComplex *work, int *lwork, 
	      int *info);
 
  int zhetrs_(char *uplo, int *n, int *nrhs, 
	      DComplex *a, int *lda, int *ipiv, DComplex *b, 
	      int *ldb, int *info);

  //BLAS routines  
  int cgemm_(char *TRANSA, char *TRANSB,int *M,int *N,int *K, 
             Complex *ALPHA,
             Complex *A, int *LDA, 
             Complex *B, int *LDB,
             Complex *BETA,
             Complex *C, int *LDC) ;

  int zgemm_(char *TRANSA, char *TRANSB,int *M,int *N,int *K, 
             DComplex *ALPHA,
             DComplex *A, int *LDA, 
             DComplex *B, int *LDB,
             DComplex *BETA,
             DComplex *C, int *LDC) ;

  int cgemv_(char *trans, int *M, int *N, Complex *alpha, 
	       Complex *A, int *lda, Complex *x, int *incx, 
	       Complex *beta, Complex *y, int *incy);

  int zgemv_(char *trans, int *M, int *N, DComplex *alpha, 
	     DComplex *A, int *lda, DComplex *x, int *incx, 
	     DComplex *beta, DComplex *y, int *incy);



}

