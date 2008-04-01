/*******************************************************************************
 *   Adapted from PRIMME: Feb 7, 2008
 *
 *   PRIMME PReconditioned Iterative MultiMethod Eigensolver
 *   Copyright (C) 2005  James R. McCombs,  Andreas Stathopoulos
 *
 *   This file is part of PRIMME.
 *
 *   PRIMME is free software; you can redistribute it and/or
 *   modify it under the terms of the GNU Lesser General Public
 *   License as published by the Free Software Foundation; either
 *   version 2.1 of the License, or (at your option) any later version.
 *
 *   PRIMME is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 *   Lesser General Public License for more details.
 *
 *   You should have received a copy of the GNU Lesser General Public
 *   License along with this library; if not, write to the Free Software
 *   Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
 *
 * File: numerical.c
 *
 * Purpose - This file contains for the most part C wrapper routines for
 *    calling various BLAS and LAPACK FORTRAN routines.
 *
 ******************************************************************************/

#include <stdio.h>
#include <stdarg.h>
#include "qdp-lapack_Complex.h"
#include "qdp-lapack_numerical_private.h"
#include "qdp-lapack_numerical.h"
#ifdef USE_QMP
#include "qmp.h"
#endif
/******************************************************************************/
/* These wrappers perform the blas functions on a parallel distributed array
 * and apply the global sum of the resulting values through a user provided
 * global sum function, eg., MPI_Allreduce().
 *
 * The dot products do not call BLAS directly but a Fortran wrapper 
 * subroutine because Fortran functions cannot return complex to a C caller.
 *
 * These functions can be replaced with functions provided in libraries with 
 * their own dot+globalsum functions
 *
 ******************************************************************************/

/******************************************************************************/
void wrap_cgemv(char *transa, int *m, int *n, Complex_C *alpha, Complex_C *a,
   int *lda, Complex_C *x, int *incx, Complex_C *beta, Complex_C *y, int *incy,
   Complex_C *work, void *params)
{

   BLAS_CGEMV(transa, m, n, alpha, a, lda, x, incx, beta, work, incy);

   globalSumDouble(work, y, n, params );

}

/******************************************************************************/
/* Perform result = x'*y where all are single precision but the summation 
 * is performed in double precision */
Complex_C wrap_zsum_cdot(int *n, Complex_C *x, int *incx, Complex_C *y, int *incy, void *params) 
{
   int i;
   Complex_C cdotc, xconj, prod;
   Complex_Z sum, gsum;

   sum.r = 0.0;sum.i = 0.0;
   for (i=0;i<*n;i++) {
       s_cnjg_primme(&xconj,&x[i]);
       c_mul_primme(&prod,&xconj,&y[i]);
       sum.r = sum.r + prod.r;
       sum.i = sum.i + prod.i;
   }
   i = 2;  /* 1 double complex = 16 bytes = 2 doubles */
   globalSumDouble(&sum, &gsum, &i, params);
   cdotc.r = gsum.r;
   cdotc.i = gsum.i;
   return(cdotc);

}
/******************************************************************************/
Complex_C wrap_cdot(int *n, Complex_C *x, int *incx, Complex_C *y, int *incy,
   void *params) 
{
   int ONE = 1;
   Complex_C cdotc_r, cdotc;
   CDOTCSUB(&cdotc_r, n, x, incx, y, incy);
   //cblas_cdotc_sub(*n, x, *incx, y, *incy, &cdotc_r);
   globalSumDouble(&cdotc_r, &cdotc, &ONE, params);
   return(cdotc);

}
/******************************************************************************/
void wrap_cgemm(char *transa, char *transb, int *m, int *n, int *k,
   Complex_C *alpha, Complex_C *a, int *lda, Complex_C *b, int *ldb,
   Complex_C *beta, Complex_C *c, int *ldc,
   Complex_C *work, void *params)
{
   int i;

   BLAS_CGEMM(transa, transb, m, n, k, alpha, a, lda, b, ldb, beta, work, m);

   if (*ldc == *m) {
      i = (*m)*(*n);
      globalSumDouble(work, c, &i, params);
   }
   else {
      /* Sum all n columns only up to m, not to ldc */
      for (i=0; i<*n; i++)
         globalSumDouble(&work[i*(*m)], &c[i*(*ldc)], m, params);
   }

}
/******************************************************************************/
/* Double precision */
/******************************************************************************/
void wrap_zgemv(char *transa, int *m, int *n, Complex_Z *alpha, Complex_Z *a,
   int *lda, Complex_Z *x, int *incx, Complex_Z *beta, Complex_Z *y, int *incy,
   Complex_C *work, void *params)
{

   BLAS_CGEMV(transa, m, n, alpha, a, lda, x, incx, beta, work, incy);

   globalSumDouble(work, y, n, params );

}

/******************************************************************************/
Complex_Z wrap_zdot(int *n, Complex_Z *x, int *incx, Complex_Z *y, int *incy,
   void *params) 
{
   int ONE=1;
   Complex_Z zdotc_r, zdotc;
   ZDOTCSUB(&zdotc_r, n, x, incx, y, incy);
   //cblas_zdotc_sub(*n, x, *incx, y, *incy, &zdotc_r);
   globalSumDouble(&zdotc_r, &zdotc, &ONE, params);
   return(zdotc);

}

/******************************************************************************
 * function 
 * primme_seq_globalSumDouble(void *sendBuf, double *recvBuf, int count) 
 *
 * This is the sequential default for the function globalSumDouble. 
 * If the program is parallel, the user must replace this with 
 * an Allreduce() function
 * 
 ******************************************************************************
 *        NOTE: The count and the copying refers to double datatypes
 ******************************************************************************/
#ifdef USE_QMP
void globalSumDouble(void *sendBuf, void *recvBuf, int *count, void *params) {
   int ONE = 1;
   QMP_status_t f;

   BLAS_DCOPY(count, (double *) sendBuf, &ONE, (double *) recvBuf, &ONE);
   f = QMP_sum_double_array((double* ) recvBuf, *count);
}
#else
void globalSumDouble(void *sendBuf, void *recvBuf, int *count, void *params) 
{
   int ONE = 1;

   BLAS_DCOPY(count, (double *) sendBuf, &ONE, (double *) recvBuf, &ONE);

}
#endif

