/*******************************************************************************
 * Subroutine restart_X - This subroutine computes X*hVecs and places 
 *    the result in X.
 *
 * INPUT ARRAYS AND PARAMETERS
 * ---------------------------
 * ldx 		Leading dimension of X
 *
 * nLocal       Number of rows of V assigned to the node
 *
 * basisSize    Current size of the basis V
 *
 * restartSize  Number of Ritz vectors V/W will be restarted with 
 *
 * rwork        Work array that must be at least of size restartSize
 *
 * rworkSize    The size availble in rwork. Matrix multiply blocks of X with 
 * 		hVecs, producing blocks of the new X of size
 * 		(AvailRows * restartSize) = rworkSize
 * 		Therefore rworkSize must be at least restartSize.
 *
 * INPUT/OUTPUT ARRAYS
 * -------------------
 * X      Holds either V or W before and after restarting
 *
 * hVecs  The eigenvectors of V'*A*V before and after restarting
 *
 ******************************************************************************/
// Commented out lines with // can be used to time the BLAS functions
//#include <stdio.h>
//#include "TEST/wtime.h"
#include "qdp-lapack_Complex.h"
#include "qdp-lapack_numerical_private.h"
#define min(a, b) (a < b ? a : b)
  
void restart_X(Complex_C *X, int ldx, Complex_C *hVecs, int nLocal, 
   int basisSize, int restartSize, Complex_C *rwork, int rworkSize) {

//		double st1, st2, ut1, ut2, wt1, wt2, uTotMM, uTotCP, wTotMM, wTotCP;
   char cN = 'N';
   int ONE = 1;
   int i, k;  /* Loop variables */
   int AvailRows = min(rworkSize/restartSize, nLocal);
   Complex_C tpone = {+1.0e+00,+0.0e00}, tzero = {+0.0e+00,+0.0e00};
   i = 0;

//   				wTotMM = 0.0;wTotCP = 0.0;
//   				uTotMM = 0.0;uTotCP = 0.0;
   while (i < nLocal) {
      /* Block matrix multiply */

//   wt1 = primme_get_wtime();
//   primme_get_time(&ut1,&st1);

      BLAS_CGEMM(&cN, &cN, &AvailRows, &restartSize, &basisSize, &tpone,
         &X[i], &ldx, hVecs, &basisSize, &tzero, rwork, &AvailRows );

//   wt2 = primme_get_wtime();
//   primme_get_time(&ut2,&st2);
//   wTotMM += (wt2-wt1);
//   uTotMM += ut2-ut1;


      /* Copy the result in the desired location of X */
//   wt1 = primme_get_wtime();
//   primme_get_time(&ut1,&st1);

      for (k=0; k < restartSize; k++) {
         BLAS_CCOPY(&AvailRows, &rwork[AvailRows*k],&ONE, &X[i+ldx*k],&ONE);
      }

//   wt2 = primme_get_wtime();
//   primme_get_time(&ut2,&st2);
//   wTotCP += wt2-wt1;
//   uTotCP += ut2-ut1;

      i = i+AvailRows;
      AvailRows = min(AvailRows, nLocal-i);
   }

//   printf("Wallclock GEMM   : %f\n", wTotMM);
//   printf("User Time GEMM   : %f seconds\n",uTotMM);
//   printf("Wallclock COPY   : %f\n", wTotCP);
//   printf("User Time COPY   : %f seconds\n", uTotCP);

}
