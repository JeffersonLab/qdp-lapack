/*******************************************************************************
 * Function IncrEigpcg -- Incremental eigpcg. 
 *
 ******************************************************************************/
#include <assert.h>
#include <limits.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "qdp-lapack_IncrEigpcg.h"

/* DEBUG **********************************************************************/
void RayleighRitz(Complex_C *evecs, int lde, int n, int numEvecs,
     Complex_C *H, int ldh,
     FILE *outputFile, void (*matvec)(void *, void *, void *),  void *params)
{
     int Eelems = lde*numEvecs;
     int Helems = numEvecs*numEvecs;
     int lwork = 2*ldh;
     int ONE = 1;
     int i,info;
     Complex_C *V;
     Complex_C *Ht;
     Complex_C *cwork;
     float *Hevals;
     float *awork;
     if (lwork < n) lwork = n;

     float rnorm;
     char cV = 'V'; char cU = 'U'; 

     V = (Complex_C *) calloc(Eelems, sizeof(Complex_C));
     Ht = (Complex_C *) calloc(Helems, sizeof(Complex_C));
     cwork = (Complex_C *) calloc(lwork, sizeof(Complex_C));
     Hevals = (float *) calloc(numEvecs, sizeof(float));
     awork = (float *) calloc(3*ldh, sizeof(float));
	     
     BLAS_CCOPY(&Eelems, evecs, &ONE, V, &ONE);
     for (i=0;i<numEvecs;i++)
        BLAS_CCOPY(&numEvecs, &H[i*ldh], &ONE, &Ht[i*numEvecs], &ONE);

     BLAS_CHEEV(&cV,&cU,&numEvecs,Ht,&numEvecs,Hevals,cwork,&lwork,awork,&info);
     restart_X(V, lde, Ht, n, numEvecs, numEvecs, cwork, lwork);

     for (i=0;i<numEvecs;i++)  {
         computeResNorm(&V[i*lde], &Hevals[i], &rnorm, cwork, n, matvec,params);
         fprintf(outputFile, "RR Eval[%d]: %22.15E rnorm: %22.15E\n", 
			                i+1, Hevals[i], rnorm);
     }
	
     free(V);
     free(Ht);
     free(cwork);
     free(Hevals);
     free(awork);
}

/* end DEBUG ******************************************************************/

/******************************************************************************/
void IncrEigpcg(int n, int lde, /* n dim of matrix A, lde leading dim of vecs */
	int nrhs, 		/* The number of right hand sides to solve */
	Complex_C *X,   	/* The nrhs solutions */
	Complex_C *B,   	/* The nrhs right hand sides */
	int *ncurEvals,		/* Num of eigenpairs already in evecs,evals */
	int ldh,  		/* Max eigenpairs that can be stored in evecs */
	Complex_C *evecs, 	/* The eigenvector basis augmented by eigpcg */
	float *evals, 		/* The eigenvalues augmented by eigpcg */
	Complex_C *H,		/* The ncurEvals^2 matrix: H=evecs'*A*evecs; */
	Complex_C *HU,  	/* The Cholesky factorization of the H = U'U */
				/* ---------------------------------------- */
	void (*matvec)(void *, void *, void *),  /* Matvec operator */
	void (*precon)(void *, void *, void *),  /* Precondition operator */
	void *params,  		/* parameters to be passed to operators */
				/* ---------------------------------------- */
			/* Work arrays will be allocated if null on input */
	Complex_C *work,        /* work array for CG etc (max(4lde,nev*ldh) */
	Complex_C *V, 		/* work array for eigenvector search basis */
	Complex_C *ework, 	/* work array. If esize < N+2*nev allocate */
	int esize, 		/* N+2*nev <= esize <= (2*nev+1)*N */
				/* ---------------------------------------- */
	float tol, 		/* CG linear system tolerance |res| < tol*b */
	int maxit, 		/* Maximum number of CG iterations per system */
	int plvl, 		/* Printing level [0..5] */
	int nev, 		/* Number of eigenvalues to target in each CG */
	int v_max,		/* Maximum number of n-dim vectors in V */
	FILE *outputFile)       /* File to write reports. Could be STDERR */
{

  /* Timing vars */
  double wt1,wt2,ut1,ut2,st1,st2;

  /* Pointers */
  Complex_C *oldEwork = ework;
  Complex_C *x, *b;
  float *rnorms, *reshist;
 
  /* Variables */
  Complex_C tempc;
  int i,j, ONE = 1;
  int zs, cs, ds, tmpsize;
  int numIts, flag, nAdded;
  int freeV = 0;
  int freeWork = 0;
  int freeEwork = 0;
  float machEps = 2e-8;
  float lambda;
  int Cholesky = 0;   % 1=> Hinv by Cholesky 0=> Solve H eigenproblem

  Complex_C tpone = {+1.0e+00,+0.0e00}, tzero = {+0.0e+00,+0.0e00};
  char cR = 'R'; char cL = 'L'; char cN ='N'; 
  char cV = 'V'; char cU = 'U'; char cC ='C';

  /* ------------------------------------------------------------------- */
  /* Work allocations */
  /* ------------------------------------------------------------------- */
  zs = sizeof(Complex_Z); 
  cs = sizeof(Complex_C); 
  ds = sizeof(double);
  if (work == NULL) {
     tmpsize = 4*lde;
     if (tmpsize < nev*ldh)
	tmpsize = nev*ldh;
     if ((work = calloc(tmpsize, cs)) == NULL) 
        fprintf(stderr, "ERROR IncrEigpcg could not allocate work\n");
     freeWork = 1;
  }
  if (V == NULL && nev > 0) {
     if ((V = calloc(v_max*n, cs)) == NULL) 
        fprintf(stderr, "ERROR IncrEigpcg could not allocate V\n");
     freeV = 1;
  }
  if (ework == NULL && nev > 0) {
     if (esize < n+2*nev)
	esize = 3*n;      /* Defaulting to 3*nev */
     if ((ework = calloc(esize, cs)) == NULL) 
        fprintf(stderr, "ERROR IncrEigpcg could not allocate ework\n");
     freeEwork = 1;
  }
  else { /* Not null but check if enough space */
     if (esize < n+2*nev && nev > 0) { /* Realloc ework but remember old pntr */
        freeEwork = 1;
	esize = 3*n;      /* Defaulting to 3*nev and allocating */
        if ((ework = calloc(esize, cs)) == NULL) 
           fprintf(stderr, "ERROR IncrEigpcg could not allocate ework\n");
     }
  }
  if (  (rnorms = calloc(nev, sizeof(float)) ) == NULL )
     fprintf(stderr, "ERROR IncrEigpcg could not allocate rnorms\n");
  if (  (reshist = calloc(maxit, sizeof(float)) ) == NULL)
     fprintf(stderr, "ERROR IncrEigpcg could not allocate reshist\n");

  /* ------------------------------------------------------------------- */
  /* end Work allocations */
  /* ------------------------------------------------------------------- */

  /* ------------------------------------------------------------------- */
  /* Solving one by one the nrhs systems with incremental init-eigpcg    */
  /* ------------------------------------------------------------------- */
  for (j=0; j<nrhs; j++) {

      /* The j-th system */
      x = &X[j*lde];
      b = &B[j*lde]; 

      if (plvl) fprintf(outputFile, "System %d\n", j);

      while ( !thisSystemConverged ) {
         /* ------------------------------------------------------------ */
	 /* Always solve A eps = residual, even when xinit = 0           */
         /* ------------------------------------------------------------ */
         wt1 = primme_get_wtime(); 
         primme_get_time(&ut1,&st1);
   
         tempc = wrap_cdot(&n, x, &ONE, x, &ONE, params); /* norm ||x||^2 */
         if (tempc.r > 0.0) { /* If initial guess */
            matvec(x, work, params);  /* work(1:n) used for residual b-A*x */
            for (i = 0; i < n; i ++) {
                work[i].r = b[i].r - work[i].r;
                work[i].i = b[i].i - work[i].i;
            }
         }
         else 
  	    BLAS_CCOPY(&n, b, &ONE, work, &ONE);

         /* ------------------------------------------------------------ */
         /* Perform init-CG on A eps = res with evecs vectors            */
         /* ------------------------------------------------------------ */
         /* xinit = evecs*Hinv*evec'*(b-Ax0) */
   
         if (*ncurEvals > 0) {
            /* work(n:n+curEvals) = evecs(n x ncurEvals)'*work(1:n) */
            wrap_cgemv(&cC,&n,ncurEvals,&tpone,evecs,&lde,work,&ONE,
	           &tzero,&work[n],&ONE,&work[2*n],params);
   
	    if (Cholesky) 
            /* HU is the upper triangular factor of the Cholesky */
            /* work(n:n+ncurEvals)=H^(-1)*work(n:n+ncurEvals)=(evecs'*work) */
	       BLAS_CPOTRS(&cU,ncurEvals,&ONE,HU,&ldh,&work[n],ncurEvals,&flag);
	    else {
            /* HU is the eigenvectors of the matrix H */
            /* work(n:n+ncurEvals) = HU*Lambda^(-1)*HU'*work(n:n+ncurEvals) */
               BLAS_CGEMV(&cC,ncurEvals,ncurEvals,&tpone,HU,&ldh,&work[n],&ONE,
		       &tzero,work,&ONE);
	       for (i = 0; i<*ncurEvals ; i++) {
	           work[i].r = work[i].r/evals[i];
	           work[i].i = work[i].i/evals[i];
	       }
               BLAS_CGEMV(&cN,ncurEvals,ncurEvals,&tpone,HU,&ldh,work,&ONE,
		       &tzero,&work[n],&ONE);
	    }
            /* work[1:n] = evecs*work(n:n+ncurEvals) */
            BLAS_CGEMV(&cN,&n,ncurEvals,&tpone,evecs,&lde,&work[n],&ONE,
		       &tzero,work,&ONE);
            /* x = x + work the new initial guess */
	    BLAS_CAXPY(&n, &tpone, work, &ONE, x, &ONE);
   
            wt2 = primme_get_wtime();
            primme_get_time(&ut2,&st2);
	    if (plvl) {
               fprintf(outputFile, "InitCG\n");
               fprintf(outputFile, "I Wallclock : %-f\n", wt2-wt1);
               fprintf(outputFile, "I User Time  : %f seconds\n", ut2-ut1);
               fprintf(outputFile, "I Syst Time  : %f seconds\n", st2-st1);
	    }
         }
         /* ------------------------------------------------------------ */
         /* end of init-CG with evecs vectors                            */
         /* ------------------------------------------------------------ */
   
         /* ------------------------------------------------------------ */
         /* Adjust nev because we can only store ldh eigenvectors        */
         /* ------------------------------------------------------------ */
         if (ldh-*ncurEvals < nev)
	    nev = ldh - *ncurEvals;
         
         /* ------------------------------------------------------------ */
         /* If simple CG (nev == 0) compute how often to restart CG      */
         /* ------------------------------------------------------------ */
          if (nev > 0) {
             RemainingTol = bnorm*tol/resNorm;
   		continue here:
          }

         /* ------------------------------------------------------------ */
         /* Solve Ax = b with x initial guess                            */
         /* ------------------------------------------------------------ */
         wt1 = primme_get_wtime(); 
         primme_get_time(&ut1,&st1);
   
         eigpcg(n, lde, x, b, tol, maxit, plvl, &numIts, reshist, 
	        &flag, work, matvec, precon, params,
                nev, &evals[*ncurEvals], rnorms, v_max, V, esize, ework);
   
         wt2 = primme_get_wtime();
         primme_get_time(&ut2,&st2);

      /* ---------- */
      /* Reporting  */
      /* ---------- */

      if (plvl) {
         fprintf(outputFile, "eigCG\n");
         fprintf(outputFile, "E Wallclock 	: %-f\n", wt2-wt1);
         fprintf(outputFile, "E User Time    : %f seconds\n", ut2-ut1);
         fprintf(outputFile, "E Syst Time    : %f seconds\n", st2-st1);
         fprintf(outputFile, "Iterations: %-d\n", numIts); 
         fprintf(outputFile, "Actual Resid of LinSys  : %e\n", reshist[numIts]);
	 if (plvl > 1) 
            for (i=0; i < nev; i++) 
               fprintf(outputFile, "Eval[%d]: %-22.15E rnorm: %-22.15E\n", 
			                i+1, evals[*ncurEvals+i], rnorms[i]); 
   
         if (flag != 0) {
            fprintf(outputFile, 
               "Error: eigpcg returned with nonzero exit status\n");
            return;
         }
      }
      /* ------------------------------------------------------------------- */

      /* ------------------------------------------------------------------- */
      /* Update the evecs and the factorization of evecs'*A*evecs            */
      /* ------------------------------------------------------------------- */
      wt1 = primme_get_wtime(); 
      primme_get_time(&ut1,&st1);

      if (nev > 0) {
         /* Append new Ritz pairs to evecs */
         for (i=0; i<nev; i++)
	     BLAS_CCOPY(&n, &V[i*n], &ONE, &evecs[(*ncurEvals+i)*lde], &ONE);
   
         /* Orthogonalize the new Ritz vectors. This is necessary only if
	  * we eigen-decompose H. With Cholesky it only reduces conditioning */
	 if (Cholesky)
            nAdded = ortho_new_vectors(evecs, lde, n, *ncurEvals, 2,
                &evecs[(*ncurEvals)*lde], nev, machEps, ework, work, params);
	 else
            nAdded = ortho_new_vectors(evecs, lde, n, *ncurEvals, 3,
                &evecs[(*ncurEvals)*lde], nev, machEps, ework, work, params);

         /* V = A*evecs(new) */
         for (i=0; i<nAdded; i++) 
	     matvec(&evecs[(*ncurEvals+i)*lde], &V[i*n], params);
   
         /* Hnew = evecs(old+new)'*V(new) = evecs(old+new)'*A*evecs(new) */
         tmpsize = (*ncurEvals)+nAdded;
         wrap_cgemm(&cC, &cN, &tmpsize, &nAdded, &n, &tpone, evecs, &lde, 
	         V, &n, &tzero, &H[(*ncurEvals)*ldh], &ldh, work, params);
   
         *ncurEvals = *ncurEvals + nAdded;
   
         /* Copy H into HU */
         for (i=0; i<*ncurEvals; i++)
     	    BLAS_CCOPY(ncurEvals, &H[i*ldh], &ONE, &HU[i*ldh], &ONE);
   
	 if (Cholesky)  /* Cholesky factorize H = HU'*HU */
            BLAS_CPOTRF(&cU, ncurEvals, HU, &ldh, &flag);
	 else {         /* Eigendecompose: H = HU Lambda HU'. Don't form V(HU)*/
	    i = 4*lde;
            BLAS_CHEEV(&cV,&cU,ncurEvals,HU,&ldh,evals,work,&i,(float *)ework,
		       &flag);
	 }

	 /* Reporting */
         wt2 = primme_get_wtime();
         primme_get_time(&ut2,&st2);
         if (plvl) {
            fprintf(outputFile, "Update\n");
            fprintf(outputFile, "Added %d vecs\n",nAdded);
            fprintf(outputFile, "U Wallclock : %-f\n", wt2-wt1);
            fprintf(outputFile, "U User Time  : %f seconds\n", ut2-ut1);
            fprintf(outputFile, "U Syst Time  : %f seconds\n", st2-st1);
	    if (plvl >= 2 ) {
	       if (Cholesky) RayleighRitz(evecs,lde,n,*ncurEvals,H,ldh,
			       		  outputFile,matvec,params);
	       else {
     		  for (i=0;i<*ncurEvals;i++)  {
                     BLAS_CGEMV(&cN,&n,ncurEvals,&tpone,evecs,&lde,&HU[i*ldh],
			       &ONE,&tzero,work,&ONE);
         	     computeResNorm(work, &lambda, rnorms,
				      &work[n], n, matvec, params);
         	     fprintf(outputFile,"RR Eval[%d]: %22.15E rnorm: %22.15E\n",
				      i+1, lambda, rnorms[0]);
		  }
	       } /* Printing if evals are known */
	    } /* plvl >=5 */
         } /* plvl > 0 */
      } /* if nev>0 */
      /* ------------------------------------------------------------------- */
      /* end of update phase */
      /* ------------------------------------------------------------------- */

   } /* end of nrhs loop */


   if (freeEwork) {
      free(ework);
      ework = oldEwork;
   }
   if (freeWork) free(work);
   if (freeV) free(V);
   free(rnorms);
   free(reshist);

   return;
}
/******************************************************************************/
/* END OF IncrEigpcg */
/******************************************************************************/
