#ifndef LIK_IPP

#include <math.h>
#include <float.h>
#include "matrix.ipp"
#include "cov.ipp"
#include "mean.ipp"

//#include "mex.h"


template <class T>
T log_marginal_likelihood(
		T* &f, T* &s2, T** &X, int &n, int &d,
		T (*covFunc) ( T* &, T* &, int &, T* &),
		T* &theta)
{
	T lik = 0;
	
	T **K = allocMatrix<T>( n, n);
	covMatrix<T>( K, X, n,d, covFunc, theta);
	//fprintf( stderr, "Covariance Matrix: K =\n");
	//printMatrix( K, n, n);
	
	// noise
	if( s2==0)
		for(int i=0; i<n; i++)
			K[i][i] += pow(10, theta[0]);
	else
		for(int i=0; i<n; i++)
			K[i][i] += s2[i] + pow(10, theta[0]);

	//T   detK = determinantCholesky( K, n);
	T   detK = determinantLU( K, n);
	//fprintf( stderr, "det(K) = %10.2lf\n", detK);
	
	T** invK = allocMatrix<T>( n, n);
	//inverseCholesky( K, invK, n);
	inverseLU( K, invK, n);
	//invertGaussJordan( K, invK, n);
	//fprintf( stderr, "Inverse of Covariance Matrix: K^-1 =\n");
	//printMatrix( invK, n, n);
	
	T tmp[n];
	T fmm[n], *pfmm = fmm;
	T* ptmp = tmp;
	for(int i=0; i<n; i++)
		fmm[i] = f[i] -  theta[3];

	vectorMatrixProduct( pfmm, invK, ptmp, n, n);
	lik = - 0.5 * dotProduct( ptmp, pfmm, n)
	      - 0.5 * log( detK)
		  - 0.5 * n * log(2*M_PI// TODO: *pow(10, theta[0])
	);
    //- 0.5 * d * log(2*M_PI);
	
	freeMatrix( K,    n);
	freeMatrix( invK, n);

	return lik;
}

//#include "mex.h"

template <class T> //  >>> M E X <<<
T log_marginal_likelihood(
		T* &f, T* &s2, T* &X, int &n, int &d,
		T (*covFunc) ( T* &, T* &, int &, T* &),
		T* &theta)
{
	T lik = 0;

	T **K = allocMatrix<T>( n, n);
	covMatrix<T>( K, X, n,d, covFunc, theta);
	//fprintf( stderr, "Covariance Matrix: K =\n");
	//printMatrix( K, n, n);

	// noise
	if( s2==0)
		for(int i=0; i<n; i++)
			K[i][i] += pow(10, theta[0]);
	else
		for(int i=0; i<n; i++)
			K[i][i] += s2[i] + pow(10, theta[0]);

	T detK;
	//T   detK = determinantCholesky( K, n);
	//T   detK = determinantLU( K, n);

	T** invK = allocMatrix<T>( n, n);
	//inverseCholesky( K, invK, n);
	//inverseLU( K, invK, n);
	//invertGaussJordan( K, invK, n);
	//fprintf( stderr, "Inverse of Covariance Matrix: K^-1 =\n");
	//mexPrintMatrix( invK, n, n);
	//printMatrix( invK, n, n);

	gsl_error_handler_t *old_handler = gsl_set_error_handler_off();

	int error=determinantAndInverseLU( K, invK, detK, n);
	if( error ){
		//mexPrintf ("error: %s (%i)\n", gsl_strerror (error), error);
		return nan("");
	}
		//mexPrintf ("info: matrix inversion succeeded\n");

	/* restore original handler */
	gsl_set_error_handler( old_handler);

	//inverseCholesky( K, invK, n);
	//determinantCholesky( K, n);

	// mean function
	double fmm[n], *pfmm = fmm;
	for(int i=0; i<n; i++)
		fmm[i] = f[i] -  theta[3];

	T tmp[n];
	T* ptmp = tmp;
	vectorMatrixProduct( pfmm, invK, ptmp, n, n);
	lik = - 0.5 * dotProduct( ptmp, pfmm, n)
	      - 0.5 * log( detK)
		  - 0.5 * n * log(2*M_PI
			//* pow(10, theta[0])
		  );
    //- 0.5 * d * log(2*M_PI);

	freeMatrix( K,    n);
	freeMatrix( invK, n);

	return lik;
}

#define LIK_IPP
#endif
