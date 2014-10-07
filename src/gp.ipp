#ifndef GP_IPP

#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "matrix.ipp"
#include "cov.ipp"
//#include "mex.h"

template <class T>
void evalVariance(
	T** &x,    T* &f,    T* &s2,    int &n,
	T** &xnew, T* fnew, T* s2new, int &nnew,
	int &d,
	T (*covFunc) ( T* &, T* &, int &, T* &), double* &theta)
{
	// covariance
	double **K = newMatrix<double>( n, n);
	covMatrix<double>( K, x, n,d, covFunc, theta);

	// noise
	if(s2==0)
		for(int i=0; i<n; i++)
			K[i][i] += pow(10,theta[0]);
	else
		for(int i=0; i<n; i++)
			K[i][i] += s2[i] + pow(10,theta[0]);

	// inverse of K
	double** invK = newMatrix<double>( n, n);
	//inverseCholesky( K, invK, n);
	inverseLU( K, invK, n);

	double tmp[n], *ptmp = tmp;
	double cov[n], *pcov = cov;

	// mean function
	double fmm[n], *pfmm = fmm;
	for(int i=0; i<n; i++)
		fmm[i] = f[i] - theta[3];

	for( int i=0; i<nnew; i++){
		// cov with new point
		covVector( pcov, xnew[i], x, n, d, covFunc, theta);
		vectorMatrixProduct( pcov, invK, ptmp, n, n);

		// mean
		fnew[i] = dotProduct<double>( ptmp, pfmm, n) + theta[3];

		// std deviation
		s2new[i] = covFunc( xnew[i], xnew[i], d, theta) - dotProduct<double>( ptmp, pcov, n);
		//fprintf(stderr, "%e - %e\n", fnew[i], s2new[i]);

	}

	deleteMatrix( K,    n);
	deleteMatrix( invK, n);
}

template <class T> // 				>>> M E X <<<
void evalVariance(
	double* &x,    double* &f,    double* &s2,    int &n,
	double* &xnew, double* &fnew, double* &s2new, int &nnew,
	int &d,
	T (*covFunc) ( T* &, T* &, int &, T* &), double* &theta)
{
	// covariance
	double **K = newMatrix<double>( n, n);
	covMatrix<double>( K, x, n,d, covFunc, theta);

	// noise
	if(s2==0)
		for(int i=0; i<n; i++)
			K[i][i] += pow(10,theta[0]);
	else
		for(int i=0; i<n; i++)
			K[i][i] += s2[i] + pow(10,theta[0]);

	// inverse of K
	double** invK = newMatrix<double>( n, n);
	//inverseCholesky( K, invK, n);

	gsl_error_handler_t *old_handler = gsl_set_error_handler_off();
	int error = inverseLU( K, invK, n);
	if( error ){
		//mexPrintf ("error: %s (%i)\n", gsl_strerror (error), error);
	}
	gsl_set_error_handler( old_handler);


	//invertGaussJordan( K, invK, n);

	double tmp[n], *ptmp = tmp;
	double cov[n], *pcov = cov;

	// mean function
	double fmm[n], *pfmm = fmm;
	for(int i=0; i<n; i++)
		fmm[i] = f[i] - theta[3];


	for( int i=0; i<nnew; i++){
		// cov with new point
		T* pRow = &xnew[i*d];
		covVector( pcov, pRow, x, n, d, covFunc, theta);
		vectorMatrixProduct( pcov, invK, ptmp, n, n);

		// mean
		fnew[i] = dotProduct<double>( ptmp, pfmm, n) + theta[3];

		// std deviation
		s2new[i] = covFunc( pRow, pRow, d, theta) - dotProduct<double>( ptmp, pcov, n);
	}

	deleteMatrix( K,    n);
	deleteMatrix( invK, n);
}


#define GP_IPP
#endif
