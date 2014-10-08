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
	double **K = allocMatrix<double>( n, n);
	covMatrix<double>( K, x, n,d, covFunc, theta);

	// noise
	if(s2==0)
		for(int i=0; i<n; i++)
			K[i][i] += pow(10,theta[0]);
	else
		for(int i=0; i<n; i++)
			K[i][i] += s2[i] + pow(10,theta[0]);

	// inverse of K
	double** invK = allocMatrix<double>( n, n);
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

	freeMatrix( K,    n);
	freeMatrix( invK, n);
}

template <class T> // 				>>> M E X <<<
void evalVariance(
	T* &x,    T* &f,    T* &s2,    int &n,
	T* &xnew, T* &fnew, T* &s2new, int &nnew,
	int &d,
	T (*covFunc) ( T* &, T* &, int &, T* &), double* &theta)
{
	// covariance
	T **K = allocMatrix<T>( n, n);
	covMatrix<T>( K, x, n,d, covFunc, theta);

	// noise
	if(s2==0)
		for(int i=0; i<n; i++)
			K[i][i] += pow(10,theta[0]);
	else
		for(int i=0; i<n; i++)
			K[i][i] += s2[i] + pow(10,theta[0]);

	// inverse of K
	T** invK = allocMatrix<T>( n, n);
	//inverseCholesky( K, invK, n);

	gsl_error_handler_t *old_handler = gsl_set_error_handler_off();
	int error = inverseLU( K, invK, n);
	if( error ){
		//mexPrintf ("error: %s (%i)\n", gsl_strerror (error), error);
	}
	gsl_set_error_handler( old_handler);


	//invertGaussJordan( K, invK, n);

	T tmp[n], *ptmp = tmp;
	T cov[n], *pcov = cov;

	// mean function
	T fmm[n], *pfmm = fmm;
	for(int i=0; i<n; i++)
		fmm[i] = f[i] - theta[3];


	for( int i=0; i<nnew; i++){
		// cov with new point
		T* pRow = &xnew[i*d];
		covVector( pcov, pRow, x, n, d, covFunc, theta);
		vectorMatrixProduct( pcov, invK, ptmp, n, n);

		// mean
		fnew[i] = dotProduct<T>( ptmp, pfmm, n) + theta[3];

		// std deviation
		s2new[i] = covFunc( pRow, pRow, d, theta) - dotProduct<T>( ptmp, pcov, n);
	}

	freeMatrix<T>( K,    n);
	freeMatrix<T>( invK, n);
}

#include "SparseMatrix.hpp"
#include "Solver.hpp"
template <class T> // 				>>> S P A R S E <<<
void evalVarianceSparse(
	T* &x,    T* &f,    T* &s2,    int &n,
	T* &xnew, T* &fnew, T* &s2new, int &nnew,
	int &d,
	T (*covFunc) ( T* &, T* &, int &, T* &), double* &theta)
{
	// Solver
	Solver<double> *S = new Solver<double>( n, Solver<double>::BiCGSTABl, 1e-30, 30);
	S->setDegree( 8);

	// covariance
	T **Kfull = allocMatrix<T>( n, n);
	covMatrix<T>( Kfull, x, n,d, covFunc, theta);
	if(s2==0)
		for(int i=0; i<n; i++)
			Kfull[i][i] += pow(10,theta[0]);
	else
		for(int i=0; i<n; i++)
			Kfull[i][i] += s2[i] + pow(10,theta[0]);

	// factorization
	T factorize[n];
	for(int i=0; i<n; i++)
		factorize[i] = 1./Kfull[i][i];

	T tmp[n], *ptmp = tmp;
	T cov[n], *pcov = cov;

	// mean function
	T fmm[n], *pfmm = fmm;
	for(int i=0; i<n; i++)
		fmm[i] = f[i] - theta[3];

	// Build System
	SparseMatrix<T> *K = new SparseMatrix<T>( n, n);
	for( int i=0; i<n; i++){
		for( int j=i; j<n; j++){
			double value = Kfull[i][j] * factorize[i];
			//if( fabs(value) > 1e-3)
				K->setLast(i,j, value);
			value = Kfull[i][j] * factorize[j];
			//if( fabs(value) > 1e-3)
				K->setLast(j,i, value);
		}
		fmm[i] *= factorize[i];
		tmp[i] = 0;
	}
	fprintf( stderr, "Sparsity = %f\n", K->countNonZero() / (double)( n*n));

	//K->printMatrix( "Kprec", "%4.2e ");
	//printVector( fmm, n, "y", "%4.2e ");

	// Solve sub-system: K^-1 * y
	//S->PreconditionJacobi(K, fmm);
	S->solve( K, fmm, ptmp);
	//printVector( ptmp, n, "K^-1 * y", "%4.2e ");

	fprintf( stderr, "Solved!\n");
	for( int i=0; i<nnew; i++){

		// cov with new point
		T* pXnewi = &xnew[i*d];
		covVector( pcov, pXnewi, x, n, d, covFunc, theta);

		// mean
		fnew[i] = dotProduct<T>( pcov, ptmp, n) + theta[3];

		// std deviation
		// construct system
		/*T k[n], *pk = k;
		T tmp2[n],*ptmp2= tmp2;
		for( int i=0; i<n; i++){
			pk[i] = pcov[i] * factorize[i];
			tmp2[i] = 0;
		}
		// Solve sub-system: K^-1 * k
		S->solve( K, pk, ptmp2);

		s2new[i] = covFunc( pXnewi, pXnewi, d, theta) - dotProduct<T>( pcov, ptmp2, n);;
*/
		s2new[i] = 0;
	}

	delete K;
	delete S;

}
template <class T> // 				>>> S P A R S E <<<
void evalVarianceSparse(
	T** x,    T* f,    T* s2,    int n,
	T** xnew, T* fnew, T* s2new, int nnew,
	int d,
	T (*covFunc) ( T* &, T* &, int &, T* &), double* theta)
{
	// Solver
	Solver<double> *S = new Solver<double>( n, Solver<double>::BiCGSTABl, 1e-30, 30);
	S->setDegree( 8);

	// covariance
	T **Kfull = allocMatrix<T>( n, n);
	covMatrix<T>( Kfull, x, n,d, covFunc, theta);
	if(s2==0)
		for(int i=0; i<n; i++)
			Kfull[i][i] += pow(10,theta[0]);
	else
		for(int i=0; i<n; i++)
			Kfull[i][i] += s2[i] + pow(10,theta[0]);

	// factorization
	T factorize[n];
	for(int i=0; i<n; i++)
		factorize[i] = 1./Kfull[i][i];

	T tmp[n], *ptmp = tmp;
	T cov[n], *pcov = cov;

	// mean function
	T fmm[n], *pfmm = fmm;
	for(int i=0; i<n; i++)
		fmm[i] = f[i] - theta[3];

	// Build System
	SparseMatrix<T> *K = new SparseMatrix<T>( n, n);
	for( int i=0; i<n; i++){
		for( int j=i; j<n; j++){
			double value = Kfull[i][j] * factorize[i];
			//if( fabs(value) > 1e-3)
				K->setLast(i,j, value);
			value = Kfull[i][j] * factorize[j];
			//if( fabs(value) > 1e-3)
				K->setLast(j,i, value);
		}
		fmm[i] *= factorize[i];
		tmp[i] = 0;
	}
	fprintf( stderr, "Sparsity = %f\n", K->countNonZero() / (double)( n*n));

	//K->printMatrix( "Kprec", "%4.2e ");
	//printVector( fmm, n, "y", "%4.2e ");

	// Solve sub-system: K^-1 * y
	//S->PreconditionJacobi(K, fmm);
	S->solve( K, fmm, ptmp);
	//printVector( ptmp, n, "K^-1 * y", "%4.2e ");

	fprintf( stderr, "Solved!\n");
	for( int i=0; i<nnew; i++){

		// cov with new point
		T* pXnewi = xnew[i];
		covVector( pcov, pXnewi, x, n, d, covFunc, theta);

		// mean
		fnew[i] = dotProduct<T>( pcov, ptmp, n) + theta[3];

		// std deviation
		// construct system
		/*T k[n], *pk = k;
		T tmp2[n],*ptmp2= tmp2;
		for( int i=0; i<n; i++){
			pk[i] = pcov[i] * factorize[i];
			tmp2[i] = 0;
		}
		// Solve sub-system: K^-1 * k
		S->solve( K, pk, ptmp2);

		s2new[i] = covFunc( pXnewi, pXnewi, d, theta) - dotProduct<T>( pcov, ptmp2, n);;
*/
		//s2new[i] = 0;
	}

	delete K;
	delete S;

}
template <class T> // 				>>> S P A R S E <<<
void evalGradientSparse(
	T** x,    T* f,    T* s2,    int n,
	T** xnew, T** dfnew, int nnew,
	int d,
	T (*covFunc) ( T* &, T* &, int &, T* &), double* theta)
{
	// Solver
	Solver<double> *S = new Solver<double>( n, Solver<double>::BiCGSTABl, 1e-30, 30);
	S->setDegree( 8);

	// covariance
	T **Kfull = allocMatrix<T>( n, n);
	covMatrix<T>( Kfull, x, n,d, covFunc, theta);
	if(s2==0)
		for(int i=0; i<n; i++)
			Kfull[i][i] += pow(10,theta[0]);
	else
		for(int i=0; i<n; i++)
			Kfull[i][i] += s2[i] + pow(10,theta[0]);

	// factorization
	T factorize[n];
	for(int i=0; i<n; i++)
		factorize[i] = 1./Kfull[i][i];

	T tmp[n], *ptmp = tmp;
	T cov[n], *pcov = cov;

	// mean function
	T fmm[n], *pfmm = fmm;
	for(int i=0; i<n; i++)
		fmm[i] = f[i] - theta[3];

	// Build System
	SparseMatrix<T> *K = new SparseMatrix<T>( n, n);
	for( int i=0; i<n; i++){
		for( int j=i; j<n; j++){
			double value = Kfull[i][j] * factorize[i];
			//if( fabs(value) > 1e-3)
				K->setLast(i,j, value);
			value = Kfull[i][j] * factorize[j];
			//if( fabs(value) > 1e-3)
				K->setLast(j,i, value);
		}
		fmm[i] *= factorize[i];
		tmp[i] = 0;
	}
	fprintf( stderr, "Sparsity = %f\n", K->countNonZero() / (double)( n*n));

	//K->printMatrix( "Kprec", "%4.2e ");
	//printVector( fmm, n, "y", "%4.2e ");

	// Solve sub-system: K^-1 * y
	//S->PreconditionJacobi(K, fmm);
	S->solve( K, fmm, ptmp);
	//printVector( ptmp, n, "K^-1 * y", "%4.2e ");

	fprintf( stderr, "Solved!\n");
	for( int i=0; i<nnew; i++){

		// cov with new point
		//T* pXnewi = &xnew[i*d];
		T* pXnewi = xnew[i];
		covVector( pcov, pXnewi, x, n, d, covFunc, theta);

		// mean
		double fnew_i = dotProduct<T>( pcov, ptmp, n);

		for( int di=0; di<d; di++)
		{
			double dx = 1e-6;
			pXnewi[di] += dx;
			covVector( pcov, pXnewi, x, n, d, covFunc, theta);
			pXnewi[di] -= dx;

			dfnew[i][di] = (dotProduct<T>( pcov, ptmp, n) - fnew_i) / dx;
		}
	}

	delete K;
	delete S;

}


#define GP_IPP
#endif
