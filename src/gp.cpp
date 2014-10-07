#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <float.h>

#include <gsl/gsl_sf_bessel.h>

#include "matrix.hpp"
#include "cov.hpp"
#include "lik.hpp"
#include "gp.hpp"

double rand01( double *x, int n)
{ return rand() / (double)RAND_MAX; }

double sum_of_squares( double *x, int n)
{ 
	double ret = 0;
	for( int j=0; j<n; j++)
		ret += (x[j]-0.5)*(x[j]-0.5);
	return ret; 
}

double sum_of_sinus( double *x, int n)
{
	double ret = 0;
	for( int j=0; j<n; j++)
		ret += sin(x[j]*(M_PI*2.));
	return ret;
}

template <class T>
void eval( 
	double** &x,    double* &f,    int &n,
	double** &xnew, double* &fnew, int &nnew,
	int &d,
	T (*covFunc) ( T* &, T* &, int &, T* &), double* &theta)
{
	// covariance
	double **K = allocMatrix<double>( n, n);
	covMatrix( K, x, n,d, covFunc, theta);
	
	// inverse of K
	double** invK = allocMatrix<double>( n, n);
	inverseCholesky( K, invK, n);
	
	double tmp[n], *ptmp = tmp;
	double cov[n], *pcov = cov;
	
	for( int i=0; i<nnew; i++){
		// cov with new point
		covVector( pcov, xnew[i], x, n, d, covFunc, theta);
		vectorMatrixProduct( pcov, invK, ptmp, n, n);
	
		// mean
		fnew[i] = dotProduct<double>( ptmp, f, n);	
	}
}


/*int main()
{
	int n=40;
	int d=2;
//	double (*covFunc) ( double* &, double* &, int &, double* &) = covSquaredExponential;//covMatern5;
	double (*covFunc) ( double* &, double* &, int &, double* &) = covMatern5;

	double **x = allocMatrix<double>( n, d);
	double *f  = (double*) malloc( n*sizeof(double));
	
	double (*pFunc) (double*, int) = sum_of_sinus; //sum_of_squares;
	
	// INIT
	fprintf( stderr, "Input Data: x | f\n");
	FILE *fp = fopen( "input.dat", "w+");
	for( int i=0; i<n; i++){
		for( int j=0; j<d; j++){
			x[i][j] = rand() / (double)RAND_MAX;
			fprintf( fp, "%10.2e ", x[i][j]);
		}
		double *px = x[i];
		f[i] = (*pFunc)(px, d)+rand01(0,0)*1;
		fprintf( fp, "  %10.2e\n", f[i]);
	}
	fclose(fp);
	

	double theta[3] = { 10,3,1.}, *ptheta = theta, maxtheta[3];

	//printMatrix( K, n,n);
	
	double lik = log_marginal_likelihood( f, x, n, d, covFunc, ptheta);
	double maxlik = -DBL_MAX;
	fprintf( stderr, "likelihood = %10.2lf\n", lik);
	
	fp = fopen( "likelihood.dat", "w+");
	//for( theta[0]=0.01; theta[0]<50; theta[0]+=0.1)
	theta[0]=-6;
	theta[0]=-6;
	for( theta[0]=-10.; theta[0]<0; theta[0]+=0.1){
	for( theta[1]=-2.; theta[1]<2; theta[1]+=0.1)
	//for( theta[2]=-2; theta[2]<2; theta[2]+=0.01)
	{
			//fprintf( stderr, "%10.2e %10.2e\n", theta[0],theta[1]);
			lik = log_marginal_likelihood( f, x, n, d, covFunc, ptheta);
			if( maxlik < lik){
				maxlik = lik;
				maxtheta[0] = theta[0];
				maxtheta[1] = theta[1];
				maxtheta[2] = theta[2];
			}
			fprintf( fp, "%10.2e %10.2e %10.2e %10.2e\n", theta[0],theta[1],theta[2],lik);
		}
		fprintf( fp, "\n");
	}
	fclose(fp);

	fprintf( stderr, "maxtheta = %10.2e %10.2e %10.2e maxlik = %10.2e\n", maxtheta[0],maxtheta[1],maxtheta[2], maxlik);
	theta[0] = maxtheta[0];
	theta[1] = maxtheta[1];
	theta[2] = maxtheta[2];

	int nnew = 1000;
	/*double **xnew = allocMatrix<double>(nnew,d);
	for( int i=0; i<nnew; i++)
		for( int j=0; j<d; j++)
			xnew[i][j] = rand() / (double)RAND_MAX;
	 * /

	int nx=100, ny=(d==1?1:100); nnew = nx*ny;
	double **xnew = allocMatrix<double>(nnew,d);
	switch(d){
	case 1:
		for( int ix=0; ix<nx; ix++)
			xnew[ix][0] = ix/(double)nx;

		break;
	case 2:
		for( int ix=0; ix<nx; ix++)
		for( int iy=0; iy<ny; iy++){
			xnew[ix+iy*nx][0] = ix/(double)nx;
			xnew[ix+iy*nx][1] = iy/(double)ny;
		}
	}
	
	double *fnew   = (double*) malloc( nnew*sizeof(double));
	double *s2new  = (double*) malloc( nnew*sizeof(double));
	
	/*eval(
		x,    f,    n,
		xnew, fnew, nnew,
		d,
		covFunc, ptheta);* /
	evalVariance( // TODO: integrate uncertainty s2
		x,    f,     0,      n,
		xnew, fnew, s2new, nnew,
		d,
		covFunc, ptheta);
	
	fp = fopen( "output.dat", "w+");
	for (int ix = 0; ix < nx; ix++) {
		for (int iy = 0; iy < ny; iy++)
		//for( int i=0; i<nnew; i++)
		{
			int i = ix + iy * nx;
			for (int j = 0; j < d; j++)
				fprintf(fp, "%10.2e ", xnew[i][j]);
			fprintf(fp, "  %10.2e  %10.2e\n", fnew[i], s2new[i]);
		}
		if (d == 2)
			fprintf(fp, "\n");
	}
	fclose(fp);
	return 0;	
}*/

#ifdef MATLAB

#include "mex.h"

void
mexFunction(int nlhs,mxArray *plhs[],int nrhs,const mxArray *prhs[])
{
	// TRAINING POINTS
	int d = (int) mxGetM(prhs[0]);
	int n = (int) mxGetN(prhs[0]);
	double *X = mxGetPr(prhs[0]);
	double *Xm[n];
	double *f = mxGetPr(prhs[1]);
	double *s2 = 0;

	mexPrintf( "%i points with dimension %i\n", n, d);
	for( int i=0; i<n; i++){
		Xm[i] = &X[i*n];
		for( int j=0; j<d; j++)
			mexPrintf( "%10.4lf ", X[i*d+j]);
		mexPrintf( "| %10.4lf\n", f[i]);
	}

	// TEST POINTS
	int ntest = (int) mxGetN(prhs[2]);
	double *Xtest= mxGetPr(prhs[2]);
	double *Xmtest[ntest];

	mexPrintf( "%i points with dimension %i\n", ntest, d);
	for( int i=0; i<ntest; i++){
		Xmtest[i] = &Xtest[i*d];
		for( int j=0; j<d; j++)
			mexPrintf( "%10.4lf ", Xtest[i*d+j]);
		mexPrintf( "\n");
	}

	// TRAIN HYPERPARAMETERS
	//double (*covFunc) ( double* &, double* &, int &, double* &) = covMatern5;
	double (*covFunc) ( double* &, double* &, int &, double* &) = covSquaredExponential;
	double theta[3], maxtheta[3], *ptheta = theta;
	double maxlik = -DBL_MAX;
	theta[0]=-6;
	for( theta[0]=-6.; theta[0]<0; theta[0]+=1)
	for( theta[1]=-2.; theta[1]<2; theta[1]+=0.1){
	//for( theta[2]=-2; theta[2]<2; theta[2]+=0.01)
	{
			//fprintf( stderr, "%10.2e %10.2e\n", theta[0],theta[1]);
			double lik = log_marginal_likelihood<double>( f, s2, X, n, d, covFunc, ptheta);
			if( maxlik < lik){
				maxlik = lik;
				maxtheta[0] = theta[0];
				maxtheta[1] = theta[1];
				maxtheta[2] = theta[2];
			}
			//fprintf( fp, "%10.2e %10.2e %10.2e %10.2e\n", theta[0],theta[1],theta[2],lik);
		}
		//fprintf( fp, "\n");
	}
	//fclose(fp);

	mexPrintf( "maxtheta = %10.2e %10.2e %10.2e maxlik = %10.2e\n", maxtheta[0],maxtheta[1],maxtheta[2], maxlik);
	theta[0] = maxtheta[0];
	theta[1] = maxtheta[1];
	theta[2] = maxtheta[2];

	// OUTPUT
	plhs[0] = mxCreateDoubleMatrix(1, ntest, mxREAL); double *ftest = mxGetPr(plhs[0]);
	plhs[1] = mxCreateDoubleMatrix(1, ntest, mxREAL); double *s2test = mxGetPr(plhs[1]);
	evalVariance(
			X,     f,     s2,     n,
			Xtest, ftest, s2test, ntest,
			d,
			covFunc, ptheta);
}
#endif
