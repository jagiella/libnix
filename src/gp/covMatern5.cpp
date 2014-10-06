#include "cov.ipp"
#include "matrix.ipp"
#include "mex.h"

extern double (*covFunc) ( double* &, double* &, int &, double* &);

void
mexFunction(int nlhs,mxArray *plhs[],int nrhs,const mxArray *prhs[])
{
	//mexPrintf( "Set covFunc to covMatern5\n");
	covFunc = covMatern5<double>;

	// READ TRAINING POINTS
	int d = (int) mxGetM(prhs[0]);
	int n = (int) mxGetN(prhs[0]);
	double *X = mxGetPr(prhs[0]);
	double *f = mxGetPr(prhs[1]);
	double *ptheta = mxGetPr(prhs[2]);


	double **K = newMatrix<double>( n, n);
	covMatrix<double>( K, X, n,d, covFunc, ptheta);

	plhs[0] = mxCreateDoubleMatrix(n, n, mxREAL);
	double *k = mxGetPr(plhs[0]);
	for( int i=0; i<n; i++)
		for( int j=0; j<n; j++)
			k[i*n+j] = K[i][j];



}
