#include "mex.h"


#include "gp.ipp"

void
mexFunction(int nlhs,mxArray *plhs[],int nrhs,const mxArray *prhs[])
{
	// READ TRAINING POINTS
	int d = (int) mxGetM(prhs[0]);
	int n = (int) mxGetN(prhs[0]);
	double *X = mxGetPr(prhs[0]);
	double *Xm[n];
	double *f = mxGetPr(prhs[1]);

	int i=0;
	double *s2 = 0;
	if( nrhs==5){
		s2 = mxGetPr(prhs[2]);
		i++;
		//mexPrintf( "USE UNCERTAINTY INFO\n");
	}

	// READ TEST POINTS
	int ntest = (int) mxGetN(prhs[i+2]);
	double *Xtest= mxGetPr(prhs[i+2]);
	double *Xmtest[ntest];

	// READ HYPERPARAMETERS
	double (*covFunc) ( double* &, double* &, int &, double* &) = covMatern5;
	//double (*covFunc) ( double* &, double* &, int &, double* &) = covSquaredExponential;

	double *ptheta = mxGetPr(prhs[i+3]);

	// OUTPUT PREDICTION OF TEST POINTS
	plhs[0] = mxCreateDoubleMatrix(1, ntest, mxREAL); double *ftest = mxGetPr(plhs[0]);
	plhs[1] = mxCreateDoubleMatrix(1, ntest, mxREAL); double *s2test = mxGetPr(plhs[1]);
	evalVariance(
			X,     f,     s2,     n,
			Xtest, ftest, s2test, ntest,
			d,
			covFunc, ptheta);
}
