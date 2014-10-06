#include "mex.h"




#include "gp.ipp"
#include "lik.ipp"
#include "cov.ipp"

void
mexFunction(int nlhs,mxArray *plhs[],int nrhs,const mxArray *prhs[])
{
	// READ TRAINING POINTS
	int d = (int) mxGetM(prhs[0]);
	int n = (int) mxGetN(prhs[0]);
	double *X = mxGetPr(prhs[0]);
	double *f = mxGetPr(prhs[1]);
	double *s2 = 0;
	double *ptheta = 0;
	if( nrhs==4){
		s2     = mxGetPr(prhs[2]);
		ptheta = mxGetPr(prhs[3]);

	}else
		ptheta = mxGetPr(prhs[2]);

	double (*covFunc) ( double* &, double* &, int &, double* &) = covMatern5;

	//mexPrintf( "alloc\n");
	plhs[0] = mxCreateDoubleMatrix(1,1, mxREAL);
	double *lik = mxGetPr(plhs[0]);
	//mexPrintf( "mean = %e\n", ptheta[3]);

	lik[0] = log_marginal_likelihood<double>( f, s2, X, n, d, covFunc, ptheta);

	if(isnan(lik[0]) || isinf(lik[0])){
		//mexPrintf( "WARNING: likelihood = %e\n", lik[0] );
		//mexPrintf( "WARNING: likelihood = %e\n -> set to %e", lik[0], -FLT_MAX );
		lik[0] = -FLT_MAX;
		//lik[0] = -DBL_MAX;
		//lik[0] = - 1.0 / 0.0;
		//lik[0] = 0x7F800000;
	}

}
