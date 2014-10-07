#include "mex.h"




#include "gp.ipp"
#include "lik.ipp"
#include "cov.ipp"

void minimizeSampling(
		double* &f, double* &X, int &n, int &d,
		double*&theta, double*&lb, double*&ub, double dtheta,
		double (*covFunc) ( double* &, double* &, int &, double* &)){

	double maxtheta[3], *ptheta = theta;
	double maxlik = -DBL_MAX;

	double* s2 = 0; // TODO

	int count = 0;

	for( theta[0]=lb[0]; theta[0]<=ub[0]; theta[0]+=dtheta)
	for( theta[1]=lb[1]; theta[1]<=ub[1]; theta[1]+=dtheta)
	for( theta[2]=lb[2]; theta[2]<=ub[2]; theta[2]+=dtheta)
	for( theta[3]=lb[3]; theta[3]<=ub[3]; theta[3]+=dtheta)
	{
		double lik = log_marginal_likelihood<double>( f, s2, X, n, d, covFunc, ptheta);
		if( maxlik < lik){
			maxlik = lik;
			maxtheta[0] = theta[0];
			maxtheta[1] = theta[1];
			maxtheta[2] = theta[2];
		}
		count++;
	}

	//mexPrintf( "sampled %i points in hyper-parameter space\n", count);

	ptheta[0] = maxtheta[0];
	ptheta[1] = maxtheta[1];
	ptheta[2] = maxtheta[2];
}

void minimizeSteepestDecent(
		double* &f, double* &X, int &n, int &d,
		double*&theta, double dtheta,
		double (*covFunc) ( double* &, double* &, int &, double* &)){

	double* s2 = 0; // TODO

	double maxtheta[3] = {theta[0],theta[1],theta[2]}, *ptheta = maxtheta;
	double maxlik = log_marginal_likelihood<double>( f, s2, X, n, d, covFunc, ptheta);

	ptheta = theta;

	int count = 0;

	bool change = true;
	int i=0;
	for(; i<200 && change; i++){
		change = false;
		for(int j=0; j<3; j++){
			ptheta[j] = maxtheta[j]+dtheta;
			double lik = log_marginal_likelihood<double>( f, s2, X, n, d, covFunc, ptheta); count++;
			//mexPrintf( "%i: maxlik=%e lik=%e\n", j, maxlik, lik);
			if( maxlik < lik){
				// go step forth
				//mexPrintf( "%i: maxlik=%e lik=%e => go step forth\n", j, maxlik, lik);
				maxlik = lik;
				maxtheta[j] = ptheta[j];
				change = true;
			}else{
				ptheta[j] = maxtheta[j]-dtheta;
				lik = log_marginal_likelihood<double>( f, s2, X, n, d, covFunc, ptheta); count++;
				if( maxlik < lik){
					// go step back
					//mexPrintf( "%i: maxlik=%e lik=%e => go step back\n", j, maxlik, lik);
					maxlik = lik;
					maxtheta[j] = ptheta[j];
					change = true;
				}else{
					// do nothing
					//mexPrintf( "%i: do nothing\n", j, maxlik, lik);
					ptheta[j] = maxtheta[j];
				}
			}
		}

		if( !change && dtheta > 0.01){
			dtheta /= 2;
			change = true;
		}
	}

	//mexPrintf( "sampled %i points in hyper-parameter space\n", count);

	ptheta[0] = maxtheta[0];
	ptheta[1] = maxtheta[1];
	ptheta[2] = maxtheta[2];
}

//extern double (*covFunc) ( double* &, double* &, int &, double* &);
//double (*covMatern5) ( double* &, double* &, int &, double* &);

void
mexFunction(int nlhs,mxArray *plhs[],int nrhs,const mxArray *prhs[])
{
	/*if(nrhs==2){
		// input & output arguments for function handle
		mxArray *rhs[2];
		mxArray *lhs[1];

		rhs[0] = const_cast<mxArray *>(prhs[0]);
		rhs[1] = (mxArray*) prhs[1];

		mexCallMATLAB(0, lhs, 1, rhs, "feval");

		double (*covFuncMatern5) ( double* &, double* &, int &, double* &) = covMatern5;
		if( covFunc==covFuncMatern5)
			mexPrintf( "covMatern5\n");
		else
			mexPrintf( "covFunc = %i\n", covFunc);

		//plhs[0] = mxCreateDoubleScalar(*mxGetPr(lhs[0]));

		return;
	}*/

	// READ TRAINING POINTS
	int d = (int) mxGetM(prhs[0]);
	int n = (int) mxGetN(prhs[0]);
	double *X = mxGetPr(prhs[0]);
	double *Xm[n];
	double *f = mxGetPr(prhs[1]);

	/*mexPrintf( "%i points with dimension %i\n", n, d);
	for( int i=0; i<n; i++){
		Xm[i] = &X[i*n];
		for( int j=0; j<d; j++)
			mexPrintf( "%10.4lf ", X[i*d+j]);
		mexPrintf( "| %10.4lf\n", f[i]);
	}*/



	// TRAIN HYPERPARAMETERS
	double (*covFunc) ( double* &, double* &, int &, double* &) = covMatern5;
	//double (*covFunc) ( double* &, double* &, int &, double* &) = covSquaredExponential;
	//fclose(fp);

	double theta[4] = {-6., 0, -6., 0}, *ptheta = theta;
	double lb[4] = {-6, -6, -6., 0};
	double ub[4] = {2, 2, 2, 0};
	double dtheta = 1;
	double *plb=lb, *pub=ub;
	switch( nrhs){
	case 5:
		dtheta = mxGetScalar(prhs[4]);
	case 4:
		plb=mxGetPr(prhs[2]);
		pub=mxGetPr(prhs[3]);
	case 2:
		// INITIAL SAMPLING
		//mexPrintf( "initial sampling of hyper parameters\n");
		minimizeSampling( f, X, n, d, ptheta, plb, pub, dtheta, covFunc);
		break;
	case 3:
		//mexPrintf( "use predefined hyper parameters\n");
		ptheta = mxGetPr(prhs[2]);
		minimizeSteepestDecent( f, X, n, d, ptheta, dtheta, covFunc);
		break;
	}

	//
	//return;

	double *s2 = 0; // TODO
	double maxlik = log_marginal_likelihood<double>( f, s2, X, n, d, covFunc, ptheta);
	/*mexPrintf( "maxtheta = \n", ptheta[0],ptheta[1],ptheta[2], maxlik);
	for( int i=0; i<4; i++)
		mexPrintf( "%10.2e ", ptheta[i]);
	mexPrintf( " maxlik = %10.2e\n", maxlik);
	 */

	// OUTPUT HYPER PARAMETERS
	plhs[0] = mxCreateDoubleMatrix(1, 4, mxREAL);
	double *maxtheta = mxGetPr(plhs[0]);
	maxtheta[0] = ptheta[0];
	maxtheta[1] = ptheta[1];
	maxtheta[2] = ptheta[2];
	maxtheta[3] = ptheta[3];

	if( nlhs==2){
		plhs[1] = mxCreateDoubleMatrix(1,1, mxREAL);

		double *lik = mxGetPr(plhs[1]);
		*lik = maxlik;
	}
}
