/*
 * matrix.cpp
 *
 *  Created on: 04.11.2013
 *      Author: jagiella
 */

#include "mex.h"
#include "matrix.ipp"

//#include "mex.h"

template <class T>
void mexPrintMatrix( T** &A, int n, int m)
{
	for( int i=0; i<n; i++){
		for( int j=0; j<m; j++){
			mexPrintf( "%10.2lf ", A[i][j]);
		}
		mexPrintf( "\n");
	}
}


void
mexFunction(int nlhs,mxArray *plhs[],int nrhs,const mxArray *prhs[])
{
	// READ TRAINING POINTS
	int m = (int) mxGetM(prhs[0]);
	int n = (int) mxGetN(prhs[0]);
	double *A = mxGetPr(prhs[0]);

	plhs[0] = mxCreateDoubleMatrix(m,n, mxREAL);
	double *invA = mxGetPr(plhs[0]);

	invertGaussJordan( A, invA, n);
}

