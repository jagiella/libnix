/*
 * optimization.cpp
 *
 *  Created on: Oct 7, 2014
 *      Author: jagiella
 */


#include <stdio.h>
#include <math.h>

#include "matrix.hpp"
#include "statistics.hpp"


void getGradient(
		int parameterSize, double *parameters, double *gradient,
		double (*f)( int, const double*, double*, void *), void *f_data,
		double dx)
{
	double yn = f( parameterSize, parameters, 0, f_data),
		   yn1;

	for( int i=0; i<parameterSize; i++){
		double delta = dx*fmax( fabs( parameters[i]), 1);

		parameters[i] += delta; yn1 = f( parameterSize, parameters, 0, f_data);
		gradient[i] = (yn1 - yn) / delta;
		parameters[i] -= delta;

	}
}
void getHessian(
		int parameterSize, double *parameters, double **H,
		double (*f)( int, const double*, double*, void *), void *f_data,
		double dx)
{
	double y = f( parameterSize, parameters, 0, f_data),
		   yi, yj, yij, yii;

	for( int i=0; i<parameterSize; i++){
		double deltai = dx*fmax( fabs( parameters[i]), 1);

		for( int j=i+1; j<parameterSize; j++){
			double deltaj = dx*fmax( fabs( parameters[j]), 1);

			parameters[i] += deltai;
			yi = f( parameterSize, parameters, 0, f_data);

			parameters[j] += deltaj;
			yij = f( parameterSize, parameters, 0, f_data);

			parameters[i] -= deltai;
			yj = f( parameterSize, parameters, 0, f_data);

			parameters[j] -= deltaj;
			H[i][j] = H[j][i] = (y - yi - yj + yij) / (deltai*deltaj);
		}

		parameters[i] += deltai;
		yi = f( parameterSize, parameters, 0, f_data);
		parameters[i] += deltai;
		yii = f( parameterSize, parameters, 0, f_data);
		parameters[i] -= 2*deltai;
		H[i][i] = (y - 2*yi + yii) / (deltai*deltai);

		//parameters[i] = xold;
	}
}


double LevenbergMarquardt(
	int parameterSize, double *parameters, double *parametersMin, double *parametersMax,
	double differentiationStepSize, double lambda, double minSquares, int maxIterations,
	double (*f)( int, const double*, double*, void *), void *f_data
	)
{
	double **H = allocMatrix<double>( parameterSize, parameterSize);
	double **A = allocMatrix<double>( parameterSize, parameterSize);
	double b[parameterSize];
	double dparameters[parameterSize];
	double gradient[parameterSize];
	double initialGuess[parameterSize];
	for (int i = 0; i < parameterSize; i++)
		initialGuess[i] = parameters[i];

	double newS = (*f)(parameterSize, parameters, 0, f_data),
		  lastS = newS+minSquares;
	bool stop=false;

	for( int iter=0; iter<maxIterations && !stop; iter++){
		lambda=0;
		//fprintf(stderr, "S(%i) = %e\n", lastS, iter);
		//fprintf(stderr, "par\n"); printVector( dparameters, parameterSize, "%5.2e ");

		// get gradient
		getGradient(parameterSize, parameters, gradient,
					f, f_data,
					differentiationStepSize);
		//fprintf( stderr, "grad: %.2e, %.2e\n", gradient[0], gradient[1]);

		// get Hessian
		getHessian(	parameterSize, parameters, H,
					f, f_data,
					differentiationStepSize);
		//printMatrix( H, parameterSize, parameterSize);
		//fprintf(stderr, "H\n"); printMatrix( H, parameterSize, parameterSize);
		//fprintf(stderr, "grad\n"); printVector( gradient, parameterSize, "%5.2e ");

		// Construct Linear System: H*dx = grad
		for (int i = 0; i < parameterSize; i++) {

			// Levenberg
			for (int j = 0; j < parameterSize; j++)
				A[i][j] = H[i][j];
				//A[i][j] = gradient[i] * gradient[j];

			// Marquardt
			A[i][i] *= (1. + lambda);

			b[i] = -gradient[i];
		}

		// Normalize
		for (int i = 0; i < parameterSize; i++) if(A[i][i]!=0){
			double Aii = A[i][i];
			for (int j = 0; j < parameterSize; j++)
				A[i][j] /= Aii;
			b[i] /= Aii;
		}

		// Solve Linear System
		solveLinearSystem(A, b, dparameters, parameterSize);

		// Update Parameters
		for (int i = 0; i < parameterSize; i++)
			parameters[i] += dparameters[i];

		lastS = newS;
		newS = (*f)(parameterSize, parameters, 0, f_data);
		if( fabs(lastS-newS) < minSquares)
			stop = true;

		if( isnan( newS) || isinf( newS)){
			fprintf( stderr, "ERROR in LevenbergMarquardt: Score is %e\n", newS);
			printVector( dparameters, parameterSize, "%5.2e ");
			//exit( 0);
			lambda += 1;
			fprintf( stderr, "retray with bigger lamda = %f\n", lambda);
			for (int i = 0; i < parameterSize; i++)
						parameters[i] =initialGuess[i];
		}

		fprintf( stderr, "x: %.2e, %.2e, %.2e ->  %e\n", parameters[0], parameters[1], parameters[2], newS);
	}

	freeMatrix( H, parameterSize);
	freeMatrix( A, parameterSize);

	return newS;
}


double LineSearch(
	int parameterSize, double *parameters, double *parametersMin, double *parametersMax,
	double differentiationStepSize, double lambda, double minSquares, int maxIterations,
	double (*f)( int, const double*, double*, void *), void *f_data
	)
{
	double **H = allocMatrix<double>( parameterSize, parameterSize);
	double **A = allocMatrix<double>( parameterSize, parameterSize);
	double b[parameterSize];
	double dparameters[parameterSize];
	double gradient[parameterSize];
	double initialGuess[parameterSize];
	for (int i = 0; i < parameterSize; i++)
		initialGuess[i] = parameters[i];

	double newS = (*f)(parameterSize, parameters, 0, f_data),
		  lastS = newS+minSquares;
	bool stop=false;

	for( int iter=0; iter<maxIterations && !stop; iter++){
		// get gradient
		getGradient(parameterSize, parameters, gradient,
					f, f_data,
					differentiationStepSize);

		// Normalize Gradient
		double abs_gradient = 0;
		for( int i=0; i<parameterSize; i++){
			abs_gradient += pow( gradient[i], 2);
		}
		abs_gradient = sqrt( abs_gradient);
		for( int i=0; i<parameterSize; i++) if(abs_gradient>0)
			gradient[i] /= abs_gradient;
		fprintf( stderr, "grad: %.2e, %.2e\n", gradient[0], gradient[1]);

		// Follow Line as long as improving
		int pot=0;
		do{
			lastS = newS;
			double distance = differentiationStepSize*pow(2, pot);
			for (int i = 0; i < parameterSize; i++)	parameters[i] -= gradient[i]*distance;
			newS = (*f)(parameterSize, parameters, 0, f_data);
			fprintf( stderr, "x(%i): %.2e, %.2e, %.2e ->  %e\n", pot, parameters[0], parameters[1], parameters[2], newS);

			for (int i = 0; i < parameterSize; i++)	parameters[i] += gradient[i]*distance;
			pot++;
		}while( newS < lastS /*&& inbound(parameters, parametersMin, parametersMax, parameterSize)*/);
		pot--;
		do{
			pot--;
			lastS = newS;
			double distance = differentiationStepSize*pow(2, pot);
			for (int i = 0; i < parameterSize; i++)	parameters[i] -= gradient[i]*distance;
			newS = (*f)(parameterSize, parameters, 0, f_data);
			fprintf( stderr, "x(%i): %.2e, %.2e, %.2e ->  %e\n", pot, parameters[0], parameters[1], parameters[2], newS);

			for (int i = 0; i < parameterSize; i++)	parameters[i] += gradient[i]*distance;
		}while( newS < lastS /*&& inbound(parameters, parametersMin, parametersMax, parameterSize)*/);
		//pot++;

		// Update Parameters
		for (int i = 0; i < parameterSize; i++)
			parameters[i] -= gradient[i]*differentiationStepSize*pow(2, pot);

		fprintf( stderr, "x: %.2e, %.2e, %.2e ->  %e\n", parameters[0], parameters[1], parameters[2], newS);
	}

	freeMatrix( H, parameterSize);
	freeMatrix( A, parameterSize);

	return newS;
}
