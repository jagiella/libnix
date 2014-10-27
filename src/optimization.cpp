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
#include "optimization.hpp"


double getGradient(
		int parameterSize, double *parameters, double *gradient,
		double (*f)( int, const double*, double*, void *), void *f_data,
		double dx)
{
	double yn = f( parameterSize, parameters, 0, f_data),
		   yn1;
	//fprintf(stderr, "y(x) = %e\n", yn);
	for( int i=0; i<parameterSize; i++){
		double delta = dx*fmax( fabs( parameters[i]), 1);

		parameters[i] += delta; yn1 = f( parameterSize, parameters, 0, f_data);
		//fprintf(stderr, "y(%i) = %e\n", i,yn1);
		gradient[i] = (yn1 - yn) / delta;
		parameters[i] -= delta;

		if( isinf( gradient[i]) || isnan(gradient[i])){
			fprintf( stderr, "WARNING: gradient[%i] = %e\n", i, gradient[i]);
		}
	}

	return yn;
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

	double lastParameters[parameterSize];
	double newS = (*f)(parameterSize, parameters, 0, f_data),
		  lastS = newS+minSquares;
	bool stop=false;

	int iter=0;
	for( ; iter<maxIterations && !stop; iter++){
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

		// Construct Linear System: H*dx = -grad
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
		//fprintf( stderr, "x: ");
		for (int i = 0; i < parameterSize; i++){
			lastParameters[i]=parameters[i];
			parameters[i] += dparameters[i];
			parameters[i] = fmin( parameters[i], parametersMax[i]);
			parameters[i] = fmax( parameters[i], parametersMin[i]);
			//fprintf( stderr, "%.2e, ", parameters[i]);
		}
		//fprintf( stderr, " ->  %e\n", newS);

		lastS = newS;
		newS = (*f)(parameterSize, parameters, 0, f_data);
		if( fabs(lastS-newS) < minSquares)
			stop = true;

		if( isnan( newS) || isinf( newS)){
			fprintf( stderr, "ERROR in LevenbergMarquardt: Score is %e\n", newS);
			printVector( dparameters, parameterSize, "%5.2e ");
			//exit( 0);
			lambda += 1;
			fprintf( stderr, "retry with bigger lamda = %f\n", lambda);
			for (int i = 0; i < parameterSize; i++)
						parameters[i] =initialGuess[i];
		}


	}

	fprintf( stderr, "STOP because: lastS=%e, newS=%e, iter=%i\n", lastS, newS, iter);
	printVector<double>( gradient, parameterSize, "gradient", "%.2e");
	printMatrix(H, parameterSize, parameterSize, "Hessian", "%.2e");
	printVector<double>( parameters, parameterSize, "parameters", "%.2e");
	printVector<double>( lastParameters, parameterSize, "lastParameters", "%.2e");

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

#include <time.h>

void gradientDecent( double *x0, double *lb, double *ub, int dim, double (*func) (int, const double*, double*, void*), void* func_data, optimoptions *options, double *sol)
{
	//srand(2);
	bool stop = false;
	double X[dim], *lastX=sol, grad[dim],
	       Y,      lastY;
	double alpha = 1;

	unsigned int seed = clock();
	if(x0)
		for( int i=0; i<dim; i++)
			lastX[i] = x0[i];
	else
		for( int i=0; i<dim; i++){
			double beta = unifrnd(double, &seed);
			lastX[i] = ub[i]*beta + lb[i]*(1-beta); //(ub[i]+lb[i])/2.;
			//lastX[i] = ub[i]*7/23. + lb[i]*16/23.; //(ub[i]+lb[i])/2.;
		}
	if( options->Display)fprintf( stderr, "<< init %e >>\n", lastX[0]);
	int evals = 0;

	//alpha = norm_diff( lb, ub, dim);
	double alpha_min = 0;
	for( int i=0; i<dim; i++)
		alpha_min += fmax( pow( lastX[i], 2), 1);
	alpha_min = options->FinDiffRelStep * sqrt(alpha_min);


	for( int iter=0; iter<options->MaxIter && !stop; iter++)
	{
		if( options->Display){
			fprintf( stderr, "Iteration %i\n", iter);
			printVector( lastX, dim, "lastX", "%.2e ");
		}

		// GET NORMALIZED GRADIENT

		if( options->Display)fprintf( stderr, "Get Gradient\n");
		if( options->GradObj){
			Y = (*func) ( dim, lastX, grad, func_data); evals++;
		}else{
			Y = getGradient( dim, lastX, grad, func, func_data, options->FinDiffRelStep); evals+=dim+1;
		}
		if( options->Display)printVector( grad, dim, "grad", "%.2e ");

		// minimal step length
		alpha=0;
		for( int i=0; i<dim; i++)
			alpha += grad[i]*grad[i];
		alpha = sqrt( alpha);

		if( alpha == 0 || evals >= options->MaxFunEvals ){
			for( int i=0; i<dim; i++)
				sol[i] = lastX[i];
			return;
		}

		// normalize gradient
		for( int i=0; i<dim; i++)
			grad[i] /= alpha;
		alpha = alpha_min;


		// DECENT

		do{
			lastY = Y;

			// double step size
			alpha *= 2;
			for( int i=0; i<dim; i++){
				X[i] = lastX[i] - grad[i]*alpha;
				X[i] = fmin( X[i], ub[i]);
				X[i] = fmax( X[i], lb[i]);
			}

			// evaluate
			Y = (*func) ( dim, X, 0, func_data); evals++;

		}while(  Y < lastY);

		for( int i=0; i<dim; i++)
			lastX[i] = X[i];
		lastY = Y;


		if( options->Display)
			fprintf( stderr, "%i: %e <- %e\n", iter, Y, X[0]);

		//alpha *= 10;
	}
	fprintf( stderr, "Max. Iterations finished (alpha = %e)\n", alpha);
	//fprintf( stderr, "<< %e >>\n", lastX[0]); //exit(0);
}
void gradientDecentOLD( double *x0, double *lb, double *ub, int dim, double (*func) (int, const double*, double*, void*), void* func_data, optimoptions *options, double *sol)
{
	//srand(2);
	bool stop = false;
	double X[dim], *lastX=sol, grad[dim],
	       Y,      lastY;
	double alpha = 1;

	unsigned int seed = clock();
	if(x0)
		for( int i=0; i<dim; i++)
			lastX[i] = x0[i];
	else
		for( int i=0; i<dim; i++){
			double beta = unifrnd(double, &seed);
			lastX[i] = ub[i]*beta + lb[i]*(1-beta); //(ub[i]+lb[i])/2.;
			//lastX[i] = ub[i]*7/23. + lb[i]*16/23.; //(ub[i]+lb[i])/2.;
		}
	if( options->Display)fprintf( stderr, "<< init %e >>\n", lastX[0]);
	int evals = 0;

	alpha = norm_diff( lb, ub, dim);

	for( int iter=0; iter<options->MaxIter && !stop; iter++)
	{
		//alpha = uniform_norm( lb, ub, dim);

		// get gradient
		if( options->GradObj){
			lastY = (*func) ( dim, lastX, grad, func_data); evals++;
		}else{
			lastY = (*func) ( dim, lastX, 0,    func_data); evals++;
			for( int i=0; i<dim; i++)
			if(ub[i]-lb[i] != 0.){
				lastX[i] += options->FinDiffRelStep;
				grad[i] = ((*func) ( dim, lastX, 0, func_data) - lastY) / options->FinDiffRelStep; evals++;
				lastX[i] -= options->FinDiffRelStep;

				if( isinf(grad[i])){
					lastX[i] -= options->FinDiffRelStep;
					grad[i] = (lastY - (*func) ( dim, lastX, 0, func_data)) / options->FinDiffRelStep; evals++;
					lastX[i] += options->FinDiffRelStep;
				}
				//fprintf(stderr, "%.3e ", grad[i]);
			}else{
				grad[i] = 0;
			}

			if( options->Display)fprintf(stderr, " <- grad, alpha -> %e\n", alpha);
		}
		double abs_grad=0;
		for( int i=0; i<dim; i++)
			abs_grad += grad[i]*grad[i];
		abs_grad = sqrt( abs_grad);

		// decent gradient
		if( norm( grad, dim) > 0){

			for( int i=0; i<dim; i++)
				X[i] = lastX[i] - grad[i]*alpha/abs_grad;

			// step size respecting boundary?
			while( !inbound( X, lb, ub, dim)){
				alpha /= 2;
				for( int i=0; i<dim; i++){
					X[i] = lastX[i] - grad[i]*alpha/abs_grad;
				}
				if( alpha < options->TolX){
					//if( options->Display)
					if( options->Display)fprintf( stderr, "%i: dX %e < TolX %e\n", iter, alpha, options->TolX);
					return;
				}
			}

			// evaluation
			Y = (*func) ( dim, X, 0, func_data); evals++;
			//fprintf( stderr, "%i: %e (%e)\n", iter, Y, alpha);
			while( Y >= lastY){

				alpha /= 2;
				for( int i=0; i<dim; i++)
					X[i] = lastX[i] - grad[i]*alpha/abs_grad;
				Y = (*func) ( dim, X, 0, func_data); evals++;
				if( options->Display)
				{
					fprintf( stderr, "%i: %e <- {", iter, Y);
					for( int i=0; i<dim; i++)
						fprintf( stderr, "%e ", X[i]);
					fprintf( stderr, "} (alpha=%e)\n", alpha);
				}
				if( evals==options->MaxFunEvals || fabs(Y-lastY) < options->TolFun || alpha < options->TolX){
					if( options->Display)fprintf( stderr, "#eval = %i/%i, |dy| = %e (< %e), dx = %e(< %e)\n", evals, options->MaxFunEvals, fabs(Y-lastY), options->TolFun, alpha, options->TolX);
					return;
				}
			}


			// breaking criterion
			if( norm_diff( lastX, X, dim) < options->TolX){
				fprintf( stderr, "change (%.3e) < TolX\n", X[0] - lastX[0]);
				stop = true;
			}

			for( int i=0; i<dim; i++)
				lastX[i] = X[i];
			lastY = Y;

		}else{
			fprintf( stderr, "Zero gradient\n");
			stop = true;
		}



		if( options->Display)
			fprintf( stderr, "%i: %e <- %e\n", iter, Y, X[0]);

		//alpha *= 10;
	}
	fprintf( stderr, "Max. Iterations finished (alpha = %e)\n", alpha);
	//fprintf( stderr, "<< %e >>\n", lastX[0]); //exit(0);
}
