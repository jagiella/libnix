/*
 * optimization.cpp
 *
 *  Created on: Oct 7, 2014
 *      Author: jagiella
 */


#include <stdio.h>
#include <math.h>
#include <time.h>
#include <float.h>

#include "matrix.hpp"
#include "statistics.hpp"
#include "optimization.hpp"
#include "gp.hpp"

#define USE_OMP

extern double EpsilonLimit;

optimoptions getoptions(){
	optimoptions options;

	options.Display = false;
	options.FinDiffRelStep = 1e-6;
	options.TolFun = 1e-6;
	options.TolX   = 1e-6;
	options.MaxIter     = 400;
	options.MaxFunEvals = 1000;
	options.GradObj     = false;
	options.ParameterScaling = Linear;
	sprintf( options.OutputFile, "output");

	// gp-opt:
	options.MaxFunEvalsAvg = 10;

	// abc
	options.PopulationSize = 100;
	options.AllowInterception = false;

	return options;
}


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
			double beta = unifrnd<double>( &seed);
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
			double beta = unifrnd<double>( &seed);
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


typedef struct {
	double **x_old;
	double **invcov;
	double   detcov;
	int      n_old;
} kdedata;
double myfunckde(int n, const double *x, double *grad, void *my_func_data)
{
    return kdepdf( x, ((kdedata*)my_func_data)->x_old, ((kdedata*)my_func_data)->invcov, ((kdedata*)my_func_data)->detcov, n, ((kdedata*)my_func_data)->n_old);
}

void log2lin( double *xlog, double *xlin, int n){
	for( int i=0; i<n; i++)
		xlin[i] = pow( 10, xlog[i]);
}
int cumsum( double *x, int n, double partial_sum)
{
	double sum = 0;

	for( int i=0; i<n; i++){
		sum += x[i];
		if( sum >= partial_sum){
			return i;
		}
	}

	fprintf( stderr, "ERROR: cumsum (%e < %e)\n", sum, partial_sum);
	printVector( x, n, "%e ");
	return -1;
}

void abc(
			// input
			double **X0, double *lb, double *ub, int d, double (*func) (int, const double*, double*, void*), void* func_data,
			// options
			optimoptions *options,
			// output
			double **&sol, int &n)
{
	srand(time(NULL));
	// ABC parameters
	//int d = 2;
	//int options->PopulationSize = 100;
	int nparallel = 10;
	//int nthreads  = 5;
	//int options->MaxIter = 100;
	//int options->MaxFunEvals= 1000000;
	unsigned int func_seed = 0;
	//int binom_data[2] = { 24, 16};
	//	int binom_data[2] = { 6, 4};
	double quantile_value = 0.5;


	int iCPU = 10;
	//options->PopulationSize = 100;

	char suffix[1024];
	char filename[1024];
	sprintf( suffix, "%s", options->OutputFile );

	//enum {Log, Lin};
	//char SCALING = Lin;








	// init variables
	fprintf( stderr, "INIT\n");
	double epsilon = DBL_MAX;
	EpsilonLimit   = DBL_MAX;
	double **x_new, **x_old, **x_tmp, **x_par, *x_lin;
	double  *f_new,  *f_old,  *f_tmp,  *f_par;
	double  *w_new,  *w_old,  *w_tmp,  *w_par;
	int      n_new,   n_old,            n_par = nparallel;
	double _mean[d], **_cov = allocMatrix<double>(d,d), **invcov = allocMatrix<double>(d,d), **covL = allocMatrix<double>(d,d), detcov;

	x_old = allocMatrix<double>( options->PopulationSize+n_par, d);	f_old = (double*) malloc((options->PopulationSize+n_par)*sizeof(double));	w_old = (double*) malloc((options->PopulationSize+n_par)*sizeof(double)); n_old = 0;
	x_new = allocMatrix<double>( options->PopulationSize+n_par, d);	f_new = (double*) malloc((options->PopulationSize+n_par)*sizeof(double));	w_new = (double*) malloc((options->PopulationSize+n_par)*sizeof(double)); n_new = 0;
	x_par = allocMatrix<double>( n_par,  d);	f_par = (double*) malloc(n_par*sizeof(double));    	w_par = (double*) malloc(n_par*sizeof(double));

	if( options->ParameterScaling == Logarithmic)
		x_lin = (double*) malloc(d*sizeof(double));;

#ifdef USE_OMP
	//int iCPU = nparallel;//= omp_get_num_procs();
	fprintf( stderr, "SET NUM THREADS: %i\n", iCPU);
	omp_set_num_threads( iCPU);
#endif
	//printMatrix( x_old, n_old, d, "%10.3e ");

	double x_best[d], f_best = DBL_MAX;

	// stats
	clock_t t, t_sampling = 0, t_evaluation = 0;
	sprintf( filename, "%s.population.dat", suffix);
	FILE *fp_population = fopen( filename, "w+");

	// run iteration
	int count_evals_per_iteration[options->MaxIter];
	int count_evals = 0;
	unsigned int sampling_seed = 0;

	fprintf( stderr, " iteration  #fun eval acceptance    epsilon\n");

	for( int it=0; it<options->MaxIter && count_evals < options->MaxFunEvals; it++){

		n_new = 0;

		// COPY OLD ONES < EPSILON
		/*for( int i=0; i<n_old; i++)
			if( f_old[i] < epsilon){
				for( int j=0; j<d; j++)
				x_new[n_new][j] = x_old[i][j];
				f_new[n_new] = f_old[i];
				w_new[n_new] = 1 / kdepdf( x_old[i], x_old, invcov, detcov, d, n_old);
				n_new ++;
			}*/


		// sample new points
		count_evals_per_iteration[it]=0;

		while( n_new<options->PopulationSize){

			// >>> SEQUENTIAL PART <<<
			if( options->Display)
				fprintf( stderr, "SAMPLING\n");
			for( int ip=0; ip<nparallel; ip++){

				t = clock();
				//if( options->Display)
				//	fprintf( stderr, "SAMPLING %i\n", omp_get_thread_num());
				if( it == 0){
					// UNIFORM SAMPLING
					for( int j=0; j<d; j++)
						x_par[ip][j] = lb[j] + (ub[j]-lb[j])*unifrnd<double>(&sampling_seed);
				}else{
					// KERNEL DENSITY ESTIMATE
					// choice point
					int idx = cumsum( w_old, n_old, unifrnd<double>(&sampling_seed));

					// sample from MVnorm dist around chosen point
					do{
						mvnrnd( x_par[ip], x_old[idx], covL, d, &sampling_seed);
						//printVector( lb, d, "%5.3f ");
						//printVector( x_par[ip], d, "%5.3f ");
						//printVector( ub, d, "%5.3f ");
						//printf( "%i\n", inbound( x_par[ip], lb, ub, d));
					}while( !inbound( x_par[ip], lb, ub, d));
				}
				t_sampling += clock() - t;
			}

			// >>> PARALLEL PART <<<

			if( options->Display)
				fprintf( stderr, "EVALUATION\n");

	#ifdef USE_OMP
	#pragma omp parallel for
	#endif
			for( int ip=0; ip<nparallel; ip++){
				// EVALUATION

				t = clock();
				//fprintf( stderr, "EVALUATION THREAD %i\n", omp_get_thread_num());
				switch( options->ParameterScaling ){
				case Linear:
					f_par[ip] = (*func)( d, x_par[ip], 0, func_data); break;
				case Logarithmic:
					log2lin( x_par[ip], x_lin, d);
					f_par[ip] = (*func)( d, x_lin, 0, func_data); break;
				}

				t_evaluation += clock() - t;


				// WEIGHTS
				//fprintf( stderr, "WEIGHTS   THREAD %i\n", omp_get_thread_num());
				w_par[ip]=0;
				if( f_par[ip] < epsilon){
					if( it == 0){
						w_par[ip] = 1;
					}else{
						w_par[ip] = 1 / kdepdf( x_par[ip], x_old, invcov, detcov, d, n_old);
					}
				}

				//fprintf( stderr, " <- End   THREAD %i\n", omp_get_thread_num());

			}
			count_evals += nparallel;
			count_evals_per_iteration[it] += nparallel;



			// >>> SEQUENTIAL PART <<<
			if( options->Display)
				fprintf( stderr, "STORE ACCEPTED\n");

			for( int ip=0; ip<nparallel; ip++){
				//fprintf( stderr, "%i = %5e \n", ip, f_par[ip]);
				//fprintf( stderr, "%s\n", (w_par[ip] ? "accepted":"refused"));
				if( w_par[ip]){
					// COPY TO NEXT GENERATION
					w_new[n_new] = w_par[ip];
					f_new[n_new] = f_par[ip];
					for( int id=0; id<d; id++)
						x_new[n_new][id] = x_par[ip][id];
					n_new++;
				}
			}

			//return 0;

			if( count_evals >= options->MaxFunEvals){
				//break;
				fprintf( stderr, "END\n");
				double f_min;
				int    i_min;
				min( f_old, n_old, f_min, i_min);
				printVector( x_old[i_min], d, "%e ");

				fprintf ( stderr,
						"TIME STATS: sampling = %e seconds per point, evaluation = %e sec / point.\n",
						( (float)t_sampling  /count_evals)/CLOCKS_PER_SEC,
						( (float)t_evaluation/count_evals)/CLOCKS_PER_SEC);

				// COPY SOLUTION
				n = options->PopulationSize;
				sol = allocMatrix<double>( n, d);
				for( int i=0; i<n; i++)
					for( int j=0; j<d; j++)
						sol[i][j] = x_old[i][j];
				return;

				// PRINT KDE
				sprintf( filename, "%s.kde.dat", suffix);
				FILE *fp_kde = fopen( filename, "w+");
				double dx = 0.1;
				double intergral_kde = 0;
				if(d==1){
					dx = 0.01;
					for( double x=lb[0]; x<ub[0]; x+=dx){
							double xtest[1] = {x};
							fprintf( fp_kde, "%e %e\n", x, kdepdf( xtest, x_old, invcov, detcov, d, n_old));
							intergral_kde += dx*kdepdf( xtest, x_old, invcov, detcov, d, n_old);
						}
				}

				if(d==2){
					int n0 = (int) (ub[0]-lb[0])/dx + 1;
					int n1 = (int) (ub[1]-lb[1])/dx + 1;
					double kde[n0][n1];
					double lik[n0][n1];
#ifdef USE_OMP
#pragma omp parallel for
#endif
					for( int i0=0; i0<n0; i0++){
						for( int i1=0; i1<n1; i1++){
							double x = lb[0] + dx*i0, y = lb[1] + dx*i1;
							double xtest[2] = {x, y};
							kde[i0][i1] = kdepdf( xtest, x_old, invcov, detcov, d, n_old);
							lik[i0][i1] = (*func)( d, xtest, 0, func_data);
						}
					}
					for( int i0=0; i0<n0; i0++){
						for( int i1=0; i1<n1; i1++){
							double x = lb[0] + dx*i0, y = lb[1] + dx*i1;
							fprintf( fp_kde, "%e %e %e %e\n", x,y, kde[i0][i1], lik[i0][i1]);
							intergral_kde += dx*dx*kde[i0][i1];
						}
						fprintf( fp_kde, "\n");
					}
				}
				fclose( fp_kde);
				fprintf ( stderr, "[ %e ]\n", intergral_kde);

			}
			//fprintf( stderr, "---------> %5i+%4i %5i / %i\n", count_evals-count_evals_per_iteration[it], count_evals_per_iteration[it], n_new, options->PopulationSize );
//			fprintf( stderr, "%10i %10i %9.3f%% %10.3e\n", it, count_evals, 100*(double)n_new/count_evals_per_iteration[it], epsilon);
			fprintf( stderr, "\r%10i %10i [%3i/%4i] %10.3e\r", it, count_evals, n_new, options->PopulationSize, epsilon);

			//if( count_evals_per_iteration[it] >= 10 && (double)n_new/count_evals_per_iteration[it] < 0.25)
			//	epsilon *= 1.1;

		}
		//if( /*count_evals >= 10 &&*/ (double)n_new/count_evals_per_iteration[it] < 0.25)
		//	quantile_value *= 1.1;

		if( options->Display)
			fprintf( stderr, "PREPARE NEXT ITERATION\n");

		normalizeVector( w_new, n_new);
//		EpsilonLimit = DBL_MAX;
		//EpsilonLimit = epsilon;
		epsilon = quantile( f_new, n_new, quantile_value);
		if(	options->AllowInterception)
			EpsilonLimit = epsilon;
		else
			EpsilonLimit = DBL_MAX;

		for( int i=0; i<n_new; i++){
			for( int j=0; j<d; j++)
				fprintf( fp_population, "%e ", x_new[i][j]);
			fprintf( fp_population, "%e %e %i %i %e\n", f_new[i], w_new[i], it, count_evals, epsilon);
		}
		fprintf( fp_population, "\n\n");

		//fprintf( stderr, "new epsilon = %e\n", epsilon);



		fprintf( stderr, "%10i %10i %9.3f%% %10.3e\n", it, count_evals, 100*(double)n_new/count_evals_per_iteration[it], epsilon);

		x_tmp = x_old; x_old = x_new; x_new = x_tmp;
		f_tmp = f_old; f_old = f_new; f_new = f_tmp;
		w_tmp = w_old; w_old = w_new; w_new = w_tmp;
		n_old = n_new;                n_new = 0;


		//fprintf( stderr, "[ ITERATION %i ]\n", it);
		count_evals_per_iteration[it]=0;

		fprintf( stderr, "mean: ");
		mean( x_old, n_old, d, _mean, 2);
		if( options->ParameterScaling == Logarithmic){
			log2lin( _mean, x_lin, d);
			printVector( x_lin, d, "%10e ");
		}else
			printVector( _mean, d, "%10e ");

		int i_temp; double f_temp;
		min( f_old, n_old, f_temp, i_temp);
		if( f_best > f_temp){
			f_best = f_temp;
			for( int i=0; i<d; i++)
				x_best[i] = x_old[i_temp][i];
		}
		fprintf( stderr, "best: ");
		printVector( x_best, d, "%10e ");

		//fprintf( stderr, "cov\n");
		cov(  x_old, n_old, d, _mean, _cov);
		//printMatrix( cov, d, d, "%10e ");

		detcov = determinantLU( _cov, d); //fprintf( stderr, "|cov| = %e\n", detcov);

		//fprintf( stderr, "L\n");
		decompCholeskyOwn(_cov, covL, d); // Cholesky-decomposition
		//printMatrix( covL, d, d, "%10e ");

		//fprintf( stderr, "L*L^T\n");
		/*double **test = allocMatrix<double>( d,d);
		for( int i=0; i<d; i++)
			for( int j=0; j<d; j++){
				test[i][j] = 0;
				for( int k=0; k<d; k++)
					test[i][j] += covL[i][k]*covL[j][k];
			}*/
		//printMatrix( test, d, d, "%10e ");


		//inverseCholesky( cov, invcov, d);
		inverseLU<double>( _cov, invcov, d);

	}


	fprintf( stderr, "END\n");
	double f_min;
	int    i_min;
	min( f_old, n_old, f_min, i_min);
	printVector( x_old[i_min], d, "%e ");

	fprintf ( stderr,
			"TIME STATS: sampling = %e seconds per point, evaluation = %e sec / point.\n",
			( (float)t_sampling  /count_evals)/CLOCKS_PER_SEC,
			( (float)t_evaluation/count_evals)/CLOCKS_PER_SEC);

}


void summary_statistics(
		double (*func) (int, const double*, double*, void*), void *func_data, int d, const double *parameters, int n,
		double &y, double &s2)
{
	y  = 0;
	s2 = 0;
	for( int i=0; i<n; i++){
		double val = (*func) (d, parameters, 0, func_data);
		y  += val;
		s2 += val*val;
	}
	y  /= (double)n;
	s2  = s2/(double)n - pow(y, 2); // s2
	s2 /= (double)n; // SEM^2

}

void transform( double **xin, double **xout, int n, int d, double *m, double **inv_sqrt_S){
	for( int i=0; i<n; i++){
		for( int di=0; di<d; di++){
			xout[i][di] = 0;
			for( int dj=0; dj<d; dj++)
				xout[i][di] += inv_sqrt_S[di][dj]*(xin[i][dj]-m[dj]);
		}
	}
}



void abcgp(
			// input
			double **X0, double *lb, double *ub, int d, double (*func) (int, const double*, double*, void*), void* func_data,
			// options
			optimoptions *options,
			// output
			double **&sol, int &nsol)
{
	srand(time(NULL));
	// ABC parameters
	//int d = 2;
	//int options->PopulationSize = 100;
	int nparallel = 10;
	//int nthreads  = 5;
	//int options->MaxIter = 100;
	//int options->MaxFunEvals= 1000000;
	unsigned int func_seed = 0;
	//int binom_data[2] = { 24, 16};
	//	int binom_data[2] = { 6, 4};


	int iCPU = 10;
	//options->PopulationSize = 100;

	char suffix[1024];
	char filename[1024];
	sprintf( suffix, "%s", options->OutputFile );

	//enum {Log, Lin};
	//char SCALING = Lin;
	bool TRANSFORM = true;


	double (*covFunc) ( double* &, double* &, int &, double* &) = covWendland;




		// init variables
	fprintf( stderr, "INIT\n");
	//double epsilon = DBL_MAX;
	//double **x_new, **x_old, **x_tmp, **x_par;
	//double  *f_new,  *f_old,  *f_tmp,  *f_par;
	//double  *w_new,  *w_old,  *w_tmp,  *w_par;
	//int      n_new,   n_old,            n_par = nparallel;
	double _mean[d], **_cov = allocMatrix<double>(d,d), **invcov = allocMatrix<double>(d,d), **invcovL = allocMatrix<double>(d,d), **covL = allocMatrix<double>(d,d), detcov;

	int n_max = 10000;
	double **x, *f, *s2, *w;
	x = allocMatrix<double>( n_max, d);
	double **x_tf = allocMatrix<double>( n_max, d);
	f = (double*) malloc(n_max*sizeof(double));
	s2= (double*) malloc(n_max*sizeof(double));
	w = (double*) malloc(n_max*sizeof(double));
	int n = 0;
	//int n_avg = 4;

	double *x_lin = 0;
	if( options->ParameterScaling == Logarithmic)
		x_lin = (double*) malloc(d*sizeof(double));;


#ifdef USE_OMP
	//int iCPU = nparallel;//= omp_get_num_procs();
	fprintf( stderr, "SET NUM THREADS: %i\n", iCPU);
	omp_set_num_threads( iCPU);
#endif
	//printMatrix( x_old, n_old, d, "%10.3e ");


	// stats
	clock_t t, t_sampling = 0, t_evaluation = 0;


	// run iteration
	int count_evals_per_iteration[options->MaxIter];
	int count_evals = 0;
	unsigned int sampling_seed = 0;

	// INITIAL SAMPLING
	sprintf( filename, "%s.population.dat", suffix);
	FILE *fp_population = fopen( filename, "w+");
	LHS( x, lb, ub, options->PopulationSize, d, &sampling_seed);
	count_evals = 0;
	count_evals_per_iteration[0] = 0;
	for( int i=0; i<options->PopulationSize; i++){
		for( int j=0; j<d; j++){
			//x[i][j] = lb[j] + (ub[j]-lb[j])*RAND01_R(&sampling_seed);
			fprintf( fp_population, "%e ", x[i][j]);
		}
		w[i] = 1./options->PopulationSize;

		switch( options->ParameterScaling){
		case Linear:
			summary_statistics( func, func_data, d, x[i], options->MaxFunEvalsAvg, f[i], s2[i]); break;
		case Logarithmic:
			for( int j=0; j<d; j++) x_lin[j] = pow( 10, x[i][j]);
			summary_statistics( func, func_data, d, x_lin, options->MaxFunEvalsAvg, f[i], s2[i]); break;
		}
		count_evals += options->MaxFunEvalsAvg;
		count_evals_per_iteration[0] += options->MaxFunEvalsAvg;

		fprintf( fp_population, "%e %e %i %i %e\n", f[i], s2[i], 0, count_evals, DBL_MAX);

		fprintf( stderr, "\r%10i %10i %9.3f%% %10.3e\r", 0, count_evals, 100*(double)count_evals_per_iteration[0]/(options->MaxFunEvalsAvg*options->PopulationSize), DBL_MAX);
	}
	fprintf( stderr, "%10i %10i %9.3f%% %10.3e\n", 0, count_evals, 100*(double)count_evals_per_iteration[0]/(options->MaxFunEvalsAvg*options->PopulationSize), DBL_MAX);

	n += options->PopulationSize;
	fprintf( fp_population, "\n\n");
	fclose(fp_population);

	sprintf( filename, "%s.prediction.dat", suffix);
	fp_population = fopen( filename, "w+");
	fclose(fp_population);

	int      n_kde = options->PopulationSize/2;
	double **x_kde = (double**)malloc(sizeof(double*)*n_kde);
	double  *f_kde = (double*)malloc(sizeof(double)*n_kde);
	double  *w_kde = (double*)malloc(sizeof(double)*n_kde);


	fprintf( stderr, " iteration  #fun eval acceptance    epsilon\n");

	for( int it=1; it<options->MaxIter && count_evals < options->MaxFunEvals; it++){

		count_evals_per_iteration[it] = 0;

		// Prepare next iteration
		double **x_new = &x[n];
		double **x_new_tf = &x_tf[n];
		double  *f_new = &f[n];
		double  *s2_new= &s2[n];
		double  *w_new = &w[n];
		int      n_new = 0;

		// KDE stuff
		double epsilon_kde = quantile( &f[n-options->PopulationSize], options->PopulationSize, n_kde/(double)options->PopulationSize);
		/*double **x_kde = &x[n-options->PopulationSize];
		double  *f_kde = &f[n-options->PopulationSize];
		double  *w_kde = &w[n-options->PopulationSize];*/
		int i_kde=0;
		for( int i=0; i<options->PopulationSize && i_kde<n_kde; i++)
		if( f[n-options->PopulationSize+i] <= epsilon_kde){
			f_kde[i_kde] = f[n-options->PopulationSize+i];
			w_kde[i_kde] = w[n-options->PopulationSize+i];
			x_kde[i_kde] = x[n-options->PopulationSize+i];
			i_kde++;
		}
		normalizeVector( w_kde, n_kde);
		//getKDESample(x,f,n, x_kde,f_kde,n_kde);
		// covariance
		mean( x_kde, n_kde, d, _mean, 2);
		if( options->ParameterScaling == Logarithmic){
			log2lin( _mean, x_lin, d);
			printVector( x_lin, d, "%10e ");
		}else
			printVector( _mean, d, "%10e ");

		cov(  x_kde, n_kde, d, _mean, _cov);
		// determinant
		detcov = determinantLU( _cov, d); //fprintf( stderr, "|cov| = %e\n", detcov);
		// sqrt of coveriance
		decompCholeskyOwn(_cov, covL, d); // Cholesky-decomposition
		// inverse of covariance
		inverseLU( _cov, invcov, d);


		// GP hyperparameters
		int      n_gp = n;//min( n, 5*options->PopulationSize);
		double **x_gp = &x[n-n_gp];
		double  *f_gp = &f[n-n_gp];
		double  *s2_gp= &s2[n-n_gp];


		double hyp[4], *phyp = hyp;
		//getHyperParameters( x_gp, f_gp, n_gp, d, phyp);
		int n_lim = n;//min( n_gp, 3*options->PopulationSize);
		if( TRANSFORM){
			inverseLU( covL, invcovL, d);
			transform( x_gp, x_tf, n_gp, d, _mean, invcovL);
			getHyperParameters( &x_tf[n_gp-n_lim], &f_gp[n_gp-n_lim], n_lim, d, phyp);
		}else
			getHyperParameters( &x_gp[n_gp-n_lim], &f_gp[n_gp-n_lim], n_lim, d, phyp);
		//printVector( hyp, 4, "%.2e ");

		// K
		double **K = allocMatrix<double>( n_gp, n_gp);
		if( TRANSFORM)
			covMatrix<double>( K, x_tf, n_gp,d, covFunc, phyp);
		else
			covMatrix<double>( K, x_gp, n_gp,d, covFunc, phyp);
		for(int i=0; i<n_gp; i++)
			K[i][i] += s2_gp[i] + pow(10,phyp[0]);

		// inverse of K
		double **invK = allocMatrix<double>( n_gp, n_gp);
		inverseLU( K, invK, n_gp);



		// Threshold
		double f_kde_gp[n_kde];
		double s2_kde_gp[n_kde];
		if( TRANSFORM)
			evalVariance<double>( x_tf, f_gp, s2_gp, n_gp, x_tf, f_kde_gp, s2_kde_gp, n_kde, d, covFunc, phyp, invK);
		else
			evalVariance<double>( x_gp, f_gp, s2_gp, n_gp, x_gp, f_kde_gp, s2_kde_gp, n_kde, d, covFunc, phyp, invK);
		double epsilon_gp = quantile( f_kde_gp, n_kde, 0.1);

		//epsilon = quantile( f_kde, n_kde, 0.5);


		while( n_new<options->PopulationSize){

			// KERNEL DENSITY ESTIMATE

			// choice point
			int idx;
			do{
				idx = cumsum( w_kde, n_kde, unifrnd<double>(&sampling_seed));
			}while( f_kde[idx] > epsilon_kde);

			// sample from MVnorm dist around chosen point
			do{
				mvnrnd( x_new[n_new], x_kde[idx], covL, d, &sampling_seed);
			}while( !inbound( x_new[n_new], lb, ub, d));


			// EVALUATION OF GP
			//double x_new_tf[d], *px_new_tf=x_new_tf;

			if( TRANSFORM){
				transform( &x_new[n_new], &x_new_tf[n_new], 1, d, _mean, invcovL);
				evalVariance<double>( x_tf, f_gp, s2_gp, n_gp, &x_new_tf[n_new], &f_new[n_new], &s2_new[n_new], 1, d, covFunc, phyp, invK);
	//			double dummy;
	//			evalVariance<double>( x_tf, f_gp, s2_gp, n_gp, &x_new_tf[n_new], &f_new[n_new], &dummy, 1, d, covFunc, phyp, invK);
	//			evalVariance<double>( x_tf, s2_gp, 0,    n_gp, &x_new_tf[n_new], &s2_new[n_new], &dummy, 1, d, covFunc, phyp, invK);
			}else
				evalVariance<double>( x_gp, f_gp, s2_gp, n_gp, &x_new[n_new], &f_new[n_new], &s2_new[n_new], 1, d, covFunc, phyp, invK);


			// SAMPLE FROM GP
			double rndval = f_new[n_new] + normrnd( &sampling_seed) * sqrt( s2_new[n_new]);


			// ADD POINT
			//fprintf( stderr, "WEIGHTS   THREAD %i\n", omp_get_thread_num());
			if( rndval < epsilon_gp){
				// EVALUATION OF MODEL
				switch( options->ParameterScaling){
				case Linear:
					summary_statistics( func, func_data, d, x_new[n_new], options->MaxFunEvalsAvg, f_new[n_new], s2_new[n_new]); break;
				case Logarithmic:
					for( int i=0; i<d; i++) x_lin[i] = pow( 10, x_new[n_new][i]);
					summary_statistics( func, func_data, d, x_lin, options->MaxFunEvalsAvg, f_new[n_new], s2_new[n_new]); break;
				}

				count_evals_per_iteration[it] += options->MaxFunEvalsAvg;
				count_evals += options->MaxFunEvalsAvg;

				// WEIGHT
				w_new[n_new] = 1 / kdepdf( x_new[n_new], x_kde, invcov, detcov, d, n_kde);
				//w_new[n_new] = 1;
				n_new++;
			}
			fprintf( stderr, "\r%10i %10i %9.3f%% %10.3e\r", it, count_evals, 100*(double)count_evals_per_iteration[it]/(options->MaxFunEvalsAvg*options->PopulationSize), epsilon_kde);

			if( count_evals >= options->MaxFunEvals){

				// TODO: STOP HERE

			}

		}
		normalizeVector( w_new, n_new);

		if( d==2){
			sprintf( filename, "%s.prediction.dat", suffix);
			fp_population = fopen( filename, "a+");
			double _x[2], *_px=_x, _f, _s2;
			for( _x[0]=lb[0]; _x[0]<=ub[0]; _x[0]+=0.1){
			for( _x[1]=lb[1]; _x[1]<=ub[1]; _x[1]+=0.1)
			{
				for( int j=0; j<d; j++){
					fprintf( fp_population, "%e ", _x[j]);
				}

				if( TRANSFORM){
					double _xt[d], *_pxt=_xt;
					transform( &_px, &_pxt, 1, d, _mean, invcovL);
					/*for( int j=0; j<d; j++){
						fprintf( fp_population, "%e ", _xt[j]);
					}*/
					evalVariance<double>( x_tf, f_gp, s2_gp, n_gp, &_pxt, &_f, &_s2, 1, d, covFunc, phyp, invK);
				}else
					evalVariance<double>( x_gp, f_gp, s2_gp, n_gp, &_px, &_f, &_s2, 1, d, covFunc, phyp, invK);

				fprintf( fp_population, "%e %e %i %e\n", _f, _s2, count_evals, epsilon_kde);
			}
			fprintf( fp_population, "\n");
						}
			fprintf( fp_population, "\n");
			fclose(fp_population);
		}


		freeMatrix( K,    n_gp);
		freeMatrix( invK, n_gp);

		n += n_new;

		fprintf( stderr, "%10i %10i %9.3f%% %10.3e\n", it, count_evals, 100*(double)n_new*options->MaxFunEvalsAvg/(count_evals_per_iteration[it]), epsilon_kde);



		sprintf( filename, "%s.population.dat", suffix);
		fp_population = fopen( filename, "a+");
		for( int i=0; i<options->PopulationSize; i++)
		if( f_new[i] < epsilon_kde){
			for( int j=0; j<d; j++){
				fprintf( fp_population, "%e ", x_new[i][j]);
			}
			fprintf( fp_population, "%e %e %i %i %e\n", f_new[i], s2_new[i], it+1, count_evals, epsilon_kde);
		}
		fprintf( fp_population, "\n\n");
		fclose(fp_population);
	}

	// COPY SOLUTION
	nsol = options->PopulationSize;
	sol = allocMatrix<double>( n, d);
	for( int i=0; i<n; i++)
		for( int j=0; j<d; j++)
			sol[i][j] = x[n-nsol+i][j];
	return;
}
