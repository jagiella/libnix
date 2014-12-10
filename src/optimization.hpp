/*
 * optimization.hpp
 *
 *  Created on: Oct 7, 2014
 *      Author: jagiella
 */

#ifndef OPTIMIZATION_HPP_
#define OPTIMIZATION_HPP_


// DATA STRUCTURES

enum ParameterScaling { Linear, Logarithmic};

typedef struct {
	bool   Display;
	double FinDiffRelStep;
	double TolFun;
	double TolX;
	int    MaxIter;
	int    MaxFunEvals;
	bool   GradObj;
	int    MaxFunEvalsAvg;
	char   ParameterScaling;
	int    PopulationSize;
	bool   AllowInterception;
	char   OutputFile[512];
} optimoptions;


// FUNCTION PROTOTYPES

optimoptions getoptions();

template <class T>
void LHS( T **x, int n, int d, unsigned int *seed){

	unsigned int own_seed = 0;
	if( !seed) seed = &own_seed;

	for(int i=0;i<n;++i)
		for(int j=0;j<d;++j)
			x[i][j]=(i + unifrnd<T>( seed)) / n;

	for (int i=n-1; i>=0; --i)
		for(int j=0;j<d;++j)
		{
			int ii = rand_r( seed) % (i+1);
		    T temp = x[i ][j];
		    x[i ][j]    = x[ii][j];
		    x[ii][j]    = temp;
		}
}

template <class T>
void LHS( T **x, T *lb, T *ub, int n, int d, unsigned int *seed){

	unsigned int own_seed = 0;
	if( !seed) seed = &own_seed;

	// uniformly distributed points
	for(int i=0;i<n;++i)
		for(int j=0;j<d;++j)
			x[i][j]=(i + unifrnd<T>( seed)) / n * (ub[j]-lb[j]) + lb[j];

	// Shuffling
	for (int i=n-1; i>=0; --i)
		for(int j=1;j<d;++j)
		{
			int ii = rand_r( seed) % (i+1);
		    T temp = x[i ][j];
		    x[i ][j]    = x[ii][j];
		    x[ii][j]    = temp;
		}
}

// >> LOCAL <<

double LineSearch(
	int parameterSize, double *parameters, double *parametersMin, double *parametersMax,
	double differentiationStepSize, double lambda, double minSquares, int maxIterations,
	double (*f)( int, const double*, double*, void *), void *f_data
	);

double LevenbergMarquardt(
	int parameterSize, double *parameters, double *parametersMin, double *parametersMax,
	double differentiationStepSize, double lambda, double minSquares, int maxIterations,
	double (*f)( int, const double*, double*, void *), void *f_data
	);

void gradientDecent(
		double *x0, double *lb, double *ub, int dim,
		double (*f) (int, const double*, double*, void*), void* f_data,
		optimoptions *options, double *sol);


// >> GLOBAL <<

void multistart(
		// input
		double **X0, double *lb, double *ub, int dim, double (*func) (int, const double*, double*, void*), void*,
		// options
		optimoptions *options,
		// methode
		void (*methode) (double*,double*,double*,int,double (*) (int, const double*, double*, void*),void*, optimoptions*,double*),
		// output
		double **sol, int n);

void gpopt(
		// input
		double *x0, double *lb, double *ub, int dim, double (*func) (int, const double*, double*, void*),
		// options
		optimoptions *options,
		// output
		double *sol);

void abc(
			// input
			double **X0, double *lb, double *ub, int dim, double (*func) (int, const double*, double*, void*), void*,
			// options
			optimoptions *options,
			// output
			double **&sol, int &n);

void abcgp(
			// input
			double **X0, double *lb, double *ub, int dim, double (*func) (int, const double*, double*, void*), void*,
			// options
			optimoptions *options,
			// output
			double **&sol, int &n);



#endif /* OPTIMIZATION_HPP_ */
