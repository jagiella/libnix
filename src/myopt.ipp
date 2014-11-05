/*
 * myopt.ipp
 *
 *  Created on: Sep 2, 2014
 *      Author: jagiella
 */

#include "statistics.hpp"
#include "matrix.hpp"
#include "optimization.hpp"



#include <math.h>

bool uniform_norm_le( double *x1, double *x2, double max, int dim)
{
	for( int i=0; i<dim; i++)
		if( fabs(x1[i]-x2[i]) > max)
			return false;
	return true;
}

bool uniform_norm_le( double *x, double max, int dim)
{
	for( int i=0; i<dim; i++)
		if( fabs(x[i]) > max)
			return false;
	return true;
}
void multistart(
		// input
		double **X0, double *lb, double *ub, int dim, double (*func) (int, const double*, double*, void*), void *func_data,
		// options
		optimoptions *options,
		// methode
		void (*methode) (double*,double*,double*,int,double (*) (int, const double*, double*, void*),void*, optimoptions*,double*),
		// output
		double **sol, int n)
{
	for( int i=0; i<n; i++){
		(*methode) ( X0[i], lb, ub, dim, func, func_data, options, sol[i]);
		//for( int d=0; d<dim; d++)
		//	fprintf(stderr, "%e ", sol[i][d]);
		//fprintf(stderr, "\n");
	}
}


typedef struct {
	double *x, *y, *w;
	char methode, weight;
	int     n;
	double (*function) ( double, double*, int);
} linregdata;

typedef struct {
	double **x, *y, *w;
	char methode, weight;
	int     n, d;
	double (*function) ( double, double*, int);
} multiregdata;

double polynom( double x, double *a, int dim){

	double y=a[0];
	for( int d=1; d<dim; d++)
		y += a[d]*pow( x, d);
	return y;
}
double dpolynom( double x, double *a, int dim){

	double dy=a[1];
	for( int d=2; d<dim; d++)
		dy += a[d]*pow( x, d-1)*d;
	return dy;
}

double plane( double *x, int size_x, double *a, int size_a){

	double y;
	//if( size_a == 3)		y = a[0] + a[1]*x[0] + a[2]*x[1]; return y;
	//if( size_a == 4)		y = a[0] + a[1]*x[0] + a[2]*x[1] + a[3]*x[0]*x[1];


	// y = a[0]*x[0] + a[1]*x[1] + a[2] + a[3]*x[0]*x[1];
	y = a[0]; // const


	for( int i=0; i<size_x; i++)
		y += a[i+1]*x[i]; // linear

	if( size_a > 1 + size_x){
		double prod_x = 1;
		for( int i=0; i<size_x; i++){
			prod_x *= x[i];
		}
		y += a[size_x+1]*prod_x;
	}

	return y;
}
void dplane( double *x, int size_x, double *a, int size_a, double *grad){

	//if( size_a == 3){ grad[0] = a[1]; grad[1] = a[2]; } return;
	//if( size_a == 4){
	//	grad[0] = a[1] + a[3]*x[1];
	//	grad[1] = a[2] + a[3]*x[1]; }

	for( int i=0; i<size_x; i++)
		grad[i] = a[i+1]; // linear

	if( size_a > 1 + size_x)
		for( int i=0; i<size_x; i++){
			double prod_x=1;
			for( int j=0; j<size_x; j++)
				if(i!=j)
					prod_x *= x[j];
			grad[i] += a[size_x+1]*prod_x;
		}
}


void linregweigths( int d, const double *g, double *grad, void *my_func_data){
	double *w = ((linregdata*)my_func_data)->w;
	double *x = ((linregdata*)my_func_data)->x;
	double *y = ((linregdata*)my_func_data)->y;
	int     n = ((linregdata*)my_func_data)->n;
	//char methode = ((linregdata*)my_func_data)->methode;
	char weight  = ((linregdata*)my_func_data)->weight;

	double k = 0;
	if( weight > 0){
		double e[n];
		for( int i=0; i<n; i++)
			e[i] = polynom( x[i], (double*)g, d) - y[i];
		k = mad( e, n) / 0.6745;
	}

	for( int i=0; i<n; i++){
		double e = polynom( x[i], (double*)g, d) - y[i];

		switch( weight){
		case 0: // LEAST SQUARES
			w[i] = 1;
			break;

		case 1: // HUBER
			if( fabs(e) <= k)
				w[i] = 1;
			else
				w[i] = k/fabs(e);
			break;

		case 2: // BISQUARE
			if( fabs(e) <= k)
				w[i] = pow(1-pow(e/k,2),2);
			else
				w[i] = 0;
			break;
		}
	}
}
void multiregweigths( int ng, const double *g, double *grad, void *my_func_data){
	double *w = ((multiregdata*)my_func_data)->w;
	double **x = ((multiregdata*)my_func_data)->x;
	double *y = ((multiregdata*)my_func_data)->y;
	int     d = ((multiregdata*)my_func_data)->d;
	int     n = ((multiregdata*)my_func_data)->n;
	//char methode = ((linregdata*)my_func_data)->methode;
	char weight  = ((multiregdata*)my_func_data)->weight;

	double k = 0;
	if( weight > 0){
		double e[n];
		for( int i=0; i<n; i++)
			e[i] = plane( x[i], d, (double*)g, ng) - y[i];
		k = mad( e, n) / 0.6745;
	}

	for( int i=0; i<n; i++){
		double e = plane( x[i], d, (double*)g, ng) - y[i];

		switch( weight){
		case 0: // LEAST SQUARES
			w[i] = 1;
			break;

		case 1: // HUBER
			if( fabs(e) <= k)
				w[i] = 1;
			else
				w[i] = k/fabs(e);
			break;

		case 2: // BISQUARE
			if( fabs(e) <= k)
				w[i] = pow(1-pow(e/k,2),2);
			else
				w[i] = 0;
			break;
		}
	}
}

void linreg_firstguess( int d, double *g, double *x, double *y, int n){

	// find center of mass
	double mean_x = median( x, n);
	double mean_y = median( y, n);

	// mean grad
	g[1] = 0;
	double gradient[n];
	for( int i=0; i<n; i++){
		//g[1] += (y[i]-mean_y) / (x[i]-mean_x);
		gradient[i] = (y[i]-mean_y) / (x[i]-mean_x);
	}
	//g[1] /= n;
	g[1] = median(gradient, n);

	g[0] = -g[1]*mean_x + mean_y;

	for( int i=2; i<d; i++)
		g[i] = 0;
}
double linreg( int d, const double *g, double *grad, void *my_func_data)
{

	double *w = ((linregdata*)my_func_data)->w;
	double *x = ((linregdata*)my_func_data)->x;
	double *y = ((linregdata*)my_func_data)->y;
	int     n = ((linregdata*)my_func_data)->n;
	char methode = ((linregdata*)my_func_data)->methode;
	char weight  = ((linregdata*)my_func_data)->weight;
	//double (*f) ( double, double*, int) = ((linregdata*)my_func_data)->function;

	double sqrErr = 0;

	double k = 0;
	if( methode > 0 || weight > 0){
		double e[n];
		for( int i=0; i<n; i++)
			e[i] = polynom( x[i], (double*)g, d) - y[i];
		k = mad( e, n) / 0.6745;
	}
	//k = 1.345;


	// objective function
	sqrErr = 0;
	for( int i=0; i<n; i++){
		double p, e = polynom( x[i], (double*)g, d) - y[i];

		switch( methode){
		case 0: // LEAST SQUARES
			p = pow( e, 2);
			break;

		case 1: // HUBER
			if( fabs(e) <= k)
				p = 0.5 * pow( e, 2);
			else
				p = k*fabs(e) - 1/2.*k*k;
			break;

		case 2: // BISQUARE
			if( fabs(e) <= k)
				p = k*k/6.*(1-pow(1-pow(e/k,2),3));
			else
				p = k*k/6.;
			break;
		}

		/*switch( weight){
		case 0: // LEAST SQUARES
			w[i] = 1;
			break;

		case 1: // HUBER
			if( fabs(e) <= k)
				w[i] = 1;
			else
				w[i] = k/fabs(e);
			break;

		case 2: // BISQUARE
			if( fabs(e) <= k)
				w[i] = pow(1-pow(e/k,2),2);
			else
				w[i] = 0;
			break;
		}*/

		sqrErr += w[i]*w[i]*p;

		//sqrErr += p;
	}

	// LASSO - Penalty
	double lassoPenalty = 0;
	for( int i=0; i<d; i++)
		lassoPenalty += fabs(g[i]);

	//sqrErr += lassoPenalty*1e0;

	//fprintf(stderr, "[ linreg: (%.2e, %.2e) -> %e ]\n", g[0], g[1], sqrErr);

	return sqrErr;
}

double linregnoise( int d, const double *g, double *grad, void *my_func_data)
{

	double *w = ((linregdata*)my_func_data)->w;
	double *x = ((linregdata*)my_func_data)->x;
	double *y = ((linregdata*)my_func_data)->y;
	int     n = ((linregdata*)my_func_data)->n;
	char methode = ((linregdata*)my_func_data)->methode;
	char weight  = ((linregdata*)my_func_data)->weight;
	//double (*f) ( double, double*, int) = ((linregdata*)my_func_data)->function;

	// L = 1/2. * sum( log(2*pi*s2) + (m - y).^2 / s2 );

	double sqrErr = 0;
	double s2 = pow(10, g[d-1]);
	for( int i=0; i<n; i++){
		double e, yavg = polynom( x[i], (double*)g, d-1);
		double p;
		switch(2){
		case 1: // normal distribution
			e = y[i] - yavg;
			p = log(2*M_PI*s2) + pow( e, 2) / s2;
			break;
		case 2: // log-normal distribution
			if( y[i] > yavg){
				e = log(y[i]-yavg);
				p = log(2*M_PI*s2*pow(y[i]-yavg, 2)) + pow( e, 2) / s2;
			}else{
				e = (y[i] - yavg);
				p = pow( e, 2) * 1e1;
				//p = pow( e, 2) * floor(g[d-1]);
			}
			break;
		}
		sqrErr += p;
	}

	return 0.5*sqrErr;
}

double multireg( int size_a, const double *a, double *grad, void *my_func_data)
{

	double *w = ((multiregdata*)my_func_data)->w;
	double **x = ((multiregdata*)my_func_data)->x;
	double *y = ((multiregdata*)my_func_data)->y;
	int     n = ((multiregdata*)my_func_data)->n;
	int     d = ((multiregdata*)my_func_data)->d;
	char methode = ((multiregdata*)my_func_data)->methode;
	char weight  = ((multiregdata*)my_func_data)->weight;

	double sqrErr = 0;

	double k = 0;
	if( methode > 0 || weight > 0){
		double e[n];
		for( int i=0; i<n; i++)
			e[i] = plane( x[i], d, (double*)a, size_a) - y[i];
		k = mad( e, n) / 0.6745;
	}
	//k = 1.345;


	// objective function
	sqrErr = 0;
	for( int i=0; i<n; i++){
		double p, e = plane( x[i], d, (double*)a, size_a) - y[i];

		switch( methode){
		case 0: // LEAST SQUARES
			p = pow( e, 2);
			break;

		case 1: // HUBER
			if( fabs(e) <= k)
				p = 0.5 * pow( e, 2);
			else
				p = k*fabs(e) - 1/2.*k*k;
			break;

		case 2: // BISQUARE
			if( fabs(e) <= k)
				p = k*k/6.*(1-pow(1-pow(e/k,2),3));
			else
				p = k*k/6.;
			break;
		}

		sqrErr += w[i]*w[i]*p;
	}

	return sqrErr;
}

double robustRegression( double *x, double *y, int n, double *theta, int ntheta){
	double w[n];

	linregdata f_data;
	f_data.methode = 0;
	f_data.w = w;
	f_data.x = x;
	f_data.y = y;
	f_data.n = n;

	f_data.weight = 0;
	/*linregweigths( ntheta, theta, 0, &f_data);
	double lastS, newS = LevenbergMarquardt(
			ntheta, theta, 0, 0,
		1e-6, 0, 1e-6, 1000,
		linreg, &f_data
		);

	//fprintf( stderr, "-> after: %.2e, %.2e\n", theta[0], theta[1]);
	if(true)for( int i=0; i<100; i++){
		f_data.methode = 0;
		f_data.weight = 1;
		linregweigths( ntheta, theta, 0, &f_data);
		lastS = newS;
		newS = LevenbergMarquardt(
				ntheta, theta, 0, 0,
			1e-6, 10, 1e-6, 100,
			linreg, &f_data
			);
		//fprintf( stderr, "-> after: %.2e, %.2e\n", theta[0], theta[1]);

		if( fabs(newS-lastS) < 1e-6)
			return newS;
	}
	return newS;*/

//	return LevenbergMarquardt(
	return LineSearch(
				ntheta, theta, 0, 0,
			1e-6, 1, 1e-16, 1000,
			linregnoise, &f_data
			);
}

double robustRegression( double **x, double *y, int n, int d, double *theta, int ntheta){
	double w[n];

	multiregdata f_data;
	f_data.methode = 0;
	f_data.w = w;
	f_data.x = x;
	f_data.y = y;
	f_data.n = n;
	f_data.d = d;

	f_data.weight = 0;
	multiregweigths( ntheta, theta, 0, &f_data);
	double lastS, newS = LevenbergMarquardt(
		ntheta, theta, 0, 0,
		1e-6, 0, 1e-6, 1000,
		multireg, &f_data
		);
	//fprintf( stderr, "-> after: %.2e, %.2e\n", theta[0], theta[1]);
	for( int i=0; i<100; i++){
		f_data.methode = 0;
		f_data.weight = 1;
		multiregweigths( ntheta, theta, 0, &f_data);
		lastS = newS;
		newS = LevenbergMarquardt(
			ntheta, theta, 0, 0,
			1e-6, 0, 1e-6, 10,
			multireg, &f_data
			);
		//fprintf( stderr, "-> after: %.2e, %.2e\n", theta[0], theta[1]);

		if( fabs(newS-lastS) < 1e-6)
			return newS;
	}
	return newS;
}

double linearRegression( double *x, double *y, int n, double *theta, int ntheta)
{
	double w[n];
	linregdata f_data;

	f_data.w = w;
	f_data.x = x;
	f_data.y = y;
	f_data.n = n;

	f_data.methode = 0;
	f_data.weight = 0;

	linregweigths( ntheta, theta, 0, &f_data);

	return LevenbergMarquardt(
			ntheta, theta, 0, 0,
			1e-6, 0, 1e-6, 1000,
			linreg, &f_data
			);
}


#include "lik.ipp"
#include "cov.ipp"
#include "gp.ipp"


enum gpmode { MEAN, VAR, LOSS};

typedef struct{
	double* f;
	double* s2;
	double**X;
	int n;
	int d;
	double (*covFunc) ( double* &, double* &, int &, double* &);
	double* hyp;
	gpmode mode;
} gplikdata;
double gplik(int n, const double *hyp, double *grad, void *my_func_data)
{
	double *phyp = (double*)hyp;
	return - log_marginal_likelihood<double>(
			((gplikdata*)my_func_data)->f,
			((gplikdata*)my_func_data)->s2,
			((gplikdata*)my_func_data)->X,
			((gplikdata*)my_func_data)->n,
			((gplikdata*)my_func_data)->d,
			((gplikdata*)my_func_data)->covFunc,
			phyp);
}
double normalCFD(double x, double m, double s2 )
{
   return 0.5 * (1 + erf((x-m)/sqrt(2*s2)) );
}


double gpeval(int d, const double *x, double *grad, void *my_func_data)
{
	double y, s2;
	double **px = (double**)&x,
			*py = &y,
			*ps2= &s2;
	int n=1;

	double min;
	int    min_idx;
	findmin( ((gplikdata*)my_func_data)->f, ((gplikdata*)my_func_data)->n, min, min_idx);

	evalVariance<double>(
			((gplikdata*)my_func_data)->X,
			((gplikdata*)my_func_data)->f,
			((gplikdata*)my_func_data)->s2,
			((gplikdata*)my_func_data)->n,
			px,	py,	ps2, n,
			d,
			((gplikdata*)my_func_data)->covFunc,
			((gplikdata*)my_func_data)->hyp);

	switch( ((gplikdata*)my_func_data)->mode ){
	case VAR: // maximize variance
		return -s2;
	case LOSS: // maximize loss
		return log( 1 - normalCFD( min, y, s2) );
	case MEAN: // minimize mean
		return y;
	}
}



/*inline double mean( double *&x, int &n){
	double _mean=0;
	for( int j=0; j<n; j++)
		_mean += x[j];
	_mean /= (double)n;
}*/



#include <float.h>
void summaryStatistics( double (*func) (int, const double*, double*, void*),
		double *x, int dim, int n,
		double &mean, double &s2
	){
	mean = 0;
	s2 = 0;
	double y[n];
	for( int j=0; j<n; j++){
		double Y = (*func) ( dim, x, 0,    0);
		y[j] = Y;
		mean += Y;
		s2 += Y*Y;
	}

	mean /= (double) n;
	s2 = fmax( s2 / (double) n - mean*mean, 0); // VARIANCE
	s2 = s2 / (double) n; // VARIANCE OF MEAN
	//s2 = sqrt( s2 / (double) n); // STANDARD ERROR OF MEAN

	//mean = median(y,n);
}

void getAlphaBounds( double *x0, double *v, double *lb, double *ub, int dim, double &alpha_lb, double &alpha_ub){
	double alpha;
	alpha_lb =-DBL_MAX;
	alpha_ub = DBL_MAX;
	for( int d=0; d<dim; d++)
	if(v[d]!=0){
		alpha =  (lb[d] - x0[d]) / -v[d];
		if( alpha < 0)
			alpha_lb = fmax( alpha_lb, alpha);
		else
			alpha_ub = fmin( alpha_ub, alpha);

		alpha =  (ub[d] - x0[d]) / -v[d];
		if( alpha < 0)
			alpha_lb = fmax( alpha_lb, alpha);
		else
			alpha_ub = fmin( alpha_ub, alpha);
	}



	fprintf( stderr, "alpha %c [%.3e, %.3e]\n", 0xE2, alpha_lb, alpha_ub);
}


void getGradient( double (*func) (int, const double*, double*, void*), double *x, int dim, int maxit, double *&Mgrad){

	//int maxit = 100;

	double    y[maxit];
	double *grad[dim]; for( int d=0; d<dim; d++) grad[d] = (double*) malloc( maxit*sizeof(double));
	double   dx[dim], dx_min = 1e-6, dx_max = 0.1;
	double SEMy, My, SEMy_dx[dim], My_dx[dim];

	//int evals = 0;
	int it = 0;

	for( int d=0; d<dim; d++){
		dx[d] = 0.001;
	}

	while( it < maxit){
		// finite differences
		y[it] = (*func)( dim, x, 0, 0);
		stdmean( y, 0,0, it+1, SEMy, My);

		for( int d=0; d<dim; d++){
			x[d] += dx[d];
			grad[d][it] = ( y[it] - (*func)( dim, x, 0, 0) ) / dx[d];
			x[d] -= dx[d];

			//stdmean( y, grad[d], dx[d], it+1, SEMy_dx[d], My_dx[d]);
			stdmean( grad[d], 0, 0, it+1, SEMy_dx[d], My_dx[d]);
		}

		//mean( grad, dim, it+1, Mgrad);

		// uncertainty
		double uncertainty[dim];
		for( int d=0; d<dim; d++){
			uncertainty[d] = SEMy_dx[d]+SEMy - fabs( My-My_dx[d]) / 10.;

			if( uncertainty > 0){
				int m = ceil( pow( (SEMy_dx[d]+SEMy) / (My-My_dx[d]), 2 )*(it+1) );
				if( m>maxit && dx[d]<dx_max)
					dx[d] = dx[d]*2;

				if( m<maxit && dx[d]>dx_min)
					dx[d] = dx[d]/2;

				//Mgrad[d] = 0;
			}
		}

		// Breaking Criterion
		for( int d=0; d<dim; d++)
			if( it > 10 && uncertainty[d] <= 0){
				median( grad, dim, maxit, Mgrad);
				 for( d=0; d<dim; d++)
					 free( grad[d]);
				 fprintf( stderr, "finish grad. estim. after %i iterations\nSEMy_dx=%.2e + SEMy=%.2e - fabs( My=%.2e-My_dx=%.2e)", it+1,
						 SEMy_dx[d], SEMy, My, My_dx[d]);

				 return;
			}


		it++;
	}

	median( grad, dim, maxit, Mgrad);
	for( int d=0; d<dim; d++) free( grad[d]);
	fprintf( stderr, "finish grad. estim. after %i iterations\nSEMy_dx=%.2e + SEMy=%.2e - fabs( My=%.2e-My_dx=%.2e)", it+1,
		SEMy_dx[0], SEMy, My, My_dx[0]);
}


#include "gp.hpp"
void gplinesearch( double *x0, double *lb, double *ub, int dim, double (*func) (int, const double*, double*, void*), optimoptions *options, double *sol)
{
/*	{
		// create noisy data
		unsigned int p_seed = 0;
		int _n=500;
		double _x[_n], _y[_n], _xtest;
		for( int i=0; i<_n; i++){
			double rnd;
			BoxMuller( &rnd, 1, &p_seed);
			_x[i] = i/(_n+1.);
			_y[i] = _x[i]*1e1+3 + exp(rnd*100);
		}

		int order = 3;
		double linregpar[order+1]; for( int i=0; i<order+1; i++)linregpar[i]=0;
		robustRegression( _x, _y, _n, linregpar, order+1);
		double Mgrad = dpolynom( _xtest, linregpar, order);

		char filename[1024];
		sprintf( filename, "robustreg.dat");
		FILE* fp = fopen( filename, "w+");
		for( int i=0; i<_n; i++){
			fprintf( fp, "%.12e %.12e %.12e\n", _x[i], _y[i], polynom( _x[i], (double*)linregpar, order));
		}
		fclose( fp);
	}
exit(0);*/

	double solY;

	// TEST END
	unsigned int seed = 0;
	double S[dim]; for( int i=0; i<dim; i++) S[i] = 1;
	double Smean = 1;

	// DECLARE VARIABLES
	double **popX = allocMatrix<double>( options->MaxFunEvals, dim);
	double  *popY = (double*)malloc( options->MaxFunEvals*sizeof(double));
	double  *popS2= (double*)malloc( options->MaxFunEvals*sizeof(double));
	int      popSize = 0;

	FILE *fp_pop = fopen( "population.out", "w+");
	FILE *fp_pred = fopen( "prediction.out", "w+");

	// INITIAL SAMPLING
	popSize = 1;
//	LHS( popX, lb, ub, popSize, dim);
//	for( int j=0; j<dim; j++)
//		sol[j] = popX[0][j];
	if( x0 == 0)
		for( int j=0; j<dim; j++)
			sol[j] = popX[0][j] = (lb[j]+ub[j]) / 2.;
	else
		for( int j=0; j<dim; j++)
			sol[j] = popX[0][j] = x0[j];


	// EVALUATE MEAN & STD
	for( int i=0; i<popSize; i++){
		summaryStatistics( func, popX[i], dim, options->MaxFunEvalsAvg, popY[i], popS2[i]);

		fprintf( fp_pop, "0 ");
		for( int j=0; j<dim; j++)
			fprintf( fp_pop, "%e ", popX [i][j]);
		fprintf( fp_pop, "%e %e\n", popY [i],popS2 [i]);

		for( int j=0; j<dim; j++)
			fprintf( stderr, "%e ", popX [i][j]);
		fprintf( stderr, "mean=%e,  s2=%e\n", popY [i], popS2[i]);
	}
	fclose(fp_pop);

	solY = popY[0];

	for( int iter=0; iter<options->MaxIter; iter ++){

		fprintf(stderr, "iteration: %i -> ", iter);

		// ESTIMATE GRADIENT
		double Mgrad[dim], *pMgrad = Mgrad;
		switch( 5){
		case 3:
			for( int i=0; i<dim; i++)
				Mgrad[i] = ( iter % dim == i ? 1. : 0) * ((iter/dim) % 2 ? -1 : 1);
			break;

		case 1:
			getGradient( func, sol, dim,
					//10+iter/2,
					100,
					pMgrad); //exit(0);
			break;

		case 2:{
			FILE *fp;

			for( int d=0; d<dim; d++){

				double dx = 0.5e-0;
				int    _n = 500;
				double _x[_n];
				double _y[_n];
				double weight = fmin( Smean/S[d], 1);
				for( int i=0; i<_n; i++){
					//double perturbe = ( i<_n/3 ? -dx : (i<2*_n/3 ? 0 : dx)) * weight ;
					double perturbe = dx*( 2.*i/(_n-1.) - 1) * weight ;
					//double perturbe = dx*(i-_n/2. < 0 ? -1 : 1)*pow( 1.1, fabs( i-_n/2.));

					sol[d] +=perturbe;
					_x[i] = sol[d];
					_y[i] = (*func)( dim, sol, 0,0);
					sol[d] -= perturbe;

					//fprintf( fp, "%e %e\n", _x[i], _y[i]);
				}
				int order = 2;
				double linregpar[order+1]; for( int i=0; i<order+1; i++)linregpar[i]=0;
				//linreg_firstguess( order, linregpar, _x,_y,_n);
				S[d] = robustRegression( _x, _y, _n, linregpar, order+1);
				Mgrad[d] = dpolynom( sol[d], linregpar, order);
				/*linreg_firstguess( order, linregpar, _x,_y,2*_n/3);
				S[d] = robustRegression( _x, _y, 2*_n/3, linregpar, order);
				Mgrad[d] = dpolynom( sol[d], linregpar, order);

				linreg_firstguess( order, linregpar, &_x[_n/3-1],&_y[_n/3-1],2*_n/3);
				S[d] = robustRegression( &_x[_n/3-1],&_y[_n/3-1],2*_n/3, linregpar, order);
				Mgrad[d] += dpolynom( sol[d], linregpar, order);*/

				char filename[1024];
				sprintf( filename, "linreg%i.dat", d);
				fp = fopen( filename, "w+");
				for( int i=0; i<_n; i++){
					fprintf( fp, "%.12e %.12e %.12e\n", _x[i], _y[i], polynom( _x[i], (double*)linregpar, order));
				}
				fclose( fp);
			}

			Smean = mean( S, dim);
			//for( int d=0; d<dim; d++)
			//	Mgrad[d] *= fmin( Smean/S[d], 1);


		}break;

		case 4:{
			// SAMPLE


			double dx = 1e-2;
			int    _n = 1000;
			double **_x = allocMatrix<double>(_n, dim);
			double _y[_n];
			for( int i=0; i<_n; i++){
				for( int d=0; d<dim; d++){
					_x[i][d] = sol[d] + dx*(unifrnd<double>( &seed)*2-1 );
				}
				_y[i] = (*func)( dim, _x[i], 0,0);
			}

			int order = dim+1;
			double multiregpar[order];//, linregpar_old[2];
			for( int i=0; i<order; i++) multiregpar[i] = 0;
			robustRegression( _x, _y, _n, dim, multiregpar, order);


			FILE* fp_grad = fopen( "multireg.dat", "w+");
			for( int i=0; i<_n; i++){
				for( int d=0; d<dim; d++){
					fprintf( fp_grad, "%e ", _x[i][d]);
				}
				fprintf( fp_grad, "%e %e\n", _y[i], plane( _x[i], dim, multiregpar, order));
			}
			fclose(fp_grad);


			freeMatrix( _x, _n);

			dplane( sol, dim, multiregpar, order, Mgrad);
		}break;

		case 5:{ 			// Hyper-cube sampling & Gaussprocess approximation
			int _n=500;
			double **_x = allocMatrix<double>( _n, dim );
			double _lb[_n];
			double _ub[_n];
			double _y[_n];
			for( int i=0; i<_n; i++){
				_lb[i] = sol[i] - 1e-1;
				_ub[i] = sol[i] + 1e-1;
			}
			LHS<double>( _x, _lb, _ub, _n, dim);
			for( int i=0; i<_n; i++){
				_y[i] = (*func)( dim, _x[i], 0,0);
			}
			double hyp[4];
			getHyperParameters( _x, _y, _n, dim, hyp);
			hyp[1] += 2;
			hyp[0] = -1;
			//hyp[2] = 1;

			double *pMgrad = Mgrad;
			evalGradientSparse<double>(
					_x, _y,   0,  _n,
					&sol, &pMgrad, 1,
					dim,
					covMatern5, hyp);

			double _yt[_n];
			evalVarianceSparse<double>(
			//evalVariance(
					_x, _y,  0,  _n,
					_x, _yt, 0,   _n,
					dim,
					covMatern5, hyp);

			char filename[1024];
			FILE *fp;
			for( int d=0; d<dim; d++){
				sprintf( filename, "linreg%i.dat", d);
				fp = fopen( filename, "w+");
				for( int i=0; i<_n; i++){
					fprintf( fp, "%.12e %.12e %.12e %.12e\n", _x[i][d], _y[i], (_x[i][d]-sol[d])*pMgrad[d] + solY, _yt[i]);
				}
				fclose( fp);
			}

			freeMatrix(_x, _n );
		}break;
		}

		// NORMALIZE GRADIENT
		double abs_grad = 0;
		for( int i=0; i<dim; i++)
			abs_grad += Mgrad[i]*Mgrad[i];
		abs_grad = sqrt(abs_grad);

		FILE *fp_grad = fopen( "grad.dat", "w+");
		fprintf( stderr, "grad [");
		for( int i=0; i<dim; i++){
			fprintf( stderr, "%.2e (%.2e) ", Mgrad[i], Mgrad[i]/abs_grad);
			Mgrad[i]/=abs_grad;
			fprintf(fp_grad, "%e %e ", sol[i], Mgrad[i]);
		}
		fprintf( stderr, "] ");
		fclose(fp_grad);

		switch(1){
		case 1:{   // GAUSSIAN PROCESS LINESEARCH

			//boundaries
			int    alpha_size= 30;
			double alpha_min;
			double alpha_max;
			getAlphaBounds( sol, Mgrad, lb, ub, dim, alpha_min, alpha_max);
			double alpha_x[alpha_size], *p_alpha_x = alpha_x;
			double 				*p_alpha_y = &popY[popSize];
			double 			    *p_alpha_s2= &popS2[popSize];
			alpha_min = 0;
			alpha_x[0] = alpha_min;
			alpha_x[1] =(alpha_min+alpha_max) / 2.;
			alpha_x[2] = alpha_max;
			// init
			for( int i=0; i<3; i++){
				// set new point
				for( int d=0; d<dim; d++){
					popX[popSize+i][d] = sol[d] - Mgrad[d]*alpha_x[i];
					//fprintf( fp_gpline, "%e ", popX[popSize+i][d]);
				}

				// evaluate summary stats
				summaryStatistics( func, popX[popSize+i], dim, options->MaxFunEvalsAvg, popY[popSize+i], popS2[popSize+i]);
			}
			//popY[popSize] = solY;

			for( int gp_it=3; gp_it<alpha_size; gp_it++){
				// test points
				int    alpha_test_size = 100;
				double alpha_test_x[alpha_test_size], *p_alpha_test_x=alpha_test_x;
				double alpha_test_y[alpha_test_size], *p_alpha_test_y=alpha_test_y;
				double alpha_test_s2[alpha_test_size], *p_alpha_test_s2=alpha_test_s2;

				for( int i=0; i<alpha_test_size; i++){
					alpha_test_x[i] = alpha_min + (alpha_max-alpha_min)*i/(alpha_test_size-1);
				}

				// predict test points with GP
				double hyp[4], *p_hyp=hyp; int dim_alpha=1;
				getHyperParameters( p_alpha_x, p_alpha_y, gp_it, hyp);
				evalVariance<double>(
							p_alpha_x,	    p_alpha_y,	    p_alpha_s2,	     gp_it,
							p_alpha_test_x,	p_alpha_test_y,	p_alpha_test_s2, alpha_test_size,
							dim_alpha,
							covMatern5,
							p_hyp);

				// find minimum of GP
				int idx;
				double min_y;
				if( gp_it%2)
					findmin( p_alpha_test_y, alpha_test_size, min_y, idx);
				else
					findmax( p_alpha_test_s2, alpha_test_size, min_y, idx);
				alpha_x[gp_it] = p_alpha_test_x[idx];
				fprintf( stderr, " - - - - > min(x=%f, i=%i) = %f\n", p_alpha_test_x[idx], idx, p_alpha_test_y[idx]);

				// evaluate model
				for( int d=0; d<dim; d++){
					popX[popSize+gp_it][d] = sol[d] - Mgrad[d]*alpha_x[gp_it];
					//fprintf( fp_gpline, "%e ", popX[popSize+gp_it][d]);
				}
				summaryStatistics( func, popX[popSize+gp_it], dim, options->MaxFunEvalsAvg, popY[popSize+gp_it], popS2[popSize+gp_it]);
			}

			char filename[512];
			//sprintf( filename, "gpline%i.dat", iter);
			sprintf( filename, "gpline.dat");
						FILE *fp_gpline = fopen( filename, "w+");
			for( int i=0; i<alpha_size; i++){
				for( int d=0; d<dim; d++)
					fprintf( fp_gpline, "%e ", popX[popSize+i][d]);
				fprintf( fp_gpline, "%e %e %e\n", alpha_x[i], p_alpha_y[i], p_alpha_s2[i]);
			}
			fclose( fp_gpline);


			int idx;
			double min_y;
			findmin( p_alpha_y, alpha_size, min_y, idx);
			fprintf( stderr, "linsearch result: x[%i] = %e -> f = %e \n", idx, p_alpha_x[idx], p_alpha_y[idx]);
			for( int d=0; d<dim; d++){
				popX[popSize][d] = sol[d] - Mgrad[d]*alpha_x[idx];
				popX[popSize][d] = fmax( popX[popSize][d], lb[d]+1e-3);
				popX[popSize][d] = fmin( popX[popSize][d], ub[d]-1e-3);
			}
			popY[popSize] = p_alpha_y[idx];
			popS2[popSize]= p_alpha_s2[idx];

		}  break;

		case 2:{   // LINESEARCH
			double alpha = 1;
			double alpha_min = 1e-6;
			// first guess
			for( int i=0; i<dim; i++)
				popX[popSize][i] = sol[i] - Mgrad[i]*alpha;

			// borders respected ?
			while( !inbound( popX[popSize], lb, ub, dim) && alpha > alpha_min){
				alpha /= 2;
				for( int i=0; i<dim; i++)
					popX[popSize][i] = sol[i] - Mgrad[i]*alpha;
			}

			// improvement ?
			summaryStatistics( func, popX[popSize], dim, options->MaxFunEvalsAvg, popY[popSize], popS2[popSize]);
			while( popY[popSize] > solY && alpha > alpha_min){
				alpha /= 2;
				for( int i=0; i<dim; i++)
					popX[popSize][i] = sol[i] - Mgrad[i]*alpha;
				summaryStatistics( func, popX[popSize], dim, options->MaxFunEvalsAvg, popY[popSize], popS2[popSize]);
			}
		}
		}

		for( int d=0; d<dim; d++){
			fprintf( stderr, "%e ", popX[popSize][d]);
		}

		fp_pop = fopen( "population.out", "a+");
		fprintf( fp_pop, "%i ", iter);
		for( int i=0; i<dim; i++)
			fprintf( fp_pop, "%e ", popX [popSize][i]);
		fprintf( fp_pop, "%e %e\n", popY [popSize],popS2 [popSize]);
		//fprintf( fp, "\n\n");
		fclose(fp_pop);



		//popY [popSize] = (*func) ( dim, popX[popSize], 0,    0);
		fprintf(stderr, " -> f=%.3e, s2=%.3e\n", popY[popSize], popS2[popSize]);

		//if( popY[popSize] < solY)
		{
			solY = popY[popSize];
			fprintf(stderr, "Solution (%e) is ", solY);
			for( int i=0; i<dim; i++){
				sol[i] = popX[popSize][i];
				fprintf(stderr, "%e ", sol[i]);
			}
			fprintf(stderr, "\n");
		}


		popSize++;

		// Contract around best
		/*int idx;
		double solf;
		findmin( popY, popSize, solf, idx);
		fprintf(stderr, "Solution (%e) is ", solf);
		for( int i=0; i<dim; i++){
			sol[i] = popX[idx][i];
			fprintf(stderr, "%e ", sol[i]);
		}
		fprintf(stderr, "\n");*/

		/*double beta=0.99;
		for( int i=0; i<dim; i++){
			lb[i] = lb[i]*beta + sol[i]*(1-beta);
			ub[i] = ub[i]*beta + sol[i]*(1-beta);
		}*/

	}

	//fclose(fp_pop);
	fclose(fp_pred);
}


void gpopt( double *x0, double *lb, double *ub, int dim, double (*func) (int, const double*, double*, void*), optimoptions *options, double *sol)
	{
	double **popX = allocMatrix<double>( options->MaxFunEvals, dim);
	double  *popY = (double*)malloc( options->MaxFunEvals*sizeof(double));
	double  *popS2= (double*)malloc( options->MaxFunEvals*sizeof(double));
	int      popSize = 0;
	double   minX[dim];
	double solf;
	unsigned int seed = 0;

	double hyp[4];
	double hlb[4] = {2,  0,  -1, 1e5};
	double hub[4] = {2,  3,  7,  1e5};

	int evals = 0;
	FILE *fp;
	FILE *fp_pred = fopen( "prediction.out", "w+");

	// INITIAL SAMPLING
	popSize = dim*3;
	LHS( popX, lb, ub, popSize, dim);
	for( int j=0; j<dim; j++)
		sol[j] = popX[0][j];

	// EVALUATE MEAN & STD
	fp = fopen( "population.out", "w+");
	for( int i=0; i<popSize; i++){
		summaryStatistics( func, popX[i], dim, options->MaxFunEvalsAvg, popY[i], popS2[i]);

		fprintf( fp, "0 ");
		for( int j=0; j<dim; j++)
			fprintf( fp, "%e ", popX [i][j]);
		fprintf( fp, "%e %e\n", popY [i],popS2 [i]);

		for( int j=0; j<dim; j++)
			fprintf( stderr, "%e ", popX [i][j]);
		fprintf( stderr, "mean=%e,  s2=%e\n", popY [i], popS2[i]);
	}
	fclose(fp);

	for( int iter=0; iter<options->MaxIter; iter ++){

		fprintf(stderr, "iteration: %i -> ", iter);

		// TRAIN HYPER PARAMETERS
		//double (*covFunc) ( double* &, double* &, int &, double* &) = covMatern5<double>;
		//log_marginal_likelihood<double>( popY, popS2, popX, popSize, dim, covFunc, phyp);
		//fprintf(stderr, "TRAIN HYPER PARAMETERS\n");
		gplikdata gpdata;
			gpdata.X = popX;
			gpdata.f = popY;
			gpdata.s2= popS2;
			gpdata.n = popSize;
			gpdata.d = dim;
			gpdata.covFunc = covMatern5<double>;
			gpdata.hyp = hyp;

		if(false){
			int nsols=10;
			int dhyp=4;
			double **sols = allocMatrix<double>( nsols, dhyp);
			double **hyp0 = allocMatrix<double>( nsols, dhyp);
			double   lik[nsols];
			LHS( hyp0, hlb, hub, nsols, dhyp);

			optimoptions hoptions = getoptions();
				hoptions.FinDiffRelStep = 1e-6;

			multistart( hyp0, hlb, hub, dhyp, gplik, (void*)&gpdata, &hoptions, gradientDecent, sols, nsols);
		//	fprintf( stderr, "hyp = [%e %e %e %e], lik=%e\n", hyp[0], hyp[1], hyp[2], hyp[3],
		//			log_marginal_likelihood<double>( popY, popS2, popX, popSize, dim, covFunc, phyp));
			for( int i=0; i<nsols; i++){
				lik[i] = -log_marginal_likelihood<double>( popY, popS2, popX, popSize, dim, covMatern5<double>, sols[i]);
				//for( int d=0; d<dhyp; d++)
				//	fprintf(stderr, "%e ", sols[i][d]);
				//fprintf(stderr, " -> %e\n", lik[i]);
			}
			double min;
			int idx;
			findmin( lik, nsols, min, idx);
			fprintf(stderr, "hyper parameters: %i -> ", idx);
			for( int i=0; i<dhyp; i++){
				hyp[i] = sols[idx][i];
				fprintf(stderr, "%8.2f ", sols[idx][i]);
			}
			fprintf(stderr, " -> liklihood: %.3e  ", min);

		}else{
			getHyperParameters( popX, popY, popSize, dim, hyp);
			fprintf(stderr, "hyper parameters -> ");
			for( int i=0; i<4; i++){
				fprintf(stderr, "%8.2e ", hyp[i]);
			}

		}


		// SAMPLE FROM GP
		//if( false)
		if( dim < 3){
			int x[2];
			int nx[2] = {100, 100};
			int      new_popSize = nx[0]*nx[1];
			double **new_popX = allocMatrix<double>( new_popSize, dim);
			double  *new_popY = (double*)malloc( new_popSize*sizeof(double));
			double  *new_popS2= (double*)malloc( new_popSize*sizeof(double));

			new_popSize=0;
			for(  x[0]=0; x[0]<nx[0]; x[0]++)
				for(  x[1]=0; x[1]<nx[1]; x[1]++){
					for( int i=0; i<dim; i++)
						new_popX[new_popSize][i] = lb[i] + x[i]/(double)(nx[i]-1)*(ub[i]-lb[i]);
					new_popY[new_popSize]=0;
					new_popS2[new_popSize]=0;
					new_popSize++;
				}
			evalVariance<double>(
					popX,    	popY,        popS2,        popSize,
					new_popX,  	new_popY,    new_popS2,    new_popSize,
					dim,
					covMatern5<double>, gpdata.hyp);

			new_popSize=0;
			for(  x[0]=0; x[0]<nx[0]; x[0]++){
				for(  x[1]=0; x[1]<nx[1]; x[1]++){
					fprintf( fp_pred, "%i ", iter );
					for( int i=0; i<dim; i++)
						fprintf( fp_pred, "%e ", new_popX[new_popSize][i]);
					fprintf( fp_pred, "%e %e %e\n", new_popY[new_popSize], new_popS2[new_popSize], normalCFD( solf, new_popY[new_popSize], new_popS2[new_popSize]));
					new_popSize++;
				}
				if( dim>1)
				fprintf( fp_pred, "\n");
			}
			fprintf( fp_pred, "\n");
			if( dim==1)
			fprintf( fp_pred, "\n");
		}

		//fprintf(stderr, "PREDICT MINIMUM\n");
		//double minX[2];
		optimoptions gpoptions = getoptions();
		gpoptions.MaxFunEvals = 1000;
		gpoptions.Display = false;

		double gp_lb[dim];
		double gp_ub[dim];
		for( int i=0; i<dim; i++){	gp_lb[i] = lb[i]; gp_ub[i] = ub[i]; }

		int mode = iter%(1+dim);
		switch( mode){
		case 0:
			gpdata.mode = LOSS;
			fprintf(stderr, "predicted loss: ");
			for( int i=0; i<dim; i++){	minX[i] = gp_lb[i] + (gp_ub[i]-gp_lb[i])*unifrnd<double>( &seed);}
			break;
		case 1:
			gpdata.mode = VAR;
			fprintf(stderr, "predicted maximum: ");
			for( int i=0; i<dim; i++){	minX[i] = gp_lb[i] + (gp_ub[i]-gp_lb[i])*unifrnd<double>( &seed);}
			break;
			//gpoptions.Display = true;
		default:
			gpdata.mode = VAR;
			double beta = 1. /  (double)mode;
			fprintf(stderr, "predicted max (%.2f): ", beta);
			for( int i=0; i<dim; i++){	gp_lb[i] = beta*gp_lb[i] + (1-beta)*sol[i]; gp_ub[i] = beta*gp_ub[i] + (1-beta)*sol[i]; }
			for( int i=0; i<dim; i++){	minX[i] = sol[i];}
			break;
		}

		//for( int i=0; i<dim; i++)
		//	minX[i] = gp_lb[i] + (gp_ub[i]-gp_lb[i])*unifrnd(double, &seed);

		//double *pminX = (double*)minX;
		//multistart( 0, lb, ub, dim, gpeval, (void*)&gpdata, &gpoptions, gradientDecent, &pminX, 1);
		gradientDecent( minX, gp_lb, gp_ub, dim, gpeval, (void*)&gpdata, &gpoptions, minX);
		// cross check
		/*double minY = gpeval( dim, minX, 0, (void*)&gpdata);
		for( int j=0; j<1000; j++){
			double tempX[dim];
			for( int i=0; i<dim; i++)
				tempX[i] = lb[i] + (ub[i]-lb[i])*unifrnd(double, &seed);
			double tempY = gpeval( dim, tempX, 0, (void*)&gpdata);
			if( tempY < minY){
				fprintf( stderr, "ERROR!!! %e < %e\n", tempY, minY);
				exit(0);
			}

		}*/

		//LevenbergMarquardt(dim, minX, gp_lb, gp_ub, 1e-3, 1, 1e-30, 100, gpeval, (void*)&gpdata);

		for( int i=0; i<dim; i++)
			fprintf(stderr, "%5.3f ", minX[i]);
		//fprintf(stderr, "\n");

		switch(2){
		case 1:
			// ADD NEW POINT
			//fprintf(stderr, "EVALUATE MINIMUM\n");
			for( int i=0; i<dim; i++)
				popX[popSize][i] = minX[i];

			summaryStatistics( func, popX[popSize], dim, options->MaxFunEvalsAvg, popY[popSize], popS2[popSize]);

			fp = fopen( "population.out", "a+");
			fprintf( fp, "%i ", iter);
			for( int i=0; i<dim; i++)
				fprintf( fp, "%e ", popX [popSize][i]);
			fprintf( fp, "%e %e\n", popY [popSize],popS2 [popSize]);
			fprintf( fp, "\n\n");
			fclose(fp);

			//popY [popSize] = (*func) ( dim, popX[popSize], 0,    0);
			fprintf(stderr, " -> f=%.3e, s2=%.3e\n", popY[popSize], popS2[popSize]);

			popSize++;

			break;

		case 2:{
			double beta = 0.1;
			int n=10;
			for( int i=0; i<dim; i++){
				gp_lb[i] = beta*lb[i] + (1-beta)* sol[i] - (ub[i]-lb[i])*beta;
				gp_ub[i] = beta*lb[i] + (1-beta)* sol[i] + (ub[i]-lb[i])*beta;
			}
			LHS( &popX[popSize], gp_lb, gp_ub, n, dim, &seed);

			fp = fopen( "population.out", "a+");
			for( int i=0; i<n; i++){
				summaryStatistics( func, popX[popSize+i], dim, options->MaxFunEvalsAvg, popY[popSize+i], popS2[popSize+i]);

				fprintf( fp, "%i ", iter);
				for( int j=0; j<dim; j++)
					fprintf( fp, "%e ", popX [popSize+i][j]);
				fprintf( fp, "%e %e\n", popY [popSize+i],popS2 [popSize+i]);
			}
			fprintf( fp, "\n\n");
			fclose(fp);
			popSize+=n;
			}break;

		}


		// Contract around best
		int idx;
		findmin( popY, popSize, solf, idx);
		for( int i=0; i<dim; i++)
			sol[i] = popX[idx][i];

		/*double beta=0.99;
		for( int i=0; i<dim; i++){
			lb[i] = lb[i]*beta + sol[i]*(1-beta);
			ub[i] = ub[i]*beta + sol[i]*(1-beta);
		}*/

		fprintf(stderr, "predicted minimum: \n");
		gpdata.mode = MEAN;
		gradientDecent( sol, gp_lb, gp_ub, dim, gpeval, (void*)&gpdata, &gpoptions, minX);
		for( int i=0; i<dim; i++)
			fprintf(stderr, "%5.3f ", minX[i]);
		fprintf(stderr, "\n");
	}

	// SOLUTION == PREDICTION
	for( int i=0; i<dim; i++)
		sol[i] = minX[i];

	fclose(fp_pred);
}

/*void LevenbergMarquardt(
	int sampleSize, double *sampleX, double *sampleY,
	int parameterSize, double *parameters, double *parametersMin, double *parametersMax,
	double differentiationStepSize, double lambda, double minSquares, int maxIterations,
	void (*f)(int, double*, double*, int, double*)
	)
{
	fprintf(stderr, "Start LevenbergMarquardt\n");

	double S, S_old = 0.;
	int it=0;
	double **jacobi = newDoubleMatrix(sampleSize, parameterSize);
	double fitY[sampleSize];

	double **A = newDoubleMatrix(parameterSize, parameterSize);
	double b[parameterSize];
	double x[parameterSize];

	// COPY PARAMETERS
	double log_parameters[parameterSize];
	double log_parametersMin[parameterSize];
	double log_parametersMax[parameterSize];
	for (int i = 0; i < parameterSize; i++) {
		setParameter( log_parameters,    i, parameters[i]);
		setParameter( log_parametersMin, i, parametersMin[i]);
		setParameter( log_parametersMax, i, parametersMax[i]);
	}


	// get function sample
	(*f)( sampleSize, sampleX, fitY, parameterSize, log_parameters);

	// compare data to fit
	S = SumOfSquares( sampleSize, sampleY, fitY);
	fprintf(stderr, "0 iterations: S = %e -> step size = %e \b", S, differentiationStepSize);



	for( ;  it<maxIterations  && S>minSquares; it++)
	{

		// Construct Jacobi-Matrix
		getJacobiMatrix(sampleSize, sampleX, sampleY, parameterSize, log_parameters, jacobi, f, differentiationStepSize);

		// Construct Linear System
		for (int i = 0; i < parameterSize; i++) {

			// Levenberg
			for (int j = 0; j < parameterSize; j++) {
				A[i][j] = 0.;
				for (int k = 0; k < sampleSize; k++)
					// A = J^T J
					A[i][j] += jacobi[k][i] * jacobi[k][j];
			}
			// Marquardt
			A[i][i] *= (1. + lambda);

			b[i] = 0.;
			for (int k = 0; k < sampleSize; k++)
				b[i] -= jacobi[k][i] * (sampleY[k] - fitY[k]);
		}

		// Solve Linear System
		solveLinearSystem(A, b, x, parameterSize);

		// Update Parameters
		for (int i = 0; i < parameterSize; i++)
			log_parameters[i] -= x[i];

		// Check upper and lower bounds
		for( int i=0; i<parameterSize; i++){
			if( log_parameters[i] < log_parametersMin[i]) log_parameters[i] = log_parametersMin[i];
			if( log_parameters[i] > log_parametersMax[i]){
				log_parameters[i] = log_parametersMax[i];
				//fprintf( stderr, "UPPER BOUND: log(%e) = %e\n", exp(log_parameters[i]), log_parameters[i]);
			}
		}


		// get function sample
		(*f)( sampleSize, sampleX, fitY, parameterSize, log_parameters);

		// squares
		S_old = S;
		S = SumOfSquares( sampleSize, sampleY, fitY);
		fprintf(stderr, "\r%i iterations: S = %e -> step size = %e \b", it, S, differentiationStepSize);
	}

	// COPY BACK PARAMETERS
	for (int i = 0; i < parameterSize; i++) {
		parameters[i]    = getParameter( log_parameters, i);
		parametersMin[i] = getParameter( log_parametersMin, i);
		parametersMax[i] = getParameter( log_parametersMax, i);
	}



	deleteDoubleMatrix( A, parameterSize);
	deleteDoubleMatrix( jacobi, sampleSize);

	fprintf(stderr, "\nParameters: ");
	for( int i=0; i<parameterSize; i++)
		fprintf(stderr, "%10.3e ", parameters[i]);
	fprintf(stderr, "\n");
}*/




