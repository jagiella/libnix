/*
 * statistics.hpp
 *
 *  Created on: Oct 7, 2014
 *      Author: jagiella
 */

#ifndef STATISTICS_HPP_
#define STATISTICS_HPP_



#define unifrnd( T, seed) ((T)rand_r(seed)/(T)(RAND_MAX))
void    normrnd( double *rnd, int n, unsigned int *p_seed);
double  normrnd( unsigned int *p_seed);



typedef struct{
	bool   H;
	double pValue;
	double KSstatistic;
} KSTestResult;

KSTestResult KolmogorovSmirnoffTest( int n1, double *x1, int n2, double *x2);

double logLikelihood( double *v1, double *v2, double *s21, double *s22, int n);



enum VectorOrientation {rowwise = 1, columnwise = 2};
double *mean( double **x, int n, int d, VectorOrientation ori);
void    mean( double **x, int d, int n, double *&mean_x);
double  mean( double  *x, int n);
double *std2(   double **x, double *m, int n, int d, VectorOrientation orientation);
void    stdmean( double *x, double *dx, double d, int n, double &SEM_x, double &mean_x);
double median( double  *x, int n);
void   median( double **x, int d, int n, double *&median_x);
double mad(    double  *x, int n);
void findmin( double *x, int n, double &o_min, int &o_idx);
void findmax( double *x, int n, double &max, int &idx);
double  min( double  *x, int n);
double  max( double  *x, int n);

#endif /* STATISTICS_HPP_ */
