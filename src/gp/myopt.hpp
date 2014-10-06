/*
 * myopt.hpp
 *
 *  Created on: Sep 2, 2014
 *      Author: jagiella
 */

#ifndef MYOPT_HPP_
#define MYOPT_HPP_

// DATA STRUCTURES

typedef struct {
	bool   Display;
	double FinDiffRelStep;
	double TolFun;
	double TolX;
	int    MaxIter;
	int    MaxFunEvals;
	bool   GradObj;
	int    MaxFunEvalsAvg;
} optimoptions;

// RANDOM SAMPLINGS

template <class T>
void LHS( T **x, int n, int d);

template <class T>
void LHS( T **x, T* lb, T* ub, int n, int d);

// OPTIMIZER

optimoptions getoptions();
void gradientDecent( double *x0, double *lb, double *ub, int dim, double (*func) (int, const double*, double*, void*), void*, optimoptions *options, double *sol);
void multistart(
		// input
		double **X0, double *lb, double *ub, int dim, double (*func) (int, const double*, double*, void*), void*,
		// options
		optimoptions *options,
		// methode
		void (*methode) (double*,double*,double*,int,double (*) (int, const double*, double*, void*),void*, optimoptions*,double*),
		// output
		double **sol, int n);
void gpopt( double *x0, double *lb, double *ub, int dim, double (*func) (int, const double*, double*, void*), optimoptions *options, double *sol);


#include "myopt.ipp"

#endif /* MYOPT_HPP_ */
