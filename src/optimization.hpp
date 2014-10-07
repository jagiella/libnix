/*
 * optimization.hpp
 *
 *  Created on: Oct 7, 2014
 *      Author: jagiella
 */

#ifndef OPTIMIZATION_HPP_
#define OPTIMIZATION_HPP_


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


#endif /* OPTIMIZATION_HPP_ */
