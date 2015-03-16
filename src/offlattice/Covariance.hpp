/*
 * Covariance.hpp
 *
 *  Created on: 07.11.2013
 *      Author: jagiella
 */

#ifndef COVARIANCE_HPP_
#define COVARIANCE_HPP_


void gsl_covariance_matrix(double **A, int m, int n, double **covA);
void gsl_covariance_matrix(double **A, int m, int n, double *covA);



#endif /* COVARIANCE_HPP_ */
