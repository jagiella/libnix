/*
 * abc.hpp
 *
 *  Created on: Oct 29, 2014
 *      Author: jagiella
 */

#ifndef ABC_HPP_
#define ABC_HPP_


void abc(
			// input
			double **X0, double *lb, double *ub, int dim, double (*func) (int, const double*, double*, void*), void*,
			// options
			optimoptions *options,
			// methode
			void (*methode) (double*,double*,double*,int,double (*) (int, const double*, double*, void*),void*, optimoptions*,double*),
			// output
			double **sol, int n);


#endif /* ABC_HPP_ */
