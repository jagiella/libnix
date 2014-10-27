/*
 * myopt.hpp
 *
 *  Created on: Sep 2, 2014
 *      Author: jagiella
 */

#ifndef MYOPT_HPP_
#define MYOPT_HPP_


// RANDOM SAMPLINGS

template <class T>
void LHS( T **x, int n, int d, unsigned int *seed = 0);

template <class T>
void LHS( T **x, T* lb, T* ub, int n, int d, unsigned int *seed = 0);

// OPTIMIZER



#include "myopt.ipp"

#endif /* MYOPT_HPP_ */
