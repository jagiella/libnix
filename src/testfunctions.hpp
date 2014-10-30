/*
 * testfunctions.hpp
 *
 *  Created on: Oct 29, 2014
 *      Author: jagiella
 */

#ifndef TESTFUNCTIONS_HPP_
#define TESTFUNCTIONS_HPP_

double myfunc(int n, const double *x, double *grad, void *my_func_data);

double myfunc2(int n, const double *x, double *grad, void *my_func_data);

double myfuncnorm(int n, const double *x, double *grad, void *my_func_data);

double myfuncnorm2(int n, const double *x, double *grad, void *my_func_data);

double myfuncbinom(int n, const double *x, double *grad, void *my_func_data);





double myfunc3(int n, const double *x, double *grad, void *my_func_data);

double externalfuncPipe(int n, const double *x, double *grad, void *my_func_data);

double externalfuncFile(int n, const double *x, double *grad, void *my_func_data);

double parabol(int n, const double *x, double *grad, void *my_func_data);

#endif /* TESTFUNCTIONS_HPP_ */
