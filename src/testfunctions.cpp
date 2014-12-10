/*
 * testfunctions.cpp
 *
 *  Created on: Oct 29, 2014
 *      Author: jagiella
 */

#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <omp.h>

#include "statistics.hpp"
#include "model.hpp"

double myfunc(int n, const double *x, double *grad, void *my_func_data)
{
	if (grad) {
        grad[0] = 0.0;
        grad[1] = 0.5 / sqrt(x[1]);
    }
    return sqrt(x[1]);
}

double myfunc2(int n, const double *x, double *grad, void *my_func_data)
{
	if (grad) {
        grad[0] = 0.0;
        grad[1] = 0.5 / sqrt(x[1]);
    }

	double dist = 0;
	for( int i=0; i<n; i++)
		dist += x[i]*x[i] + cos(x[i])*10;

	unsigned int seed;
    return dist + unifrnd<double>( &seed) * 0;
}

double myfuncnorm(int n, const double *x, double *grad, void *my_func_data)
{

	double dist = 0;
	for( int i=0; i<n; i++)
		dist += pow( x[i] - 2, 2);

	double rnd = 0;
	normrnd( &rnd, 1, (unsigned int *)my_func_data);

	double s2 = 10;
//	return sqrt(dist) + rnd*2*sqrt(s2);
	return dist + rnd*4*s2;
}

double myfuncnorm2(int n, const double *x, double *grad, void *my_func_data)
{

	double result = 1;

	double m=0, s=1;

	for( int i=0; i<n; i++)
		result *= exp( -0.5 * pow( (x[i] - m) / s, 2) ) / sqrt(2*M_PI) / s;

	return result; //+ RAND01*0.01;
}

double parabol(int n, const double *x, double *grad, void *my_func_data)
{

	double result = 0;

	for( int i=0; i<n; i++)
		result += x[i]*x[i];

	return result + unifrnd<double>();
}


double myfuncbinom(int n, const double *x, double *grad, void *my_func_data)
{
    int nsamples = ((int*)my_func_data)[0]+((int*)my_func_data)[1];
    int sample[2] = { 0, 0};
    unsigned int seed = time(0);

    // flip coins
    for( int i=0; i<nsamples; i++){
    	if( x[0] > unifrnd<double>( &seed)){
    		sample[0]++;
    	}else
    		sample[1]++;
    }

    return pow( ((int*)my_func_data)[0]-sample[0], 2) + pow( ((int*)my_func_data)[1]-sample[1], 2);
}





double myfunc3(int n, const double *x, double *grad, void *my_func_data)
{
	double parv1[] = {0, 1/24., 100,   50, 0.25,   0.0033, 0.0005, 0.003};
	double parv2[n+1];

	parv2[0] = 0; // random seed
	for( int i=0; i<n; i++){
		parv2[i+1] = pow( 10., x[i]);
	}
	//printVector( parv2, n+1, " %10.3e");

    return log10( compare( n+1, parv1, parv2, 5, 1));
}

double externalfuncPipe(int n, const double *x, double *grad, void *my_func_data){
	char command[1024];
	sprintf( command, "%s", (const char *)my_func_data); // executable
	for( int i=0; i<n; i++)
		sprintf( command, "%s %20e", command, x[i]);
	//fprintf( stderr, "COMMAND: %s\n", command);

	FILE *lsofFile_p = popen( command, "r");
	char buffer[1024];
	double rtn = atof( fgets(buffer, sizeof(buffer), lsofFile_p));
	pclose(lsofFile_p);

	return rtn;
	//return RAND01;
}

extern double EpsilonLimit;

double externalfuncFile(int n, const double *x, double *grad, void *my_func_data){
	char command[1024];
	char rnd_filename[1024];

	sprintf( rnd_filename, "tmp_%i_%i.out", rand(), omp_get_thread_num());

	sprintf( command, "%s -e%e ", (const char *)my_func_data, EpsilonLimit); // executable
	for( int i=0; i<n; i++)
		sprintf( command, "%s %20e", command, x[i]);
	sprintf( command, "%s > %s", command, rnd_filename);
	//fprintf(stderr, "%s\n", command);
	system( command);

	FILE *lsofFile_p = fopen( rnd_filename, "r");
	char buffer[1024];
	double rtn = atof( fgets(buffer, sizeof(buffer), lsofFile_p));
	fclose(lsofFile_p); remove(rnd_filename);

	return rtn;
	//return RAND01;
}
