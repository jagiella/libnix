/*
 * bin_abc.cpp
 *
 *  Created on: Oct 29, 2014
 *      Author: jagiella
 */

#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <omp.h>
#include <unistd.h>
#include <time.h>
#include <float.h>


#include "../model.hpp"
#include "../statistics.hpp"
#include "../optimization.hpp"
#include "../testfunctions.hpp"



void usage(){

	fprintf( stderr, "usage: nix-abc [OPTION(S)]\nOPTIONS:\n");
	fprintf( stderr, " -m METHODE  is the function to optimize\n");
	fprintf( stderr, "             either builtin (tumor, norm, binom)\n");
	fprintf( stderr, "             or custom command\n");
	fprintf( stderr, " -l LIST     lower bound(s) as comma divided list, e.g. -l 1,-2,1\n");
	fprintf( stderr, " -u LIST     upper bound(s) as comma divided list, e.g. -u 3,0,10\n");
	fprintf( stderr, " -n INT      population size\n");
	fprintf( stderr, " -o FILENAME filename for output\n");
	fprintf( stderr, " -L          perform optimization in logarithmic scaling\n");
	fprintf( stderr, " -G          use speed-up by Gauss process approximation\n");
}

int main( int argc, char **argv){

	if(argc<=1){
		usage();
		return 0;
	}

	srand(time(0));
	// ABC parameters
	int d = 2;
	//int n = 100;
	int nparallel = 10;
	//int nthreads  = 5;
	//int iterations = 100;
	//int evaluations= 1000000;
	double (*func) (int, const double*, double*, void*) = parabol;
	void *func_data = 0;
	unsigned int func_seed = 0;
	int binom_data[2] = { 24, 16};
	//	int binom_data[2] = { 6, 4};

	double *lb=0;
	double *ub=0;

	int iCPU = 10;


//	char suffix[1024];
//	char filename[1024];
//	sprintf( suffix, "%s", "output" );

	optimoptions options = getoptions();
	options.MaxFunEvalsAvg = 4;
	options.MaxFunEvals    = 100000;
	options.MaxIter		   = 100;

	bool GP_Approximation = false;


	char c;
	while ((c = getopt (argc, argv, "m:l:u:n:o:a:LGIv")) != -1)
	switch (c){
	      case 'm': { // function to minimize
	  		if( strcmp( optarg, "tumor") == 0){
	  			func = myfunc3;
	  		}else if( strcmp( optarg, "norm") == 0){
	  			func = myfuncnorm;
	  			func_data = &func_seed;
	  			d=2;
	  		}else if( strcmp( optarg, "norm2") == 0){
	  			func = myfuncnorm2;
	  			func_data = 0;
	  			d=2;
	  		}else if( strcmp( optarg, "binom") == 0){
	  			func = myfuncbinom;
	  			func_data = &binom_data;
	  			d=1;
	  		}else{
	  			fprintf( stderr, "Use user-provided command [ %s ]\n", optarg);
	  			func_data = optarg;
	  #if __linux
	  			func = externalfuncPipe;
	  #else
	  			func = externalfuncFile;
	  #endif
	  		}
	      }break;

	      case 'l':{ // lower bound
	    	  	d=0;
				char *str = strtok( optarg,",");
				while (str != NULL){
					lb = (double*) realloc( lb, sizeof(double)*(d+1));
					lb[d++] = atof( str);
					str = strtok (NULL, ",");
				}
				fprintf( stderr, "Detected dimensionality to be d=%i\n", d);
	      }break;

	      case 'u':{ // upper bound
				d=0;
				char *str = strtok( optarg,",");
				while (str != NULL){
					ub = (double*) realloc( ub, sizeof(double)*(d+1));
					ub[d++] = atof( str);
					str = strtok (NULL, ",");
				}
	      }break;

	      case 'n': options.PopulationSize = atoi( optarg); break;
	      case 'a': options.MaxFunEvalsAvg = atoi( optarg); break;
	      case 'o': sprintf( options.OutputFile, "%s", optarg ); break;
	      case 'L': options.ParameterScaling = Logarithmic; break;
	      case 'G': GP_Approximation = true; break;
	      case 'I': options.AllowInterception = true; break;
	      case 'v': options.Display = true; break;
	}


	if( lb==0){
		lb = (double*) malloc( sizeof(double)*d);
		for(int i=0; i<d; i++) lb[i] = -FLT_MAX;
		fprintf( stderr, "Assume dimensionality to be d=%i\n", d);
	}


	if( ub == 0){
		ub = (double*) malloc( sizeof(double)*d);
		for(int i=0; i<d; i++) ub[i] = FLT_MAX;
	}


	double **sol = 0;
	int n=0;

	if( GP_Approximation)
		abcgp(
				// input
				0, lb, ub, d, func, func_data,
				// options
				&options,
				// output
				sol, n);
	else
		abc(
				// input
				0, lb, ub, d, func, func_data,
				// options
				&options,
				// output
				sol, n);

	return 0;
}
