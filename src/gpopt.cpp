/*
 * mymincon.cpp
 *
 *  Created on: Sep 2, 2014
 *      Author: jagiella
 */


#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "myopt.hpp"







double parabol(int n, const double *x, double *grad, void *my_func_data)
{
	if (grad)
	for( int i=0; i<n; i++)
        grad[i] = 2*x[i]
              //- sin(x[i])*10
                ;

	double dist = 0;
	for( int i=0; i<n; i++)
		dist += x[i]*x[i]
		   // + cos(x[i])*10
			  ;

    return dist + 10.*((double)rand()/(double)RAND_MAX - 0.5);
}

#include "model.hpp"
double myfunc3(int n, const double *x, double *grad, void *my_func_data)
{
	double parv1[] = {0, 1/24., 100,   50, 0.25,   0.0033, 0.0005, 0.003};
	double parv2[n+1];

	parv2[0] = 0; // random seed
	for( int i=0; i<n; i++){
		parv2[i+1] = pow( 10., x[i]);
	}
	//printVector( parv2, n+1, " %10.3e");

    return log10( compare( n+1, parv1, parv2, 50, 1) );
}
double tumorModel2x3(int n, const double *x, double *grad, void *my_func_data)
{
	//if(n==3) n=4;
	int nex = (n==2 ? n+2 : n+1);
	double parv1[] = {0, 1/24., 100,   50, 0.25,   0.0033, 0.0005, 0.003};
	double parv2[nex];

	parv2[0] = 0; // random seed
	for( int i=0; i<n; i++){
		parv2[i+1] = pow( 10., x[i]);
	}

	if( n>2){
		double tmp = parv2[3];
		parv2[3] = parv2[2];
		parv2[2] = tmp;
	}
	if( n==2){
		parv2[3] = parv2[2];
		parv2[2] = parv1[2];
	}
	//printVector( parv2, nex, "%5.2e | ");

    return log10( compare( nex, parv1, parv2, 50, 1) );
}




int main( int argc, char **argv){

	int dim = 10;
	double x0[dim], lb[dim], ub[dim], sol[dim], sol_true[dim];
	optimoptions options = getoptions();

	options.GradObj = true;
	options.MaxIter = 100;
	options.MaxFunEvalsAvg = 10;
	double (*func) (int, const double*, double*, void*);

	switch( 2){
	case 1:
		dim = 3;
		for( int i=0; i<dim; i++){
			x0[i] = -3.12312414;
			lb[i] = -100;
			ub[i] = 100;
			sol_true[i] = 0;
		}
		func = parabol;



		//gpopt( 0, lb, ub, dim, parabol, &options, sol);
		break;
	case 2:{
		dim = 3;
		double tlb[] = {-3,1,0,-3,-5,-5,-5};
		double tub[] = {-1,3,3,0,0,0,0};
		double tx0[] = {-2,1.5, log10(100), log10(0.75) };
		double tsol_true[] = {log10(1/24.),2, log10(50), log10(0.75) };
		for( int i=0; i<dim; i++){
			lb[i]=tlb[i];ub[i]=tub[i];x0[i]=tx0[i];sol_true[i]=tsol_true[i];
		}
		func = myfunc3;
	}break;
	case 3:{
		dim = 2;
		double tlb[] = {-3,	0,	1,	-3,-5,-5,-5};
		double tub[] = {-1,	3,	3,	0,0,0,0};
		double tx0[] = {-2,	1.5, log10(100), log10(0.75) };
		double tsol_true[] = {log10(1/24.), log10(50), 2, log10(0.75) };
		for( int i=0; i<dim; i++){
			lb[i]=tlb[i];ub[i]=tub[i];x0[i]=tx0[i];sol_true[i]=tsol_true[i];
		}
		func = tumorModel2x3;
	}
	}

	// test sampling
/*	int n=10000;
	int nhist = 1000;

	double minhist=0, maxhist=1e1, dhist=(maxhist-minhist)/(nhist-1);
	double hist[nhist];
	for( int i=0; i<nhist; i++)
		hist[i] = 0;
	for( int i=0; i<n; i++){
		double value = (*func)( dim, x0, 0, 0);
		if( minhist <= value && value <= maxhist){
			int idx = (int)((value-minhist)/dhist);
			hist[idx]++;
		}else{
			hist[nhist-1]++;
		}
	}
	for( int i=0; i<nhist; i++)
		fprintf(stdout, "%e %e\n", (i+0.5)*dhist+minhist, hist[i]);
	exit(0);
*/

	// CREATE GNUPLOT SCRIPTS
	FILE *fp = fopen( "gnuplot.glt", "w+");
	fprintf( fp, "set multiplot; set size %f,%f\n", 1./dim, 1./dim);
	for( int i=0; i<dim; i++){		double offset_i = (double)(i)/dim;
		for( int j=0; j<dim; j++){	double offset_j = (double)(dim-1-j)/dim;
			// DIAGONAL
			if( i==j){
				//fprintf( fp, "set logscale y; \n");
				fprintf( fp, "set origin %f,%f; plot [][:]\"linreg%i.dat\"u 1:2, \"linreg%i.dat\"u 1:3w l ls -1; unset logscale\n", offset_i, offset_j, i, i);
				//fprintf( fp, "unset logscale\n");
			}if( i<j)
				fprintf( fp, "set origin %f,%f; plot \"population.out\"u %i:%iw lp, '+' using (%e):(%e) lc 3 pt 6 ps 5 t'', \"grad.dat\"u %i:%i:(-$%i):(-$%i) notitle with vectors arrowstyle 1\n",
						offset_i, offset_j, i+2, j+2, sol_true[i], sol_true[j],
						i*2+1, j*2+1, i*2+2, j*2+2);
			if( i>j)
				fprintf( fp, "set origin %f,%f; splot \"gpline.dat\"u %i:%i:%iw p\n", offset_i, offset_j, i+1, j+1, dim+2, dim+3);

		}
	}
	fprintf( fp, "unset multiplot\n");
	fclose(fp);

	gplinesearch( x0, lb, ub, dim, func, &options, sol);


	for( int i=0; i<dim; i++)
		fprintf( stderr, "%.3e ", sol[i]);
	fprintf( stderr, " <- solution \n");


	return 0;
}
