/*
 * compare.cpp
 *
 *  Created on: Oct 27, 2014
 *      Author: jagiella
 */

#include <stdlib.h>
#include <stdio.h>
#include <string.h> // strcpy

#include "Montecarlo.h"
#include "../matrix.hpp"
#include "../statistics.hpp"
#include <gsl/gsl_interp.h>
#include <gsl/gsl_spline.h>
// DATA: X-Y

typedef struct {
	int 	dim, size;
	double *x; // x
	double **y;// y's
	double *m; // mean of y
	double *s; // standard deviation of y
} comparison_t;

comparison_t create_comparison(){
	comparison_t c;
	c.size = c.dim = 0;
	c.x = c.m = c.s = 0;
	c.y = 0;
	return c;
}
enum { mean_vs_mean, mean_vs_single};

#define inrange( a, x, b) ( a<=x && x<=b )

void interpolate( double *xi, double *yi, int ni, double *xo, double *yo, int no){

	// init
	gsl_interp_accel *acc = gsl_interp_accel_alloc ();
	gsl_spline *spline = gsl_spline_alloc (gsl_interp_cspline, ni);
	gsl_spline_init (spline, xi, yi, ni);

	// interpolate
	for( int i=0; i<no; i++ ) // for all data points
	if( xo[i] >= xi[0] && xo[i] <= xi[ni-1] )
		// interpolate
		yo[i] = gsl_spline_eval (spline, xo[i], acc);
	else
		yo[i] = 1./0.;

	// free
	gsl_spline_free (spline);
	gsl_interp_accel_free (acc);

}

double compare( comparison_t d1, comparison_t d2, char mode ){
	double negloglik = 0;

	if( d1.m == 0){
		d1.m = (double*)malloc( d1.dim * sizeof(double));
		mean( d1.y, d1.size, d1.dim, d1.m, 2);
	}

	switch( mode){
	case mean_vs_mean:{

		if( d2.m == 0){
			d2.m = (double*)malloc( d2.dim * sizeof(double));
			mean( d2.y, d2.size, d2.dim, d2.m, 2);
		}

		double d2_m[ d1.dim];
		interpolate(
				d2.x, d2.m, d2.dim, // map d2.m onto d1.x  |
				d1.x, d2_m, d1.dim);//                     v

		for( int i=0; i<d1.dim; i++ ) // for all data points
			if( inrange( d2.x[0], d1.x[i], d2.x[d2.dim-1]) && d1.s[i]>0)
			negloglik += 0.5 * pow( (d1.m[i] - d2_m[i])/d1.s[i], 2);
	}break;

	case mean_vs_single:{
		double d2_y[ d1.dim];


		for( int j=0; j<d2.size; j++){
		interpolate(
				d2.x, d2.y[j], d2.dim, // map d2.y onto d1.x  |
				d1.x, d2_y,    d1.dim);//                     v


		for( int i=0; i<d1.dim; i++ ) // for all data points
			if( inrange( d2.x[0], d1.x[i], d2.x[d2.dim-1])  && d1.s[i]>0)
			negloglik += 0.5 * pow( (d1.m[i] - d2_y[i])/d1.s[i], 2);
		}
	}break;
	}

	return negloglik;
}


int readFile( const char* filename, double **&cols, int ncol)
{
	//printf( "[ %i ]\n", sizeof(cols));
	if( !cols){
		cols = (double**)realloc( cols, sizeof(double*)*ncol);
		for( int i=0; i<ncol; i++)
			cols[i] = 0;
	}
	int nrow = 0;

	FILE *fp = fopen( filename, "r");
	char buffer[1024], *ptr;

	while( fgets( buffer, 1024, fp ))if(buffer[0] != '#'){
		ptr = buffer;
		nrow++;

		for( int i=0; i<ncol; i++){
			cols[i] = (double*)realloc( cols[i], sizeof(double)*nrow);
			cols[i][nrow-1]  = strtof( ptr, &ptr) ;
		}
	}
	fclose(fp);

	return nrow;
}
int readFileColumn( const char* filename, double *&col, int icol)
{
	int nrow = 0;

	FILE *fp = fopen( filename, "r");
	char buffer[1024], *ptr;

	while( fgets( buffer, 1024, fp ))
	if(buffer[0] != '#'){
		ptr = buffer;
		nrow++;

		for( int i=0; i<icol; i++)
			strtof( ptr, &ptr) ;

		col = (double*)realloc( col, sizeof(double)*nrow);
		col[nrow-1]  = strtof( ptr, &ptr) ;
	}
	fclose(fp);

	return nrow;
}

int readFileColumn( char** filename, int n, double **&col, int icol)
{
	int nrow = 0;
	double **old=col;
	col = (double**)realloc( col, sizeof(double*)*n);
	for( int i=0; i<n; i++){
		if( old != col) col[i] = 0;
		nrow = readFileColumn( filename[i], col[i], icol);
	}

	return nrow;
}

void addOption( int &n, char **&options, const char *new_option)
{
	n++;
	options = (char**) realloc( options, sizeof(char*)*n);
	options[n-1] = (char*)malloc(sizeof(char) * 512);
	//strcpy( options[n-1], new_option);
	sprintf( options[n-1], "%s",new_option);
}
#include <unistd.h>
int main( int argc, char **argv)
{
	// READ DATA

	char data_filename[1024];
	comparison_t data_growthcurve = create_comparison();
	sprintf( data_filename, "/Users/jagiella/Dropbox/Work/PUBLICATIONS/Paper-wp2-lungsys/data/SK-MES-1/SK-MES_Radius_O0.28_G25.dat");
	data_growthcurve.dim = readFileColumn( data_filename, data_growthcurve.x, 0);
	                        readFileColumn( data_filename, data_growthcurve.m, 1);
	                        readFileColumn( data_filename, data_growthcurve.s, 2);
	comparison_t data_KI67 = create_comparison();
	sprintf( data_filename, "/Users/jagiella/Dropbox/Work/PUBLICATIONS/Paper-wp2-lungsys/data/SK-MES-1-CryoSections-MACBOOK/ConditionIII/T3/Ki67_MEDIAN_histDiv_StandardDerivation.dat");
	data_KI67.dim = readFileColumn( data_filename, data_KI67.x, 0);
					readFileColumn( data_filename, data_KI67.m, 3);
					readFileColumn( data_filename, data_KI67.s, 11);
	/*for( int i=0; i<data_KI67.dim; i++){
			fprintf( stderr, "%e %e\n", data_KI67.x[i], data_KI67.m[i]);
		}*/

	// SIMULATE

	bool writeRawDataToFile = false;
	int n = 1;
	int parc = 0;
	char **parv = 0;
	char simdir[512];
	pid_t pid = getpid();
	sprintf( simdir, "TEST%i", pid);
	char option[512];
	addOption( parc, parv, "nix-tumor3d");
	addOption( parc, parv, "-x1");
	addOption( parc, parv, "-y500");
	addOption( parc, parv, "-RNoRadialProfiles");
	addOption( parc, parv, "-RRadialProfilesTime17");
	addOption( parc, parv, "-RNoSliceOutput");
	addOption( parc, parv, "-k20");
	addOption( parc, parv, "-RExponentialReentranceProbability");
	addOption( parc, parv, "-M10");
	sprintf( option, "-d%s", simdir);
	addOption( parc, parv, option);

	for( int i=1; i<argc; i++){
		if( argv[i][0] == 'R'){
			writeRawDataToFile = true;
		}else{
			switch(i){
			case 1: sprintf( option, "-v%s", argv[i]); break;
			case 2: sprintf( option, "-RReentranceProbabilityLength%s", argv[i]); break;
			case 3: sprintf( option, "-RInitialRadius%s", argv[i]); break;
			case 4: sprintf( option, "-RInitialQuiescentFraction%s", argv[i]); break;
			case 5: sprintf( option, "-RECMProductionRate%s", argv[i]); break;
			case 6: sprintf( option, "-RECMDegradationRate%s", argv[i]); break;
			case 7: sprintf( option, "-RECMThresholdQuiescence%s", argv[i]); break;
			}
			addOption( parc, parv, option);
		}
	}

	//fprintf(stderr, " SIM\n");
	montecarlo( parc, parv);
	//fprintf(stderr, "EMD SIM\n");

	// read sim data
	char **sim_filename = allocMatrix<char>( n, 1024);
	for( int i=0; i<n; i++)
		sprintf( sim_filename[i], "%s/rel%i.data.dat", simdir, i);
	comparison_t sim_growthcurve = create_comparison();
	sim_growthcurve.size = n;
	sim_growthcurve.dim = readFileColumn( sim_filename[0], sim_growthcurve.x, 0);
	sim_growthcurve.dim = readFileColumn( sim_filename, n, sim_growthcurve.y, 2);
	for( int j=0; j<sim_growthcurve.dim; j++)
		sim_growthcurve.x[j] /= 24.;
	//printVector<double>( sim_growthcurve.x, sim_growthcurve.dim, "%.2e ");

	comparison_t sim_KI67 = create_comparison();
	for( int i=0; i<n; i++)
		sprintf( sim_filename[i], "%s/rel%i.radialProfiles_day16-17.pov.dat", simdir, i);
	sim_KI67.size = n;
	sim_KI67.dim = readFileColumn( sim_filename[0], sim_KI67.x, 0);
	sim_KI67.dim = readFileColumn( sim_filename, n, sim_KI67.y, 3);
	for(int i=0; i<sim_KI67.size; i++) for(int j=0; j<sim_KI67.dim; j++) if( isnan( sim_KI67.y[i][j]) ) sim_KI67.y[i][j] = 0;

	// COMPARE

	//double negLogLik = compare( data, sim, mean_vs_mean);
	double negLogLik
			= compare( data_growthcurve, sim_growthcurve, mean_vs_single) / data_growthcurve.dim
			+ compare( data_KI67, sim_KI67, mean_vs_single) / data_KI67.dim;

	fprintf( stdout, "%e\n", log10( negLogLik) );

	//
	if( writeRawDataToFile){
		FILE *fp;

		fp = fopen( "raw.data.dat", "w+");
		for( int i=0; i<data_growthcurve.dim; i++)
			fprintf( fp, "%e %e %e\n", data_growthcurve.x[i], data_growthcurve.m[i], data_growthcurve.s[i]);
		fclose( fp);
		fp = fopen( "raw.sim.dat", "w+");
		for( int i=0; i<sim_growthcurve.size; i++){
			for( int j=0; j<sim_growthcurve.dim; j++)
				fprintf( fp, "%e %e\n", sim_growthcurve.x[j], sim_growthcurve.y[i][j]);
			fprintf( fp, "\n");
		}
		fclose( fp);

		fp = fopen( "raw.KI67.data.dat", "w+");
		for( int i=0; i<data_KI67.dim; i++)
			fprintf( fp, "%e %e %e\n", data_KI67.x[i], data_KI67.m[i], data_KI67.s[i]);
		fclose( fp);
		fp = fopen( "raw.KI67.sim.dat", "w+");
		for( int i=0; i<sim_KI67.size; i++){
			for( int j=0; j<sim_KI67.dim; j++)
				fprintf( fp, "%e %e\n", sim_KI67.x[j], sim_KI67.y[i][j]);
			fprintf( fp, "\n");
		}
		fclose( fp);

	}

	char command[512];
	sprintf( command, "rm -rf %s", simdir);
	system( command);

}