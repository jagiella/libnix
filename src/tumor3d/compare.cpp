/*
 * compare.cpp
 *
 *  Created on: Oct 27, 2014
 *      Author: jagiella
 */

#include <stdlib.h>
#include <stdio.h>
#include <string.h> // strcpy
#include <sys/utsname.h> // uname
#include <time.h> // time

#include "AsciiIO.h"
#include "Montecarlo.h"

#include "../matrix.hpp"
// DATA: X-Y












void addOption( int &n, char **&options, const char *new_option)
{
	n++;
	options = (char**) realloc( options, sizeof(char*)*n);
	options[n-1] = (char*)malloc(sizeof(char) * 512);
	//strcpy( options[n-1], new_option);
	sprintf( options[n-1], "%s",new_option);
}
#include <unistd.h>
#include <float.h>
int main( int argc, char **argv)
{
	// READ OPTIONS
	/*double maxEpsilon = DBL_MAX;
	char c;
	while ((c = getopt (argc, argv, "e:")) != -1){
		switch (c){
			case 'e': {
				//maxEpsilon = atof( optarg);
			}break;

		}
	}
	fprintf(stderr, "(( %i ))\n", optind);*/

	// READ DATA
	FILE *fp;
	char data_filename_gc[1024];
	sprintf( data_filename_gc, "/Users/jagiella/Dropbox/Work/PUBLICATIONS/Paper-wp2-lungsys/data/SK-MES-1/SK-MES_Radius_O0.28_G25.dat");
/*	comparison_t data_growthcurve = create_comparison();
	data_growthcurve.dim = readFileColumn( data_filename_gc, data_growthcurve.x, 0);
	                       readFileColumn( data_filename_gc, data_growthcurve.m, 1);
	                       readFileColumn( data_filename_gc, data_growthcurve.s, 2);
*/
	char data_filename_ki67[1024];
	sprintf( data_filename_ki67, "/Users/jagiella/Dropbox/Work/PUBLICATIONS/Paper-wp2-lungsys/data/SK-MES-1-CryoSections-MACBOOK/ConditionIII/T3/Ki67_MEDIAN_histDiv_StandardDerivation.dat");
/*	comparison_t data_KI67 = create_comparison();
	data_KI67.dim = readFileColumn( data_filename_ki67, data_KI67.x, 0);
					readFileColumn( data_filename_ki67, data_KI67.m, 3);
					readFileColumn( data_filename_ki67, data_KI67.s, 11);
*/

	// SIMULATE

	bool writeRawDataToFile = false;
	int n = 1;
	int parc = 0;
	char **parv = 0;
	char simdir[512];
	pid_t pid = getpid();
	struct utsname unameData;
	uname(&unameData); // Might check return value here (non-0 = failure)
	//printf("%s", unameData.nodename);

	sprintf( simdir, "%s_PID%i", unameData.nodename, pid); //fprintf( stderr, "[OUTPUTDIR: %s]\n", simdir);
	bool remove_simulation_directory = true;

	char option[512];
	addOption( parc, parv, "nix-tumor3d");
	addOption( parc, parv, "-x1");
	addOption( parc, parv, "-y500");
	addOption( parc, parv, "-RNoRadialProfiles");
	addOption( parc, parv, "-RRadialProfilesTime17");
	addOption( parc, parv, "-RNoSliceOutput");
	addOption( parc, parv, "-k10");
	addOption( parc, parv, "-RExponentialReentranceProbability");
	addOption( parc, parv, "-M10");
	srand (time(NULL));
	sprintf( option, "-s%i", rand() % 1000);
	addOption( parc, parv, option);
	/*if (fp = fopen(data_filename_gc, "r")){
		fclose(fp);
		sprintf( option, "-1%s", data_filename_gc);
		addOption( parc, parv, option);
	}
	if (fp = fopen(data_filename_ki67, "r")){
		fclose(fp);
		sprintf( option, "-2%s", data_filename_ki67);
		addOption( parc, parv, option);
	}*/
	/*if( maxEpsilon != DBL_MAX){
	sprintf( option, "-3%e", maxEpsilon);
	addOption( parc, parv, option);
	}*/
	//sprintf( option, "-i /Users/jagiella/tmp/100pow%i", DIMENSIONS);
	//addOption( parc, parv, option);


	int par_count=0;
	for( int i=1; i<argc; i++){
		if( argv[i][0] == '-' && argv[i][1] == 'm'){
			sprintf( option, "-5%e", atof( &argv[i][2]));
			addOption( parc, parv, option);
		}else if( argv[i][0] == '-' && argv[i][1] == 'e'){
			sprintf( option, "-3%e", atof( &argv[i][2]));
			addOption( parc, parv, option);
		}else if( argv[i][0] == '-' && argv[i][1] == 'g'){
			sprintf( option, "-1%s", &argv[i][2]);
			addOption( parc, parv, option);
		}else if( argv[i][0] == '-' && argv[i][1] == 'k'){
			sprintf( option, "-2%s", &argv[i][2]);
			addOption( parc, parv, option);
		}else if( argv[i][0] == '-' && argv[i][1] == 'E'){
			sprintf( option, "-4%s", &argv[i][2]);
			addOption( parc, parv, option);
		}else if( argv[i][0] == '-' && argv[i][1] == 'd'){
			sprintf( simdir, "%s", &argv[i][2]);
			remove_simulation_directory = false;
			addOption( parc, parv, option);
		}else{
			switch(par_count){
			case 0: sprintf( option, "-v%s", argv[i]); break;
			case 1: sprintf( option, "-RReentranceProbabilityLength%s", argv[i]); break;
			case 2: sprintf( option, "-RInitialRadius%s", argv[i]); break;
			case 3: sprintf( option, "-RInitialQuiescentFraction%s", argv[i]); break;
			case 4: sprintf( option, "-RECMProductionRate%s", argv[i]); break;
			case 5: sprintf( option, "-RECMDegradationRate%s", argv[i]); break;
			case 6: sprintf( option, "-RECMThresholdQuiescence%s", argv[i]); break;
			}
			addOption( parc, parv, option);
			par_count++;
			//fprintf( stderr, "%i: %s\n", parc, parv[parc-1]);
		}
	}

	sprintf( option, "-d%s", simdir);
	addOption( parc, parv, option);


	//fprintf(stderr, " parc = %i\n", parc);

	//fprintf(stderr, " SIM\n");
	//parc--;
	double negLogLik = montecarlo( parc, parv);
	//fprintf(stderr, "Epsilon = %e\n", negLogLik);
	//fprintf(stderr, "EMD SIM\n");

	// read sim data
/*	char **sim_filename = allocMatrix<char>( n, 1024);
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
			= compare( data_growthcurve, sim_growthcurve, mean_vs_single) // data_growthcurve.dim
			+ compare( data_KI67, sim_KI67, mean_vs_single); // data_KI67.dim;
*/
	fprintf( stdout, "%e\n", ( negLogLik) );

	//
	/*if( writeRawDataToFile){
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

	}*/

	if( remove_simulation_directory){
		char command[512];
		sprintf( command, "rm -rf %s", simdir);
		system( command);
	}

}
