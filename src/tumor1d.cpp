#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>
#include <math.h>

#include "model.hpp"

int main( int argc, char **argv)
{
/*	fprintf( stderr, "n1=%s, n2=%s, bias=%s\n", argv[1], argv[2], argv[3]);
	// TEST
	int n1=atoi(argv[1]), n2=atoi(argv[2]);
	double x1[n1], x2[n2];
	for( int i=0; i<n1; i++)
		x1[i] = rand()/(double)RAND_MAX;
	for( int i=0; i<n2; i++)
		x2[i] = rand()/(double)RAND_MAX + atof(argv[3]);

	KSTestResult rslt = KolmogorovSmirnoffTest(  n1, x1, n2, x2);
	fprintf( stderr, "H=%i, p=%e, D=%e => %s\n", rslt.H, rslt.pValue, rslt.KSstatistic, (rslt.H ? "different" :  "similar") );
*/
	// TEST END

	// OPTION PART
	int n=1;
	char c;
	while ((c = getopt (argc, argv, "n:")) != -1)
	    switch (c){
	    case 'n':
	    	n = atoi(optarg);
	    	break;
	    }

	// PARAMETER PART

	int    parn = argc-optind;
	double parv[parn];
	//parv[0] = 0; // random seed
	for( int i=optind; i<argc; i++){
		parv[i-optind] = atof(argv[i]);
		fprintf( stderr, "%e, ", parv[i-optind]);
	}
	fprintf( stderr, "optind=%i\n", optind);
	double **mout = (double **) malloc(n*sizeof(double*));
	for( int i=0; i<n; i++){
		parv[0]++;
		mout[i] = model( parn, parv);
	}

	FILE *fp = fopen( "tmp.dat", "w+");
	for( int j=0; j<600+300+300+300+300+300+300; j++){
		double mean = 0;
		double mean2= 0;
		for( int i=0; i<n; i++){
			mean += mout[i][j];
			mean2+= mout[i][j]*mout[i][j];
		}
		fprintf( fp, "%f %f\n ", mean/(double)n, sqrt( mean2/(double)n - pow(mean/(double)n,2)));
	}
	fclose( fp);

	for( int i=0; i<n; i++)
		free( mout[i]);
	free( mout);

	// GNUPLOT
	fp = fopen( "tumor1d.glt", "w+");
	char plot[] = "\"tmp.dat\"u 0:($1-$2):($1+$2) w filledcu, \"tmp.dat\"u 0:1w l lc -1";
	fprintf( fp, "set multiplot;\n");
	fprintf( fp, "set size %f,%f\n", 1/3., 1/2.);
	fprintf( fp, "set origin %f,%f; set xlabel ''; set ylabel 'growth curve'; plot [   0: 600] %s\n", 0/3., 1/2., plot);
	fprintf( fp, "set size %f,%f\n", 1/3., 1/3.);
	fprintf( fp, "set origin %f,%f; set xlabel ''; set ylabel 'KI67';  plot [ 601: 620] %s\n", 1/3., 2/3., plot);
	fprintf( fp, "set origin %f,%f; set xlabel ''; set ylabel 'KI67';  plot [ 901: 920] %s\n", 2/3., 2/3., plot);
	fprintf( fp, "set origin %f,%f; set xlabel ''; set ylabel 'TUNEL'; plot [1201:1220] %s\n", 1/3., 1/3., plot);
	fprintf( fp, "set origin %f,%f; set xlabel ''; set ylabel 'TUNEL'; plot [1501:1520] %s\n", 2/3., 1/3., plot);
	fprintf( fp, "set origin %f,%f; set xlabel ''; set ylabel 'ECM';   plot [1801:1820] %s\n", 1/3., 0/3., plot);
	fprintf( fp, "set origin %f,%f; set xlabel ''; set ylabel 'ECM';   plot [2101:2120] %s\n", 2/3., 0/3., plot);
/*	fprintf( fp, "set origin %f,%f; set xlabel ''; set ylabel 'KI67';  plot [ 601: 900] %s\n", 1/3., 2/3., plot);
	fprintf( fp, "set origin %f,%f; set xlabel ''; set ylabel 'KI67';  plot [ 901:1200] %s\n", 2/3., 2/3., plot);
	fprintf( fp, "set origin %f,%f; set xlabel ''; set ylabel 'TUNEL'; plot [1201:1500] %s\n", 1/3., 1/3., plot);
	fprintf( fp, "set origin %f,%f; set xlabel ''; set ylabel 'TUNEL'; plot [1501:1800] %s\n", 2/3., 1/3., plot);
	fprintf( fp, "set origin %f,%f; set xlabel ''; set ylabel 'ECM';   plot [1801:2100] %s\n", 1/3., 0/3., plot);
	fprintf( fp, "set origin %f,%f; set xlabel ''; set ylabel 'ECM';   plot [2101:2400] %s\n", 2/3., 0/3., plot);*/
	fprintf( fp, "unset multiplot\n");
	fclose(fp);

}
