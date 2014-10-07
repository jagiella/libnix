#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <ctype.h>

#define RAND01 ((double)rand()/(double)(RAND_MAX+1.))

int  main(int argc, char **argv)
{
	// INITIALIZE SETTINGS
	double (*func) (double) = sin;
	double randomContribution = 0;

	// ADAPT SETTINGS
	char c;
	bool stop = false;
	int last_optind = optind;
	while (!stop && (c = getopt (argc, argv, "m:r:0::1::2::3::4::5::6::7::8::9::")) != -1){
		switch (c){
		case 'm':
			if( strcmp( optarg, "cos") == 0)
				func = cos;
			break;
		case 'r':
			randomContribution = atof(optarg);
			break;
		default:
			//fprintf(stderr, "optopt='%i'\n", (int)optopt);
			if( isdigit( optopt) || optopt==0){
				//fprintf(stderr, "ISDIGIT: optind=%i\n", optind);
				optind = last_optind;
				stop = true;
			}
			break;
		}
		last_optind = optind;
	}
//fprintf(stderr, "optind=%i\n", optind);
	// READ INPUT VALUES
	int    n = argc - optind;
	double x[n];
	int    ix = 0;
	for (int i = optind; i < argc; i++){
		//fprintf( stderr, "[ %s ]", argv[i]);
		x[ix++]=atof( argv[i]);
	}

	// CALCULATE & PRINT OUTPUT
	double dist = 0;
	for( int i=0; i<n; i++)
		dist += x[i]*x[i] + (*func)(x[i])*10;
    printf( "%e", dist + RAND01*randomContribution);

    return 0;
}
