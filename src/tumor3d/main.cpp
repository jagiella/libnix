#include "Montecarlo.h"

#include <stdio.h>

#define __MAIN__

int main( int argc, char **argv)
	{

	double epsilon = montecarlo( argc, argv);
	fprintf(stdout, "%e\n", epsilon);

	return 0;
	}


