template <int s>
struct butchertableau{
	int n;
	double c[s];
	double a[s][s];
	double b[s];
	double e[s];
} ;

struct butchertableau<4> rungekutta4thorder = {
		4,
		{  0., 0.5,   0.5, 1.0},
		{{  0.,   0.,  0.,  0.},
		{1/2.,   0.,  0.,  0.},
		{  0., 1/2.,  0.,  0.},
		{  0.,   0., 1.0,  0.}},
		{1/6., 1/3., 1/3., 1/6.},
		{1/6., 2/3., 1/6., 0.0}
};

struct butchertableau<4> rungekutta38rule = {
	4,
	{  0., 1/3.,   2/3., 1.0},
	{
			  {  0.,   0.,  0.,  0.},
			  {1/3.,   0.,  0.,  0.},
			  {-1/3., 1.0,  0.,  0.},
			  {  1., -1., 1.0,  0.},
	},
	{1/8., 3/8., 3/8., 1/8.},
	{1/6., 2/3., 1/6., 0.0}
} ;

/*struct {
	int n=3;
	double c[3] = {0.0, 0.5, 1.0};
	double a[3][3] = {
			  {0.0, 0.0, 0.0},
			  {1/2., 0.0, 0.0},
			  {-1.,    2., 0.0}
	};
	double b[4] = {1/6., 2/3., 1/6.};
} rungekutta3rdorder;*/


struct butchertableau<6> rungekuttafehlberg = {
		6,
		{0, 1/4., 3/8., 12/13., 1., 1/2.},
		{{0,			0,				0,			0,			0,		0},
		 {1/4.,			0,				0,			0,			0,		0},
		 {3/32.,		9/32.,			0,			0,			0,		0},
		 {1932/2197.,	-7200/2197.,	7296/2197.,	0,			0,		0},
		 {439/216.,		-8.,			3680/513.,	-845/4104.,	0,		0},
		 {-8/27.,		2.,				-3544/2565.,1859/4104.,	-11/40.	,0}
		},
		{16/135.,	0.,	6656/12825.,	28561/56430.,	-9/50.,	2/55.},
		{25/216.,	0.,	1408/2565.,		2197/4104.,		-1/5.,	0.}
};


#define FIX_TIME_STEP 0

template <int order >
int solveRungeKutta4( gsl_odeiv2_system sys, double *t, double tmax, double *h, double *y, double error, struct butchertableau<order> tableau,
		double** &k, double* &dydt, double* &yi, int &length // temporary data
	)
{
	//double k   [sys.dimension][tableau.n];
	//double dydt[sys.dimension];
	//double yi  [sys.dimension];
/*	double **k   = (double **) malloc(sys.dimension*sizeof(double*));
	for( int i=0; i<sys.dimension; i++) k[i] = (double * ) malloc(tableau.n*sizeof(double));
	double *dydt = (double * ) malloc(sys.dimension*sizeof(double));
	double *yi   = (double * ) malloc(sys.dimension*sizeof(double));*/
	//double e   [sys.dimension];
	double ti, abs_e=0;

	double h_prop = *h;

	do{
		if( *t + *h > tmax)
			*h = tmax-*t;

		// calculate k's
		for( int i=0; i<tableau.n; i++){
			// intermediate time point
			ti = (*t) + tableau.c[i]* (*h);
	//		fprintf( stderr, "t[%i] = %e\n", i, ti);

			// intermediate y
	//		fprintf( stderr, "y[d] = [");
			for( unsigned int d=0; d<sys.dimension; d++){
				//if(*t >= 270.00)fprintf( stderr, "y[%i] = \n", d);
				yi[d] = y[d];
				for( int j=0; j<i; j++)
					yi[d] += tableau.a[i][j]*k[d][j];
	//			fprintf( stderr, "%5.2e ", yi[d]);
			}
	//		fprintf( stderr, "]\n");

			(*sys.function)( ti, yi, dydt, sys.params);
	//		fprintf( stderr, "k[d][%i] = [", i);
			for( unsigned int d=0; d<sys.dimension; d++){
				k[d][i] = (*h) * dydt[d];
	//			fprintf( stderr, "%5.2e ", k[d][i]);
			}
	//		fprintf( stderr, "]\n");
		}

		// update y's
		abs_e=0;
		for( unsigned int d=0; d<sys.dimension; d++){
			//e[d]=0;
			double e=0;
			for( int i=0; i<tableau.n; i++){
				//y[d] += tableau.b[i]*k[d][i];
				//e[d] += (tableau.b[i] - tableau.e[i])*k[d][i];
				e += (tableau.b[i] - tableau.e[i])*k[d][i];
			}
			//abs_e += fabs(e[d]);
			abs_e += fabs(e);
		}
		if( isnan( abs_e)){
			//fprintf( stderr, "%e %e %e %e\n", (*h), error, abs_e, 1./order);
			//exit(0);
		}

		if( !FIX_TIME_STEP){
			// PROPOSE ADAPTATION
			h_prop = (*h) * 0.9 *pow(error / abs_e, 1./order);
			//fprintf( stderr, "err: %e, h:%e\n", abs_e, (*h));

			if( abs_e > error){
				//fprintf( stderr, "err: %e, h:%e (-> h:%e, t:%e - %e)\n", abs_e, (*h), (*h) * 0.9 *pow(error / abs_e, 1./order), *t, tmax);
				(*h) = h_prop;
				//fprintf( stderr, "adapt h: %e\n", (*h));
			}

			if( (*h)<=1e-7){
				(*h) =1e-7;
				//fprintf( stderr, "adapt h: %e\n", (*h));
			}

			if( isinf(*h)){
				(*h) = tmax - *t;
				//fprintf( stderr, "adapt h: %e\n", (*h));
			}
		}




	}while( !FIX_TIME_STEP && abs_e > error);

	// update y's
	for( unsigned int d=0; d<sys.dimension; d++){
		for( int i=0; i<tableau.n; i++){
			y[d] += tableau.b[i]*k[d][i];// + normrnd( 0, 0.001)*sqrt(*h);
		}
	}

	// update t


	/*if(!FIX_TIME_STEP){
		(*h) = (*h) * 0.9 *pow(error / abs_e, 1./order);
		if( error == 0) (*h) = tmax - *t;
	}else*/
		(*t) += (*h);
		*h = h_prop;

		//exit(0);

	//fprintf( stderr, "FINISHED\n");
//		free(yi);
//		for( int i=0; i<sys.dimension; i++) free(k[i]);
//		free(k);
//		free(dydt);

	return 0;
}
