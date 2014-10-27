#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <float.h>
#include <unistd.h>
#include <sys/stat.h>

//#if __linux
#define USE_OMP
//#endif

#ifdef USE_OMP
#include <omp.h>
#endif

#include <gsl/gsl_matrix.h>
#include <gsl/gsl_permutation.h>
#include <gsl/gsl_linalg.h>

#include "matrix.hpp"
#include "statistics.hpp"


#define RAND01         ((double)rand()      /((double)RAND_MAX+1.))
#define RAND01_R(seed) ((double)rand_r(seed)/(double)(RAND_MAX))
#define DET2( A) (A[0][0]*A[1][1] - A[0][1]*A[1][0])

double bla( double *a, double **B, double *c, int n){

	double rtn=0, tmpi;
	for( int i=0; i<n; i++){
		tmpi = 0;
		for( int j=0; j<n; j++)
			tmpi += a[j]*B[j][i];
		rtn += tmpi*c[i];
	}

	return rtn;
}

void minimum( double *x, int n, double &min_x, int &min_i)
{
	min_x = x[0];
	min_i = 0;

	for( int i=1; i<n; i++){
		if( min_x > x[i]){
			min_x = x[i];
			min_i = i;
		}
	}
}

#include <gsl/gsl_sort.h>
#include <gsl/gsl_statistics_double.h>
double quantile( double *x, int n, double f)
{
	double tmp[n];
	for( int i=0; i<n; i++) tmp[i] = x[i];
	gsl_sort( tmp, 1, n);
	return gsl_stats_quantile_from_sorted_data ( tmp, 1, n, f);
}

void _mean( double **A, int n, int m, double *mean){
	// mean
	for( int j=0; j<m; j++)
		mean[j] = A[0][j];
	for( int i=1; i<n; i++)
	for( int j=0; j<m; j++)
		mean[j]+= A[i][j];
	for( int j=0; j<m; j++)
		mean[j]/= (double)n;
}

void _cov( double **A, int n, int m, double *mean, double **cov){
	// covariance
	for( int j=0;  j<m;  j++)
		for( int jj=0; jj<m; jj++){
			cov[j][jj] = 0;
			for( int i=0;  i<n;  i++)
				cov[j][jj] += (A[i][j] - mean[j])*(A[i][jj] - mean[jj]);
			cov[j][jj] /= (double)(n-1);
		}
}

/*double **allocMatrix<double>( int n, int m){
	double **x_old = (double**)malloc( n*sizeof(double**));
	for( int i=0; i<n; i++)
		x_old[i] = (double*)malloc( m*sizeof(double*));
	return x_old;
}
void printVector( double *x, int m, const char *fmt = "%5.2e "){

	for( int j=0; j<m; j++)
		fprintf( stderr, fmt, x[j]);
	fprintf( stderr, "\n");
}*/
void normalizeVector( double *x, int m){

	double sum = 0;
	for( int j=0; j<m; j++)
		sum += x[j];
	for( int j=0; j<m; j++)
		x[j] /= sum;
}
/*void printMatrix( double **A, int n, int m, const char *fmt){

	for( int i=0; i<n; i++)
		printVector( A[i], m, fmt);
}*/
int cumsum( double *x, int n, double partial_sum)
{
	double sum = 0;

	for( int i=0; i<n; i++){
		sum += x[i];
		if( sum >= partial_sum){
			return i;
		}
	}

	fprintf( stderr, "ERROR: cumsum (%e < %e)\n", sum, partial_sum);
	printVector( x, n, "%e ");
	return -1;
}

void matrixMultiplication( double **A, double **B, int n, int m, int l, double **C){
	for( int i=0; i<n; i++)
		for( int j=0; j<l; j++){
			C[i][j] = 0;
			for( int k=0; k<m; k++)
				C[i][j] += A[i][k]*B[k][j];
		}
}

double mvnpdf( const double *x, double *mean, double **inv_s2, double det_s2, int n)
{
	double x0[n];
	for( int i=0; i<n; i++)
		x0[i] = x[i] - mean[i];
	return 1. / sqrt( pow(2*M_PI, n) * det_s2) * exp(-0.5*bla( x0, inv_s2, x0,n));
}

double kdepdf( const double *x, double **samples, double **inv_s2, double det_s2, int n, int nsamples)
{
	double pdf = 0;
	for( int i=0; i<nsamples; i++)
		pdf += mvnpdf( x, samples[i], inv_s2, det_s2, n);
	return pdf / (double)nsamples;
}

void mvnrnd( double *x, double *mean, double **s2, int n, unsigned int *p_seed)
{
	// sample from MVnorm dist around chosen point
	double z[n];
	normrnd( z, n, p_seed);

	//fprintf( stderr, "Sample: ");
	for( int i=0; i<n; i++){
		x[i] = mean[i];
		double tmp = 0;
		for( int j=0; j<n; j++){
			x[i] += s2[i][j] * z[j];
			tmp += s2[i][j] * z[j];
		}
		//fprintf( stderr, "%e (%e) ", tmp, x[i]);
	}
	//fprintf( stderr, "\n");


}

bool inbound( double *x, double *lb, double *ub, int d)
{
	for( int i=0; i<d; i++)
		if( x[i] < lb[i] || x[i] > ub[i])
			return false;
	return true;
}

template <class T>
void decompCholesky( T** &A, T** &L, int n) {

	gsl_matrix *gsl_L = gsl_matrix_alloc( n, n);
	for( int i=0; i<n; i++)
		for( int j=0; j<n; j++){
			gsl_L->data[i * gsl_L->tda + j] = A[i][j];
		}

	gsl_linalg_cholesky_decomp ( gsl_L);
	//gsl_linalg_cholesky_invert ( gsl_L);

	for( int i=0; i<n; i++){
		for( int j=0; j<=i; j++){
			L[i][j] = gsl_L->data[i * gsl_L->tda + j];
		}
		for( int j=i+1; j<n; j++)
			L[i][j] = 0;
	}

	gsl_matrix_free( gsl_L);
}

template <class T>
void decompCholeskyOwn( T** &A, T** &L, int n) {

	for( int i=0; i<n; i++)
		for( int j=0; j<=i; j++){
			double Summe = A[i][j];
			for( int k=0; k<=j-1; k++)
				Summe = Summe - L[i][k] * L[j][k];
			if( i > j){
				L[i][j] = Summe / L[j][j];   // Untere Dreiecksmatrix
				L[j][i] = 0.;                // Obere Dreiecksmatrix
				//L[j][i] = Summe / A[j][j]; // Obere Dreiecksmatrix
			}
			else if( Summe > 0)          // Diagonalelement
				L[i][i] = sqrt( Summe);       // ... ist immer groesser Null
			else
	            fprintf( stderr, "Matrix not positive definite\n");   // ERROR
		}
}

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

    return dist + RAND01*0;
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

double myfuncbinom(int n, const double *x, double *grad, void *my_func_data)
{
    int nsamples = ((int*)my_func_data)[0]+((int*)my_func_data)[1];
    int sample[2] = { 0, 0};

    // flip coins
    for( int i=0; i<nsamples; i++){
    	if( x[0] > RAND01){
    		sample[0]++;
    	}else
    		sample[1]++;
    }

    return pow( ((int*)my_func_data)[0]-sample[0], 2) + pow( ((int*)my_func_data)[1]-sample[1], 2);
}

typedef struct {
	double **x_old;
	double **invcov;
	double   detcov;
	int      n_old;
} kdedata;
double myfunckde(int n, const double *x, double *grad, void *my_func_data)
{
    return kdepdf( x, ((kdedata*)my_func_data)->x_old, ((kdedata*)my_func_data)->invcov, ((kdedata*)my_func_data)->detcov, n, ((kdedata*)my_func_data)->n_old);
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

    return log10( compare( n+1, parv1, parv2, 5, 1) );
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

double externalfuncFile(int n, const double *x, double *grad, void *my_func_data){
	char command[1024];
	char rnd_filename[1024];

	sprintf( rnd_filename, "tmp_%i_%i.out", rand(), omp_get_thread_num());

	sprintf( command, "%s", (const char *)my_func_data); // executable
	for( int i=0; i<n; i++)
		sprintf( command, "%s %20e", command, x[i]);
	sprintf( command, "%s > %s", command, rnd_filename);

	system( command);

	FILE *lsofFile_p = fopen( rnd_filename, "r");
	char buffer[1024];
	double rtn = atof( fgets(buffer, sizeof(buffer), lsofFile_p));
	fclose(lsofFile_p); remove(rnd_filename);

	return rtn;
	//return RAND01;
}

//void read2array( char *str, const char delimiter, ){

//}

double rekursiveIntegration( double *startingPoint, int startingDimension, int excludeDimension,
		double *lowerBounds, double *upperBounds, double *stepSize, int dimensions,
		double (*func) (int, const double*, double*, void*), void *my_func_data){

	if( startingDimension == excludeDimension )
		startingDimension ++;


	if( startingDimension >= dimensions){
		// Evaluate
		return (*func) ( dimensions, startingPoint, 0, my_func_data);
	}else{
		double newPoint[dimensions];
		double newDimension = startingDimension + ( startingDimension + 1 != excludeDimension ? 1 : 2);
		double result = 0;
		for( int d=0; d<dimensions; d++)
			newPoint[d] = startingPoint[d];
		for( newPoint[startingDimension] = lowerBounds[startingDimension]; newPoint[startingDimension] <= upperBounds[startingDimension]; newPoint[startingDimension] += stepSize[startingDimension]){
			result += rekursiveIntegration( newPoint, newDimension, excludeDimension,
					lowerBounds, upperBounds, stepSize, dimensions,
					func, my_func_data);
		}
		return result*stepSize[startingDimension];
	}
}

void summary_statistics(
		double (*func) (int, const double*, double*, void*), void *func_data, int d, const double *parameters, int n,
		double &y, double &s2)
{
	y  = 0;
	s2 = 0;
	for( int i=0; i<n; i++){
		double val = (*func) (d, parameters, 0, func_data);
		y  += val;
		s2 += val*val;
	}
	y  /= (double)n;
	s2  = s2/(double)n - pow(y, 2); // s2
	s2 /= (double)n; // SEM^2

}

#include <stdlib.h>
int cmp( const void * a, const void * b){
	if( *(double*)a < *(double*)b ) return -1;
	if( *(double*)a ==*(double*)b ) return  0;
	if( *(double*)a > *(double*)b ) return +1;
}

void getKDESample(
	double **x, 	double *f, 		int n,		// in
	double **x_kde, double *f_kde, 	int n_kde	// out
){
	double f_crit = quantile( f, n, (double)n_kde/(double)n);

	n_kde = 0;
	for( int i=0; i<n; i++){
		if( f[i] <= f_crit){
			f_kde[n_kde] = f[i];
			x_kde[n_kde] = x[i];
			n_kde++;
		}
	}
}

void transform( double **xin, double **xout, int n, int d, double *m, double **inv_sqrt_S){
	for( int i=0; i<n; i++){
		for( int di=0; di<d; di++){
			xout[i][di] = 0;
			for( int dj=0; dj<d; dj++)
				xout[i][di] += inv_sqrt_S[di][dj]*(xin[i][dj]-m[dj]);
		}
	}
}

#include "gp.hpp"
#include "cov.hpp"
#include "optimization.hpp"
int main( int argc, char **argv){
	srand(time(NULL));
	//srand( 0);

	bool TRANSFORM = true;

/*	double startingPoint[] = { 0, 0, 0};
	double lowerBounds[] = { -10, -10, -10};
	double upperBounds[] = {  10,  10,  10};
	double stepSize[] = { 0.3, 0.3, 0.3};
	int startingDimension = 0;
	for( startingPoint[startingDimension]=-1; startingPoint[startingDimension]<=1; startingPoint[startingDimension]+=0.1)
		fprintf(stderr, "%e\n", rekursiveIntegration( startingPoint, 0, startingDimension,
				lowerBounds, upperBounds, stepSize, 2,
				myfuncnorm2, 0));
	fprintf(stderr, "%e\n", rekursiveIntegration( startingPoint, 0, -1,
			lowerBounds, upperBounds, stepSize, 3,
			myfuncnorm2, 0));
	exit(0);
*/
	// ABC parameters
	int d = 2;
	int n = 100;
	int nparallel = 10;
	//int nthreads  = 5;
	int iterations = 100;
	int evaluations= 1000000;
	double (*func) (int, const double*, double*, void*) = myfunc2;
	void *func_data = 0;
	unsigned int func_seed = 0;
	int binom_data[2] = { 24, 16};
	//	int binom_data[2] = { 6, 4};


	// USER-PROVIDED PARAMETERS
	if(argc>=2){ // function to minimize
		if( strcmp( argv[1], "tumor") == 0){
			func = myfunc3;
		}else if( strcmp( argv[1], "norm") == 0){
			func = myfuncnorm;
			func_data = &func_seed;
			d=2;
		}else if( strcmp( argv[1], "norm2") == 0){
			func = myfuncnorm2;
			func_data = 0;
			d=2;
		}else if( strcmp( argv[1], "binom") == 0){
			func = myfuncbinom;
			func_data = &binom_data;
			d=1;
		}else{
			fprintf( stderr, "Use user-provided command [ %s ]\n", argv[1]);
			func_data = argv[1];
#if __linux
			func = externalfuncPipe;
#else
			func = externalfuncFile;
#endif
		}
	}

	double *lb=0;
	if(argc>2){ // lower bound
		d=0;
		char *str = strtok( argv[2],",");
		while (str != NULL){
		    lb = (double*) realloc( lb, sizeof(double)*(d+1));
		    lb[d++] = atof( str);
		    str = strtok (NULL, ",");
		}
		fprintf( stderr, "Detected dimensionality to be d=%i\n", d);
	}else{
		lb = (double*) malloc( sizeof(double)*d);
		for(int i=0; i<d; i++) lb[i] = -FLT_MAX;
		fprintf( stderr, "Assume dimensionality to be d=%i\n", d);
	}

	double *ub=0;
	if(argc>3){ // upper bound
		d=0;
		char *str = strtok( argv[3],",");
		while (str != NULL){
		    ub = (double*) realloc( ub, sizeof(double)*(d+1));
		    ub[d++] = atof( str);
		    str = strtok (NULL, ",");
		}
	}else{
		ub = (double*) malloc( sizeof(double)*d);
		for(int i=0; i<d; i++) ub[i] = FLT_MAX;
	}


	int iCPU = 10;
	int n_accepted = ( argc>4 ? atoi( argv[4]) : 100 );

	char suffix[1024];
	char filename[1024];
	sprintf( suffix, "%s", ( argc>5 ? argv[5] : "output" ));


	// init variables
	fprintf( stderr, "INIT\n");
	//double epsilon = DBL_MAX;
	//double **x_new, **x_old, **x_tmp, **x_par;
	//double  *f_new,  *f_old,  *f_tmp,  *f_par;
	//double  *w_new,  *w_old,  *w_tmp,  *w_par;
	//int      n_new,   n_old,            n_par = nparallel;
	double mean[d], **cov = allocMatrix<double>(d,d), **invcov = allocMatrix<double>(d,d), **invcovL = allocMatrix<double>(d,d), **covL = allocMatrix<double>(d,d), detcov;

	int n_max = 10000;
	double **x, *f, *s2, *w;
	x = allocMatrix<double>( n_max, d);
	double **x_tf = allocMatrix<double>( n_max, d);
	f = (double*) malloc(n_max*sizeof(double));
	s2= (double*) malloc(n_max*sizeof(double));
	w = (double*) malloc(n_max*sizeof(double));
	n = 0;
	int n_avg = 4;


#ifdef USE_OMP
	//int iCPU = nparallel;//= omp_get_num_procs();
	fprintf( stderr, "SET NUM THREADS: %i\n", iCPU);
	omp_set_num_threads( iCPU);
#endif
	//printMatrix( x_old, n_old, d, "%10.3e ");


	// stats
	clock_t t, t_sampling = 0, t_evaluation = 0;


	// run iteration
	int count_evals_per_iteration[iterations];
	int count_evals = 0;
	unsigned int sampling_seed = 0;

	// INITIAL SAMPLING
	sprintf( filename, "%s.population.dat", suffix);
	FILE *fp_population = fopen( filename, "w+");
	LHS( x, lb, ub, n_accepted, d, &sampling_seed);
	count_evals = n_accepted*n_avg;
	for( int i=0; i<n_accepted; i++){
		for( int j=0; j<d; j++){
			//x[i][j] = lb[j] + (ub[j]-lb[j])*RAND01_R(&sampling_seed);
			fprintf( fp_population, "%e ", x[i][j]);
		}
		w[i] = 1./n_accepted;
		summary_statistics( func, func_data, d, x[i], n_avg, f[i], s2[i]);
		fprintf( fp_population, "%e %e %i %i %e\n", f[i], s2[i], 0, count_evals, DBL_MAX);
	}

	n += n_accepted;
	fprintf( fp_population, "\n\n");
	fclose(fp_population);

	sprintf( filename, "%s.prediction.dat", suffix);
	fp_population = fopen( filename, "w+");
	fclose(fp_population);

	int      n_kde = n_accepted/2;
	double **x_kde = (double**)malloc(sizeof(double*)*n_kde);
	double  *f_kde = (double*)malloc(sizeof(double)*n_kde);
	double  *w_kde = (double*)malloc(sizeof(double)*n_kde);


	fprintf( stderr, " iteration  #fun eval acceptance    epsilon\n");

	for( int it=0; it<iterations; it++){

		count_evals_per_iteration[it] = 0;

		// Prepare next iteration
		double **x_new = &x[n];
		double **x_new_tf = &x_tf[n];
		double  *f_new = &f[n];
		double  *s2_new= &s2[n];
		double  *w_new = &w[n];
		int      n_new = 0;

		// KDE stuff
		double epsilon_kde = quantile( &f[n-n_accepted], n_accepted, n_kde/(double)n_accepted);
		/*double **x_kde = &x[n-n_accepted];
		double  *f_kde = &f[n-n_accepted];
		double  *w_kde = &w[n-n_accepted];*/
		int i_kde=0;
		for( int i=0; i<n_accepted && i_kde<n_kde; i++)
		if( f[n-n_accepted+i] <= epsilon_kde){
			f_kde[i_kde] = f[n-n_accepted+i];
			w_kde[i_kde] = w[n-n_accepted+i];
			x_kde[i_kde] = x[n-n_accepted+i];
			i_kde++;
		}
		normalizeVector( w_kde, n_kde);
		//getKDESample(x,f,n, x_kde,f_kde,n_kde);
		// covariance
		_mean( x_kde, n_kde, d, mean);
		_cov(  x_kde, n_kde, d, mean, cov);
		// determinant
		detcov = determinantLU( cov, d); //fprintf( stderr, "|cov| = %e\n", detcov);
		// sqrt of coveriance
		decompCholeskyOwn(cov, covL, d); // Cholesky-decomposition
		// inverse of covariance
		inverseLU( cov, invcov, d);


		// GP hyperparameters
		int      n_gp = n;//min( n, 5*n_accepted);
		double **x_gp = &x[n-n_gp];
		double  *f_gp = &f[n-n_gp];
		double  *s2_gp= &s2[n-n_gp];


		double hyp[4], *phyp = hyp;
		//getHyperParameters( x_gp, f_gp, n_gp, d, phyp);
		int n_lim = n;//min( n_gp, 3*n_accepted);
		if( TRANSFORM){
			inverseLU( covL, invcovL, d);
			transform( x_gp, x_tf, n_gp, d, mean, invcovL);
			getHyperParameters( &x_tf[n_gp-n_lim], &f_gp[n_gp-n_lim], n_lim, d, phyp);
		}else
			getHyperParameters( &x_gp[n_gp-n_lim], &f_gp[n_gp-n_lim], n_lim, d, phyp);
		printVector( hyp, 4, "%.2e ");

		// K
		double **K = allocMatrix<double>( n_gp, n_gp);
		if( TRANSFORM)
			covMatrix<double>( K, x_tf, n_gp,d, covMatern5, phyp);
		else
			covMatrix<double>( K, x_gp, n_gp,d, covMatern5, phyp);
		for(int i=0; i<n_gp; i++)
			K[i][i] += s2_gp[i] + pow(10,phyp[0]);

		// inverse of K
		double **invK = allocMatrix<double>( n_gp, n_gp);
		inverseLU( K, invK, n_gp);



		// Threshold
		double f_kde_gp[n_kde];
		double s2_kde_gp[n_kde];
		if( TRANSFORM)
			evalVariance<double>( x_tf, f_gp, s2_gp, n_gp, x_tf, f_kde_gp, s2_kde_gp, n_kde, d, covMatern5, phyp, invK);
		else
			evalVariance<double>( x_gp, f_gp, s2_gp, n_gp, x_gp, f_kde_gp, s2_kde_gp, n_kde, d, covMatern5, phyp, invK);
		double epsilon_gp = quantile( f_kde_gp, n_kde, 0.1);

		//epsilon = quantile( f_kde, n_kde, 0.5);


		while( n_new<n_accepted){

			// KERNEL DENSITY ESTIMATE

			// choice point
			int idx;
			do{
				idx = cumsum( w_kde, n_kde, RAND01_R(&sampling_seed));
			}while( f_kde[idx] > epsilon_kde);

			// sample from MVnorm dist around chosen point
			do{
				mvnrnd( x_new[n_new], x_kde[idx], covL, d, &sampling_seed);
			}while( !inbound( x_new[n_new], lb, ub, d));


			// EVALUATION OF GP
			//double x_new_tf[d], *px_new_tf=x_new_tf;

			if( TRANSFORM){
				transform( &x_new[n_new], &x_new_tf[n_new], 1, d, mean, invcovL);
				//evalVariance<double>( x_tf, f_gp, s2_gp, n_gp, &x_new_tf[n_new], &f_new[n_new], &s2_new[n_new], 1, d, covMatern5, phyp, invK);
				double dummy;
				evalVariance<double>( x_tf, f_gp, s2_gp, n_gp, &x_new_tf[n_new], &f_new[n_new], &dummy, 1, d, covMatern5, phyp, invK);
				evalVariance<double>( x_tf, s2_gp, 0,    n_gp, &x_new_tf[n_new], &s2_new[n_new], &dummy, 1, d, covMatern5, phyp, invK);
			}else
				evalVariance<double>( x_gp, f_gp, s2_gp, n_gp, &x_new[n_new], &f_new[n_new], &s2_new[n_new], 1, d, covMatern5, phyp, invK);


			// SAMPLE FROM GP
			double rndval = f_new[n_new] + normrnd( &sampling_seed) * sqrt( s2_new[n_new]);


			// ADD POINT
			//fprintf( stderr, "WEIGHTS   THREAD %i\n", omp_get_thread_num());
			if( rndval < epsilon_gp){
				// EVALUATION OF MODEL
				summary_statistics( func, func_data, d, x_new[n_new], n_avg, f_new[n_new], s2_new[n_new]);
				count_evals_per_iteration[it] += n_avg;
				count_evals += n_avg;

				// WEIGHT
				w_new[n_new] = 1 / kdepdf( x_new[n_new], x_kde, invcov, detcov, d, n_kde);
				//w_new[n_new] = 1;
				n_new++;
			}
		}
		normalizeVector( w_new, n_new);

		if( d==2){
			sprintf( filename, "%s.prediction.dat", suffix);
			fp_population = fopen( filename, "a+");
			double _x[2], *_px=_x, _f, _s2;
			for( _x[0]=lb[0]; _x[0]<=ub[0]; _x[0]+=0.1){
			for( _x[1]=lb[1]; _x[1]<=ub[1]; _x[1]+=0.1)
			{
				for( int j=0; j<d; j++){
					fprintf( fp_population, "%e ", _x[j]);
				}

				if( TRANSFORM){
					double _xt[d], *_pxt=_xt;
					transform( &_px, &_pxt, 1, d, mean, invcovL);
					/*for( int j=0; j<d; j++){
						fprintf( fp_population, "%e ", _xt[j]);
					}*/
					evalVariance<double>( x_tf, f_gp, s2_gp, n_gp, &_pxt, &_f, &_s2, 1, d, covMatern5, phyp, invK);
				}else
					evalVariance<double>( x_gp, f_gp, s2_gp, n_gp, &_px, &_f, &_s2, 1, d, covMatern5, phyp, invK);

				fprintf( fp_population, "%e %e %i %e\n", _f, _s2, count_evals, epsilon_kde);
			}
			fprintf( fp_population, "\n");
						}
			fprintf( fp_population, "\n");
			fclose(fp_population);
		}


		freeMatrix( K,    n_gp);
		freeMatrix( invK, n_gp);

		n += n_new;

		fprintf( stderr, "%10i %10i %9.3f%% %10.3e\n", it, count_evals, 100*(double)n_new/count_evals_per_iteration[it], epsilon_kde);



		sprintf( filename, "%s.population.dat", suffix);
		fp_population = fopen( filename, "a+");
		for( int i=0; i<n_accepted; i++)
		if( f_new[i] < epsilon_kde){
			for( int j=0; j<d; j++){
				fprintf( fp_population, "%e ", x_new[i][j]);
			}
			fprintf( fp_population, "%e %e %i %i %e\n", f_new[i], s2_new[i], it+1, count_evals, epsilon_kde);
		}
		fprintf( fp_population, "\n\n");
		fclose(fp_population);


	}


	return 0;
}
