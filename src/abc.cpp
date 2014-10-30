#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <float.h>
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



bool inbound( double *x, double *lb, double *ub, int d)
{
	for( int i=0; i<d; i++)
		if( x[i] < lb[i] || x[i] > ub[i])
			return false;
	return true;
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

void lin2log( double *xlin, double *xlog, int n){
	for( int i=0; i<n; i++)
		xlog[i] = log10(xlin[i]);
}

void log2lin( double *xlog, double *xlin, int n){
	for( int i=0; i<n; i++)
		xlin[i] = pow( 10, xlog[i]);
}

#include <unistd.h> // getopt

int main( int argc, char **argv){
	srand(time(NULL));
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

	double *lb=0;
	double *ub=0;

	int iCPU = 10;
	n = 100;

	char suffix[1024];
	char filename[1024];
	sprintf( suffix, "%s", "output" );

	enum {Log, Lin};
	char SCALING = Lin;

	char c;
	while ((c = getopt (argc, argv, "m:l:u:n:o:L")) != -1)
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

	      case 'n': n = atoi( optarg); break;

	      case 'o': sprintf( suffix, "%s", optarg ); break;

	      case 'L': SCALING = Log;
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





	// init variables
	fprintf( stderr, "INIT\n");
	double epsilon = DBL_MAX;
	double **x_new, **x_old, **x_tmp, **x_par, *x_lin;
	double  *f_new,  *f_old,  *f_tmp,  *f_par;
	double  *w_new,  *w_old,  *w_tmp,  *w_par;
	int      n_new,   n_old,            n_par = nparallel;
	double mean[d], **cov = allocMatrix<double>(d,d), **invcov = allocMatrix<double>(d,d), **covL = allocMatrix<double>(d,d), detcov;

	x_old = allocMatrix<double>( n+n_par, d);	f_old = (double*) malloc((n+n_par)*sizeof(double));	w_old = (double*) malloc((n+n_par)*sizeof(double)); n_old = 0;
	x_new = allocMatrix<double>( n+n_par, d);	f_new = (double*) malloc((n+n_par)*sizeof(double));	w_new = (double*) malloc((n+n_par)*sizeof(double)); n_new = 0;
	x_par = allocMatrix<double>( n_par,  d);	f_par = (double*) malloc(n_par*sizeof(double));    	w_par = (double*) malloc(n_par*sizeof(double));

	if(SCALING == Log)
		x_lin = (double*) malloc(d*sizeof(double));;

#ifdef USE_OMP
	//int iCPU = nparallel;//= omp_get_num_procs();
	fprintf( stderr, "SET NUM THREADS: %i\n", iCPU);
	omp_set_num_threads( iCPU);
#endif
	//printMatrix( x_old, n_old, d, "%10.3e ");


	// stats
	clock_t t, t_sampling = 0, t_evaluation = 0;
	sprintf( filename, "%s.population.dat", suffix);
	FILE *fp_population = fopen( filename, "w+");

	// run iteration
	int count_evals_per_iteration[iterations];
	int count_evals = 0;
	unsigned int sampling_seed = 0;

	fprintf( stderr, " iteration  #fun eval acceptance    epsilon\n");

	for( int it=0; it<iterations; it++){

		n_new = 0;

		// COPY OLD ONES < EPSILON
		/*for( int i=0; i<n_old; i++)
			if( f_old[i] < epsilon){
				for( int j=0; j<d; j++)
				x_new[n_new][j] = x_old[i][j];
				f_new[n_new] = f_old[i];
				w_new[n_new] = 1 / kdepdf( x_old[i], x_old, invcov, detcov, d, n_old);
				n_new ++;
			}*/


		// sample new points
		count_evals_per_iteration[it]=0;

		while( n_new<n){

			// >>> SEQUENTIAL PART <<<

			for( int ip=0; ip<nparallel; ip++){

				t = clock();
				//fprintf( stderr, "SAMPLING %i\n", omp_get_thread_num());
				if( it == 0){
					// UNIFORM SAMPLING
					for( int j=0; j<d; j++)
						x_par[ip][j] = lb[j] + (ub[j]-lb[j])*RAND01_R(&sampling_seed);
				}else{
					// KERNEL DENSITY ESTIMATE
					// choice point
					int idx = cumsum( w_old, n_old, RAND01_R(&sampling_seed));

					// sample from MVnorm dist around chosen point
					do{
						mvnrnd( x_par[ip], x_old[idx], covL, d, &sampling_seed);
						//printVector( lb, d, "%5.3f ");
						//printVector( x_par[ip], d, "%5.3f ");
						//printVector( ub, d, "%5.3f ");
						//printf( "%i\n", inbound( x_par[ip], lb, ub, d));
					}while( !inbound( x_par[ip], lb, ub, d));
				}
				t_sampling += clock() - t;
			}

			// >>> PARALLEL PART <<<

	#ifdef USE_OMP
	#pragma omp parallel for
	#endif
			for( int ip=0; ip<nparallel; ip++){
				// EVALUATION

				t = clock();
				//fprintf( stderr, "EVALUATION THREAD %i\n", omp_get_thread_num());
				switch( SCALING){
				case Lin:
					f_par[ip] = (*func)( d, x_par[ip], 0, func_data); break;
				case Log:
					log2lin( x_par[ip], x_lin, d);
					f_par[ip] = (*func)( d, x_lin, 0, func_data); break;
				}

				t_evaluation += clock() - t;


				// WEIGHTS
				//fprintf( stderr, "WEIGHTS   THREAD %i\n", omp_get_thread_num());
				w_par[ip]=0;
				if( f_par[ip] < epsilon){
					if( it == 0){
						w_par[ip] = 1;
					}else{
						w_par[ip] = 1 / kdepdf( x_par[ip], x_old, invcov, detcov, d, n_old);
					}
				}

				//fprintf( stderr, " <- End   THREAD %i\n", omp_get_thread_num());

			}
			count_evals += nparallel;
			count_evals_per_iteration[it] += nparallel;



			// >>> SEQUENTIAL PART <<<
			for( int ip=0; ip<nparallel; ip++){
				//fprintf( stderr, "%i = %5e \n", ip, f_par[ip]);
				//fprintf( stderr, "%s\n", (w_par[ip] ? "accepted":"refused"));
				if( w_par[ip]){
					// COPY TO NEXT GENERATION
					w_new[n_new] = w_par[ip];
					f_new[n_new] = f_par[ip];
					for( int id=0; id<d; id++)
						x_new[n_new][id] = x_par[ip][id];
					n_new++;
				}
			}

			//return 0;

			if( count_evals >= evaluations){
				//break;
				fprintf( stderr, "END\n");
				double f_min;
				int    i_min;
				minimum( f_old, n_old, f_min, i_min);
				printVector( x_old[i_min], d, "%e ");

				fprintf ( stderr,
						"TIME STATS: sampling = %e seconds per point, evaluation = %e sec / point.\n",
						( (float)t_sampling  /count_evals)/CLOCKS_PER_SEC,
						( (float)t_evaluation/count_evals)/CLOCKS_PER_SEC);

				// PRINT KDE
				sprintf( filename, "%s.kde.dat", suffix);
				FILE *fp_kde = fopen( filename, "w+");
				double dx = 0.1;
				double intergral_kde = 0;
				if(d==1){
					dx = 0.01;
					for( double x=lb[0]; x<ub[0]; x+=dx){
							double xtest[1] = {x};
							fprintf( fp_kde, "%e %e\n", x, kdepdf( xtest, x_old, invcov, detcov, d, n_old));
							intergral_kde += dx*kdepdf( xtest, x_old, invcov, detcov, d, n_old);
						}
				}

				if(d==2){
					int n0 = (int) (ub[0]-lb[0])/dx + 1;
					int n1 = (int) (ub[1]-lb[1])/dx + 1;
					double kde[n0][n1];
					double lik[n0][n1];
#ifdef USE_OMP
#pragma omp parallel for
#endif
					for( int i0=0; i0<n0; i0++){
						for( int i1=0; i1<n1; i1++){
							double x = lb[0] + dx*i0, y = lb[1] + dx*i1;
							double xtest[2] = {x, y};
							kde[i0][i1] = kdepdf( xtest, x_old, invcov, detcov, d, n_old);
							lik[i0][i1] = (*func)( d, xtest, 0, func_data);
						}
					}
					for( int i0=0; i0<n0; i0++){
						for( int i1=0; i1<n1; i1++){
							double x = lb[0] + dx*i0, y = lb[1] + dx*i1;
							fprintf( fp_kde, "%e %e %e %e\n", x,y, kde[i0][i1], lik[i0][i1]);
							intergral_kde += dx*dx*kde[i0][i1];
						}
						fprintf( fp_kde, "\n");
					}
				}
				fclose( fp_kde);
				fprintf ( stderr, "[ %e ]\n", intergral_kde);

				{
					double startingPoint[d];
					double stepSize[d];
					int startingDimension = 0;
					for( int i=0; i<d; i++)
						stepSize[i] = (ub[i]-lb[i]) / 10.;
					kdedata _kdedata;
					_kdedata.detcov = detcov;
					_kdedata.invcov = invcov;
					_kdedata.n_old  = n_old;
					_kdedata.x_old  = x_old;
					#pragma omp parallel for
					for( startingDimension=0; startingDimension<d; startingDimension++){
						sprintf( filename, "%s.kde%i.dat", suffix, startingDimension);
						FILE *fp_kde_d = fopen( filename, "w+");
						double dx = (ub[startingDimension]-lb[startingDimension]) / 10.;
						for( startingPoint[startingDimension]=lb[startingDimension]; startingPoint[startingDimension]<=ub[startingDimension]; startingPoint[startingDimension]+=dx)
							fprintf( fp_kde_d, "%e %e\n", startingPoint[startingDimension], rekursiveIntegration( startingPoint, 0, startingDimension,
									lb, ub, stepSize, d,
									myfunckde, &_kdedata));
						fclose( fp_kde_d);
					}
				}

				return 0;
			}
			//fprintf( stderr, "---------> %5i+%4i %5i / %i\n", count_evals-count_evals_per_iteration[it], count_evals_per_iteration[it], n_new, n );
//			fprintf( stderr, "%10i %10i %9.3f%% %10.3e\n", it, count_evals, 100*(double)n_new/count_evals_per_iteration[it], epsilon);
			fprintf( stderr, "\r%10i %10i [%3i/%4i] %10.3e\r", it, count_evals, n_new, n, epsilon);

		}

		normalizeVector( w_new, n_new);
		epsilon = quantile( f_new, n_new, 0.5);

		for( int i=0; i<n_new; i++){
			for( int j=0; j<d; j++)
				fprintf( fp_population, "%e ", x_new[i][j]);
			fprintf( fp_population, "%e %e %i %i %e\n", f_new[i], w_new[i], it, count_evals, epsilon);
		}
		fprintf( fp_population, "\n\n");

		//fprintf( stderr, "new epsilon = %e\n", epsilon);



		fprintf( stderr, "%10i %10i %9.3f%% %10.3e\n", it, count_evals, 100*(double)n_new/count_evals_per_iteration[it], epsilon);

		x_tmp = x_old; x_old = x_new; x_new = x_tmp;
		f_tmp = f_old; f_old = f_new; f_new = f_tmp;
		w_tmp = w_old; w_old = w_new; w_new = w_tmp;
		n_old = n_new;                n_new = 0;


		//fprintf( stderr, "[ ITERATION %i ]\n", it);
		count_evals_per_iteration[it]=0;

		//fprintf( stderr, "mean\n");
		_mean( x_old, n_old, d, mean);
		//printVector( mean,d,"%10e ");

		//fprintf( stderr, "cov\n");
		_cov(  x_old, n_old, d, mean, cov);
		//printMatrix( cov, d, d, "%10e ");

		detcov = determinantLU( cov, d); //fprintf( stderr, "|cov| = %e\n", detcov);

		//fprintf( stderr, "L\n");
		decompCholeskyOwn(cov, covL, d); // Cholesky-decomposition
		//printMatrix( covL, d, d, "%10e ");

		//fprintf( stderr, "L*L^T\n");
		/*double **test = allocMatrix<double>( d,d);
		for( int i=0; i<d; i++)
			for( int j=0; j<d; j++){
				test[i][j] = 0;
				for( int k=0; k<d; k++)
					test[i][j] += covL[i][k]*covL[j][k];
			}*/
		//printMatrix( test, d, d, "%10e ");


		//inverseCholesky( cov, invcov, d);
		inverseLU( cov, invcov, d);

	}


	fprintf( stderr, "END\n");
	double f_min;
	int    i_min;
	minimum( f_old, n_old, f_min, i_min);
	printVector( x_old[i_min], d, "%e ");

	fprintf ( stderr,
			"TIME STATS: sampling = %e seconds per point, evaluation = %e sec / point.\n",
			( (float)t_sampling  /count_evals)/CLOCKS_PER_SEC,
			( (float)t_evaluation/count_evals)/CLOCKS_PER_SEC);

	return 0;
}
