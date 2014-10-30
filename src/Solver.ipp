#include <stdlib.h>
#include <stdio.h>

//#ifdef OMP
  #include <omp.h>
//#else
//  #warning "no openmp support"
//#endif


#include "SparseMatrix.hpp"
#include "Solver.hpp"
#include "Mathematix.hpp"

#define ABS( x) (x<0?-x:x)

/*template <class T> Solver<T>::Solver( int matrix_size) : Solver(matrix_size, BiCGSTAB)
{
	//this( matrix_size, Solver::BiCGSTAB);
}*/

template <class T> Solver<T>::Solver( int matrix_size, char type, T error, int iterations)
{
	this->type = type;
	this->size = matrix_size;
	this->maxit = iterations;
	this->maxerr = error;

	switch( type){
	case Solver::SOR:
		this->omega = 0.5;
		vectorCount = 1;
		break;
	case Solver::GaussSeidel:
	case Solver::Jacobi:
		vectorCount = 1;
		break;

	case Solver::CG:
		vectorCount = 4;
		break;
	case Solver::BiCGSTAB:
		vectorCount = 7;
		break;
	case Solver::BiCGSTABl:
		l=2;
		vectorCount = 3+2*l;
		break;
	case Solver::ILU:
		vectorCount = 1;
		break;
	default:
		fprintf( stderr, "Solver of type %i still not defined!\n", type);
		break;
	}

	vectors = (T**)malloc( vectorCount*sizeof(T*));
	for( int v=0; v<vectorCount; v++)
		vectors[v] = (T*)malloc( size*sizeof(T));
}

template <class T> Solver<T>::~Solver()
{
	for( int v=0; v<vectorCount; v++)
		free( vectors[v]);
	free(vectors);
}

template <class T> int Solver<T>::solve( SparseMatrix<T> *A, T *b, T *x)
{
	switch( type){
	case Solver::GaussSeidel:
		this->SolveGaussSeidelMethod( A, b, x);
		break;
	case Solver::Jacobi:
		this->SolveJacobiMethod( A, b, x);
		break;
	case Solver::SOR:
		this->SolveSuccessiveOverRelaxationMethod( A, b, x, omega);
		break;

	case Solver::BiCGSTAB:
		return this->SolveBiCGSTAB( A, b, x);
		break;

	case Solver::BiCGSTABl:
		return this->SolveBiCGSTABl( A, b, x);
		break;

	case Solver::CG:
		return this->SolveConjugateGradient( A, b, x);
		break;
	case Solver::ILU:
		IncompleteLUfactorization(A);
		return this->SolveLU( A, b, vectors[0], x);
	default:
		fprintf( stderr, "Solver of type %i still not defined!\n", type);
		break;
	}

	return 0;
}



// o = s1*v1 + s2*v2 + s3*v3
template <class T> void Solver<T>::calculateScaledSum( T *o, T s1, T *v1, T s2, T *v2, T s3, T *v3, int N)
{
	for( int i=0; i<N; i++)
		o[i] = s1*v1[i] + s2*v2[i] + s3*v3[i];
}

// o = s1*v1 + s2*v2
template <class T> void Solver<T>::calculateScaledSum( T *o, T s1, T *v1, T s2, T *v2, int N)
{
	for( int i=0; i<N; i++)
		o[i] = s1*v1[i] + s2*v2[i];
}

// x = x - (omega*w + alpha*p)
template <class T> void Solver<T>::calculate_x( T *x, T omega, T *w, T alpha, T *p, int N)
{
	for( int i=0; i<N; i++)
		x[i] -=	(omega*w[i] + alpha*p[i]);
}

// r = r - (omega*v + alpha*s )
template <class T> void Solver<T>::calculate_r( T *r, T omega, T *v, T alpha, T *s, int N)
{
	for( int i=0; i<N; i++)
		r[i] -=	(omega*v[i] + alpha*s[i]);
}

// p = r - beta*(p - omega*s)
template <class T> void Solver<T>::calculate_p( T *x, T *y, T alpha, T *a, T beta, T*b, int N)
{
	for( int i=0; i<N; i++)
		x[i] = y[i] + alpha*a[i] + beta*b[i];
}

// p = r - beta*(p - omega*s)
template <class T> void Solver<T>::calculate_p1( T *x, T *y, T alpha, T *a, int N)
{
	for( int i=0; i<N; i++)
		x[i] = y[i] + alpha*a[i];
}

#include <stdlib.h>
#include <time.h>
#include <math.h>

template <class T> int Solver<T>::SolveBiCGSTAB( SparseMatrix<T> *A, T *b, T *x)
// If A is symmetric positive matrix
{
	/*char filename[512];
	FILE *fp;

	//sprintf( filename, "convergenceBiCGSTAB_%iD_%iPoints.dat", DIMENSIONS, A->dimI);
	sprintf( filename, "convergenceBiCGSTAB%iPoints.dat", A->dimI);
	fp = fopen( filename, "w+");
	long passedTime = clock();
	 */
	int N = A->columns();
	//fprintf( stderr, "test\n");
	// vectors
	T *r0 = vectors[0];
	T *r  = vectors[1];

	T *p  = vectors[2];

	T *s  = vectors[3];

	T *v  = vectors[4];
	T *w  = vectors[5];

	T *temp = vectors[6];
	//fprintf( stderr, "test end\n");

	// scalars
	T alpha;
	T beta;
	T omega;
	T rho;
	T sigma;


	// initialization of residual vector: r = b-Ax
	//fprintf( stderr, "r = Ax\n");
	SparseMatrix<T>::MatrixVectorProduct( A, x, temp);
	//fprintf( stderr, "...finished\n");
	vectorDifference( temp, b, r, N);

	// y, p, v
	for( int i=0; i<N; i++){
		r0[i] = r[i];
		p[i]  = r[i];
		//v[i] = 0.;
	}

	// rho
	rho = dotProduct( r, r0, N);

	// omega
	//omega = 1;

	// alpha
	//alpha = 1;

	for( int i=0; i<maxit; i++){
		//fprintf( stderr, "Inner CG-Iteration %i\n", i+1);

		// v = Ap
		//fprintf( stderr, "s = Ap\n");
		SparseMatrix<T>::MatrixVectorProduct( A, p, s);
		//fprintf( stderr, "...finished\n");

		// sigma
		//beta = rho;
		sigma = dotProduct( s, r0, N);
		if(sigma == 0){
			fprintf( stderr, "sigma = %lf, it = %i\n", sigma, i+1);
			return i;
		}

		// alpha
		alpha = rho / sigma;

		// w = r - alpha*s
		Solver<T>::calculate_p1( w, r, -alpha, s, N);

		// v = Ap
		//fprintf( stderr, "v = Aw\n");
		SparseMatrix<T>::MatrixVectorProduct( A, w, v);
		//fprintf( stderr, "...finished\n");

		// omega
		omega = dotProduct( v, w, N) / dotProduct( v, v, N);

		// x = x - (omega*w + alpha*p)
		Solver<T>::calculate_x( x, omega, w, alpha, p, N);

		// r = r - (omega*v + alpha*s )
		Solver<T>::calculate_r( r, omega, v, alpha, s, N);

		// rho
		beta = rho; // backup old rho -> beta
		rho = dotProduct( r, r0, N);

		// quadratic error
		//float myError = sqrt( dotProduct( r, r, N));

		// maximal error
		T myError = vectorNormInfinity( r, N);;
		//T myError = vectorNorm( r, N, 2);
		//double myError2 = 0;
		/*for(int m=0; m<N; m++){
			if( myError < ABS(r[m]))
				myError = ABS(r[m]);
#if COMPARE_WITH_EXACT_SOLUTION
			if( this->useExactSolution)
			if( myError2 < fabs(exactSolution[m]-x[m]))
				myError2 = fabs(exactSolution[m]-x[m]);
#endif
		}*/

//		fprintf( fp, "%i %e %e %e \n", i, (float)(clock() - passedTime) / (float)CLOCKS_PER_SEC, myError, myError2);

		//fprintf( stderr, "\r\t\t\t\t\t\t\t\t\t\t  CG: error = %e \b", myError);
		if( myError < maxerr){
			//fprintf( stderr, "ConjugateGradientSparse: reached error %e < %e after %i iterations\n", myError, maxerr, i+1);
			return i+1;// i;
		}

		// new search direction
		// beta
		beta = rho / beta * alpha / omega;

		// p = r + beta*(p - omega*s)
		Solver<T>::calculate_p( p, r, beta, p, -beta*omega, s, N);

		//fprintf( stderr, "\rCG: %i its, error = %10.3e, alpha:%10.3e, beta:%10.3e, omega:%10.3e, rho:%10.3e, sigma:%10.3e \n", i, myError, alpha, beta, omega, rho, sigma);
		//if( myError > 10e5)
		//	exit(0);
	}

	//fprintf( stderr, "error is still %e (=sqrt(%e)) < %e\n", myError, dotProduct( r, r, N), err);

//	fclose(fp);
	//return iterations;
	return maxit;
}

template <class T> void Solver<T>::setIterations( int maximalIterations)
{	this->maxit = maximalIterations;

}
template <class T> void Solver<T>::setError( T maximalResidual)
{	this->maxerr = maximalResidual;}

template <class T> void Solver<T>::setDegree( int newl)
{
	assert( type==Solver::BiCGSTABl); // only for BiCGSTAB(l)
	assert( 0==	newl%2); // must be even
	assert( newl>0); // must be larger than zero


	if( newl > l){
		// allocate new vectors
		vectorCount = 3+2*newl;
		vectors = (T**)realloc( vectors, vectorCount*sizeof(T*));
		for( int v=3+2*l; v<vectorCount; v++)
			vectors[v] = (T*)malloc( size*sizeof(T));
		l = newl;
	}

	if( newl < l){
		// free memory
		vectorCount = 3+2*newl;
		for( int v=vectorCount; v<3+2*l; v++)
			free( vectors[v]);
		vectors = (T**)realloc( vectors, vectorCount*sizeof(T*));
		l = newl;
	}
}

//#define RAND01 ((double)rand()/(RAND_MAX+1.))
template <class T> int Solver<T>::SolveBiCGSTABl( SparseMatrix<T> *A, T *b, T *x0)
// If A is symmetric positive matrix
{
	int N = A->columns();

	T error = this->maxerr;

	// VECTORS: 5+2*l
	// 3 (b inclusive)
	T *r_tilde0 = vectors[0];
	T *x_hat0 = x0;
	// 2*l+2
	T *r_hat[l+1];
	T *u_hat[l+1];
	for( int i=0; i<l+1; i++){
		r_hat[i] = vectors[2*i+1];
		u_hat[i] = vectors[2*i+2];
	}

	// SCALARS
	T omega=1, alpha=0, beta, gamma, rho0=1, rho1;
	T Sigma[l+1], Gamma[l+1], Gamma_prime[l+1], Gamma_prime2[l+1];
	T Tau[l][l+1];

	// INIT
	SparseMatrix<T>::MatrixVectorProduct( A, x0, r_hat[0]);
	vectorDifference( b, r_hat[0], r_hat[0], N);

	// u & r_tilde0
	for( int i=0; i<N; i++){
		r_tilde0[i] = unifrnd<double>();//r_hat[0][i];
		u_hat[0][i]=0;
	}

	int k=-l;

	int it=0;
	for( it=0; error >= maxerr && it<maxit; it++){
		k=k+l;
		//u_hat[0] = u;
		//r_hat[0] = r;
		//x_hat0 = x;

		rho0 = -omega*rho0;

		// Bi-CG Part
		for( int j=0; j<=l-1; j++){
			rho1 = dotProduct( r_hat[j], r_tilde0, N);
			beta = alpha*rho1/rho0;
			rho0 = rho1;

			for( int i=0; i<=j; i++){
				calculateScaledSum( u_hat[i], 1., r_hat[i], -beta, u_hat[i], N);
			}

			SparseMatrix<T>::MatrixVectorProduct( A, u_hat[j], u_hat[j+1]);
			gamma = dotProduct( u_hat[j+1], r_tilde0, N);
			alpha = rho0/gamma;

			for( int i=0; i<=j; i++){
				calculateScaledSum( r_hat[i], 1., r_hat[i], -alpha, u_hat[i+1], N);
			}
			SparseMatrix<T>::MatrixVectorProduct( A, r_hat[j], r_hat[j+1]);
			calculateScaledSum( x_hat0, 1., x_hat0, +alpha, u_hat[0], N);
			//fprintf( stderr, "omega, alpha, beta, gamma, rho0, rho1\n%e %e %e %e %e %e\n", omega, alpha, beta, gamma, rho0, rho1);
		}

		// (mod.G-S) (MR part)
		for( int j=1; j<=l; j++){
			for( int i=1; i<=j-1; i++){
				Tau[i][j] = 1/Sigma[i] * dotProduct( r_hat[j], r_hat[i], N);
				calculateScaledSum( r_hat[j], 1., r_hat[j], -Tau[i][j], r_hat[i], N);
			}
			Sigma[j] = dotProduct( r_hat[j], r_hat[j], N);
			Gamma_prime[j] = 1/Sigma[j] * dotProduct( r_hat[0], r_hat[j], N);
		}

		Gamma[l] = Gamma_prime[l];
		omega = Gamma[l];
		for( int j=l-1; j>=1; j--){
			Gamma[j] = Gamma_prime[j];
			for( int i=j+1; i<=l; i++)
				Gamma[j] -= Tau[j][i] * Gamma[i];
		}
		for( int j=1; j<=l-1; j++){
			Gamma_prime2[j] = Gamma[j+1];
			for( int i=j+1; i<=l-1; i++)
				Gamma_prime2[j] += Tau[j][i] * Gamma[i+1];
		}

		// Update
		calculateScaledSum( x_hat0, 1., x_hat0,  Gamma[1],       r_hat[0], N);
		calculateScaledSum( r_hat[0], 1., r_hat[0], -Gamma_prime[l], r_hat[l], N);
		calculateScaledSum( u_hat[0], 1., u_hat[0], -Gamma[l],       u_hat[l], N);
		for( int j=1; j<=l-1; j++){
			//fprintf( stderr, "Gamma: %e %e %e\n", Gamma[j], Gamma_prime[j], Gamma_prime2[j]);
			calculateScaledSum( u_hat[0],    1.,u_hat[0],    - Gamma[j],       u_hat[j],    N);
			calculateScaledSum( x_hat0,      1.,x_hat0,      + Gamma_prime2[j],r_hat[j],    N);
			calculateScaledSum( r_hat[0],    1.,r_hat[0],    - Gamma_prime[j], r_hat[j],    N);
		}
//		exit(0);

		//u = u_hat[0];
		//r = r_hat[0];
		//x = x_hat0;
		//error = dotProduct( r_hat[0], r_hat[0], N);
		error = vectorNormInfinity( r_hat[0], N);
		//error = vectorNorm( r_hat[0], N, 2);

	}

	//x0=x_hat0;

	return it;
}


template <class T> void Solver<T>::PreconditionJacobi( SparseMatrix<T> *A, T *b)
{
#pragma omp parallel for
	for( int i=0; i<A->rows(); i++){
		T Aii = A->get(i,i);
		for( int jj=0; jj<A->columnsNonZero(i)/*sizeA[i]*/; jj++){
			//int j = A->JA[i][jj];
			A->getNonZero( i, jj) /= Aii;/*A[i][jj] /= Aii;*/
		}
		b[i] /= Aii;
	}
}
template <class T> void Solver<T>::PreconditionSOR( SparseMatrix<T> *A, T *b, T omega)
{
	T _omega = 1.-omega;
	for( int i=0; i<A->dimI; i++){
		double Aii = 1.;
		double sigma = 0.;
		for( int jj=0; jj<A->sizeA[i]; jj++){
			int j=A->JA[i][jj];
			if( j!=i)
				sigma += A->A[i][jj];
			else
				Aii = A->A[i][jj];
		}

		for( int jj=0; jj<A->sizeA[i]; jj++){
			A->A[i][jj] *= (_omega*sigma + (b[i]-sigma) * omega / Aii);
		}
		b[i] *= (_omega*sigma + (b[i]-sigma) * omega / Aii);
	}
}


template <class T> int Solver<T>::SolveJacobiMethod( SparseMatrix<T> *A, T *b, T *x0)
{
	T *x1 = vectors[0];
	T *temp;
	T err = this->maxerr;

	int m=0;
	for( ; m<maxit && err >= maxerr; m++){
		for( int i=0; i<A->rows(); i++){
			x1[i] = 0.;
			for( int jj=0; jj<A->columnsNonZero(i); jj++){
				int j=A->getNonZeroIndex(i,jj);
				if( j!=i)
					x1[i] += A->getNonZero(i,jj) * x0[j];
			}
			x1[i] = (b[i]-x1[i]) / A->get(i,i);
			// max error
			err = ( err>fabs(x1[i]-x0[i]) ? err : fabs(x1[i]-x0[i]));
		}
		temp = x0;
		x0 = x1;
		x1 = temp;
	}

	return m;
}


template <class T> int Solver<T>::SolveGaussSeidelMethod( SparseMatrix<T> *A, T *b, T *x0)
{
	T *x1 = vectors[0];
	T *temp;
	T err = maxerr;

	int m=0;
	for( ; m<maxit && err >= maxerr; m++){
		for( int k=0; k<A->rows(); k++){
			//x0[k] = x1[k];

			x1[k] = 0.;
			for( int ii=0; ii<A->columnsNonZero(k); ii++){
				int i=A->getNonZeroIndex(k,ii);
				if( i<k)
					x1[k] += A->getNonZero(k,ii) * x1[i];
				if( i>k)
					x1[k] += A->getNonZero(k,ii) * x0[i];
			}

			x1[k] = (b[k]-x1[k]) / A->get(k,k);
			err = ( err>fabs(x1[k]-x0[k]) ? err : fabs(x1[k]-x0[k]));
		}
		temp = x0;
		x0 = x1;
		x1 = temp;
	}

	return m;
}


template <class T> int Solver<T>::SolveSuccessiveOverRelaxationMethod( SparseMatrix<T> *A, T *b, T *x0, T omega)
{
	//char filename[512];
	//FILE *fp;

	/*sprintf( filename, "convergenceSOR_%iD_%iPoints.dat", DIMENSIONS, A->dimI);
	fp = fopen( filename, "w+");
	long passedTime = clock();
*/
	T *x1 = vectors[0];
	T *temp;
	T err = maxerr;

	T _omega = 1.-omega;

	int m=0;
	for( ; m<maxit && err >= maxerr; m++){
		for( int i=0; i<A->rows(); i++){
			//x0[k] = x1[k];
			T Aii = 1.;
			x1[i] = 0.;
			for( int jj=0; jj<A->columnsNonZero(i); jj++){
				int j=A->getNonZeroIndex(i,jj);
				if( j<i)
					x1[i] += A->getNonZero(i,jj) * x1[j];
				else if( j==i)
					Aii    = A->getNonZero(i,jj);
				else if( j>i)
					x1[i] += A->getNonZero(i,jj) * x0[j];
			}

			//x1[k] = (1.-omega)*x0[k] + (b[k]-x1[k]) * omega / A->get(k,k);
			x1[i] = _omega*x0[i] + (b[i]-x1[i]) * omega / Aii;
			err = ( err>fabs(x1[i]-x0[i]) ? err : fabs(x1[i]-x0[i]));
		}
		temp = x0;
		x0 = x1;
		x1 = temp;

		// maximal error
		//double myError = 0;
		//double myError2 = 0;
		for(int n=0; n<A->rows(); n++){
			//if( myError < fabs(r[m]))
			//	myError = fabs(r[m]);
#if COMPARE_WITH_EXACT_SOLUTION
			if( myError2 < fabs(exactSolution[n]-x0[n]))
				myError2 = fabs(exactSolution[n]-x0[n]);
#endif
		}

	/*	fprintf( fp, "%i %e %e %e \n", m, (float)(clock() - passedTime) / (float)CLOCKS_PER_SEC, myError, myError2);
*/	}

	//fclose(fp);

	return m;
}


template <class T> int Solver<T>::SolveExplicit( SparseMatrix<T> *A, T *b, T *x)
{
#pragma omp parallel for
	for( int i=0; i<A->dimI; i++)
		x[i] = b[i] / A->A[i][0];

	return A->dimI;
}

template <class T> int Solver<T>::SolveConjugateGradient( SparseMatrix<T> *A, T *b, T *x)
// If A is symmetric positive matrix
{
	/*char filename[512];
	FILE *fp;

	//sprintf( filename, "convergencePreCondCG_%iD_%iPoints.dat", DIMENSIONS, A->dimI);
	sprintf( filename, "convergencePreCondCG_%iPoints.dat", A->dimI);
	fp = fopen( filename, "w+");
	long passedTime = clock();
	*/
	int N = A->columns();

	// vectors
	T *r = vectors[0];
	T *w = vectors[1];
	T *z = vectors[2];

	// scalars
	T alpha;
	T beta;

	T *temp = vectors[3];

	//float err = 1e-10;
	T myError = 0.;
	//fprintf( stderr, "1\n");
	// initialization of residual vector: r
	SparseMatrix<T>::MatrixVectorProduct( A, x, temp);
	vectorDifference( b, temp, r, N);

	// w
	vectorScale( r, -1, w, N);
	//fprintf( stderr, "2\n");
	// z
	SparseMatrix<T>::MatrixVectorProduct( A, w, z);

	// alpha
	alpha = dotProduct( r, w, N) / dotProduct( w, z, N);

	// beta
	beta = 0.;

	// x
	vectorScale( w, alpha, temp, N);
	vectorSum( x, temp, x, N);

	int i=0;
	for( ; i<maxit; i++){
		//fprintf( stderr, "Inner CG-Iteration %i\n", i+1);

		// r = r - alpha*z
		vectorScale( z, alpha, temp, N);
		vectorDifference( r, temp, r, N);

		/*float errSqr = dotProduct( r, r, N);
		if( isnan(errSqr)){
			fprintf( stderr, "ConjugateGradientSparse: ERROR is nan!\n");
			exit( 0);
		}
		myError = sqrt( errSqr);*/

		//fprintf( stderr, "\r\t\t\t\t\t\t\t\t\t\t  CG: error = %e \b", myError);

		//myError = 0;
		myError = vectorNormInfinity( r, N);
		//myError = vectorNorm( r, N, 2);

		//T myError2 = 0;
		/*for(int m=0; m<N; m++){
			if( myError < ABS(r[m]))
				myError = ABS(r[m]);
#if COMPARE_WITH_EXACT_SOLUTION
			if( this->useExactSolution)
			if( myError2 < fabs(exactSolution[m]-x[m]))
				myError2 = fabs(exactSolution[m]-x[m]);
#endif
		}*/


//		fprintf( fp, "%i %e %e %e\n", i, (float)(clock() - passedTime) / (float)CLOCKS_PER_SEC, myError, myError2);

		if( myError < maxerr){
			//fprintf( stderr, "ConjugateGradientSparse: reached error %e < %e after %i iterations\n", myError, maxerr, i+1);
			return i;// i;
		}

		// B = (r'*z)/(w'*z);
		beta = dotProduct( r, z, N) / dotProduct( w, z, N);
		//if( beta<0.) beta = 0.;

		// w = -r + B*w;
		vectorScale( w, beta, w, N);
		vectorDifference( w, r, w, N);
		//fprintf( stderr, "3\n");
		// z = A*w;
		SparseMatrix<T>::MatrixVectorProduct( A, w, z);

		// a = (r'*w)/(w'*z);
		alpha = dotProduct( r, w, N) / dotProduct( w, z, N);

		// x = x + a*w;
		vectorScale( w, alpha, temp, N);
		vectorSum( x, temp, x, N);
	}
/*
	fprintf( stderr, "error is still %e (=sqrt(%e)) < %e\n", myError, dotProduct( r, r, N), minerr);
*/
//	fclose(fp);
	return i;
}


template <class T> int Solver<T>::SolverPreconditionedConjugateGradient( SparseMatrix<T> *A, SparseMatrix<T> *M, T *b, T *x)
{
	// vectors
	//float r[N];
	//float z[N];
	//float p[N];
	int N = A->dimI;

	T *r = vectors[0];
	T *z = vectors[1];
	T *p = vectors[2];

	// scalars
	T alpha;
	T beta;
	T rho;

	//float temp[N];
	T *temp = vectors[3];

	//float err = 1e-4;
	T myError = 0.;

	// initialization:

	// r = b - A*x
	/*for( int i=0; i<N; i++)
		if( isnan(x[i])){
			fprintf( stderr, "x[%i] = %lf\n", i, x[i]);
			exit(0);
		}*/
	//fprintf( stderr, "r = b - A*x\n");
	SparseMatrix<T>::MatrixVectorProduct( A, x, temp);
	vectorDifference( b, temp, r, N);

	// z = M*r
	//fprintf( stderr, "z = M*r\n");
	SparseMatrix<T>::MatrixVectorProduct( M, r, z);

	// p = z
	vectorScale( z, 1, p, N);


	for( int i=0; i<maxit; i++){
		// rho = r*z
		rho = dotProduct( r, z, N);

		// alpha = r*z / p*A*p
		SparseMatrix<T>::vectorSparseMatrixProduct( p, A, temp);
		alpha = rho / dotProduct( temp, p, N);

		// x = x + alpha*p
		vectorScale( p, alpha, temp, N);
		vectorSum( x, temp, x, N);

		// r = r - alpha*Ap
		SparseMatrix<T>::MatrixVectorProduct( A, p, temp);
		vectorScale( temp, alpha, temp, N);
		vectorDifference( r, temp, r, N);

		myError = sqrt( dotProduct( r, r, N));
		//if( sqrt( dotProduct( r, r, N)) < err){
		fprintf( stderr, "\r\t\t\t\t\t\t\t\t\t\t  CG: error = %e \b", myError);

		if( myError < maxerr){
			//fprintf( stderr, "ConjugateGradientSparse: reached error %e < %e after %i iterations\n", myError, err, i+1);
			//fprintf( stderr, "\r\t\t\t\t\t\t\t\t\t\t  CG: error = %e \b", myError);
			return myError;
		}else if( myError > 1.){
			//fprintf( stderr, "WARNING: ConjugateGradientSparse reached error %e > 1 after %i iterations\n", myError, i+1);
		}

		// z = M^-1*r
		// Solve Mz = r
		SparseMatrix<T>::MatrixVectorProduct( M, r, z);

		// beta = r*z/rho
		beta = dotProduct( r, z, N) / rho;

		// p = z + beta*p
		vectorScale( p, beta, p, N);
		vectorSum( z, p, p, N);
	}

	fprintf( stderr, "error is still %e < %e\n", myError, maxerr);
	return myError;
}


template <class T> void Solver<T>::SuccessiveOverRelaxationPreconditioner( SparseMatrix<T> *A, SparseMatrix<T> *M)
{
	// (ixj) = (ixn) x (nxj)
	double omega = .01;

	// row
	for( int i=0; i<A->dimI; i++){
		//double factor = 0.;
		//column
		//for( int j=0; j<A->dimI; j++){
		for( int jj=0; jj<A->sizeA[i]; jj++){
			int j = A->JA[i][jj];
			double value = 0.;
			double ab = 0.;
			//double b = 0.;
			double c = 0.;

			// matrix product
			for( int pp=0; pp<A->sizeA[i]; pp++){
				int p = A->JA[i][pp];
				ab=0.;
				c=0.;
				// A*B
				if( i>p){
					// A*B = L * w/(2-w)D^-1
					ab = A->A[i][pp] * omega/(2.-omega) / A->get(p,p);
				}
				if( p==i){
					// A*B = D/w * w/(2-w)D^-1
					ab = A->A[i][pp]/(2.-omega) / A->get(p,p);
				}

				// C
				//if( ab!=0.){
				if( p<j){
					// C = U
					c = A->get(p,j) * omega/(2.-omega) / A->get(p,p);
				}
				if( p==j){
					// C = D/w
					c = A->get(p,p)/omega;
				}
				//}

				value += ab	* c;
			}

			if(value!=0.)
				M->set( i,j, value);
		}

		//M->set( i,i, A->get(i,i));
	}
//	fprintf( stderr, "SuccessiveOverRelaxationPreconditioner still not implimented\n");
//	exit( 0);
}

template <class T> void Solver<T>::IncompleteLUfactorization( SparseMatrix<T> *A)
{
	// NAIVE
	/*int n=A->dimI;
	for( int k=0; k<n-1; k++)
		for( int i=k+1; i<n; i++)
		if((*A)(i,k)!=0)
		{
			A->set( i,k, (*A)(i,k) / (*A)(k,k));
			for( int j=k+1; j<n; j++)
			if((*A)(i,j)!=0)
			{
				A->set( i,j, (*A)(i,j) - (*A)(i,k) * (*A)(k,j));
			}
		}*/


	int n=A->rows();
	for( int k=0; k<n-1; k++){
		fprintf( stderr, "\rILU: %5.2lf\%%  \b", (double)100.*k/n);
		for( int i=k+1; i<n; i++)
		{
			int jj, kk;
			for( kk=0; kk<A->columnsNonZero(i) && A->getNonZeroIndex(i,kk)<k; kk++) ;

			// if((*A)(i,k)!=0)
			if( kk<A->columnsNonZero(i) && A->getNonZeroIndex(i,kk)==k){
				// A->set( i,k, (*A)(i,k) / (*A)(k,k));
				A->getNonZero(i,kk) /= (*A)(k,k);

				for( jj=kk+1; jj<A->columnsNonZero(i); jj++){
					int j = A->getNonZeroIndex(i,jj);
					A->getNonZero(i,jj) -= A->getNonZero(i,kk) * (*A)(k,j);
				}
			}
		}
	}
	fprintf( stderr, "\rILU Factorization complete\n");
}

template <class T> void Solver<T>::IncompleteLUfactorization( SparseMatrix<T> *original, SparseMatrix<T> *copy)
{
	(*copy) = original;

	IncompleteLUfactorization( copy);
}

template <class T> int Solver<T>::SolveLowerTriangular( SparseMatrix<T> *A, T *b, T *x)
{
	// L00  0   0
	// L01 L11  0
	// ... ... L22
	// Lij = 1

	int n=A->dimI;
	x[0] = b[0] / (*A)(0,0);
	for( int i=1; i<n; i++){
		fprintf( stderr, "\rSolve Ly=b: %5.2lf\%%  \b", (double)100.*i/n);
		x[i] = b[i];
		//for( int j=0; j<i; j++)
		//	x[i] -= x[j] * (*A)(i,j);
		for( int jj=0; /*jj<A->sizeA[i] &&*/ A->JA[i][jj]<i; jj++){
			int j = A->JA[i][jj];
			x[i] -= x[j] * A->A[i][jj];
		}
		x[i] /= (*A)(i,i);
	}

	return 1;
}

template <class T> int Solver<T>::SolveUpperTriangular( SparseMatrix<T> *A, T *b, T *x)
{
	int n=A->dimI;
	x[n-1] = b[n-1] / (*A)(n-1,n-1);
	for( int i=n-2; i>=0; i--){
		fprintf( stderr, "\rSolve Ux=y: %5.2lf\%%  \b", (double)100.*(n-i)/n);
		x[i] = b[i];
		//for( int j=n-1; j>i; j--)
		//	x[i] -= x[j] * (*A)(i,j);
		for( int jj=A->sizeA[i]-1; /*jj>=0 && */A->JA[i][jj]>i; jj--){
			int j = A->JA[i][jj];
			x[i] -= x[j] * A->A[i][jj];
		}
		x[i] /= (*A)(i,i);
	}

	return 1;
}

template <class T> int Solver<T>::SolveLU( SparseMatrix<T> *LU, T *b, T *y, T *x)
{
	//  1   0   0
	// L10  1   0
	// L20 L21  1
	// ...

	// SOLVE Ly=b

	int m=LU->rows();
	y[0] = b[0];
	for( int i=1; i<m; i++){
		fprintf( stderr, "\rSolve Ly=b: %5.2lf\%%  \b", (double)100.*i/m);
		y[i] = b[i];
		for( int jj=0; LU->getNonZeroIndex( i, jj)<i; jj++){
			int j = LU->getNonZeroIndex( i, jj);
			y[i] -= y[j] * LU->getNonZero( i, jj);
		}
	}

	// U00 U01 U02 ...
	//  0  U11 U12 ...
	//  0   0  U22 ...
	// ...

	// SOLVE Ux=y

	x[m-1] = y[m-1] / (*LU)(m-1,m-1);
	for( int i=m-2; i>=0; i--){
		fprintf( stderr, "\rSolve Ux=y: %5.2lf\%%  \b", (double)100.*(m-i)/m);
		x[i] = y[i];
		for( int jj=LU->columnsNonZero(i)-1; LU->getNonZeroIndex(i,jj)>i; jj--){
			int j = LU->getNonZeroIndex(i,jj);
			x[i] -= x[j] * LU->getNonZero(i,jj);
		}
		x[i] /= (*LU)(i,i);
	}

	return 1;
}
