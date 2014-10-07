#ifndef COV_IPP

#include <math.h>
#include <assert.h>

//#include "mex.h"

template <class T>
T covSquaredExponential( T* &x1, T* &x2, int &d, T* &theta)
{
			// distance
			double dist = 0;
			for( int k=0; k<d; k++)
				dist += pow( x1[k] - x2[k], 2);
			dist = sqrt(dist);

			// covariance
			T l  = exp(theta[1]);
			T sf2= exp(theta[2]);

			assert( l>0);

			T temp = -0.5*(dist*dist) / (l*l);

			if(dist>0){
				return sf2*exp(temp);
			}else{
				return sf2;
			}
}

template <class T>
T covMatern1( T* &x1, T* &x2, int &d, T* &theta)
{
			// distance
			double dist = 0;
			for( int k=0; k<d; k++)
				dist += pow( x1[k] - x2[k], 2);
			dist = sqrt(dist);

			// covariance
			T nu = 1./2.;
			T l  = theta[1];
			T sf2= theta[2];

			assert( nu>0);
			assert( l>0);

			T temp = sqrt(2*nu)*dist / l;

			if(dist>0){
				return (1)*exp(-temp);
			}else{
				return 1.;
			}
}

template <class T>
T covMatern3( T* &x1, T* &x2, int &d, T* &theta)
{
			// distance
			double dist = 0;
			for( int k=0; k<d; k++)
				dist += pow( x1[k] - x2[k], 2);
			dist = sqrt(dist);

			// covariance
			T nu = 3./2.;
			T l  = theta[1];
			T sf2= theta[2];

			assert( nu>0);
			assert( l>0);

			T temp = sqrt(2*nu)*dist / l;

			if(dist>0){
				return (1+temp)*exp(-temp);
			}else{
				return 1.;
			}
}

template <class T>
T covMatern5( T* &x1, T* &x2, int &d, T* &theta)
{
			// distance
			double dist = 0;
			for( int k=0; k<d; k++)
				dist += pow( x1[k] - x2[k], 2);
			dist = sqrt(dist);

			// covariance
			T nu = 5./2.;
			T l  = pow(10,theta[1]);
			T sf2= pow(10,theta[2]);

			assert( nu>0);
			assert( l>0);

			T t = sqrt(2*nu)*dist / l;
			//fprintf( stderr, "covMatern5 = %e, l=%e, sf2=%e, dist=%e\n", sf2*(1 + t*(1+t/3))*exp(-t), l, sf2, dist);
			if(dist>0 && !isnan(t) && !isinf(t)){
				return sf2*(1 + t*(1+t/3))*exp(-t);
			}else{
				return sf2*1.;
			}
}

template <class T>
T covMatern( T* &x1, T* &x2, int &d, T* &theta)
{
			// distance
			double dist = 0;
			for( int k=0; k<d; k++)
				dist += pow( x1[k] - x2[k], 2);
			dist = sqrt(dist);

			// covariance
			T nu = theta[0];
			T l  = theta[1];
			T sf2= theta[2];

			assert( nu>0);
			assert( l>0);

			T temp = sqrt(2*nu*dist) / l;
			//
			/*fprintf(stderr, "%e %e %e\n",
					pow( 2, 1-nu),
					gamma(nu),
					pow(temp, nu));*/
			if(dist>0){
				return sf2*pow( 2, 1-nu) / gamma(nu) * pow(temp, nu) * gsl_sf_bessel_Knu( nu, temp);
			}else{
				return sf2*1.;
			}
}

template <class T>
void covVector( T* &Ki, T* &x, T** &X, int &n, int &d, T (*covFunc) ( T* &, T* &, int &, T* &), T* &theta)
{
	for( int j=0; j<n; j++){
		Ki[j] = covFunc(X[j], x, d, theta );
	}
}

template <class T>
void covVector( T* &Ki, T* &x, T* &X, int &n, int &d, T (*covFunc) ( T* &, T* &, int &, T* &), T* &theta)
{
	for( int j=0; j<n; j++){
		T* pRow = &X[j*d];
		Ki[j] = covFunc(pRow, x, d, theta );
	}
}

template <class T>
void covMatrix( T** &K, T** &X, int &n, int &d, T (*covFunc) ( T* &, T* &, int &, T* &), T* &theta)
{
	for( int i=0; i<n; i++){
		// covariance of row i with all other rows
		covVector<T>( K[i], X[i], X, n, d, covFunc, theta);
	}
}

template <class T>
void covMatrix( T** &K, T* &X, int &n, int &d, T (*covFunc) ( T* &, T* &, int &, T* &), T* &theta)
{
	for( int i=0; i<n; i++){
		// covariance of row i with all other rows
		T* pRow = &X[i*d];
		covVector<T>( K[i], pRow, X, n, d, covFunc, theta);
	}
}



#define COV_IPP
#endif
