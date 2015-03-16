/*
 * Macros.hpp
 *
 *  Created on: 07.11.2013
 *      Author: jagiella
 */

#ifndef MACROS_HPP_
#define MACROS_HPP_


#ifdef NOGUI
#define RAND01  (rand()/(RAND_MAX + 1.))
#define RANDE0E1 ((rand() + 1.)/(RAND_MAX + 1.))
#else
#include <QtGlobal>
#define RAND01   (qrand()/(RAND_MAX + 1.))
#define RANDE0E1 ((qrand() + 1.)/(RAND_MAX + 1.))
#endif
#define RANDNORM( mean, variance) ( mean + variance*sqrt(-2 * log(RANDE0E1)) * cos( 2 * M_PI* RANDE0E1) )
#define MIN(a,b) (a<b?a:b)
#define MAX(a,b) (a>b?a:b)
#define MINMAX(a,b,c) (b<a ? a : (b>c ? c : b))

#define UNUSED(x) (void)(x)

template<class T>
T **allocMatrix(int m, int n) {
	T **A = (T**) malloc(sizeof(T*) * m);
	for (int i = 0; i < m; i++)
		A[i] = (T*) malloc(sizeof(T) * n);
	return A;
}

template<class T>
void freeMatrix(int m, T **A) {
	for (int i = 0; i < m; i++)
		free( A[i]);
	free( A);
}

template<class T>
void initMatrix(T **&A, int m, int n, T value) {
	for (int i = 0; i < m; i++)
		for (int j = 0; j < n; j++)
			A[i][j] = value;
}

template<class T>
void printMatrix(T **A, int m, int n) {
	for (int i = 0; i < m; i++) {
		for (int j = 0; j < n; j++)
			printf("%10.2lf", A[i][j]);
		printf("\n");
	}
}

template<class T>
void printMatrix(T *A, int m, int n) {
	for (int i = 0; i < m; i++) {
		for (int j = 0; j < n; j++)
			printf("%10.2lf", A[i * m + j]);
		printf("\n");
	}
}



#endif /* MACROS_HPP_ */

