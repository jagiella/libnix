/*
 * SparseMatrix.tcc
 *
 *  Created on: Oct 27, 2010
 *      Author: jagiella
 */

#include <stdlib.h>
#include <stdio.h>
#include <assert.h>

//#include "Mathematix.h"
#include "SparseMatrix.hpp"

#ifndef NULL
	#define NULL 0
#endif

template <class T> SparseMatrix<T>::SparseMatrix( int m, int n)
{
	A = (SparseVector<T>**) malloc( m * sizeof(SparseVector<T>**)); // values
	assert( A);

	for( int i=0; i<m; i++){
		A[i] = new SparseVector<T>( n); // values
		assert( A[i]);
	}

	this->m = m;
	this->n = n;
}


template <class T> SparseMatrix<T>::~SparseMatrix()
{
	for( int i=0; i<m; i++){
		delete A[i];
	}

	free( A);
}

template <class T> T SparseMatrix<T>::get( int i, int j)
{
	return A[i]->get( j);
}

template <class T> T& SparseMatrix<T>::getNonZero( int i, int jj)
{
	return A[i]->getNonZero( jj);
}

template <class T> int SparseMatrix<T>::getNonZeroIndex( int i, int jj)
{
	return A[i]->getNonZeroIndex( jj);
}

template <class T> void SparseMatrix<T>::resetRow( int i)
{
	A[i]->reset();
}


template <class T> void SparseMatrix<T>::setLast( int i, int j, T value)
{
	if( value == 0.) return;

	A[i]->set( j, value);
}


template <class T> void SparseMatrix<T>::set( int i, int j, T value)
{
	A[i]->set(j, value);
}

template <class T> void SparseMatrix<T>::setNonZero( int i, int jj, T value)
{
	A[i]->setNonZero(jj, value);
}


template <class T> void SparseMatrix<T>::add( int i, int j, T value)
{
	A[i]->add( j, value);
}


template <class T> void SparseMatrix<T>::MatrixVectorProduct( SparseMatrix<T> *A, T *b, T *x)
{
	//x = A*b
	// (n x m) * (m x p) = (n x p)
	// (i x j) * (j x 1) = (i x 1)
	for( int i=0; i<A->m; i++){
		x[i] = 0.;
		for( int jj=0; jj<A->A[i]->lengthNonZero(); jj++){
			int j = A->A[i]->getNonZeroIndex(jj);
			x[i] += A->A[i]->getNonZero(jj) * b[j];
		}
	}
}

template <class T> void SparseMatrix<T>::VectorMatrixProduct( T *b, SparseMatrix<T> *A, T *x)
{
	//x = b*A
	// (n x m) * (m x p) = (n x p)
	// (1 x i) * (i x j) = (1 x j)
	//int jj = 0;
	for( int j=0; j<A->dimI; j++)
		x[j] = 0.;

	for( int i=0; i<A->dimI; i++){
		for( int jj=0; jj<A->sizeA[i]; jj++){
			int j = A->JA[i][jj];
			x[j] += b[i] * A->A[i][jj];
#if DEBUG > 0
			if( isnan(x[j])){
				fprintf( stderr, "nan occures in vectorsparseMatrixProduct\n");
				fprintf( stderr, "vector b[%i] = %lf\n", i, b[i]);
				for( int jj=0; jj<sA->sizeA[i]; jj++){
					j = sA->JA[i][jj];
					fprintf( stderr, "matrix A[%i][%i] = %lf\n", i, j, sA->A[i][jj]);
				}
				exit( 0);
			}
#endif
		}
	}
}


template <class T> void SparseMatrix<T>::printMatrix( const char* name, const char* format)
{
	fprintf( stderr, "\"%s\" = [\n", name);
	for( int i=0; i<m; i++)
	{
		for( int j=0; j<n; j++)
		{
			fprintf( stderr, format, this->get(i,j));
		}
		fprintf( stderr, "\n");
	}
	fprintf( stderr, "]\n");
}


template <class T> int SparseMatrix<T>::columns()
{
	return n;
}

template <class T> int SparseMatrix<T>::columnsNonZero( int i)
{
	return A[i]->lengthNonZero();
}


template <class T> int SparseMatrix<T>::rows()
{
	return m;
}

template <class T>
T& SparseMatrix<T>::operator() ( int i, int j)
{
	return (*A[i])(j);//this->get( i, j);
}

template <class T> void SparseMatrix<T>::operator() ( int i, int j, T value)
{
	this->set( i, j, value);
}

template <class T> SparseMatrix<T>* SparseMatrix<T>::copy()
{
	SparseMatrix<T>* copy = new SparseMatrix<T>( m, n);
	for( int i=0; i<m; i++)
		for( int jj=0; jj<A[i]->lengthNonZero(); jj++){
			copy->setLast( i, A[i]->getNonZeroIndex(jj), A[i]->getNonZero(jj));
		}

	return copy;
}
