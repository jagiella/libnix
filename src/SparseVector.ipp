/*
 * SparseVector.ipp
 *
 *  Created on: 20.04.2012
 *      Author: jagiella
 */

#ifndef SPARSEVECTOR_IPP_
#define SPARSEVECTOR_IPP_

#define ROW_EXTENSION_SIZE 10

#include <assert.h>
#include <stdlib.h>

//#include "Mathematix.hpp"
#include "SparseVector.hpp"

template <class T> SparseVector<T>::SparseVector( int m)
{
	A        = (T*) malloc( ROW_EXTENSION_SIZE * sizeof(T)); // values
	JA       = (int*)   malloc( ROW_EXTENSION_SIZE * sizeof(int));  // column index of values
	size = 0;
	sizeAllocated = ROW_EXTENSION_SIZE;

	this->m = m;
	value = 0;
}

template <class T> SparseVector<T>::~SparseVector()
{
	free(A);
	free(JA);
}

template <class T> int SparseVector<T>::length()
{	return m;}

template <class T> int SparseVector<T>::lengthNonZero()
{	return size;}

template <class T> void SparseVector<T>::setNonZero( int jj, T value)
{	A[jj] = value;}

template <class T> void SparseVector<T>::set( int j, T value)
{
	if( value == 0.) return;

	int ji = 0;
	int jj = size-1;
	int k  = 0;
	// look for position
	if( size>0){
		ji = 0;
		int x = j;

		do{
			k = (ji + jj) / 2;
			if( x > JA[k])
				ji = k + 1;
			else
				jj = k - 1;
		}while( (JA[k] != x) && (ji <= jj));
	}

	if( size==0 || (JA[k] != j && jj<ji)){
		if( jj<ji)
			k = ji;
		// realloc memory
		if( size == sizeAllocated){
			sizeAllocated += ROW_EXTENSION_SIZE;
			JA = (int*)   realloc( JA, sizeAllocated * sizeof(int));
			assert( JA);
			A  = (T*) realloc( A, sizeAllocated * sizeof(T));
			assert( A);
		}
		// shift following elements
		for( int jjj=size; jjj>k; jjj--){
			JA[jjj] = JA[jjj-1];
			A[jjj] = A[jjj-1];
		}
		JA[k] = j;
		A[k] = value;
		size ++;

	}else{
		A[k] = value;
	}
	return;
}

template <class T> void SparseVector<T>::add( int j, T value)
{
	if( value == 0) return;

	int ji = 0;
	int jj = this->size-1;
	int k  = 0;
	// look for position
	if( this->size>0){
		ji = 0;
		//int ji = 0;
		//jj = this->sizeA[i]-1;
		int x = j;

		do{
			k = (ji + jj) / 2;
			//fprintf( stderr, "i=%i, j=%i, k=%i\n", ji, jj, k);
			if( x > this->JA[k])
				ji = k + 1;
			else
				jj = k - 1;
		}while( (this->JA[k] != x) && (ji <= jj));
		// jj > ji =>
	}
	//fprintf( stderr, "(i, j) = (%i, %i), ji=%i, jj=%i, k=%i, this->sizeA[%i]=%i\n", i, j, ji, jj, k, i, this->sizeA[i]);

//	if( this->sizeA[i]==0 || this->JA[i][ji] != j /*jj==-1*/){
	if( this->size==0 || (this->JA[k] != j && jj<ji)){
		if( jj<ji)
			k = ji;
		//fprintf( stderr, "new element at the beginning (%i, %i)\n", i, k);
		// realloc memory
		if( size == sizeAllocated){
			sizeAllocated += ROW_EXTENSION_SIZE;
			JA = (int*) realloc( JA, sizeAllocated * sizeof(int));
			assert( JA);
			A  = (T*) realloc( A, sizeAllocated * sizeof(T));
			assert( A);
		}
		// shift following elements
		for( int jjj=this->size; jjj>k; jjj--){
			this->JA[jjj] = this->JA[jjj-1];
			this->A[jjj] = this->A[jjj-1];
		}
		JA[k] = j;
		A[k] = value; // set first time
		size ++;

	}else{
		//if( ji>jj){
			//fprintf( stderr, "replace (%i, %i)\n", i, k);
			A[k] += value; // add to old value
		//}else{
		//	fprintf( stderr, "don't know\n");
		//}
	}
	return;
}

template <class T> void SparseVector<T>::setLast( int j, T value)
{
	if( value == 0.) return;

	// allocate memory
	if( size == sizeAllocated){
		sizeAllocated += ROW_EXTENSION_SIZE;
		JA = (int*)   realloc( JA, sizeAllocated * sizeof(int));
		assert( JA);
		A  = (T*) realloc( A, sizeAllocated * sizeof(T));
		assert( A);
	}

	// add element
	JA[size] = j;
	A[size] = value;
	size ++;

	return;

}


template <class T> T& SparseVector<T>::getNonZero( int j)
{	return A[j];}

template <class T> int SparseVector<T>::getNonZeroIndex( int j)
{	return JA[j];}

template <class T> T SparseVector<T>::get( int j)
{
	//householding();
	for( int jj=0; jj<size; jj++){
		if( JA[jj] == j){
			return A[jj];
		}else if( JA[jj] > j){
			return 0.;
		}
	}

	return 0.;
}

template <class T> void SparseVector<T>::reset()
{
	size = 0;
}

template <class T>
void SparseVector<T>::householding()
{
	if( value != 0){
		set( column, value);
		value = 0;
	}
}

template <class T>
T& SparseVector<T>::operator()( int index)
{
	householding();
	int ii=0;
	for( ; ii<size && JA[ii]<index; ii++) ;

	if( ii<size && JA[ii]==index)
		return A[ii];
	else{
		column= index;
		return value;
	}
}

#endif /* SPARSEVECTOR_IPP_ */
