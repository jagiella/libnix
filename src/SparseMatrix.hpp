#ifndef __SPARSE_MATRIX_H
#define __SPARSE_MATRIX_H

//#include "VoronoiDiagramExtended.h"

//#define DEBUG 0

//#define ROW_EXTENSION_SIZE 10

//#define COMPARE_WITH_EXACT_SOLUTION 0

//typedef struct _SparseVector SparseVector;

#include <iostream>
#include <ostream>
#include <string>

#include "SparseVector.hpp"

/** \addtogroup LinearAlgebra
 *  @{
 */

// voronoi cell type
template <class T> class SparseMatrix {
	

private:
	SparseVector<T> **A; // values
	int n, m; // dimensions
	
public:
	/** Initializes a sparse matrix object of a certain size.
	 *  \param n is the number of rows.
	 *  \param n is the number of columns.
	 */
	SparseMatrix( int n, int m);
	~SparseMatrix();
	
	/** \brief returns the number of columns
	 * \return the number of columns. */
	int columns();
	/** \brief returns the number of non-zero elements in row i
	 * \return the number of columns. */
	int columnsNonZero( int i);
	/** \brief returns the number of rows
	 * \return the number of rows. */
	int rows();

	/** \brief returns the element in row i, column j
	 * \param i is the row index.
	 * \param j is the column index.
	 * \return the value of member (i,j). */
	T	get( int i, int j);
	/** \brief returns a reference (value can be modified) to the element in row i, column j
	 * \param i is the row index.
	 * \param j is the column index.
	 * \return the value of element (i,j). */
	T&	getNonZero( int i, int j);
	/** \brief returns the row index of the j-th non-zero element in row i
	 * \param i is the row index.
	 * \param j is the column index.
	 * \return the value of element (i,j). */
	int	getNonZeroIndex( int i, int sparsej);

	/** Sets (or updates) value of member (i,j).
	 * \param value is the value to set
	 * \param i is the row index.
	 * \param j is the column index.
	 * \sa setLast(), operator()()*/
	void set( int i, int j, T value);
	void setNonZero( int i, int jj, T value);

	/** Sets value of member (i,j). Is faster than set(), but (i,j) is required to be the last non-zero element of row i!
	 * \param value is the value to set
	 * \param i is the row index.
	 * \param j is the column index.
	 * \sa set()*/
	void setLast( int i, int j, T value);

	/** Adds up value to member (i,j).
	 * \param value is the summand to add
	 * \param i is the row index.
	 * \param j is the column index.
	 * \sa setLast()*/
	void add( int i, int j, T value);

	/** Sets (or updates) value of member (i,j).
	 * \param value is the value to set
	 * \param i is the row index.
	 * \param j is the column index.
	 * \sa set(), setLast()*/
	void operator () ( int i, int j, T value);

	/** Grants access to member (i,j) for write and read.
	 * \param i is the row index.
	 * \param j is the column index.
	 * \sa get(), set()*/
	T& operator () ( int i, int j);

	/** Creates an exact copy of another sparse matrix object
	 * \param A is the sparse matrix to copy.
	 * \sa copy() */
	void operator = (SparseMatrix<T> *A);

	/**
	 * \param i is the row index.
	 * \param j is the column index.
	 * \return a copy of the sparse matrix object
	 * \sa operator=() */
	SparseMatrix<T>* copy();

	/** Resets the i-th row.
	 * \param i is the row index. */
	void resetRow( int i);

	/** Counts all non-zero elements.
	 * \return number of non-zero elements */
	int countNonZero(){
		int count = 0;
		for( int i=0; i<this->m; i++)
			count += columnsNonZero( i);
		return count;
	}


	void printMatrix( const char* name, const char* format);
	friend std::ostream& operator<<(std::ostream& os, SparseMatrix<T>& A)
	{
	for( int i=0; i<A.m; i++)
		os << (*A.A[i]) << std::endl;
	return os;
	}

	// static methodes
	static void MatrixVectorProduct( SparseMatrix<T> *A, T *b, T *x);
	static void VectorMatrixProduct( T *b, SparseMatrix<T> *A, T *x);
};

/** @}*/

#include "SparseMatrix.ipp"

#endif
