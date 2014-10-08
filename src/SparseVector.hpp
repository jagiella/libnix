/*
 * SparseVector.hpp
 *
 *  Created on: 20.04.2012
 *      Author: jagiella
 */

#ifndef SPARSEVECTOR_HPP_
#define SPARSEVECTOR_HPP_

#include <ostream>

/** \addtogroup LinearAlgebra
 *  @{
 */


/** Brief description which ends at this dot. Details follow
 *  here.
 */

template <class T> class SparseVector {
private:
	T 	*A;   // values of members
	int *JA;  // column index of members
	int m;    // size of vector

	int size; // non-zero members
	int sizeAllocated;

	// only for operator usage
	T	value;
	int	column;

public:
	/** Initializes a sparse vector object of a certain size.
	 *  \param n is the size of the vector.
	 */
	SparseVector( int n);
	/** Destructs the sparse vector object. */
	~SparseVector();

	/** \return the number elements. */
	int length();
	/** \return the number non-zero elements. */
	int lengthNonZero();

	/** Updates the value of a vector member.
	 *  \param i is the position inside the vector.
	 *  \param value is the new value.
	 */
	void set( int i, T value);
	void setNonZero( int ii, T value);

	/** Sets the value of a vector member. Member is required to be the last non-zero element of the vector.
	 *  \param i is the position inside the vector.
	 *  \param value is the new value.
	 */
	void setLast( int j, T value);

	void add( int j, T value);

	/** Returns the value of vector member i.
	 *  \param i is the position inside the vector.
	 *  \return the value of the member.
	 */
	T get( int i);


	/** Returns the value of the i-th non-zero element.
	 *  \param i is the index in respect to all non-zero elements.
	 *  \return the value of the member.
	 */
	T& getNonZero( int i);


	/** Returns the the i-th non-zero element.
	 *  \param i is the index in respect to all non-zero elements.
	 *  \return the index of the member.
	 */
	int getNonZeroIndex( int i);


	/** Resets the sparse vector. */
	void reset();


	/** Grants access to member i.
	 * ATTENTION: If non-zero elements are manipulated,
	 * it needs a householding() call after the last
	 * manipulation to insure the matrix to be uptodate!)
	 *  \param i is the index of the requested member.
	 *  \return the value of the member.
	 */
	T& operator()( int i);

	/** All manipulations done with operator()() needs a
	 * housholding to insure the matrix to be uptodate. */
	void householding();

	friend std::ostream& operator<<(std::ostream& os, SparseVector<T>& A)
	{
		A.householding();
		int lasti=0;
		for( int ii=0; ii<A.size; ii++){
			int newi=A.JA[ii];
			for( int i=lasti; i<newi; i++)
				os << "0 ";
			os << A.A[ii] << " ";
			lasti = newi+1;
		}
		for( int i=lasti; i<A.m; i++)
			os << "0 ";
		return os;
	}

};

/** @}*/


#include "SparseVector.ipp"

#endif /* SPARSEVECTOR_HPP_ */
