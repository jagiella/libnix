/*
 * Matrix.hpp
 *
 *  Created on: 27.04.2012
 *      Author: jagiella
 */

#ifndef MATRIX_HPP_
#define MATRIX_HPP_

/** \addtogroup LinearAlgebra
 *  @{
 */


/** Matrix operations*/
template <class T>
class Matrix
{
private:
	Matrix();
public:
	static T** newMatrix( int dimi, int dimj);
	static void deleteMatrix( T ** matrix, int dimi);
	static T getDeterministe( T **A, int dim);
	static void matrixInversion( T **Ai, T **Ao, int n);
	static void matrixTranspose( T **Ai, T **Ao, int n, int m);
	static void solveLinearSystem( T **A, T *b, T *x, int dim);
	static void solveLinearSystem( T **A, T *b, T *x, int dim, T **B);
	static void matrixVectorProduct( T **A, T *b, T *x, int dim);
};

/** Vector operations*/
template <class T>
class Vector
{
private:
	Vector();
public:
	static T dotProduct( T *vectorA, T *vectorB, int dim);
	static T dotProduct3D( T *vectorA, T *vectorB);
	static void crossProduct3D( T *vectorA, T *vectorB, T *vectorAxB);
	static void vectorScale3D( T *vector, T scalar, T *vectorXscalar);
};

/** @}*/

#include "Matrix.ipp"

#endif /* MATRIX_HPP_ */
