/*
 * testMatrix.cpp
 *
 *  Created on: 16.04.2014
 *      Author: jagiella
 */

#include "SparseMatrix.hpp"


int main()
{
	SparseMatrix<int> *A = new SparseMatrix<int>(5,5);

	for( int i=0; i<5; i++)
		A->add(i,i, 1);
	A->printMatrix( "A", "  %3i");


	for( int i=0; i<4; i++)
		A->add(i,4-i, -3);
	A->printMatrix( "A", "  %3i");


	for( int i=0; i<5; i++)
		A->add(i,i, 1);
	A->printMatrix( "A", "  %3i");

	return 0;
}
