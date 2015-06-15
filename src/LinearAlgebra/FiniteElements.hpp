/*
 * FiniteElements.h
 *
 *  Created on: Jun 4, 2010
 *      Author: jagiella
 */

#ifndef FINITEELEMENTS_H_
#define FINITEELEMENTS_H_

#include "SparseMatrix.hpp"
#include "../Geometry/Triangulation.hpp"

template <int Dimensions>
void assembleGlobalStiffnessMatrix( Triangulation<Dimensions> *vd, SparseMatrix<double> *B);

template <int Dimensions>
void assembleGlobalMassMatrix( Triangulation<Dimensions>  *vd, SparseMatrix<double> *B);

template <int Dimensions>
void assembleGlobalMassMatrix( Triangulation<Dimensions>  *vd, SparseMatrix<double> *B, double *q);

template <int Dimensions>
void assembleGlobalMassMatrix( Triangulation<Dimensions>  *vd, int mOrder, SparseMatrix<double> *B, double *q, int qOrder);

template <int Dimensions>
void getJacobi( Simplex<Dimensions> *tet, double **jacobi);


#include "FiniteElements.ipp"

#endif /* FINITEELEMENTS_H_ */
