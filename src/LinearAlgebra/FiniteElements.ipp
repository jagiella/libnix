/*
 * FinitsElements.tcc
 *
 *  Created on: 21.03.2012
 *      Author: jagiella
 */

//#include "../Geometry/Triangulation.hpp"

#include <stdio.h>

template <int Dimensions>
void assembleGlobalMassMatrix( Triangulation<Dimensions> *vd, SparseMatrix<double> *B)
{
	// Reset Mass Matrix
	for( int i=0; i<vd->countVertices; i++)
		B->resetRow(i);

	double **jacobi = newDoubleMatrix( Dimensions, Dimensions);


	// for each tet/triangle
	for( int t=0; t<vd->countTetrahedra; t++)
	if( !vd->SimplexContainsFramePoints( vd->tetrahedra[t])){
		double detJ;
		//int vv=0;

		//do{
		// Jacobian of Tet
		for( int p=0; p<Dimensions; p++)
			for( int v=0; v<Dimensions; v++)
				jacobi[p][v] = vd->tetrahedra[t]->vertices[v]->position[p] - vd->tetrahedra[t]->vertices[Dimensions]->position[p];
			//for( int v=0; v<Dimensions+1; v++)
			//	if( v!=vv)
			//	jacobi[p][(v<vv?v:v-1)] = vd->tetrahedra[t]->vertices[v]->position[p] - vd->tetrahedra[t]->vertices[vv]->position[p];
		//printMatrix( jacobi, Dimensions, Dimensions, "Jacobian J", "%10.5lf ");

		// Determinant
		detJ = getDeterministe( jacobi, Dimensions);
		//fprintf( stderr, "det(J) = %lf\n", detJ);
		//	vv++;
		//}while( detJ<0.);

		//matrixScale( mass, mass, detJ, 4, 4);
		//printMatrix( mass, 4, "det(J)*B", "%10.5lf ");
		for( int i=0; i<Dimensions+1; i++)
			for( int j=0; j<Dimensions+1; j++){
				//fprintf( stderr, "add (%i, %i)\n", vd->tetrahedra[t]->vertices[i]->index, vd->tetrahedra[t]->vertices[j]->index);
				switch( Dimensions){
				case 1:
					if( i==j)
						B->add(vd->tetrahedra[t]->vertices[i]->index, vd->tetrahedra[t]->vertices[j]->index, 1./3.  * fabs(detJ));
					else
						B->add(vd->tetrahedra[t]->vertices[i]->index, vd->tetrahedra[t]->vertices[j]->index, 1./6. * fabs(detJ));
					break;

				case 2:
					if( i==j)
						B->add(vd->tetrahedra[t]->vertices[i]->index, vd->tetrahedra[t]->vertices[j]->index, 1./12.  * fabs(detJ));
					else
						B->add(vd->tetrahedra[t]->vertices[i]->index, vd->tetrahedra[t]->vertices[j]->index, 1./24. * fabs(detJ));
					break;

				case 3:
					if( i==j)
						B->add(vd->tetrahedra[t]->vertices[i]->index, vd->tetrahedra[t]->vertices[j]->index, 1./60.  * fabs(detJ));
					else
						B->add(vd->tetrahedra[t]->vertices[i]->index, vd->tetrahedra[t]->vertices[j]->index, 1./120. * fabs(detJ));
					break;
				}
			}
	}
}

template <int Dimensions>
void assembleGlobalMassMatrix( Triangulation<Dimensions> *vd, SparseMatrix<double> *B, double *q)
{
	// Reset Mass Matrix
	for( int i=0; i<vd->countVoronoiCells; i++)
		B->resetRow(i);

	// Local Mass Matrix
	double M[4][4];

	// Jacobi Matrix
	double **jacobi = newDoubleMatrix( Dimensions, Dimensions);


	// for each tet/triangle
	for( int t=0; t<vd->countTetrahedra; t++)
	if( !vd->tetrahedronContainsFramePoints( vd->tetrahedra[t])){
		double detJ;

		// Jacobian of Tet
		getJacobi( vd->tetrahedra[t], jacobi);
		//printMatrix( jacobi, Dimensions, Dimensions, "Jacobian J", "%10.5lf ");

		// Determinant
		detJ = getDeterministe( jacobi, Dimensions);
		//fprintf( stderr, "det(J) = %lf\n", detJ);

		for( int i=0; i<Dimensions+1; i++)
			for( int j=0; j<Dimensions+1; j++)
				M[i][j]=0;

		for( int k=0; k<Dimensions+1; k++){
			//fprintf( stderr, "\"%s%i\" = [\n", "I",k);
			for( int i=0; i<Dimensions+1; i++){
				for( int j=0; j<Dimensions+1; j++){
					double Iijk;
					switch(Dimensions){
					case 3:
					if( i==j && j==k)
						Iijk = 1./120.;
					else if(i!=j && i!=k && j!=k)
						Iijk = 1./720.;
					else
						Iijk = 1./360.;
					break;

					case 2:
					if( i==j && j==k)
						Iijk = 1./20.;
					else if(i!=j && i!=k && j!=k)
						Iijk = 1./120.;
					else
						Iijk = 1./60.;
					break;

					case 1:
					if( i==j && j==k)
						Iijk = 1./4.;
					else
						Iijk = 1./12.;
					break;
					}
					M[i][j] += q[vd->tetrahedra[t]->vertices[k]->index]*detJ*Iijk;
					//fprintf( stderr, "1/%4.0lf  ", 1/Iijk);
				}
				//fprintf( stderr, "\n");
			}
			//fprintf( stderr, "]\n");
		}


		for( int i=0; i<Dimensions+1; i++)
			for( int j=0; j<Dimensions+1; j++)
				//fprintf( stderr, "add (%i, %i)\n", vd->tetrahedra[t]->vertices[i]->index, vd->tetrahedra[t]->vertices[j]->index);
					B->add(vd->tetrahedra[t]->vertices[i]->index, vd->tetrahedra[t]->vertices[j]->index, M[i][j]);
	}
}

template <int Dimensions>
void assembleGlobalStiffnessMatrix( Triangulation<Dimensions> *vd, SparseMatrix<double> *B)
{
	// Reset Stiffness Matrix
	for( int i=0; i<vd->countVertices; i++)
		B->resetRow(i);

	//VoronoiCell *points[Dimensions];

	double **jacobi = newDoubleMatrix( Dimensions, Dimensions);

	double **refGradient = newDoubleMatrix( Dimensions, Dimensions+1);
	for( int p=0; p<Dimensions; p++){
		for( int v=0; v<Dimensions; v++)
			if( v==p)
				refGradient[p][v+1] = 1.;
			else
				refGradient[p][v+1] = 0.;
		refGradient[p][0] = -1.;
	}
	//printMatrix( refGradient, Dimensions, Dimensions+1, "Reference Gradient", "%10.5lf ");

	double **refGradientTransposed = newDoubleMatrix( Dimensions+1, Dimensions);
	matrixTranspose( refGradient, refGradientTransposed, Dimensions, Dimensions+1);
	//printMatrix( refGradientTransposed, Dimensions+1, Dimensions, "Transposed Reference Gradient", "%10.5lf ");

	//jacobi[0][0] = 1; jacobi[0][1] = 2; //jacobi[0][2] = 0;
	//jacobi[1][0] = 2.5; jacobi[1][1] = 3; //jacobi[1][2] = 0;
	//jacobi[2][0] = 3; jacobi[2][1] = 4; jacobi[2][2] = 1;
	//printMatrix( jacobi, Dimensions, Dimensions, "Jacobi", "%10.5lf ");
	double **jacobiInverse = newDoubleMatrix( Dimensions, Dimensions);
	//matrixInversion( jacobi, jacobiInverse, Dimensions);
	//printMatrix( jacobiInverse, Dimensions, Dimensions, "Inversed Jacobi", "%10.5lf ");

	double **temp = newDoubleMatrix( Dimensions+1, Dimensions);
	//matrixProduct( refGradientTransposed, jacobiInverse, temp, Dimensions+1, Dimensions, Dimensions);
	//printMatrix( temp, Dimensions+1, Dimensions, "Grad^T * J^-1", "%10.5lf ");

	double **jacobiInverseTransposed = newDoubleMatrix( Dimensions, Dimensions);
	//matrixTranspose( jacobiInverse, jacobiInverseTransposed, Dimensions, Dimensions);

	double **temp2 = newDoubleMatrix( Dimensions+1, Dimensions);
	//matrixProduct( temp, jacobiInverseTransposed, temp2, Dimensions+1, Dimensions, Dimensions);
	//printMatrix( temp2, Dimensions+1, Dimensions, "Grad^T * J^-1 * (J^-1)^T", "%10.5lf ");

	double **temp3 = newDoubleMatrix( Dimensions+1, Dimensions+1);
	//matrixProduct( temp2, refGradient, temp3, Dimensions+1, Dimensions, Dimensions+1);
	//printMatrix( temp3, Dimensions+1, Dimensions+1, "Grad^T * J^-1 * (J^-1)^T * Grad", "%10.5lf ");

	//matrixScale( double **Ai, double **Ao, double c, Dimensions+1, Dimensions+1);

	//exit(0);




	// for each tet/triangle
	for( int t=0; t<vd->countTetrahedra; t++)
	if( !vd->SimplexContainsFramePoints( vd->tetrahedra[t])){
		double detJ;

		int vv=0;

		//do{
			/*if(vv){
				// Permutate
				VoronoiCell *point = vd->tetrahedra[t]->vertices[2];
				//vd->tetrahedra[t]->vertices[4]=vd->tetrahedra[t]->vertices[3];
				vd->tetrahedra[t]->vertices[2]=vd->tetrahedra[t]->vertices[1];
				//vd->tetrahedra[t]->vertices[1]=vd->tetrahedra[t]->vertices[0];
				vd->tetrahedra[t]->vertices[1]=point;
			//	fprintf( stderr, "%i %i %i\n", vd->tetrahedra[t]->vertices[0]->index, vd->tetrahedra[t]->vertices[1]->index, vd->tetrahedra[t]->vertices[2]->index);
			}*/

			// Jacobi Matrix
			getJacobi( vd->tetrahedra[t], jacobi);
			while( jacobi[0][0]==0){
				fprintf( stderr, "Jacobi Matrix has to be permutated!\n");
				// Permutate
				Vertex<Dimensions> *point = vd->tetrahedra[t]->vertices[Dimensions-vv];
				//vd->tetrahedra[t]->vertices[4]=vd->tetrahedra[t]->vertices[3];
				vd->tetrahedra[t]->vertices[Dimensions-vv]=vd->tetrahedra[t]->vertices[0];
				//vd->tetrahedra[t]->vertices[1]=vd->tetrahedra[t]->vertices[0];
				vd->tetrahedra[t]->vertices[0]=point;
				getJacobi( vd->tetrahedra[t], jacobi);
				vv++;
			}
			// Determinant
			detJ = getDeterministe( jacobi, Dimensions);
			//fprintf( stderr, "vv=%i: det(J) = %lf\n", vv, detJ);

			//fprintf( stderr, "%i: det(J) = %lf\n", vv, detJ);

			vv++;
			if(vv==4)
				exit(0);
		//}while( detJ<0.);

		//
		matrixInversion( jacobi, jacobiInverse, Dimensions);
		//}while( jacobiInverse);
		//printMatrix( jacobi, Dimensions, Dimensions, "Jacobi", "%10.5lf ");
		//printMatrix( jacobiInverse, Dimensions, Dimensions, "Inversed Jacobi", "%10.5lf ");

		//double **temp = newDoubleMatrix( Dimensions+1, Dimensions);
		matrixProduct( refGradientTransposed, jacobiInverse, temp, Dimensions+1, Dimensions, Dimensions);
		//printMatrix( refGradientTransposed, Dimensions+1, Dimensions, "Grad^T", "%10.5lf ");
		//printMatrix( temp, Dimensions+1, Dimensions, "Grad^T * J^-1", "%10.5lf ");

		//double **jacobiInverseTransposed = newDoubleMatrix( Dimensions, Dimensions);
		matrixTranspose( jacobiInverse, jacobiInverseTransposed, Dimensions, Dimensions);

		//double **temp2 = newDoubleMatrix( Dimensions+1, Dimensions);
		matrixProduct( temp, jacobiInverseTransposed, temp2, Dimensions+1, Dimensions, Dimensions);
		//printMatrix( temp2, Dimensions+1, Dimensions, "Grad^T * J^-1 * (J^-1)^T", "%10.5lf ");

		//double **temp3 = newDoubleMatrix( Dimensions+1, Dimensions+1);
		matrixProduct( temp2, refGradient, temp3, Dimensions+1, Dimensions, Dimensions+1);
		//printMatrix( temp3, Dimensions+1, Dimensions+1, "Grad^T * J^-1 * (J^-1)^T * Grad", "%10.5lf ");

		double refTetVolume;
		switch(Dimensions){
		case 1:
			refTetVolume = 1.;
		break;
		case 2:
			refTetVolume = 1./2.;
		break;
		case 3:
			refTetVolume = 1./6.;
		break;
		}
		matrixScale( temp3, temp3, fabs(detJ)*refTetVolume, Dimensions+1, Dimensions+1);
		//printMatrix( temp3, Dimensions+1, Dimensions+1, "Grad^T * J^-1 * (J^-1)^T * Grad * |Omega|*det(J)", "%10.5lf ");

		//matrixScale( mass, mass, detJ, 4, 4);
		//printMatrix( mass, 4, "det(J)*B", "%10.5lf ");
		for( int i=0; i<Dimensions+1; i++)
			for( int j=0; j<Dimensions+1; j++){
#if DEBUG > 0
					if( isnan(temp3[i][j])){
						fprintf( stderr, "tet: [%i %i %i %i]\n", vd->tetrahedra[t]->vertices[0]->index, vd->tetrahedra[t]->vertices[1]->index, vd->tetrahedra[t]->vertices[2]->index, vd->tetrahedra[t]->vertices[3]->index);
						fprintf( stderr, "  x: [%lf %lf %lf %lf]\n", vd->tetrahedra[t]->vertices[0]->position[0], vd->tetrahedra[t]->vertices[1]->position[0], vd->tetrahedra[t]->vertices[2]->position[0], vd->tetrahedra[t]->vertices[3]->position[0]);
						fprintf( stderr, "  y: [%lf %lf %lf %lf]\n", vd->tetrahedra[t]->vertices[0]->position[1], vd->tetrahedra[t]->vertices[1]->position[1], vd->tetrahedra[t]->vertices[2]->position[1], vd->tetrahedra[t]->vertices[3]->position[1]);
						fprintf( stderr, "  z: [%lf %lf %lf %lf]\n", vd->tetrahedra[t]->vertices[0]->position[2], vd->tetrahedra[t]->vertices[1]->position[2], vd->tetrahedra[t]->vertices[2]->position[2], vd->tetrahedra[t]->vertices[3]->position[2]);
						printMatrix( jacobi, Dimensions, Dimensions, "Jacobi", "%e ");
						fprintf( stderr, "det(J) = %lf\n", detJ);
						printMatrix( jacobiInverse, Dimensions, Dimensions, "Inversed Jacobi", "%10.5lf ");

						printMatrix( refGradientTransposed, Dimensions+1, Dimensions, "Grad^T", "%10.5lf ");
						printMatrix( temp, Dimensions+1, Dimensions, "Grad^T * J^-1", "%10.5lf ");

						printMatrix( temp2, Dimensions+1, Dimensions, "Grad^T * J^-1 * (J^-1)^T", "%10.5lf ");

						printMatrix( temp3, Dimensions+1, Dimensions+1, "Grad^T * J^-1 * (J^-1)^T * Grad * |Omega|*det(J)", "%10.5lf ");
						exit( 0);
					}
#endif
					B->add( vd->tetrahedra[t]->vertices[i]->index, vd->tetrahedra[t]->vertices[j]->index, temp3[i][j]);
			}
	}
	//printMatrix( B, vd->countVoronoiCells, vd->countVoronoiCells, "Global Stiffness", "%10.5lf ");
}

template <int Dimensions>
void getJacobi( Simplex<Dimensions> *tet, double **jacobi)
{
	switch(0){
	case 0:
	// [ p1-p0, p2-p0, p3-p0 ]
	for( int p=0; p<Dimensions; p++)
		for( int v=0; v<Dimensions; v++)
			jacobi[p][v] = tet->vertices[v+1]->position[p] - tet->vertices[0]->position[p];
	break;

	case 1:
	// [ p0-p3, p1-p3, p2-p3 ]
	for( int p=0; p<Dimensions; p++)
		for( int v=0; v<Dimensions; v++)
			jacobi[p][v] = tet->vertices[v]->position[p] - tet->vertices[Dimensions]->position[p];
	break;

	case 2:
	// [ p3-p0, p2-p0, p1-p0 ]
	for( int p=0; p<Dimensions; p++)
		for( int v=0; v<Dimensions; v++)
			jacobi[p][v] = tet->vertices[Dimensions-v]->position[p] - tet->vertices[0]->position[p];
	break;

	case 3:
	// [ p2-p3, p1-p3, p0-p3 ]
	for( int p=0; p<Dimensions; p++)
		for( int v=0; v<Dimensions; v++)
			jacobi[p][v] = tet->vertices[Dimensions-v-1]->position[p] - tet->vertices[Dimensions]->position[p];
	break;
	}
}

template <int Dimensions>
void assembleGlobalStiffnessMatrix( Triangulation<Dimensions> *vd, double **B)
{
	// Reset Stiffness Matrix
	for( int i=0; i<vd->countVertices; i++)
		for( int j=0; j<vd->countVertices; j++)
			B[i][j] = 0.;

//	VoronoiCell *points[Dimensions];

	double **jacobi = newDoubleMatrix( Dimensions, Dimensions);

	double **refGradient = newDoubleMatrix( Dimensions, Dimensions+1);
	for( int p=0; p<Dimensions; p++){
		for( int v=0; v<Dimensions; v++)
			if( v==p)
				refGradient[p][v+1] = 1.;
			else
				refGradient[p][v+1] = 0.;
		refGradient[p][0] = -1.;
	}
	//printMatrix( refGradient, Dimensions, Dimensions+1, "Reference Gradient", "%10.5lf ");

	double **refGradientTransposed = newDoubleMatrix( Dimensions+1, Dimensions);
	matrixTranspose( refGradient, refGradientTransposed, Dimensions, Dimensions+1);
	//printMatrix( refGradientTransposed, Dimensions+1, Dimensions, "Transposed Reference Gradient", "%10.5lf ");

	//jacobi[0][0] = 1; jacobi[0][1] = 2; //jacobi[0][2] = 0;
	//jacobi[1][0] = 2.5; jacobi[1][1] = 3; //jacobi[1][2] = 0;
	//jacobi[2][0] = 3; jacobi[2][1] = 4; jacobi[2][2] = 1;
	//printMatrix( jacobi, Dimensions, Dimensions, "Jacobi", "%10.5lf ");
	double **jacobiInverse = newDoubleMatrix( Dimensions, Dimensions);
	//matrixInversion( jacobi, jacobiInverse, Dimensions);
	//printMatrix( jacobiInverse, Dimensions, Dimensions, "Inversed Jacobi", "%10.5lf ");

	double **temp = newDoubleMatrix( Dimensions+1, Dimensions);
	//matrixProduct( refGradientTransposed, jacobiInverse, temp, Dimensions+1, Dimensions, Dimensions);
	//printMatrix( temp, Dimensions+1, Dimensions, "Grad^T * J^-1", "%10.5lf ");

	double **jacobiInverseTransposed = newDoubleMatrix( Dimensions, Dimensions);
	//matrixTranspose( jacobiInverse, jacobiInverseTransposed, Dimensions, Dimensions);

	double **temp2 = newDoubleMatrix( Dimensions+1, Dimensions);
	//matrixProduct( temp, jacobiInverseTransposed, temp2, Dimensions+1, Dimensions, Dimensions);
	//printMatrix( temp2, Dimensions+1, Dimensions, "Grad^T * J^-1 * (J^-1)^T", "%10.5lf ");

	double **temp3 = newDoubleMatrix( Dimensions+1, Dimensions+1);
	//matrixProduct( temp2, refGradient, temp3, Dimensions+1, Dimensions, Dimensions+1);
	//printMatrix( temp3, Dimensions+1, Dimensions+1, "Grad^T * J^-1 * (J^-1)^T * Grad", "%10.5lf ");

	//matrixScale( double **Ai, double **Ao, double c, Dimensions+1, Dimensions+1);

	//exit(0);




	// for each tet/triangle
	for( int t=0; t<vd->countTetrahedra; t++)
	if( !vd->SimplexContainsFramePoints( vd->tetrahedra[t])){
		double detJ;

		int vv=0;

		//do{
			if(vv){
				// Permutate
				Vertex<Dimensions> *point = vd->tetrahedra[t]->vertices[2];
				//vd->tetrahedra[t]->vertices[4]=vd->tetrahedra[t]->vertices[3];
				vd->tetrahedra[t]->vertices[2]=vd->tetrahedra[t]->vertices[1];
				//vd->tetrahedra[t]->vertices[1]=vd->tetrahedra[t]->vertices[0];
				vd->tetrahedra[t]->vertices[1]=point;
			//	fprintf( stderr, "%i %i %i\n", vd->tetrahedra[t]->vertices[0]->index, vd->tetrahedra[t]->vertices[1]->index, vd->tetrahedra[t]->vertices[2]->index);
			}

			// Jacobi Matrix
			getJacobi<Dimensions>( vd->tetrahedra[t], jacobi);

			// Determinant
			detJ = getDeterministe( jacobi, Dimensions);
			//fprintf( stderr, "vv=%i: det(J) = %lf\n", vv, detJ);

			//fprintf( stderr, "%i: det(J) = %lf\n", vv, detJ);

			vv++;
			if(vv==4)
				exit(0);
		//}while( detJ<0.);

		//
		matrixInversion( jacobi, jacobiInverse, Dimensions);
		//printMatrix( jacobi, Dimensions, Dimensions, "Jacobi", "%10.5lf ");
		//printMatrix( jacobiInverse, Dimensions, Dimensions, "Inversed Jacobi", "%10.5lf ");

		//double **temp = newDoubleMatrix( Dimensions+1, Dimensions);
		matrixProduct( refGradientTransposed, jacobiInverse, temp, Dimensions+1, Dimensions, Dimensions);
		//printMatrix( refGradientTransposed, Dimensions+1, Dimensions, "Grad^T", "%10.5lf ");
		//printMatrix( temp, Dimensions+1, Dimensions, "Grad^T * J^-1", "%10.5lf ");

		//double **jacobiInverseTransposed = newDoubleMatrix( Dimensions, Dimensions);
		matrixTranspose( jacobiInverse, jacobiInverseTransposed, Dimensions, Dimensions);

		//double **temp2 = newDoubleMatrix( Dimensions+1, Dimensions);
		matrixProduct( temp, jacobiInverseTransposed, temp2, Dimensions+1, Dimensions, Dimensions);
		//printMatrix( temp2, Dimensions+1, Dimensions, "Grad^T * J^-1 * (J^-1)^T", "%10.5lf ");

		//double **temp3 = newDoubleMatrix( Dimensions+1, Dimensions+1);
		matrixProduct( temp2, refGradient, temp3, Dimensions+1, Dimensions, Dimensions+1);
		//printMatrix( temp3, Dimensions+1, Dimensions+1, "Grad^T * J^-1 * (J^-1)^T * Grad", "%10.5lf ");

		double refTetVolume;
		switch(Dimensions){
		case 1:
			refTetVolume = 1.;
			break;
		case 2:
			 refTetVolume = 1./2.;
			break;
		case 3:
			refTetVolume = 1./6.;
			break;
		}
		matrixScale( temp3, temp3, fabs(detJ)*refTetVolume, Dimensions+1, Dimensions+1);
		//printMatrix( temp3, Dimensions+1, Dimensions+1, "Grad^T * J^-1 * (J^-1)^T * Grad * |Omega|*det(J)", "%10.5lf ");

		//matrixScale( mass, mass, detJ, 4, 4);
		//printMatrix( mass, 4, "det(J)*B", "%10.5lf ");
		for( int i=0; i<Dimensions+1; i++)
			for( int j=0; j<Dimensions+1; j++){
					B[vd->tetrahedra[t]->vertices[i]->index][vd->tetrahedra[t]->vertices[j]->index] += temp3[i][j];
			}
	}
	//printMatrix( B, vd->countVoronoiCells, vd->countVoronoiCells, "Global Stiffness", "%10.5lf ");
}

template <int Dimensions>
void assembleGlobalMassMatrix( Triangulation<Dimensions> *vd, double **B)
{
	// Reset Mass Matrix
	for( int i=0; i<vd->countVertices; i++)
		for( int j=0; j<vd->countVertices; j++)
			B[i][j] = 0.;

	double **jacobi = newDoubleMatrix( Dimensions, Dimensions);


	// for each tet/triangle
	for( int t=0; t<vd->countTetrahedra; t++)
	if( !vd->SimplexContainsFramePoints( vd->tetrahedra[t])){
		double detJ;
		//int vv=0;

		//do{
		// Jacobian of Tet
		for( int p=0; p<Dimensions; p++)
			for( int v=0; v<Dimensions; v++)
				jacobi[p][v] = vd->tetrahedra[t]->vertices[v]->position[p] - vd->tetrahedra[t]->vertices[Dimensions]->position[p];
			//for( int v=0; v<DIMENSIONS+1; v++)
			//	if( v!=vv)
			//	jacobi[p][(v<vv?v:v-1)] = vd->tetrahedra[t]->vertices[v]->position[p] - vd->tetrahedra[t]->vertices[vv]->position[p];
		//printMatrix( jacobi, DIMENSIONS, DIMENSIONS, "Jacobian J", "%10.5lf ");

		// Determinant
		detJ = getDeterministe( jacobi, Dimensions);
		//fprintf( stderr, "det(J) = %lf\n", detJ);
		//	vv++;
		//}while( detJ<0.);

		//matrixScale( mass, mass, detJ, 4, 4);
		//printMatrix( mass, 4, "det(J)*B", "%10.5lf ");
		for( int i=0; i<Dimensions+1; i++)
			for( int j=0; j<Dimensions+1; j++){
				switch(Dimensions){
				case 1:
					if( i==j)
						B[vd->tetrahedra[t]->vertices[i]->index][vd->tetrahedra[t]->vertices[j]->index] += 1./3.  * fabs(detJ);
					else
						B[vd->tetrahedra[t]->vertices[i]->index][vd->tetrahedra[t]->vertices[j]->index] += 1./6. * fabs(detJ);
					break;
				case 2:
					if( i==j)
						B[vd->tetrahedra[t]->vertices[i]->index][vd->tetrahedra[t]->vertices[j]->index] += 1./12.  * fabs(detJ);
					else
						B[vd->tetrahedra[t]->vertices[i]->index][vd->tetrahedra[t]->vertices[j]->index] += 1./24. * fabs(detJ);
					break;
				case 3:
					if( i==j)
						B[vd->tetrahedra[t]->vertices[i]->index][vd->tetrahedra[t]->vertices[j]->index] += 1./60.  * fabs(detJ);
					else
						B[vd->tetrahedra[t]->vertices[i]->index][vd->tetrahedra[t]->vertices[j]->index] += 1./120. * fabs(detJ);
					break;
				}
			}
	}
}
