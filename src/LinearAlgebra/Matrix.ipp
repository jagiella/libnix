/*
 * Mtrix.ipp
 *
 *  Created on: 27.04.2012
 *      Author: jagiella
 */

#ifndef MATRIX_IPP_
#define MATRIX_IPP_

#include <math.h>
#include <stdlib.h>


#include "Mathematix.hpp"


	template <class T> T ** Matrix<T>::newMatrix( int dimi, int dimj)
	{
		int i;
		T **A;

		A = (T**) calloc( sizeof(T*), dimi );
		if( A == NULL) memoryAllocationError( (char*)"A");

		for( i=0; i<dimi; i++){
			A[i] = (T*) calloc( sizeof(T), dimj);
			if( A[i] == NULL) memoryAllocationError( (char*)"A[i]");
		}

		return A;
	}
	/*****************************************************************************/

	template <class T> void Matrix<T>::deleteMatrix( T ** matrix, int dimi)
	{
		int i;

		for( i=0; i<dimi; i++)
			free( matrix[i]);
		free( matrix);
	}
	/*****************************************************************************/

	template <class T> T Matrix<T>::getDeterministe( T **A, int dim)
	{
		if(dim==1)
			return A[0][0];
		if(dim==2)
			return A[0][0]*A[1][1] - A[0][1]*A[1][0];
		if(dim==3)
			return A[0][0]*A[1][1]*A[2][2] + A[0][1]*A[1][2]*A[2][0] + A[0][2]*A[1][0]*A[2][1] - A[0][2]*A[1][1]*A[2][0] - A[0][1]*A[1][0]*A[2][2] - A[0][0]*A[1][2]*A[2][1];

		T **B = Matrix<T>::newMatrix( dim-1, dim-1);
		if( B == NULL) exit( 0);

		T det = 0.;
		for( int i=0; i<dim; i++)
			for( int j=0; j<dim; j++){
				for( int ii=0; ii<dim; ii++)
					for( int jj=0; jj<dim; jj++){
						B[(ii<=i?ii:ii-1)][(jj<=j?jj:jj-1)] = A[i][j];
					}
				det += pow(-1, i+j) * A[i][j] * getDeterministe( B, dim-1);
			}
		deleteDoubleMatrix( B, dim-1);

		return det;
	}
	/*****************************************************************************/

	template <class T> void Matrix<T>::matrixInversion( T **Ai, T **Ao, int n)
	{
		//double **A = newDoubleMatrix( n, n);//[n][n];
		T A[n][n];

		// Init Ao & A
		for( int i=0; i<n; i++)
			for( int j=0; j<n; j++){
				if(i==j)
					Ao[i][j] = 1;
				else
					Ao[i][j] = 0;
				A[i][j] = Ai[i][j];
			}
		//printMatrix( Ao, n, n, "Init Ao", "%10.5lf ");

		// Gauss-Jordan-Algorithm
		//forward
		//actual diagonal element
		for( int d=0; d<n; d++){
			// rescale first line
			for( int j=d+1; j<n; j++)
				A[d][j] /= A[d][d];
			for( int j=0; j<n; j++)
				Ao[d][j] /= A[d][d];
			A[d][d] = 1.;

			//
			for( int i=d+1; i<n; i++)
				if(A[i][d]!=0){
					// adapt folowing lines
					//double c = (A[i][d]/A[d][d]);
					for( int j=d+1; j<n; j++)
						//B[i][j]=B[i][j]-B[i][k]*B[k][j];
						A[i][j] = A[i][j] -A[d][j] *A[i][d];
					for( int j=0; j<n; j++)
						Ao[i][j]= Ao[i][j]-Ao[d][j]*A[i][d];
					A[i][d] = 0.; // NOT NECESSARY!
				}
		}
		//printMatrix( Ao, n, n, "forward Ao", "%10.5lf ");
		//printMatrix( A, n, n, "forward A", "%10.5lf ");

		// backward
		for( int d=n-1; d>0; d--){
			for( int i=0; i<d; i++)
				if(A[i][d]!=0){
					for( int j=0; j<n; j++)
						Ao[i][j]= Ao[i][j]-Ao[d][j]*A[i][d];
					//fprintf( stderr , "test\n");
					A[i][d] = 0.; // NOT NECESSARY!
				}
				//else
					//fprintf( stderr , "test2 \n");
		}
		//printMatrix( Ao, n, n, "backward Ao", "%10.5lf ");
		//printMatrix( A, n, n, "backward A", "%10.5lf ");
	}
	/*****************************************************************************/

	template <class T> void Matrix<T>::matrixTranspose( T **Ai, T **Ao, int n, int m)
	{
	//#pragma omp parallel for
		for( int i=0; i<n; i++)
			for( int j=0; j<m; j++)
				Ao[j][i] = Ai[i][j];

	}
	/*****************************************************************************/
	template <class T> void Matrix<T>::solveLinearSystem( T **A, T *b, T *x, int dim)
	{
		T **B = newMatrix( dim, dim);
		if( B == NULL){
			exit( 0);
		}
		solveLinearSystem( A, b, x, dim, B);
		deleteMatrix( B, dim);
	}
	/*****************************************************************************/
	// PIVOT GAUSSE //
	template <class T> void Matrix<T>::solveLinearSystem( T **A, T *b, T *x, int dim, T **B)
	{
		int i, // index of equation
		    j; // index of column
		int k;
		//double B[dim][dim];
		//float **B = newMatrix( dim, dim);

		// copy matrix
		for( i=0; i<dim; i++){
			for( j=0; j<dim; j++){
				B[i][j] = A[i][j];
		//		printf("%lf  ", B[i][j]);
			}
		//	printf("| %lf\n", b[i]);
		}

		// solving the linear system

		// forward reduction
		for( k=0; k<dim; k++){
			//printf("forward reduction: line %i/%i!\n", k+1, dim);
			if(B[k][k]==0){
				// find better line
	//			printf("looking for better line!\n");
				for( i=k+1; i<dim && B[i][k]==0; i++);

				// resort lines
				if(i<dim && B[i][k]!=0){
	//				printf("resort!\n");
					T temp;
					for( j=k; j<dim; j++){
						temp = B[i][j];
						B[i][j] = B[k][j];
						B[k][j] = temp;
					}
					temp = b[i];
					b[i] = b[k];
					b[k] = temp;
				}
			}
			if(B[k][k]!=0){
				// normalize first row
				for( j=k+1; j<dim; j++)
					B[k][j]=B[k][j]/B[k][k];
				b[k]=b[k]/B[k][k];
				B[k][k]=1.;

				// reduce following rows
				for( i=k+1; i<dim; i++){
					for( j=k+1; j<dim; j++)
						B[i][j]=B[i][j]-B[i][k]*B[k][j];
					b[i]=b[i]+b[k]*-B[i][k];
					B[i][k]=0;
				}
			}
		}

		/*printf("----------------------------\n");
		for( i=0; i<dim; i++){
			for( j=0; j<dim; j++){
				printf("%lf  ", B[i][j]);
			}
			printf("| %lf\n", b[i]);
		}//*/

		// backward reduction
		for( k=dim-1; k>=0; k--){
			if( B[k][k]!=0)
			for( i=0; i<k; i++){
				b[i]=b[i]+b[k]*-B[i][k];
				B[i][k]=0;
			}
		}

		/*printf("----------------------------\n");
		for( i=0; i<dim; i++){
			for( j=0; j<dim; j++){
				printf("%lf  ", B[i][j]);
			}
			printf("| %lf\n", b[i]);
		}//*/

		// copy solution
		for( i=0; i<dim; i++)
			x[i] = b[i];
			//x[i] = B[i][i];
	}
	/*****************************************************************************/

	template <class T> void Matrix<T>::matrixVectorProduct( T **A, T *b, T *x, int dim)
	{
	//#pragma omp parallel for
		for( int m=0; m<dim; m++){
			x[m] = 0.;
			for( int n=0; n<dim; n++){
				x[m] += A[m][n] * b[n];
			}
		}
	}
	/*****************************************************************************/


	template <class T> T Vector<T>::dotProduct( T *vectorA, T *vectorB, int dim)
	{
		T ret = 0.;

	#pragma omp parallel for reduction(+:ret)
		for( int i=0; i<dim; i++){
			ret += vectorA[i]*vectorB[i];
		}

		return ret;
	}
	/*****************************************************************************/


	template <class T> T Vector<T>::dotProduct3D( T *vectorA, T *vectorB)
	{
		T ret = 0.;

	//#pragma omp parallel for reduction(+:ret)
		for( int i=0; i<3; i++){
			ret += vectorA[i]*vectorB[i];
		}

		return ret;
	}
	/*****************************************************************************/


	template <class T> void Vector<T>::vectorScale3D( T *vector, T scalar, T *vectorxscalar)
	{
		for( int i=0; i<3; i++){
			vectorxscalar[i] = vector[i] * scalar;
		}
	}
	/*****************************************************************************/

	template <class T> void Vector<T>::crossProduct3D( T *vectorA, T *vectorB, T *vectorAxB)
	{
		//(a2b3 − a3b2, a3b1 − a1b3, a1b2 − a2b1)
		vectorAxB[0] = vectorA[1]*vectorB[2] - vectorA[2]*vectorB[1];
		vectorAxB[1] = vectorA[2]*vectorB[0] - vectorA[0]*vectorB[2];
		vectorAxB[2] = vectorA[0]*vectorB[1] - vectorA[1]*vectorB[0];
	}
	/*****************************************************************************/

	template <int Dimensions, class T> inline T triangleInterface( T *A, T *B, T *C)
	{
		T a=0;
		T b=0;
		T c=0;
		for( int i=0; i<Dimensions; i++){
			a += pow(B[i]-C[i],2.);
			b += pow(A[i]-C[i],2.);
			c += pow(A[i]-B[i],2.);
		}
		a=sqrt(a);
		b=sqrt(b);
		c=sqrt(c);

		T s=(a+b+c)/2.;
		return sqrt(s*(s-a)*(s-b)*(s-c));
	}

	template <class T> inline T tetrahedronHeight( T *A, T *B, T *C, T *D)
	{
		// HEIGHT
		T a[3] = {A[0]-B[0],A[1]-B[1],A[2]-B[2]};
		T b[3] = {C[0]-B[0],C[1]-B[1],C[2]-B[2]};
		T norm[3] ={
				a[1]*b[2] - a[2]*b[1],
				a[2]*b[0] - a[0]*b[2],
				a[0]*b[1] - a[1]*b[0]
		};
		T norm_abs=sqrt(pow(norm[0],2)+pow(norm[1],2)+pow(norm[2],2));
		norm[0]/=norm_abs;
		norm[1]/=norm_abs;
		norm[2]/=norm_abs;

		return norm[0]*(D[0] - A[0])
		     + norm[1]*(D[1] - A[1])
		     + norm[2]*(D[2] - A[2]);
	}

	template <class T> inline T tetrahedronVolume( T *A, T *B, T *C, T *D)
	{
		// INTERFACE ABC
		T interfaceABC = triangleInterface<3,T>( A, B, C);

		// HEIGHT
		T height = tetrahedronHeight<T>( A, B, C, D);

		return fabs(height)*interfaceABC/3.;
	}

#endif /* MATRIX_IPP_ */
