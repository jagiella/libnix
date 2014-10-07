#ifndef MATRIX_IPP

#include <stdio.h>
#include <math.h>

// Definition of a singular matrix
#define SINGULAR 1e-16

template <class T>
T **newMatrix( int n, int m)
{
	T **A = (T**) malloc( n*sizeof(T*));
	for( int i=0; i<n; i++){
		A[i] = (T*) malloc( m*sizeof(T));
		for( int j=0; j<m; j++)
			A[i][j] = 0.;
	}
	return A;
}

template <class T>
void deleteMatrix( T** &A, int n)
{

	for( int i=0; i<n; i++)
		free(A[i]);
	free(A);
}

template <class T>
void printMatrix( T** &A, int n, int m)
{
	for( int i=0; i<n; i++){
		for( int j=0; j<m; j++){
			fprintf( stderr, "%10.2lf ", A[i][j]);
		}
		fprintf( stderr, "\n");
	}
}


template <class T>
void vectorMatrixProduct( T* &v, T** &A, T* &u, int m, int n)
{
	for( int j=0; j<n; j++){
		u[j] = 0;
		for( int i=0; i<m; i++)
			u[j] += v[i] * A[i][j];
	}
}
template <class T>
void matrixVectorProduct( T** &A, T* &v, T* &u, int m, int n)
{
	for( int i=0; i<m; i++){
		u[i] = 0;
		for( int j=0; j<n; j++)
			u[i] += A[i][j] * v[j];
	}
}
template <class T>
T dotProduct( T* &v1, T* &v2, int n)
{
	T dotProd = 0;
	for( int i=0; i<n; i++){
		dotProd += v1[i] * v2[i];
	}
	return dotProd;
}


#include <gsl/gsl_matrix.h>
#include <gsl/gsl_permutation.h>
#include <gsl/gsl_linalg.h>

template <class T>
double determinantLU( T** &A, int n) {

	gsl_matrix *invA = gsl_matrix_alloc( n, n);
	for( int i=0; i<n; i++)
		for( int j=0; j<n; j++)
			invA->data[i * invA->tda + j] = A[i][j];

 	double det;
 	int signum;
 	gsl_permutation *p = gsl_permutation_alloc( invA->size1);

	gsl_linalg_LU_decomp(invA , p , &signum);
 	det = gsl_linalg_LU_det(invA , signum);

 	gsl_permutation_free(p);
 	gsl_matrix_free(invA);


  	return det;
}

template <class T>
double determinantCholesky( T** &A, int n) {

	gsl_matrix *invA = gsl_matrix_alloc( n, n);
	for( int i=0; i<n; i++)
		for( int j=0; j<n; j++)
			invA->data[i * invA->tda + j] = A[i][j];

 	double det = 1;

	gsl_linalg_cholesky_decomp(invA);

	for( int i=0; i<n; i++)
		det = pow( invA->data[i * invA->tda + i], 2);

 	gsl_matrix_free(invA);

   
  	return det;
}

template <class T>
void inverseCholesky( T** &A, T** &invA, int n) {

	gsl_matrix *invCholesky = gsl_matrix_alloc( n, n);
	for( int i=0; i<n; i++)
		for( int j=0; j<n; j++){
			invCholesky->data[i * invCholesky->tda + j] = A[i][j];
		}

	gsl_linalg_cholesky_decomp ( invCholesky);
	gsl_linalg_cholesky_invert ( invCholesky);

	for( int i=0; i<n; i++)
		for( int j=0; j<n; j++){
			invA[i][j] = invCholesky->data[i * invCholesky->tda + j];
		}


	gsl_matrix_free( invCholesky);
}

template <class T>
int inverseLU( T** &A, T** &invA, int n) {
	
 	int signum;
 	gsl_permutation *p = gsl_permutation_alloc( n);
	gsl_matrix *LU    = gsl_matrix_alloc( n, n);
	gsl_matrix *invLU = gsl_matrix_alloc( n, n);
	for( int i=0; i<n; i++)
		for( int j=0; j<n; j++){
			LU->data[i * LU->tda + j] = A[i][j];
			invLU->data[i * LU->tda + j] = A[i][j];
		}

	int error = 0;
	if( (error=gsl_linalg_LU_decomp(LU, p, &signum)) ) return error;
	if( (error=gsl_linalg_LU_invert(LU, p, invLU)) ) return error;
	
	for( int i=0; i<n; i++)
		for( int j=0; j<n; j++){
			invA[i][j] = invLU->data[i * LU->tda + j];
		}

	gsl_permutation_free(p);
	gsl_matrix_free(    LU);
	gsl_matrix_free( invLU);

	return 0;
}

template <class T>
int determinantAndInverseLU( T** &A, T** &invA, T &det, int n) {

	bool success = true;
 	int signum;
 	gsl_permutation *p = gsl_permutation_alloc( n);
	gsl_matrix *LU    = gsl_matrix_alloc( n, n);
	gsl_matrix *invLU = gsl_matrix_alloc( n, n);
	for( int i=0; i<n; i++)
		for( int j=0; j<n; j++){
			LU->data[i * LU->tda + j] = A[i][j];
			invLU->data[i * LU->tda + j] = A[i][j];
		}

	// LU Decomposition
	gsl_linalg_LU_decomp(LU, p, &signum);

	// Determinant
	det = gsl_linalg_LU_det(LU , signum);

	if( det < SINGULAR)
		success = 1;
	else{

		// Inversion
		int error=gsl_linalg_LU_invert(LU, p, invLU);
		if(error /*&& error != GSL_EDOM*/)
			return error;

		for( int i=0; i<n; i++)
			for( int j=0; j<n; j++){
				invA[i][j] = invLU->data[i * LU->tda + j];
			}
	}

	gsl_permutation_free(p);
	gsl_matrix_free(    LU);
	gsl_matrix_free( invLU);

	return 0;
}



template <class T>
void decomposeGaussJordan( T** &A, T** &LR, int n) {
	// Gauss Elimination
	for( int i=0; i<n; n++){
	       // Bestimmen von R
	       for( int j=i; j<n; j++){
	           for( int k=0; k<i-1; k++){
	               A[i][j] -= A[i][k] * A[k][j];
	           }
	       }
	       // Bestimmen von L
	       for( int j=i+1; j<n; j++){
	           for( int k=0; k<i-1; k++){
	               A[j][i] -= A[j][k] * A[k][i];
	           }
	           A[j][i] /= A[i][i];
	       }
	}
}

template <class T>
void invertGaussJordan( T** &a, T** &ainv, int n) {
	int i, j;                    // Zeile, Spalte
	int s;                       // Elimininationsschritt
	int pzeile;                  // Pivotzeile
	int fehler = 0;              // Fehlerflag
	double f;                      // Multiplikationsfaktor
	const double Epsilon = 0.01;   // Genauigkeit
	double Maximum;                // Zeilenpivotisierung
	FILE *fout = stdout;
	int pivot = 1;

	// ergänze die Matrix a um eine Einheitsmatrix (rechts anhängen)
	/*for (i = 0; i < n; i++) {
		for (j = 0; j < n; j++) {
			a[i][n + j] = 0.0;
			if (i == j)
				a[i][n + j] = 1.0;
		}
	}*/
	for (i = 0; i < n; i++) {
		for (j = 0; j < n; j++) {
			ainv[i][j] = 0.0;
			if (i == j)
				ainv[i][j] = 1.0;
		}
	}
#if DEBUG
	MatOut (stdout, a, n, 2*n);
#endif

	// die einzelnen Eliminationsschritte
	s = 0;
	do {
		// Pivotisierung vermeidet unnötigen Abbruch bei einer Null in der Diagnonalen und
		// erhöht die Rechengenauigkeit
		Maximum = fabs(a[s][s]);
		if (pivot) {
			pzeile = s;
			for (i = s + 1; i < n; i++)
				if (fabs(a[i][s]) > Maximum) {
					Maximum = fabs(a[i][s]);
					pzeile = i;
				}
		}
		fehler = (Maximum < Epsilon);

		if (fehler)
			break;           // nicht lösbar

		if (pivot) {
			if (pzeile != s)  // falls erforderlich, Zeilen tauschen
					{
				double h;
				/*for (j = s; j < 2 * n; j++) {
					h = a[s][j];
					a[s][j] = a[pzeile][j];
					a[pzeile][j] = h;
				}*/
				for (j = s; j < n; j++) {
					h = a[s][j];
					a[s][j] = a[pzeile][j];
					a[pzeile][j] = h;
				}
				for (j = 0; j < n; j++) {
					h = ainv[s][j];
					ainv[s][j] = ainv[pzeile][j];
					ainv[pzeile][j] = h;
				}
			}
		}

		// Eliminationszeile durch Pivot-Koeffizienten f = a[s][s] dividieren
		f = a[s][s];
		/*for (j = s; j < 2 * n; j++)
			a[s][j] = a[s][j] / f;*/
		for (j = s; j < n; j++)
			a[s][j] = a[s][j] / f;
		for (j = 0; j < n; j++)
			ainv[s][j] = ainv[s][j] / f;

		// Elimination --> Nullen in Spalte s oberhalb und unterhalb der Diagonalen
		// durch Addition der mit f multiplizierten Zeile s zur jeweiligen Zeile i
		for (i = 0; i < n; i++) {
			if (i != s) {
				f = -a[i][s];                 // Multiplikationsfaktor
				/*for (j = s; j < 2 * n; j++)    // die einzelnen Spalten
					a[i][j] += f * a[s][j];*/       // Addition der Zeilen i, s
				for (j = s; j < n; j++)    // die einzelnen Spalten
					a[i][j] += f * a[s][j];       // Addition der Zeilen i, s
				for (j = 0; j < n; j++)    // die einzelnen Spalten
					ainv[i][j] += f * ainv[s][j];       // Addition der Zeilen i, s
			}
		}
#if DEBUG
		fprintf(stdout, "Nach %1i-tem Eliminationschritt:\n", s+1);
		MatOut (stdout, a, n, 2*n);
#endif
		s++;
	} while (s < n);

	if (fehler) {
		fprintf(fout, "Inverse: Matrix ist singulaer\n");
		//return 0;
	}
	// Die angehängte Einheitsmatrix Matrix hat sich jetzt in die inverse Matrix umgewandelt
	// Umkopieren auf die Zielmatrix
	/*for (i = 0; i < n; i++) {
		for (j = 0; j < n; j++) {
			ainv[i][j] = a[i][n + j];
		}
	}*/
}


template <class T>
void invertGaussJordan( T* &a, T* &ainv, int n) {
	int i, j;                    // Zeile, Spalte
	int s;                       // Elimininationsschritt
	int pzeile;                  // Pivotzeile
	int fehler = 0;              // Fehlerflag
	double f;                      // Multiplikationsfaktor
	const double Epsilon = 0.01;   // Genauigkeit
	double Maximum;                // Zeilenpivotisierung
	FILE *fout = stdout;
	int pivot = 1;

	// ergänze die Matrix a um eine Einheitsmatrix (rechts anhängen)
	/*for (i = 0; i < n; i++) {
		for (j = 0; j < n; j++) {
			a[i][n + j] = 0.0;
			if (i == j)
				a[i][n + j] = 1.0;
		}
	}*/
	for (i = 0; i < n; i++) {
		for (j = 0; j < n; j++) {
			ainv[i*n+j] = 0.0;
			if (i == j)
				ainv[i*n+j] = 1.0;
		}
	}
#if DEBUG
	MatOut (stdout, a, n, 2*n);
#endif

	// die einzelnen Eliminationsschritte
	s = 0;
	do {
		// Pivotisierung vermeidet unnötigen Abbruch bei einer Null in der Diagnonalen und
		// erhöht die Rechengenauigkeit
		Maximum = fabs(a[s*n+s]);
		if (pivot) {
			pzeile = s;
			for (i = s + 1; i < n; i++)
				if (fabs(a[i*n+s]) > Maximum) {
					Maximum = fabs(a[i*n+s]);
					pzeile = i;
				}
		}
		fehler = (Maximum < Epsilon);

		if (fehler)
			break;           // nicht lösbar

		if (pivot) {
			if (pzeile != s)  // falls erforderlich, Zeilen tauschen
					{
				double h;
				/*for (j = s; j < 2 * n; j++) {
					h = a[s][j];
					a[s][j] = a[pzeile][j];
					a[pzeile][j] = h;
				}*/
				for (j = s; j < n; j++) {
					h = a[s*n+j];
					a[s*n+j] = a[pzeile*n+j];
					a[pzeile*n+j] = h;
				}
				for (j = 0; j < n; j++) {
					h = ainv[s*n+j];
					ainv[s*n+j] = ainv[pzeile*n+j];
					ainv[pzeile*n+j] = h;
				}
			}
		}

		// Eliminationszeile durch Pivot-Koeffizienten f = a[s][s] dividieren
		f = a[s*n+s];
		/*for (j = s; j < 2 * n; j++)
			a[s][j] = a[s][j] / f;*/
		for (j = s; j < n; j++)
			a[s*n+j] = a[s*n+j] / f;
		for (j = 0; j < n; j++)
			ainv[s*n+j] = ainv[s*n+j] / f;

		// Elimination --> Nullen in Spalte s oberhalb und unterhalb der Diagonalen
		// durch Addition der mit f multiplizierten Zeile s zur jeweiligen Zeile i
		for (i = 0; i < n; i++) {
			if (i != s) {
				f = -a[i*n+s];                 // Multiplikationsfaktor
				/*for (j = s; j < 2 * n; j++)    // die einzelnen Spalten
					a[i][j] += f * a[s][j];*/       // Addition der Zeilen i, s
				for (j = s; j < n; j++)    // die einzelnen Spalten
					a[i*n+j] += f * a[s*n+j];       // Addition der Zeilen i, s
				for (j = 0; j < n; j++)    // die einzelnen Spalten
					ainv[i*n+j] += f * ainv[s*n+j];       // Addition der Zeilen i, s
			}
		}
#if DEBUG
		fprintf(stdout, "Nach %1i-tem Eliminationschritt:\n", s+1);
		MatOut (stdout, a, n, 2*n);
#endif
		s++;
	} while (s < n);

	if (fehler) {
		fprintf(fout, "Inverse: Matrix ist singulaer\n");
		//return 0;
	}
	// Die angehängte Einheitsmatrix Matrix hat sich jetzt in die inverse Matrix umgewandelt
	// Umkopieren auf die Zielmatrix
	/*for (i = 0; i < n; i++) {
		for (j = 0; j < n; j++) {
			ainv[i][j] = a[i][n + j];
		}
	}*/
}

#define MATRIX_IPP
#endif
