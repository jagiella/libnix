/*
 * Covariance.cpp
 *
 *  Created on: 07.11.2013
 *      Author: jagiella
 */


#include <gsl/gsl_math.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_statistics_double.h>


void gsl_covariance_matrix(double **A, int m, int n, double **covA) {
	for (int i = 0; i < m; i++) {
		for (int j = i; j < m; j++) {
			covA[i][j] = covA[j][i] = gsl_stats_covariance(A[i], 1, A[j], 1, n);
		}
	}
}

void gsl_covariance_matrix(double **A, int m, int n, double *covA) {
	for (int i = 0; i < m; i++) {
		for (int j = i; j < m; j++) {
			covA[i * m + j] = covA[j * m + i] = gsl_stats_covariance(A[i], 1,
					A[j], 1, n);
		}
	}
}

