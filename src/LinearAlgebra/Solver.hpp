/*
 * Solver.h
 *
 *  Created on: 21.08.2010
 *      Author: jagiella
 */

#ifndef SOLVER_H_
#define SOLVER_H_

#include "SparseMatrix.hpp"

/** \addtogroup LinearAlgebra
 *  @{
 */

template <class T> class Solver
{
private:
	// Vectors Allocated
	int		size;
	T **	vectors;
	int 	vectorCount;

	// Solver Type & Convergence Behavior
	char 	type;
	int  	maxit;
	T	 	maxerr;

	// Solver-Specific Parameters
	T		omega;	// SOR
	int		l;		// BiCGSTAB(l)

	// Exact Solution
	bool 	useExactSolution;
	T*		exactSolution;

	void calculate_x( T *x, T omega, T *w, T alpha, T *p, int N);
	void calculate_r( T *r, T omega, T *v, T alpha, T *s, int N);
	void calculate_p( T *x, T *y, T alpha, T *a, T beta, T *b, int N);
	void calculate_p1( T *x, T *y, T alpha, T *a, int N);
	void calculateScaledSum( T *o, T s1, T *v1, T s2, T *v2, T s3, T *v3, int N);
	void calculateScaledSum( T *o, T s1, T *v1, T s2, T *v2, int N);

	int SolveConjugateGradient( SparseMatrix<T> *A, T *b, T *x);
	int SolverPreconditionedConjugateGradient( SparseMatrix<T> *A, SparseMatrix<T> *M, T *b, T *x);
	int SolveBiCGSTAB(     SparseMatrix<T> *A, T *b, T *x);
	int SolveBiCGSTABl(     SparseMatrix<T> *A, T *b, T *x);
	int SolveBiCGSTAB2( SparseMatrix<T> *A, T *b, T *x);
	int SolveJacobiMethod( SparseMatrix<T> *A, T *b, T *x);
	int SolveGaussSeidelMethod( SparseMatrix<T> *A, T *b, T *x);
	int SolveSuccessiveOverRelaxationMethod( SparseMatrix<T> *A, T *b, T *x, T omega);
	int SolveExplicit( SparseMatrix<T> *A, T *b, T *x);
	int SolveLowerTriangular( SparseMatrix<T> *A, T *b, T *x);
	int SolveUpperTriangular( SparseMatrix<T> *A, T *b, T *x);
	int SolveLU( SparseMatrix<T> *LU, T *b, T *y, T *x);

public:
	/** All types the Solver class can use to iteratively solve a system. */
	typedef  const char solverType ;

	/** \brief Conjugate Gradient
	 *
	 * Hestenes, Magnus R.; Stiefel, Eduard (December 1952),
	 * "Methods of Conjugate Gradients for Solving Linear Systems",
	 * Journal of Research of the National Bureau of Standards 49 (6). */
	static solverType CG 			= 1;

	/** \brief Bi-Conjugate Gradient-Stabilized
	 *
	 * H.A. Van der Vorst,
	 * "Bi-CGSTAB: a fast and smoothly converging variant of Bi-CG for the solution of nonsymmetric linear systems",
	 * SIAM J. Sci. Statist. Comput., 13 (1992), pp. 631-644 */
	static solverType BiCGSTAB		= 2;

	/** \brief Jacobi Method */
	static solverType Jacobi		= 3;

	/** \brief Gauss-Seidel
	 *
	 * L.  Seidel,
	 * "Ueber  ein  Verfahren  die  Gleichungen,  auf  welche  die  Methode  der  kleinsten
	 * Quadrate  fuehrt,  sowie  lineare  Gleichungen  ueberhaupt  durch  successive  Annaeherung  aufzuloesen",
	 * Abhandlungen  der  Bayerischen  Akademie,  vol. 11, Dritte  Abteilung  (1873),  pp.  81-108.*/
	static solverType GaussSeidel	= 4;

	/** \brief Successive Over-Relaxation
	 *
	 * S.  Frankel,
	 * "Convergence  rates  of  iterative  treatments  of  partial  differential equations",
	 * Mathematical  Tables  and  Other  Aids  to  Computation,  vol.  4  (1950), pp.  65-75.
	 *
	 * David Young,
	 * "Iterative Methods for Solving Partial Difference Equations of Elliptic Type",
	 * Transactions of the American Mathematical Society Vol. 76, No. 1 (Jan., 1954), pp. 92-111 */
	static solverType SOR			= 5;

	/** \brief Incomplete LU-Factorization */
	static solverType ILU			= 6;

	/** \brief Bi-Conjugate Gradient-Stabilized (l)
	 *
	 * The algorithm generalizes the Bi-CGSTAB algorithm and overcomes some shortcomings of BiCGStab2.
	 * In some sense, the new algorithm combines GMRES(l) and Bi-CG. The parameter default value of l=2
	 * can be modified with setBiCGSTABdegree().
	 *
	 * It was published in:
	 * G.L.G. Sleijpen, D.R. Fokkema.
	 * "BICGSTAB(L) for linear equations involving unsymmetric matrices with complex spectrum".
	 * Electron. Trans. Numer. Anal., 1 (1993), pp. 11-32 */
	static solverType BiCGSTABl		= 7;

	// Solver Type (Preconditioner)
	static const char Diagonal			= 6;
	static const char Tridiagonal		= 7;
	static const char LowerTriangular	= 8;
	static const char UpperTriangular	= 9;

	/** The constructor
	 * \param matrix_size is the number of degrees of liberty.
	 * \param type is the solver type. (Default: BiCGSTAB) */
	Solver( int matrix_size, solverType type=BiCGSTAB, T error=1e-8, int iterations=1000);
	/** The destructor */
	~Solver();

	/** Solves the given system x=b/A
	 * \param A is the matrix.
	 * \param b is the vector.
	 * \param x is the vector with the solution of the system.
	 * \return the number of iterations needed to solve the system. */
	int solve( SparseMatrix<T> *A, T *b, T *x);

	/** modifies the degree of the BiCGSTABl method
	 *
	 * \param l has to be even and larger than 0
	 * \sa BiCGSTABl*/
	void setDegree( int l);
	void setIterations( int iterations);
	void setError( T error);

	/** modifies the relaxation parameter of SOR method
	 * \param omega should be a value between 0 and 2 (convergence proved).
	 * Values of omega>1 are used to speedup convergence of a slow-converging process,
	 * while values of  are often used to help establish convergence of a diverging iterative process.
	 * \sa SOR*/
	void setRelaxationParameter( T omega);

	void PreconditionJacobi( SparseMatrix<T> *A, T *b);
	void PreconditionSOR( SparseMatrix<T> *A, T *b, T omega);
	void SuccessiveOverRelaxationPreconditioner( SparseMatrix<T> *A, SparseMatrix<T> *M);
	void IncompleteLUfactorization( SparseMatrix<T> *A);
	void IncompleteLUfactorization( SparseMatrix<T> *A, SparseMatrix<T> *M);
};

/** @}*/

#include "Solver.ipp"


#endif /* SOLVER_H_ */
