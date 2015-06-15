#include "../LinearAlgebra.hpp"

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

#define ABS( x) (x<0. ? -x : x)

float *b;//[20];
float *x;//[20];
double *exactSolution;//[20];



void setGramianMatrix( double **M, double *x, int N)
{
	// Integrate "hat functions" over whole domain
	for( int i=0; i<N; i++)
		for( int j=0; j<N; j++)
		{
			// border
			if( i==0 || i==N-1){
				M[i][j] = 0.;
				if(i==j)
					M[i][j] = 0.; // border value

			}else

			// i and j are equal?
			if( i==j){
				//M[i][j] = 2./3.;
				M[i][j] = 1./3.*( (i>0 ? x[i]-x[i-1] : 0.) + (i+1<N ? x[i+1]-x[i] : 0.) );
			}else
			
			// i and j are neighbors?
			if( i==j+1 || i==j-1){//abs(i-j) == 1){
				M[i][j] = 1./6. * (x[i]-x[j]) * (x[i]-x[j]);//4./3.;
			}else
			
				M[i][j] = 0.;
		}
	
}

/*void setGramianMatrix2D( double **M, VoronoiDiagram *vd, int N)
{
	// Reset Matrix
	for( int i=0; i<N; i++)
		for( int j=0; j<N; j++)
			M[i][j] = 0.;

	// Integrate "hat functions" over whole domain
	for( int t=0; t<vd->countTetrahedra; t++){
		for( int v=0; v<NR_TETRAHEDRON_POINTS; v++){
			for( int vv=v; vv<NR_TETRAHEDRON_POINTS; vv++){
				vd->tetrahedra[t];

			}
		}
	}


	for( int i=0; i<N; i++)
		for( int j=0; j<N; j++)
		{
			// border
			if( i==0 || i==N-1){
				M[i][j] = 0.;
				if(i==j)
					M[i][j] = 0.; // border value

			}else

			// i and j are equal?
			if( i==j){
				//M[i][j] = 2./3.;
				M[i][j] = 1./3.*( (i>0 ? x[i]-x[i-1] : 0.) + (i+1<N ? x[i+1]-x[i] : 0.) );
			}else

			// i and j are neighbors?
			if( i==j+1 || i==j-1){//abs(i-j) == 1){
				M[i][j] = 1./6. * (x[i]-x[j]) * (x[i]-x[j]);//4./3.;
			}else

				M[i][j] = 0.;
		}

}*/

void setDerivedGramianMatrix( double **M, double *x, int N)
{
	// Integrate "hat functions" over whole domain
	for( int i=0; i<N; i++)
		for( int j=0; j<N; j++)
		{
			// border
			if( i==0 || i==N-1){
				M[i][j] = 0.;
				if(i==j)
					M[i][j] = 1.; // border condition (DON'T TOUCHE!)
			}else

			// i and j are equal?
			if( i==j){
				//M[i][j] = 2.;
				M[i][j] = (i>0 ? 1./(x[i]-x[i-1]) : 0.) + (i+1<N ? 1./(x[i+1]-x[i]) : 0.);
				//M[i][j] = 0.;
				//if(i>0)   M[i][j] += 1./(x[i]-x[i-1]);
				//if(i+1<N) M[i][j] += 1./(x[i+1]-x[i]);

			}else
			
			// i and j are neighbors?
			if( i==j+1 || i==j-1){//abs(i-j) == 1){
				//M[i][j] = -1.;
				M[i][j] = -1./(i<j ? (x[j]-x[i]) : (x[i]-x[j]));
			}else
			
				M[i][j] = 0.;
		}
	
}

void main1D( int argc, char** argv)
{
	// INIT DOMAIN
	int N = (argc>1 ? atoi(argv[1]) : 11);
	double pos[N];
	double u[N];
	double f[N];
	double **M = newDoubleMatrix( N, N),
	       **L = newDoubleMatrix( N, N);

	// init values x
	for( int i=0; i<N; i++)
		pos[i] = 1.*(double)i/(double)(N-1);
	//printVector( pos, N);

	// init values u
	for( int i=0; i<N; i++)
		u[i] = 0.;

	// set Gramian Matrix M: v*v
	setGramianMatrix( M, pos, N);
	//printMatrix( M, N);

	// set Load Matrix L: v'*v'
	setDerivedGramianMatrix( L, pos, N);
	//printMatrix( L, N);

	// set vector f
	for( int i=0; i<N; i++){
		/*if( i==N*1/3)
			f[i] = 1.;
		else
		if( i==N*2/3)
			f[i] = 0.5;*/
		//if( i==N/2)
		//	f[i] = 1./N;
		//else
			f[i] = 1.;
	}
	//printVector( f, N);

	// multiply M*f
	double temp[N];
	matrixVectorProduct( M, f, temp, N);
	//printVector( temp, N);

	// solve system
	solveLinearSystem( L, temp, u, N);
	//printVector( u, N);

	for( int i=0; i<N; i++){
		printf( "%lf %lf\n", pos[i], u[i]);
	}
}







double getSqDist( double *p0, double *p1){
#if DIMENSIONS == 1
	return pow(p0[0]-p1[0], 2) + pow(p0[1]-p1[1], 2);
#elif DIMENSIONS == 2
	return pow(p0[0]-p1[0], 2) + pow(p0[1]-p1[1], 2) + pow(p0[2]-p1[2], 2);
#elif DIMENSIONS == 3
	return pow(p0[0]-p1[0], 2) + pow(p0[1]-p1[1], 2) + pow(p0[2]-p1[2], 2) + pow(p0[3]-p1[3], 2);
#endif
}


double getHeightNEW( double *p0, double *p1, double *p2, double *p3){
	double vector1[3], vector2[3], normal[3];

	vector1[0] = p2[0] - p1[0];
	vector1[1] = p2[1] - p1[1];
	vector1[2] = p2[2] - p1[2];

	vector2[0] = p3[0] - p1[0];
	vector2[1] = p3[1] - p1[1];
	vector2[2] = p3[2] - p1[2];

	// normal vector of plane p1-p2-p3
	normal[0] = vector1[1]*vector2[2] - vector1[2]*vector2[1];
	normal[1] = vector1[2]*vector2[0] - vector1[0]*vector2[2];
	normal[2] = vector1[0]*vector2[1] - vector1[1]*vector2[0];

	// normalize normal vector
	double length = sqrt( pow(normal[0],2.) + pow(normal[1],2.) + pow(normal[2],2.));
	if( length==0.){
		//std::cerr << "WARNING: Couldn't determine normal to plane:" << std::endl << p1 << std::endl << p2 << std::endl << p3 << std::endl;
		return 0.;
	}
	normal[0] /= length;
	normal[1] /= length;
	normal[2] /= length;

	// distance: D = n*(v-pi), i=1,2 or 3
	return fabs( normal[0]*(p0[0]-p1[0]) + normal[1]*(p0[1]-p1[1]) + normal[2]*(p0[2]-p1[2]));
//	return fabs(a*v->point().x() + b*v->point().y() + c*v->point().z() + d) /
//			sqrt( a*a + b*b + c*c);
}

double getVolumeIntegral( double teti[4][3], int i, int j)
{
	//printVector( teti[0], 3);
	//printVector( teti[1], 3);
	//printVector( teti[2], 3);
	//printVector( teti[3], 3);

	int A, B, C=i, D=j;

	FILE *fp;

		// copy & translate
		double tet[4][3];

		// index for A
		A=0;
		for(; A==i || A==j; A++) ;

		// index for B
		B=A+1;
		for(; B==i || B==j; B++) ;

		// index for C
		if( i==j){
			C=B+1;
			for(; C==i; C++) ;
		}

		//printVector( teti[A], 3);
		//printVector( teti[B], 3);
		//printVector( teti[C], 3);
		//printVector( teti[D], 3);

		// Translate
		for( int p=0; p<3; p++){
			tet[D][p] = teti[D][p] - teti[A][p];
			tet[C][p] = teti[C][p] - teti[A][p];
			tet[B][p] = teti[B][p] - teti[A][p];
			tet[A][p] = 0.;
		}

		//printVector( tet[A], 3);
		//printVector( tet[B], 3);
		//printVector( tet[C], 3);
		//printVector( tet[D], 3);

		// Rotate around xy-axis
		double alpha = (tet[D][1]==tet[C][1] ?
				0. :
				atan((tet[D][1]-tet[C][1])/(tet[C][2]-tet[D][2])));
		double beta = atan(
				(tet[D][0]-tet[C][0])
				/( sin(alpha)*(tet[C][1]-tet[D][1]) + cos(alpha)*(tet[D][2]-tet[C][2]) )
			);
		//fprintf( stderr, "alpha=%lf (%lf), beta=%lf\n", alpha, atan(0), beta);
//		for( int v=0; v<4; v++)
//			rotateZYX( tet[v][0], tet[v][1], tet[v][2], tet[v][0], tet[v][1], tet[v][2], alpha, beta, 0.);
		rotateZYX( tet[B][0], tet[B][1], tet[B][2], tet[B][0], tet[B][1], tet[B][2], alpha, beta, 0.);
		rotateZYX( tet[C][0], tet[C][1], tet[C][2], tet[C][0], tet[C][1], tet[C][2], alpha, beta, 0.);
		rotateZYX( tet[D][0], tet[D][1], tet[D][2], tet[D][0], tet[D][1], tet[D][2], alpha, beta, 0.);
		fp = fopen( "tetRotationYX.dat", "w+");
		for( int f=0; f<4; f++) for( int v=0; v<4; v++)	fprintf( fp, "%lf %lf %lf \n", tet[((v==3?0:v)+f)%4][0], tet[((v==3?0:v)+f)%4][1], tet[((v==3?0:v)+f)%4][2]);
		fclose(fp);

		// Rotate around z-axis
		double delta = atan((tet[B][0]-tet[C][0])/(tet[B][1]-tet[C][1]));
		//fprintf( stderr, "delta=%lf\n", delta);
		//for( int v=0; v<4; v++)
		//	rotateZ( tet[v][0], tet[v][1], tet[v][2], tet[v][0], tet[v][1], tet[v][2], delta);
		rotateZ( tet[B][0], tet[B][1], tet[B][2], tet[B][0], tet[B][1], tet[B][2], delta);
		rotateZ( tet[C][0], tet[C][1], tet[C][2], tet[C][0], tet[C][1], tet[C][2], delta);
		rotateZ( tet[D][0], tet[D][1], tet[D][2], tet[D][0], tet[D][1], tet[D][2], delta);
		fp = fopen( "tetRotationZ.dat", "w+");
		for( int f=0; f<4; f++) for( int v=0; v<4; v++)	fprintf( fp, "%lf %lf %lf \n", tet[((v==3?0:v)+f)%4][0], tet[((v==3?0:v)+f)%4][1], tet[((v==3?0:v)+f)%4][2]);
		fclose(fp);

		//fprintf( stderr, "Finite Element Integral: %lf\n", fabs(1./120. * tet[B][0]*(tet[C][1] - tet[B][1])*(tet[D][2] - tet[C][2])));



	if( i==j){
		// 1/60 * xB * (yC-yB) * (zD-zC)
		return fabs(1./60. * tet[B][0]*(tet[C][1] - tet[B][1])*(tet[D][2] - tet[C][2]));
	}else{
		// 1/120 * xB * (yC-yB) * (zD-zC)
		return fabs(1./120. * tet[B][0]*(tet[C][1] - tet[B][1])*(tet[D][2] - tet[C][2]));
	}
}

/*int main( int argc, char** argv)
{
	srand(1);
#if DIMENSIONS == 3
	main3DSparse( argc, argv);
	//main3D( argc, argv);
#elif DIMENSIONS == 2
	main2DSparse( argc, argv);
	//main2D( argc, argv);
#endif
	FILE *fp;


	// tet
	fp = fopen( "tetOriginal.dat", "w+");
	double tet[4][3]; // = {{0,0,0},{1,1,1},{1,2,2},{1,2,3}};
	for( int v=0; v<4; v++){
		for( int p=0; p<3; p++){
			tet[v][p] = 1.*myRand();
			//fprintf( fp, "%lf ", tet[v][p]);
		}
		//fprintf( fp, "\n");
	}
	for( int f=0; f<4; f++) for( int v=0; v<4; v++)	fprintf( fp, "%lf %lf %lf \n", tet[((v==3?0:v)+f)%4][0], tet[((v==3?0:v)+f)%4][1], tet[((v==3?0:v)+f)%4][2]);
	fclose(fp);

	double **mass = newDoubleMatrix( 4, 4);
	for( int i=0; i<4; i++){
		for( int j=0; j<4; j++)
			mass[i][j] = getVolumeIntegral( tet, i, j);
	}
	printMatrix( mass, 4, 4, "Mass Matrix B", "%10.5lf ");

	// TET OF REFERENCE
	double refTet[4][3] = {{0,0,0},{1,0,0},{0,1,0},{0,0,1}};

	for( int i=0; i<4; i++){
		for( int j=0; j<4; j++)
			mass[i][j] = getVolumeIntegral( refTet, i, j);
			//fprintf( stderr, " %10.3lf", 120.*getVolumeIntegral( tet, i, j));
		//fprintf( stderr, "\n");
	}
	//fprintf( stderr, "]\n");
	printMatrix( mass, 4, 4, "Reference Mass Matrix B", "%10.5lf ");

	// Jacobian of Tet
	double **jacobi = newDoubleMatrix( 3, 3);
	for( int p=0; p<3; p++)
		for( int v=0; v<3; v++)
			jacobi[p][v] = tet[3][p] - tet[v][p];
	printMatrix( jacobi, 3, 3, "Jacobian J", "%10.5lf ");

	// Determinant
	double detJ = getDeterministe( jacobi, 3);
	fprintf( stderr, "det(J) = %lf\n", getDeterministe( jacobi, 3));
	//fprintf( stderr, "det(J) = %lf\n", solveDeterminant( jacobi, 3));

	matrixScale( mass, mass, detJ, 4, 4);
	printMatrix( mass, 4, 4, "det(J)*B", "%10.5lf ");

	exit( 0);

	// Translate
	for( int p=0; p<3; p++){
		tet[3][p] -= tet[0][p];
		tet[2][p] -= tet[0][p];
		tet[1][p] -= tet[0][p];
		tet[0][p] -= tet[0][p];
	}
	fp = fopen( "tetTranslation.dat", "w+");
	//for( int v=0; v<4; v++)	fprintf( fp, "%lf %lf %lf \n", tet[v][0], tet[v][1], tet[v][2]);
	for( int f=0; f<4; f++) for( int v=0; v<4; v++)	fprintf( fp, "%lf %lf %lf \n", tet[((v==3?0:v)+f)%4][0], tet[((v==3?0:v)+f)%4][1], tet[((v==3?0:v)+f)%4][2]);
	fclose(fp);
	
	double angle;
	// Rotate around x-axis
	/*double angle = atan((tet[2][1]-tet[3][1])/(tet[2][2]-tet[3][2]));
	for( int v=0; v<4; v++)
		rotateX( tet[v][0], tet[v][1], tet[v][2], tet[v][0], tet[v][1], tet[v][2], angle);
		//rotateZYX( tet[v][0], tet[v][1], tet[v][2], tet[v][0], tet[v][1], tet[v][2], angle, 0., 0.);
	fp = fopen( "tetRotationX.dat", "w+");
	for( int f=0; f<4; f++) for( int v=0; v<4; v++)	fprintf( fp, "%lf %lf %lf \n", tet[((v==3?0:v)+f)%4][0], tet[((v==3?0:v)+f)%4][1], tet[((v==3?0:v)+f)%4][2]);
	fclose(fp);
	
	// Rotate around y-axis
	angle = atan((tet[3][0]-tet[2][0])/(tet[2][2]-tet[3][2]));
	for( int v=0; v<4; v++)
		rotateY( tet[v][0], tet[v][1], tet[v][2], tet[v][0], tet[v][1], tet[v][2], angle);
		//rotateZYX( tet[v][0], tet[v][1], tet[v][2], tet[v][0], tet[v][1], tet[v][2], 0., angle, 0.);
	fp = fopen( "tetRotationY.dat", "w+");
	for( int f=0; f<4; f++) for( int v=0; v<4; v++)	fprintf( fp, "%lf %lf %lf \n", tet[((v==3?0:v)+f)%4][0], tet[((v==3?0:v)+f)%4][1], tet[((v==3?0:v)+f)%4][2]);
	fclose(fp);* /

	// Rotate around xy-axis
	double angleA = atan((tet[3][1]-tet[2][1])/(tet[2][2]-tet[3][2]));
	double angleB = atan(
			(tet[3][0]-tet[2][0])
			/( sin(angleA)*(tet[2][1]-tet[3][1]) + cos(angleA)*(tet[3][2]-tet[2][2]) )
		);
	for( int v=0; v<4; v++)
		//rotateY( tet[v][0], tet[v][1], tet[v][2], tet[v][0], tet[v][1], tet[v][2], angle);
		rotateZYX( tet[v][0], tet[v][1], tet[v][2], tet[v][0], tet[v][1], tet[v][2], angleA, angleB, 0.);
	fp = fopen( "tetRotationYX.dat", "w+");
	for( int f=0; f<4; f++) for( int v=0; v<4; v++)	fprintf( fp, "%lf %lf %lf \n", tet[((v==3?0:v)+f)%4][0], tet[((v==3?0:v)+f)%4][1], tet[((v==3?0:v)+f)%4][2]);
	fclose(fp);
	
	// Rotate around z-axis
	angle = atan((tet[1][0]-tet[2][0])/(tet[1][1]-tet[2][1]));
	for( int v=0; v<4; v++)
		rotateZ( tet[v][0], tet[v][1], tet[v][2], tet[v][0], tet[v][1], tet[v][2], angle);
		//rotateZYX( tet[v][0], tet[v][1], tet[v][2], tet[v][0], tet[v][1], tet[v][2], 0., 0., angle);
	fp = fopen( "tetRotationZ.dat", "w+");
	for( int f=0; f<4; f++) for( int v=0; v<4; v++)	fprintf( fp, "%lf %lf %lf \n", tet[((v==3?0:v)+f)%4][0], tet[((v==3?0:v)+f)%4][1], tet[((v==3?0:v)+f)%4][2]);
	fclose(fp);

	
	fprintf( stderr, "Finite Element Integral: %lf\n", fabs(1./120. * tet[1][0]*(tet[2][1] - tet[1][1])*(tet[3][2] - tet[2][2])));
	
	fprintf( stderr, "xA = yA = zA = %lf\n", 0.);
	
	fprintf( stderr, "xB = xC = xD = %lf\n", getHeightNEW( tet[0], tet[1], tet[2], tet[3]));
	
	exit( 0);
	
	// tet
	//double tet[4][3] = {{0,0,0},{1,1,1},{1,2,2},{1,2,3}};
	for( int v=0; v<4; v++){
		for( int p=0; p<3; p++)
			printf( "%lf ", tet[v][p]);
		printf( "\n");
	}

	float xB = tet[1][0],
	      yB = tet[1][1],
	      zB = tet[1][2],
	      zC = tet[2][2],
	      yD = tet[3][1],
	      zD = tet[3][2];
	      

	float dx = 0.05;
	for( float x=0; x<=xB; x+=dx){
		for( float y=yB/xB*x; y<=yD/xB*x; y+=dx){
			for( float z=(zC-zB)/(yD-yB)*y+(zB-yB*(zC-zB)/(yD-yB))*x/xB; z<=(zD-zB)/(yD-yB)*y+(zB-yB*(zD-zB)/(yD-yB))*x/xB; z+=dx)
			{
				printf( "%lf %lf %lf \n", x, y, z);
			}	
		}	
	}


	exit(0);
	main1D( argc, argv);
}*/



double canonicalBaseFunctionGradient( double *x, int dimensions, int dimension, int index, int order)
{
	double ret = 0.;
	if(index == 0){
		for( int d=0; d<dimensions; d++)
			ret += x[d];
		ret = -order*pow( 1. - ret, order-1);
	}else
		if(index-1==dimension)
		ret = order*pow( x[index-1], order);
	return ret;
}

double canonicalBaseFunction( double *x, int dimensions, int index, int order)
{
	double ret = 0.;
	if(index == 0){
		for( int d=0; d<dimensions; d++)
			ret += x[d];
		ret = pow( 1. - ret, order);
	}else
		ret = pow( x[index-1], order);
	return ret;
}


double gaussQuadratureCanonical1D( int orderI, int orderJ, int orderK, int I, int J, int K)
{
	// Exact Integration for polynomials upto order 2n-1
	double a = 0;
	double b = 1;
	int o = orderI + orderJ + orderK;
	int n = (o+1)/2 + (o+1)%2;
	double w[n];
	double p[n];
	//fprintf( stderr, "order=%i -> %i-point Gauss quadrature\n", o, n);

	switch(n){
	case 1:
		// upto order 1
		w[0] = 2;
		p[0] = 0.;
		break;

	case 2:
		// upto order 3
		w[0] = w[1] = 1;
		p[0] = -sqrt(1./3.);
		p[1] =  sqrt(1./3.);
		break;

	case 3:
		// upto order 5
		w[0] = 8./9.;
		p[0] = 0.;
		w[1] = w[2] = 5./9.;
		p[1] = -sqrt(3./5.);//-0.774596669241483;
		p[2] =  sqrt(3./5.);//0.774596669241483;
		break;

	case 4:
		// upto order 7
		w[0] = w[1] = 0.652145154862546;
		p[0] = -0.339981043584856;
		p[1] =  0.339981043584856;

		w[2] = w[3] = 0.347854845137454;
		p[2] = -0.861136311594053;
		p[3] =  0.861136311594053;
		break;
	default:
		fprintf( stderr, "ERROR: No quadrature rule defined for n=%i!\n", n);
		exit(0);
	}

	// map evaluation points from interval [-1,1] to interval [0,1]
	for( int k=0; k<n; k++){
		p[k] = (b-a)/2.*(p[k]+1.);
	}

	double ret = 0.;
	double x[1];
	for( int k=0; k<n; k++){
		//fprintf( stderr, "w[%i]=%lf, x[%i]=%lf\n", k, w[k], k, p[k]);
		//ret += w[k]*f(p[k]);
		x[0] = p[k];
		ret += w[k]
		         * canonicalBaseFunction( x, 1, I, orderI)
		         * canonicalBaseFunction( x, 1, J, orderJ)
		         * canonicalBaseFunction( x, 1, K, orderK);
	}

	/*{
		char filename[512];
		sprintf( filename, "baseFunction%i%i%i.dat", I,J,K);
		FILE *fp = fopen( filename, "w+");
		for( double x[1] = {0}; x[0]<=1.; x[0]+=0.01)
			fprintf( fp, "%lf %lf\n", x[0],
					   canonicalBaseFunction( x, 1, I, orderI)
			         * canonicalBaseFunction( x, 1, J, orderJ)
			         * canonicalBaseFunction( x, 1, K, orderK));
	}*/

	// det(J_e) * I_0 + E
	return (b-a)/2.*ret;
}
/*****************************************************************************/



