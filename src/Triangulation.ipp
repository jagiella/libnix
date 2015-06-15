// C HEADER
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#ifdef OMP
  #include <omp.h>
#else
  #warning "no openmp support"
#endif

// C++ HEADER
#include <fstream>

// CUSTOM HEADER
#include "Triangulation.hpp"
#include "LinearAlgebra.hpp"
#include "IO.hpp"

#define signof( a) (a==0 ? 0 : (a<0 ? -1 : 1))

template <int Dimensions, class VertexType>
inline void Simplex<Dimensions, VertexType>::getCircumCenter( POSITION_T center[Dimensions])
{
	switch( Dimensions){
	case 2:
		getCircumCenter2D( center);
		return;
	case 3:
		getCircumCenter3D( center);
		return;
	}
}

template <int Dimensions, class VertexType>
inline void Simplex<Dimensions, VertexType>::getCircumCenter3D( POSITION_T center[Dimensions])
{
	//Simplex<Dimensions, VertexType> *tet = this;

//        |d-a|^2 [(b-a)x(c-a)] + |c-a|^2 [(d-a)x(b-a)] + |b-a|^2 [(c-a)x(d-a)]
//m = a + ---------------------------------------------------------------------.
//                                | bx-ax  by-ay  bz-az |
//                              2 | cx-ax  cy-ay  cz-az |
//                                | dx-ax  dy-ay  dz-az |
	POSITION_T *a = this->vertices[0]->pos(),
	       *b = this->vertices[1]->pos(),
	       *c = this->vertices[2]->pos(),
	       *d = this->vertices[3]->pos();
	
	POSITION_T ba[3], ca[3], da[3];
	vectorDifference3D( b, a, ba);
	vectorDifference3D( c, a, ca);
	vectorDifference3D( d, a, da);
	
	POSITION_T bc[3], cd[3], db[3];
	Vector<POSITION_T>::crossProduct3D( (POSITION_T*)ba, (POSITION_T*)ca, (POSITION_T*)bc);
	Vector<POSITION_T>::crossProduct3D( (POSITION_T*)ca, (POSITION_T*)da, (POSITION_T*)cd);
	Vector<POSITION_T>::crossProduct3D( (POSITION_T*)da, (POSITION_T*)ba, (POSITION_T*)db);
	
	POSITION_T bcd = 2. * Vector<POSITION_T>::dotProduct3D( ba, cd);

	Vector<POSITION_T>::vectorScale3D( bc, Vector<POSITION_T>::dotProduct3D( da, da), bc);
	Vector<POSITION_T>::vectorScale3D( cd, Vector<POSITION_T>::dotProduct3D( ba, ba), cd);
	Vector<POSITION_T>::vectorScale3D( db, Vector<POSITION_T>::dotProduct3D( ca, ca), db);
		
	center[0] = a[0] + (bc[0] + cd[0] + db[0])/bcd;
	center[1] = a[1] + (bc[1] + cd[1] + db[1])/bcd;
	center[2] = a[2] + (bc[2] + cd[2] + db[2])/bcd;
}

template <int Dimensions, class VertexType>
inline void Simplex<Dimensions, VertexType>::getCircumCenter2D( POSITION_T center[Dimensions])
{
	//Simplex<Dimensions, VertexType> *tet = this;

	POSITION_T *a = this->vertices[0]->pos(),
	           *b = this->vertices[1]->pos(),
	           *c = this->vertices[2]->pos();

	POSITION_T D = 2.*( a[0] * (b[1] - c[1]) + b[0] * ( c[1] - a[1] ) + c[0] * ( a[1] - b[1] ) );

	    center[0]  = ( a[1]*a[1] + a[0]*a[0] )*(b[1] - c[1]);
	    center[0] += ( b[1]*b[1] + b[0]*b[0] )*(c[1] - a[1]);
	    center[0] += ( c[1]*c[1] + c[0]*c[0] )*(a[1] - b[1]);
	    center[0] /= D;

	    center[1]  = ( a[1]*a[1] + a[0]*a[0] )*(c[0] - b[0]);
	    center[1] += ( b[1]*b[1] + b[0]*b[0] )*(a[0] - c[0]);
	    center[1] += ( c[1]*c[1] + c[0]*c[0] )*(b[0] - a[0]);
	    center[1] /= D;

}

/*template <int Dimensions>
Vertex<Dimensions> *newGridPoint()
{
	return (Vertex<Dimensions> *) malloc( sizeof(Vertex<Dimensions>));
}*/

template <int Dimensions, class VertexType, class SimplexType>
void Triangulation<Dimensions, VertexType, SimplexType>::setVoronoiGrid()
{
	fprintf( stderr, "Triangulation<Dimensions>::setVoronoiGrid()\n");

	// init points
	this->voronoiGridPoints = (VertexType *) malloc( this->countTetrahedra * sizeof(VertexType));
	
	// init links
	this->voronoiCellToVoronoiGridPoints = (VertexType ***) malloc( this->countVertices * sizeof(VertexType**));
	this->voronoiCellCountVoronoiGridPoints = (int *) malloc( this->countVertices * sizeof(int));
	for( int i=0; i<this->countVertices; i++){
		this->voronoiCellToVoronoiGridPoints[i] = NULL;
		this->voronoiCellCountVoronoiGridPoints[i] = 0;
	}

	// set points
	for( int i=0; i<this->countTetrahedra; i++){
		// set index
		this->voronoiGridPoints[i].index = i;
		
		// set point position
		this->tetrahedra[i]->getCircumCenter( this->voronoiGridPoints[i].position);

		// set neighborhoods of points
		// ...
	}
	
	// set links
	for( int t=0; t<this->countTetrahedra; t++){
		for( int v=0; v<Dimensions+1; v++){
			// for each tet vertice
			int i = this->tetrahedra[t]->vertices[v]->index;
			int l;
			for( l=0; l<this->voronoiCellCountVoronoiGridPoints[i] && this->voronoiCellToVoronoiGridPoints[i][l]!=&(this->voronoiGridPoints[t]); l++);
			if( l==this->voronoiCellCountVoronoiGridPoints[i]){
				this->voronoiCellCountVoronoiGridPoints[i]++;
				this->voronoiCellToVoronoiGridPoints[i] = (VertexType **)realloc( this->voronoiCellToVoronoiGridPoints[i], this->voronoiCellCountVoronoiGridPoints[i]*sizeof(VertexType*));
				this->voronoiCellToVoronoiGridPoints[i][l] = &(this->voronoiGridPoints[t]);
			}
		}
	}

	// test output of links
	/*for( int i=0; i<this->countVertices; i++){
		fprintf(stderr, "%i: ", i);
		for( int l=0; l<this->voronoiCellCountVoronoiGridPoints[i]; l++)
			fprintf(stderr, "%i ", this->voronoiCellToVoronoiGridPoints[i][l]->index);
		fprintf(stderr, "\n");
	}*/

	
}
/*template <int Dimensions>
void Vertex<Dimensions>::validate()
{
	int i;
	
	int countFreeNeighbors = 0;
	for( i=0; i<this->countNeighbors; i++){
		if( this->neighbors[i]->isFree())
			countFreeNeighbors++;
	}
	if( countFreeNeighbors!=this->countFreeNeighborCells){
		fprintf( stderr, "ERROR: Vertex %i has wrong number of free neighbors: %i (correct: %i)\n", this->index, this->countFreeNeighborCells, countFreeNeighbors);
		exit( 0);
	}

	countFreeNeighbors = 0;
	for( i=0; i<this->countExtendedNeighborCells; i++){
		if( this->extendedNeighborhood[i]->isFree())
			countFreeNeighbors++;
	}
	if( countFreeNeighbors!=this->countFreeExtendedNeighborCells){
		if( !this->isFree())
			fprintf( stderr, "ERROR: Vertex %i (agent %i) has wrong number of free extended neighbors: %i/%i (correct: %i)\n",
			         this->index, GetAgent( this)->index, this->countFreeExtendedNeighborCells, this->countExtendedNeighborCells, countFreeNeighbors);
		else{
			fprintf( stderr, "ERROR: Vertex %i has wrong number of free extended neighbors: %i/%i (correct: %i)\n",
			         this->index, this->countFreeExtendedNeighborCells, this->countExtendedNeighborCells, countFreeNeighbors);			
		}

			for( i=0; i<this->countExtendedNeighborCells; i++){
				int found = FALSE;
				int ii;
				for( ii=0; ii<this->extendedNeighborhood[i]->countExtendedNeighborCells; ii++){
					if( this->extendedNeighborhood[i]->extendedNeighborhood[ii] == this){
						found = TRUE;
					//fprintf( stderr, "ERROR: Vertex %i -> extended neighbor %i -> extended neighbor %i :-)\n",
					//         this->index, this->extendedNeighborhood[i]->index, this->extendedNeighborhood[i]->extendedNeighborhood[ii]->index);
					}
				}			
				if( !found){
					fprintf( stderr, "ERROR: Vertex %i has a extended neighbor %i which doesn't contain %i as extended neighbor!\n", this->index, this->extendedNeighborhood[i]->index, this->index);
				}
			}

		exit( 0);
	}


}*/
/*****************************************************************************/

#define FREE	0


#ifdef USE_AGENT
int Vertex<Dimensions>::getState()
{
	return (this->agent==NULL ? FREE : GetAgent( this)->state);	
}
/*****************************************************************************/


int Vertex<Dimensions>::isFree()
{
	return (GetAgent( this)==NULL ? TRUE : GetAgent( this)->isFree());	
}
#endif
/*****************************************************************************/


template <typename T>
int isElementOf( T **pointList, int nrOfPoints, T *point){
	int i;

	//printf("{to search point:%i,list(%i): ", point->nr, nrOfPoints);

	for( i=0; i<nrOfPoints; i++){
		//printf("%i, ", pointNrList[i]);
		if( pointList[i] == point) {
			//printf("[found]}");
			return 1;
		}
	}
	//printf("}");

	return 0;
}
/*****************************************************************************/



template <int Dimensions, class VertexType, class SimplexType>
VertexType* Triangulation<Dimensions,VertexType,SimplexType>::searchForVoronoiCell(
		int countIgnoredPoints,
		VertexType** ignoredPoints,
		int countSourcePoints, VertexType** sourcePoints,
		int &countExploredPoints,
		VertexType** exploredPoints){
	int i, j;

	/*fprintf( stderr, "Ignored Points: [");
	for( j=0; j<countIgnoredPoints; j++)
		fprintf( stderr, "%i ", ignoredPoints[j]->index);
	fprintf( stderr, "\n");

	fprintf( stderr, "Source Points: [");
	for( j=0; j<countSourcePoints; j++)
		fprintf( stderr, "%i ", sourcePoints[j]->index);
	fprintf( stderr, "\n");
	*/

	for( i = 0; i < countSourcePoints; i++){ // for every element in the previous list
		for( j = 0; j < sourcePoints[i]->countNeighbors; j++){ // for every neighbor of a element in the previous list
        //printf("  candidates:%i\n", (voronoiDiagram[ point->reachableNeighbors[ base_index + i ] ].neighbors[j])->nr);
			VertexType *candidate = sourcePoints[i]->neighbors[j];
			if( !isElementOf( ignoredPoints,  // pointer to list of all generations
			                  countIgnoredPoints,     // number of points in the list
			                  candidate) &&
			    !isElementOf( sourcePoints,  // pointer to list of all generations
			                  countSourcePoints,     // number of points in the list
			                  candidate) &&
			    !isElementOf( exploredPoints,  // pointer to list of all generations
			                  countExploredPoints,     // number of points in the list
			                  candidate)
			    ){                                                                       // actual neighbor of neighbor 
			    exploredPoints[countExploredPoints] = candidate;
				countExploredPoints++;
			}
		}
	}
  
	return 0;
}
/*****************************************************************************/


inline char *long2string( char *str, int number)
{
	int base = 128 - 3;
	int i = 0;
	// 4 = [100]
	int mod;
	
	do{
		mod = number%base;
		number = number/base;
		str[i] = mod+3;
		//fprintf( stderr, "%i|", mod+3);
		i++;
	}while(number>0);
	str[i] = '\0';
	
	// resort
	int ii;
	int temp;
	for( ii=0; ii<i/2; ii++){
		temp = str[ii];
		str[ii] = str[i-ii-1];
		str[i-ii-1] = temp;
	}
	
	return str;
}

inline int string2long( char *str)
{
	int base = 128 - 3;
	int i=0;
	int number = 0;
	
	while( str[i]!='\0'){
		number *= base;
		number += (str[i]-3);
		i++;
	}
	
	return number;

}

inline int stringlength( char *str)
{
	int i=0;
	
	while( str[i]!='\0'){
		i++;
	}
	
	return i;
}
/*****************************************************************************/


template <int Dimensions, class VertexType, class SimplexType>
 VertexType *Triangulation<Dimensions, VertexType, SimplexType>::searchClosestFreeVoronoiCell( VertexType *explorationCenter){

/*	int j;
	Vertex<Dimensions>* tempExplorationList[this->countVertices];// = ( Vertex<Dimensions>**) calloc ( this->countVertices, sizeof( Vertex<Dimensions>*));
	//Vertex<Dimensions>** tempExplorationList = ( Vertex<Dimensions>**) calloc ( this->countVertices, sizeof( Vertex<Dimensions>*));
	int countExploredCells=0, countIgnoredPoints=0, countSourcePoints=0, countExploredPoints=0;
	int indexIgnoredPoints, indexSourcePoints, indexExploredPoints;

	// 0st generation (Set point itself)
	tempExplorationList[0]  = explorationCenter;        // neighbors in 0th generation
	countExploredCells = 1;                   // number of neighbors in 0th generation
	indexIgnoredPoints = 0;
	countIgnoredPoints = 1;
	//printf("Add: %i to Ignored\n", tempExplorationList[0]->index);


	// 1st generation (Set direct neighbors)
	for(j = 0; j < explorationCenter->countNeighbors; j++){
		tempExplorationList[countExploredCells++]  = explorationCenter->neighbors[j];        // neighbors in 0th generation
		if(explorationCenter->countExtendedNeighborCells>0){
			//printf("Add: %i to Ignored\n", tempExplorationList[countExploredCells-1]->index);
			countIgnoredPoints++;
		}else{
			//printf("Add: %i to Source\n", tempExplorationList[countExploredCells-1]->index);
			countSourcePoints++;
		}
	}

	// following generations (Set extended neighbors)
	for(j = 0; j < explorationCenter->countExtendedNeighborCells; j++){
		tempExplorationList[countExploredCells++]  = explorationCenter->extendedNeighborhood[j];        // neighbors in 0th generation
		//printf("Add: %i to Source\n", tempExplorationList[countExploredCells-1]->index);
		countSourcePoints++;
	}
	indexSourcePoints = countIgnoredPoints;
	
	// get all point below a certain radius
	countExploredPoints = 0;
	indexExploredPoints = countIgnoredPoints + countSourcePoints;	

	do{
		//printf("Ignored: count %i, index %i; Source: count %i, index %i;\n", countIgnoredPoints, indexIgnoredPoints, countSourcePoints, indexSourcePoints);
		//fprintf( stderr, ". \b");
		
		// explore neighborhood
		Triangulation<Dimensions>::searchForVoronoiCell( countIgnoredPoints,  &tempExplorationList[indexIgnoredPoints],
		                            countSourcePoints,   &tempExplorationList[indexSourcePoints], 
		                            countExploredPoints, &tempExplorationList[indexExploredPoints]);
		
		// search closest free point
		Vertex<Dimensions>* closestPoint = NULL;
		double distanceToClosestPoint = 0.;
		double tempDistance;
		for( j=0; j<countExploredPoints; j++){
			if( tempExplorationList[indexExploredPoints+j]->getState()==FREE){
				tempDistance = explorationCenter->getDistanceTo( tempExplorationList[indexExploredPoints+j]);
				if( tempDistance < distanceToClosestPoint || closestPoint==NULL){
					distanceToClosestPoint = tempDistance;
					closestPoint = tempExplorationList[indexExploredPoints+j];
				}
			}
		}
		if( closestPoint!=NULL){
			//fprintf( stderr, "clostest cell to cell (%i) is %i\n", explorationCenter->index, closestPoint->index);
			return closestPoint;
		}
		
		// initialize next
		indexIgnoredPoints+= countIgnoredPoints;
		countIgnoredPoints = countSourcePoints;
		
		indexSourcePoints+= countSourcePoints;
		countSourcePoints = countExploredPoints;
		
		indexExploredPoints+= countExploredPoints;
		countExploredPoints = 0;
	}while( countSourcePoints > 0);

	
	//free( tempExplorationList);
	*/
	fprintf( stderr, "WARNING: No free cell was found which is close to cell %i\n", explorationCenter->index);
	return NULL;
}
/*****************************************************************************/



template <int Dimensions, class VertexType, class SimplexType>
 VertexType *Triangulation<Dimensions, VertexType, SimplexType>::searchClosestVoronoiCell( VertexType *explorationCenter, double targetPosition[Dimensions]){

	double min_dist = explorationCenter->getDistanceTo( targetPosition);
	VertexType *candidate = explorationCenter;

	do{
		explorationCenter = candidate;
		
		for( int i=0; i<explorationCenter->countNeighbors; i++){
			double dist = explorationCenter->neighbors[i]->getDistanceTo( targetPosition);
			if( dist < min_dist){
				min_dist = dist;
				candidate = explorationCenter->neighbors[i];
			}
		}	
	}while( candidate!=explorationCenter);
	
	return candidate;
	
	//free( tempExplorationList);
	fprintf( stderr, "WARNING: No free cell was found which is close to cell %i\n", explorationCenter->index);
	return NULL;
}
/*****************************************************************************/


template <int Dimensions>
double Vertex<Dimensions>::getDistanceTo( Vertex<Dimensions>* cell)
{
	double distance = 0.;
	
	for( int i=0; i<Dimensions; i++)
		distance += pow( this->position[i] - cell->position[i], 2);
	
	return sqrt( distance);
}

template <int Dimensions, class VertexType, class SimplexType>
void Triangulation<Dimensions,VertexType,SimplexType>::setDomain()
{
	if( !domainSet){
		POSITION_T xMin[Dimensions];
		POSITION_T xMax[Dimensions];
		//int xN[Dimensions];

		if( countVertices == 0){
			for( int dim = 0; dim<Dimensions; dim++){
				xMin[dim] = 0;
				xMax[dim] = 0;
			}
			return;
		}

		int i = 0, dim;
	
		// domain borders
		for( dim = 0; dim<Dimensions; dim++){
			xMin[dim] = vertices[i]->x(dim);
			xMax[dim] = vertices[i]->x(dim);
		}

		for( i = 1; i < countVertices; i++){
			for( dim = 0; dim<Dimensions; dim++){
				if( xMin[dim] > vertices[i]->x(dim))
					xMin[dim] = vertices[i]->x(dim);
				if( xMax[dim] < vertices[i]->x(dim))
					xMax[dim] = vertices[i]->x(dim);
			}
		}
	
		// boundary thickness
		boundaryThickness = 1.;
		for( dim = 0; dim<Dimensions; dim++)
			boundaryThickness *= (xMax[dim]-xMin[dim]); 
		boundaryThickness = pow( boundaryThickness / countVertices ,1./Dimensions);
					
		// elements in each direction
		/*for( dim = 0; dim<Dimensions; dim++)
			//xN[dim] = (int)(xMax[dim]-xMin[dim] + 1.);
			xN[dim] = ceil(xMax[dim])-ceil(xMin[dim])+1;*/
		
		domain = new Domain<Dimensions>(xMin,xMax);
		domainSet = TRUE;
	}
}
/*****************************************************************************/
  

template <int Dimensions, class VertexType, class SimplexType>
inline void Triangulation<Dimensions, VertexType, SimplexType>::setDomain( double xmin, double xmax, double ymin, double ymax, double zmin, double zmax)
{
	/*xMin[0] = xmin;
	xMax[0] = xmax;
	xMin[1] = ymin;
	xMax[1] = ymax;
	xMin[2] = zmin;
	xMax[2] = zmax;*/
	domain = new Domain<3>(xmin,  xmax,  ymin,  ymax,  zmin,  zmax);

	domainSet = TRUE;

}
/*****************************************************************************/


template <int Dimensions, class VertexType, class SimplexType>
VertexType * Triangulation<Dimensions,VertexType,SimplexType>::getCentralCell(){
	Triangulation<Dimensions,VertexType,SimplexType> *voronoiDiagram = this;
	double min[Dimensions], max[Dimensions], mid[Dimensions];
	
	// minima, maxima and average of x and y
	for( int ii=0; ii<Dimensions; ii++)
		min[ii] = max[ii] = voronoiDiagram->vertices[0]->x(ii);
  
 	for( int i = 1; i < voronoiDiagram->countVertices; i++){
		for( int ii=0; ii<Dimensions; ii++){
			if( max[ii] < voronoiDiagram->vertices[i]->x(ii))
				max[ii] = voronoiDiagram->vertices[i]->x(ii);
			if( min[ii] > voronoiDiagram->vertices[i]->x(ii) )
				min[ii] = voronoiDiagram->vertices[i]->x(ii);
		}
	}
  
	for( int ii=0; ii<Dimensions; ii++)
		mid[ii] = (min[ii] + max[ii]) / 2;
  
	// central cell
	VertexType* central_cell = voronoiDiagram->vertices[0];
	double distance = 0.;
	for( int ii=0; ii<Dimensions; ii++)
		distance += pow( mid[ii] - voronoiDiagram->vertices[0]->x(ii), 2.);
	distance = sqrt( distance);
	
	for( int i = 1; i < voronoiDiagram->countVertices; i++){
		double temp_distance = 0.;
		for( int ii=0; ii<Dimensions; ii++)
			temp_distance += pow( mid[ii] - voronoiDiagram->vertices[i]->x(ii), 2.);
		temp_distance = sqrt( temp_distance);
		if( distance > temp_distance){
			central_cell = voronoiDiagram->vertices[i];
			distance = temp_distance;
		}
	}
  
	//fprintf( stderr, "central cell: %i, %lf\n", central_cell->index, central_cell->position[0]);
  
	return central_cell;
}
/*****************************************************************************/


template <int Dimensions, class VertexType>
void Simplex<Dimensions,VertexType>::initCircumsphereRadius()
{
	Simplex<Dimensions,VertexType>* tet = this;
	int k, l;
	int NR_Simplex_POINTS = Dimensions+1;
	//double A[NR_Simplex_POINTS][NR_Simplex_POINTS]/*={{0., 0., 0., 0.},{0., 0., 0., 0.},{0., 0., 0., 0.},{0., 0., 0., 0.}}*/;
	//double **A = newDoubleMatrix( NR_Simplex_POINTS, NR_Simplex_POINTS);
	//double x[NR_Simplex_POINTS];

	// initialize matrix and vector of linear system
	for( k=0; k<NR_Simplex_POINTS; k++){
		tet->circumsphere[k]=0.;
		for( l=0; l<Dimensions; l++){
			tet->circumsphere[k] -= pow( tet->getVertex(k)->x(l), 2);
			//tet->circumsphere[k] -= myPow( tet->vertices[k]->position[l], 2);
			Simplex<Dimensions,VertexType>::_A[k][l] = tet->getVertex(k)->x(l);
			//A[k][l] = tet->vertices[k]->position[l];
		}
		Simplex<Dimensions,VertexType>::_A[k][l] = 1.;
		//A[k][l] = 1.;
	}
	Matrix<double>::solveLinearSystem(
			(double**)Simplex<Dimensions,VertexType>::_A,
			tet->circumsphere,
			(double*)Simplex<Dimensions,VertexType>::_x, NR_Simplex_POINTS,
			(double**)Simplex<Dimensions,VertexType>::_B);
	//solveLinearSystem( (double**)A, tet->circumsphere, x, NR_Simplex_POINTS);
	tet->radius = 0.;
	for( l=0; l<Dimensions; l++)
		tet->radius += tet->circumsphere[l]*tet->circumsphere[l];
	tet->radius = tet->radius/4. - tet->circumsphere[l];
	tet->circumsphereInitialized = 1;

	//deleteDoubleMatrix( A, NR_Simplex_POINTS);

	// NUMERIC ERROR

	/*double R  = 0.;
	for( l=0; l<Dimensions; l++)
		//R += pow(tet->vertices[0]->position[l] + tet->circumsphere[l]/2., 2);
		R += myPow(tet->vertices[0]->position[l] + tet->circumsphere[l]/2., 2);
	//double maxNumericError = R - tet->radius;
	double avNumericError = R - tet->radius;
	for( k=1; k<NR_Simplex_POINTS; k++){
		R  = 0.;
		for( l=0; l<Dimensions; l++)
			//R += pow(tet->vertices[k]->position[l] + tet->circumsphere[l]/2., 2);
			R += myPow(tet->vertices[k]->position[l] + tet->circumsphere[l]/2., 2);
		tet->numericError = R - tet->radius;
		fprintf( stderr, "INFO: numerical error of tet %i for vertice %i: %e\n", tet->index, tet->vertices[k]->index, tet->numericError );
		//if( maxNumericError < tet->numericError)
		//	maxNumericError = tet->numericError;
		avNumericError += R - tet->radius;
	}
	//tet->numericError = maxNumericError;
	tet->numericError = avNumericError/NR_Simplex_POINTS;
*/
	// END NUMERIC ERROR
}
/*****************************************************************************/


template <int Dimensions,class VertexType>
void Simplex<Dimensions,VertexType>::getCircumsphereRadius( long double *circumsphere, long double &radius)
{
	Simplex<Dimensions,VertexType>* tet = this;

	int k, l;
	int NR_Simplex_POINTS = Dimensions+1;

	//long double A[NR_Simplex_POINTS][NR_Simplex_POINTS]/*={{0., 0., 0., 0.},{0., 0., 0., 0.},{0., 0., 0., 0.},{0., 0., 0., 0.}}*/;
	//long double B[NR_Simplex_POINTS][NR_Simplex_POINTS]/*={{0., 0., 0., 0.},{0., 0., 0., 0.},{0., 0., 0., 0.},{0., 0., 0., 0.}}*/;
	long double **A = newLongDoubleMatrix(NR_Simplex_POINTS,NR_Simplex_POINTS);
	long double **B = newLongDoubleMatrix(NR_Simplex_POINTS,NR_Simplex_POINTS);
	long double x[NR_Simplex_POINTS];

	// initialize matrix and vector of linear system
	for( k=0; k<NR_Simplex_POINTS; k++){
		circumsphere[k]=0.;
		for( l=0; l<Dimensions; l++){
			circumsphere[k] -= pow( tet->getVertex(k)->x(l), 2);
			A[k][l] = tet->getVertex(k)->x(l);
		}
		A[k][l] = 1.;
	}
	solveLinearSystemB( (long double**)A, circumsphere, (long double*)x, NR_Simplex_POINTS, (long double**)B);

	for( l=0; l<Dimensions; l++)
		radius += circumsphere[l]*circumsphere[l];
	radius = radius/4. - circumsphere[l];

	deleteLongDoubleMatrix( A, NR_Simplex_POINTS);
	deleteLongDoubleMatrix( B, NR_Simplex_POINTS);
}
/*****************************************************************************/


template <int Dimensions,class VertexType>
double Simplex<Dimensions,VertexType>::getCircumsphereRadius()
{
	Simplex<Dimensions,VertexType>* tet = this;
	// initialize circumsphere
	if( tet->circumsphereInitialized == 0)
		tet->initCircumsphereRadius( );

	return tet->radius;
}
/*****************************************************************************/


template <int Dimensions> template <class SimplexType>
double Vertex<Dimensions>::getDistanceToCircumsphereCenter( SimplexType* tet)
{
	Vertex<Dimensions>* point = this;
	int l;

	// initialize circumsphere
	if( tet->circumsphereInitialized == 0)
		tet->initCircumsphereRadius();

	double R  = 0.;
	for( l=0; l<Dimensions; l++)
		R += pow(point->position[l] + tet->circumsphere[l]/2., 2);
		//R += myPow(point->position[l] + tet->circumsphere[l]/2., 2);

	return R;
}
/*****************************************************************************/


template <int Dimensions> template <class SimplexType>
double Vertex<Dimensions>::getDistanceToCircumsphere( SimplexType* tet)
{
	Vertex<Dimensions>* point = this;
	int l;

	// initialize circumsphere
	if( tet->circumsphereInitialized == 0)
		tet->initCircumsphereRadius( );

	DISTANCE_T R  = 0.;
	for( l=0; l<Dimensions; l++)
		//R += pow((float)point->position[l] + (float)tet->circumsphere[l]/2., 2);
		R += myPow(point->position[l] + tet->circumsphere[l]/2., 2);

	if( isnan( R - tet->radius )){
		fprintf( stderr, "dist = %e\n", (R - tet->radius));
		fprintf( stderr, "R = %e\n", R);
		fprintf( stderr, "tet->radius = %e\n", tet->radius);
		for( l=0; l<Dimensions; l++){
			//R += myPow(point->position[l] + tet->circumsphere[l]/2., 2);
			fprintf( stderr, "point->position[%i] = %e\n", l, point->position[l]);
			fprintf( stderr, "tet->circumsphere[%i] = %e\n", l, tet->circumsphere[l]);
		}

		exit( 0);
	}

	if( fabs(R - tet->radius)<1e-7)
	//if( fabs(R - tet->radius)<1e-9)
	{

		HIGH_PRECISION_DISTANCE_T circumsphere[Dimensions+1];
		HIGH_PRECISION_DISTANCE_T radius = 0.;
		tet->getCircumsphereRadius( circumsphere, radius);
		HIGH_PRECISION_DISTANCE_T ldR  = 0.;
		for( l=0; l<Dimensions; l++)
			ldR += ((HIGH_PRECISION_DISTANCE_T)point->position[l] + circumsphere[l]/2.)*((HIGH_PRECISION_DISTANCE_T)point->position[l] + circumsphere[l]/2.);
		//fprintf( stderr, "WARNING: dist = %.10e ", (R - tet->radius));
		//fprintf( stderr, "----> ld dist = %.10e ", (double)(ldR - radius));
		//fprintf( stderr, "------> error = %.10e\n", fabs((double)(ldR - radius)-(R - tet->radius)));
		//if( signof(ldR - radius) != signof(R - tet->radius))
		if( (ldR - radius)*(R - tet->radius) < 0.)
		{
			fprintf( stderr, "ERROR: wrong sign!\n");
			fprintf( stderr, "WARNING: dist = %.10e\n", (R - tet->radius));
			fprintf( stderr, "----> ld dist = %.10e ", (double)(ldR - radius));
			fprintf( stderr, "------> error = %.10e\n", fabs((double)(ldR - radius)-(R - tet->radius)));
			return (DISTANCE_T)(ldR - radius);
			//exit(0);
		}

		/*long double Rl = 0.;
		for( l=0; l<Dimensions; l++)
			//R += pow((float)point->position[l] + (float)tet->circumsphere[l]/2., 2);
			R += ((long double)point->position[l] + (long double)tet->circumsphere[l]/2.) * ((long double)point->position[l] + (long double)tet->circumsphere[l]/2.);
		fprintf( stderr, "-------> dist = %.10e\n", (Rl - (long double)tet->radius));*/
	}
	return (R - tet->radius)/* - tet->numericError*/;
}
/*****************************************************************************/


template <int Dimensions, class VertexType, class SimplexType>
SimplexType* Triangulation<Dimensions,VertexType,SimplexType>::getSimplexContainingPointInCircumSphere( VertexType* point)
{
	Triangulation<Dimensions,VertexType> *voronoiDiagram = this;

	int i;

	// look for first tet containing new point
	//fprintf( stderr, "look for first tet containing new point\n");
	SimplexType *actualTet = NULL,
	            *candidateTet = voronoiDiagram->tetrahedra[(int)(myRandE(1.)*(double)voronoiDiagram->countTetrahedra)];
	double distanceToCircumsphere;
	double minDistanceToCircumsphere = point->getDistanceToCircumsphere( candidateTet);

	// chose closest random point as start point
	for( i=0; i<5; i++){
		actualTet = voronoiDiagram->tetrahedra[(int)(myRandE(1.)*(double)voronoiDiagram->countTetrahedra)];
		distanceToCircumsphere = point->getDistanceToCircumsphere( actualTet);
		if( minDistanceToCircumsphere > distanceToCircumsphere){
			minDistanceToCircumsphere = distanceToCircumsphere;
			candidateTet = actualTet;
		}
	}
	actualTet = NULL;

	// follow shortest way to first Simplex containing new point
	while( minDistanceToCircumsphere > 0.){
	//	fprintf( stderr, "INFO: actualTet!=candidateTet ti=%i\n", ti);
		actualTet=candidateTet;
		//fprintf( stderr, "INFO: actualTet(%i)! ti=%i\n", actualTet->index, ti);
		//fprintf( stderr, "INFO: actualTet(%i)==candidateTet(%i)! ti=%i\n", actualTet->index, candidateTet->index, ti);
		for( i=0; i<actualTet->countNeighborTetrahedra; i++){
			// circumsphere contains point i
			distanceToCircumsphere = point->getDistanceToCircumsphere( actualTet->neighborTetrahedra[i]);
			if( minDistanceToCircumsphere > distanceToCircumsphere){
				minDistanceToCircumsphere = distanceToCircumsphere;
				candidateTet = actualTet->neighborTetrahedra[i];
			}
		}
		if( actualTet==candidateTet){
			//fprintf( stderr, "INFO: actualTet(%i)==candidateTet(%i)! dist=%e\n", actualTet->index, candidateTet->index, minDistanceToCircumsphere);
			//exit(0);
			candidateTet = actualTet->neighborTetrahedra[(int)(myRand()*actualTet->countNeighborTetrahedra)];
			//this->checkDelaunayCondition();
		}

	}

	return candidateTet;
}

/*****************************************************************************/


template <int Dimensions, class VertexType>
int Simplex<Dimensions,VertexType>::countCommonPointOfTetrahedra( Simplex<Dimensions,VertexType>* tetB)
{
	Simplex<Dimensions,VertexType>* tetA = this;

	int l, m;
	int nrCommonPoints = 0;
		
	for( l=0; l<Dimensions+1; l++) // for all points of j
		for( m=0; m<Dimensions+1; m++) // for all points of k
			if( tetA->vertices[l] == tetB->vertices[m])
				nrCommonPoints++;
				
	return nrCommonPoints;
}
/*****************************************************************************/


template <int Dimensions, class VertexType>
int Simplex<Dimensions,VertexType>::removeSimplexNeighbor( Simplex<Dimensions,VertexType>* neighborTet)
{
	Simplex<Dimensions,VertexType>* tet = this;
	int i;

	//fprintf( stderr, "INFO: removeNeighbor(): remove %i as neighbor of %i!\n", neighborCell->index, voronoiCell->index);
	
	for( i=0; i<tet->countNeighborTetrahedra; i++){
		if( tet->neighborTetrahedra[i] == neighborTet){
			tet->countNeighborTetrahedra--;
			tet->neighborTetrahedra[i] = tet->neighborTetrahedra[tet->countNeighborTetrahedra];
			//tet->neighbors = (Vertex<Dimensions>**) realloc( tet->neighbors, tet->countNeighbors*sizeof(Vertex<Dimensions>*));
			return 1;
		}	
	}

	fprintf( stderr, "ERROR: removeSimplexNeighbor() can not find the neighbor to delete!\n");
	for( i=0; i<tet->countNeighborTetrahedra; i++){
		fprintf( stderr, "neighbor %i\n", tet->neighborTetrahedra[i]->index);
	}
	exit( 0);
	return 0;
}
/*****************************************************************************/


template <int Dimensions,class VertexType>
int Simplex<Dimensions,VertexType>::addSimplexNeighbor( Simplex<Dimensions,VertexType>* neighborTet)
{
	this->neighborTetrahedra[this->countNeighborTetrahedra] = neighborTet;
	this->countNeighborTetrahedra++;	

	return 1;
}
/*****************************************************************************/


template <int Dimensions,class VertexType>
int Simplex<Dimensions,VertexType>::addOrReplaceSimplexNeighbor( Simplex<Dimensions,VertexType>* neighborTet)
{
	Simplex<Dimensions,VertexType>* tet = this;
	int i;
	for( i=0; i<tet->countNeighborTetrahedra; i++){
		if( tet->neighborTetrahedra[i] == neighborTet){
			//fprintf( stderr, "ERROR: addSimplexNeighbor() found neighbor being already included!\n");
			//exit( 0);
			return 0;
		}	
	}
	

	//printf("realloc: %i -> %i\n", voronoiCell->countNeighbors, voronoiCell->countNeighbors+1);
	//if( voronoiCell->countNeighbors != 0)
		//tet->neighbors = (Vertex<Dimensions>**) realloc( voronoiCell->neighbors, (voronoiCell->countNeighbors+1)*sizeof(Vertex<Dimensions>*));
	//else
	//	voronoiCell->neighbors = (Vertex<Dimensions>**) malloc( sizeof());
	//printf("set new neighbor at the end of neighbor list\n");
	tet->neighborTetrahedra[tet->countNeighborTetrahedra] = neighborTet;
	//printf("increase neighbor counter\n");
	tet->countNeighborTetrahedra++;	
	//printf("done!\n");

	return 1;
}
/*****************************************************************************/


template <int Dimensions,class VertexType>
int Simplex<Dimensions,VertexType>::SimplexContainsFace( VertexType** face)
{
	VertexType** tet = this->vertices;
	int i, j = 0;
	
	for( i=0; i<Dimensions && j!=Dimensions+1; i++)
		for( j=0; j<Dimensions+1 && face[i]!=tet[j]; j++);

	if( i==Dimensions && j!=Dimensions+1){
		//fprintf( stderr, "( %i, %i, %i, %i) contains ( %i, %i, %i)\n", tet[0]->index, tet[1]->index, tet[2]->index, tet[3]->index, face[0]->index, face[1]->index, face[2]->index);
		return 1;
	}
	else
		return 0;
}
/*****************************************************************************/


template <int Dimensions, class VertexType>
inline VertexType*& Simplex<Dimensions,VertexType>::getVertex( int i)
{
	return vertices[i];
}
/*****************************************************************************/


template <int Dimensions, class VertexType>
inline void Simplex<Dimensions,VertexType>::setVertex( VertexType* &vertex, int i)
{
	vertices[i] = vertex;
}
/*****************************************************************************/


template <int Dimensions,class VertexType>
inline void Simplex<Dimensions,VertexType>::setIndex( int index)
{
	this->index = index;
}
/*****************************************************************************/


template <int Dimensions, class VertexType, class SimplexType>
void Triangulation<Dimensions,VertexType,SimplexType>::setSimplexNeighbors()
{
	Triangulation<Dimensions,VertexType,SimplexType>* voronoiDiagram = this;
	int i, j;

	for( i=0; i<voronoiDiagram->countTetrahedra-1; i++){
		// find neighbors of Simplex i
		
		for( j=i+1; j<voronoiDiagram->countTetrahedra; j++){
			// check if Simplex j is a neighbor of Simplex i

			if( voronoiDiagram->tetrahedra[i]->countCommonPointOfTetrahedra( voronoiDiagram->tetrahedra[j]) == Dimensions){
				voronoiDiagram->tetrahedra[i]->addSimplexNeighbor( voronoiDiagram->tetrahedra[j]);
				voronoiDiagram->tetrahedra[j]->addSimplexNeighbor( voronoiDiagram->tetrahedra[i]);
			}
		}
	}

	// TEST: check if all neighborhood relations are reversible
	/*int k;
	for( i=0; i<voronoiDiagram->countTetrahedra; i++){
		for( j=0; j<voronoiDiagram->tetrahedra[i]->countNeighborTetrahedra; j++){
			for( k=0; k<voronoiDiagram->tetrahedra[i]->neighborTetrahedra[j]->countNeighborTetrahedra &&
				voronoiDiagram->tetrahedra[i] != voronoiDiagram->tetrahedra[i]->neighborTetrahedra[j]->neighborTetrahedra[k]; k++);
			
			if( k==voronoiDiagram->tetrahedra[i]->neighborTetrahedra[j]->countNeighborTetrahedra){
				fprintf( stderr, "ERROR: neighborship relation between tets %p -> %p!!\n", voronoiDiagram->tetrahedra[i], voronoiDiagram->tetrahedra[i]->neighborTetrahedra[j]->neighborTetrahedra[k]);
				exit( 0);
			}
		}
	}*/
}
/*****************************************************************************/


template <int Dimensions, class VertexType>
void Simplex<Dimensions,VertexType>::printSimplexNeighbors()
{
	Simplex<Dimensions,VertexType>* s = this;
	int j;
	
	/*for( i=0; i<voronoiDiagram->countTetrahedra; i++){
		fprintf( stderr, "Simplex %p has %i neighbors: ", voronoiDiagram->tetrahedra[i], voronoiDiagram->tetrahedra[i]->countNeighborTetrahedra);
		for( j=0; j<voronoiDiagram->tetrahedra[i]->countNeighborTetrahedra; j++){
			fprintf( stderr, "%p ", voronoiDiagram->tetrahedra[i]->neighborTetrahedra[j]);
		}
		fprintf( stderr, "\n");
	}*/
	fprintf( stderr, "Simplex %p has %i neighbors: ", s, s->countNeighborTetrahedra);
	for( j=0; j<s->countNeighborTetrahedra; j++){
		fprintf( stderr, "%p ", s->neighborTetrahedra[j]);
	}
	fprintf( stderr, "\n");		
}
/*****************************************************************************/


template <int Dimensions, class VertexType, class SimplexType>
void Triangulation<Dimensions,VertexType,SimplexType>::printTetrahedraNeighbors()
{
	int i;
	
	for( i=0; i<countTetrahedra; i++){
		tetrahedra[i]->printSimplexNeighbors( );
	}		
}
/*****************************************************************************/


template < int Dimensions >
int Vertex<Dimensions>::removeNeighbor( Vertex<Dimensions>* neighborCell)
{
	Vertex<Dimensions>* voronoiCell = this;

	int i;

	//fprintf( stderr, "INFO: removeNeighbor(): remove %i as neighbor of %i!\n", neighborCell->index, voronoiCell->index);
	
	for( i=0; i<voronoiCell->countNeighbors; i++){
		if( voronoiCell->neighbors[i] == neighborCell){
			voronoiCell->countNeighbors--;
			voronoiCell->neighbors[i] = voronoiCell->neighbors[voronoiCell->countNeighbors];
			voronoiCell->neighbors = (Vertex<Dimensions>**) realloc( voronoiCell->neighbors, voronoiCell->countNeighbors*sizeof(Vertex<Dimensions>*));
			assert(voronoiCell->countNeighbors==0 || voronoiCell->neighbors);
			return 1;
		}	
	}

	//fprintf( stderr, "ERROR: removeNeighbor() can not find the neighbor to delete!\n");
	//exit( 0);
	return 0;
}
/*****************************************************************************/


template <int Dimensions>
int Vertex<Dimensions>::addNeighbor( Vertex<Dimensions>* neighborCell)
{
	Vertex<Dimensions>* voronoiCell = this;
	int i;
	for( i=0; i<voronoiCell->countNeighbors; i++){
		if( voronoiCell->neighbors[i] == neighborCell){
	//		fprintf( stderr, "ERROR: addNeighbor() found neighbor being already included!\n");
			//exit( 0);
			return 0;
		}	
	}
	

	//printf("realloc: %i -> %i\n", voronoiCell->countNeighbors, voronoiCell->countNeighbors+1);
	//if( voronoiCell->countNeighbors != 0)
		voronoiCell->neighbors = (Vertex<Dimensions>**) realloc( voronoiCell->neighbors, (voronoiCell->countNeighbors+1)*sizeof(Vertex<Dimensions>*));
		assert(voronoiCell->neighbors);
	//else
	//	voronoiCell->neighbors = (Vertex<Dimensions>**) malloc( sizeof());
	//printf("set new neighbor at the end of neighbor list\n");
	voronoiCell->neighbors[voronoiCell->countNeighbors] = neighborCell;
	//printf("increase neighbor counter\n");
	voronoiCell->countNeighbors++;
	//printf("done!\n");

	return 1;
}
/*****************************************************************************/



template <int Dimensions, class VertexType, class SimplexType>
void Triangulation<Dimensions,VertexType,SimplexType>::detriangulate( VertexType* removedVoronoiCell)
{
	Triangulation<Dimensions,VertexType,SimplexType>* voronoiDiagram = this;

	int NR_TETRAHEDRON_POINTS = Dimensions+1;
	int NR_FACE_POINTS = Dimensions;
	//int DIMENSIONS = Dimensions;

	int i, j, k, l, m;

	/*if(removedVoronoiCell->index==4674){
		double volume = 0.;
		for( int i=0; i<voronoiDiagram->countTetrahedra; i++){
			volume += getVolume(
					voronoiDiagram->tetrahedra[i]->vertices[0]->position,
					voronoiDiagram->tetrahedra[i]->vertices[1]->position,
					voronoiDiagram->tetrahedra[i]->vertices[2]->position,
					voronoiDiagram->tetrahedra[i]->vertices[3]->position);
		}
		fprintf(stderr, "volume of all tets: %.10lf\n", volume);
	/ *FILE* fp = fopen("removedCell.dat","w+");
	fprintf(fp, "\n%lf %lf %lf %i\n", removedVoronoiCell->position[0], removedVoronoiCell->position[1], removedVoronoiCell->position[2], removedVoronoiCell->index);
	for(int i=0; i<removedVoronoiCell->countNeighborCells; i++){
		fprintf(fp, "%lf %lf %lf %i\n", removedVoronoiCell->neighborCells[i]->position[0], removedVoronoiCell->neighborCells[i]->position[1], removedVoronoiCell->neighborCells[i]->position[2], removedVoronoiCell->neighborCells[i]->index);
	}
	fclose(fp);
	checkDelaunayCondition(voronoiDiagram, 0, 0);* /
	}*/
	// tetrahedra
	int ti = voronoiDiagram->countTetrahedra;
	int timax = voronoiDiagram->maxTetrahedra;
	Simplex<Dimensions> **tetrahedra = voronoiDiagram->tetrahedra;

	// faces
	int fi = 0;			// face counter
	int fimax = voronoiDiagram->maxFaces;			// maximal number of faces in list
	VertexType ***faces = voronoiDiagram->faces;		// list of faces

	// points
	int pi = 0;			// face counter
	int api = 0;			// face counter
	int pimax = INIT_FACES_ARRAY_SIZE;	// maximal number of faces in list
	VertexType **points;	// list of faces
	points = (VertexType**) calloc( sizeof(VertexType*), INIT_FACES_ARRAY_SIZE);
	assert(points);
	int *pointsInvolved;
	pointsInvolved = (int*) calloc( sizeof(int), INIT_FACES_ARRAY_SIZE);
	assert(pointsInvolved);
	VertexType **allPoints;	// list of faces



	// DELETE POINT
	//voronoiDiagram->countVoronoiCells--;



	// SEARCH INVOLVED TETRAHEDRA

#if TRACK_DELAUNAY_REGION_NEIGHBORHOOD

	SimplexType **correspondingTetrahedron = (SimplexType**) calloc( sizeof(SimplexType*), fimax);

	//fprintf( stderr, "look for first tetrahedra containing deleted point\n");

	SimplexType *actualTet = voronoiDiagram->getSimplexContainingPoint( removedVoronoiCell);
	//fprintf( stderr, "END!!!!\n");

	// look for all tetrahedra deleted by the point insertion	
	int maxTetStackLength = TETRAHEDRA_ARRAY_EXTENSION_SIZE;		
	SimplexType **tetStack = (SimplexType**) calloc( sizeof(SimplexType*), TETRAHEDRA_ARRAY_EXTENSION_SIZE);
	assert(tetStack);
	int maxUsedTetStackLength = TETRAHEDRA_ARRAY_EXTENSION_SIZE;		
	SimplexType **usedTetStack = (SimplexType**) calloc( sizeof(SimplexType*), TETRAHEDRA_ARRAY_EXTENSION_SIZE);
	assert(usedTetStack);

	tetStack[0] = actualTet; //tetrahedra[0];
	int tetStackLength = 1;	
	int usedTetStackLength = 0;	


	// TEST:
	/*int test_count_tetrahedra_deleted = 0;
	for( j=0; j<ti; j++){
		if( getDistanceOfPointToCircumsphereCenter( newVoronoiCell, tetrahedra[j]) < getCircumsphereRadius( tetrahedra[j]))
			test_count_tetrahedra_deleted++;
	}*/
	int countIts = 0;
	// END TEST

	//fprintf( stderr, "look for all tetrahedra deleted by the point deletion\n");

	//fprintf( stderr, "TEST1!!!!\n");

	do{
		// TAKE NEXT TETRAHEDRON FROM STACK
		tetStackLength--;
		actualTet = tetStack[tetStackLength];
		if( usedTetStackLength == maxUsedTetStackLength){
			maxUsedTetStackLength += TETRAHEDRA_ARRAY_EXTENSION_SIZE;
			usedTetStack = (SimplexType**) realloc( usedTetStack, sizeof(SimplexType*) * maxUsedTetStackLength);
			assert(usedTetStack);
		}
		usedTetStack[usedTetStackLength] = actualTet;
		usedTetStackLength++;

		//fprintf( stderr, "INFO: Delete tet %i (%i, %i, %i, %i)!\n", actualTet->index, actualTet->vertices[0]->index, actualTet->vertices[1]->index, actualTet->vertices[2]->index, actualTet->vertices[3]->index);

		// ADD NEIGHBOR TETRAHEDRA TO STACK
		//int temp_tetStackLength = tetStackLength;
		for( j=0; j<actualTet->countNeighborTetrahedra /*&& temp_tetStackLength == tetStackLength*/; j++){
			if( tetStackLength + actualTet->countNeighborTetrahedra < maxTetStackLength){
				maxTetStackLength += TETRAHEDRA_ARRAY_EXTENSION_SIZE;
				tetStack = (SimplexType**) realloc( tetStack, sizeof(SimplexType*) * maxTetStackLength);
				assert(tetStack);
			}
			//distanceToCircumsphereCenter = getDistanceOfPointToCircumsphereCenter( newVoronoiCell, actualTet->neighborTetrahedra[j]);
			for( k=0; k<usedTetStackLength && usedTetStack[k]!=actualTet->neighborTetrahedra[j]; k++);
			for( l=0; l<tetStackLength     && tetStack[l]!=actualTet->neighborTetrahedra[j];     l++);
			// circumsphere contains new point ?
			if( k==usedTetStackLength && l==tetStackLength && 
			    ( actualTet->neighborTetrahedra[j]->vertices[0]==removedVoronoiCell
			    		|| actualTet->neighborTetrahedra[j]->vertices[1]==removedVoronoiCell
			    		|| actualTet->neighborTetrahedra[j]->vertices[2]==removedVoronoiCell
			    		|| (Dimensions==3 && actualTet->neighborTetrahedra[j]->vertices[3]==removedVoronoiCell)
			    		)){
				tetStack[tetStackLength]=actualTet->neighborTetrahedra[j];
				tetStackLength++;
			}
		}

		// SET j TO POSITION OF VERTEX EQUAL TO REMOVED POINT
		for( j=0; j<NR_TETRAHEDRON_POINTS && actualTet->vertices[j] != removedVoronoiCell; j++);
		j++;

#else //TRACK_DELAUNAY_REGION_NEIGHBORHOOD

	// for all tetrahedra
	for( i=0; i<ti; i++)
	{
		Tetrahedron *actualTet = tetrahedra[i];

		// for all vertices
		int containsRemovedPoint = 0;
		for( j=0; j<NR_TETRAHEDRON_POINTS && !containsRemovedPoint; j++){
			// tetrehedron contains removed point 
			if( actualTet->vertices[j] == removedVoronoiCell){
				//fprintf( stderr, "Tetrahedra %i contains removed point %i\n", i, removedVoronoiCell->index);
				containsRemovedPoint++;
			}
		}

		if( containsRemovedPoint){
			//fprintf( stderr, "INFO: Delete tet %i (%i, %i, %i, %i)!\n", actualTet->index, actualTet->vertices[0]->index, actualTet->vertices[1]->index, actualTet->vertices[2]->index, actualTet->vertices[3]->index);

#endif //TRACK_DELAUNAY_REGION_NEIGHBORHOOD

			// ADD FACES TO LIST
			if( fi+NR_TETRAHEDRON_POINTS >= fimax){
				//fprintf( stderr, "INFO: realloc: %i -> %i!\n", fimax, fimax + FACES_ARRAY_EXTENSION_SIZE);
				faces = (VertexType***) realloc( faces, sizeof(VertexType**) * (fimax + FACES_ARRAY_EXTENSION_SIZE));
				assert(faces);
				for( l=fimax; l<fimax + FACES_ARRAY_EXTENSION_SIZE; l++){
					faces[l] = (VertexType**) calloc( sizeof(VertexType*), NR_FACE_POINTS);
					assert(faces[l]);
				}
				fimax += FACES_ARRAY_EXTENSION_SIZE;
			}
			//k=j;
			//fprintf( stderr, "j=%i, Add faces %i (", j, fi);
			for( l=0; l<NR_FACE_POINTS; l++){
				k=0;
				for( k=0; k<pi && points[k]!=actualTet->vertices[(j+l)%(NR_TETRAHEDRON_POINTS)]; k++);
				//if( pi==0 || points[k]!=actualTet->vertices[(j+l)%(NR_TETRAHEDRON_POINTS)]){
				if( k==pi){
					// realloc
					if( pi==pimax){
						pimax += FACES_ARRAY_EXTENSION_SIZE;
						points = (VertexType**) realloc( points, sizeof(VertexType*) * pimax);
						assert(points);
						pointsInvolved = (int*)    realloc( pointsInvolved,    sizeof(int) * pimax);
						assert(pointsInvolved);
					}
					//j=
					points[pi] = actualTet->vertices[(j+l)%(NR_TETRAHEDRON_POINTS)];
					pointsInvolved[pi] = 1;
					pi++;
				}else
					pointsInvolved[k]++;
				faces[fi][l%NR_FACE_POINTS]=actualTet->vertices[(j+l)%(NR_TETRAHEDRON_POINTS)];
				//fprintf( stderr, " %i", faces[fi][l]->index);

				//pointsInvolved[k]++;


				// REMOVE NEIGHBORSHIP RELATIONS TO REMOVED POINT
				removedVoronoiCell->removeNeighbor( faces[fi][l]);
				faces[fi][l]->removeNeighbor( removedVoronoiCell);
			}
			//fprintf( stderr, ")\n");
#if TRACK_DELAUNAY_REGION_NEIGHBORHOOD

			// DELETE TETRAHEDRA NEIGHBORSHIP RELATIONS

			for( l=actualTet->countNeighborTetrahedra-1; l>=0; l--){
				if( actualTet->neighborTetrahedra[l]->SimplexContainsFace( faces[fi]))
					correspondingTetrahedron[fi] = actualTet->neighborTetrahedra[l];
				actualTet->neighborTetrahedra[l]->removeSimplexNeighbor( actualTet);
				actualTet->removeSimplexNeighbor( actualTet->neighborTetrahedra[l]);
			}

#endif	
			fi++;


			// DELETE TETRAHEDRON
				
			//fprintf( stderr, "INFO: delete tetrahedron!\n");
			//actualTet->countNeighborTetrahedra = 0;

			ti--;
			int tempIndex = actualTet->index;
			SimplexType * temp = tetrahedra[tempIndex];
	
			tetrahedra[tempIndex] = tetrahedra[ti];
			tetrahedra[tempIndex]->index = tempIndex;
	
			tetrahedra[ti] = temp;
			tetrahedra[ti]->index = ti;
				
			//printPoints( "Points: ", points, pi );


#if TRACK_DELAUNAY_REGION_NEIGHBORHOOD

	}while( tetStackLength != 0);

	free( tetStack);
	free( usedTetStack);

#else // TRACK_DELAUNAY_REGION_NEIGHBORHOOD
			
			if( i>=tempIndex)
				i--;
		}
	}

#endif // TRACK_DELAUNAY_REGION_NEIGHBORHOOD

	api = pi;
	allPoints = (VertexType**) calloc( sizeof(VertexType*), api);
	assert(allPoints);
	for( i=0; i<pi; i++)
		allPoints[i] = points[i];
		//pointsInvolved[i] = 2;

	//fprintf( stderr, "TEST!!!!\n");
	//printPoints( "Points: ", allPoints, api );




	// CREATE NEW TETRAHEDRA FROM FACES
		
	int countAddedTetrahedra = 0;
	//int last_minCountValidTetrahedra = 1;
	//int minCountValidTetrahedra = fi;
	//int countNotAddedFacesInAChain = 0;

	for( i=0; pi!=0; i = (fi==0 ? 0 : (i+1)%fi)){ // for all faces

		// STOP ENDLESS LOOP
		if( countIts==2000){
			fprintf( stderr, "WARNING: Apparently there is a Problem in removeVoronoiCell()\n");
			//fprintf( stderr, "Algorithm imprisoned in endless loop while trying to remove point %i!\n", removedVoronoiCell->index);
			fprintf( stderr, "Print Point coordinates to file\n");
			/*FILE *fp;
			fp = fopen( "errorPointSet.nodes", "w+");
			fprintf( fp, "%i %i\n", voronoiDiagram->countVoronoiCells, DIMENSIONS);
			int index = 0;
			for( int v=0; v<voronoiDiagram->countVoronoiCells; v++)
			if(voronoiDiagram->vertices[v]->countNeighborCells>0 || voronoiDiagram->vertices[v] == removedVoronoiCell){
				if( voronoiDiagram->vertices[v] == removedVoronoiCell)
					fprintf( stderr, "Algorithm imprisoned in endless loop while trying to remove point %i!\n", index);

				for( int p=0; p<DIMENSIONS; p++)
					fprintf( fp, "%.20lf ", voronoiDiagram->vertices[v]->position[p]);
				fprintf( fp, "\n");
				index++;
			}
			fclose( fp);*/

			FILE *pFile = fopen( "errorPointSet.nodes", "w+");
			int index = 0;
			for( int v=0; v<voronoiDiagram->countVertices; v++)
			if(voronoiDiagram->vertices[v]->numberOfNeighbors()>0 || voronoiDiagram->vertices[v] == removedVoronoiCell)
			{
				if( voronoiDiagram->vertices[v] == removedVoronoiCell){
					fprintf( stderr, "Algorithm imprisoned in endless loop while trying to remove point %i (in output file: %i)!\n", removedVoronoiCell->getIndex(), index);
					for( int p=0; p<Dimensions; p++)
						fprintf( stderr, "%lf ", voronoiDiagram->vertices[v]->pos()[p]);
					fprintf( stderr, "\n");
				}
				fwrite( voronoiDiagram->vertices[v]->pos() , sizeof(voronoiDiagram->vertices[v]->pos()), 1, pFile );
				index++;
			}
			fclose(pFile);

			exit( 0);
		}
		countIts++;
		// STOP ENDLESS LOOP


		if( ti == timax){
			//fprintf( stderr, "INFO: realloc: %i -> %i!\n", timax, timax + TETRAHEDRA_ARRAY_EXTENSION_SIZE);
			tetrahedra = (SimplexType **) realloc( tetrahedra, sizeof(SimplexType*)*(timax + TETRAHEDRA_ARRAY_EXTENSION_SIZE));
			assert(tetrahedra);
			for( k=timax; k<timax + TETRAHEDRA_ARRAY_EXTENSION_SIZE;k++)
				tetrahedra[k] = new SimplexType();
			timax += TETRAHEDRA_ARRAY_EXTENSION_SIZE;
		}
		//fprintf( stderr, "Complete new tetrahedron %i: (", ti);
		// set face as points of new tetrahedron
		for( k=0; k<NR_FACE_POINTS; k++){
			//fprintf( stderr, " %i", faces[i][k]->index);	
			tetrahedra[ti]->vertices[k] = faces[i][k];
		}
		//fprintf( stderr, " ?)\n");	

		// look for point to complete tetrahedron
		int tetrahedronValid = 0;
		int countValidTetrahedra = 0;
		VertexType * validPoint = NULL;


		for( j=0; j<pi /*TEST && !tetrahedronValid*/; j++){ // for all points
			// check if point is not contained yet by face
			for( k=0; k<NR_FACE_POINTS && points[j]!=faces[i][k]; k++);
			if( k==NR_FACE_POINTS){
				//fprintf( stderr, "point %i is candidate to form a tetrahedron with face %i\n", points[j]->index, i);

				// set point as point of new tetrahedron
				tetrahedra[ti]->vertices[k] = points[j];
				tetrahedra[ti]->circumsphereInitialized = 0;
				
				// check validy of new tetrahedron
				tetrahedronValid = 1;

				// removed point lies in the circumsphere?
				if( removedVoronoiCell->getDistanceToCircumsphere( tetrahedra[ti]) < 0.){
					// other points lie in the circumsphere?
					/*for( k=0; k<pi; k++){ // for all points
						// check if point is contained tetrahedron
						for( l=0; l<NR_TETRAHEDRON_POINTS && points[k]!=tetrahedra[ti]->vertices[l]; l++);
						if( l==NR_TETRAHEDRON_POINTS){
							// is point in the new tetrahedrons circumsphere ???
							if( getDistanceOfPointToCircumsphere( points[k], tetrahedra[ti]) < 0.){
								tetrahedronValid = 0;
							}					
						}
					}*/
					for( k=0; k<api; k++){ // for all points
						// check if point is contained tetrahedron
						for( l=0; l<NR_TETRAHEDRON_POINTS && allPoints[k]!=tetrahedra[ti]->vertices[l]; l++);
						if( l==NR_TETRAHEDRON_POINTS){
							// is point in the new tetrahedrons circumsphere ???
							if( allPoints[k]->getDistanceToCircumsphere( tetrahedra[ti]) < 0.){
								tetrahedronValid = 0;
							}					
						}
					}

//#if TRACK_DELAUNAY_REGION_NEIGHBORHOOD
					// tet exists already?
					if( correspondingTetrahedron[i]!=NULL){
						for( k=0; k<NR_TETRAHEDRON_POINTS && correspondingTetrahedron[i]->vertices[k]!=points[j]; k++);
						if( k<NR_TETRAHEDRON_POINTS)
							tetrahedronValid = 0;
					}
//#endif
					
				}else{
					//fprintf( stderr, "WARNING: new tetrahedron %i is not valid: removed point is not within the circumsphere\n", ti);
					//TEST 
					tetrahedronValid = 0;
				}

				if( tetrahedronValid){
					countValidTetrahedra++;
					validPoint = points[j];
				}
				//printPoints( "test1 points:      ", allPoints, api );
			}
		}

		//fprintf( stderr, "INFO: fi=%i, i=%i, ti=%i:  countValidTetrahedra=%i   < last_minCountValidTetrahedra=%i ?\n", fi, i, ti, countValidTetrahedra, last_minCountValidTetrahedra);

		/*if( countValidTetrahedra > 1){
			fprintf( stderr, "INFO: countValidTetrahedra = %i\n", countValidTetrahedra);
			for( k=0; k<NR_FACE_POINTS; k++){
				fprintf( stderr, " %i", faces[i][k]->index);
				//tetrahedra[ti]->vertices[k] = faces[i][k];
			}
			fprintf( stderr, "\n");
			if( countIts>1000)
				countValidTetrahedra=1;
		}
		if( countValidTetrahedra < 1){
			fprintf( stderr, "INFO: countValidTetrahedra = %i\n", countValidTetrahedra);
			for( k=0; k<NR_FACE_POINTS; k++){
				fprintf( stderr, " %i", faces[i][k]->index);
				//tetrahedra[ti]->vertices[k] = faces[i][k];
			}
			fprintf( stderr, "\n");
			exit( 0);
		}*/

		if( countValidTetrahedra == 1){
		//if( countValidTetrahedra && countValidTetrahedra == last_minCountValidTetrahedra){
			tetrahedra[ti]->vertices[Dimensions] = validPoint;
			tetrahedra[ti]->circumsphereInitialized = 0;

			countAddedTetrahedra++;

			// print new tetrahedron
			//fprintf( stderr, "new tetrahedron %i is (", ti);
			for( k=0; k<NR_FACE_POINTS; k++){
				tetrahedra[ti]->vertices[k]->addNeighbor( validPoint);
				validPoint->addNeighbor( tetrahedra[ti]->vertices[k]);
				//fprintf( stderr, " %i", tetrahedra[ti]->vertices[k]->index);	
			}
			//fprintf( stderr, ")\n");	
			

			// DELETE USED FACE FROM LIST
			fi--;
			//fprintf( stderr, "Remove Face %i: (", i);	
			for( k=0; k<NR_FACE_POINTS; k++){
				//fprintf( stderr, " %i", faces[i][k]->index);	
				for( l=0; faces[i][k]!=points[l]; l++);
				pointsInvolved[l]--;
				faces[i][k] = faces[fi][k];
			}

#if TRACK_DELAUNAY_REGION_NEIGHBORHOOD

			if(correspondingTetrahedron[i]!=NULL){
				tetrahedra[ti]->addOrReplaceSimplexNeighbor( correspondingTetrahedron[i]);
				correspondingTetrahedron[i]->addOrReplaceSimplexNeighbor( tetrahedra[ti]);
			}
			correspondingTetrahedron[i] = correspondingTetrahedron[fi];

#endif
			
			//fprintf( stderr, ")\n");	
			i--;
			// ADD OTHER FACES TO LIST
			for( j=1; j<NR_TETRAHEDRON_POINTS; j++){ // for all new faces
				int commonVertices = 0; //TODO: optimize: replace faceAlreadyInList by commonVertices!!!
				for( k=0; commonVertices<NR_FACE_POINTS && k<fi; k++){ // check all faces for a duplicate
					commonVertices = 0;
					for( l=0; l<NR_FACE_POINTS; l++){ // for all vertices of new face
						//commonVertices = 0;
						for( m=0; m<NR_FACE_POINTS; m++){ // for all vertices of face in list
							if( tetrahedra[ti]->vertices[(j+l)%NR_TETRAHEDRON_POINTS] == faces[k][m])
								commonVertices++;
						}
					}
				}
				// face already in the list ?
				if( commonVertices>=NR_FACE_POINTS){
					// already in the list => delete face
					k--;
					fi--;
					//fprintf( stderr, "Remove Face %i: (", k);	
					for( l=0; l<NR_FACE_POINTS; l++){
						//fprintf( stderr, " %i", faces[k][l]->index);
						// deleted face has common points with previous deleted faces?	
						for( m=0; faces[k][l]!=points[m]; m++);
						/*if(m>=pi){
							fprintf( stderr, "\n ERROR: face point not in point list!!!\n");	
							exit( 0);
						}*/
						pointsInvolved[m]--;
						if( pointsInvolved[m] == 0){
							// delete from point list
							pi--;
							points[m] = points[pi];
							pointsInvolved[m] = pointsInvolved[pi];
								
						}

						faces[k][l]=faces[fi][l];
					}
					//fprintf( stderr, ")\n");

#if TRACK_DELAUNAY_REGION_NEIGHBORHOOD

					if(correspondingTetrahedron[k]!=NULL){
						tetrahedra[ti]->addOrReplaceSimplexNeighbor( correspondingTetrahedron[k]);
						correspondingTetrahedron[k]->addOrReplaceSimplexNeighbor( tetrahedra[ti]);
					}
					correspondingTetrahedron[k] = correspondingTetrahedron[fi];

#endif

				}else{
					// not yet in the list => add face	
					//fprintf( stderr, "Add Face %i: (", fi);	
					for( l=0; l<NR_FACE_POINTS; l++){
						faces[fi][l]=tetrahedra[ti]->vertices[(j+l)%NR_TETRAHEDRON_POINTS];
						//fprintf( stderr, " %i", faces[fi][l]->index);	
						for( m=0; faces[fi][l]!=points[m]; m++);
						/*if(m>=pi){
							fprintf( stderr, "\n ERROR: face point not in point list!!!\n");	
							exit( 0);
						}*/
						pointsInvolved[m]++;
					}
					//fprintf( stderr, ")\n");

#if TRACK_DELAUNAY_REGION_NEIGHBORHOOD

					correspondingTetrahedron[fi] = tetrahedra[ti];

#endif

					fi++;
				}
	
			}
			//fprintf( stderr, "\n");	
			tetrahedra[ti]->index = ti;
			ti++;

			//countNotAddedFacesInAChain = 0;				
			//last_minCountValidTetrahedra = 1;
		}
		/*else{
			countNotAddedFacesInAChain++;
			if( countValidTetrahedra == 0){
				fprintf( stderr, "ERROR: countValidTetrahedra == 0!!!\n");
				exit( 0);
			}
		}
		if( fi>0 && countNotAddedFacesInAChain == fi){
			last_minCountValidTetrahedra = minCountValidTetrahedra;
			minCountValidTetrahedra = fi;
			countNotAddedFacesInAChain = 0;
			fprintf( stderr, "\n WARNING: More than one posibility to include new tetrahedra: %i!\n", last_minCountValidTetrahedra);
			fprintf( stderr, "INFO: ti=%i, fi=%i\n", ti, fi);
			printPoints( "all points:  ", allPoints, api );
			printPoints( "points left: ", points, pi );
			voronoiDiagram->countTetrahedra = ti;
			//checkDelaunayCondition( voronoiDiagram, NULL, 0);
			//exit( 0);
			//occurenceOfAmbiguity++;
		}*/

	}


	// DELETE POINT
	//free( removedVoronoiCell);
	

	// FREE ALLOCATED MEMORY

	//deleteMatrix( A, NR_TETRAHEDRON_POINTS);
	free( points);
	free( pointsInvolved);
	free( allPoints);

#if TRACK_DELAUNAY_REGION_NEIGHBORHOOD

	free( correspondingTetrahedron);

#endif	

	// RESET VALUES OF VORONOI DIAGRAM

	voronoiDiagram->countFaces = fi;			// face counter
	voronoiDiagram->maxFaces = fimax;			// maximal number of faces in list
	voronoiDiagram->faces = faces;		// list of faces

	voronoiDiagram->countTetrahedra = ti; 
	voronoiDiagram->maxTetrahedra = timax;
	voronoiDiagram->tetrahedra = tetrahedra;



	// FINAL TESTING
	/*if( occurenceOfAmbiguity)
		checkDelaunayCondition( voronoiDiagram, NULL, 0);
	*/
	
}
/*****************************************************************************/


template <int Dimensions, class VertexType, class SimplexType>
inline void Triangulation<Dimensions,VertexType,SimplexType>::add( VertexType* newVoronoiCell)
{
	newVoronoiCell->setIndex( countVertices);

	// ADD NEW POINT TO POINT ARRAY
	if( countVertices==maxVertices){
		maxVertices += POINTS_ARRAY_EXTENSION_SIZE;
		vertices = (VertexType **) realloc( vertices, sizeof(VertexType*) * maxVertices);
		assert(vertices);
	}
	vertices[countVertices] = newVoronoiCell;
	countVertices++;
}
/*****************************************************************************/


template <int Dimensions, class VertexType,class SimplexType>
void Triangulation<Dimensions,VertexType,SimplexType>::remove( VertexType* v)
{
	countVertices--;
	vertices[v->index] = vertices[countVertices];
}
/*****************************************************************************/


template <int Dimensions, class VertexType, class SimplexType>
void Triangulation<Dimensions,VertexType,SimplexType>::removeAndDetriangulate( VertexType* v)
{
	this->detriangulate( v);
	countVertices--;
	vertices[v->getIndex()] = vertices[countVertices];
	vertices[v->getIndex()]->setIndex( v->getIndex());
}
/*****************************************************************************/


template <int Dimensions, class VertexType, class SimplexType>
void Triangulation<Dimensions,VertexType,SimplexType>::triangulate( VertexType* newVoronoiCell)
{
	Triangulation<Dimensions,VertexType,SimplexType>* voronoiDiagram = this;

	int j, k, l;
	int NR_Simplex_POINTS=Dimensions+1;
	int NR_Face_POINTS=Dimensions;


	//double **A = newMatrix( Dimensions + 1, Dimensions + 1);

	int count_faces_added = 0;
	int count_common_faces_deleted = 0;
	int count_tetrahedra_deleted = 0;
	int count_tetrahedra_added = 0;
	//int npr = 0;

	// faces
	int fi = 0;			// face counter
	int fimax = voronoiDiagram->maxFaces;			// maximal number of faces in list
	VertexType*** faces = voronoiDiagram->faces;		// list of faces


	// tetrahedra
	int ti = voronoiDiagram->countTetrahedra;
	int timax = voronoiDiagram->maxTetrahedra;
	SimplexType** tetrahedra = voronoiDiagram->tetrahedra;

	//newVoronoiCell->index = voronoiDiagram->countVertices;


	// DELAUNAY TRIANGULATION (WATSON)

#if TRACK_DELAUNAY_REGION_NEIGHBORHOOD

	SimplexType **correspondingSimplex = (SimplexType**) calloc( sizeof(SimplexType*), fimax);
	assert(correspondingSimplex);

	SimplexType *actualTet = voronoiDiagram->getSimplexContainingPointInCircumSphere( newVoronoiCell);

	// look for all tetrahedra deleted by the point insertion
	int maxTetStackLength = TETRAHEDRA_ARRAY_EXTENSION_SIZE;
	SimplexType **tetStack = (SimplexType**) calloc( sizeof(SimplexType*), TETRAHEDRA_ARRAY_EXTENSION_SIZE);
	assert(tetStack);
	int maxUsedTetStackLength = TETRAHEDRA_ARRAY_EXTENSION_SIZE;
	SimplexType **usedTetStack = (SimplexType**) calloc( sizeof(SimplexType*), TETRAHEDRA_ARRAY_EXTENSION_SIZE);
	assert(usedTetStack);

	tetStack[0] = actualTet; //tetrahedra[0];
	int tetStackLength = 1;
	int usedTetStackLength = 0;


	// TEST:
	/*int test_count_tetrahedra_deleted = 0;
	for( j=0; j<ti; j++){
		if( getDistanceOfPointToCircumsphereCenter( newVoronoiCell, tetrahedra[j]) < getCircumsphereRadius( tetrahedra[j]))
			test_count_tetrahedra_deleted++;
	}*/
	// END TEST

	//fprintf( stderr, "look for all tetrahedra deleted by the point insertion\n");


	do{

		// END TEST

		// TAKE NEXT Simplex FROM STACK
		tetStackLength--;
		actualTet = tetStack[tetStackLength];
		if( usedTetStackLength == maxUsedTetStackLength){
			maxUsedTetStackLength += TETRAHEDRA_ARRAY_EXTENSION_SIZE;
			usedTetStack = (SimplexType**) realloc( usedTetStack, sizeof(SimplexType*) * maxUsedTetStackLength);
			assert(usedTetStack);
		}
		usedTetStack[usedTetStackLength] = actualTet;
		usedTetStackLength++;


		// TEST
		/*int i;
		for( i=0; i<tetStackLength; i++)
		{
			for( j=0; j<actualTet->countNeighborTetrahedra; j++){
				if( actualTet->neighborTetrahedra[j] == NULL)
					exit( 0);

				if( actualTet->neighborTetrahedra[j]->circumsphereInitialized == 1){
					if( actualTet->neighborTetrahedra[j]->circumsphere[0] == 0.123456789)
						exit( 0);
				}
			}

			if( tetStack[i]->circumsphereInitialized == 1){
				if( tetStack[i]->circumsphere[0] == 0.123456789)
					exit( 0);
			}
		}*/


		//fprintf( stderr, "INFO: Delete tet %i (%i, %i, %i, %i)!\n", actualTet->index, actualTet->vertices[0]->index, actualTet->vertices[1]->index, actualTet->vertices[2]->index, actualTet->vertices[3]->index);

		// ADD NEIGHBOR TETRAHEDRA TO STACK
		//int temp_tetStackLength = tetStackLength;
		for( j=0; j<actualTet->countNeighborTetrahedra /*&& temp_tetStackLength == tetStackLength*/; j++){
			if( tetStackLength + actualTet->countNeighborTetrahedra < maxTetStackLength){
				maxTetStackLength += TETRAHEDRA_ARRAY_EXTENSION_SIZE;
				tetStack = (SimplexType**) realloc( tetStack, sizeof(SimplexType*) * maxTetStackLength);
				assert(tetStack);
			}
			//distanceToCircumsphereCenter = getDistanceOfPointToCircumsphereCenter( newVoronoiCell, actualTet->neighborTetrahedra[j]);
			for( k=0; k<usedTetStackLength && usedTetStack[k]!=actualTet->neighborTetrahedra[j]; k++);
			for( l=0; l<tetStackLength     && tetStack[l]!=actualTet->neighborTetrahedra[j];     l++);
			// circumsphere contains new point ?
			if( k==usedTetStackLength && l==tetStackLength &&
					newVoronoiCell->getDistanceToCircumsphere( actualTet->neighborTetrahedra[j]) < 0./*getCircumsphereRadius( actualTet->neighborTetrahedra[j])*/){
				tetStack[tetStackLength] = actualTet->neighborTetrahedra[j];
				tetStackLength++;
				//fprintf( stderr, "INFO: Add neighbor tet %i (%i, %i, %i, %i) as candidate\n", actualTet->neighborTetrahedra[j]->index, actualTet->neighborTetrahedra[j]->vertices[0]->index, actualTet->neighborTetrahedra[j]->vertices[1]->index, actualTet->neighborTetrahedra[j]->vertices[2]->index, actualTet->neighborTetrahedra[j]->vertices[3]->index);

			}//else
				//fprintf( stderr, "INFO: Neighbor tet %i (%i, %i, %i, %i) doesn't contain new point\n", actualTet->neighborTetrahedra[j]->index, actualTet->neighborTetrahedra[j]->vertices[0]->index, actualTet->neighborTetrahedra[j]->vertices[1]->index, actualTet->neighborTetrahedra[j]->vertices[2]->index, actualTet->neighborTetrahedra[j]->vertices[3]->index);

		}

#else //TRACK_DELAUNAY_REGION_NEIGHBORHOOD

	for( j=0; j<ti; j++){

		Simplex<Dimensions> *actualTet = tetrahedra[j];
		if( getDistanceOfPointToCircumsphere( newVoronoiCell, actualTet) < 0. /*getCircumsphereRadius( actualTet)*/){
			//fprintf( stderr, "INFO: Delete tet %i (%i, %i, %i, %i)!\n", actualTet->index, actualTet->vertices[0]->index, actualTet->vertices[1]->index, actualTet->vertices[2]->index, actualTet->vertices[3]->index);

#endif //TRACK_DELAUNAY_REGION_NEIGHBORHOOD

		// ADD FACES

		// reallocation if nessessary
		if( fi+NR_Simplex_POINTS >= fimax){
			//fprintf( stderr, "INFO: realloc: %i -> %i!\n", fimax, fimax + FACES_ARRAY_EXTENSION_SIZE);
#if TRACK_DELAUNAY_REGION_NEIGHBORHOOD
			correspondingSimplex = (SimplexType**) realloc( correspondingSimplex, sizeof(SimplexType*) * (fimax + FACES_ARRAY_EXTENSION_SIZE));
#endif
			faces = (VertexType***) realloc( faces, sizeof(VertexType**) * (fimax + FACES_ARRAY_EXTENSION_SIZE));
			assert(faces);
			for( l=fimax; l<fimax + FACES_ARRAY_EXTENSION_SIZE; l++){
				faces[l] = (VertexType**) calloc( sizeof(VertexType*), NR_Face_POINTS);
				assert(faces[l]);
			}
			fimax += FACES_ARRAY_EXTENSION_SIZE;
		}

		for( k=0; k<NR_Simplex_POINTS; k++){ // for all faces of the Simplex
			// add face
			for( l=0; l<NR_Face_POINTS; l++){ // for all points of one face
				faces[fi][l]=actualTet->vertices[(l+k)%(NR_Simplex_POINTS)];
			}

#if TRACK_DELAUNAY_REGION_NEIGHBORHOOD

			//fprintf( stderr, "INFO: save neighbor tet!\n");
			for( l=0; l<actualTet->countNeighborTetrahedra &&
				!actualTet->neighborTetrahedra[l]->SimplexContainsFace( faces[fi]); l++);
			if( l!=actualTet->countNeighborTetrahedra){
				correspondingSimplex[fi] = actualTet->neighborTetrahedra[l];
				actualTet->neighborTetrahedra[l]->removeSimplexNeighbor( actualTet);
				actualTet->removeSimplexNeighbor( actualTet->neighborTetrahedra[l]);
			}else{
				correspondingSimplex[fi] = NULL;
			}

#endif

			// sort point labels of face
			VertexType * temp;
			int dim;
			for( dim=NR_Face_POINTS; dim>1; dim--){
				for( l=1; l<dim; l++){
					// (memory address)
					//if( faces[fi][l-1] > faces[fi][l])
					// (point index)
					if( faces[fi][l-1]->getIndex() > faces[fi][l]->getIndex())
					{
						temp = faces[fi][l-1];
						faces[fi][l-1] = faces[fi][l];
						faces[fi][l] = temp;
					}
				}
			}


			fi++;
			count_faces_added++;
		}


		// DELETE TETRAHEDRA FROM LIST

		ti--;
		int tempIndex = actualTet->index;
		SimplexType * temp = tetrahedra[tempIndex];

		tetrahedra[tempIndex] = tetrahedra[ti];
		tetrahedra[tempIndex]->index = tempIndex;

		tetrahedra[ti] = temp;
		tetrahedra[ti]->index = ti;

		count_tetrahedra_deleted++;

		/*for( k=0; k<ti; k++){
			if( k != tetrahedra[k]->index){
				fprintf( stderr, "ERROR: index of tet %i wrong: %i ! ti=%i\n", k, tetrahedra[k]->index, ti);
				exit( 0);
			}
		}*/

#if TRACK_DELAUNAY_REGION_NEIGHBORHOOD

	}while( tetStackLength != 0);

	free( tetStack);
	free( usedTetStack);

#else // TRACK_DELAUNAY_REGION_NEIGHBORHOOD

		j--;
	}}

#endif // TRACK_DELAUNAY_REGION_NEIGHBORHOOD



		// DELETE COMMON FACES FROM LIST

		//fprintf( stderr, "Number of Faces: fi=%d", fi);
		for( j=0; j<fi-1; j++){
			//fprintf( stderr, "Duplicate Check for Face %d\n", j);

			// check for copy of face j
			int duplicateFound = 0;
			for( l=j+1; l<fi && !duplicateFound; l++){ // for all following faces
				// check if faces j and l are equal
				for( k=0; k<Dimensions && faces[j][k] == faces[l][k]; k++); // for all points of the faces
				if( k==Dimensions && faces[j][k-1] == faces[l][k-1]){
					duplicateFound++;
					//fprintf( stderr, "%d. Face Duplicate found!!! Face %d (%d %d %d) and %d (%d %d %d) are equal. \n", duplicateFound, j, faces[j][0], faces[j][1], faces[j][2],  l, faces[l][0], faces[l][1], faces[l][2]);
					//exit(0);

					// remove neighborships
					for( k=0; k<NR_Face_POINTS; k++){
						faces[l][k]->removeNeighbor( faces[l][(k+1)%NR_Face_POINTS]);
						faces[l][(k+1)%NR_Face_POINTS]->removeNeighbor( faces[l][k]);
					}

					// remove face l
					//int faceContainsInsertedPoint = 0;
					for( k=0; k<NR_Face_POINTS; k++){
						faces[l][k]=faces[fi-1][k];
					}
#if TRACK_DELAUNAY_REGION_NEIGHBORHOOD
					correspondingSimplex[l]=correspondingSimplex[fi-1];
#endif
					//if
					//facesSimplexIndex[l]=facesSimplexIndex[fi-1];
					fi--;
					count_common_faces_deleted++;

					// remove face j
					for( k=0; k<NR_Face_POINTS; k++){
						faces[j][k]=faces[fi-1][k];
					}
#if TRACK_DELAUNAY_REGION_NEIGHBORHOOD
					correspondingSimplex[j]=correspondingSimplex[fi-1];
#endif
					fi--;
					count_common_faces_deleted++;
					j--;
				}
			}

		}
		//fprintf( stderr, "( after deleting common faces: fi=%d\n", fi);



		// ADD NEW TETRAHEDRA TO LIST

		for( j=0; j<fi; j++){ // for all faces
			if( ti == timax){
				//fprintf( stderr, "INFO: realloc: %i -> %i!\n", timax, timax + TETRAHEDRA_ARRAY_EXTENSION_SIZE);
				tetrahedra = (SimplexType **) realloc( tetrahedra, sizeof(SimplexType*)*(timax + TETRAHEDRA_ARRAY_EXTENSION_SIZE));
				assert(tetrahedra);
				for( k=timax; k<timax + TETRAHEDRA_ARRAY_EXTENSION_SIZE;k++)
					tetrahedra[k] = new SimplexType();
				timax += TETRAHEDRA_ARRAY_EXTENSION_SIZE;
			}

			//tetrahedra[ti] = new Simplex();
			//fprintf( stderr, "Add tetrahedra %d: (", ti);
			for( k=0; k<NR_Face_POINTS; k++){ // for all points of face
				tetrahedra[ti]->vertices[k] = faces[j][k]; // set point of face
				//fprintf( stderr, "%d ", tetrahedra[ti]->vertices[k]->index);
				newVoronoiCell->addNeighbor( faces[j][k]);
				faces[j][k]->addNeighbor( newVoronoiCell);
				faces[j][k]->addNeighbor( faces[j][(k+1)%NR_Face_POINTS]);
				faces[j][(k+1)%NR_Face_POINTS]->addNeighbor( faces[j][k]);
			}
			tetrahedra[ti]->vertices[k] = newVoronoiCell; // set added point
			//fprintf( stderr, "%d )\n", tetrahedra[ti]->vertices[k]->index);
			//tetrahedra[ti]->included=1;
			tetrahedra[ti]->countNeighborTetrahedra=0;
			tetrahedra[ti]->circumsphereInitialized=0;
			tetrahedra[ti]->index=ti;

			count_tetrahedra_added++;
			ti++;

		}
		//fprintf( stderr, ", add %d Tetrehedra\n", fi);
		//fprintf( stderr, "to %d already existing Tetrehedra\n", ti);


#if TRACK_DELAUNAY_REGION_NEIGHBORHOOD

		// DESTINE NEIGHBORSHIPS OF NEW TETRAHEDRA
		for( j=ti-fi/*-1*/; j<ti-1; j++){ // for each new tetrahedra j
		//for( j=ti-fi-1; j<ti-1; j++){ // for each new tetrahedra j


			for( k=j+1; k<ti; k++){ // for all other new tetrahedra k
				//fprintf( stderr, "Compare tetrahedra %i vs %i, max %i\n", j, k, ti);
				if( tetrahedra[j]->countCommonPointOfTetrahedra( tetrahedra[k]) == Dimensions){
					// j and k are neighbors
					tetrahedra[j]->addSimplexNeighbor( tetrahedra[k]);
					tetrahedra[k]->addSimplexNeighbor( tetrahedra[j]);
				}
			}
		}

		for( j=ti-fi/*-1*/; j<ti; j++){ // for each new tetrahedra j
			for( k=0; k<fi; k++){ // for all surrounding tetrahedra k
				//fprintf( stderr, "INFO: (%i)\n", k);
				if( correspondingSimplex[k]!=NULL && tetrahedra[j]->countCommonPointOfTetrahedra( correspondingSimplex[k]) == Dimensions){
					// j and k are neighbors
					tetrahedra[j]->addOrReplaceSimplexNeighbor( correspondingSimplex[k]);
					correspondingSimplex[k]->addOrReplaceSimplexNeighbor( tetrahedra[j]);
				}
			}
		}

#endif

	// ADD NEW POINT TO POINT ARRAY
	/*if( voronoiDiagram->countVertices==voronoiDiagram->maxVertices){
		voronoiDiagram->maxVertices += POINTS_ARRAY_EXTENSION_SIZE;
		voronoiDiagram->vertices = (Vertex<Dimensions> **) realloc( voronoiDiagram->vertices, sizeof(Vertex<Dimensions>*) * voronoiDiagram->maxVertices);
	}
	voronoiDiagram->vertices[voronoiDiagram->countVertices] = newVoronoiCell;
	voronoiDiagram->countVertices++;
	*/

	// FREE MEMORY

	//deleteMatrix( A, Dimensions+1);

#if TRACK_DELAUNAY_REGION_NEIGHBORHOOD
	free( correspondingSimplex);
#endif

	// RESET VALUES OF VORONOI DIAGRAM

	voronoiDiagram->faces = faces;
	voronoiDiagram->maxFaces = fimax;
	voronoiDiagram->countFaces = fi;

	voronoiDiagram->tetrahedra = tetrahedra;
	voronoiDiagram->maxTetrahedra = timax;
	voronoiDiagram->countTetrahedra = ti;
}


template <int Dimensions, class VertexType, class SimplexType>
void Triangulation<Dimensions,VertexType,SimplexType>::addAndTriangulate( VertexType* newVoronoiCell)
{
	Triangulation<Dimensions,VertexType,SimplexType>* voronoiDiagram = this;
	int j, k, l;
	int NR_Simplex_POINTS=Dimensions+1;
	int NR_Face_POINTS=Dimensions;


	//double **A = newMatrix( Dimensions + 1, Dimensions + 1);

	int count_faces_added = 0;
	int count_common_faces_deleted = 0;
	int count_tetrahedra_deleted = 0;
	int count_tetrahedra_added = 0;
	//int npr = 0;

	// faces
	int fi = 0;			// face counter
	int fimax = voronoiDiagram->maxFaces;			// maximal number of faces in list
	VertexType ***faces = voronoiDiagram->faces;		// list of faces


	// tetrahedra
	int ti = voronoiDiagram->countTetrahedra;
	int timax = voronoiDiagram->maxTetrahedra;
	SimplexType **tetrahedra = voronoiDiagram->tetrahedra;



	// DELAUNAY TRIANGULATION (WATSON)

#if TRACK_DELAUNAY_REGION_NEIGHBORHOOD

	SimplexType **correspondingSimplex = (SimplexType**) calloc( sizeof(SimplexType*), fimax);
	assert(correspondingSimplex);

	SimplexType *actualTet = voronoiDiagram->getSimplexContainingPointInCircumSphere( newVoronoiCell);

	// look for all tetrahedra deleted by the point insertion	
	int maxTetStackLength = TETRAHEDRA_ARRAY_EXTENSION_SIZE;		
	SimplexType **tetStack = (SimplexType**) calloc( sizeof(SimplexType*), TETRAHEDRA_ARRAY_EXTENSION_SIZE);
	assert(tetStack);
	int maxUsedTetStackLength = TETRAHEDRA_ARRAY_EXTENSION_SIZE;		
	SimplexType **usedTetStack = (SimplexType**) calloc( sizeof(SimplexType*), TETRAHEDRA_ARRAY_EXTENSION_SIZE);
	assert(usedTetStack);

	tetStack[0] = actualTet; //tetrahedra[0];
	int tetStackLength = 1;	
	int usedTetStackLength = 0;	


	// TEST:
	/*int test_count_tetrahedra_deleted = 0;
	for( j=0; j<ti; j++){
		if( getDistanceOfPointToCircumsphereCenter( newVoronoiCell, tetrahedra[j]) < getCircumsphereRadius( tetrahedra[j]))
			test_count_tetrahedra_deleted++;
	}*/
	// END TEST

	//fprintf( stderr, "look for all tetrahedra deleted by the point insertion\n");


	do{
		// TAKE NEXT Simplex FROM STACK
		tetStackLength--;
		actualTet = tetStack[tetStackLength];
		if( usedTetStackLength == maxUsedTetStackLength){
			maxUsedTetStackLength += TETRAHEDRA_ARRAY_EXTENSION_SIZE;
			usedTetStack = (SimplexType**) realloc( usedTetStack, sizeof(SimplexType*) * maxUsedTetStackLength);
			assert(usedTetStack);
		}
		usedTetStack[usedTetStackLength] = actualTet;
		usedTetStackLength++;

		//fprintf( stderr, "INFO: Delete tet %i (%i, %i, %i, %i)!\n", actualTet->index, actualTet->vertices[0]->index, actualTet->vertices[1]->index, actualTet->vertices[2]->index, actualTet->vertices[3]->index);

		// ADD NEIGHBOR TETRAHEDRA TO STACK
		//int temp_tetStackLength = tetStackLength;
		for( j=0; j<actualTet->countNeighborTetrahedra /*&& temp_tetStackLength == tetStackLength*/; j++){
			if( tetStackLength + actualTet->countNeighborTetrahedra < maxTetStackLength){
				maxTetStackLength += TETRAHEDRA_ARRAY_EXTENSION_SIZE;
				tetStack = (SimplexType**) realloc( tetStack, sizeof(SimplexType*) * maxTetStackLength);
				assert(tetStack);
			}
			//distanceToCircumsphereCenter = getDistanceOfPointToCircumsphereCenter( newVoronoiCell, actualTet->neighborTetrahedra[j]);
			for( k=0; k<usedTetStackLength && usedTetStack[k]!=actualTet->neighborTetrahedra[j]; k++);
			for( l=0; l<tetStackLength     && tetStack[l]!=actualTet->neighborTetrahedra[j];     l++);
			// circumsphere contains new point ?
			if( k==usedTetStackLength && l==tetStackLength &&
					newVoronoiCell->getDistanceOfPointToCircumsphere( actualTet->neighborTetrahedra[j]) < 0./*getCircumsphereRadius( actualTet->neighborTetrahedra[j])*/){
				tetStack[tetStackLength] = actualTet->neighborTetrahedra[j];
				tetStackLength++;
				//fprintf( stderr, "INFO: Add neighbor tet %i (%i, %i, %i, %i) as candidate\n", actualTet->neighborTetrahedra[j]->index, actualTet->neighborTetrahedra[j]->vertices[0]->index, actualTet->neighborTetrahedra[j]->vertices[1]->index, actualTet->neighborTetrahedra[j]->vertices[2]->index, actualTet->neighborTetrahedra[j]->vertices[3]->index);

			}//else
				//fprintf( stderr, "INFO: Neighbor tet %i (%i, %i, %i, %i) doesn't contain new point\n", actualTet->neighborTetrahedra[j]->index, actualTet->neighborTetrahedra[j]->vertices[0]->index, actualTet->neighborTetrahedra[j]->vertices[1]->index, actualTet->neighborTetrahedra[j]->vertices[2]->index, actualTet->neighborTetrahedra[j]->vertices[3]->index);

		}

#else //TRACK_DELAUNAY_REGION_NEIGHBORHOOD

	for( j=0; j<ti; j++){

		SimplexType *actualTet = tetrahedra[j];
		if( getDistanceOfPointToCircumsphere( newVoronoiCell, actualTet) < 0. /*getCircumsphereRadius( actualTet)*/){
			//fprintf( stderr, "INFO: Delete tet %i (%i, %i, %i, %i)!\n", actualTet->index, actualTet->vertices[0]->index, actualTet->vertices[1]->index, actualTet->vertices[2]->index, actualTet->vertices[3]->index);

#endif //TRACK_DELAUNAY_REGION_NEIGHBORHOOD

		// ADD FACES

		// reallocation if nessessary 
		if( fi+NR_Simplex_POINTS >= fimax){
			//fprintf( stderr, "INFO: realloc: %i -> %i!\n", fimax, fimax + FACES_ARRAY_EXTENSION_SIZE);
#if TRACK_DELAUNAY_REGION_NEIGHBORHOOD
			correspondingSimplex = (SimplexType**) realloc( correspondingSimplex, sizeof(SimplexType*) * (fimax + FACES_ARRAY_EXTENSION_SIZE));
			assert(correspondingSimplex);
#endif
			faces = (VertexType***) realloc( faces, sizeof(Vertex<Dimensions>**) * (fimax + FACES_ARRAY_EXTENSION_SIZE));
			assert(faces);
			for( l=fimax; l<fimax + FACES_ARRAY_EXTENSION_SIZE; l++){
				faces[l] = (VertexType**) calloc( sizeof(Vertex<Dimensions>*), NR_Face_POINTS);
				assert(faces[l]);
			}
			fimax += FACES_ARRAY_EXTENSION_SIZE;
		}

		for( k=0; k<NR_Simplex_POINTS; k++){ // for all faces of the Simplex
			// add face
			for( l=0; l<NR_Face_POINTS; l++){ // for all points of one face
				faces[fi][l]=actualTet->vertices[(l+k)%(NR_Simplex_POINTS)];
			}

#if TRACK_DELAUNAY_REGION_NEIGHBORHOOD

			//fprintf( stderr, "INFO: save neighbor tet!\n");
			for( l=0; l<actualTet->countNeighborTetrahedra && 
				!actualTet->neighborTetrahedra[l]->SimplexContainsFace( faces[fi]); l++);
			if( l!=actualTet->countNeighborTetrahedra){
				correspondingSimplex[fi] = actualTet->neighborTetrahedra[l];
				actualTet->neighborTetrahedra[l]->removeSimplexNeighbor( actualTet);
				actualTet->removeSimplexNeighbor( actualTet->neighborTetrahedra[l]);
			}else{
				correspondingSimplex[fi] = NULL;
			}

#endif	

			// sort point labels of face (memory address)
			VertexType * temp;
			int dim;
			for( dim=NR_Face_POINTS; dim>1; dim--){
				for( l=1; l<dim; l++){
					// (memory address)
					//if( faces[fi][l-1] > faces[fi][l])
					// (point index)
					if( faces[fi][l-1]->index > faces[fi][l]->index)
					{
						temp = faces[fi][l-1];
						faces[fi][l-1] = faces[fi][l];
						faces[fi][l] = temp;
					}
				}
			}

			fi++;			
			count_faces_added++;
		}


		// DELETE TETRAHEDRA FROM LIST

		ti--;
		int tempIndex = actualTet->index;
		SimplexType * temp = tetrahedra[tempIndex];

		tetrahedra[tempIndex] = tetrahedra[ti];
		tetrahedra[tempIndex]->index = tempIndex;

		tetrahedra[ti] = temp;
		tetrahedra[ti]->index = ti;
				
		count_tetrahedra_deleted++;

		/*for( k=0; k<ti; k++){
			if( k != tetrahedra[k]->index){
				fprintf( stderr, "ERROR: index of tet %i wrong: %i ! ti=%i\n", k, tetrahedra[k]->index, ti);
				exit( 0);
			}
		}*/

#if TRACK_DELAUNAY_REGION_NEIGHBORHOOD

	}while( tetStackLength != 0);

	free( tetStack);
	free( usedTetStack);

#else // TRACK_DELAUNAY_REGION_NEIGHBORHOOD

		j--;
	}}

#endif // TRACK_DELAUNAY_REGION_NEIGHBORHOOD



		// DELETE COMMON FACES FROM LIST

		//fprintf( stderr, "Number of Faces: fi=%d", fi);
		for( j=0; j<fi-1; j++){
			//fprintf( stderr, "Duplicate Check for Face %d\n", j);

			// check for copy of face j
			int duplicateFound = 0;
			for( l=j+1; l<fi && !duplicateFound; l++){ // for all following faces
				// check if faces j and l are equal
				for( k=0; k<Dimensions && faces[j][k] == faces[l][k]; k++); // for all points of the faces
				if( k==Dimensions && faces[j][k-1] == faces[l][k-1]){
					duplicateFound++;
					//fprintf( stderr, "%d. Face Duplicate found!!! Face %d (%d %d %d) and %d (%d %d %d) are equal. \n", duplicateFound, j, faces[j][0], faces[j][1], faces[j][2],  l, faces[l][0], faces[l][1], faces[l][2]);
					//exit(0);

					// remove neighborships
					for( k=0; k<NR_Face_POINTS; k++){
						faces[l][k]->removeNeighbor( faces[l][(k+1)%NR_Face_POINTS]);
						faces[l][(k+1)%NR_Face_POINTS]->removeNeighbor( faces[l][k]);
					}
					
					// remove face l
					//int faceContainsInsertedPoint = 0;
					for( k=0; k<NR_Face_POINTS; k++){
						faces[l][k]=faces[fi-1][k];
					}
#if TRACK_DELAUNAY_REGION_NEIGHBORHOOD
					correspondingSimplex[l]=correspondingSimplex[fi-1];
#endif
					//if
					//facesSimplexIndex[l]=facesSimplexIndex[fi-1];
					fi--;
					count_common_faces_deleted++;
					
					// remove face j
					for( k=0; k<NR_Face_POINTS; k++){
						faces[j][k]=faces[fi-1][k];
					}
#if TRACK_DELAUNAY_REGION_NEIGHBORHOOD
					correspondingSimplex[j]=correspondingSimplex[fi-1];
#endif
					fi--;
					count_common_faces_deleted++;
					j--;
				}
			}		
			
		}
		//fprintf( stderr, "( after deleting common faces: fi=%d\n", fi);
			


		// ADD NEW TETRAHEDRA TO LIST

		for( j=0; j<fi; j++){ // for all faces
			if( ti == timax){
				//fprintf( stderr, "INFO: realloc: %i -> %i!\n", timax, timax + TETRAHEDRA_ARRAY_EXTENSION_SIZE);
				tetrahedra = (SimplexType **) realloc( tetrahedra, sizeof(SimplexType*)*(timax + TETRAHEDRA_ARRAY_EXTENSION_SIZE));
				if( tetrahedra==NULL){
					fprintf( stderr, "ERROR reallocating memory for Simplex array!\n");
					exit(0);
				}
				for( k=timax; k<timax + TETRAHEDRA_ARRAY_EXTENSION_SIZE;k++)
					tetrahedra[k] = new SimplexType();
				timax += TETRAHEDRA_ARRAY_EXTENSION_SIZE;
			}

			//tetrahedra[ti] = new Simplex();
			//fprintf( stderr, "Add tetrahedra %d: (", ti);
			for( k=0; k<NR_Face_POINTS; k++){ // for all points of face
				tetrahedra[ti]->vertices[k] = faces[j][k]; // set point of face
				//fprintf( stderr, "%d ", tetrahedra[ti]->vertices[k]->index);
				newVoronoiCell->addNeighbor( faces[j][k]);
				faces[j][k]->addNeighbor( newVoronoiCell);
				faces[j][k]->addNeighbor( faces[j][(k+1)%NR_Face_POINTS]);
				faces[j][(k+1)%NR_Face_POINTS]->addNeighbor( faces[j][k]);
			}
			tetrahedra[ti]->vertices[k] = newVoronoiCell; // set added point
			//fprintf( stderr, "%d )\n", tetrahedra[ti]->vertices[k]->index);
			//tetrahedra[ti]->included=1;
			tetrahedra[ti]->countNeighborTetrahedra=0;
			tetrahedra[ti]->circumsphereInitialized=0;
			tetrahedra[ti]->index=ti;
			
			count_tetrahedra_added++;
			ti++;

		}
		//fprintf( stderr, ", add %d Tetrehedra\n", fi);


#if TRACK_DELAUNAY_REGION_NEIGHBORHOOD
           		
		// DESTINE NEIGHBORSHIPS OF NEW TETRAHEDRA
		for( j=ti-fi; j<ti; j++){ // for each new tetrahedra j
		//for( j=ti-fi-1; j<ti-1; j++){ // for each new tetrahedra j

           		
			for( k=j+1; k<ti; k++){ // for all other new tetrahedra k
				if( tetrahedra[j]->countCommonPointOfTetrahedra( tetrahedra[k]) == Dimensions){
					// j and k are neighbors
					tetrahedra[j]->addSimplexNeighbor( tetrahedra[k]);
					tetrahedra[k]->addSimplexNeighbor( tetrahedra[j]);
				}
			}
		}
           		
		for( j=ti-fi-1; j<ti; j++){ // for each new tetrahedra j
			for( k=0; k<fi; k++){ // for all surrounding tetrahedra k
				//fprintf( stderr, "INFO: (%i)\n", k);
				if( correspondingSimplex[k]!=NULL && tetrahedra[j]->countCommonPointOfTetrahedra( correspondingSimplex[k]) == Dimensions){
					// j and k are neighbors
					tetrahedra[j]->addOrReplaceSimplexNeighbor( correspondingSimplex[k]);
					correspondingSimplex[k]->addOrReplaceSimplexNeighbor( tetrahedra[j]);
				}
			}
		}

#endif

	// ADD NEW POINT TO POINT ARRAY
	if( voronoiDiagram->countVertices==voronoiDiagram->maxVertices){
		voronoiDiagram->maxVertices += POINTS_ARRAY_EXTENSION_SIZE;
		voronoiDiagram->vertices = (Vertex<Dimensions> **) realloc( voronoiDiagram->vertices, sizeof(Vertex<Dimensions>*) * voronoiDiagram->maxVertices);
		assert(voronoiDiagram->vertices);
	}
	voronoiDiagram->vertices[voronoiDiagram->countVertices] = newVoronoiCell;
	voronoiDiagram->countVertices++;

	// FREE MEMORY

	//deleteMatrix( A, Dimensions+1);

#if TRACK_DELAUNAY_REGION_NEIGHBORHOOD
	free( correspondingSimplex);
#endif

	// RESET VALUES OF VORONOI DIAGRAM

	voronoiDiagram->faces = faces;
	voronoiDiagram->maxFaces = fimax;
	voronoiDiagram->countFaces = fi; 

	voronoiDiagram->tetrahedra = tetrahedra;
	voronoiDiagram->maxTetrahedra = timax;
	voronoiDiagram->countTetrahedra = ti; 
}
/*****************************************************************************/


template <int Dimensions, class VertexType>
Simplex<Dimensions,VertexType>::Simplex()
{
	// initialize values
	this->circumsphereInitialized = 0;
	
	this->countNeighborTetrahedra = 0;
}
/*****************************************************************************/


template <int Dimensions>
Vertex<Dimensions>::~Vertex()
{
	free(this->neighbors);
}
/*****************************************************************************/


template <int Dimensions>
void Vertex<Dimensions>::Initialize()
{
	// neighborhood
	this->index = 0;
	this->neighbors = 0;
	this->countNeighbors = 0;
	this->neighborsInitialized = FALSE;
}
/*****************************************************************************/


template <>
inline Vertex<2>::Vertex(double x, double y)
{
	// initialize values
	this->position[0] = x;
	this->position[1] = y;

	Initialize();
}
/*****************************************************************************/


template <>
inline Vertex<3>::Vertex(double x, double y, double z)
{
	// initialize values
	this->position[0] = x;
	this->position[1] = y;
	this->position[2] = z;

	Initialize();
}
/*****************************************************************************/


template <>
inline Vertex<1>::Vertex(double x)
{
	// initialize values
	this->position[0] = x;

	Initialize();
}
/*****************************************************************************/


template <int Dimensions>
Vertex<Dimensions>::Vertex()
{
	// initialize values
	for( int d=0; d<Dimensions; d++)
		this->position[d] = 0;

	Initialize();
}
/*****************************************************************************/


template <int Dimensions>
Vertex<Dimensions>::Vertex( double x, double y, double z)
{
	// initialize values
	this->position[0] = x;
	this->position[1] = y;
	this->position[2] = z;

	Initialize();
}
template <int Dimensions>
int Vertex<Dimensions>::numberOfNeighbors()
{	return countNeighbors;}

template <int Dimensions>
Vertex<Dimensions>* Vertex<Dimensions>::getNeighbor( int i)
{	return neighbors[i];}

/*****************************************************************************/
template <int Dimensions, class VertexType, class SimplexType>
inline int Triangulation<Dimensions,VertexType,SimplexType>::numberOfVertices()
{
	return this->countVertices;
}

/*****************************************************************************/

template <int Dimensions, class VertexType, class SimplexType>
Triangulation<Dimensions,VertexType,SimplexType>::~Triangulation()
{
	int NR_Simplex_POINTS=Dimensions+1;

	//fprintf( stderr, "Triangulation<Dimensions>::~Triangulation()\n");
	int i;
	
	// free tetrahedra
	for( i=0; i<this->maxTetrahedra; i++)
		delete( (this->tetrahedra[i]));
	free( this->tetrahedra);

	// free faces
	for( i=0; i<this->maxFaces; i++)
		free( (this->faces[i]));
	free( this->faces);

	// free points
	for( i=0; i<this->countVertices; i++){
		delete this->vertices[i];
		//free( this->vertices[i]->neighbors);
		//free( this->vertices[i]);
	}
	free( this->vertices);
	
	// free frame points
	for( i=0; i<this->countFramePoints; i++){
		delete this->framePoints[i];
		//free( this->framePoints[i]->neighbors);
		//free( this->framePoints[i]);
	}
	free( this->framePoints);

	// free voronoi grid
	/*if( voronoiGridSet == TRUE){
		free( this->voronoiGridPoints);
		free( this->voronoiCellCountVoronoiGridPoints);
		for( int i=0; i<this->countVertices; i++)
			free( this->voronoiCellToVoronoiGridPoints[i]);
		free( this->voronoiCellToVoronoiGridPoints);
	}*/
	
	SimplexType::countUsers--;
	if( SimplexType::countUsers==0){
		Matrix<double>::deleteMatrix(SimplexType::_A, NR_Simplex_POINTS);
		Matrix<double>::deleteMatrix(SimplexType::_B, NR_Simplex_POINTS);
		free(SimplexType::_x);
	}
}
/*****************************************************************************/


template <int Dimensions, class VertexType, class SimplexType>
Triangulation<Dimensions,VertexType,SimplexType>::Triangulation()
{
	Initialize();
}
/*****************************************************************************/


template <int Dimensions,class VertexType,class SimplexType>
Triangulation<Dimensions,VertexType,SimplexType>::Triangulation( const char* filename){

	Initialize();

	/****************************************************************************
	 * Definieren der Variablen                                                 *
	 ****************************************************************************/

		//int Dimensions = 3;

		//Triangulation<Dimensions>	*newVoronoiDiagram = new Triangulation<Dimensions>();
	Triangulation<Dimensions,VertexType,SimplexType>	*newVoronoiDiagram = this;

	    char fin[ FILENAMESIZE ];
	    char f_node[ FILENAMESIZE ];
		sprintf( fin,	"%s.qele", filename );
		sprintf( f_node,"%s.qnode", filename );
	    FILE*	pin;
	    FILE*	fp_node;

	    int		i, j, l, new_point,
				t, ta, tx,
				tetrahedra	= 0,
				last		= 0,
				lastpoint	= 0;
	    long	triangle[ Dimensions+2 ];

	    char	buffer[ READBUFFERSIZE ],
	            cTri[ READBUFFERSIZE ],
			    buffer_node[ READBUFFERSIZE ];
	    char*	ptr = buffer;

	/****************************************************************************
	 * Einlesen der Dreiecke			                                        *
	 ****************************************************************************/



		/* open input file (*.qele)*/
	    pin = fopen( fin, "r" );
	    if( pin == NULL ){
	      fprintf( stderr, "Error opening file %s\n", fin);
	      exit( 1);
	    }
	    else fprintf( stderr, ">%s eingelesen\n", fin);

		/* open input file (*.qnode)*/
	    fp_node = fopen( f_node, "r" );
	    if( fp_node == NULL ){
	      fprintf( stderr, "Error opening file %s\n", f_node);
	      exit( 1 );
	    }
	    else fprintf( stderr, ">%s eingelesen\n", f_node);

		/* erste Zeile einlesen und Anzahl der Dreiecke auslesen */
	    if( fgets( cTri, READBUFFERSIZE, pin ) == NULL )
	      perror( "Error reading first line" );

	    /* maximale Punkte = 3*Dreiecke;
	     * jeder Punkt kann maximal doppelt so viele Nachbarn haben,
	     * wie es Dreiecke gibt */
	    tetrahedra = atoi( cTri );
	    newVoronoiDiagram->countTetrahedra = tetrahedra;
	    if (tetrahedra > MAX_TRIANGLES) {
	       printf("\n Too many triangle (%d)", tetrahedra);
	       exit(0);
	       }
	    //static int point[MAX_TRIANGLES][MAX_NNS]; /* * 2; * 3 */
	    	fgets( buffer_node, READBUFFERSIZE, fp_node );

		//printf("Allocate Memory!\n");
		lastpoint = atoi( buffer_node);
		//printf("%i Points to initialize!\n", lastpoint);
		int **point;
		int *pointAllocSize;
		point = (int**) calloc( sizeof(int*), lastpoint);
		assert(point);
		pointAllocSize = (int*) calloc( sizeof(int), lastpoint);
		assert(pointAllocSize);
		for( i=0; i<lastpoint; i++){
			point[i] = (int*) calloc( sizeof(int), MAX_NNS);
			assert(point[i]);
			pointAllocSize[i] = MAX_NNS;
		}
		//printf("Finished: Allocate Memory!\n");

	/****************************************************************************
	 * Dreieck- in Punktbeziehung wandeln                                       *
	 ****************************************************************************/
		newVoronoiDiagram->tetrahedra = ( SimplexType **) calloc ( tetrahedra, sizeof(SimplexType *));
		assert(newVoronoiDiagram->tetrahedra);
		for( t = 0; t < tetrahedra; t++ ){
			newVoronoiDiagram->tetrahedra[t] = new SimplexType(); //( Simplex<Dimensions> *) calloc ( Dimensions+1, sizeof( Simplex));
			newVoronoiDiagram->tetrahedra[t]->setIndex(t);
		}
		newVoronoiDiagram->maxTetrahedra =      tetrahedra;
		newVoronoiDiagram->countTetrahedra =    tetrahedra;

		/* Anzahl der Punkte speichern */
		newVoronoiDiagram->countVertices = lastpoint;            // Anzahl der Punkte speichern
		newVoronoiDiagram->maxVertices = lastpoint;            // Anzahl der Punkte speichern
		newVoronoiDiagram->vertices = ( VertexType**) malloc ( newVoronoiDiagram->countVertices * sizeof( VertexType*));
		assert(newVoronoiDiagram->vertices);
		for( i=0; i<newVoronoiDiagram->countVertices; i++){
			newVoronoiDiagram->vertices[i] = new VertexType();
			assert(newVoronoiDiagram->vertices[i]);
		}

		/* alle Dreiecke verarbeiten */
		for( t = 0; t < tetrahedra; t++ ){
	#if _COMMENTS_ > 2
			fprintf( stderr, "\rRead and Process Tetrahedra Data %.3lf%%                         \b", 100.*(t+1.)/tetrahedra);
	#endif
			if( fgets( buffer, READBUFFERSIZE, pin ) == NULL )
				perror( "Error reading line1" );

			ptr = buffer;
			triangle[ 0 ] = t;

			/* alle Punkte eines Dreiecks einlesen */
			for( tx = 1; tx < Dimensions+2; tx++){
				triangle[ tx ]  = strtol( ptr, &ptr, 0 );
				//newVoronoiDiagram->tetrahedra[t]->vertices[tx-1] = newVoronoiDiagram->vertices[triangle[ tx ]];
				newVoronoiDiagram->tetrahedra[t]->setVertex( newVoronoiDiagram->vertices[triangle[ tx ]], tx-1);

			}

			/* zu jedem Punkt die Nachbarn dieses Dreiecks sichern */
			for( ta = 1; ta < Dimensions+2; ta++ ){
				if( point[ triangle[ ta ] ][ 0 ] == 0 ){	// Punkt neu?
					point[ triangle[ ta ] ][ 0 ] = 1;
					last = 1;
				if( lastpoint < triangle[ ta ] )
					lastpoint = triangle[ ta ];
				}else
					last = point[ triangle[ ta ] ][ 0 ];	// Speicherposition
			for( tx = 1; tx < Dimensions+2; tx++ ){
				new_point = 1;
				if( triangle[ tx ] == triangle[ ta ] )	// identischer Punkt
					new_point = 0;
				else
					for( l = 1; l < last; l++ )		// alle Nachbarn pruefen
						if( point[ triangle[ ta ] ][ l ] == triangle[ tx ] )
							new_point = 0;	// Punkt ist schon Nachbar
					if( new_point != 0 ){
						if( pointAllocSize[triangle[ ta ]] == last){
							pointAllocSize[triangle[ ta ]] += MAX_NNS;
							point[triangle[ ta ]] = (int*) realloc( point[triangle[ ta ]], sizeof(int) * pointAllocSize[triangle[ ta ]]);
							assert(point[triangle[ ta ]]);
						}
						point[ triangle[ ta ] ][ last++ ] = triangle[ tx ];
					}
				}
				point[ triangle[ ta ] ][ 0 ] = last;
			}
		}

	/****************************************************************************
	 * Dreieckbeziehung einlesen                                                *
	 ****************************************************************************/
		if( fgets( buffer, READBUFFERSIZE, pin ) != NULL )
		for( t = 0; t < tetrahedra; t++ ){
	#if _COMMENTS_ > 2
			fprintf( stderr, "\rDetermine Neighborhood Relations between vertices %.3lf%%                    \b", 100.*(t+1.)/tetrahedra);
	#endif
			if( fgets( buffer, READBUFFERSIZE, pin ) == NULL )
				perror( "Error reading line1" );

			ptr = buffer;

			/* alle Punkte eines Dreiecks einlesen */
			strtol( ptr, &ptr, 0 );
			for( tx = 0; tx < Dimensions+1; tx++){
				int index  = strtol( ptr, &ptr, 0 );
				if( index >= 0)
					newVoronoiDiagram->tetrahedra[t]->addSimplexNeighbor( newVoronoiDiagram->tetrahedra[index]);
				//newVoronoiDiagram->tetrahedra[t]->vertices[tx-1] = &newVoronoiDiagram->vertices[triangle[ tx ]];
			}
		}else{
			perror( "Error reading line1" );
		}

	/*	for( i=0; i<newVoronoiDiagram->countTetrahedra-1; i++){
			// find neighbors of Simplex i

			for( j=i+1; j<newVoronoiDiagram->countTetrahedra; j++){
				// check if Simplex j is a neighbor of Simplex i

				if( countCommonPointOfTetrahedra( newVoronoiDiagram->tetrahedra[i], newVoronoiDiagram->tetrahedra[j]) == Dimensions){
					addSimplexNeighbor( newVoronoiDiagram->tetrahedra[i], newVoronoiDiagram->tetrahedra[j]);
					addSimplexNeighbor( newVoronoiDiagram->tetrahedra[j], newVoronoiDiagram->tetrahedra[i]);
				}
			}
		}
		int k;
		for( i=0; i<newVoronoiDiagram->countTetrahedra; i++){
			for( j=0; j<newVoronoiDiagram->tetrahedra[i]->countNeighborTetrahedra; j++){
				for( k=0; k<newVoronoiDiagram->tetrahedra[i]->neighborTetrahedra[j]->countNeighborTetrahedra &&
					newVoronoiDiagram->tetrahedra[i] != newVoronoiDiagram->tetrahedra[i]->neighborTetrahedra[j]->neighborTetrahedra[k]; k++);

				if( k==newVoronoiDiagram->tetrahedra[i]->neighborTetrahedra[j]->countNeighborTetrahedra){
					fprintf( stderr, "ERROR: neighborship relation between tets %p -> %p!!\n", newVoronoiDiagram->tetrahedra[i], newVoronoiDiagram->tetrahedra[i]->neighborTetrahedra[j]->neighborTetrahedra[k]);
					exit( 0);
				}
			}
		}*/


	/****************************************************************************
	 *  alle Punkte in Datei schreiben                                          *
	 ****************************************************************************/

		/* Gitter initialisieren*/


		/* Punkte mit Nachbarn in Datei sichern */
		for( i = 0; i<lastpoint; i++ ){
	#if _COMMENTS_ > 2
			fprintf( stderr, "\rInitialize Triangulation %.3lf%%                              \b", 100.*(i+1.)/newVoronoiDiagram->countVertices);
	#endif
			if( fgets( buffer_node, READBUFFERSIZE, fp_node ) == NULL )
				perror( "Error reading line4" );
			ptr = buffer_node;

		  //printf("initialize point %i\n", i);
			/* set number of voronoi cell */
			newVoronoiDiagram->vertices[i]->setIndex(i);

			//newVoronoiDiagram->vertices[i]->agent = NULL;

			//newVoronoiDiagram->vertices[i]->refined = false;

			/* set coordinates */
			int dim;
			for( dim=0; dim<Dimensions; dim++)
				newVoronoiDiagram->vertices[i]->setx( strtod( ptr, &ptr ), dim);

			/* number of neighbours */
			//newVoronoiDiagram->vertices[i]->countNeighbors = point[ i ][ 0 ] - 1;
			//newVoronoiDiagram->vertices[i]->countFreeNeighborCells = point[ i ][ 0 ] - 1;

			/* set neighbors of voronoi cell */
			//newVoronoiDiagram->vertices[i]->neighbors = ( Vertex<Dimensions>**) calloc ( point[ i ][ 0 ] - 1, sizeof( Vertex<Dimensions>*)); // allocate memory for list of neighbors
			//assert(newVoronoiDiagram->vertices[i]->neighbors);
			for( j = 0; j < point[i][0] - 1; j++) {
				newVoronoiDiagram->vertices[i]->addNeighbor( newVoronoiDiagram->vertices[point[ i ][ j + 1]]);
				//newVoronoiDiagram->vertices[i]->neighbors[j] = newVoronoiDiagram->vertices[point[ i ][ j + 1]]; // Nachbarn in Liste eifimaxuegen
			}
			//newVoronoiDiagram->vertices[i]->neighborsInitialized = TRUE;

			// extended neighbors
			//newVoronoiDiagram->vertices[i]->countExtendedNeighborCells = 0;
			//newVoronoiDiagram->vertices[i]->countFreeExtendedNeighborCells = 0;
			//newVoronoiDiagram->vertices[i]->extendedNeighborhood = NULL;
			//newVoronoiDiagram->vertices[i]->extendedNeighborCellsInitialized = FALSE;

		}
		fprintf( stderr, "\n");

	#ifdef LIMITED_NEIGHBOR_DISTANCE
		/*for( int i=0; i<newVoronoiDiagram->countVertices; i++){
			int countReachable = 0;
			for( int ii=0; ii<newVoronoiDiagram->vertices[i]->countNeighbors; ii++){
				if( abs( (int)(newVoronoiDiagram->vertices[i]->position[0]) - (int)(newVoronoiDiagram->vertices[i]->neighbors[ii]->position[0])) < LIMITED_NEIGHBOR_DISTANCE
				    &&(Dimensions < 2 || abs( (int)(newVoronoiDiagram->vertices[i]->position[1]) - (int)(newVoronoiDiagram->vertices[i]->neighbors[ii]->position[1])) < LIMITED_NEIGHBOR_DISTANCE)
				    &&(Dimensions < 3 || abs( (int)(newVoronoiDiagram->vertices[i]->position[2]) - (int)(newVoronoiDiagram->vertices[i]->neighbors[ii]->position[2])) < LIMITED_NEIGHBOR_DISTANCE)
				){

					//countReachable++;
					newVoronoiDiagram->vertices[i]->neighbors[countReachable++] = newVoronoiDiagram->vertices[i]->neighbors[ii];
				}
			}

			newVoronoiDiagram->vertices[i]->countNeighbors     = countReachable;
			//newVoronoiDiagram->vertices[i]->countFreeNeighborCells = countReachable;
		}*/
	#endif

		/*for( i = 0; i<lastpoint; i++ ){
			printf( "point %i: ", newVoronoiDiagram->vertices[i].index);
			for( j = 0; j < newVoronoiDiagram->vertices[i].countNeighbors; j++)
				printf( " %i", newVoronoiDiagram->vertices[i].neighbors[j]->index);
			printf( "\n");
		}*/

		/* Ein- und Ausgabedatei schließen */
	    if( fclose( pin ) )
	      printf( "can't close %s\n", fin );
	    if( fclose( fp_node ) )
	      printf( "can't close %s\n", f_node );

		for( i=0; i<newVoronoiDiagram->countVertices; i++)
			free( point[i]);
		free( point);
		free( pointAllocSize);

		this->setDomain();
}
/*****************************************************************************/

template <int Dimensions,class VertexType, class SimplexType>
void Triangulation<Dimensions, VertexType, SimplexType>::Initialize()
{
#ifndef NDEBUG
	fprintf(stderr, "Triangulation used in Debug Mode\nIf you want to speed up recompile with \"-D NDEBUG\"\n");
#endif

	int NR_Simplex_POINTS=Dimensions+1;

	// init vertices
	vertices = 0;
	this->countVertices = 0;
	this->maxVertices = 0;

	// init tets
	this->tetrahedra = 0;
	this->countTetrahedra = 0;
	this->maxTetrahedra = 0;

	// voronoiDiagram->faces;
	this->faces = 0;
	this->maxFaces = 0;
	this->countFaces = 0;

	// init Voronoi diagram
	//voronoiGridSet = false;

	// init domain information
	domainSet = false;

	SimplexType::countUsers++;
	if( SimplexType::countUsers==1){
		SimplexType::_A = newDoubleMatrix( NR_Simplex_POINTS, NR_Simplex_POINTS);
		SimplexType::_B = newDoubleMatrix( NR_Simplex_POINTS, NR_Simplex_POINTS);
		SimplexType::_x = (double*) malloc( sizeof(double) * NR_Simplex_POINTS);
		assert( SimplexType::_x);
		for( int i=0; i<NR_Simplex_POINTS; i++){
			for( int ii=0; ii<NR_Simplex_POINTS; ii++){
				SimplexType::_A[i][ii]=0;
				SimplexType::_B[i][ii]=0;
			}
			SimplexType::_x[i]=0;
		}
	}
}

template <int Dimensions,class VertexType, class SimplexType>
Triangulation<Dimensions, VertexType, SimplexType>::Triangulation( int x, int y, int z){

	Initialize();

	Triangulation<Dimensions, VertexType, SimplexType>	*newVoronoiDiagram = this;


	newVoronoiDiagram->countVertices = x*y*z;            // Anzahl der Punkte speichern
	newVoronoiDiagram->maxVertices = x*y*z;            // Anzahl der Punkte speichern
	newVoronoiDiagram->vertices = ( VertexType**) malloc ( newVoronoiDiagram->countVertices * sizeof( VertexType*));
	assert(	newVoronoiDiagram->vertices);
	//for( i=0; i<newVoronoiDiagram->countVertices; i++)
	//	newVoronoiDiagram->vertices[i] = ( Vertex<Dimensions>*) malloc ( sizeof( Vertex));



	/* Punkte mit Nachbarn in Datei sichern */
	int i=0;
	// Regular Mesh
	//double alpha = 0.0001;
	// Irregular Mesh
	double alpha = 1.;

	int ix[3] = {0,0,0};
	for( ix[0] = 0; ix[0]<x; ix[0]++ )
	for( ix[1] = 0; ix[1]<y; ix[1]++ )
	for( ix[2] = 0; ix[2]<z; ix[2]++ )
		{
#if _COMMENTS_ > 2
		fprintf( stderr, "\rInitialize Triangulation %.3lf%%                              \b", 100.*(i+1.)/newVoronoiDiagram->countVertices);
#endif

		newVoronoiDiagram->vertices[i] = new VertexType();
		for( int d=0; d<Dimensions; d++)
			newVoronoiDiagram->vertices[i]->setx( (float)(ix[d]+alpha*myRandE(1.)), d);

		/* set number of voronoi cell */
		newVoronoiDiagram->vertices[i]->setIndex(i);
		i++;
	}
	fprintf( stderr, "\n");

	newVoronoiDiagram->setFramePoints();
	newVoronoiDiagram->triangulate();

	// Post-Triangulation Initialization
	//for( int i=0; i<newVoronoiDiagram->countVertices; i++ )
	//	newVoronoiDiagram->vertices[i]->countFreeNeighborCells = newVoronoiDiagram->vertices[i]->countNeighbors;

}
/*****************************************************************************/

template <int Dimensions, class VertexType, class SimplexType>
void Triangulation<Dimensions, VertexType, SimplexType>::printPoints( const char * title, VertexType ** allPoints, int api )
{
	int j;
	int negativeValues = 0;
		// print points
		fprintf( stderr, "%s ", title);	
		for( j=0; j<api; j++){
			fprintf( stderr, "%i ", allPoints[j]->index);
			if( allPoints[j]->index<0)
				negativeValues++;	
		}
		fprintf( stderr, "\n");
		if( negativeValues)
			exit( 0);	
}
/*****************************************************************************/


template <int Dimensions, class VertexType, class SimplexType>
void Triangulation<Dimensions,VertexType,SimplexType>::checkDelaunayCondition()
{
	Triangulation<Dimensions,VertexType,SimplexType> *voronoiDiagram = this;
	int ii, l;
	int NR_Simplex_POINTS=Dimensions+1;

	fprintf( stderr, "Check Delaunay Condition of Voronoi Diagram...\n");

	for( int i=0; i<voronoiDiagram->countTetrahedra; i++)
		if( voronoiDiagram->tetrahedra[i]->circumsphereInitialized == 0)
			 voronoiDiagram->tetrahedra[i]->initCircumsphereRadius();

#pragma omp parallel for private(l,ii)
	//for( int i=voronoiDiagram->countTetrahedra-1; i>=0; i--){
	for( int i=0; i<voronoiDiagram->countTetrahedra; i++){

		for( ii=0; ii<voronoiDiagram->countVertices; ii++){
			if( voronoiDiagram->vertices[ii]->countNeighbors > 0){

				for( l=0; l<NR_Simplex_POINTS && voronoiDiagram->vertices[ii]!=voronoiDiagram->tetrahedra[i]->vertices[l]; l++);
				if( 	l==NR_Simplex_POINTS &&
						voronoiDiagram->vertices[ii]->getDistanceOfPointToCircumsphere( voronoiDiagram->tetrahedra[i]) < 0.){
					fprintf( stderr, "... Delaunay Condition is not valid!!! :-(\n");
					/*fprintf( stderr, "%i\n", voronoiDiagram->tetrahedra[i]->index);
					fprintf( stderr, "%i\n", voronoiDiagram->tetrahedra[i]->vertices[0]->index);
					fprintf( stderr, "%i\n", voronoiDiagram->tetrahedra[i]->vertices[1]->index);
					fprintf( stderr, "%i\n", voronoiDiagram->tetrahedra[i]->vertices[2]->index);
					//fprintf( stderr, "%i\n", voronoiDiagram->tetrahedra[i]->vertices[3]->index);
					fprintf( stderr, "%i\n", voronoiDiagram->vertices[ii]->index);
					fprintf( stderr, "%f\n", voronoiDiagram->vertices[ii]->position[0]);
					fprintf( stderr, "%f\n", voronoiDiagram->vertices[ii]->position[1]);*/
					//fprintf( stderr, "%f\n", voronoiDiagram->vertices[ii]->position[2]);

					/*fprintf( stderr, "tet %i (%i [%p] %i [%p] %i [%p]) should not contain point %i [%p] (%lf %lf)\n",
										voronoiDiagram->tetrahedra[i]->index,
										voronoiDiagram->tetrahedra[i]->vertices[0]->index,
										voronoiDiagram->tetrahedra[i]->vertices[0],
										voronoiDiagram->tetrahedra[i]->vertices[1]->index,
										voronoiDiagram->tetrahedra[i]->vertices[1],
										voronoiDiagram->tetrahedra[i]->vertices[2]->index,
										voronoiDiagram->tetrahedra[i]->vertices[2],
										voronoiDiagram->vertices[ii]->index,
										voronoiDiagram->vertices[ii],
										voronoiDiagram->vertices[ii]->position[0],
										voronoiDiagram->vertices[ii]->position[1]);
*/
					fprintf( stderr, "tet %i (%i %i %i %i) should not contain point %i (%lf %lf %lf)\n",
						voronoiDiagram->tetrahedra[i]->index,
						voronoiDiagram->tetrahedra[i]->vertices[0]->index,
						voronoiDiagram->tetrahedra[i]->vertices[1]->index,
						voronoiDiagram->tetrahedra[i]->vertices[2]->index,
						voronoiDiagram->tetrahedra[i]->vertices[3]->index,
						voronoiDiagram->vertices[ii]->index,
						voronoiDiagram->vertices[ii]->position[0],
						voronoiDiagram->vertices[ii]->position[1],
						voronoiDiagram->vertices[ii]->position[2]);
					fprintf( stderr, "vertex %i (%lf %lf %lf)\n",
							voronoiDiagram->tetrahedra[i]->vertices[0]->index,
							voronoiDiagram->tetrahedra[i]->vertices[0]->position[0],
							voronoiDiagram->tetrahedra[i]->vertices[0]->position[1],
							voronoiDiagram->tetrahedra[i]->vertices[0]->position[2]);
					fprintf( stderr, "vertex %i (%lf %lf %lf)\n",
							voronoiDiagram->tetrahedra[i]->vertices[1]->index,
							voronoiDiagram->tetrahedra[i]->vertices[1]->position[0],
							voronoiDiagram->tetrahedra[i]->vertices[1]->position[1],
							voronoiDiagram->tetrahedra[i]->vertices[1]->position[2]);
					fprintf( stderr, "vertex %i (%lf %lf %lf)\n",
							voronoiDiagram->tetrahedra[i]->vertices[2]->index,
							voronoiDiagram->tetrahedra[i]->vertices[2]->position[0],
							voronoiDiagram->tetrahedra[i]->vertices[2]->position[1],
							voronoiDiagram->tetrahedra[i]->vertices[2]->position[2]);
					fprintf( stderr, "vertex %i (%lf %lf %lf)\n",
							voronoiDiagram->tetrahedra[i]->vertices[3]->index,
							voronoiDiagram->tetrahedra[i]->vertices[3]->position[0],
							voronoiDiagram->tetrahedra[i]->vertices[3]->position[1],
							voronoiDiagram->tetrahedra[i]->vertices[3]->position[2]);
					/*fprintf( stderr, "distance of point to circumsphere of tet = %e\n", getDistanceOfPointToCircumsphere( voronoiDiagram->vertices[ii], voronoiDiagram->tetrahedra[i]));
					fprintf( stderr, "distance of point to circumsphere of tet = %e\n", getDistanceOfPointToCircumsphere( voronoiDiagram->tetrahedra[i]->vertices[0], voronoiDiagram->tetrahedra[i]));
					fprintf( stderr, "distance of point to circumsphere of tet = %e\n", getDistanceOfPointToCircumsphere( voronoiDiagram->tetrahedra[i]->vertices[1], voronoiDiagram->tetrahedra[i]));
					fprintf( stderr, "distance of point to circumsphere of tet = %e\n", getDistanceOfPointToCircumsphere( voronoiDiagram->tetrahedra[i]->vertices[2], voronoiDiagram->tetrahedra[i]));
					fprintf( stderr, "distance of point to circumsphere of tet = %e\n", getDistanceOfPointToCircumsphere( voronoiDiagram->tetrahedra[i]->vertices[3], voronoiDiagram->tetrahedra[i]));*/

					exit( 0);
				}
			}else{
				//fprintf( stderr, "Point %i is not part of the triangulation and thus will be ignored!\n", voronoiDiagram->vertices[ii]->index);
			}
			/*for( ii=0; ii<countNewCells; ii++){
				for( l=0; l<NR_Simplex_POINTS && newCell[ii]!=voronoiDiagram->tetrahedra[i]->vertices[l]; l++);
				if( l==NR_Simplex_POINTS && getDistanceOfPointToCircumsphere( newCell[ii], voronoiDiagram->tetrahedra[i]) < 0){
					fprintf( stderr, "... Delaunay Condition is not valid!!! :-(\n");
					exit( 0);
				}
			}*/
		}
	}

	fprintf( stderr, "... Delaunay Condition is valid :-)\n");
}

/*
Hey salut Luna!
tu as encore reussi hier de prendre le dernier metro? Pour nous c'etait un peu chaud...en plus avec des trajets prolonges de la forme sinusoidaire a cause la biere ;) Je t'ecrit juste pour te proposer a boire quelque chose ce soir/weekend dans le cas que il te tombe le plafond sur la tete dans ta maison model ikea :) Il aura meme Chiara a Paris ce weekend je croit.

 */

/*****************************************************************************/


template <int Dimensions, class VertexType, class SimplexType>
void Triangulation<Dimensions,VertexType,SimplexType>::getConvexHull( double thresholdDistance)
{
	int NR_Simplex_POINTS=Dimensions+1;
	int NR_Face_POINTS=Dimensions;

	this->convexFaces = (VertexType ***)malloc( sizeof(VertexType **) * this->countVertices);
	this->countConvexFaces = 0;
	this->maxConvexFaces = this->countVertices;
	//int countConvexFaces;
	//int maxConvexFaces;
	//Vertex<Dimensions> ***convexFaces;

	SimplexType ** innerTets = (SimplexType **)malloc( sizeof(SimplexType *) * this->countVertices);
	int countInnerTets = 0;
	int maxInnerTets = this->countVertices;

	SimplexType ** usedTets = (SimplexType **)malloc( sizeof(SimplexType *) * this->countVertices);
	int countUsedTets = 0;
	int maxUsedTets = this->countVertices;

	int i, t, f;

	// ADD ALL FACES CONNECTED TO FRAME

	for( i=0; i<this->countTetrahedra; i++){ // for each existing tet
		// check if tet contains one of the frame points
		f=this->countFramePoints;
		int countFramePointsInTet = 0;
		int tt = 0, ff = 0;
		for( t=0; t<NR_Simplex_POINTS; t++){
			//fprintf( stderr, "%i ", this->tetrahedra[i]->vertices[t]->index);
			for( f=0; f<this->countFramePoints && this->framePoints[f]!=this->tetrahedra[i]->vertices[t]; f++);
			if( f<this->countFramePoints){
				countFramePointsInTet++;
				tt = t;
				ff = f;
				//fprintf( stderr, "(f) ");
			}
		}
		//fprintf( stderr, "\n");
		//fprintf( stderr, "INFO: tet %i contains %i frame points\n", i, countFramePoints);


		int j, jj;
		double dist = 0.;
		double maxDist = 0.;
		if( countFramePointsInTet==1){
//		if( countFramePointsInTet>=1){
			/*fprintf( stderr, "INFO: vertice %i of tet %i (%i %i %i %i) contains frame point %i\n",
			         this->tetrahedra[i]->vertices[tt]->index, i,
			         this->tetrahedra[i]->vertices[0]->index,
			         this->tetrahedra[i]->vertices[1]->index,
			         this->tetrahedra[i]->vertices[2]->index,
			         this->tetrahedra[i]->vertices[3]->index,
			         this->framePoints[ff]->index);
			*/
			for( j=0; j<NR_Simplex_POINTS-1; j++){
				for( jj=1; jj<NR_Simplex_POINTS; jj++){
					if( j!=tt && jj!=tt){
						dist = pow( this->tetrahedra[i]->vertices[j]->position[0] - this->tetrahedra[i]->vertices[jj]->position[0], 2)
						     + pow( this->tetrahedra[i]->vertices[j]->position[1] - this->tetrahedra[i]->vertices[jj]->position[1], 2)
						     + pow( this->tetrahedra[i]->vertices[j]->position[2] - this->tetrahedra[i]->vertices[jj]->position[2], 2);
						dist = sqrt( dist);
						//fprintf( stderr, "INFO: distance of vertice %i and %i is %lf\n", this->tetrahedra[i]->vertices[j]->index, this->tetrahedra[i]->vertices[jj]->index, dist);
						if( maxDist<dist)
							maxDist = dist;
					}
				}
			}
		//}

		//if( countFramePoints==1 && maxDist<thresholdDistance){
			if( maxDist<thresholdDistance){
				if( this->countConvexFaces == this->maxConvexFaces){
					this->maxConvexFaces += FACES_ARRAY_EXTENSION_SIZE;
					this->convexFaces = (VertexType ***)realloc( this->convexFaces, sizeof(VertexType **) * this->maxConvexFaces);
				}

				t = tt+1;
				f = ff;

				//tet contains one of the frame points
				//fprintf( stderr, "INFO: vertice %i of tet %i contains frame point %i\n", this->tetrahedra[i]->vertices[t-1]->index, i, this->framePoints[f]->index);

				// add face
				this->convexFaces[this->countConvexFaces] = (VertexType **)malloc( sizeof(VertexType *) * NR_Face_POINTS);
				//fprintf( stderr, "Add %ith Convex Face: ", this->countConvexFaces);
				for( f=0; f<NR_Face_POINTS; f++){
					this->convexFaces[this->countConvexFaces][f] = this->tetrahedra[i]->vertices[t];
					//fprintf( stderr, " %i", this->convexFaces[this->countConvexFaces][f]->index);

					t++;
					t = t%NR_Simplex_POINTS;
					//t = (++t)%NR_Simplex_POINTS;
				}
				//fprintf( stderr, "\n");
				this->countConvexFaces++;

			}else{
				// find inner neighbors of actual tet
				//fprintf( stderr, "Add inner tets: ");
				for( j=0; j<this->tetrahedra[i]->countNeighborTetrahedra; j++){
					for( jj=0; jj<NR_Simplex_POINTS && this->tetrahedra[i]->neighborTetrahedra[j]->vertices[jj]!=this->framePoints[ff]; jj++);
					if(jj==NR_Simplex_POINTS){
						// inner tet
						for( jj=0; jj<countInnerTets && innerTets[jj]!=this->tetrahedra[i]->neighborTetrahedra[j]; jj++);
						if(jj==countInnerTets){
							if( countInnerTets == maxInnerTets){
								maxInnerTets += FACES_ARRAY_EXTENSION_SIZE;
								innerTets = (SimplexType **)realloc( innerTets, sizeof(SimplexType *) * maxInnerTets);
							}
							if( countUsedTets == maxUsedTets){
								maxUsedTets += FACES_ARRAY_EXTENSION_SIZE;
								usedTets = (SimplexType **)realloc( usedTets, sizeof(SimplexType *) * maxUsedTets);
							}
							innerTets[countInnerTets++] = this->tetrahedra[i]->neighborTetrahedra[j];
							usedTets[countUsedTets++] = this->tetrahedra[i];

							//fprintf( stderr, "%i (used:%i)", this->tetrahedra[i]->neighborTetrahedra[j]->index, this->tetrahedra[i]->index);
						}//else
					}
				}
				//fprintf( stderr, "\n");


			}
		}else{
			//fprintf( stderr, "INFO: tet %i doesn't contain any frame point\n", i);
		}

	}

	// ADD INNER FACES

	for( i=0; i<countInnerTets; i++){ // for each inner tet
		int j, jj, ii;
		VertexType* tempFace[Dimensions];

		for( ii=0; ii<NR_Simplex_POINTS; ii++){ // for each face

			// actual face
			for( j=0; j<NR_Face_POINTS; j++)
				tempFace[j] = innerTets[i]->vertices[(ii+j)%NR_Simplex_POINTS];

			// distance of face points
			double dist = 0.;
			double maxDist = 0.;

			for( j=0; j<NR_Face_POINTS-1; j++){
				for( jj=j+1; jj<NR_Face_POINTS; jj++){
					dist = pow( tempFace[j]->position[0] - tempFace[jj]->position[0], 2)
					     + pow( tempFace[j]->position[1] - tempFace[jj]->position[1], 2)
					     + pow( tempFace[j]->position[2] - tempFace[jj]->position[2], 2);
					dist = sqrt( dist);
					//fprintf( stderr, "INFO: distance of face points %i and %i is %lf (%i-%i)\n", tempFace[j]->index, tempFace[jj]->index, dist, j, jj);
					if( maxDist<dist)
						maxDist = dist;
				}
			}
			if( maxDist<thresholdDistance){
				if( this->countConvexFaces == this->maxConvexFaces){
					this->maxConvexFaces += FACES_ARRAY_EXTENSION_SIZE;
					this->convexFaces = (VertexType ***)realloc( this->convexFaces, sizeof(VertexType **) * this->maxConvexFaces);
				}

				// add face
				this->convexFaces[this->countConvexFaces] = (VertexType **)malloc( sizeof(VertexType *) * NR_Face_POINTS);
				//fprintf( stderr, "Add %ith Convex Face: ", this->countConvexFaces);
				for( j=0; j<NR_Face_POINTS; j++){
					this->convexFaces[this->countConvexFaces][j] = tempFace[j];
					//fprintf( stderr, " %i", this->convexFaces[this->countConvexFaces][j]->index);
				}
				//fprintf( stderr, "\n");
				this->countConvexFaces++;
				//fprintf( stderr, "INFO: distance of vertice %i and %i is %lf\n", this->tetrahedra[i]->vertices[j]->index, this->tetrahedra[i]->vertices[jj]->index, dist);
			}else{
				// add corresponding neighbor tet to inner tet list
				for( j=0; j<innerTets[i]->countNeighborTetrahedra && !innerTets[i]->neighborTetrahedra[j]->SimplexContainsFace( tempFace); j++);
				if( j<innerTets[i]->countNeighborTetrahedra){
					//fprintf( stderr, "Add %ith Neighbor tet: %i\n", j, innerTets[i]->neighborTetrahedra[j]->index);
					// is already used?
					for( jj=0; jj<countUsedTets && usedTets[jj]!=innerTets[i]->neighborTetrahedra[j]; jj++);
					if( jj<countUsedTets){
						//fprintf( stderr, "Neighbor tet %i is already used!!!\n", innerTets[i]->neighborTetrahedra[j]->index);
						//exit( 0);
					}else{
						// insert?
						//fprintf( stderr, "Add inner tets: ");
						for( jj=0; jj<countInnerTets && innerTets[jj]!=innerTets[i]->neighborTetrahedra[j]; jj++);
						if(jj==countInnerTets){
							if( countInnerTets == maxInnerTets){
								maxInnerTets += FACES_ARRAY_EXTENSION_SIZE;
								innerTets = (SimplexType **)realloc( innerTets, sizeof(SimplexType *) * maxInnerTets);
							}
							if( countUsedTets == maxUsedTets){
								maxUsedTets += FACES_ARRAY_EXTENSION_SIZE;
								usedTets = (SimplexType **)realloc( usedTets, sizeof(SimplexType *) * maxUsedTets);
							}
							innerTets[countInnerTets++] = innerTets[i]->neighborTetrahedra[j];
							usedTets[countUsedTets++] = innerTets[i];

							//fprintf( stderr, "%i (used:%i)", innerTets[i]->neighborTetrahedra[j]->index, innerTets[i]->index);
						}
						//fprintf( stderr, "\n");
					}
				}
			}
		}
	}
}
/*****************************************************************************/

template <int Dimensions, class VertexType, class SimplexType>
void Triangulation<Dimensions,VertexType,SimplexType>::setConvexHullNeighborhood()
{
	int NR_Face_POINTS=Dimensions;
	int i, ii;
	//fprintf( stderr, "countVertices=%i\n", countVertices);
	for( i=0; i<this->countVertices; i++){
		this->vertices[i]->countNeighbors = 0;
	}

	for( i=0; i<this->countConvexFaces; i++){
		for( ii=0; ii<NR_Face_POINTS; ii++){
			this->convexFaces[i][ii]->addNeighbor( this->convexFaces[i][(ii+1)%NR_Face_POINTS]);
			this->convexFaces[i][(ii+1)%NR_Face_POINTS]->addNeighbor( this->convexFaces[i][ii]);
		}
	}


}
/*****************************************************************************/

template <int Dimensions, class VertexType,class SimplexType>
void Triangulation<Dimensions, VertexType, SimplexType>::printToVTK( const char *filename)
{
	FILE *fp = fopen( filename, "w+");

	// HEADER
	fprintf( fp, "# vtk DataFile Version 2.0\n");
	fprintf( fp, "Really cool data\n");
	fprintf( fp, "ASCII\n");
	fprintf( fp, "DATASET UNSTRUCTURED_GRID\n");

	// POINTS
	fprintf( fp, "POINTS %i float\n", this->countVertices+this->countFramePoints);
	for( int i=0; i<countVertices; i++){
		for( int d=0; d<Dimensions; d++){
			fprintf( fp, "%f ", this->vertices[i]->pos()[d]);
		}
		fprintf( fp, "\n");
	}
	for( int i=0; i<countFramePoints; i++){
		for( int d=0; d<Dimensions; d++){
			fprintf( fp, "%f ", this->framePoints[i]->pos()[d]);
		}
		fprintf( fp, "\n");
	}


	// TETRAHEDRA
	fprintf( fp, "CELLS %i %i\n", this->countTetrahedra, countTetrahedra*(1+1+Dimensions));
	for( int t=0; t<countTetrahedra; t++){
		fprintf( fp, "4 ");
		for( int v=0; v<Dimensions+1; v++){
			if( tetrahedra[t]->vertices[v]->getIndex()<0)
				fprintf( fp, "%i ", countVertices - 1 - tetrahedra[t]->vertices[v]->getIndex());
			else
				fprintf( fp, "%i ", this->tetrahedra[t]->vertices[v]->getIndex());
		}
		fprintf( fp, "\n");
	}

	fprintf( fp, "CELL_TYPES %i\n", countTetrahedra);
	for( int t=0; t<countTetrahedra; t++)
		fprintf( fp, "10\n");

	fclose(fp);
}

/*****************************************************************************/

template <int Dimensions, class VertexType, class SimplexType>
void Triangulation<Dimensions, VertexType, SimplexType>::printToPovray( const char *filename,
		bool printPoints, bool printFramePoints, bool printNeighborship, bool printTetrahedra, bool printConvexHull)
{
	int NR_Simplex_POINTS=Dimensions+1;
	int NR_Face_POINTS=Dimensions;

	Triangulation<Dimensions,VertexType,SimplexType> *voronoiDiagram=this;
	VertexType **newCell=0;
	int countNewCells=0;

	//int NR_Simplex_POINTS = Dimensions+1;
	int i, j, k, l;

	FILE* fp;

	fp = fopen( filename, "w+" );

	// PRINT GRID TO POV-FILE
	fprintf( fp, "#include \"finish.inc\"\n"\
				"#include \"colors.inc\"\n"\
				"#include \"textures.inc\"\n"\
				"background { color White }\n"\
				"camera {location <%lf, %lf, %lf> look_at  <%lf, %lf, 0>}\n"\
				"light_source {<5, 5, 12> color rgb <1, 1, 1>}\n",
				0.5*(domain->min(0)+domain->max(0)), 0.5*(domain->min(1)+domain->max(1)), MAX(domain->max(0)-domain->min(0), domain->max(1)-domain->min(1)),
				0.5*(domain->min(0)+domain->max(0)), 0.5*(domain->min(1)+domain->max(1)));

	// all points (GRAY)
	if( printPoints)
	for( i=0; i<voronoiDiagram->countVertices; i++){
		{
		fprintf( fp, "sphere {<");
		for( j=0; j<Dimensions; j++){
			fprintf( fp, " %.3lf", voronoiDiagram->vertices[i]->pos()[j]);
			if( j+1!=Dimensions)
				fprintf( fp, ",");
		}
		fprintf( fp, ">,0.1 texture {pigment { color rgbf <0.95, 0.95, 0.95, 0.2> }}}\n");
		}

		// all neighborships (YELLOW)
		if( printNeighborship){
		for( j=0; j<voronoiDiagram->vertices[i]->numberOfNeighbors(); j++)
			if(printFramePoints || voronoiDiagram->vertices[i]->getNeighbor(j)->getIndex() >= 0){
			fprintf( fp, "  cylinder {<%lf, %lf, %lf>, <%lf, %lf, %lf>, %lf pigment { color %s }}\n",
	    		voronoiDiagram->vertices[i]->pos()[0],
				voronoiDiagram->vertices[i]->pos()[1],
				voronoiDiagram->vertices[i]->pos()[2],
	       		voronoiDiagram->vertices[i]->pos()[0]
				+ 0.4 * (voronoiDiagram->vertices[i]->getNeighbor(j)->pos()[0] - voronoiDiagram->vertices[i]->pos()[0]),
	       		voronoiDiagram->vertices[i]->pos()[1]
				+ 0.4 * (voronoiDiagram->vertices[i]->getNeighbor(j)->pos()[1] - voronoiDiagram->vertices[i]->pos()[1]),
				voronoiDiagram->vertices[i]->pos()[2]
				+ 0.4 * (voronoiDiagram->vertices[i]->getNeighbor(j)->pos()[2] - voronoiDiagram->vertices[i]->pos()[2]),
				0.05,
	           	"Yellow");
			fprintf( fp, "\n");
		}
		}
	}

	// all frame points (GRAY)
	if( printFramePoints)
	for( i=0; i<voronoiDiagram->countFramePoints; i++){
		{
		fprintf( fp, "sphere {<");
		for( j=0; j<Dimensions; j++){
			fprintf( fp, " %.3lf", voronoiDiagram->framePoints[i]->pos()[j]);
			if( j+1!=Dimensions)
				fprintf( fp, ",");
		}
		fprintf( fp, ">,0.1 texture {pigment { color rgbf <0.95, 0.95, 0.95, 0.2> }}}\n");
		}

		// all neighborships (YELLOW)
		if( printNeighborship){
		for( j=0; j<voronoiDiagram->framePoints[i]->numberOfNeighbors(); j++)
			if(printPoints || voronoiDiagram->framePoints[i]->getNeighbor(j)->getIndex() < 0){
			fprintf( fp, "  cylinder {<%lf, %lf, %lf>, <%lf, %lf, %lf>, %lf pigment { color %s }}\n",
	    		voronoiDiagram->framePoints[i]->pos()[0],
				voronoiDiagram->framePoints[i]->pos()[1],
				voronoiDiagram->framePoints[i]->pos()[2],
	       		voronoiDiagram->framePoints[i]->pos()[0]
				+ 0.4 * (voronoiDiagram->framePoints[i]->getNeighbor(j)->pos()[0] - voronoiDiagram->framePoints[i]->pos()[0]),
	       		voronoiDiagram->framePoints[i]->pos()[1]
				+ 0.4 * (voronoiDiagram->framePoints[i]->getNeighbor(j)->pos()[1] - voronoiDiagram->framePoints[i]->pos()[1]),
				voronoiDiagram->framePoints[i]->pos()[2]
				+ 0.4 * (voronoiDiagram->framePoints[i]->getNeighbor(j)->pos()[2] - voronoiDiagram->framePoints[i]->pos()[2]),
				0.05,
	           	"Yellow");
			fprintf( fp, "\n");
		}
		}
	}

	// convex hull (BLUE)
	if( printConvexHull)
	for( j=0; j<voronoiDiagram->countConvexFaces; j++){
		for( k=0; k<NR_Face_POINTS; k++)
		for( l=k+1; l<NR_Face_POINTS; l++){ // TEST
			fprintf( fp, "  cylinder {<%lf, %lf, %lf>, <%lf, %lf, %lf>, %lf pigment { color %s }}\n",
	    			voronoiDiagram->convexFaces[j][k]->pos()[0],
				voronoiDiagram->convexFaces[j][k]->pos()[1],
				voronoiDiagram->convexFaces[j][k]->pos()[2],
				//voronoiDiagram->tetrahedra[j]->vertices[(k+1)%(Dimensions+1)]->position[0],
				//voronoiDiagram->tetrahedra[j]->vertices[(k+1)%(Dimensions+1)]->position[1],
				//voronoiDiagram->tetrahedra[j]->vertices[(k+1)%(Dimensions+1)]->position[2],
				voronoiDiagram->convexFaces[j][l]->pos()[0],
				voronoiDiagram->convexFaces[j][l]->pos()[1],
				voronoiDiagram->convexFaces[j][l]->pos()[2],
				0.04,
           		"rgbf <0, 0, 1, 0.2>");//"Blue");
			// zugehoerigkeitsmarkierung
			/*double centerX, centerY, centerZ;
			centerX = 0.5 * (voronoiDiagram->tetrahedra[j]->vertices[k]->position[0] + voronoiDiagram->tetrahedra[j]->vertices[(k+1)%(Dimensions+1)]->position[0]);
			centerY = 0.5 * (voronoiDiagram->tetrahedra[j]->vertices[k]->position[1] + voronoiDiagram->tetrahedra[j]->vertices[(k+1)%(Dimensions+1)]->position[1]);
			centerZ = 0.5 * (voronoiDiagram->tetrahedra[j]->vertices[k]->position[2] + voronoiDiagram->tetrahedra[j]->vertices[(k+1)%(Dimensions+1)]->position[2]);
			fprintf( fp, "  cylinder {<%lf, %lf, %lf>, <%lf, %lf, %lf>, %lf pigment { color %s }}\n",
	    		centerX,
				centerY,
				centerZ,
				0.8 * centerX + 0.2 * voronoiDiagram->tetrahedra[j]->vertices[(k+2)%(Dimensions+1)]->position[0],
	       		0.8 * centerY + 0.2 * voronoiDiagram->tetrahedra[j]->vertices[(k+2)%(Dimensions+1)]->position[1],
				0.8 * centerZ + 0.2 * voronoiDiagram->tetrahedra[j]->vertices[(k+2)%(Dimensions+1)]->position[2],
				0.02,
           		"rgbf <0, 0, 1, 0.2>");//"Blue");*/
		}
		fprintf( fp, "\n");
       		fprintf( fp, "triangle {<%lf, %lf, %lf>, <%lf, %lf, %lf>, <%lf, %lf, %lf> pigment { color %s }}\n",
       		        voronoiDiagram->convexFaces[j][0]->pos()[0],
			voronoiDiagram->convexFaces[j][0]->pos()[1],
			voronoiDiagram->convexFaces[j][0]->pos()[2],
			voronoiDiagram->convexFaces[j][1]->pos()[0],
			voronoiDiagram->convexFaces[j][1]->pos()[1],
			voronoiDiagram->convexFaces[j][1]->pos()[2],
			voronoiDiagram->convexFaces[j][2]->pos()[0],
			voronoiDiagram->convexFaces[j][2]->pos()[1],
			voronoiDiagram->convexFaces[j][2]->pos()[2],
			//"rgbf <0, 0, 1, 0.2>"
			"rgbf <0.95, 0.95, 0.95, 0.2> "
			);
	}


	// all tetrahedra (BLUE)
	if( printTetrahedra)
	for( j=0; j<voronoiDiagram->countTetrahedra; j++){
		for( k=0; k<NR_Simplex_POINTS; k++)
		for( l=k+1; l<NR_Simplex_POINTS; l++){ // TEST
			fprintf( fp, "  cylinder {<%lf, %lf, %lf>, <%lf, %lf, %lf>, %lf pigment { color %s }}\n",
	    			voronoiDiagram->tetrahedra[j]->vertices[k]->pos()[0],
				voronoiDiagram->tetrahedra[j]->vertices[k]->pos()[1],
				voronoiDiagram->tetrahedra[j]->vertices[k]->pos()[2],
				//voronoiDiagram->tetrahedra[j]->vertices[(k+1)%(Dimensions+1)]->position[0],
				//voronoiDiagram->tetrahedra[j]->vertices[(k+1)%(Dimensions+1)]->position[1],
				//voronoiDiagram->tetrahedra[j]->vertices[(k+1)%(Dimensions+1)]->position[2],
				voronoiDiagram->tetrahedra[j]->vertices[l]->pos()[0],
				voronoiDiagram->tetrahedra[j]->vertices[l]->pos()[1],
				voronoiDiagram->tetrahedra[j]->vertices[l]->pos()[2],
				0.04,
           		"rgbf <0, 0, 1, 0.2>");//"Blue");
			// zugehoerigkeitsmarkierung
			/*double centerX, centerY, centerZ;
			centerX = 0.5 * (voronoiDiagram->tetrahedra[j]->vertices[k]->position[0] + voronoiDiagram->tetrahedra[j]->vertices[(k+1)%(Dimensions+1)]->position[0]);
			centerY = 0.5 * (voronoiDiagram->tetrahedra[j]->vertices[k]->position[1] + voronoiDiagram->tetrahedra[j]->vertices[(k+1)%(Dimensions+1)]->position[1]);
			centerZ = 0.5 * (voronoiDiagram->tetrahedra[j]->vertices[k]->position[2] + voronoiDiagram->tetrahedra[j]->vertices[(k+1)%(Dimensions+1)]->position[2]);
			fprintf( fp, "  cylinder {<%lf, %lf, %lf>, <%lf, %lf, %lf>, %lf pigment { color %s }}\n",
	    		centerX,
				centerY,
				centerZ,
				0.8 * centerX + 0.2 * voronoiDiagram->tetrahedra[j]->vertices[(k+2)%(Dimensions+1)]->position[0],
	       		0.8 * centerY + 0.2 * voronoiDiagram->tetrahedra[j]->vertices[(k+2)%(Dimensions+1)]->position[1],
				0.8 * centerZ + 0.2 * voronoiDiagram->tetrahedra[j]->vertices[(k+2)%(Dimensions+1)]->position[2],
				0.02,
           		"rgbf <0, 0, 1, 0.2>");//"Blue");*/
		}
		fprintf( fp, "\n");
	}

	for( i=0; i<countNewCells; i++){
	// new point
	fprintf( fp, "sphere {<");
	for( j=0; j<Dimensions; j++){
		fprintf( fp, " %.3lf", newCell[i]->pos()[j]);
		if( j+1!=Dimensions)
			fprintf( fp, ",");
	}
	fprintf( fp, ">,0.11 texture {pigment { color rgbf <0.95, 0., 0.95, 0.2> }}}\n");

	// all neighborships
	for( j=0; j<newCell[i]->numberOfNeighbors(); j++){
		//for( k=0; k<NR_Simplex_POINTS; k++){ // TEST
			fprintf( fp, "  cylinder {<%lf, %lf, %lf>, <%lf, %lf, %lf>, %lf pigment { color %s }}\n",
	    		newCell[i]->pos()[0],
				newCell[i]->pos()[1],
				newCell[i]->pos()[2],
				newCell[i]->pos()[0]
				+ 0.4 * (newCell[i]->getNeighbor(j)->pos()[0] - newCell[i]->pos()[0]),
	       		newCell[i]->pos()[1]
				+ 0.4 * (newCell[i]->getNeighbor(j)->pos()[1] - newCell[i]->pos()[1]),
				newCell[i]->pos()[2]
				+ 0.4 * (newCell[i]->getNeighbor(j)->pos()[2] - newCell[i]->pos()[2]),
				0.051,
           		"rgbf <0.95, 0., 0.95, 0.2>");
		//}
		fprintf( fp, "\n");
	}
	}


	fclose( fp);
}
/*****************************************************************************/


template <int Dimensions, class VertexType, class SimplexType>
void Triangulation<Dimensions,VertexType,SimplexType>::setPointsOnSphericalSurface( double radius, int N, double variance, double center[Dimensions])
{
	int k;
	double height, theta, phi, alpha, var_angle;
	double x, y, z, x_moved, y_moved, z_moved;

	// extend Voronoi diagram
	this->vertices = (VertexType**) realloc( this->vertices, sizeof(VertexType*)*(this->countVertices+N));
	assert(this->vertices);

	// polar coordinates of point k of N points
	//alpha = 0.5*2*acos(1-1/(2*PI*radius*radius*density)); // wenn Flächeninhalt der Kaloppe gleich 1/densitiy
	alpha = acos(1.-2./(double)N)*variance;

	/*for( countNewCells=1; countNewCells<=maxNewCells; countNewCells++){
		srand ( countNewCells+randomGeneratorSeed);
		//sprintf( filename, "test%i.pov", countNewCells);
		//testCell = newVoronoiCell( myRand()*8.+1., myRand()*8.+1., myRand()*8.+1.);
		//testCell = newVoronoiCell( myRand()*60., myRand()*60., myRand()*60.);
		testCell = newVoronoiCell( myRand()*2., myRand()*2., myRand()*2.);
		//printf("Insert %i. Voronoi Cell ( %lf, %lf, %lf)\n", countNewCells, testCell->position[0], testCell->position[1], testCell->position[2]);
		testCell->index = testGrid->countVertices;
		insertVoronoiCell( testGrid, testCell);
		//newCells[countNewCells-1] = testCell;
		//printToPovray( filename, testGrid, newCells, countNewCells);

	}*/
	srand(0);

	for(k=1; k<=N; k++){

		// spherical coordinates
		if( k==1){
			height = -1.;
			theta = acos(height);
			phi = 0.;
		}
		if( 1<k && k<N){
			height = -1. + 2.*((double)k-1.)/((double)N-1.);
			theta = acos(height);
			phi = (phi + 3.6/(sqrt((double)N) * sqrt(1.-height*height)));
		}
		if( k==N ){
			height = 1.;
			theta = acos(height);
			phi = 0.;
		}

		/*x = 1.;
		y = 0.;
		z = 0.;

		//rotateXYZ( x, y, z, x_moved, y_moved, z_moved, phi, 0, theta);
		rotateZ( x, y, z, x, y, z, theta);
		rotateX( x, y, z, x_moved, y_moved, z_moved, phi);*/


		// cartesian coordinates
		x = radius * sin(theta) * cos(phi);
		y = radius * sin(theta) * sin(phi);
		z = radius * cos(theta);

		// random angle of variance
		var_angle = 0.;//asin( sin(alpha) * myRand());

		x_moved = radius * sin(theta - var_angle) * cos(phi);
		y_moved = radius * sin(theta - var_angle) * sin(phi);
		z_moved = radius * cos(theta - var_angle);

		rotate( &x_moved, &y_moved, &z_moved, x, y, z, 2*PI*myRand()-PI);



		// set new voronoi cell
		//this->vertices[this->countVertices] = newVoronoiCell( x_moved, y_moved, z_moved);
		//this->countVertices++;
		VertexType* testCell = new VertexType( center[0]+x_moved+myRand()*POINT_POSITION_PRECISION, center[1]+y_moved+myRand()*POINT_POSITION_PRECISION, center[2]+z_moved+myRand()*POINT_POSITION_PRECISION);
		testCell->index = this->countVertices;
		this->add( testCell);

		// set spherical coordinates
		//testCell->sphericalCoordinates.phi = phi;
		//testCell->sphericalCoordinates.theta = theta;

		//printf("%lf %lf %lf\n", x_moved, y_moved, z_moved);
		//printf("%lf %lf %lf\n", x, y, z);
	}
}
/*****************************************************************************/


template <int Dimensions, class VertexType, class SimplexType>
void Triangulation<Dimensions,VertexType,SimplexType>::triangulate()
{
	int i;

	for( i=0; i<this->countVertices; i++){
		//printf("\rAdd Cell %i \b", i);
		//fprintf( stderr, "Add Cell %i \n", i);
		fprintf( stderr, "\r%.3lf%%\t\t\b", 100.*(i+1.)/this->countVertices);
		triangulate( this->vertices[i]);
		//checkDelaunayCondition( this, 0, 0);
//#pragma omp parallel for
		/*for( int t=0; t<this->countTetrahedra; t++){
			if(this->tetrahedra[t]->countNeighborTetrahedra != 4 && !this->SimplexContainsFramePoints( this->tetrahedra[t])){
				fprintf(stderr, "[%i] => %i\n",this->tetrahedra[t]->index, this->tetrahedra[t]->countNeighborTetrahedra);
				exit(0);
			}
		}*/
	}
	fprintf( stderr, "\n");
}
/*****************************************************************************/


template <int Dimensions, class VertexType, class SimplexType>
void Triangulation<Dimensions,VertexType,SimplexType>::sortTetrahedra()
{
	int i,ii;

	// sort vertex order of tets
	for( i=0; i<this->countTetrahedra; i++){
		for( int v=0; v<Dimensions; v++)
			for( int vv=0; vv<Dimensions-v; vv++)
				if( abs(this->tetrahedra[i]->vertices[vv]->index) > abs(this->tetrahedra[i]->vertices[vv+1]->index))
				{
					VertexType *temp = this->tetrahedra[i]->vertices[vv];
					this->tetrahedra[i]->vertices[vv] = this->tetrahedra[i]->vertices[vv+1];
					this->tetrahedra[i]->vertices[vv+1] = temp;
				}

		/*for( int v=0; v<Dimensions+1; v++)
			fprintf( stderr, "%i ", this->tetrahedra[i]->vertices[v]->index);
		fprintf( stderr, "\n");*/
	}

	// sort tets
	for( i=0; i<this->countTetrahedra-1; i++)
		for( ii=0; ii<this->countTetrahedra-i-1; ii++)
			if( this->tetrahedra[ii]->vertices[0]->index <= this->tetrahedra[ii+1]->vertices[0]->index &&
				this->tetrahedra[ii]->vertices[1]->index <= this->tetrahedra[ii+1]->vertices[1]->index &&
				this->tetrahedra[ii]->vertices[2]->index <= this->tetrahedra[ii+1]->vertices[2]->index &&
				this->tetrahedra[ii]->vertices[3]->index <= this->tetrahedra[ii+1]->vertices[3]->index ){

				SimplexType *temp      = this->tetrahedra[ii];
				this->tetrahedra[ii]   = this->tetrahedra[ii+1];
				this->tetrahedra[ii+1] = temp;
			}

	for( i=0; i<this->countTetrahedra; i++){
		for( int v=0; v<Dimensions+1; v++)
			fprintf( stderr, "%i ", this->tetrahedra[i]->vertices[v]->index);
		fprintf( stderr, "\n");
	}

}
/*****************************************************************************/


template <int Dimensions, class VertexType, class SimplexType>
bool Triangulation<Dimensions,VertexType,SimplexType>::SimplexContainsFramePoints( SimplexType* tet)
{
	int NR_Simplex_POINTS=Dimensions+1;

	for( int v=0; v<NR_Simplex_POINTS; v++){
		if( isElementOf( this->framePoints, this->countFramePoints, tet->vertices[v]))
			return true;
	}
	return false;
}
/*****************************************************************************/


template <int Dimensions, class VertexType,class SimplexType>
inline void Triangulation<Dimensions,VertexType,SimplexType>::setFramePoints3D()
{
	Triangulation<Dimensions,VertexType,SimplexType> *triangulation = this;
	int NR_Simplex_POINTS = 3+1;

	//double dimMin[3], dimMax[3];
	int i, dim;

	if( triangulation->countVertices==0)
		return;

	/*for( dim=0; dim<3; dim++){
		dimMin[dim] = dimMax[dim] = triangulation->vertices[0]->position[dim] ;
	}

	// find border
	for( i=1; i<triangulation->countVertices; i++){
		for( dim=0; dim<3; dim++){
			if( dimMin[dim] > triangulation->vertices[i]->position[dim])
				dimMin[dim] = triangulation->vertices[i]->position[dim] ;
			if( dimMax[dim] < triangulation->vertices[i]->position[dim])
				dimMax[dim] = triangulation->vertices[i]->position[dim] ;
		}
	}*/
	triangulation->setDomain();

	// BACKUP CELL LIST
	//Vertex<3> ** tempList = triangulation->vertices;
	//int tempCountCells = triangulation->countVertices;

	// NEW EMPTY & EXTENDED CELL LIST
	triangulation->framePoints = (VertexType**) malloc( sizeof(VertexType*) * ((int)pow( 2, 3)));
	assert(triangulation->framePoints);
	triangulation->countFramePoints = 0;
	//triangulation->maxFramePoints   = tempCountCells + (int)pow( 2, 3);

	// SET INITIAL FRAME POINTS

	double strech = MAX( triangulation->domain->max(0) - triangulation->domain->min(0),
			MAX( triangulation->domain->max(1) - triangulation->domain->min(1),
					triangulation->domain->max(2) - triangulation->domain->min(2)));

	for( i=0; i<(int)pow( 2, 3); i++){
		int temp = i;
		int countOnes = 0;
		triangulation->framePoints[triangulation->countFramePoints] = new VertexType( 0., 0., 0.);
		for( dim=0; dim<3; dim++){
			countOnes += temp%2;
			//printf( "%lf ", (1 - temp%2)*xMin[dim] + (temp%2)*xMax[dim] + ((temp%2 - 0.5)));
			//printf( "%lf ", (1 - temp%2)*xMin[dim] + (temp%2)*xMax[dim] + (2*(temp%2)-1)*myRand()*1.);

			triangulation->framePoints[triangulation->countFramePoints]->pos()[dim] = (1 - temp%2)*triangulation->domain->min(dim) + (temp%2)*triangulation->domain->max(dim) + (2*(temp%2)-1)*myRand()*strech;
			temp = temp/2;
		}
		//fprintf( stderr, " (%i)", countOnes);
		triangulation->framePoints[triangulation->countFramePoints]->setIndex(-1-i);

		temp = i;
		if( countOnes%2 == 1)
		for( dim=0; dim<3; dim++){
			triangulation->framePoints[triangulation->countFramePoints]->pos()[dim] += (myRand()*strech*0.1+strech)*(2*(temp%2) - 1);
			temp = temp/2;
		}

		//for( dim=0; dim<3; dim++)
			//printf( "%lf ", triangulation->framePoints[triangulation->countFramePoints]->position[dim]);
		//printf( "\n");

		triangulation->countFramePoints++;
	}

	// SET INITIAL TETS
	triangulation->countTetrahedra = 0;
	triangulation->maxTetrahedra = 1 + 4 + 6;
	triangulation->tetrahedra = (SimplexType**) realloc( triangulation->tetrahedra, sizeof(SimplexType*)*triangulation->maxTetrahedra);
	assert(triangulation->tetrahedra);

	// central tet
	triangulation->tetrahedra[triangulation->countTetrahedra] = new SimplexType();
	for( i=0; i<NR_Simplex_POINTS; i++){
		int temp = 0;
		for( dim=0; dim<3; dim++){
			temp = temp*2;
			if(dim==i || 3==i)
				temp++;
		}
		triangulation->tetrahedra[triangulation->countTetrahedra]->vertices[i] = triangulation->framePoints[triangulation->countFramePoints-temp-1];
		//fprintf( stderr, "%i ", 7-temp);
	}
	triangulation->tetrahedra[triangulation->countTetrahedra]->index = triangulation->countTetrahedra;
	//fprintf( stderr, "\n");

	// 4 surrounding tets
	for( i=0; i<NR_Simplex_POINTS; i++){

		// first tet vertex
		int temp = 0;
		int temp_others[3];
		for( dim=0; dim<3; dim++){
			temp = temp*2;
			if(dim!=i && 3!=i){
				temp++;
				temp_others[dim]=1;
			}else{
				temp_others[dim]=0;
			}
		}

		triangulation->tetrahedra[triangulation->countTetrahedra+1+i] = new SimplexType();

		//triangulation->tetrahedra[triangulation->countTetrahedra+1+i]->addSimplexNeighbor( triangulation->tetrahedra[triangulation->countTetrahedra]    );
		//triangulation->tetrahedra[triangulation->countTetrahedra    ]->addSimplexNeighbor( triangulation->tetrahedra[triangulation->countTetrahedra+1+i]);

		triangulation->tetrahedra[triangulation->countTetrahedra+1+i]->vertices[0] = triangulation->framePoints[triangulation->countFramePoints-temp-1];
		//fprintf( stderr, "%i ", 7-temp);

		// following 3 vertices
		for( dim=0; dim<3; dim++){
			int dim2;
			int temp = 0;
			for( dim2=0; dim2<3; dim2++){
				temp = temp*2;
				if(dim!=dim2)
					temp += temp_others[dim2];
				else
					temp += 1 - temp_others[dim2];
			}
			triangulation->tetrahedra[triangulation->countTetrahedra+1+i]->vertices[1+dim] = triangulation->framePoints[triangulation->countFramePoints-temp-1];
			//fprintf( stderr, "%i ", triangulation->tetrahedra[triangulation->countTetrahedra+1+i]->vertices[1+dim]->index); //7-temp);
		}
		//fprintf( stderr, "\n");
		triangulation->tetrahedra[triangulation->countTetrahedra+1+i]->index = triangulation->countTetrahedra+1+i;
	}
	triangulation->countTetrahedra += 5;

	// 6 peripheric tets
	for( dim=0; dim<3; dim++){ // for all 3
		for( i=0; i<2; i++){        // for each direction
			//fprintf( stderr, "tet %i:", triangulation->countTetrahedra+dim*2+i);
			int ii;
			triangulation->tetrahedra[triangulation->countTetrahedra+dim*2+i] = new SimplexType();
			triangulation->tetrahedra[triangulation->countTetrahedra+dim*2+i]->index = triangulation->countTetrahedra+dim*2+i;
			for( ii=0; ii<4; ii++){ // for all tet points
				int temp2 = 0;
				int temp = ii;
				int dim2;
				for( dim2=0; dim2<3; dim2++){
					temp2 = temp2*2;
					if(dim==dim2)
						temp2 += i;
					else{
						temp2 += temp%2;
						temp = temp/2;
					}
				}
				triangulation->tetrahedra[triangulation->countTetrahedra+dim*2+i]->vertices[ii] = triangulation->framePoints[temp2];
				//fprintf( stderr, "%i ", triangulation->tetrahedra[triangulation->countTetrahedra+dim*2+i]->vertices[ii]->index);//temp2);

			}
			//fprintf( stderr, "\n");

		}

	}
	triangulation->countTetrahedra += 6;
	triangulation->setSimplexNeighbors( );

	// SET NEIGHBORHOOD RELATIONS
	int ii, iii;
	for( i=0; i<triangulation->countTetrahedra; i++){
		//fprintf( stderr, "%i \n", triangulation->tetrahedra[i]->index);

		for( ii=0; ii<NR_Simplex_POINTS; ii++){
			for( iii=ii+1; iii<NR_Simplex_POINTS; iii++){
				//fprintf( stderr, "add %i & %i \n", triangulation->tetrahedra[i]->vertices[ii]->index, triangulation->tetrahedra[i]->vertices[iii]->index);
				triangulation->tetrahedra[i]->vertices[ii]->addNeighbor( triangulation->tetrahedra[i]->vertices[iii]);
				triangulation->tetrahedra[i]->vertices[iii]->addNeighbor( triangulation->tetrahedra[i]->vertices[ii]);
			}
		}
	}

	//addNeighbor(

	// TEST
	//for( i=0; i<tempCountCells; i++)

	// PRINT ALL

	/*fprintf( stderr, "countVertices:%i, maxVertices:%i\n", countVertices, maxVertices);
	for( i=0; i<this->countVertices; i++)
		fprintf( stderr, "%i: index=%i, (%lf %lf %lf)\n",
		         i, this->vertices[i]->index,
		         this->vertices[i]->position[0], this->vertices[i]->position[1], this->vertices[i]->position[2]);

	fprintf( stderr, "countFramePoints:%i, maxFramePoints:%i\n", countFramePoints, countFramePoints);
	for( i=0; i<this->countFramePoints; i++)
		fprintf( stderr, "%i: index=%i, (%lf %lf %lf)\n",
		         i, this->framePoints[i]->index,
		         this->framePoints[i]->position[0], this->framePoints[i]->position[1], this->framePoints[i]->position[2]);

	fprintf( stderr, "countTetrahedra:%i, maxTetrahedra:%i\n", countTetrahedra, maxTetrahedra);
	for( i=0; i<this->countTetrahedra; i++){
		fprintf( stderr, "%i: index=%i, ", i, this->tetrahedra[i]->index);
		fprintf( stderr, "(%i %i %i %i)",
		         this->tetrahedra[i]->vertices[0]->index, this->tetrahedra[i]->vertices[1]->index, this->tetrahedra[i]->vertices[2]->index, this->tetrahedra[i]->vertices[3]->index);
		fprintf( stderr, " => %i neighbors: (", this->tetrahedra[i]->countNeighborTetrahedra);
		for( int ii=0; ii<this->tetrahedra[i]->countNeighborTetrahedra; ii++)
			fprintf( stderr, " %i", this->tetrahedra[i]->neighborTetrahedra[ii]->index);
		fprintf( stderr, " )\n");
	}*/
	// TRIANGULATION
	//Vertex<3> ** tempList = this->vertices;
	//int tempCountCells = this->countVertices;
	/*printToPovray( "test0.pov", this, NULL, 0);
	checkDelaunayCondition( this, NULL, 0);

	Triangulation<3> *testDiagram = Triangulation<3>::newVoronoiDiagramFromFile( "data");
	fprintf( stderr, "countVertices:%i, maxVertices:%i\n", countVertices, maxVertices);
	for( i=0; i<testDiagram->countVertices; i++)
		fprintf( stderr, "%i: index=%i, (%lf %lf %lf)\n",
		         i, testDiagram->vertices[i]->index,
		         testDiagram->vertices[i]->position[0], testDiagram->vertices[i]->position[1], testDiagram->vertices[i]->position[2]);

	fprintf( stderr, "countTetrahedra:%i, maxTetrahedra:%i\n", countTetrahedra, maxTetrahedra);
	for( i=0; i<testDiagram->countTetrahedra; i++){
		fprintf( stderr, "%i: index=%i, ", i, testDiagram->tetrahedra[i]->index);
		fprintf( stderr, "(%i %i %i %i)",
		         testDiagram->tetrahedra[i]->vertices[0]->index, testDiagram->tetrahedra[i]->vertices[1]->index, testDiagram->tetrahedra[i]->vertices[2]->index, testDiagram->tetrahedra[i]->vertices[3]->index);
		fprintf( stderr, " => %i neighbors: (", testDiagram->tetrahedra[i]->countNeighborTetrahedra);
		for( int ii=0; ii<testDiagram->tetrahedra[i]->countNeighborTetrahedra; ii++)
			fprintf( stderr, " %i", testDiagram->tetrahedra[i]->neighborTetrahedra[ii]->index);
		fprintf( stderr, " )\n");
	}

	char filename[100];
	for( i=0; i<tempCountCells; i++){
		fprintf( stderr, "Add point:%i  => (%lf %lf %lf)\n", tempList[i]->index = testDiagram->countVertices, tempList[i]->position[0], tempList[i]->position[1], tempList[i]->position[2]);
		insertVoronoiCell( testDiagram, tempList[i]);
		sprintf( filename, "anim_%i.pov", i+1);
		printToPovray( filename, testDiagram, tempList, i+1);
		//checkDelaunayCondition( testDiagram, NULL, 0);
	}	*/

	// SET INDEX TO NEGATIVE VALUES
	//for( int i=0; i<this->countFramePoints; i++)
	//	this->framePoints[i]->index = -i-1;
}

template <int Dimensions, class VertexType,class SimplexType>
inline void Triangulation<Dimensions,VertexType,SimplexType>::setFramePoints2D()
{
	//Triangulation<Dimensions,VertexType,SimplexType> *triangulation = this;
	int NR_Simplex_POINTS = 2+1;

	this->setDomain();

	// BACKUP CELL LIST
	//Vertex<2> ** tempList = this->vertices;
	//int tempCountCells = this->countVertices;

	// NEW EMPTY & EXTENDED CELL LIST
	this->framePoints = (VertexType**) malloc( sizeof(VertexType*) * ((int)pow( 2, 2)));
	assert(this->framePoints);
	this->countFramePoints = 0;
	//this->maxFramePoints   = tempCountCells + (int)pow( 2, 2);

	//fprintf( stderr, " min frame:(%f, %f) - (%f, %f)\n", xMin[0], xMin[1], xMax[0], xMax[1]);

	// SET INITIAL FRAME POINTS
	for( int i=0; i<(int)pow( 2, 2); i++){
		int temp = i;
		int countOnes = 0;
		this->framePoints[this->countFramePoints] = new VertexType();
		this->framePoints[this->countFramePoints]->setx( 0.);
		this->framePoints[this->countFramePoints]->sety( 0.);
		//fprintf( stderr, " (%i, %i) [%p]", i%2, i/2, this->framePoints[this->countFramePoints]);
		for( int dim=0; dim<2; dim++){
			countOnes += temp%2;
			//fprintf( stderr, " %i,", temp%2);
			//fprintf( stderr, " %f+%f,",
			//		(1 - temp%2)*xMin[dim] + (temp%2)*xMax[dim],
			//		( i%2 == i/2 ? 2.*(temp%2)-1. : (2.*(temp%2)-1.)*(1+1) ));
			//printf( "%lf ", (1 - temp%2)*xMin[dim] + (temp%2)*xMax[dim] + ((temp%2 - 0.5)));
			//printf( "%lf ", (1 - temp%2)*xMin[dim] + (temp%2)*xMax[dim] + (2*(temp%2)-1)*myRand()*1.);
			this->framePoints[this->countFramePoints]->pos()[dim] = (1 - temp%2)*domain->min(dim) + (temp%2)*domain->max(dim) + ( i%2 != i/2 ? (2.*(temp%2)-1.)*myRand() : (2.*(temp%2)-1.)*(myRand()+1));
			temp = temp/2;
		}
		//fprintf( stderr, " (%i)", countOnes);
		//fprintf( stderr, " (%i, %i)", countOnes%2, temp);
		this->framePoints[this->countFramePoints]->setIndex(-1-i);

		/*printf( " -> point %i: ", i);
		for( int dim=0; dim<2; dim++)
			printf( "%lf ", this->framePoints[this->countFramePoints]->position[dim]);
		printf( "\n");*/
		//printf( "%i %i\n", i%2, i%3);

		this->countFramePoints++;
	}



	// SET INITIAL TRIANGLES
	this->countTetrahedra = 0;
	this->maxTetrahedra = 2;
	this->tetrahedra = (SimplexType**) realloc( this->tetrahedra, sizeof(SimplexType*)*this->maxTetrahedra);
	assert(this->tetrahedra);

	for( int t=0; t<2; t++){
		// alloc triangle
		this->tetrahedra[this->countTetrahedra] = new SimplexType();

		// set vertices
		for( int i=0; i<NR_Simplex_POINTS; i++){
			this->tetrahedra[this->countTetrahedra]->setVertex( this->framePoints[t+i], i);
			//fprintf( stderr, "%i ", 7-temp);
		}

		// set index
		this->tetrahedra[this->countTetrahedra]->setIndex( this->countTetrahedra);
		//fprintf( stderr, "\n");

		this->countTetrahedra++;
	}

	this->setSimplexNeighbors( );

	// SET NEIGHBORHOOD RELATIONS
	int ii, iii;
	for( int i=0; i<this->countTetrahedra; i++){
		//fprintf( stderr, "%i \n", this->tetrahedra[i]->index);

		for( ii=0; ii<NR_Simplex_POINTS; ii++){
			for( iii=ii+1; iii<NR_Simplex_POINTS; iii++){
				//fprintf( stderr, "add %i & %i \n", this->tetrahedra[i]->vertices[ii]->index, this->tetrahedra[i]->vertices[iii]->index);
				this->tetrahedra[i]->vertices[ii]->addNeighbor( this->tetrahedra[i]->vertices[iii]);
				this->tetrahedra[i]->vertices[iii]->addNeighbor( this->tetrahedra[i]->vertices[ii]);
			}
		}
	}

	/*this->countTetrahedra = 0;
	this->maxTetrahedra = 2;
	this->tetrahedra = (Simplex<2>**) realloc( this->tetrahedra, sizeof(Simplex<2>*)*this->maxTetrahedra);

	// central tet
	this->tetrahedra[this->countTetrahedra] = new Simplex<2>();
	for( i=0; i<NR_Simplex_POINTS; i++){
		int temp = 0;
		for( dim=0; dim<2; dim++){
			temp = temp*2;
			if(dim==i || 2==i)
				temp++;
		}
		this->tetrahedra[this->countTetrahedra]->vertices[i] = this->framePoints[this->countFramePoints-temp-1];
		//fprintf( stderr, "%i ", 7-temp);
	}
	this->tetrahedra[this->countTetrahedra]->index = this->countTetrahedra;
	//fprintf( stderr, "\n");

	// 4 surrounding tets
	for( i=0; i<NR_Simplex_POINTS; i++){

		// first tet vertex
		int temp = 0;
		int temp_others[2];
		for( dim=0; dim<2; dim++){
			temp = temp*2;
			if(dim!=i && 2!=i){
				temp++;
				temp_others[dim]=1;
			}else{
				temp_others[dim]=0;
			}
		}

		this->tetrahedra[this->countTetrahedra+1+i] = new Simplex<2>();

		//this->tetrahedra[this->countTetrahedra+1+i]->addSimplexNeighbor( this->tetrahedra[this->countTetrahedra]    );
		//this->tetrahedra[this->countTetrahedra    ]->addSimplexNeighbor( this->tetrahedra[this->countTetrahedra+1+i]);

		this->tetrahedra[this->countTetrahedra+1+i]->vertices[0] = this->framePoints[this->countFramePoints-temp-1];
		//fprintf( stderr, "%i ", 7-temp);

		// following 3 vertices
		for( dim=0; dim<2; dim++){
			int dim2;
			int temp = 0;
			for( dim2=0; dim2<2; dim2++){
				temp = temp*2;
				if(dim!=dim2)
					temp += temp_others[dim2];
				else
					temp += 1 - temp_others[dim2];
			}
			this->tetrahedra[this->countTetrahedra+1+i]->vertices[1+dim] = this->framePoints[this->countFramePoints-temp-1];
			//fprintf( stderr, "%i ", this->tetrahedra[this->countTetrahedra+1+i]->vertices[1+dim]->index); //7-temp);
		}
		//fprintf( stderr, "\n");
		this->tetrahedra[this->countTetrahedra+1+i]->index = this->countTetrahedra+1+i;
	}
	this->countTetrahedra += 5;

	// 6 peripheric tets
	for( dim=0; dim<2; dim++){ // for all 2
		for( i=0; i<2; i++){        // for each direction
			//fprintf( stderr, "tet %i:", this->countTetrahedra+dim*2+i);
			int ii;
			this->tetrahedra[this->countTetrahedra+dim*2+i] = new Simplex<2>();
			this->tetrahedra[this->countTetrahedra+dim*2+i]->index = this->countTetrahedra+dim*2+i;
			for( ii=0; ii<4; ii++){ // for all tet points
				int temp2 = 0;
				int temp = ii;
				int dim2;
				for( dim2=0; dim2<2; dim2++){
					temp2 = temp2*2;
					if(dim==dim2)
						temp2 += i;
					else{
						temp2 += temp%2;
						temp = temp/2;
					}
				}
				this->tetrahedra[this->countTetrahedra+dim*2+i]->vertices[ii] = this->framePoints[temp2];
				//fprintf( stderr, "%i ", this->tetrahedra[this->countTetrahedra+dim*2+i]->vertices[ii]->index);//temp2);

			}
			//fprintf( stderr, "\n");

		}

	}
	this->countTetrahedra += 6;
	setSimplexNeighbors( this);

	// SET NEIGHBORHOOD RELATIONS
	int ii, iii;
	for( i=0; i<this->countTetrahedra; i++){
		//fprintf( stderr, "%i \n", this->tetrahedra[i]->index);

		for( ii=0; ii<NR_Simplex_POINTS; ii++){
			for( iii=ii+1; iii<NR_Simplex_POINTS; iii++){
				//fprintf( stderr, "add %i & %i \n", this->tetrahedra[i]->vertices[ii]->index, this->tetrahedra[i]->vertices[iii]->index);
				addNeighbor( this->tetrahedra[i]->vertices[ii], this->tetrahedra[i]->vertices[iii]);
				addNeighbor( this->tetrahedra[i]->vertices[iii], this->tetrahedra[i]->vertices[ii]);
			}
		}
	}*/

}

template <int Dimensions,class VertexType,class SimplexType>
inline void Triangulation<Dimensions,VertexType,SimplexType>::setFramePoints1D()
{
	Triangulation<Dimensions,VertexType,SimplexType> *triangulation = this;
	//int NR_Simplex_POINTS = 1+1;

	triangulation->setDomain();

	// NEW EMPTY & EXTENDED CELL LIST
	triangulation->framePoints = (VertexType**) malloc( sizeof(VertexType*) * 2);
	triangulation->countFramePoints = 2;

	// SET INITIAL FRAME POINTS
	for( int i=0; i<1+1; i++){
		triangulation->framePoints[i] = new VertexType();
		triangulation->framePoints[i]->x(0);
		triangulation->framePoints[i]->setIndex(-i-1);
		triangulation->framePoints[i]->pos()[0] = (i? triangulation->domain->min(0)-myRand() : triangulation->domain->max(0)+myRand());
	}

	// SET FRAME POINTS NEIGHBORHOOD RELATIONS
	triangulation->framePoints[0]->addNeighbor( triangulation->framePoints[1]);
	triangulation->framePoints[1]->addNeighbor( triangulation->framePoints[0]);

	// SET INITIAL TRIANGLES
	triangulation->countTetrahedra = 1;
	triangulation->maxTetrahedra = 1;
	triangulation->tetrahedra = (SimplexType**) realloc( triangulation->tetrahedra, sizeof(SimplexType*)*triangulation->maxTetrahedra);
	triangulation->tetrahedra[0] = new SimplexType();
	triangulation->tetrahedra[0]->setVertex( triangulation->framePoints[0], 0);
	triangulation->tetrahedra[0]->setVertex( triangulation->framePoints[1], 1);
	triangulation->tetrahedra[0]->setIndex(0);

	// SET TRIANGLES NEIGHBORHOOD RELATIONS
	// no neighbors in 1D

}
/*****************************************************************************/


template <int Dimensions, class VertexType,class SimplexType>
void Triangulation<Dimensions,VertexType,SimplexType>::setFramePoints()
{
	switch( Dimensions){
	case 1:
		this->setFramePoints1D();
		break;
	case 2:
		this->setFramePoints2D();
		break;
	case 3:
		this->setFramePoints3D();
		break;
	}
}



/*void Triangulation<Dimensions>::coarsen( Vertex<Dimensions> *vc, int scale)
{
	double base_x = floor(vc->position[0]);
	double base_y = floor(vc->position[1]);
#if Dimensions == 3
	double base_z = floor(vc->position[2]);
#endif

	fprintf( stderr, "COARSEN( cell: %i, scale: %i)\n", vc->index, scale);
	assert( vc->refined);
	assert( vc->countNeighbors>0);


	// SEARCH & DELETE CONSIDERED POINTS
	int actualStackElement = 0;
	int stackSize = 0;
	Vertex<Dimensions> *stack[200];
	stack[stackSize] = vc; stackSize++;
	do{
		for( int n=0; n<stack[actualStackElement]->countNeighbors; n++){
			if( stack[actualStackElement]->neighbors[n]->position[0] > base_x &&
				stack[actualStackElement]->neighbors[n]->position[0] < base_x+1 &&
				stack[actualStackElement]->neighbors[n]->position[1] > base_y &&
				stack[actualStackElement]->neighbors[n]->position[1] < base_y+1 &&
				!isElementOf( stack, stackSize, stack[actualStackElement]->neighbors[n])){

				stack[stackSize] = stack[actualStackElement]->neighbors[n];
				stackSize++;
			}
		}

		removeVoronoiCell( this, stack[actualStackElement]);

		actualStackElement++;
	}while( actualStackElement<stackSize);
	//fprintf( stderr, "found %i cells to fuse\n", stackSize);


	// REINSERT COARSE POINT INTO TRIANGULATION
#if Dimensions == 3
	int index = base_x*this->xN[1]*this->xN[2] + base_y*this->xN[2] + base_z;
#elif Dimensions == 2
	int index = base_x*this->xN[1] + base_y;
#endif
	triangulateVoronoiCell( this, this->vertices[index]);
	this->vertices[index]->refined = false;


	// ACTUALIZE ALL CONSIDERED VERTICES
	//actualStackElement = 0;
	//stackSize = 0;
	//stack[stackSize] = this->vertices[index]; stackSize++;

/ *
	// mark new cell to be actualized
	this->vertices[index]->countFreeNeighborCells = -1;

	// mark direct neighbors to be actualized
	for( int n=0; n<this->vertices[index]->countNeighbors; n++ )
		this->vertices[index]->neighbors[n]->countFreeNeighborCells = -1;

	// actualize indirect neighbors
	for( int n=0; n<this->vertices[index]->countNeighbors; n++ )
		// direct neighbors
		for( int nn=0; nn<this->vertices[index]->neighbors[n]->countNeighbors; nn++ )
			// indirect neighbors
			if( this->vertices[index]->neighbors[n]->neighbors[nn]->countFreeNeighborCells != -1)
			{
				int countFree = 0;
				for( int nnn=0; nnn<this->vertices[index]->neighbors[n]->neighbors[nn]->countNeighbors; nnn++ )
					if( this->vertices[index]->neighbors[n]->neighbors[nn]->neighbors[nnn]->isFree())
						countFree++;
				this->vertices[index]->neighbors[n]->neighbors[nn]->countFreeNeighborCells = countFree;
				fprintf( stderr, "-> update free neighbors of %ith point\n", this->vertices[index]->neighbors[n]->neighbors[nn]->index);
			}

	// actualize direct neighbors
	for( int n=0; n<this->vertices[index]->countNeighbors; n++ ){
		// direct neighbors
		int countFree = 0;
		for( int nn=0; nn<this->vertices[index]->neighbors[n]->countNeighbors; nn++ )
			if( this->vertices[index]->neighbors[n]->neighbors[nn]->isFree())
				countFree++;
		this->vertices[index]->neighbors[n]->countFreeNeighborCells = countFree;
		fprintf( stderr, "-> update free neighbors of %ith point\n", this->vertices[index]->neighbors[n]->index);
	}

	// actualize new cell
	int countFree = 0;
	for( int n=0; n<this->vertices[index]->countNeighbors; n++ )
		if( this->vertices[index]->neighbors[n]->isFree())
			countFree++;
	this->vertices[index]->countFreeNeighborCells = countFree;
	fprintf( stderr, "-> update free neighbors of %ith point\n", this->vertices[index]->index);
* /

	// mark considered cells to be actualized
	for( int n=0; n<this->vertices[index]->countNeighbors; n++ )
		// direct neighbors
		for( int nn=0; nn<this->vertices[index]->neighbors[n]->countNeighbors; nn++ )
			// indirect neighbors
			this->vertices[index]->neighbors[n]->neighbors[nn]->countFreeNeighborCells = -1;

	// actualize marked cells
	for( int n=0; n<this->vertices[index]->countNeighbors; n++ )
		// direct neighbors
		for( int nn=0; nn<this->vertices[index]->neighbors[n]->countNeighbors; nn++ )
			// indirect neighbors
			if( this->vertices[index]->neighbors[n]->neighbors[nn]->countFreeNeighborCells == -1)
			{
				int countFree = 0;
				for( int nnn=0; nnn<this->vertices[index]->neighbors[n]->neighbors[nn]->countNeighbors; nnn++ )
					if( this->vertices[index]->neighbors[n]->neighbors[nn]->neighbors[nnn]->isFree())
						countFree++;
				this->vertices[index]->neighbors[n]->neighbors[nn]->countFreeNeighborCells = countFree;
				fprintf( stderr, "-> update free neighbors of %ith point\n", this->vertices[index]->neighbors[n]->neighbors[nn]->index);
			}else{
				fprintf( stderr, "-> keep free neighbors of %ith point\n", this->vertices[index]->neighbors[n]->neighbors[nn]->index);

			}



	/ *for( int c=0; c<this->vertices[index]->countNeighbors; c++ ){
		// mark
		int countFree = 0;
		for( int n=0; n<this->vertices[index]->neighbors[c]->countNeighbors; n++)
			if( this->vertices[index]->neighbors[c]->neighbors[n]->isFree())
				countFree++;

		//fprintf( stderr, "-> %s free neighbors of %ith point in stack\n",
//				(stack[actualStackElement]->countFreeNeighborCells == countFree ? "keep" : "update"), actualStackElement+1);
		this->vertices[index]->neighbors[c]->countFreeNeighborCells = countFree;
	}* /

}*/

template <int Dimensions>
bool Vertex<Dimensions>::isDomainBorder( Triangulation<Dimensions> *vd)
{
	return this->position[0] < vd->xMin[0] + vd->boundaryThickness
	    && this->position[0] > vd->xMax[0] - vd->boundaryThickness
#if Dimensions >= 2
		&& this->position[1] < vd->xMin[1] + vd->boundaryThickness
		&& this->position[1] > vd->xMax[1] - vd->boundaryThickness
#endif
#if Dimensions >= 3
		&& this->position[2] < vd->xMin[2] + vd->boundaryThickness
		&& this->position[2] > vd->xMax[2] - vd->boundaryThickness
#endif
		;
}


/*template <int Dimensions, class VertexType, class SimplexType>
POSITION_T Triangulation<Dimensions,VertexType,SimplexType>::getVoronoiCellVolume(VertexType *v) {

	Simplex<Dimensions> *first = getSimplexContainingPoint( v);
	Simplex<Dimensions> *tet = first;
	Simplex<Dimensions> *last = first;

	double array1[Dimensions];
	double array2[Dimensions];
	double *cc1 = array1, *cc2 = array2, *cctmp;
	tet->getCircumCenter( cc2);

	double volume = 0.;

	do {
		//fprintf( stderr, "%i -> ", tet->index);

		int n=0;
		for( n=0; n<tet->countNeighborTetrahedra
			&& (!tet->neighborTetrahedra[n]->contains( v) || tet->neighborTetrahedra[n]==last); n++) ;

		last = tet;
		tet = tet->neighborTetrahedra[n];

		cctmp = cc1;
		cc1 = cc2;
		cc2 = cctmp;

		tet->getCircumCenter( cc2);

		// A = a*h_a
		double a = sqrt( pow(cc1[0] - cc2[0], 2)
		               + pow(cc1[1] - cc2[1], 2));
		double b = sqrt( pow(v->position[0] - cc2[0], 2)
		               + pow(v->position[1] - cc2[1], 2));
		double c = sqrt( pow(cc1[0] - v->position[0], 2)
				       + pow(cc1[1] - v->position[1], 2));
		double s = (a+b+c)/2.;
		volume += sqrt(s*(s-a)*(s-b)*(s-c));

	} while (tet != first);

	//fprintf( stderr, "%lf -> ", volume);
	return volume;
}*/
/*****************************************************************************/


template <int Dimensions, class VertexType>
bool Simplex<Dimensions,VertexType>::contains( Vertex<Dimensions> * vc) {
	for( int i=0; i<Dimensions+1; i++)
		if( this->vertices[i]==vc)
			return true;
	return false;
}
/*****************************************************************************/


/*template <int Dimensions, class VertexType>
bool Simplex<Dimensions,VertexType>::contains( VertexType * vc) {
	for( int i=0; i<Dimensions+1; i++)
		if( this->vertices[i]==vc)
			return true;
	return false;
}*/
/*****************************************************************************/


template <int Dimensions, class VertexType, class SimplexType>
void Triangulation<Dimensions, VertexType, SimplexType>::printVoronoiDiagramToPovray( const char * filename, Domain<Dimensions>* subDomain)
{
	switch( Dimensions)
	{
	case 2:
		this->printVoronoiDiagramToPovray2D(filename, subDomain);
		break;
	case 3:
		this->printVoronoiDiagramToPovray3D(filename, subDomain);
		break;
	}
}
/*****************************************************************************/


template <int Dimensions, class VertexType, class SimplexType>
void Triangulation<Dimensions, VertexType, SimplexType>::printVoronoiDiagramToPovray3D( const char * filename, Domain<Dimensions>* subDomain)
{
	if( subDomain==0)
		subDomain=domain;
	Triangulation<Dimensions, VertexType, SimplexType>* vd = this;
	std::fstream fs;
	fs.open( filename, std::fstream::out);
	PovrayIO::writePovrayHeader( &fs, domain->min(0), domain->min(1), domain->max(0), domain->max(1), 2*(domain->max(0)-domain->min(0)));

	//PovrayIO::beginIntersectionWithBox( &fs, vd->xMin[0], vd->xMin[1], vd->xMin[2], vd->xMax[0], vd->xMax[1], vd->xMax[2]);
	// Construction Points
/*	for( int v=0; v<vd->countVertices; v++){
		PovrayIO::writeSphere( &fs, vd->vertices[v]->position[0], vd->vertices[v]->position[1], vd->vertices[v]->position[2], 0.04, "rgb <1,0,0>");
	}

	// Voronoi Cells & Neighborships
	for( int t=0; t<vd->countTetrahedra; t++)
	//if(!vd->tetrahedronContainsFramePoints( vd->tetrahedra[t]))
	{
		// Voronoi Cell Vertice
		double pos[DIMENSIONS];
		getCircumCenter( vd->tetrahedra[t], pos);
		PovrayIO::writeSphere( &fs, pos[0], pos[1], pos[2], 0.02, "rgb <0,1,0>");

		// Neighbors
		for( int nt=0; nt<vd->tetrahedra[t]->countNeighborTetrahedra; nt++)
		if(	//!vd->tetrahedronContainsFramePoints( vd->tetrahedra[t]->neighborTetrahedra[nt]) &&
				vd->tetrahedra[t]->index < vd->tetrahedra[t]->neighborTetrahedra[nt]->index){
			double npos[DIMENSIONS];
			getCircumCenter( vd->tetrahedra[t]->neighborTetrahedra[nt], npos);
			PovrayIO::writeCylinder( &fs, pos[0], pos[1], pos[2], npos[0], npos[1], npos[2], 0.01, "rgb <0,0,1>");
		}
	}
*/
	// Faces
	//fprintf( stderr, "Faces!\n");
	for( int v=0; v<vd->countVertices; v++)
	//int v = 13+27;
	if( vd->vertices[v] != NULL &&
	    vd->vertices[v]->numberOfNeighbors()>0 &&
	    subDomain->contains( vd->vertices[v]) &&
	    // only occupied cells
	    //vd->vertices[v]->getState() != FREE &&
	    // cut
		(	(// !vd->vertices[v]->refined &&
				true//vd->vertices[v]->position[2] < 0.5*(vd->xMin[2]+vd->xMax[2])
				//&&
				//vd->vertices[v]->position[2] > 0.5*(vd->xMin[2]+vd->xMax[2])-1.
			)
			/*||
			(vd->vertices[v]->refined &&
				vd->vertices[v]->position[2] < 0.5*(vd->xMin[2]+vd->xMax[2]) &&
				vd->vertices[v]->position[2] > 0.5*(vd->xMin[2]+vd->xMax[2])-1./(pow(GetAgent(vd->vertices[v])->maxCellCount, 1./3.)))*/
		)

	)
	{
		char color[100];
		/*switch (vd->vertices[v]->getState()) {
		case FREE:
			sprintf(color, "rgb <1,1,1>");
			break;
		case ACTIVE:*/
			sprintf(color, "rgb <1,1,1>  transmit 0.5");
			//sprintf(color, "rgb <%lf,%lf,%lf> ", RAND01, RAND01, RAND01);
		/*	break;
		case NONACTIVE:
			sprintf(color, "rgb <0,1,0>");
			break;
		case COMPARTMENT:
			sprintf(color, "rgb <0,1,0>");
			break;
		}*/

		PovrayIO::writeSphere( &fs,
				vd->vertices[v]->pos()[0], vd->vertices[v]->pos()[1], vd->vertices[v]->pos()[2],
				0.2, "rgb <1,0,0>");

		// get tet containing Point v
		SimplexType *tet = vd->getSimplexContainingPoint( vd->vertices[v]);
		//fprintf( stderr, "Point %i -> tet %i (%i, %i, %i, %i)\n", v, tet->index, tet->vertices[0]->index, tet->vertices[1]->index, tet->vertices[2]->index, tet->vertices[3]->index);

		for( int nv=0; nv<vd->vertices[v]->numberOfNeighbors(); nv++)
		if( vd->vertices[v]->getNeighbor(nv)->getIndex() >= 0 /* no framepoints */	)
		//		&& vd->vertices[v]->index > vd->vertices[v]->neighbors[nv]->index)
		{
			//fprintf( stderr, "Construct Circulator for Point %i - Point %i\n", v, vd->vertices[v]->neighbors[nv]->index);

			// Circulator of Tets around Point-NeighborPoint-Axis
			SimplexType *tetCirculator[100];
			int tetCirculatorLength = 0;

			// get (first) tet containing Points v & nv
			SimplexType *first = tet;
			while( !first->contains( (VertexType*)vd->vertices[v]->getNeighbor(nv))){
				// chose random neighbor tet containing point v
				int temp;
				do{
					// chose random neighbor tet
					temp = (int)(myRandE((double)first->countNeighborTetrahedra));
				}while( ! first->neighborTetrahedra[temp]->contains( vd->vertices[v]));
				first = first->neighborTetrahedra[temp];
			}
			tetCirculator[0] = first;
			POSITION_T posFirst[3];
			POSITION_T x[3];
			POSITION_T y[3];
			POSITION_T z[3];
			tetCirculator[0]->getCircumCenter( posFirst); x[0]=posFirst[0]; y[0]=posFirst[1]; z[0]=posFirst[2];


			// get second
			int nt=0;
			for( ; nt<first->countNeighborTetrahedra &&
				!( first->neighborTetrahedra[nt]->contains( (VertexType*)vd->vertices[v]) &&  first->neighborTetrahedra[nt]->contains( (VertexType*)vd->vertices[v]->getNeighbor(nv))); nt++) ;
			tetCirculator[1] = first->neighborTetrahedra[nt];
			tetCirculatorLength = 2;

			POSITION_T pos[3];
			tetCirculator[1]->getCircumCenter( pos);
			//PovrayIO::writeCylinder( &fs, posFirst[0], posFirst[1], posFirst[2], pos[0], pos[1], pos[2], 0.01, "rgb <0,0,1>");


			// Circulate
			do{
				int nt=0;
				SimplexType *actual = tetCirculator[tetCirculatorLength-1];
				for( ; nt<actual->countNeighborTetrahedra &&
					!(  actual->neighborTetrahedra[nt] != tetCirculator[tetCirculatorLength-2] &&
						actual->neighborTetrahedra[nt]->contains( (VertexType*)vd->vertices[v]) &&
						actual->neighborTetrahedra[nt]->contains( (VertexType*)vd->vertices[v]->getNeighbor(nv))); nt++) ;
				tetCirculator[tetCirculatorLength] = actual->neighborTetrahedra[nt];
				tetCirculatorLength++;

				POSITION_T pos[3];
				tetCirculator[tetCirculatorLength-1]->getCircumCenter( pos);
				POSITION_T posLast[3];
				tetCirculator[tetCirculatorLength-2]->getCircumCenter( posLast);

				// WRITE EDGE
				//if( pow(pos[0]-posLast[0],2) + pow(pos[1]-posLast[1],2) + pow(pos[2]-posLast[2],2) > 0.01)
				//PovrayIO::writeCylinder( &fs, pos[0], pos[1], pos[2], posLast[0], posLast[1], posLast[2], 0.01, "rgb <0,0,1>");

				if( first != tetCirculator[tetCirculatorLength-1]){
					x[1]=pos[0]; y[1]=pos[1]; z[1]=pos[2];
					x[2]=posLast[0]; y[2]=posLast[1]; z[2]=posLast[2];

					// WRITE FACE
					POSITION_T interface = triangleInterface<Dimensions,POSITION_T>( pos, posLast, posFirst);
					if(	isnan(interface	) || isinf(interface)){
						fprintf( stderr, "ERROR: Polygon Interface = %lf\n", interface);
						//exit(0);
					}
					if( interface > 0.001)
					PovrayIO::writePolygon<POSITION_T>( &fs, 3, x, y, z, color);

					//PovrayIO::writePolygon( &fs, 3, x, y, z, "rgb <1,0,1>");
					//PovrayIO::writePolygon( &fs, 3, x, y, z, "rgb <1,0,1,0.7>");
					//PovrayIO::writePrism( &fs, 0.1, 0.1, 3, x, y, z, "rgb <1,0,1,0.7>");
				}
			}while( first != tetCirculator[tetCirculatorLength-1]);
			tetCirculatorLength--;

			/*double x[tetCirculatorLength];
			double y[tetCirculatorLength];
			double z[tetCirculatorLength];
			for( int i=0; i< tetCirculatorLength; i++){
				fprintf( stderr, "-> tet %i (%i, %i, %i, %i)\n", i, tetCirculator[i]->vertices[0]->index, tetCirculator[i]->vertices[1]->index, tetCirculator[i]->vertices[2]->index, tetCirculator[i]->vertices[3]->index);
				double pos[DIMENSIONS];
				getCircumCenter( tetCirculator[i], pos);
				x[i] = pos[0];
				y[i] = pos[1];
				z[i] = pos[2];
			}
			PovrayIO::writePolygon( &fs, tetCirculatorLength, x, y, z, "rgb <1,0,1>");*/
			//fs.close();
			//exit(0);
		}
	}
	/*for( int t=0; t<vd->countTetrahedra; t++){
		// Voronoi Cell Vertice
		double pos[DIMENSIONS];
		getCircumCenter( vd->tetrahedra[t], pos);
		PovrayIO::writeSphere( &fs, pos[0], pos[1], pos[2], 0.02, "rgb <0,1,0>");

		// Neighbors
		for( int nt=0; nt<vd->tetrahedra[t]->countNeighborTetrahedra; nt++)
		if(vd->tetrahedra[t]->index < vd->tetrahedra[t]->neighborTetrahedra[nt]->index){
			double npos[DIMENSIONS];
			getCircumCenter( vd->tetrahedra[t]->neighborTetrahedra[nt], npos);
			PovrayIO::writeCylinder( &fs, pos[0], pos[1], pos[2], npos[0], npos[1], npos[2], 0.01, "rgb <0,0,1>");
		}
	}*/
	//PovrayIO::endIntersectionWithBox( &fs);

	fs.close();
}
/*****************************************************************************/
template <int Dimensions, class VertexType>
inline Simplex<Dimensions,VertexType>* Simplex<Dimensions,VertexType>::getNeighbor( int i)
{
	return neighborTetrahedra[i];
}

template <int Dimensions, class VertexType, class SimplexType>
void Triangulation<Dimensions,VertexType,SimplexType>::printVoronoiDiagramToVTK( const char * filename, Domain<Dimensions>* subDomain)
{
	switch(Dimensions){
	case 3:
		printVoronoiDiagramToVTK3D( filename, subDomain);
		return;
	}
}
/*****************************************************************************/

template <int Dimensions, class VertexType, class SimplexType>
inline void Triangulation<Dimensions,VertexType,SimplexType>::printVoronoiDiagramToVTK3D( const char * filename, Domain<Dimensions>* subDomain)
{
	if( subDomain==0)
		subDomain=domain;
	Triangulation<Dimensions,VertexType,SimplexType>* vd = this;

	FILE *fp = fopen( filename, "w+");

	fprintf( fp, "# vtk DataFile Version 2.0\n");
	fprintf( fp, "Really cool data\n");
	fprintf( fp, "ASCII\n");
	fprintf( fp, "DATASET POLYDATA\n");

	//char pointString[2][100000] = {"\0", "\0"};
	int  pointCount=0;
	int  polygonCount=0;

	// Faces
	fprintf( stderr, "Faces!\n");
	for( int v=0; v<vd->countVertices; v++)
	if( subDomain->contains( vd->vertices[v]))
	{
		//char color[100];
		//sprintf(color, "rgb <1,1,1>  transmit 0.5");

		SimplexType *tet = vd->getSimplexContainingPoint( vd->vertices[v]);

		for( int nv=0; nv<vd->vertices[v]->numberOfNeighbors(); nv++)

		//if(vd->vertices[v]->getNeighbor(nv)->index >= 0 /* no framepoints */)
		//if(vd->vertices[v]->index < vd->vertices[v]->getNeighbor(nv)->index)
		//if( subDomain->contains( vd->vertices[v]))
		{
			fprintf( stderr, "Construct Circulator for Point %i - Point %i\n", v, vd->vertices[v]->getNeighbor(nv)->getIndex());
			polygonCount++;
			// Circulator of Tets around Point-NeighborPoint-Axis
			SimplexType *tetCirculator[100];
			int tetCirculatorLength = 0;

			// get (first) tet containing Points v & nv
			SimplexType *first = tet;
			while( !first->contains( (VertexType*)vd->vertices[v]->getNeighbor(nv))){
				// chose random neighbor tet containing point v
				int temp;
				do{
					// chose random neighbor tet
					temp = (int)(myRandE((double)first->countNeighborTetrahedra));
				}while( ! first->neighborTetrahedra[temp]->contains( vd->vertices[v]));
				first = first->neighborTetrahedra[temp];
			}
			tetCirculator[0] = first;
			POSITION_T posFirst[3];
			//POSITION_T x[3];
			//POSITION_T y[3];
			//POSITION_T z[3];
			tetCirculator[0]->getCircumCenter( posFirst); //x[0]=posFirst[0]; y[0]=posFirst[1]; z[0]=posFirst[2];


			// get second
			int nt=0;
			for( ; nt<first->countNeighborTetrahedra &&
				!( first->neighborTetrahedra[nt]->contains( (VertexType*)vd->vertices[v]) &&  first->neighborTetrahedra[nt]->contains( (VertexType*)vd->vertices[v]->getNeighbor(nv))); nt++) ;
			tetCirculator[1] = first->neighborTetrahedra[nt];
			tetCirculatorLength = 2;

			POSITION_T pos[3];
			tetCirculator[1]->getCircumCenter( pos);

			// WRITE POINT
			pointCount++;

			// Circulate
			do{
				int nt=0;
				SimplexType *actual = tetCirculator[tetCirculatorLength-1];
				for( ; nt<actual->countNeighborTetrahedra &&
					!(  actual->neighborTetrahedra[nt] != tetCirculator[tetCirculatorLength-2] &&
						actual->neighborTetrahedra[nt]->contains( (VertexType*)vd->vertices[v]) &&
						actual->neighborTetrahedra[nt]->contains( (VertexType*)vd->vertices[v]->getNeighbor(nv))); nt++) ;
				tetCirculator[tetCirculatorLength] = actual->neighborTetrahedra[nt];
				tetCirculatorLength++;

				POSITION_T pos[3];
				tetCirculator[tetCirculatorLength-1]->getCircumCenter( pos);
				POSITION_T posLast[3];
				tetCirculator[tetCirculatorLength-2]->getCircumCenter( posLast);

				// WRITE POINT
				pointCount++;

				// WRITE EDGE

				/*if( first != tetCirculator[tetCirculatorLength-1]){
					x[1]=pos[0]; y[1]=pos[1]; z[1]=pos[2];
					x[2]=posLast[0]; y[2]=posLast[1]; z[2]=posLast[2];

					// WRITE FACE

				}*/
			}while( first != tetCirculator[tetCirculatorLength-1]);
			tetCirculatorLength--;
		}
	}

	fprintf( fp, "POINTS %i float\n", pointCount);
	pointCount=0;

	// Faces
	fprintf( stderr, "Faces!\n");
	for( int v=0; v<vd->countVertices; v++)
	if( subDomain->contains( vd->vertices[v]))
	{
		//char color[100];
		//sprintf(color, "rgb <1,1,1>  transmit 0.5");

		SimplexType *tet = vd->getSimplexContainingPoint( vd->vertices[v]);

		for( int nv=0; nv<vd->vertices[v]->numberOfNeighbors(); nv++)
		//if(vd->vertices[v]->getNeighbor(nv)->index >= 0 /* no framepoints */)
		//if(vd->vertices[v]->index < vd->vertices[v]->getNeighbor(nv)->index)
		{
			//fprintf( stderr, "Construct Circulator for Point %i - Point %i\n", v, vd->vertices[v]->getNeighbor(nv)->index);

			// Circulator of Tets around Point-NeighborPoint-Axis
			SimplexType *tetCirculator[100];
			int tetCirculatorLength = 0;

			// get (first) tet containing Points v & nv
			SimplexType *first = tet;
			while( !first->contains( (VertexType*)vd->vertices[v]->getNeighbor(nv))){
				// chose random neighbor tet containing point v
				int temp;
				do{
					// chose random neighbor tet
					temp = (int)(myRandE((double)first->countNeighborTetrahedra));
				}while( ! first->neighborTetrahedra[temp]->contains( vd->vertices[v]));
				first = first->neighborTetrahedra[temp];
			}
			tetCirculator[0] = first;
			POSITION_T posFirst[3];
			//POSITION_T x[3];
			//POSITION_T y[3];
			//POSITION_T z[3];
			tetCirculator[0]->getCircumCenter( posFirst); //x[0]=posFirst[0]; y[0]=posFirst[1]; z[0]=posFirst[2];


			// get second
			int nt=0;
			for( ; nt<first->countNeighborTetrahedra &&
				!( first->neighborTetrahedra[nt]->contains( (VertexType*)vd->vertices[v]) &&  first->neighborTetrahedra[nt]->contains( (VertexType*)vd->vertices[v]->getNeighbor(nv))); nt++) ;
			tetCirculator[1] = first->neighborTetrahedra[nt];
			tetCirculatorLength = 2;

			POSITION_T pos[3];
			tetCirculator[1]->getCircumCenter( pos);

			// WRITE POINT
			fprintf( fp, "%f %f %f\n", pos[0], pos[1], pos[2]);
			pointCount++;

			// Circulate
			do{
				int nt=0;
				SimplexType *actual = tetCirculator[tetCirculatorLength-1];
				for( ; nt<actual->countNeighborTetrahedra &&
					!(  actual->neighborTetrahedra[nt] != tetCirculator[tetCirculatorLength-2] &&
						actual->neighborTetrahedra[nt]->contains( (VertexType*)vd->vertices[v]) &&
						actual->neighborTetrahedra[nt]->contains( (VertexType*)vd->vertices[v]->getNeighbor(nv))); nt++) ;
				tetCirculator[tetCirculatorLength] = actual->neighborTetrahedra[nt];
				tetCirculatorLength++;

				POSITION_T pos[3];
				tetCirculator[tetCirculatorLength-1]->getCircumCenter( pos);
				POSITION_T posLast[3];
				tetCirculator[tetCirculatorLength-2]->getCircumCenter( posLast);

				// WRITE POINT
				fprintf( fp, "%f %f %f\n", pos[0], pos[1], pos[2]);
				pointCount++;

				// WRITE EDGE

				/*if( first != tetCirculator[tetCirculatorLength-1]){
					x[1]=pos[0]; y[1]=pos[1]; z[1]=pos[2];
					x[2]=posLast[0]; y[2]=posLast[1]; z[2]=posLast[2];

					// WRITE FACE

				}*/
			}while( first != tetCirculator[tetCirculatorLength-1]);
			tetCirculatorLength--;
		}
	}

	fprintf( fp, "POLYGONS %i %i\n", polygonCount, polygonCount + pointCount);
	pointCount=0;

	// Faces
	fprintf( stderr, "Faces!\n");
	for( int v=0; v<vd->countVertices; v++)
	if( subDomain->contains( vd->vertices[v]))
	{
		//char color[100];
		//sprintf(color, "rgb <1,1,1>  transmit 0.5");

		SimplexType *tet = vd->getSimplexContainingPoint( vd->vertices[v]);

		for( int nv=0; nv<vd->vertices[v]->numberOfNeighbors(); nv++)
		//if(vd->vertices[v]->getNeighbor(nv)->index >= 0 /* no framepoints */)
		//if(vd->vertices[v]->index < vd->vertices[v]->getNeighbor(nv)->index)
		{
			//fprintf( stderr, "Construct Circulator for Point %i - Point %i\n", v, vd->vertices[v]->getNeighbor(nv)->index);

			char polygonString[2][1000] = {"\0","\0"};
			int countPointsPerPolygon=0;

			// Circulator of Tets around Point-NeighborPoint-Axis
			SimplexType *tetCirculator[100];
			int tetCirculatorLength = 0;

			// get (first) tet containing Points v & nv
			SimplexType *first = tet;
			while( !first->contains( (VertexType*)vd->vertices[v]->getNeighbor(nv))){
				// chose random neighbor tet containing point v
				int temp;
				do{
					// chose random neighbor tet
					temp = (int)(myRandE((double)first->countNeighborTetrahedra));
				}while( ! first->neighborTetrahedra[temp]->contains( vd->vertices[v]));
				first = first->neighborTetrahedra[temp];
			}
			tetCirculator[0] = first;
			POSITION_T posFirst[3];
			//POSITION_T x[3];
			//POSITION_T y[3];
			//POSITION_T z[3];
			tetCirculator[0]->getCircumCenter( posFirst); //x[0]=posFirst[0]; y[0]=posFirst[1]; z[0]=posFirst[2];


			// get second
			int nt=0;
			for( ; nt<first->countNeighborTetrahedra &&
				!( first->neighborTetrahedra[nt]->contains( (VertexType*)vd->vertices[v]) &&  first->neighborTetrahedra[nt]->contains( (VertexType*)vd->vertices[v]->getNeighbor(nv))); nt++) ;
			tetCirculator[1] = first->neighborTetrahedra[nt];
			tetCirculatorLength = 2;

			POSITION_T pos[3];
			tetCirculator[1]->getCircumCenter( pos);

			// WRITE POINT
			sprintf( polygonString[countPointsPerPolygon%2], "%s %i",
					polygonString[(countPointsPerPolygon+1)%2],
					pointCount);
			pointCount++;
			countPointsPerPolygon++;

			// Circulate
			do{
				int nt=0;
				SimplexType *actual = tetCirculator[tetCirculatorLength-1];
				for( ; nt<actual->countNeighborTetrahedra &&
					!(  actual->neighborTetrahedra[nt] != tetCirculator[tetCirculatorLength-2] &&
						actual->neighborTetrahedra[nt]->contains( (VertexType*)vd->vertices[v]) &&
						actual->neighborTetrahedra[nt]->contains( (VertexType*)vd->vertices[v]->getNeighbor(nv))); nt++) ;
				tetCirculator[tetCirculatorLength] = actual->neighborTetrahedra[nt];
				tetCirculatorLength++;

				POSITION_T pos[3];
				tetCirculator[tetCirculatorLength-1]->getCircumCenter( pos);
				POSITION_T posLast[3];
				tetCirculator[tetCirculatorLength-2]->getCircumCenter( posLast);

				// WRITE POINT
				sprintf( polygonString[countPointsPerPolygon%2], "%s %i",
						polygonString[(countPointsPerPolygon+1)%2],
						pointCount);
				pointCount++;
				countPointsPerPolygon++;

				// WRITE EDGE

				/*if( first != tetCirculator[tetCirculatorLength-1]){
					x[1]=pos[0]; y[1]=pos[1]; z[1]=pos[2];
					x[2]=posLast[0]; y[2]=posLast[1]; z[2]=posLast[2];

					// WRITE FACE

				}*/
			}while( first != tetCirculator[tetCirculatorLength-1]);
			tetCirculatorLength--;

			fprintf(fp, "%i%s \n", countPointsPerPolygon, polygonString[(countPointsPerPolygon+1)%2]);
		}
	}


	fclose(fp);
}
/*****************************************************************************/


template <int Dimensions, class VertexType, class SimplexType>
SimplexType* Triangulation<Dimensions,VertexType,SimplexType>::getSimplexContainingPoint( VertexType* point)
{
	Triangulation<Dimensions, VertexType, SimplexType> *voronoiDiagram = this;
	int i;

	// look for first tet containing new point
	//fprintf( stderr, "look for first tet containing new point\n");
	i=(int)(myRandE((double)voronoiDiagram->countTetrahedra));
	if(i>=voronoiDiagram->countTetrahedra){
		fprintf( stderr, "WRONG1! ;)\n");
		exit(0);
	}
	SimplexType *actualTet = NULL,
            *candidateTet = voronoiDiagram->tetrahedra[i];
			//*candidateTet = voronoiDiagram->tetrahedra[(int)(myRand()*(double)voronoiDiagram->countTetrahedra)];

	double distance;
	double candidateDistance = point->getDistanceTo( (Simplex<Dimensions>*)candidateTet);

	// chose closest random point as start point
	for( i=0; i<5; i++){
		int temp = (int)(myRandE((double)voronoiDiagram->countTetrahedra));
		if(temp>=voronoiDiagram->countTetrahedra){
			fprintf( stderr, "WRONG2! ;)\n");
			exit(0);
		}
		//actualTet = voronoiDiagram->tetrahedra[(int)(myRand()*(double)voronoiDiagram->countTetrahedra)];
		actualTet = voronoiDiagram->tetrahedra[temp];
		distance = point->getDistanceTo( (Simplex<Dimensions>*)actualTet);
		if( candidateDistance > distance){
			candidateDistance = distance;
			candidateTet = actualTet;
		}
	}


	// follow shortest way to first tetrahedron containing new point
	int stuck=0;
	while( candidateTet->vertices[0]!=point
			&& candidateTet->getVertex(1)!=point
			&& candidateTet->getVertex(2)!=point
			&& (Dimensions<3 || candidateTet->getVertex(3)!=point)
			){
		//if(point->getIndex() == 26259)
		//	fprintf( stderr, "INFO: candidateTet=%i\n", candidateTet->index);
		actualTet=candidateTet;
		//if(point->getIndex() == 26259)
		//	fprintf( stderr, "INFO: actualTet(%i), (%i %i %i %i), dist=%lf (%lf) to point %i ( %lf, %lf, %lf)\n", actualTet->index, actualTet->vertices[0]->getIndex(), actualTet->vertices[1]->getIndex(), actualTet->vertices[2]->getIndex(), actualTet->vertices[3]->getIndex(), /*minDistance*/0, point->getDistanceTo( (Simplex<Dimensions>*)actualTet), point->getIndex(), point->pos()[0], point->pos()[1], point->pos()[2]);


		for( i=0; i<actualTet->countNeighborTetrahedra; i++){
			distance = point->getDistanceTo( (Simplex<Dimensions>*)actualTet->neighborTetrahedra[i]);
			/*if(point->getIndex() == 26259){
				fprintf( stderr, "INFO: %i. neighbor tet(%i), (%i %i %i %i), dist=%lf\n", i+1, actualTet->neighborTetrahedra[i]->getIndex(), actualTet->neighborTetrahedra[i]->vertices[0]->getIndex(), actualTet->neighborTetrahedra[i]->vertices[1]->getIndex(), actualTet->neighborTetrahedra[i]->vertices[2]->getIndex(), actualTet->neighborTetrahedra[i]->vertices[3]->getIndex(), distance);
				fprintf( stderr, "INFO: dist=%lf\n", distance);
			}*/
			if( candidateDistance > distance){
			//if( minDistance > distance = getDistanceOfPointToTetrahedron( point, actualTet->neighborTetrahedra[i])){
				candidateTet      = actualTet->neighborTetrahedra[i];
				candidateDistance = distance;
			}
		}

		if( actualTet==candidateTet){
			//if(point->getIndex() == 26259)
			//	fprintf( stderr, "RANDOM CHOICE!!!\n");

			// TODO: IMPROVE!!! NOT AT ALL OPTIMAL SOLUTION! Mayhaps predicates would help

			if(stuck<1000){
				int temp = (int)(myRandE((double)actualTet->countNeighborTetrahedra));
				if(temp>=actualTet->countNeighborTetrahedra){
					fprintf( stderr, "WRONG3: %i / %i ;)\n", temp, actualTet->countNeighborTetrahedra);
					exit(0);
				}
				candidateTet      = actualTet->neighborTetrahedra[temp];
				stuck++;
			}
			else{
				//fprintf( stderr, "WARNING in getSimplexContainingPoint( Vertex %i): Algorithm stuck! Have to restart from new randomly chosen simplex!\n", point->getIndex());
				candidateTet      = this->tetrahedra[(int)(RAND01*(double)this->countTetrahedra)];
				stuck=0;
			}

			//candidateTet      = actualTet->neighborTetrahedra[(int)(myRand()*(double)actualTet->countNeighborTetrahedra)];
			candidateDistance = point->getDistanceTo( (Simplex<Dimensions>*)candidateTet);


		}

	}

	return candidateTet;
}
/*****************************************************************************/

template <int Dimensions>
double Vertex<Dimensions>::getDistanceTo( Simplex<Dimensions>* tet)
{
	Vertex<Dimensions> *point = this;
	int NR_Simplex_POINTS = Dimensions+1;
	int i=0, ii;

	double distance;
	double minDistance = 0;
	if( tet==NULL)
		fprintf(stderr, "Tet == NULL\n");
	if( point==NULL)
		fprintf(stderr, "Point == NULL\n");
	//fprintf(stderr, "Point %i (%.3lf, %.3lf, %.3lf)\n", point->index, point->position[0], point->position[1], point->position[2]);
	for( ii=0; ii<Dimensions; ii++)
		minDistance += pow( tet->getVertex(i)->x(ii) - point->x(ii), 2);
		//minDistance += myPow( tet->vertices[i]->position[ii] - point->position[ii], 2);

	for( i=1; i<NR_Simplex_POINTS; i++){
		distance = 0;
		for( ii=0; ii<Dimensions; ii++)
			distance += pow( tet->getVertex(i)->x(ii) - point->x(ii), 2);
			//distance += myPow( tet->vertices[i]->position[ii] - point->position[ii], 2);

		if( minDistance > distance)
			minDistance = distance;
	}

	// USE BARYCENTER
	/*int NR_Simplex_POINTS = Dimensions+1;
	double barycenter[Dimensions];
	for( int d=0; d<Dimensions; d++)
		barycenter[d]=tet->vertices[0]->position[d];

	for( int i=1; i<NR_Simplex_POINTS; i++)
		for( int d=0; d<Dimensions; d++)
			barycenter[d] += tet->vertices[i]->position[d];

	for( int d=0; d<Dimensions; d++)
		barycenter[d] /= (double)NR_Simplex_POINTS;


	double minDistance = 0;
	for( int d=0; d<Dimensions; d++)
		minDistance += pow( barycenter[d]-this->x(d), 2);
*/
	return minDistance;
}
/*****************************************************************************/

template <int Dimensions, class VertexType, class SimplexType>
POSITION_T Triangulation<Dimensions,VertexType,SimplexType>::getVoronoiCellVolume(VertexType *&v, POSITION_T *&interfaces) {

	switch(Dimensions){
	case 2: return getVoronoiCellVolume2D(v, interfaces);
	case 3: return getVoronoiCellVolume3D(v, interfaces);
	default: fprintf(stderr, "ERROR: getVoronoiCellVolume() not defined for dimension %i\n", Dimensions); exit(0);
	}
}
/*****************************************************************************/

template <int Dimensions, class VertexType, class SimplexType>
POSITION_T Triangulation<Dimensions,VertexType,SimplexType>::getVoronoiCellVolume(VertexType *&v) {

	POSITION_T _interfaces[v->numberOfNeighbors()];
	switch(Dimensions){
	case 2: return getVoronoiCellVolume2D(v, _interfaces);
	case 3: return getVoronoiCellVolume3D(v, _interfaces);
	default: fprintf(stderr, "ERROR: getVoronoiCellVolume() not defined for dimension %i\n", Dimensions); exit(0);
	}
}
/*****************************************************************************/

template <int Dimensions, class VertexType, class SimplexType>
POSITION_T Triangulation<Dimensions,VertexType,SimplexType>::getVoronoiCellVolume2D(VertexType *&v, POSITION_T *&interfaces) {

	assert(Dimensions==2);

	SimplexType *first = getSimplexContainingPoint( v);
	SimplexType *tet = first;
	SimplexType *last = first;

	POSITION_T array1[Dimensions];
	POSITION_T array2[Dimensions];
	POSITION_T *cc1 = array1, *cc2 = array2, *cctmp;
	tet->getCircumCenter( cc2);

	POSITION_T volume = 0.;

	do {
		//fprintf( stderr, "%i -> ", tet->index);

		int n=0;
		for( n=0; n<tet->countNeighborTetrahedra
			&& (!tet->neighborTetrahedra[n]->contains( v) || tet->neighborTetrahedra[n]==last); n++) ;

		last = tet;
		tet = tet->neighborTetrahedra[n];

		cctmp = cc1;
		cc1 = cc2;
		cc2 = cctmp;

		tet->getCircumCenter( cc2);

		// A = a*h_a
		/*POSITION_T a = sqrt( pow(cc1[0] - cc2[0], 2)
		               + pow(cc1[1] - cc2[1], 2));
		POSITION_T b = sqrt( pow(v->pos()[0] - cc2[0], 2)
		               + pow(v->pos()[1] - cc2[1], 2));
		POSITION_T c = sqrt( pow(cc1[0] - v->pos()[0], 2)
				       + pow(cc1[1] - v->pos()[1], 2));
		POSITION_T s = (a+b+c)/2.;
		volume += sqrt(s*(s-a)*(s-b)*(s-c));*/

		volume += triangleInterface<Dimensions,POSITION_T>( cc1, cc2, v->pos());

		// DETERMINE NEIGHBOR VERTEX CORRESPONDING TO EDGE
		//int j=0;
		int i=0;
		for( ; i<v->numberOfNeighbors() && (!tet->contains(v->getNeighbor(i)) || !last->contains(v->getNeighbor(i))); i++) ;
		interfaces[i] = sqrt( pow(cc1[0] - cc2[0], 2) + pow(cc1[1] - cc2[1], 2));
		//for( int i=0; i<=Dimensions && last->vertices[i]!=tet->vertices[j]; i++)
		//	for( j=0; j<=Dimensions && last->vertices[i]!=tet->vertices[j]; j++)

	} while (tet != first);

	//fprintf( stderr, "%lf -> ", volume);
	return volume;
}
/*****************************************************************************/



template <int Dimensions, class VertexType, class SimplexType>
POSITION_T Triangulation<Dimensions,VertexType,SimplexType>::getVoronoiCellVolume3D(VertexType *&vertex, POSITION_T *&interfaces) {

	assert(Dimensions==3);

	POSITION_T volume = 0.;

	/*for( int v=0; v<this->numberOfVertices(); v++)
	if( vertices[v] != NULL &&
	    vertices[v]->numberOfNeighbors()>0
	)*/
	int v=vertex->getIndex();
	{
		// get tet containing Point v
		SimplexType *tet = getSimplexContainingPoint( vertices[v]);
		//fprintf( stderr, "Point %i -> tet %i (%i, %i, %i, %i)\n", v, tet->index, tet->vertices[0]->index, tet->vertices[1]->index, tet->vertices[2]->index, tet->vertices[3]->index);

		for( int nv=0; nv<vertices[v]->numberOfNeighbors(); nv++){
			//fprintf( stderr, "Construct Circulator for Point %i - Point %i\n", v, vertices[v]->getNeighbor(nv)->getIndex());

			interfaces[nv] = 0.;
			POSITION_T height = 0.;

			// Circulator of Tets around Point-NeighborPoint-Axis
			SimplexType *tetCirculator[100];
			int tetCirculatorLength = 0;

			// get (first) tet containing Points v & nv
			SimplexType *first = tet;

			while( !first->contains( vertices[v]->getNeighbor(nv))){
				// chose random neighbor tet containing point v
				//if( v==26259)
				//	fprintf( stderr, "looking for tet containing vertex %i\n", vertices[v]->getNeighbor(nv)->getIndex());
				int temp;
				do{
					// chose random neighbor tet
					temp = (int)(myRandE((double)first->countNeighborTetrahedra));
				}while( ! first->neighborTetrahedra[temp]->contains( vertices[v]));
				first = first->neighborTetrahedra[temp];
				//if( v==26259)
				//	fprintf( stderr, "found tet containing vertex %i\n", vertices[v]->getNeighbor(nv)->getIndex());
			}
			tetCirculator[0] = first;
			POSITION_T posFirst[3];
			POSITION_T x[3];
			POSITION_T y[3];
			POSITION_T z[3];
			tetCirculator[0]->getCircumCenter( posFirst); x[0]=posFirst[0]; y[0]=posFirst[1]; z[0]=posFirst[2];


			// get second
			{
				int nt=0;
				//if( v==26259)
				//	fprintf( stderr, "looking for 2nd tet containing vertex %i\n", vertices[v]->getNeighbor(nv)->getIndex());
				for( ; nt<first->countNeighbors() &&
					!( first->neighborTetrahedra[nt]->contains( vertices[v]) &&  first->neighborTetrahedra[nt]->contains( vertices[v]->getNeighbor(nv))); nt++) ;
				//if( v==26259)
				//	fprintf( stderr, "found 2nd tet containing vertex %i\n", vertices[v]->getNeighbor(nv)->getIndex());
				tetCirculator[1] = first->neighborTetrahedra[nt];
				tetCirculatorLength = 2;

				POSITION_T pos[3];
				tetCirculator[1]->getCircumCenter( pos);
			}

			// Circulate
			POSITION_T pos[3];
			POSITION_T posLast[3];

			do{
				//if( v==26259)
				//	fprintf( stderr, "circulating\n");
				int nt=0;
				SimplexType *actual = tetCirculator[tetCirculatorLength-1];
				for( ; nt<actual->countNeighbors() &&
					!(  actual->neighborTetrahedra[nt] != tetCirculator[tetCirculatorLength-2] &&
						actual->neighborTetrahedra[nt]->contains(  (Vertex<Dimensions>*)vertices[v]) &&
						actual->neighborTetrahedra[nt]->contains( ((Vertex<Dimensions>*)vertices[v])->getNeighbor(nv))); nt++) ;
				tetCirculator[tetCirculatorLength] = actual->neighborTetrahedra[nt];
				tetCirculatorLength++;


				if( first != tetCirculator[tetCirculatorLength-1]){
					tetCirculator[tetCirculatorLength-1]->getCircumCenter( pos);
					tetCirculator[tetCirculatorLength-2]->getCircumCenter( posLast);
					x[1]=pos[0]; y[1]=pos[1]; z[1]=pos[2];
					x[2]=posLast[0]; y[2]=posLast[1]; z[2]=posLast[2];

					{
						// INTERFACE
						/*POSITION_T a=sqrt(pow(x[0]-x[1],2) + pow(y[0]-y[1],2) + pow(z[0]-z[1],2));
						POSITION_T b=sqrt(pow(x[0]-x[2],2) + pow(y[0]-y[2],2) + pow(z[0]-z[2],2));
						POSITION_T c=sqrt(pow(x[2]-x[1],2) + pow(y[2]-y[1],2) + pow(z[2]-z[1],2));
						POSITION_T s=(a+b+c)/2.;
						POSITION_T interface = sqrt(s*(s-a)*(s-b)*(s-c));*/

						POSITION_T interface = triangleInterface<Dimensions, POSITION_T>( posFirst, pos, posLast);
						if( !isnan(interface)){
							interfaces[nv] += interface;
							if( height==0 ){
								height = fabs(tetrahedronHeight<POSITION_T>( posFirst, pos, posLast, vertices[v]->pos()));
								if(isnan(height))
									height = 0;
							}
						}
						if(isnan(interfaces[nv])){
							//fprintf(stderr, "a=%lf, b=%lf, c=%lf, s=%lf\n",a,b,c,s);
							exit(0);
						}
					}
				}else{
					{
						// HEIGHT
						/*POSITION_T a[3] = {x[0]-x[1],y[0]-y[1],z[0]-z[1]};
						POSITION_T b[3] = {x[2]-x[1],y[2]-y[1],z[2]-z[1]};
						POSITION_T norm[3] ={
								a[1]*b[2] - a[2]*b[1],
								a[2]*b[0] - a[0]*b[2],
								a[0]*b[1] - a[1]*b[0]
						};
						POSITION_T norm_abs=sqrt(pow(norm[0],2)+pow(norm[1],2)+pow(norm[2],2));
						norm[0]/=norm_abs;
						norm[1]/=norm_abs;
						norm[2]/=norm_abs;
						POSITION_T height = norm[0]*(vertices[v]->x() - x[0])
									 + norm[1]*(vertices[v]->y() - y[0])
									 + norm[2]*(vertices[v]->z() - z[0]);
						//fprintf(stderr, "%lf %lf\n",height, area);
						volume += fabs(height)*interfaces[nv]/3.;*/

						/*tetCirculator[0]->getCircumCenter( posFirst);
						tetCirculator[1]->getCircumCenter( pos);
						tetCirculator[2]->getCircumCenter( posLast);*/
						// TODO:
						/*POSITION_T height = fabs(tetrahedronHeight<POSITION_T>( posFirst, pos, posLast, vertices[v]->pos()));
						for(int i=0; isnan(height) && i<tetCirculatorLength-3; i++)
						{
							tetCirculator[i  ]->getCircumCenter( posFirst);
							tetCirculator[i+1]->getCircumCenter( pos);
							tetCirculator[i+2]->getCircumCenter( posLast);
							height = fabs(tetrahedronHeight<POSITION_T>( posFirst, pos, posLast, vertices[v]->pos()));
							fprintf(stderr, "[%i]: (%.2lf,%.2lf,%.2lf)recalculate heigth = %.1lf => A=%.2lf =>\n",
									i, vertices[v]->x(),vertices[v]->y(),vertices[v]->z(),height, interfaces[nv]);
							//exit( 0);
						}*/
						volume += height * interfaces[nv]/3.;
						/*fprintf(stderr, "[h=%.1lf A=%.1lf]",
								tetrahedronHeight<POSITION_T>( posFirst, pos, posLast, vertices[v]->pos()),
								interfaces[nv]);*/

						//fprintf(stderr, "[h=%.2lf] ", height);
					}
				}
				/*if( tetCirculatorLength>100){
					fprintf( stderr, "ERROR in getVoronoiCellVolume3D(): tetCirculatorLength>1000\n");
					exit(0);
				}*/
			}while( first != tetCirculator[tetCirculatorLength-1]);
			tetCirculatorLength--;

			/*double x[tetCirculatorLength];
			double y[tetCirculatorLength];
			double z[tetCirculatorLength];
			for( int i=0; i< tetCirculatorLength; i++){
				fprintf( stderr, "-> tet %i (%i, %i, %i, %i)\n", i, tetCirculator[i]->vertices[0]->index, tetCirculator[i]->vertices[1]->index, tetCirculator[i]->vertices[2]->index, tetCirculator[i]->vertices[3]->index);
				double pos[DIMENSIONS];
				getCircumCenter( tetCirculator[i], pos);
				x[i] = pos[0];
				y[i] = pos[1];
				z[i] = pos[2];
			}
			PovrayIO::writePolygon( &fs, tetCirculatorLength, x, y, z, "rgb <1,0,1>");*/
			//fs.close();
			//exit(0);
		}
	}
	//fprintf(stderr, "\n");
	return volume;
}
/*****************************************************************************/



template <int Dimensions, class VertexType, class SimplexType>
void Triangulation<Dimensions,VertexType,SimplexType>::getVoronoiCellInterfaces(VertexType *&vertex, POSITION_T *&interfaces ) {

	switch(Dimensions){
	case 3: getVoronoiCellInterfaces3D(vertex, interfaces ); break;
	default: fprintf(stderr, "ERROR: getVoronoiCellInterfaces() not defined for dimension %i\n", Dimensions); exit(0);
	}
}
/*****************************************************************************/



template <int Dimensions, class VertexType, class SimplexType>
void Triangulation<Dimensions,VertexType,SimplexType>::getVoronoiCellInterfaces3D(VertexType *&vertex, POSITION_T *&_interfaces ) {

	assert(Dimensions==3);

	//POSITION_T volume = 0.;

	/*for( int v=0; v<this->numberOfVertices(); v++)
	if( vertices[v] != NULL &&
	    vertices[v]->numberOfNeighbors()>0
	)*/
	int v=vertex->getIndex();
	{
		// get tet containing Point v
		SimplexType *tet = getSimplexContainingPoint( vertices[v]);
		//fprintf( stderr, "Point %i -> tet %i (%i, %i, %i, %i)\n", v, tet->index, tet->vertices[0]->index, tet->vertices[1]->index, tet->vertices[2]->index, tet->vertices[3]->index);

		for( int nv=0; nv<vertices[v]->numberOfNeighbors(); nv++){
			//fprintf( stderr, "Construct Circulator for Point %i - Point %i\n", v, vertices[v]->getNeighbor(nv)->getIndex());

			_interfaces[nv] = 0.;

			// Circulator of Tets around Point-NeighborPoint-Axis
			SimplexType *tetCirculator[100];
			int tetCirculatorLength = 0;

			// get (first) tet containing Points v & nv
			SimplexType *first = tet;

			while( !first->contains( vertices[v]->getNeighbor(nv))){
				// chose random neighbor tet containing point v
				int temp;
				do{
					// chose random neighbor tet
					temp = (int)(myRandE((double)first->countNeighborTetrahedra));
				}while( ! first->neighborTetrahedra[temp]->contains( vertices[v]));
				first = first->neighborTetrahedra[temp];
			}
			tetCirculator[0] = first;
			POSITION_T posFirst[3];
			POSITION_T x[3];
			POSITION_T y[3];
			POSITION_T z[3];
			tetCirculator[0]->getCircumCenter( posFirst); x[0]=posFirst[0]; y[0]=posFirst[1]; z[0]=posFirst[2];


			// get second
			{
				int nt=0;
				for( ; nt<first->countNeighbors() &&
					!( first->neighborTetrahedra[nt]->contains( vertices[v]) &&  first->neighborTetrahedra[nt]->contains( vertices[v]->getNeighbor(nv))); nt++) ;
				tetCirculator[1] = first->neighborTetrahedra[nt];
				tetCirculatorLength = 2;

				POSITION_T pos[3];
				tetCirculator[1]->getCircumCenter( pos);
			}

			// Circulate
			do{
				int nt=0;
				SimplexType *actual = tetCirculator[tetCirculatorLength-1];
				for( ; nt<actual->countNeighbors() &&
					!(  actual->neighborTetrahedra[nt] != tetCirculator[tetCirculatorLength-2] &&
						actual->neighborTetrahedra[nt]->contains(  (Vertex<Dimensions>*)vertices[v]) &&
						actual->neighborTetrahedra[nt]->contains( ((Vertex<Dimensions>*)vertices[v])->getNeighbor(nv))); nt++) ;
				tetCirculator[tetCirculatorLength] = actual->neighborTetrahedra[nt];
				tetCirculatorLength++;

				POSITION_T pos[3];
				tetCirculator[tetCirculatorLength-1]->getCircumCenter( pos);
				POSITION_T posLast[3];
				tetCirculator[tetCirculatorLength-2]->getCircumCenter( posLast);

				if( first != tetCirculator[tetCirculatorLength-1]){
					x[1]=pos[0]; y[1]=pos[1]; z[1]=pos[2];
					x[2]=posLast[0]; y[2]=posLast[1]; z[2]=posLast[2];

					{
						// INTERFACE
						POSITION_T a=sqrt(pow(x[0]-x[1],2.) + pow(y[0]-y[1],2.) + pow(z[0]-z[1],2.));
						POSITION_T b=sqrt(pow(x[1]-x[2],2.) + pow(y[1]-y[2],2.) + pow(z[1]-z[2],2.));
						POSITION_T c=sqrt(pow(x[2]-x[0],2.) + pow(y[2]-y[0],2.) + pow(z[2]-z[0],2.));
						POSITION_T s=(a+b+c)/2.;
						POSITION_T interface = sqrt(s*(s-a)*(s-b)*(s-c));
						if( !isnan(interface))
							_interfaces[nv] += interface;
						if(isnan(_interfaces[nv])){
							fprintf(stderr, "a=%lf, b=%lf, c=%lf, s=%lf\n",a,b,c,s);
							exit(0);
						}

					}
				}/*else{
					{
						// HEIGHT
						POSITION_T a[3] = {x[0]-x[1],y[0]-y[1],z[0]-z[1]};
						POSITION_T b[3] = {x[2]-x[1],y[2]-y[1],z[2]-z[1]};
						POSITION_T norm[3] ={
								a[1]*b[2] - a[2]*b[1],
								a[2]*b[0] - a[0]*b[2],
								a[0]*b[1] - a[1]*b[0]
						};
						POSITION_T norm_abs=sqrt(pow(norm[0],2)+pow(norm[1],2)+pow(norm[2],2));
						norm[0]/=norm_abs;
						norm[1]/=norm_abs;
						norm[2]/=norm_abs;
						POSITION_T height = norm[0]*(vertices[v]->x() - x[0])
									 + norm[1]*(vertices[v]->y() - y[0])
									 + norm[2]*(vertices[v]->z() - z[0]);
						//fprintf(stderr, "%lf %lf\n",height, area);
						volume += fabs(height)*area/3.;

					}
				}*/
			}while( first != tetCirculator[tetCirculatorLength-1]);
			tetCirculatorLength--;

			if(isnan(_interfaces[nv])){
				fprintf(stderr, "interfaces[%i]=%lf\n",nv,_interfaces[nv]);
				exit(0);
			}


			/*double x[tetCirculatorLength];
			double y[tetCirculatorLength];
			double z[tetCirculatorLength];
			for( int i=0; i< tetCirculatorLength; i++){
				fprintf( stderr, "-> tet %i (%i, %i, %i, %i)\n", i, tetCirculator[i]->vertices[0]->index, tetCirculator[i]->vertices[1]->index, tetCirculator[i]->vertices[2]->index, tetCirculator[i]->vertices[3]->index);
				double pos[DIMENSIONS];
				getCircumCenter( tetCirculator[i], pos);
				x[i] = pos[0];
				y[i] = pos[1];
				z[i] = pos[2];
			}
			PovrayIO::writePolygon( &fs, tetCirculatorLength, x, y, z, "rgb <1,0,1>");*/
			//fs.close();
			//exit(0);
		}
	}

	//return volume;
}
/*****************************************************************************/



template <int Dimensions, class VertexType, class SimplexType>
inline void Triangulation<Dimensions, VertexType, SimplexType>::printVoronoiDiagramToEps( const char * filename)
{
	switch( Dimensions){
	case 2:
		printVoronoiDiagramToEps2D( filename);
		break;
	}
}

template <int Dimensions, class VertexType, class SimplexType>
inline void Triangulation<Dimensions, VertexType, SimplexType>::printVoronoiDiagramToEps2D( const char * filename)
{
	float scale = 1.;
	Triangulation<Dimensions, VertexType, SimplexType> *voronoiDiagram = this;

	std::fstream fs;
	fs.open(filename, std::fstream::out);
	EPS::PSwriteHeader(&fs, domain->min(0) * scale,
			domain->max(0) * scale, domain->min(1) * scale,
			domain->max(1) * scale);
	for (int t = 0; t < voronoiDiagram->countTetrahedra; t++) {
		POSITION_T vp1[2], vp2[2];
		voronoiDiagram->tetrahedra[t]->getCircumCenter( vp1);
		for (int n = 0; n
				< voronoiDiagram->tetrahedra[t]->countNeighborTetrahedra; n++) {
			voronoiDiagram->tetrahedra[t]->neighborTetrahedra[n]->getCircumCenter( vp2);
			EPS::PSwriteLine(&fs, vp1[0] * scale,
					vp1[1] * scale,
					//voronoiDiagram->vertices[c]->neighborCells[n]->position[0]*scale,
					//voronoiDiagram->vertices[c]->neighborCells[n]->position[1]*scale,
					(vp1[0] + vp2[0]) / 2. * scale, (vp1[1] + vp2[1]) / 2. * scale,
					0.1, (char*) (" 0 0 0 setrgbcolor"));
			/*if( voronoiDiagram->vertices[c]==centralCell)
			 EpsIO::PSwriteCircle( &fs,
			 voronoiDiagram->vertices[c]->position[0]*scale,
			 voronoiDiagram->vertices[c]->position[1]*scale,
			 1,
			 1,
			 " 0 0 1 setrgbcolor");
			 if( c==voronoiDiagram->countVertices-1)
			 EpsIO::PSwriteCircle( &fs,
			 voronoiDiagram->vertices[c]->position[0]*scale,
			 voronoiDiagram->vertices[c]->position[1]*scale,
			 1,
			 1,
			 " 1 0 0 setrgbcolor");*/
		}
		// /Times-Roman 270 selectfont
		/*if( plotIndex && voronoiDiagram->vertices[c]->countNeighbors>0){
		 fs << "/Times-Roman 4 selectfont" << endl;
		 fs << "0.5 setgray " << voronoiDiagram->vertices[c]->position[0]*scale << " " << voronoiDiagram->vertices[c]->position[1]*scale << " moveto ("<< voronoiDiagram->vertices[c]->index <<") show" << endl;
		 }*/
	}

	if(false)
	for (int c = 0; c < voronoiDiagram->countVertices; c++) {
		if (voronoiDiagram->vertices[c]->numberOfNeighbors() > 0) {
			char color[200];
			float radius;

			/*if (voronoiDiagram->vertices[c]->refined) {
				radius = 0.5;
			} else {
				radius = 2.;
			}

			switch (voronoiDiagram->vertices[c]->getState()) {
			case FREE:
				sprintf(color, " 0 0 0 setrgbcolor");
				break;
			case ACTIVE:
				sprintf(color, " 1 0 0 setrgbcolor");
				break;
			case NONACTIVE:
				sprintf(color, " 0 1 0 setrgbcolor");
				break;
			case COMPARTMENT:
				if (GetAgent(voronoiDiagram->vertices[c])->countFree
						== GetAgent(voronoiDiagram->vertices[c])->maxCellCount)
					sprintf(color, " 0 0 0 setrgbcolor");
				else if (GetAgent(voronoiDiagram->vertices[c])->countActive
						== GetAgent(voronoiDiagram->vertices[c])->maxCellCount)
					sprintf(color, " 1 0 0 setrgbcolor");
				else if (GetAgent(voronoiDiagram->vertices[c])->countNonactive
						== GetAgent(voronoiDiagram->vertices[c])->maxCellCount)
					sprintf(color, " 0 1 0 setrgbcolor");
				else
					sprintf(color, " 1 0 1 setrgbcolor");
				radius = 2.;
				break;
			}*/

			radius = 2.;
			sprintf(color, " 1 0 0 setrgbcolor");

			//if(voronoiDiagram->vertices[c]->refined)
			//if (voronoiDiagram->vertices[c]->getState() != FREE)
			EPS::PSwriteCircle(&fs,
						voronoiDiagram->vertices[c]->pos()[0] * scale,
						voronoiDiagram->vertices[c]->pos()[1] * scale,
						radius, 0., color);//" 1 0 0 setrgbcolor");
			/*else
			 EpsIO::PSwriteCircle( &fs,
			 voronoiDiagram->vertices[c]->position[0]*scale,
			 voronoiDiagram->vertices[c]->position[1]*scale,
			 2,
			 2,
			 " 0 0 1 setrgbcolor");
			 */
		}
	}

	fs.close();

}
/*****************************************************************************/


template <int Dimensions, class VertexType, class SimplexType>
inline void Triangulation<Dimensions, VertexType, SimplexType>::printVoronoiDiagramToPovray2D( const char * filename, Domain<Dimensions>* subDomain)
{
	Triangulation<Dimensions, VertexType, SimplexType> *voronoiDiagram = this;

	std::fstream fs;
	fs.open(filename, std::fstream::out);
	PovrayIO::writePovrayHeader(&fs,
			domain->min(0) * 10,
			domain->min(1) * 10,
			domain->max(0) * 10,
			domain->max(1) * 10);
	for (int t = 0; t < voronoiDiagram->countTetrahedra; t++) {
		POSITION_T vp1[2], vp2[2];
		voronoiDiagram->tetrahedra[t]->getCircumCenter( vp1);
		for (int n = 0; n < voronoiDiagram->tetrahedra[t]->countNeighborTetrahedra; n++) {
			voronoiDiagram->tetrahedra[t]->neighborTetrahedra[n]->getCircumCenter( vp2);
			if( pow(vp1[0]-vp2[0],2) + pow(vp1[1]-vp2[1],2) > 0.0001)
			PovrayIO::writeCylinder(&fs, vp1[0] * 10, vp1[1] * 10,
					//voronoiDiagram->vertices[c]->neighborCells[n]->position[0]*10,
					//voronoiDiagram->vertices[c]->neighborCells[n]->position[1]*10,
					(vp1[0] + vp2[0]) / 2. * 10, (vp1[1] + vp2[1]) / 2. * 10,
					0.1, (char*) (" Black"));
			/*if( voronoiDiagram->vertices[c]==centralCell)
			 EpsIO::PSwriteCircle( &fs,
			 voronoiDiagram->vertices[c]->position[0]*10,
			 voronoiDiagram->vertices[c]->position[1]*10,
			 1,
			 1,
			 " 0 0 1 setrgbcolor");
			 if( c==voronoiDiagram->countVertices-1)
			 EpsIO::PSwriteCircle( &fs,
			 voronoiDiagram->vertices[c]->position[0]*10,
			 voronoiDiagram->vertices[c]->position[1]*10,
			 1,
			 1,
			 " 1 0 0 setrgbcolor");*/
		}
		// /Times-Roman 270 selectfont
		/*if( plotIndex && voronoiDiagram->vertices[c]->countNeighbors>0){
		 fs << "/Times-Roman 4 selectfont" << endl;
		 fs << "0.5 setgray " << voronoiDiagram->vertices[c]->position[0]*10 << " " << voronoiDiagram->vertices[c]->position[1]*10 << " moveto ("<< voronoiDiagram->vertices[c]->index <<") show" << endl;
		 }*/
	}

	for (int c = 0; c < voronoiDiagram->countVertices; c++) {
		if (voronoiDiagram->vertices[c]->numberOfNeighbors() > 0) {
			char color[200];
			float radius;

			/*if (voronoiDiagram->vertices[c]->refined) {
				radius = 0.5;
			} else {
				radius = 2.;
			}

			switch (voronoiDiagram->vertices[c]->getState()) {
			case FREE:
				sprintf(color, " 0 0 0 setrgbcolor");
				break;
			case ACTIVE:
				sprintf(color, " 1 0 0 setrgbcolor");
				break;
			case NONACTIVE:
				sprintf(color, " 0 1 0 setrgbcolor");
				break;
			case COMPARTMENT:
				if (GetAgent(voronoiDiagram->vertices[c])->countFree
						== GetAgent(voronoiDiagram->vertices[c])->maxCellCount)
					sprintf(color, " 0 0 0 setrgbcolor");
				else if (GetAgent(voronoiDiagram->vertices[c])->countActive
						== GetAgent(voronoiDiagram->vertices[c])->maxCellCount)
					sprintf(color, " 1 0 0 setrgbcolor");
				else if (GetAgent(voronoiDiagram->vertices[c])->countNonactive
						== GetAgent(voronoiDiagram->vertices[c])->maxCellCount)
					sprintf(color, " 0 1 0 setrgbcolor");
				else
					sprintf(color, " 1 0 1 setrgbcolor");
				radius = 2.;
				break;
			}*/

			radius = 2.;
			sprintf(color, " Red");

			//if(voronoiDiagram->vertices[c]->refined)
			//if (voronoiDiagram->vertices[c]->getState() != FREE)
				/*PovrayIO::writeCircle(&fs,
						voronoiDiagram->vertices[c]->position[0] * 10,
						voronoiDiagram->vertices[c]->position[1] * 10,
						radius, radius/4., color);//" 1 0 0 setrgbcolor");*/

				PovrayIO::writeSphere(&fs,
						voronoiDiagram->vertices[c]->pos()[0] * 10,
						voronoiDiagram->vertices[c]->pos()[1] * 10,
						0.,
						radius, color);//" 1 0 0 setrgbcolor");

			/*else
			 EpsIO::PSwriteCircle( &fs,
			 voronoiDiagram->vertices[c]->position[0]*10,
			 voronoiDiagram->vertices[c]->position[1]*10,
			 2,
			 2,
			 " 0 0 1 setrgbcolor");
			 */
		}
	}

	fs.close();

}

template <int Dimensions>
inline void Vertex<Dimensions>::setx( POSITION_T x, int dimension)
{	position[dimension] = x;}

template <int Dimensions>
inline void Vertex<Dimensions>::setx( POSITION_T x)
{	position[0] = x;}

template <int Dimensions>
inline void Vertex<Dimensions>::sety( POSITION_T y)
{	position[1] = y;}

template <int Dimensions>
inline void Vertex<Dimensions>::setz( POSITION_T z)
{	position[2] = z;}

/*template <int Dimensions>
void Vertex<Dimensions>::setIndex( int &index)
{	this->index = index;}
*/
template <int Dimensions>
inline void Vertex<Dimensions>::setIndex( int index)
{	this->index = index;}

template <int Dimensions>
inline int Vertex<Dimensions>::getIndex()
{	return this->index;}

template <int Dimensions>
inline POSITION_T Vertex<Dimensions>::x( int dimension)
{	return position[dimension];}

template <int Dimensions>
inline POSITION_T* Vertex<Dimensions>::pos()
{	return position;}

template <int Dimensions>
inline POSITION_T Vertex<Dimensions>::x()
{	return position[0];}

template <int Dimensions>
inline POSITION_T Vertex<Dimensions>::y()
{	return position[1];}

template <int Dimensions>
inline POSITION_T Vertex<Dimensions>::z()
{	return position[2];}

template <int Dimensions>
Domain<Dimensions>::Domain( POSITION_T min[Dimensions], POSITION_T max[Dimensions])
{
	for( int d=0; d<Dimensions; d++){
		xMin[d] = min[d];
		xMax[d] = max[d];
	}
}

template <int Dimensions>
Domain<Dimensions>::Domain( POSITION_T xmin, POSITION_T xmax, POSITION_T ymin, POSITION_T ymax, POSITION_T zmin, POSITION_T zmax)
{
	xMin[0] = xmin;
	xMax[0] = xmax;
	if( Dimensions>=2){
		xMin[1] = ymin;
		xMax[1] = ymax;
	}
	if( Dimensions>=3){
		xMin[2] = zmin;
		xMax[2] = zmax;
	}
}

template <int Dimensions>
POSITION_T Domain<Dimensions>::min( int dimension)
{	return xMin[dimension];}

template <int Dimensions>
POSITION_T Domain<Dimensions>::max( int dimension)
{	return xMax[dimension];}

template <int Dimensions>
bool Domain<Dimensions>::contains( Vertex<Dimensions> *v)
{
	for( int d=0; d< Dimensions; d++)
		if( v->x(d) < min(d) || v->x(d) > max(d))
			return false;
	return true;
}

template <int Dimensions, class VertexType>
double **Simplex<Dimensions,VertexType>::_A;
template <int Dimensions, class VertexType>
double **Simplex<Dimensions,VertexType>::_B;
template <int Dimensions, class VertexType>
double *Simplex<Dimensions,VertexType>::_x;
template <int Dimensions, class VertexType>
int Simplex<Dimensions,VertexType>::countUsers=0;


#if false// COMMENT_OLD_CODE
#ifdef USE_AGENT
void Vertex<Dimensions>::actualizeFreeNeighbors()
{
	this->countFreeNeighborCells = 0;
	for( int n=0; n<this->countNeighbors; n++)
		if( this->neighbors[n]->isFree())
			this->countFreeNeighborCells++;
}


void Vertex<Dimensions>::checkFreeNeighbors()
{
	int countFree = 0;
	for( int n=0; n<this->countNeighbors; n++)
		if( this->neighbors[n]->isFree())
			countFree++;

	if( this->countFreeNeighborCells != countFree)
	{
		fprintf( stderr, "Vertex %i has a wrong number of free neighbors: %i (!= %i found)\n", this->index, this->countFreeNeighborCells, countFree);
		exit( 0);
	}
}


void Triangulation<Dimensions>::refine( Vertex<Dimensions> *vc, int scale, ActionTree *actionTree)
{
	int base_x = vc->position[0];
	int base_y = vc->position[1];
#if Dimensions >= 3
	int base_z = vc->position[2];
#endif

	//fprintf( stderr, "REFINE( cell: %i, scale: %i)\n", vc->index, scale);

	// LIST OF NEW VERTICES
	//int countNewVertices = scale*scale*scale;
	//fprintf( stderr, "countNewVertices = %i -> %i\n", countNewVertices, (this->countVertices+countNewVertices));

	// LIST OF NEIGHBORS
	Vertex<Dimensions> *neighbors[vc->countNeighbors];
	//int countNeighbors = vc->countNeighbors;
	for( int n=0; n<vc->countNeighbors; n++)
		neighbors[n] = vc->neighbors[n];

	// INIT STACK
	int actualStackElement = 0;
	int stackSize = 0;
	Vertex<Dimensions> *stack[2000];
	for( stackSize=0; stackSize<vc->countNeighbors; stackSize++){
		// add neighbors of refined cells
		stack[stackSize] = vc->neighbors[stackSize];
		//stack[stackSize]->countFreeNeighborCells = -1;
		//stackSize++;
		//fprintf( stderr, "-> add point %ith point to stack (neighbor of original cells)\n", stackSize+1);
	}

	// ADD NEW VERTICES & POINTS
	int i=this->countVertices;//-countNewVertices;
	int firstIndexNew = this->countVertices;
	Vertex<Dimensions> *newVC;
	for( int x=0; x<scale; x++)
		for( int y=0; y<scale; y++)
#if Dimensions >= 3
			for( int z=0; z<scale; z++)
#endif
			{

			// allocate Voronoi Cell
			double rx = (float)(base_x + (x+myRandE(1.))/(double)scale);
			double ry = (float)(base_y + (y+myRandE(1.))/(double)scale);
#if Dimensions >= 3
			double rz = (float)(base_z + (z+myRandE(1.))/(double)scale);
#endif
			//this->vertices[i] = newVoronoiCell( rx, ry, rz);
#if Dimensions >= 3
			newVC = newVoronoiCell( rx, ry, rz);
#elif Dimensions == 2
			newVC = newVoronoiCell( rx, ry, 0);
#endif

			newVC->index = i;
			newVC->position[0] = rx;
			newVC->position[1] = ry;
#if Dimensions >= 3
			newVC->position[2] = rz;
#endif
			newVC->glucose = vc->glucose;
			newVC->oxygen  = vc->oxygen;

			newVC->refined = true;

			newVC->coarseParent = vc;

			insertVoronoiCell( this, newVC);
			//fprintf( stderr, "-> add point %ith point to stack (refined cell)\n", stackSize+1);
			//checkDelaunayCondition( this, 0, 0);

			//fprintf( stderr, "-> add point %i (%lf, %lf, %lf)\n", i, rx, ry, (Dimensions>=3?rz:0));
			//fprintf( stderr, "-> add point %i (%lf, %lf)\n", i, rx, ry);

			// add refined cells
			stack[stackSize] = newVC;
			//stack[stackSize]->countFreeNeighborCells = -1;
			stackSize++;

			i++;
		}

	// DELETE OLD POINT FROM TRIANGULATION
	//fprintf( stderr, "remove point %i\n", vc->index);
	removeVoronoiCell( this, vc);
	vc->refined = true;
	vc->refinement = firstIndexNew;

	//vc = newVC;

	// ADD NEIGHBORS OF REFINED CELLS
	for( int c=firstIndexNew; c<this->countVertices; c++){
		// add neighbors
		for( int n=0; n<this->vertices[c]->countNeighbors; n++){
			if( !isElementOf( stack, stackSize, this->vertices[c]->neighbors[n])){
				stack[stackSize] = this->vertices[c]->neighbors[n];
				//stack[stackSize]->countFreeNeighborCells = -1;
				stackSize++;
				//fprintf( stderr, "-> add point %ith point to stack (neighbor of refined cell)\n", stackSize);
			}

			//fprintf( stderr, "STACK to actualize free neighbor cells in actualized lattice is not implemented yet!\n");
			//exit( 0);
		}
	}

	// ACTUALIZE ALL CONSIDERED VERTICES
	for( actualStackElement=0; actualStackElement<stackSize; actualStackElement++ ){
		// actualize number of free neighbor Voronoi cells
		stack[actualStackElement]->actualizeFreeNeighbors();

		// actualize attached agent
		if(	stack[actualStackElement]->agent!=NULL){
			//fprintf( stderr, "Actualize Agent on Voronoi cell %i\n", stack[actualStackElement]->index);
			GetAgent(stack[actualStackElement])->actualize( actionTree);
		}
	}
	//checkDelaunayCondition( this, 0,0);
	//fprintf( stderr, "...finished\n");
}

Vertex<Dimensions> *Triangulation<Dimensions>::coarsen( Vertex<Dimensions> *vc, int scale, ActionTree *actionTree, AgentList *agentList)
{
	double base_x = floor(vc->position[0]);
	double base_y = floor(vc->position[1]);
#if Dimensions == 3
	double base_z = floor(vc->position[2]);
#endif

#if Dimensions == 3
	fprintf( stderr, "COARSEN( cell: %i (%i), scale: %i)\n", vc->index, (int)base_x*this->xN[1]*this->xN[2] + (int)base_y*this->xN[2] + (int)base_z, scale);
#else
	//fprintf( stderr, "COARSEN( cell: %i (%i), scale: %i)\n", vc->index, (int)base_x*this->xN[1] + (int)base_y, scale);
#endif
	assert( vc->refined);
	assert( vc->countNeighbors>0);


	// SEARCH & DELETE CONSIDERED POINTS
	int actualStackElement = 0;
	int stackSize = 0;
	Vertex<Dimensions> *stack[2000];
	stack[stackSize] = vc; stackSize++;

	int stack2Size = 0;
	Vertex<Dimensions> *stack2[2000];

	/* do{ // OLD
		for( int n=0; n<stack[actualStackElement]->countNeighbors; n++){
			if( stack[actualStackElement]->neighbors[n]->position[0] > base_x &&
				stack[actualStackElement]->neighbors[n]->position[0] < base_x+1 &&
				stack[actualStackElement]->neighbors[n]->position[1] > base_y &&
				stack[actualStackElement]->neighbors[n]->position[1] < base_y+1 &&
#if Dimensions == 3
				stack[actualStackElement]->neighbors[n]->position[2] > base_z &&
				stack[actualStackElement]->neighbors[n]->position[2] < base_z+1 &&
#endif
				!isElementOf( stack, stackSize, stack[actualStackElement]->neighbors[n])){

				stack[stackSize] = stack[actualStackElement]->neighbors[n];
				stackSize++;
			}else{
				// mark neighbors to be actualized
				//stack[actualStackElement]->neighbors[n]->countFreeNeighborCells = -1;
				stack2[stack2Size] = stack[actualStackElement]->neighbors[n];
				stack2Size++;
			}
		}

//		removeVoronoiCell( this, stack[actualStackElement]);

		actualStackElement++;
	}while( actualStackElement<stackSize);*/
	//fprintf( stderr, "found %i cells to fuse\n", stackSize);

	stackSize=0;
	int i_base = vc->coarseParent->refinement;
	for( int i=0; i<scale*scale
#if Dimensions == 3
								*scale
#endif
								; i++){
		//fprintf( stderr, "->coarsen %i\n", i_base+i);
		stack[stackSize] = this->vertices[i_base+i];
		stackSize++;
	}
	for( int i=0; i<stackSize; i++){
		for( int n=0; n<stack[i]->countNeighbors; n++)
			if( !isElementOf( stack, stackSize, stack[i]->neighbors[n]) &&
				!isElementOf( stack2, stack2Size, stack[i]->neighbors[n])	){

				stack2[stack2Size] = stack[i]->neighbors[n];
				stack2Size++;

			}
	}


	// REMOVE FINE POINTS FROM TRIANGULATION
	//fprintf( stderr, "removeVoronoiCell\n");
	for( int i=0; i<stackSize; i++){
		//fprintf( stderr, "removeVoronoiCell: %i \b", stack[i]->index);
		removeVoronoiCell( this, stack[i]);
		//fprintf( stderr, "...finished\n");

		// Deactivate Agent
		agentList->deactivateAgent( GetAgent(stack[i]));

		// Delete Actions from List
		actionTree->deleteAllActions( GetAgent(stack[i]));

		// Detach from Grid
		GetAgent(stack[i])->detach();

		// Free Memory of coarsened Lattice Point
		this->vertices[stack[i]->index]=NULL;
		free( stack[i]);
		stack[i]=NULL;
	}


	// REINSERT COARSE POINT INTO TRIANGULATION
	//fprintf( stderr, "triangulateVoronoiCell\n");
#if Dimensions == 3
	int index = base_x*this->xN[1]*this->xN[2] + base_y*this->xN[2] + base_z;
#elif Dimensions == 2
	int index = base_x*this->xN[1] + base_y;
#elif Dimensions == 1
	int index = base_x;
#endif
	triangulateVoronoiCell( this, this->vertices[index]);
	this->vertices[index]->refined = false;


	// ACTUALIZE ALL CONSIDERED VERTICES

	// mark direct
	//fprintf( stderr, "mark direct\n");
	for( int n=0; n<this->vertices[index]->countNeighbors; n++ )
		if( this->vertices[index]->neighbors[n]->countFreeNeighborCells != -1){
			this->vertices[index]->neighbors[n]->countFreeNeighborCells = -1;
			stack2[stack2Size] = this->vertices[index]->neighbors[n];
			stack2Size++;
		}

	// mark indirect
	//fprintf( stderr, "mark indirect\n");
	for( int i=0; i<stack2Size; i++ )
		for( int n=0; n<stack2[i]->countNeighbors; n++ ){
		//if( stack2[n]->countFreeNeighborCells != -1){
			stack2[i]->neighbors[n]->countFreeNeighborCells = -1;
		}

	// actualize all
	//fprintf( stderr, "actualize all\n");
	int ii=0;
	for( int i=0; i<stack2Size; i++ )
		for( int n=0; n<stack2[i]->countNeighbors; n++ ){
			if( stack2[i]->neighbors[n]->countFreeNeighborCells == -1){

				// actualize number of free neighbor Voronoi cells
				stack2[i]->neighbors[n]->actualizeFreeNeighbors();
				//stack[ii] = stack2[i]->neighbors[n];
				//ii++;

				// actualize attached agent
				if(	stack2[i]->neighbors[n]->agent!=NULL)
					GetAgent(stack2[i]->neighbors[n])->actualize( actionTree);
			}
		}

//	for( int i=0; i<ii; i++ )
//		if(	stack[i]->agent!=NULL)
//			GetAgent(stack[i])->actualize( actionTree);
	fprintf( stderr, "...finished\n");

	return this->vertices[index];
}


bool Triangulation<Dimensions>::isHomogen( Vertex<Dimensions> *vc, int scale, int &type)
{
	double base_x = floor(vc->position[0]);
	double base_y = floor(vc->position[1]);
#if Dimensions == 3
	double base_z = floor(vc->position[2]);
#endif

	//fprintf( stderr, "HOMOGEN? cell: %i (%i)\n", vc->index, (int)base_x*this->xN[1] + (int)base_y);
	//assert( vc->refined);
	//assert( vc->countNeighbors>0);

	if( GetAgent(vc->coarseParent)->countFree == GetAgent(vc->coarseParent)->maxCellCount){
		type = FREE;
		//fprintf( stderr, "FREE\n"); //exit(0);
		return true;
	}else
	if(	GetAgent(vc->coarseParent)->countActive == GetAgent(vc->coarseParent)->maxCellCount){
		type = ACTIVE;
		//fprintf( stderr, "ACTIVE\n"); //exit(0);
		return true;
	}else
	if(	GetAgent(vc->coarseParent)->countNonactive == GetAgent(vc->coarseParent)->maxCellCount ){
		type = NONACTIVE;
		//fprintf( stderr, "NONACTIVE\n"); exit(0);
		return true;
	}else{
		//fprintf( stderr, "NOT HOMOGENOUS!\n");
		return false;
	}
	// TEST HOMOGENITY over CONSIDERED POINTS
	/*int actualStackElement = 0;
	int stackSize = 0;
	Vertex<Dimensions> *stack[200];
	stack[stackSize] = vc; stackSize++;

	type = vc->getState();

	do{
		for( int n=0; n<stack[actualStackElement]->countNeighbors; n++){
			if( stack[actualStackElement]->neighbors[n]->position[0] > base_x &&
				stack[actualStackElement]->neighbors[n]->position[0] < base_x+1 &&
				stack[actualStackElement]->neighbors[n]->position[1] > base_y &&
				stack[actualStackElement]->neighbors[n]->position[1] < base_y+1 &&
#if Dimensions == 3
				stack[actualStackElement]->neighbors[n]->position[2] > base_z &&
				stack[actualStackElement]->neighbors[n]->position[2] < base_z+1 &&
#endif
				!isElementOf( stack, stackSize, stack[actualStackElement]->neighbors[n])){
				//fprintf( stderr, "> neighbor: %i\n", stack[actualStackElement]->neighbors[n]->index);
				stack[stackSize] = stack[actualStackElement]->neighbors[n];
				stackSize++;

				if( stack[actualStackElement]->neighbors[n]->getState() != type){
					type =-1;
					return false;
				}
			}
		}


//		removeVoronoiCell( this, stack[actualStackElement]);

		actualStackElement++;
	}while( actualStackElement<stackSize);*/

#if Dimensions == 3
	int index = base_x*this->xN[1]*this->xN[2] + base_y*this->xN[2] + base_z;
#elif Dimensions == 2
	int index = base_x*this->xN[1] + base_y;
#elif Dimensions == 1
	int index = base_x;
#endif

	Vertex<Dimensions> **refinement = &this->vertices[this->vertices[index]->refinement];
	/*fprintf( stderr, "fine point:   (%lf, %lf %lf)\n", vc->position[0], vc->position[1], vc->position[2]);
	fprintf( stderr, "coarse point: (%lf, %lf %lf)\n", this->vertices[index]->position[0], this->vertices[index]->position[1], this->vertices[index]->position[2]);
	if( this->vertices[index]->refinement==NULL)
		exit( 0);
	if( this->vertices[index]->refinement[0]==NULL)
		exit( 0);
	fprintf( stderr, "1st refined point: (%lf, %lf %lf)\n", this->vertices[index]->refinement[0]->position[0], this->vertices[index]->refinement[0]->position[1], this->vertices[index]->refinement[0]->position[2]);
*/
	/*type = refinement[0]->getState();
	for( int i=1; i<scale*scale
#if Dimensions == 3
								*scale
#endif
								; i++){
		if(	refinement[i]->getState() != type)
			return false;
	}*/
	//fprintf( stderr, "homogeniously %i (%i cells)\n", type, actualStackElement);

	return true;
}

#endif
#endif // COMMENT_OLD_CODE
