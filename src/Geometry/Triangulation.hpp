#ifndef __TRIANGULATION_H
#define __TRIANGULATION_H


//#include "CellProcesses.h"

// DEFINES
#define FILENAMESIZE 512
#define READBUFFERSIZE 20000

#define MAX_TRIANGLES 6727210
#define MAX_NNS 30

#define INIT_FACES_ARRAY_SIZE 10

#define FACES_ARRAY_EXTENSION_SIZE 10
#define TETRAHEDRA_ARRAY_EXTENSION_SIZE 10
#define POINTS_ARRAY_EXTENSION_SIZE 100
/*
//#define Dimensions 3
#define NR_FACE_POINTS        Dimensions
#define NR_Simplex_POINTS (Dimensions + 1)
*/
#define TRACK_DELAUNAY_REGION_NEIGHBORHOOD 1

#define LIMITED_NEIGHBOR_DISTANCE	3

#define POINT_POSITION_PRECISION	(1e-9)

//#define USE_AGENT

/** \addtogroup Geometry
 *  @{
 */

// TYPES
typedef float 		POSITION_T;
typedef double 		DISTANCE_T;
typedef long double	HIGH_PRECISION_DISTANCE_T;

// DATA STRUCTURES
template <int Dimensions>
class Vertex;
template <int Dimensions, class VertexType = Vertex<Dimensions> >
class Simplex;
template <int Dimensions>
class Domain;


/**
 * \class Triangulation
 * \brief Delaunay triangulation of a set of vertices
 *
 * \author Nick Jagiella <nick.jagiella@googlemail.com>
 * \date April 2012
 *
 * This class template represents an Delaunay triangulation of a given set
 * of points. It provides functions to add and to remove points from an
 * initial triangulation of a fixed size.
 *
 * \tparam Dimensions The number of spatial dimensions. This is 2 by
 * default, but can also be set to another value <b>between 1 and 3</b>.
 */
template <int Dimensions = 2, class VertexType = Vertex<Dimensions>, class SimplexType = Simplex<Dimensions, VertexType> >
class Triangulation {
private:
	// methodes
	void sortTetrahedra();
	void setConvexHullNeighborhood();
	
	void setFramePoints1D();
	void setFramePoints2D();
	void setFramePoints3D();
	POSITION_T getVoronoiCellVolume2D(VertexType *&vertex, POSITION_T *&interfaces);
	POSITION_T getVoronoiCellVolume3D(VertexType *&vertex, POSITION_T *&interfaces);
	void getVoronoiCellInterfaces3D(VertexType *&vertex, POSITION_T *&interfaces );
	void printVoronoiDiagramToPovray2D(const char*, Domain<Dimensions>*);
	void printVoronoiDiagramToEps2D(const char*);
	void printVoronoiDiagramToPovray3D(const char*, Domain<Dimensions>*);
	void printVoronoiDiagramToVTK3D(const char*, Domain<Dimensions>* = 0);

	static Triangulation* newVoronoiDiagram();
	static Triangulation* newVoronoiDiagram( int x, int y, int z);
	static Triangulation* newVoronoiDiagramFromFile( char* filename);

	void Initialize();

	void setPointsOnSphericalSurface( double radius, int N, double variance, double center[Dimensions]);

	void getConvexHull( double thresholdDistance);

	void setSimplexNeighbors();

	int readExtendedNeighborhoodToFile( char* filename, int radius);
	int writeExtendedNeighborhoodToFile( char* filename, int radius);

	bool SimplexContainsFramePoints( SimplexType* tet);
	//double getVoronoiCellVolume(VertexType *v);
	//typename Triangulation<Dimensions>::Simplex*
	//getSimplexContainingPointInCircumSphere( Vertex* point);

	
	// frame
	int countFramePoints;
	VertexType **framePoints;

	// tetrahedra
	int countTetrahedra;
	int maxTetrahedra;
	SimplexType **tetrahedra;

	// voronoi grid point
	char voronoiGridSet;
	VertexType *voronoiGridPoints;

	// link between voronoi grid and voronoi cells
	VertexType ***voronoiCellToVoronoiGridPoints;
	int        *voronoiCellCountVoronoiGridPoints;

	// faces
	int countFaces;
	int maxFaces;
	VertexType ***faces;
	
	// convex faces
	int countConvexFaces;
	int maxConvexFaces;
	VertexType ***convexFaces;

protected:
	// voronoi cell
	int countVertices;
	int maxVertices;
	VertexType **vertices;

	// domain
	char domainSet;
	float boundaryThickness;
	Domain<Dimensions> *domain;
	//float xMin[Dimensions];
	//float xMax[Dimensions];
	//int xN[Dimensions];
	
public:
	void printPoints( const char * title, VertexType ** allPoints, int api );
	void printTetrahedraNeighbors();
	static VertexType *searchClosestVoronoiCell( VertexType *explorationCenter, double targetPosition[Dimensions]);
	VertexType *searchClosestFreeVoronoiCell( VertexType *explorationCenter);
	//static Vertex *searchForVoronoiCell( int countIgnoredPoints, Vertex** ignoredPoints, int countSourcePoints, Vertex** sourcePoints, int &countExploredPoints, Vertex** exploredPoints);

	/** \brief creates an empty triangulation */
	Triangulation();
	Triangulation(const char*);
	/** \brief creates an empty 1-dimensional triangulation for a fix domain size*/
	Triangulation( double xmin, double xmax);
	/** \brief creates an empty 2-dimensional triangulation for a fix domain size*/
	Triangulation( double xmin, double xmax, double ymin, double ymax);
	/** \brief creates an empty 3-dimensional triangulation for a fix domain size*/
	Triangulation( double xmin, double xmax, double ymin, double ymax, double zmin, double zmax);
	/** \brief creates an 3-dimensional triangulation of points spread homogeneously among and randomly within the cubes of an underlying square lattice*/
	Triangulation( int x, int y, int z);
	~Triangulation();

	/** \brief sets the minimal extension of the domain containing all added vertices */
	void setDomain();
	/** \brief sets a domain of a defined size */
	void setDomain( double xmin, double xmax, double ymin, double ymax, double zmin, double zmax);
	/** \brief sets a domain of a defined size */
	void setDomain( double min[Dimensions], double max[Dimensions]);
	/** \brief sets the frame points and super-simplices containing the domain */
	void setFramePoints();

	/** \brief addes a vertex (without triangulating it) */
	void add( VertexType* v);
	/** \brief addes and triangulates a vertex */
	void addAndTriangulate( VertexType* v);
	/** \brief triangulates vertex */
	void triangulate( VertexType* v);
	/** \brief triangulates all added vertices */
	void triangulate();

	/** \brief removes a vertex (without triangulating it) */
	void remove( VertexType* v);
	/** \brief removes and detriangulates a vertex */
	void removeAndDetriangulate( VertexType* v);
	/** \brief detriangulates a vertex */
	void detriangulate( VertexType* v);

	/** \brief returns the vertex with index i */
	VertexType* get( int i) { return vertices[i]; };

	int numberOfVertices();

	void setVoronoiGrid();

	/** \brief exports the triangulation to a file in Povray format*/
	void printToPovray( const char *filename, bool printPoints = true,  bool printFramePoints = true, bool printNeighborship = true, bool printTetrahedra = true, bool printConvexHull = false);
	void printToVTK( const char *filename);
	/** \brief exports the triangulation to a plain text file*/
	void printVoronoiDiagramToPovray( const char * filename, Domain<Dimensions>* subDomain = 0);
	/** \brief exports the triangulation to a file in EPS format*/
	void printVoronoiDiagramToEps( const char * filename);
	/** \brief exports the triangulation to a file in VTK format*/
	void printVoronoiDiagramToVTK( const char * filename, Domain<Dimensions>* subDomain = 0);

	void checkDelaunayCondition();
	VertexType *getCentralCell();
	SimplexType *getSimplexContainingPointInCircumSphere( VertexType* point);
	SimplexType *getSimplexContainingPoint( VertexType* point);
	VertexType *searchForVoronoiCell(
			int countIgnoredPoints,
			VertexType** ignoredPoints,
			int countSourcePoints, VertexType** sourcePoints,
			int &countExploredPoints,
			VertexType** exploredPoints);

	POSITION_T getVoronoiCellVolume(VertexType *&v);
	POSITION_T getVoronoiCellVolume(VertexType *&v, POSITION_T *&interfaces);
	void getVoronoiCellInterfaces(VertexType *&vertex, POSITION_T *&interfaces );

};

class Node{

};

/**
 * \class Vertex
 * \brief Vertex of the triangulation.
 *
 * \author Nick Jagiella <nick.jagiella@googlemail.com>
 * \date April 2012
 *
 * This class template represents an Delaunay triangulation of a given set
 * of points. It provides functions to add and to remove points from an
 * initial triangulation of a fixed size.
 *
 * \tparam Dimensions The number of spatial dimensions. This is 2 by
 * default, but can also be set to another value <b>between 1 and 3</b>.
 */
template <int Dimensions = 2>
class Vertex : public Node{
	friend class Simplex<Dimensions, Vertex<Dimensions> >; // class BinaryTree can now access data directly
	friend class Triangulation<Dimensions>; // class BinaryTree can now access data directly

protected:
	// spatial position
	POSITION_T position[Dimensions];

private:
	int index;

	// neighborhood
	bool neighborsInitialized;
	int countNeighbors;
	Vertex<Dimensions> **neighbors;

	void Initialize();

public:
	// methodes
	Vertex();
	Vertex(double x);
	Vertex(double x, double y);
	Vertex(double x, double y, double z);
	~Vertex();

	int addNeighbor( Vertex* neighborCell);
	int removeNeighbor( Vertex* neighborCell);
	int numberOfNeighbors();
	//int numberOfNeighbors() { return countNeighbors;};
	Vertex<Dimensions>* getNeighbor( int i);

	POSITION_T* pos();

	POSITION_T x( int dimension);
	POSITION_T x();
	POSITION_T y();
	POSITION_T z();

	void setx(POSITION_T, int);
	void setx(POSITION_T);
	void sety(POSITION_T);
	void setz(POSITION_T);

	int getIndex();
	void setIndex(int);

	double getDistanceTo( double[Dimensions]);
	double getDistanceTo( Vertex<Dimensions>* cell);
	double getDistanceTo( Simplex<Dimensions>* simplex);
	template <class SimplexType>
	double getDistanceToCircumsphere( SimplexType* tet);
	template <class SimplexType>
	double getDistanceToCircumsphereCenter( SimplexType* tet);

	bool isDomainBorder( Triangulation<Dimensions> *vd);

};


// Simplex type
template <int Dimensions, class VertexType >
class Simplex {
	friend class Triangulation<Dimensions, VertexType>; // class BinaryTree can now access data directly
	friend class Vertex<Dimensions>; // class BinaryTree can now access data directly
private:
	int index;

	// vertices
	VertexType *vertices[Dimensions+1];

	// neighbor tetrahedra
	int countNeighborTetrahedra;
	Simplex *neighborTetrahedra[Dimensions+1];

	// circumsphere
	char circumsphereInitialized;
	double circumsphere[Dimensions+1];
	double radius;
public:
	static int countUsers;
	static double **_A;
	static double **_B;
	static double *_x;

public:
	// methodes
	Simplex();

	int addSimplexNeighbor( Simplex<Dimensions,VertexType>* neighborTet);
	int addOrReplaceSimplexNeighbor( Simplex<Dimensions,VertexType>* neighborTet);
	int removeSimplexNeighbor( Simplex<Dimensions,VertexType>* neighborTet);
	Simplex<Dimensions,VertexType>* getNeighbor( int i);
	int countNeighbors() {return countNeighborTetrahedra;};

	void getCircumCenter( POSITION_T center[Dimensions]);
	double getCircumsphereRadius();
	//bool contains( VertexType * vc);
	bool contains( Vertex<Dimensions> * vc);
	void initCircumsphereRadius();
	int countCommonPointOfTetrahedra( Simplex<Dimensions,VertexType>* tetB);
	void getCircumsphereRadius( long double *circumsphere, long double &radius);
	int SimplexContainsFace( VertexType** face);
	void printSimplexNeighbors();

	VertexType*& getVertex( int i);
	void setVertex(VertexType*&, int);
	void setIndex(int);
	int getIndex() { return index;};

private:
	void getCircumCenter2D( POSITION_T center[2]);
	void getCircumCenter3D( POSITION_T center[3]);
};

template <int Dimensions>
class Domain {
private:
	POSITION_T xMin[Dimensions];
	POSITION_T xMax[Dimensions];
	int xN[Dimensions];
public:
	Domain( POSITION_T xmin, POSITION_T xmax, POSITION_T ymin=0, POSITION_T ymax=0, POSITION_T zmin=0, POSITION_T zmax=0);
	Domain( POSITION_T min[Dimensions], POSITION_T max[Dimensions]);
	//~Domain();

	POSITION_T min( int dimension);
	POSITION_T max( int dimension);

	bool contains( POSITION_T pos[Dimensions]);
	bool contains( Vertex<Dimensions>* v);
};

/** @}*/

#include "Triangulation.ipp"

#endif
