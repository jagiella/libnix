/*
 * Boxes.hpp
 *
 *  Created on: Nov 27, 2012
 *      Author: jagiella
 */

#ifndef BOXES_HPP_
#define BOXES_HPP_

#include <stdio.h>
#include <assert.h>

template <class T>
class Boxes{
public:
	T*  ***boxes; // array of lists of all elements in box
	int ***N;     // array of number of elements in box
	int ***NMAX;  // array of max. number of elements fitting into box
	float latticeConstant;
	int X;
	int Y;
	int Z;

	double *molecules;
	int   nbmolecules;

	Boxes(float X, float Y, float Z, float l)
	{
		this->X=(int)ceil(X/l);
		this->Y=(int)ceil(Y/l);
		this->Z=(int)ceil(Z/l);
		latticeConstant = l;
		init();
	};

	Boxes(int nX, int nY, int nZ, float l){
		this->X=nX;
		this->Y=nY;
		this->Z=nZ;
		latticeConstant = l;
		init();
	};

	void init(){
		boxes = (T* ***)  malloc( X*sizeof(T* **));
		N     = (int***)  malloc( X*sizeof(int**));
		NMAX  = (int***)  malloc( X*sizeof(int**));
		for( int x=0; x<X; x++){
			boxes[x] = (T* **) malloc( Y*sizeof(T* *));
			N[x]     = (int**) malloc( Y*sizeof(int*));
			NMAX[x]  = (int**) malloc( Y*sizeof(int*));
			for( int y=0; y<Y; y++){
				boxes[x][y] = (T* *) malloc( Y*sizeof(T* ));
				N[x][y]     = (int*) malloc( Y*sizeof(int));
				NMAX[x][y]  = (int*) malloc( Y*sizeof(int));
				for( int z=0; z<Z; z++){
					N[x][y][z]=0;
					NMAX[x][y][z]=10;
					boxes[x][y][z] = (T*) malloc( NMAX[x][y][z]*sizeof(T));
				}
			}
		}
	}

	bool inside(float* &pos){
		return ( getIndex( pos[0]) >=0 && getIndex( pos[0]) < X
			  && getIndex( pos[1]) >=0 && getIndex( pos[1]) < Y
			  && getIndex( pos[2]) >=0 && getIndex( pos[2]) < Z);
	}

	/*Boxes(int X, int Y, float l){
		this->X=X;
		this->Y=Y;
		latticeConstant = l;
		boxes = (int***) malloc( X*sizeof(int**));
		N     = (int**)  malloc( X*sizeof(int*));
		NMAX  = (int**)  malloc( X*sizeof(int*));
		for( int x=0; x<X; x++){
			boxes[x] = (int**) malloc( Y*sizeof(int*));
			N[x]     = (int*)  malloc( Y*sizeof(int));
			NMAX[x]  = (int*)  malloc( Y*sizeof(int));
			for( int y=0; y<Y; y++){
				N[x][y]=0;
				NMAX[x][y]=10;
				boxes[x][y] = (int*) malloc( NMAX[x][y]*sizeof(int));
			}
		}
	};*/

	int getIndex( float x)
	{ return (int)floor(x/latticeConstant);}

	bool add( float *x, T index){
		return add( x[0], x[1], x[2], index);
	}

	bool add( float x, float y, float z, T index){
		int xi = (int)getIndex( x);
		int yi = (int)getIndex( y);
		int zi = (int)getIndex( z);

		assert(xi<X);
		assert(yi<Y);
		assert(zi<Z);
		assert(xi>=0);
		assert(yi>=0);
		assert(zi>=0);

		if( xi<0 || yi<0 || zi<0 || xi>=X || xi>=Y || zi>=Z){
			fprintf(stderr, "Error adding: element at (%i,%i) or (%f, %f) lies outside domain!\n", xi,yi, x,y);
			//exit(0);
			return false;
		}

		//fprintf(stderr, "indices: %i,%i\n", xi,yi);
		if( N[xi][yi][zi]==NMAX[xi][yi][zi]){
			//fprintf(stderr, "Error adding: %i at (%i,%i)\n", index, xi,yi);
			NMAX[xi][yi][zi]+=10;
			boxes[xi][yi][zi]     = (T*)  realloc(boxes[xi][yi][zi], NMAX[xi][yi][zi]*sizeof(T));
		}

		boxes[xi][yi][zi][N[xi][yi][zi]] = index;
		N[xi][yi][zi]++;

		return true;
	};

	void remove( float *x, T index){
		remove( x[0], x[1], x[2], index);
	};

	void remove( float x, float y, float z, T index){
		int xi = (int)getIndex( x);
		int yi = (int)getIndex( y);
		int zi = (int)getIndex( z);

		assert(xi<X);
		assert(yi<Y);
		assert(zi<Z);
		assert(xi>=0);
		assert(yi>=0);
		assert(zi>=0);

		for( int n=0; n<N[xi][yi][zi]; n++){
			if( boxes[xi][yi][zi][n]==index ){
				N[xi][yi][zi]--;
				//fprintf(stderr, "Removing at (%i,%i,%i)\n", xi,yi,zi);
				boxes[xi][yi][zi][n] = boxes[xi][yi][zi][ N[xi][yi][zi] ];
				return;
			}
		}
		fprintf(stderr, "Error removing at (%i,%i,%i): element not found!\n", xi,yi,zi);
		//exit(0);
		return;
	};


	T* get( float x, float y, float z){
		return boxes[getIndex(x)][getIndex(y)][getIndex(z)];
	};

	int length( float x, float y, float z){
		return N[getIndex(x)][getIndex(y)][getIndex(z)];
	};
	int length( int xi, int yi, int zi){
		return N[xi][yi][zi];
	};
};

template <class T>
class BoxIterator{
	// setting of iterator
	Boxes<T> *boxes;
	int xc, yc, zc;
	int depth;

	// state
	int xi, yi, zi; // actual box
	int i;		// actual item in box

	inline int MIN( int a, int b) { return (a<b ? a : b);}
	inline int MAX( int a, int b) { return (a>b ? a : b);}
public:
	BoxIterator( Boxes<T>* boxes, float *x, int depth=0) : boxes(boxes), depth(depth){
		init( x);
	};

	BoxIterator( Boxes<T>* boxes, int depth=0) :
		boxes(boxes), xc(0), yc(0), zc(0), depth(depth), xi(0), yi(0), zi(0), i(0)
	{};

	void init( float *x){
		// central box
		xc = (int)boxes->getIndex( x[0]);
		yc = (int)boxes->getIndex( x[1]);
		zc = (int)boxes->getIndex( x[2]);

		// starting box
		xi = MAX(xc-depth,0);
		yi = MAX(yc-depth,0);
		zi = MAX(zc-depth,0);

		// starting element
		i = 0;

		// look for first non empty box (if necessary)
		if( boxes->N[xi][yi][zi] == 0)
			operator++(1);
	};

	bool end()
	{ 	/*fprintf( stderr, "ROI boundary: %i (%i + %i)\n", xc+depth, xc, depth);
		fprintf( stderr, "box boundary: %i\n", boxes->X-1);
		fprintf( stderr, "index:        %i\n", xi);*/

			return ( xi > MIN(xc+depth,boxes->X-1));
	}

	bool operator++(int){
		i++;
		if( i >= boxes->N[xi][yi][zi]){ // no more elements in box?
			i = 0;

			// change box in z-direction
			zi++;
			if( zi > MIN(zc+depth,boxes->Z-1)){ // no more boxes in z-direction?
				zi = MAX(zc-depth,0);

				// change box in y-direction
				yi++;
				if( yi > MIN(yc+depth,boxes->Y-1)){ // no more boxes in y-direction?
					yi = MAX(yc-depth,0);

					// change box in x-direction
					xi++;
					if( xi > MIN(xc+depth,boxes->X-1)) // no more boxes ?
						return false;
				}
			}
		}


		// empty box?
		if( boxes->N[xi][yi][zi] == 0)
			return operator++(1);

		//fprintf( stderr, "box(%i,%i), element=%i/%i\n", xi, yi, i, boxes->N[xi][yi]);

		return true;//boxes->boxes[xi][yi][i];
	}
	//int operator==( int index);
	T operator*() { return boxes->boxes[xi][yi][zi][i]; }
};

#endif /* BOXES_HPP_ */
