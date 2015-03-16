/*
 * Cells.hpp
 *
 *  Created on: Nov 26, 2012
 *      Author: jagiella
 */

#ifndef CELLS_HPP_
#define CELLS_HPP_

#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#define PI 3.142

class Cells{
public://private:
	enum types { cell, barrier};

    int Nmax;
    int N;
    float **position;
    float **speed;
    float *radius;
    float *cycleTime;
    float *maxCycleTime;
    bool  *overlapp;
    int   *type;

    float minRadius, maxRadius;
    float E;//=450;     // 480 Pa
    float eta;//=0.0; // 0...0.5
    float pressureThreshold;//=10; // Pa
    int   m;//=1; // 1=exponential, 2...99=erlang, 100=dirac


public:
    Cells( int N, float minRadius, float maxRadius, float E=450, float eta=0.0, float pressureThr=50, int m=1)
	: Nmax( N) , N(0), E(E), eta(eta), pressureThreshold(pressureThr), m(m) {
    	this->minRadius = minRadius;
    	this->maxRadius = maxRadius;

        position = (float**)malloc( Nmax*sizeof(float*));
        for( int i=0; i<Nmax; i++)
          	 position[i] = (float*)malloc( 3*sizeof(float));

		speed = (float**)malloc( Nmax*sizeof(float*));
		for( int i=0; i<Nmax; i++)
			speed[i] = (float*)malloc( 3*sizeof(float));

		radius = (float*)malloc( Nmax*sizeof(float));
		cycleTime = (float*)malloc( Nmax*sizeof(float));
		maxCycleTime = (float*)malloc( Nmax*sizeof(float));
		overlapp = (bool*)malloc( Nmax*sizeof(bool));
		type = (int*)malloc( Nmax*sizeof(int));
    }

    void add( float x, float y){
    	type[N] = cell;

    	position[N][0] = x;
    	position[N][1] = y;

    	radius[N] = 0.;
		cycleTime[N] = 0.;
		//maxCycleTime[N] = exp(- RAND01) * 24.;

    	N++;
    }

    inline float HertzForce( int i, int j){
    	//float E=450;     // 480 Pa
    	//float eta=0.4; // 0...0.5

    	float Estar = 1. / ((1. - eta*eta)/E + (1. - eta*eta)/E);
    	float Rstar = 1. / ((1./radius[i])   + (1./radius[j])  );

    	float distance = sqrt( pow(position[i][0]-position[j][0],2) + pow(position[i][1]-position[j][1],2));
    	float overlapp = radius[i] + radius[j] - distance;
    	if(overlapp<0){
       		return 0;
    	}else{
    		float Frep = 4./3. *Estar * sqrt(Rstar) * pow(overlapp, 3./2.);
    		//float Fadh = 0;//-2*PI*Rstar*10;
    		float Fadh = -2*PI*Rstar*10;
    		if( type[i]==barrier && type[j]==barrier){
    			Fadh = 0;
    			Frep = 0;
    		}else if(type[i]==barrier || type[j]==barrier){
    			Fadh = 0;
    		}

    		return Frep + Fadh;
    	}
    }

    void Acceleration( float *acceleration, int i, int j)
    {
    	float F = HertzForce( i, j);
       	float distance = sqrt( pow(position[i][0]-position[j][0],2) + pow(position[i][1]-position[j][1],2));
       	float mass = 1;

       	//if( type[i] != barrier)
       	acceleration[0] = - F/mass * (position[j][0]-position[i][0]) / distance;
       	//if( type[j] != barrier)
       	acceleration[1] = - F/mass * (position[j][1]-position[i][1]) / distance;

       	if( isnan(acceleration[0]) || isnan( acceleration[1]) || isinf(acceleration[0]) || isinf( acceleration[1])){
       		fprintf( stderr, "positions %i,%i = (%e,%e), (%e,%e) \n", i, j, position[i][0],position[i][1],position[j][0],position[j][1]);
       		fprintf( stderr, "distance |%i,%i| = %e\n", i, j, distance);
       		fprintf( stderr, "(Hertz) force = %e\n", F);
       		fprintf( stderr, "acceleration = (%e, %e)\n", acceleration[0], acceleration[1]);
       	}
    	return ;
    }


    void Acceleration( float *acceleration, int i, int j, float F)
    {
    	;
       	float distance = sqrt( pow(position[i][0]-position[j][0],2) + pow(position[i][1]-position[j][1],2));
       	float mass = 1;

       	//if( type[i] != barrier)
       	acceleration[0] = - F/mass * (position[j][0]-position[i][0]) / distance;
       	//if( type[j] != barrier)
       	acceleration[1] = - F/mass * (position[j][1]-position[i][1]) / distance;

       	if( isnan(acceleration[0]) || isnan( acceleration[1]) || isinf(acceleration[0]) || isinf( acceleration[1])){
       		fprintf( stderr, "distance |%i,%i| = %e\n", i, j, distance);
       		fprintf( stderr, "force = %e\n", F);
       		fprintf( stderr, "acceleration = (%e, %e)\n", acceleration[0], acceleration[1]);
       	}

    	return ;
    }


     float Pressure( int i, int j)
     {

    	 float Estar = 1. / ((1. - eta*eta)/E + (1. - eta*eta)/E);
    	 float Rstar = 1. / ((1./radius[i])   + (1./radius[j])  );

    	 float distance = sqrt( pow(position[i][0]-position[j][0],2) + pow(position[i][1]-position[j][1],2));
    	 float d = radius[i] + radius[j] - distance;

    	 float F = 4./3. *Estar * sqrt(Rstar) * pow(d, 3./2.);

    	 return 1./PI * pow( F*Estar*Estar/Rstar/Rstar, 1./3.);
     }
 };


#endif /* CELLS_HPP_ */
