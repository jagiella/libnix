/*
 * SimulationThread.hpp
 *
 *  Created on: Nov 26, 2012
 *      Author: jagiella
 */

#include "Cells.hpp"
#include "Agent.hpp"
#include "Boxes.hpp"
#include "Spring.hpp"
#include "Kinetics.hpp"
#include "Macros.hpp"
#include "Covariance.hpp"

#include <stdio.h>
#include <math.h>

#include <time.h>
#include <clocale>

#include <gsl/gsl_math.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_statistics_double.h>



#define mexPrintf(format, args...) fprintf(stderr, format , ## args)
#define eprintf(format, args...) fprintf(stderr, format , ## args)
#define OVERDAMPED
#define HERTZ_KINETICS

#define r0 60.
#define Z0 150.
#define l1 0.25
#define l2 0.1
#define CRYPT_RADIUS( z ) (( r0*(1.+l1*cos(2*M_PI*(Z0-(z))/Z0+M_PI/2.)) * sqrt(2./M_PI*acos(pow((Z0-(z))/Z0, l2)))  ))

//float m=1;


float enhance( float x, float critical_x, int sharpness){

	return 1. / (1. + pow(critical_x/x,sharpness));
}

float inhibit( float x, float critical_x, int sharpness){

	return 1. / (1. + pow(x/critical_x,sharpness));
}


int DeltaNotchModel (double t, const double y[], double dydt[], void * params)
{
	AgentList<> *al    = (AgentList<>*) ((int**)params)[0];
	Boxes<Agent*> *box = (Boxes<Agent*>*) ((int**)params)[1];

	//fprintf( stderr, "[N=%i]\n", al->size());

	fprintf( stderr, "[N=%i]\n", box->X*box->Y*box->Z);
	int delta = box->X*box->Y*box->Z+0;
	int notch = box->X*box->Y*box->Z+1;


	// Agent-Based
	for( int i=0; i<al->size(); i++){ // FOR ALL AGENTS
		float delta_neighbors = 0;
		int   count_neighbors = 0;

		BoxIterator<Agent*> boxIt( box, al->at(i)->position(), 1);
		for( ; !boxIt.end(); boxIt++){
			Agent *a = *boxIt;
			int j = a->index();
			if( i!=j && al->inContact(i,j)){
				delta_neighbors += y[2*j+delta];
				count_neighbors ++;
			}
		}

		if( count_neighbors)
			delta_neighbors /=  count_neighbors;
		//fprintf(stderr, "<<< nei: %e (%i) >>>\n", delta_neighbors, count_neighbors);

		float k=0.1, h=2;
		dydt[i*2+notch] = enhance( delta_neighbors, k, h) - y[i*2+notch]; // notch
		dydt[i*2+delta] = inhibit( y[i*2+notch],    k, h) - y[i*2+delta]; // delta
	}

	// Diffusion
	//2*al->size() + box->X*box->Y*box->Z
	//int M = 2*al->size();
	int M=0;
	int dx=box->Y*box->Z;
	int dy=box->Z;
	int dz=1;
	int j;
	for( int x=0; x<box->X; x++ )
		for( int _y=0; _y<box->Y; _y++ )
			for( int z=0; z<box->Z; z++ )
			{
				j = M + x*dx + _y*dy + z*dz;
				dydt[j] = 0;

				// BOUNDARY?
				if( x==0 || _y==0 ||  x==box->X-1 || _y==box->Y-1 ){
					// influx
					dydt[j] += 1-y[j];


				}else
				{
					// diffusion
					if( x-1>=0) dydt[j] += -y[j]+y[j-dx];
					if(_y-1>=0) dydt[j] += -y[j]+y[j-dy];
					//if( z-1>=0) dydt[j] += -y[j]+y[j-dz];
					if( x+1<box->X) dydt[j] += -y[j]+y[j+dx];
					if(_y+1<box->Y) dydt[j] += -y[j]+y[j+dy];
					//if( z+1<box->Z) dydt[j] += -y[j]+y[j+dz];
					dydt[j] *= 0.1;
				}
					// reaction
					//Agent **a = box->get( x*box->latticeConstant,_y*box->latticeConstant,z*box->latticeConstant);
					int count = box->length( x*box->latticeConstant,_y*box->latticeConstant,z*box->latticeConstant);
					dydt[j] += -count*y[j]; // uptake




			}
	return 1;
}

void init( double *array, int n, const char* string)
{
	char* str = (char*) string;
	for( int i=0; i<n; i++)
		array[i] = strtof( str, &str);
}


int Caspase8Model (double t, const double y[], double dydt[], void * params)
{
	AgentList<> *al    = (AgentList<>*) ((int**)params)[0];
	Boxes<Agent*> *box = (Boxes<Agent*>*) ((int**)params)[1];

	//fprintf( stderr, "[ N=%i ]\n", al->size());

	int M = box->nbmolecules*box->X*box->Y*box->Z;

	double k[12]; for( int i=0; i<11; i++) k[i]=1;
	double v[14];

	k[  1] =        8.12e-4; // kon,FADD
	k[  2] =        0.00567; // koff,FADD
	k[  3] =        0.000492; // kon,p55
	k[  4] =        0.0114; // kcl,D216
	k[  5] =        4.47e-4; // kcl,D374,trans,p55
	k[  6] =        0.00344; // kcl,D374,trans,p43
	k[  7] =        0.0950; // kp18,inactive
	k[  8] =        0.000529; // kcl,BID
	k[  9] =        0.00152; // kcl,probe
	k[ 10] =        8.98; // KD,R
	k[ 11] =        15.4; // KD,L

	// Agent-Based
	for( int i=0; i<al->size(); i++){
		// FOR ALL AGENTS
		double *dxdt =       &dydt[M + i*al->nbmolecules];
		double *x = (double *)&y  [M + i*al->nbmolecules];
		double a[2], u[2];


		/* Driving input */
		//u[1]     = 16.6;//(t<500 ? 0 : t);
		{
			int _x= box->getIndex(al->at(i)->position()[0]);
			int _y= box->getIndex(al->at(i)->position()[1]);
			int z= box->getIndex(al->at(i)->position()[2]);
			int dx=box->Y*box->Z;
			int dy=box->Z;
			int dz=1;
			int j =_x*dx + _y*dy + z*dz;

			double conc = box->molecules[ j ];
			u[1]     = conc;//(t<500 ? 0 : t);
		}
		//x[1]     = 12;

		/* Explicit algebraic equations */
		a[ 1]    = ((x[1])*(x[1])*(x[1]))*((k[11])*(k[11]))*u[1]/((u[1] + k[11])*(((x[1])*(x[1]))*((k[11])*(k[11])) + k[10]*((u[1])*(u[1])) + 2.0*k[10]*k[11]*u[1] + k[10]*((k[11])*(k[11]))));
		//a[ 1]    = 0;

		/* Rates */
		v[ 1]    = k[1]*a[1]*x[2];
		v[ 2]    = k[2]*x[7];
		v[ 3]    = k[3]*x[3]*x[7];
		v[ 4]    = k[4]*x[8];
		v[ 5]    = k[4]*x[10];
		v[ 6]    = k[5]*x[8]*(x[8] + x[9]);
		v[ 7]    = k[6]*x[8]*x[10];
		v[ 8]    = k[5]*x[9]*(x[8] + x[9]);
		v[ 9]    = k[6]*x[9]*x[10];
		v[10]    = k[7]*x[11];
		v[11]    = k[8]*x[4]*(x[10] + x[11]);
		v[12]    = k[9]*x[5]*(x[10] + x[11]);
		v[13]    = k[9]*x[6]*x[11];

		/* ODE */
		dxdt[ 0] = 0.0; // ignore -> start index = 1
		dxdt[ 1] = 0.0; // CD95
		dxdt[ 2] =  - v[1] + v[2]; // FADD
		dxdt[ 3] =  - v[3]; // p55free
		dxdt[ 4] =  - v[11]; // Bid
		dxdt[ 5] =  - v[12]; // PrNES-mCherry
		dxdt[ 6] =  - v[13]; // PrER-mGFP
		dxdt[ 7] = v[1] - v[2] - v[3] + v[5] + v[8] + v[9]; // DISC
		dxdt[ 8] = v[3] - v[4] - v[6] - v[7]; // DISCp55
		dxdt[ 9] = v[4] - v[8] - v[9]; // p30
		dxdt[10] =  - v[5] + v[6] + v[7]; // p43
		dxdt[11] = v[5] + v[8] + v[9] - v[10]; // p18
		dxdt[12] = v[10]; // p18inactive
		dxdt[13] = v[11]; // tBid
		dxdt[14] = v[12]; // PrNES
		dxdt[15] = v[12]; // mCherry
		dxdt[16] = v[13]; // PrER
		dxdt[17] = v[13]; // mGFP
	}

	// Diffusion
	int dx=box->Y*box->Z;
	int dy=box->Z;
	int dz=1;
	int j;
	for( int x=0; x<box->X; x++ )
		for( int _y=0; _y<box->Y; _y++ )
			for( int z=0; z<box->Z; z++ )
			{
				j = x*dx + _y*dy + z*dz;
				dydt[j] = 0;

				// BOUNDARY?
				if( x==0 || _y==0 ||  x==box->X-1 || _y==box->Y-1 ){
					// influx
					if( t>24*7)
						dydt[j] += 16.6-y[j];
					//else
					//	dydt[j] += 0-y[j];
				}else
				{
					// diffusion
					if( x-1>=0) dydt[j] += -y[j]+y[j-dx]*(!box->length( x,_y,z) && !box->length( x-1,_y,z) ? 10 : 1);
					if(_y-1>=0) dydt[j] += -y[j]+y[j-dy]*(!box->length( x,_y,z) && !box->length( x,_y-1,z) ? 10 : 1);
					//if( z-1>=0) dydt[j] += -y[j]+y[j-dz];
					if( x+1<box->X) dydt[j] += -y[j]+y[j+dx]*(!box->length( x,_y,z) && !box->length( x+1,_y,z) ? 10 : 1);
					if(_y+1<box->Y) dydt[j] += -y[j]+y[j+dy]*(!box->length( x,_y,z) && !box->length( x,_y+1,z) ? 10 : 1);
					//if( z+1<box->Z) dydt[j] += -y[j]+y[j+dz];
					dydt[j] *= 0.001;
				}
					// reaction
					//Agent **a = box->get( x*box->latticeConstant,_y*box->latticeConstant,z*box->latticeConstant);
					int count = box->length( x*box->latticeConstant,_y*box->latticeConstant,z*box->latticeConstant);
					dydt[j] += -count*y[j]; // uptake




			}
	return 1;
}




float erl(float mean, float m) {
	if (m >= 100)
		return mean;

	float value = 1;
	for (int i = 0; i < m; i++)
		value *= (1. - RAND01);
	return -log(value) * mean / m;
}


#ifndef NOGUI
#include <QtCore/QThread>

void simCrypt( AgentList<> *al, Boxes<Agent*> *box, float *time, bool verbosity, bool snapshots, int num_snapshots, double *parameters, int num_parameters);

class MyThread: public QThread {
private:
	Cells *cells;
	AgentList<> *al;
	Spring<int> **springs;
	int countSprings;
	//int boxes[1000][1000];
	Boxes<Agent*> *box;
	float *time;
	//static const float minRadius = 5;
	//static const float maxRadius = 7;

public:
	MyThread(AgentList<> *c, Boxes<Agent*> *b, float *t) {
		cells = 0;
		al = c;
		box = b;
		time = t;
		countSprings = 0;
		UNUSED(springs);


		// ADD BARRIER (LINE)
		/*for (int i = -100; i < 100; i++) {
			Agent *a = new Agent(i * 1 + b->X * b->latticeConstant / 2.,
					-0. + b->X * b->latticeConstant / 2.);
			a->set_type(2);
			a->radius() = 1.;

			box->add(a->position(), al->size());
			al->insert(a);
		}*/

		// ADD BARRIER (HALF-CIRCLE)
		/*float mR = 20*5;
		float my = b->Y * b->latticeConstant / 2. - mR, mx = b->X * b->latticeConstant / 2.;
		for (float angle = 3.142; angle > 0; angle -= 0.05) {
			Agent *a = new Agent(mx + cos(angle) * mR, my + sin(angle) * mR);
			a->set_type(2);
			a->radius() = 5.;

			box->add(a->position(), al->size());
			al->insert(a);
		}*/


		exit(0);
		//this->stop();
		return;

		// Add to Boxes
		//for( int a=0; a<al->size(); a++){
		//	box->add( al->at(a)->position(), a);
		//}
	}

/*	MyThread(Cells *c, Boxes<Agent*> *b, float *t) //: N(N)
			{
		cells = c;
		box = b;
		time = t;

		// ADD FIRST CELL
		float mx = (int) sqrt(cells->Nmax) * box->latticeConstant / 2;
		float my = (int) sqrt(cells->Nmax) * box->latticeConstant / 2;

		cells->add(mx, my);
		cells->radius[0] = cells->minRadius;
		cells->maxCycleTime[0] = erl(24, cells->m);
		box->add(cells->position[0][0], cells->position[0][1], 0);

		// ADD BARRIERS
		int previous = -1;
		int first;
		springs = (Spring<int>**) malloc(1000 * sizeof(Spring<int>*));
		countSprings = 0;
		for (float angle = 0; angle < 3.142; angle += 0.05) {
			printf("%lf %lf\n", sin(angle), cos(angle));
			cells->add(mx + sin(angle) * 200, my + cos(angle) * 200);
			cells->radius[cells->N - 1] = 5;
			cells->type[cells->N - 1] = Cells::barrier;
			box->add(cells->position[cells->N - 1][0],
					cells->position[cells->N - 1][1], cells->N - 1);

			if (previous >= 0) {
				springs[countSprings] = new Spring<int>(previous, cells->N - 1,
						1.);
				countSprings++;
			} else
				first = cells->N - 1;
			previous = cells->N - 1;
		}
		for (float angle = 3.142; angle > 0; angle -= 0.1) {
			printf("%lf %lf\n", sin(angle), cos(angle));
			cells->add(mx - sin(angle) * 200, my + cos(angle) * 200);
			cells->radius[cells->N - 1] = 5;
			cells->type[cells->N - 1] = Cells::barrier;
			box->add(cells->position[cells->N - 1][0],
					cells->position[cells->N - 1][1], cells->N - 1);

			if (previous >= 0) {
				springs[countSprings] = new Spring<int>(previous, cells->N - 1,
						1.);
				countSprings++;
			}
			previous = cells->N - 1;

		}

		springs[countSprings] = new Spring<int>(previous, first, 1.);
		countSprings++;

		fprintf(stderr, "Added first cell.\n");
	}
	;
*/
	void run() {
		fprintf(stderr, "START SIM\n");

		simCrypt( al, box, time, true, true, 0, 0);
		//simMonolayer();
		//simCylinder();
		fprintf(stderr, "STOPPED SIM\n");
	};

	//void simCylinder();
	//void simMonolayer();
	//void simCrypt() { simCrypt( al, box, time);};
};

#endif //NOGUI

/*
void MyThread::sim() {
	float end = 10000;
	float dt = .0125; //0.01;

	float randVelocity = 0.1; //0.1;
	//float repulsiveVelocity = 10.;

	qsrand(0);
	srand(0);

	//float latticeConstant = 5;
	//boxes[(int)(cells->position[0][0]/latticeConstant)][(int)(cells->position[0][1]/latticeConstant)] = 0;

	FILE *fp = fopen("dataOnTheFly.dat", "w+");

	for (float t = 0; t < end && fp != 0; t += dt)
	// if(go)
			{
		*time = t;
		float maxDisplacement = 0;
		//float maxAcceleration = 0;
		//float maxAcceleration = 0;

		// random movement
		for (int i = 0; i < cells->N; i++)
		//if( cells->type[i] != Cells::barrier)
				{
			cells->speed[i][0] += randVelocity * (RAND01 - 0.5);
			cells->speed[i][1] += randVelocity * (RAND01 - 0.5);

			cells->overlapp[i] = false;
		}

		// interaction
//#pragma omp parallel
		for (int i = 0; i < cells->N; i++)
		//if( cells->type[i] != Cells::barrier)
				{
			//for( int j=0; j<cells->N; j++)
			//for( int j=i+1; j<cells->N; j++)

			/// BOXES
			int xi = box->getIndex(cells->position[i][0]);
			int yi = box->getIndex(cells->position[i][1]);
			int l = 0;
			for (int x = MAX(xi-1,0); x <= MIN(xi+1,box->X-1); x++)
				for (int y = MAX(yi-1,0); y <= MIN(yi+1,box->Y-1); y++)
					for (l = 0; l < box->N[x][y]; l++) {

						Agent *a = box->boxes[x][y][l];
						int j = a->index();

						if (j > i)
						//if( j!=i)
						/// BOXES
								{
							/ *float dx = cells->position[i][0]
									- cells->position[j][0];
							float dy = cells->position[i][1]
									- cells->position[j][1];* /

							// acceleration
							float acceleration[2];
							cells->Acceleration(acceleration, i, j);
							//if( cells->type[i] != Cells::barrier)
							{
								cells->speed[i][0] += dt * acceleration[0];
								cells->speed[i][1] += dt * acceleration[1];
							}
							//if( cells->type[j] != Cells::barrier)
							{
								cells->speed[j][0] -= dt * acceleration[0];
								cells->speed[j][1] -= dt * acceleration[1];
							}
							//fprintf(stderr, "(%e, %e) vs (%e, %e)\n", dx,dy, acceleration[0], acceleration[1]);

							{
								//float dist = sqrt(pow( dx, 2) + pow( dy, 2));
								//if( dist < cells->radius[i]+cells->radius[j] )
								if (cells->Pressure(i, j)
										> cells->pressureThreshold) {
									cells->overlapp[i] = true;
									cells->overlapp[j] = true;
								}
							}
						}
					}
		}

		// SPRING FORCES
		//if(false)
		for (int s = 0; s < countSprings; s++) {

			int i = springs[s]->get(0);
			int j = springs[s]->get(1);

			float dist0 = springs[s]->distance();
			float dist1 = sqrt(
					pow(cells->position[i][0] - cells->position[j][0], 2)
							+ pow(cells->position[i][1] - cells->position[j][1],
									2));

			// acceleration
			float acceleration[2];
			/ *float springForce[2] = { (cells->position[j][0]
					- cells->position[i][0]) * (dist1 - dist0),
					(cells->position[j][1] - cells->position[i][1])
							* (dist1 - dist0) };* /
			cells->Acceleration(acceleration, i, j, -(dist1 - dist0) * 1000.);
			cells->speed[i][0] += dt * acceleration[0];
			cells->speed[i][1] += dt * acceleration[1];
			cells->speed[j][0] -= dt * acceleration[0];
			cells->speed[j][1] -= dt * acceleration[1];

			//fprintf( stderr, "(%i, %i) d = %lf -> a (%lf, %lf)\n", i, j, dist1, acceleration[0], acceleration[1]);

			if (dist1 > cells->radius[i] + cells->radius[j] && t > 10) {
				countSprings--;
				springs[s] = springs[countSprings];
				s--;
			}
		}

		// ESTIMATE TIME STEP
		maxDisplacement = 0.09;
		dt = 1;
		for (int i = 0; i < cells->N; i++)
		//if( cells->type[i] != Cells::barrier)
				{
			if (dt > fabs(maxDisplacement / cells->speed[i][0]))
				dt = fabs(maxDisplacement / cells->speed[i][0]);
			if (dt > fabs(maxDisplacement / cells->speed[i][1]))
				dt = fabs(maxDisplacement / cells->speed[i][1]);
		}

		// update
		bool leftDomain = false;
//#pragma omp parallel
		for (int i = 0; i < cells->N && !leftDomain; i++)
		//if( cells->type[i] != Cells::barrier)
				{
			// friction
			//cells->speed[i][0] -= cells->speed[i][0]*100*dt;
			//cells->speed[i][1] -= cells->speed[i][1]*100*dt;
			if (cells->type[i] != Cells::barrier) {
				cells->speed[i][0] -= cells->speed[i][0] * 0.1;
				cells->speed[i][1] -= cells->speed[i][1] * 0.1;
			} else {
				cells->speed[i][0] -= cells->speed[i][0] * 0.01;
				cells->speed[i][1] -= cells->speed[i][1] * 0.01;
			}

			// Move
			box->remove(cells->position[i][0], cells->position[i][1], i);

			cells->position[i][0] += dt * cells->speed[i][0];
			cells->position[i][1] += dt * cells->speed[i][1];

			if (!box->add(cells->position[i][0], cells->position[i][1], i)) {
				leftDomain = true;
				fprintf(stderr, "dt = %e\n", dt);
				fprintf(stderr, "speed = (%e, %e)\n", cells->speed[i][0],
						cells->speed[i][1]);
				fclose(fp);
				//return;
			}

			if (maxDisplacement < fabs(dt * cells->speed[i][0]))
				maxDisplacement = fabs(dt * cells->speed[i][0]);
			if (maxDisplacement < fabs(dt * cells->speed[i][1]))
				maxDisplacement = fabs(dt * cells->speed[i][1]);

			if (cells->type[i] != Cells::barrier) {
				// Cell cycle progression?
				if (cells->cycleTime[i] > 0 || !cells->overlapp[i]) {
					cells->cycleTime[i] += dt;
				}

				// Radius Grow?
				if (cells->cycleTime[i] > 0
						&& cells->radius[i] < cells->maxRadius
						|| !cells->overlapp[i]) {
					cells->radius[i] += dt
							* (cells->maxRadius - cells->minRadius)
							/ cells->maxCycleTime[i];
				}
			}
			//cells->position[i][0]
		}

		if (leftDomain)
			return;

		// division
		//if(false)
		{

			int Nold = cells->N;
			for (int i = 0; i < Nold; i++)
				if (cells->cycleTime[i] >= cells->maxCycleTime[i]
						&& cells->radius[i] >= cells->maxRadius) {

					if (cells->N == cells->Nmax) {
						fclose(fp);
						//break;
						return;
					}

					float direction[2] = { RAND01 - 0.5, RAND01 - 0.5 };
					float length = sqrt(
							pow(direction[0], 2.) + pow(direction[1], 2.));
					direction[0] /= length;
					direction[1] /= length;
					cells->position[cells->N][0] = cells->position[i][0]
							+ direction[0] * cells->maxRadius;
					cells->position[cells->N][1] = cells->position[i][1]
							+ direction[1] * cells->maxRadius;

					cells->speed[cells->N][0] = cells->speed[i][0];
					cells->speed[cells->N][1] = cells->speed[i][1];

					if (!box->add(cells->position[cells->N][0],
							cells->position[cells->N][1], cells->N)) {
						fclose(fp);
						return;
					}

					cells->radius[i] = cells->minRadius;
					cells->radius[cells->N] = cells->minRadius;
					cells->cycleTime[i] = 0;
					cells->cycleTime[cells->N] = 0;
					cells->maxCycleTime[i] = erl(24, cells->m);
					cells->maxCycleTime[cells->N] = erl(24, cells->m);

					cells->N++;
				}
		}

		fprintf(stderr, "\r%f => %i cells (dt=%e) maxDisplacement=%e \b", t,
				cells->N, dt, maxDisplacement);
		fprintf(fp, "%f %i\n", t, cells->N);

		/ *if( maxDisplacement>0.)
		 dt = cells->minRadius / (maxDisplacement / dt) / 2.;

		 if( dt>0.05)
		 dt=0.05;
		 * /
	}
	fclose(fp);
}
*/


void simCrypt( AgentList<> *al, Boxes<Agent*> *box, float *time, bool verbosity, bool write_snapshots, int num_snapshots, double *parameters, int num_parameters) {

	// INPUT
	float data_mean[box->Y][Agent::TypeSize];
	float data_std [box->Y][Agent::TypeSize];
	if( !write_snapshots){
		for( int pos=0; pos<box->Y; pos++)
			for( int typ=0; typ<Agent::TypeSize; typ++)
				data_mean[pos][typ] = data_std [pos][typ] = 0;
		FILE *fp_in;
		char line[1024], *p_word;
		int size =0;

		for( int ti=0; ti<num_snapshots; ti++){ // time points
			char filename[1024];
			sprintf( filename, "snapshot_t%i.dat", ti);
			fp_in = fopen( filename, "r");
			while( fgets(line, 1024, fp_in) != NULL){ // position
				 /* get a line, up to 80 chars from fr.  done if NULL */
				int pos = (int) strtof( line, &p_word);
				//printf ("pos = %i\n", pos);
				for( int typ=0; typ<Agent::TypeSize; typ++){ // cell types
					float value = strtof( p_word, &p_word);
					data_mean[pos][typ] += value;
					data_std [pos][typ] += value*value;
					//printf ("%f\n", data_mean[j][i]);
				}
				//printf ("\n");
			}
			fclose( fp_in);

			for( int pos=0; pos<box->Y; pos++){
				for( int typ=0; typ<Agent::TypeSize; typ++){
					data_mean[pos][typ] /= (float)num_snapshots;
					data_std [pos][typ] = 1e-6 + sqrt(data_std [pos][typ] / (float)num_snapshots - data_mean[pos][typ]*data_mean[pos][typ]);
					//printf ("%f\n", data_mean[pos][typ]);
				}
				//printf ("\n");
			}
		}
	}


	// PARAMETERS
	for( int i=0; i<num_parameters; i++)
		switch(i){
		case 0:
			Agent::TPwnt   = parameters[i]; break;
		case 1:
			Agent::TPnotch = parameters[i]; break;
		case 2:
			Agent::TDwnt = parameters[i]; break;
		}



	// OUTPUT SETTINGS
	int output_offset   = 200;
	int output_interval = 24;
	int output_duration = 24*num_snapshots;


	setlocale(LC_ALL,"C");

	al->setDimensions(3);

	float x0 = box->X * box->latticeConstant / 2.;
	float z0 = box->Z * box->latticeConstant / 2.;
	float y0 = 30;//box->Y * box->latticeConstant / 2.;//30;

	// ADD CELLS
	for (int i = 0; i < 1; i++) {
		Agent *a = new Agent(
				box->X * box->latticeConstant / 2. - RAND01,
				y0 + RAND01,
				z0);
		if(verbosity) printf("x=%e, y=%e, z=%e\n", a->position()[0], a->position()[1], a->position()[2]);

		a->set_type(Agent::StemCell);
		a->radius() = 5.1;

		box->add(a->position(), a);
		al->insert(a);
	}


	Agent::wntmin = y0+150;
	Agent::wntmax = y0;

	float dy = 0;
	float dy0 = 5;

	if(true)
	for (float y = 0; y<=150; y += dy) {

		// Crypt Radius
		float R = CRYPT_RADIUS(y);

		// Change in Crypt Radius
		float dRdy = ( CRYPT_RADIUS(y+0.0001)-CRYPT_RADIUS(y) )/0.0001;

		/*fprintf( stderr, "z=%e, R=%e, dRdz=(%e*%e=)%e=%e=%e\n", z,R,
				r0*(1.+l1*cos(2*M_PI*(Z0-z-1)/Z0+M_PI/2.)),
				sqrt(2./M_PI*acos(pow((Z0-z-1)/Z0, l2))),
				r0*(1.+l1*cos(2*M_PI*(Z0-z-1)/Z0+M_PI/2.))*
				sqrt(2./M_PI*acos(pow((Z0-z-1)/Z0, l2))),
				CRYPT_RADIUS((z+1.)),
				dRdz);*/

		// Step in z-Direction
		if( !isnan(dRdy)) dy = dy0 / sqrt(fabs( dRdy ) + 1);
		if( y<150 && y+dy >150) dy = 150 - y;

		//fprintf( stderr, "y=%e\n", box->Y * box->latticeConstant / 2 + z - 150/2.);

		for( float angle=0; angle<2*M_PI; angle+=dy0/R){
			//fprintf( stderr, "dangle=%e\n",dangle);
			      // angular position
			float cell_x = x0 + sin(angle)*R,
				  // position along crypt axis
				  cell_y = y0 + y,
				  // angular position
				  cell_z = z0 + cos(angle)*R;

			Agent *a = new Agent(
					cell_x,
					cell_y,
					cell_z);
			a->set_type( Agent::substrate);
			a->radius() = 5.;
			//eprintf( "BEFORE\n");
			box->add(a->position(), a);
			//eprintf( "AFTER\n");
			al->insert(a);
			//eprintf( "AFTER2\n");
		}


		//eprintf( "z=%e, dz=%e, r(z)=%e === dRdz=%e\n", z,dz, CRYPT_RADIUS(z), dRdz);

	}

	if(true)
	for (float x = 0; x<box->X * box->latticeConstant; x += 5)
	for (float z = 0; z<box->Z * box->latticeConstant; z += 5)
	{
		float y = y0 + 150;
		float r = sqrt(	pow(x - x0, 2) + pow(z - z0, 2)	);
		if( r > CRYPT_RADIUS(150)){
			Agent *a = new Agent(
					x,
					y,
					z);
			a->set_type( Agent::substrate);
			a->radius() = 5.;
			//eprintf( "BEFORE\n");
			box->add(a->position(), a);
			//eprintf( "AFTER\n");
			al->insert(a);

		}
	}


	// SIMULATE
	srand(0);

	char MODEL = 3;

	if(verbosity)
		mexPrintf("SIMULATE (%i)\n", rand());
	float dt = 0.1;
	float tend = output_offset + output_duration;

	float randVelocity[Agent::TypeSize], *p_randVelocity = randVelocity;
	float baseVelocity = 1000;
	for( int typ=0; typ<Agent::TypeSize; typ++) randVelocity[typ] = baseVelocity;

	int dimy = 3;
	float minR = 5;
	float maxR = 7;

	float cycleTime[Agent::TypeSize];
	for( int typ=0; typ<Agent::TypeSize; typ++) cycleTime[typ] = 0.;

	switch( MODEL){
	case 1:
		cycleTime[Agent::StemCell]	 = 24.;
		cycleTime[Agent::Absorptive] = 24.;
		cycleTime[Agent::Enterocyte] = 24.;
		break;

	case 2:
		cycleTime[Agent::StemCell]	 = 24.;
		cycleTime[Agent::Absorptive] = 24.;
		cycleTime[Agent::Secretory]	 = 24.;
		break;

	case 3:
		cycleTime[Agent::StemCell]	 = 24.;
		cycleTime[Agent::Absorptive] = 12.;
		break;

	case 4:
		cycleTime[Agent::StemCell]	 = 24.;
		cycleTime[Agent::Absorptive] = 12.;
		break;

	}

	float **adhesionEnergies = allocMatrix<float>(Agent::TypeSize, Agent::TypeSize);
	initMatrix<float>( adhesionEnergies, Agent::TypeSize, Agent::TypeSize, 0);
	float sc = 2e1;//2e1;
	for( int typ=0; typ<Agent::TypeSize; typ++)
		adhesionEnergies[typ][Agent::substrate] = adhesionEnergies[Agent::substrate][typ] = sc;
	/*adhesionEnergies[Agent::StemCell  ][Agent::substrate] = adhesionEnergies[Agent::substrate][Agent::StemCell  ] = 1*sc;
	adhesionEnergies[Agent::Paneth    ][Agent::substrate] = adhesionEnergies[Agent::substrate][Agent::Paneth    ] = 1*sc;
	adhesionEnergies[Agent::Absorptive][Agent::substrate] = adhesionEnergies[Agent::substrate][Agent::Absorptive] = sc;
	adhesionEnergies[Agent::Secretory ][Agent::substrate] = adhesionEnergies[Agent::substrate][Agent::Secretory ] = sc;
*/

	al->pressureThreshold = 1e3;//1e3;

	//for (float t = 0; t < tend; t += dt)
	float t = 0;
	dt = 1e-10;

	// COMPUTATION TIME
	clock_t clck = clock();
	clock_t clck_rand = 0;
	clock_t clck_updpos = 0;
	clock_t clck_div = 0;
	clock_t clck_forces = 0;
	clock_t clck_timestep = 0;

	// OUTPUT
	int g_stats[box->Y][Agent::TypeSize];
	int g_count[box->Y];
	for( int pos=0; pos<box->Y; pos++){
		for( int t=0; t<Agent::TypeSize; t++)
			g_stats[pos][t] = 0;
		g_count[pos] = 0;
	}
	{
		FILE *fp = fopen( "profiles.dat", "w+");
		fclose(fp);
	}

	while( t<tend)
	{
		//mexPrintf("t=%e, dt=%e, #cells=%i\n", t, dt, al->size());

		// RESET ACCELERATION
		//printf("RESET ACCELERATION\n");
		resetAcceleration( al, dimy);


		// CELL-CELL-INTERACTION
		//printf("CELL-CELL-INTERACTION\n");
		addCellCellInteractionForces( al, dimy, box, adhesionEnergies);
		//updatePressure( al, dimy, box, adhesionEnergies);
		clck_forces += clock() - clck; clck = clock();

		//printf("REMOVE DETACHED\n");
		for (int i = 0; i < al->size(); i++)
			if(al->at(i)->type() != Agent::substrate)
			if( al->at(i)->_substrateContacts == 0)
			{
				//fprintf( stderr, "cell %i lost contact to substrate: [%f, %f, %f]\n", i, al->at(i)->position()[0], al->at(i)->position()[1], al->at(i)->position()[2]);
				//t=tend;
				box->remove( al->at(i)->position(), al->at(i));
				//fprintf( stderr, "Erase\n");
				al->erase( al->at(i));
				//fprintf( stderr, "Erased\n");
			}

		//if(t>180)
		//	addChemotacticForces( al);


		// BROWNIAN MOTION
		//mexPrintf("BROWNIAN MOTION\n");
		//addBrownianMotion( al, dimy, p_randVelocity, dt);
		//printf("BROWNIAN MOTION\n");
		addBrownianMotion( al, dimy, p_randVelocity, dt);
		clck_rand += clock() - clck; clck = clock();

		dt = getTimeStep( al, dimy, OverDamped ); dt = MIN(dt, 1e-1);if (t + dt > tend) dt = tend - t;
		clck_timestep += clock() - clck; clck = clock();

		// UPDATE POSITIONS
		//printf("UPDATE POSITIONS\n");
		updatePositions(  al, dimy, box, dt, OverDamped);
		clck_updpos += clock() - clck; clck = clock();

		// DIVISION along SUBSTRATE
		//printf("DIVISION\n");
		bool WNTDependence = false;
		//performCellDivision( al, dimy, box, dt, minR, maxR, cycleTime, WNTDependence);
		//if(false)
		performCellDivisionAlongSubstrate( al, dimy, box, dt, minR, maxR, cycleTime, WNTDependence);
		clck_div += clock() - clck; clck = clock();

		// DIFFERENTIATION
		//printf("DIFFERENTIATION\n");
		if(t>180)
		for (int i = 0; i < al->size(); i++)
			if(al->at(i)->type() != Agent::substrate)
				//if( al->at(i)->radius() == minR)
				//al->at(i)->differentiate2( al, box);
				al->at(i)->differentiate( al, box);
				//al->at(i)->differentiateProbabilistic( al, box);


		t += dt;
		time[0] = t;




		if(verbosity)
			mexPrintf("\r time = %.2lf hours, clock (%.2lf, %.2lf, %.2lf, %.2lf, %.2lf),  dt =  %e \b",
				t,
				clck_rand/(float)CLOCKS_PER_SEC,
				clck_forces/(float)CLOCKS_PER_SEC,
				clck_updpos/(float)CLOCKS_PER_SEC,
				clck_div/(float)CLOCKS_PER_SEC,
				clck_timestep/(float)CLOCKS_PER_SEC,
				dt);




		//usleep(10000);

		// OUTPUT
		if( t >= output_offset &&
		    t <= output_offset + output_duration &&
			//ceil(t-dt) == floor(t)
			floor((t-dt)/output_interval) != floor(t/output_interval)
		)
		{
			int stats[box->Y][Agent::TypeSize];
			int count[box->Y];
			for( int pos=0; pos<box->Y; pos++){
				for( int typ=0; typ<Agent::TypeSize; typ++)
					stats[pos][typ] = 0;
				count[pos] = 0;
			}

			for (int i = 0; i < al->size(); i++)
				if(al->at(i)->type() != Agent::substrate)
				{
					int pos = (int) (al->at(i)->position()[1] / box->latticeConstant);
					int typ = (int) (al->at(i)->type());

					stats[pos][typ] ++;
					count[pos] ++;

					if( t>output_offset){
						g_stats[pos][typ] ++;
						g_count[pos] ++;
					}
				}

			FILE *fp;

			if(write_snapshots && floor((t-dt)/output_interval) != floor(t/output_interval)){
				// snapshot
				char filename[1024];
				sprintf( filename, "snapshot_t%i.dat", (int)((t-output_offset)/output_interval));
				fp = fopen( filename, "w+");
				for( int pos=0; pos<box->Y; pos++){
					fprintf( fp, "%i ", pos);
					for( int typ=0; typ<Agent::TypeSize; typ++)
						fprintf( fp, "%e ", (count[pos] ? (float)stats[pos][typ] / (float)count[pos] : 0));
					fprintf( fp, "\n");
				}
				fclose(fp);
			}

			{ // profiles
				fp = fopen( "profiles.dat", "a+");
				for( int pos=0; pos<box->Y; pos++)
				if( count[pos]){
					fprintf( fp, "%e %i ", floor(t), pos);
					for( int typ=0; typ<Agent::TypeSize; typ++)
						fprintf( fp, "%e ", (float)stats[pos][typ] / (float)count[pos]);
					fprintf( fp, "\n");
				}
				fprintf( fp, "\n\n");
				fclose(fp);
			}

			{ // plot cell positions
				fp = fopen( "cells.dat", "w+");
				for( int it=0; it<Agent::TypeSize; it++){
					for( int a=0; a<al->size(); a++)
						if( al->at(a)->type() == it){
						fprintf( fp, "%f %f %f %i \n", al->at(a)->position()[0], al->at(a)->position()[1], al->at(a)->position()[2], al->at(a)->type());
					}
					fprintf( fp, "\n");
				}
				fclose(fp);
			}

			// plot average profiles
			{

				fp = fopen( "profiles_avg.dat", "w+");
				for( int pos=0; pos<box->Y; pos++)
				if( g_count[pos]){
					fprintf( fp, "%i ", pos);
					for( int typ=0; typ<Agent::TypeSize; typ++)
						fprintf( fp, "%e ", (float)g_stats[pos][typ] / (float)g_count[pos]);

					fprintf( fp, "\n");
				}
				fclose(fp);


			}

		}
	}

	double lik = 0;
	for( int pos=0; pos<box->Y; pos++){
		for( int typ=0; typ<Agent::TypeSize; typ++)
		if( data_std[pos][typ] > 0.)
			lik += pow( (g_count[pos] ? (float)g_stats[pos][typ] / (float)g_count[pos] : 0) - data_mean[pos][typ], 2) / pow( data_std[pos][typ], 2);
	}
	//fprintf( stderr, "lik = %e\n", lik);
	fprintf( stdout, "%e\n", lik);


}

void simCylinder( AgentList<> *al, Boxes<Agent*> *box, float *time) {

	setlocale(LC_ALL,"C");

	al->setDimensions(2);

	mexPrintf("ADD BOUNDARY\n");
	float dz = 5;
	float maxx = 0;

	for (float y = 0; y<=1000; y += dz) {

		float x = 0;

		for( float dx = 200; dx>=5; dx /= 2){
			Agent *a;

			a = new Agent(
					x,
					y);
			a->set_type( Agent::substrate);
			a->radius() = 5.;
			box->add(a->position(), a);
			al->insert(a);

			x = x + dx;
		}

		maxx = x;
	}

	// ADD CELLS
	mexPrintf("ADD CELLS\n");
	for (int x = 5; x <= maxx-5; x+=5) {
		//printf("x=%lf\n", b->X*b->latticeConstant/2.);
		Agent *a;
		a = new Agent( x+RAND01, 5+RAND01 );
		a->set_type(Agent::StemCell);
		a->radius() = 5.;
		box->add(a->position(), a);
		al->insert(a);

		a = new Agent( x, 0 );
		a->set_type(Agent::substrate);
		a->radius() = 5.;
		box->add(a->position(), a);
		al->insert(a);

	}



	// SIMULATE
	srand(0);


	mexPrintf("SIMULATE (%i)\n", rand());
	float dt = 0.1;
	float tend = 1500;
	float randVelocity[Agent::TypeSize], *p_randVelocity = randVelocity;
	randVelocity[Agent::StemCell]   = 10;//10;
	//randVelocity[Agent::Goblet]   = 1000;
	//randVelocity[Agent::Absorptive] = 1000;
	//randVelocity[Agent::Secretory] = 1000;
	int dimy = 2;

	float minR = 5;
	float maxR = 7;
	float cycleTime[Agent::TypeSize];
	cycleTime[Agent::StemCell]	 = 24.;
	cycleTime[Agent::Absorptive] = 24.;
	cycleTime[Agent::Secretory]	 = 24.;
	cycleTime[Agent::substrate]	 = 100000000.;

	float **adhesionEnergies = allocMatrix<float>(Agent::TypeSize, Agent::TypeSize);
	initMatrix<float>( adhesionEnergies, Agent::TypeSize, Agent::TypeSize, 0) ;
	adhesionEnergies[Agent::StemCell][Agent::substrate] = adhesionEnergies[Agent::substrate][Agent::StemCell] = 5e1*0;
	adhesionEnergies[Agent::StemCell][Agent::StemCell] = 1e1*0;

	FILE *fp = fopen( "height.dat", "w+");
	int hn = 10000;
	float height[hn];
	float pressure[hn];
	float velocity[hn][2];
	for( int h=0; h<hn; h++){
		height[h] = -1e10;
		pressure[h] = -1e10;
		velocity[h][0] = 0;
		velocity[h][1] = 0;
	}



	//for (float t = 0; t < tend; t += dt)
	float t = 0;
	dt = 1e-3;
	while( t<tend)
	{

		resetAcceleration( al, dimy);

		// External Forces
		addBrownianMotion( al, dimy, p_randVelocity, dt);
		//addBrownianMotion( al, dimy, p_randVelocity);
		addCellCellInteractionForces( al, dimy, box, adhesionEnergies);
		//addSpringForces( al, dimy, springs, countSprings);
		//addCellMediumFriction( al, dimy, dt);


		// Update Positions
		updatePositions( al, dimy, box, dt, OverDamped);

		// Cell Divisions
		//performCellDivision( al, dimy, box, dt, minR, maxR, cycleTime, false);
		performCellDivisionAlongSubstrate( al, dimy, box, dt, minR, maxR, cycleTime, false);

		t += dt;
		time[0] = t;

		// Time Step
		dt = getTimeStep( al, dimy, OverDamped);
		if (t + dt > tend)
			dt = tend - t;
		eprintf("\r time = %.2lf hours,  dt =  %e \b", t, dt);


		//usleep(5000000/t);

		// Statistics

		{
			for( int a=0; a<al->size(); a++){
				if(al->at(a)->position()[2] != 0)
					t = tend;
				if(al->at(a)->type() != Agent::substrate)
				for( int col = (int)((al->at(a)->position()[0] - al->at(a)->radius()/2.) / 7);
						 col<= (int)((al->at(a)->position()[0] + al->at(a)->radius()/2.) / 7);
						 col++)
				{
					height[col] = MAX( height[col], al->at(a)->position()[1]);
					pressure[col] = MAX( pressure[col], al->at(a)->pressure());
					velocity[col][0] += al->at(a)->acceleration()[0] * dt;
					velocity[col][1] += al->at(a)->acceleration()[1] * dt;

				}
			}
		}

		if((int)(t/10) != (int)(t/10-dt/10)){

			for( int h=0; h<hn; h++)
				if( height[h] != -1e10)
				fprintf(fp, "%i %i %e %e %e %e\n",
						(int)(t-dt),
						h,
						height[h],
						pressure[h],
						velocity[h][0],
						velocity[h][1]);
			fprintf(fp, "\n\n");

			for( int h=0; h<hn; h++){
				height[h] = -1e10;
				pressure[h] = -1e10;
				velocity[h][0] = 0;
				velocity[h][1] = 0;
			}

		}

		//usleep(10000);
	}

	fclose(fp);

}

#include <gsl/gsl_odeiv2.h>
#include <gsl/gsl_errno.h>
#include "rungekutta.hpp"


void simMonolayer( AgentList<> *al, Boxes<Agent*> *box, float *time) {

	setlocale(LC_ALL,"C");

	al->setDimensions(2);

	// Define System
	//enum CellTypes{ Normal } cellTypes;

	// ADD CELLS
	mexPrintf("ADD CELLS\n");
	{
		//printf("x=%lf\n", b->X*b->latticeConstant/2.);
		Agent *a;
		a = new Agent( box->X*box->latticeConstant/2, box->Y*box->latticeConstant/2 );
		a->set_type( Agent::StemCell);
		a->radius() = 5.;
		box->add(a->position(), a);
		al->insert(a);

		a = new Agent( box->X*box->latticeConstant/2+5, box->Y*box->latticeConstant/2 );
		a->set_type( Agent::StemCell);
		a->radius() = 5.;
		box->add(a->position(), a);
		al->insert(a);
	}



	// SIMULATE
	srand(0);


	mexPrintf("SIMULATE (%i)\n", rand());
	float randVelocity[Agent::TypeSize], *p_randVelocity = randVelocity;
	randVelocity[Agent::StemCell]   = 10;//10;
	int dimy = 2;

	float minR = 5;
	float maxR = 7;
	float cycleTime[Agent::TypeSize];
	cycleTime[Agent::StemCell]	 = 24.;

	float **adhesionEnergies = allocMatrix<float>(Agent::TypeSize, Agent::TypeSize);
	initMatrix<float>( adhesionEnergies, Agent::TypeSize, Agent::TypeSize, 0) ;
	adhesionEnergies[Agent::StemCell][Agent::StemCell] = 5e1;

	float t = 0;
	float dt = 1e-3;
	double h  = dt;
	float tend = 1500;

	al->pressureThreshold  = 5e2;


	double *y = (double*)malloc(100000*sizeof(double)); for( int i=0; i<100000; i++) y[i]=0;//RAND01*1;
	normrnd( y, 100000, 100, 1);
	int *param[2] = {(int*)al, (int*)box};
	gsl_odeiv2_system sys = {0,0,0,param};

	switch( 2){
	case 1: // DELTA-NOTCH


		box->nbmolecules = 1;
		al->nbmolecules  = 2;
		sys.function = DeltaNotchModel;


		break;

	case 2: // CASPASE 8

		box->nbmolecules = 1;
		al->nbmolecules  = 18;
		sys.function = Caspase8Model;


		break;
	}
	box->molecules = y;
	al->molecules = &y[box->X*box->Y*box->Z*box->nbmolecules];

	int length = 100000;
	double **k   = (double **) malloc(length*sizeof(double*));
	for( int i=0; i<length; i++) k[i] = (double * ) malloc(6*sizeof(double));
	double *dydt = (double * ) malloc(length*sizeof(double));
	double *yi   = (double * ) malloc(length*sizeof(double));

	for( int i=0; i<al->size(); i++){
		double *x = (double *)&y  [box->X*box->Y*box->Z*box->nbmolecules + i*al->nbmolecules];
		// INITIALIZE
		fprintf( stderr, "INITIALIZE\n");
		// CD95-HeLa
		init( x, 18, "0 116 93 155 236 973 5178");

		// HeLa wt
		//init( x, 18, "0 12 90 127 224 1909 3316");
	}

	int apoptosis = 0;
	FILE *fp = fopen( "countApt.dat", "a+");

	while( t<tend)
	{

		resetAcceleration( al, dimy);

		// External Forces
		addBrownianMotion( al, dimy, p_randVelocity, dt);
		//addBrownianMotion( al, dimy, p_randVelocity);
		addCellCellInteractionForces( al, dimy, box, adhesionEnergies);
		//addSpringForces( al, dimy, springs, countSprings);
		//addCellMediumFriction( al, dimy, dt);


		// Update Positions
		updatePositions( al, dimy, box, dt, OverDamped);

		// Cell Divisions
		performCellDivision( al, dimy, box, dt, minR, maxR, cycleTime, false);
		//performCellDivisionAlongSubstrate( al, dimy, box, dt, minR, maxR, cycleTime, false);


		if( (int)t != (int)(t+dt)){
			sys.dimension = box->nbmolecules*box->X*box->Y*box->Z + al->nbmolecules*al->size();
			//fprintf(stderr, "sys.dimension = %i\n", sys.dimension);

			double tcont = t;
			double tmax  = t + 1;
			do{
				solveRungeKutta4<6>( sys, &tcont, tmax, &h, y, 1e-6, rungekuttafehlberg,
						k, dydt, yi, length);
			}while( tcont < tmax);

			fprintf( stdout, "%e ", t);
			for( int m=0; m<al->nbmolecules; m++)
				fprintf( stdout, "%e ", al->molecules[m]);
			fprintf( stdout, "\n");
		}

		t += dt;
		time[0] = t;

		// Time Step
		dt = getTimeStep( al, dimy, OverDamped);
		eprintf("\r time = %.2lf hours,  dt =  %e, h = %e \b", t, dt, h);
		dt = fmin( dt, 0.1);
		dt = fmax( dt, 0.0001);
		if (t + dt > tend)
			dt = tend - t;

		// Cell Death
		if(false)
		for( int i=0; i<al->size(); i++){
			// local conc
			/*int x= box->getIndex(al->at(i)->position()[0]);
			int _y= box->getIndex(al->at(i)->position()[1]);
			int z= box->getIndex(al->at(i)->position()[2]);
			int dx=box->Y*box->Z;
			int dy=box->Z;
			int dz=1;
			int j =x*dx + _y*dy + z*dz;

			double conc = box->molecules[ j ];*/

			//if( RAND01 < conc*dt*0.08 && t>100)

			int idx = box->X*box->Y*box->Z * box->nbmolecules + i * al->nbmolecules + 13;
			fprintf(stderr, "%e\n", y[idx]);
			if( 50 < y[idx] && t>24*7)
			{
				box->remove(al->at(i)->position(), al->at(i));
				al->erase(al->at(i));
				i--;


				apoptosis++;
				fprintf( stderr, "%e %e\n", t, apoptosis);

			}
		}


		//usleep(5000000/t);

		//if( t>24*7)
		//	usleep(100);
	}

	fclose(fp);

	free(y);

}


