/*
 * Kinetics.cpp
 *
 *  Created on: 07.11.2013
 *      Author: jagiella
 */

#include "Agent.hpp"
#include "Spring.hpp"
#include "Boxes.hpp"
#include "Kinetics.hpp"
#include "Macros.hpp"
#include "Covariance.hpp"

#include <gsl/gsl_math.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_statistics_double.h>


void randSphere3( float *vector){
	// Marsaglia (1972)
	float a;
	float b;
	do{
		a = RAND01*2-1;
		b = RAND01*2-1;
	}while(a*a+b*b >= 1);

	vector[0] = 2*a*sqrt(1-a*a-b*b);
	vector[1] = 2*b*sqrt(1-a*a-b*b);
	vector[2] = 1-2*( a*a+b*b);
}
void randSphere3( float *vector, float length){
	// Marsaglia (1972)
	float a;
	float b;
	do{
		a = RAND01*2-1;
		b = RAND01*2-1;
	}while(a*a+b*b >= 1);

	vector[0] = length*2*a*sqrt(1-a*a-b*b);
	vector[1] = length*2*b*sqrt(1-a*a-b*b);
	vector[2] = length*(1-2*( a*a+b*b));
}

void normrnd( float *vector, int dimensions, float mean, float variance)
{
	// Box-Muller
	/*float u1 = rnd();
	float u2 = rnd();

	return mean + variance*sqrt(-2 * log(u1)) * cos( 2 * M_PI* u2);
	 */

	float x1, x2, w;
	for( int d=0; d<dimensions; d+=2){
		// Alternative Box-Muller
		do {
			x1 = 2.0 * RAND01 - 1.0;
			x2 = 2.0 * RAND01 - 1.0;
			w = x1 * x1 + x2 * x2;
		} while ( w >= 1.0 );

		w = sqrt( (-2.0 * log( w ) ) / w );

		vector[d]   = mean + variance * x1 * w;
		if( d+1 < dimensions)
		vector[d+1] = mean + variance * x2 * w;
	}
}
void normrnd( double *vector, int dimensions, double mean, double variance)
{
	// Box-Muller
	/*float u1 = rnd();
	float u2 = rnd();

	return mean + variance*sqrt(-2 * log(u1)) * cos( 2 * M_PI* u2);
	 */

	double x1, x2, w;
	for( int d=0; d<dimensions; d+=2){
		// Alternative Box-Muller
		do {
			x1 = 2.0 * RAND01 - 1.0;
			x2 = 2.0 * RAND01 - 1.0;
			w = x1 * x1 + x2 * x2;
		} while ( w >= 1.0 );

		w = sqrt( (-2.0 * log( w ) ) / w );

		vector[d]   = mean + variance * x1 * w;
		if( d+1 < dimensions)
		vector[d+1] = mean + variance * x2 * w;
	}
}


float getTimeStep( AgentList<Agent> *al, int dimy, int kineticType ){
		// ESTIMATE TIME STEP
		float maxDisplacement = 0.09;
		float dt = 1;

		switch( kineticType)
		{
		case OverDamped:
			for (int i = 0; i < al->size(); i++)
			if (al->at(i)->type() != Agent::substrate)
				for (int d = 0; d < dimy; d++)
					if (al->at(i)->acceleration()[d] != 0.)
						dt = MIN( dt, sqrt(2*maxDisplacement / fabs(al->at(i)->acceleration()[d])));
			break;

		case Exact:

			for (int i = 0; i < al->size(); i++)
			if (al->at(i)->type() != Agent::substrate){
				for (int d = 0; d < dimy; d++) {
					if (al->at(i)->acceleration()[d] != 0.) {
						float sign = (al->at(i)->speed()[d] >= 0 ? 1. : -1.);

						// displacement in direction of velocity
						// if velocity & acceleration have different signes
						float temp1 = -al->at(i)->speed()[d] / al->at(i)->acceleration()[d]- sqrt(	pow(al->at(i)->speed()[d]/ al->at(i)->acceleration()[d],2)+ 2 * maxDisplacement * sign/ al->at(i)->acceleration()[d]);
						// if velocity & acceleration have same signes
						float temp2 =-al->at(i)->speed()[d]	/ al->at(i)->acceleration()[d]+ sqrt(	pow(al->at(i)->speed()[d]/ al->at(i)->acceleration()[d],2)+ 2 * maxDisplacement * sign/ al->at(i)->acceleration()[d]);

						// displacement in opposite direction of velocity
						//float temp3 = - al->at(i)->speed()[d]/al->at(i)->acceleration()[d] - sqrt( pow( al->at(i)->speed()[d]/al->at(i)->acceleration()[d], 2) - 2*maxDisplacement*sign/al->at(i)->acceleration()[d]);
						float temp4 =-al->at(i)->speed()[d]/ al->at(i)->acceleration()[d]+ sqrt(pow(al->at(i)->speed()[d]/ al->at(i)->acceleration()[d],2)- 2 * maxDisplacement * sign/ al->at(i)->acceleration()[d]);

						if (!isnan(temp1)) {
							if (temp1 > 0)
								dt = MIN( dt, temp1);
							else
								dt = MIN( dt, temp2);
						} else
							dt = MIN( dt, temp4);
						if (isnan(dt) || dt < 0) {
							fprintf( stderr, "a=%e, v=%e\n", al->at(i)->acceleration()[d],al->at(i)->speed()[d]);
							exit(0);
						}
					}
				}
			}

			break;

		}

		return dt;
}


void resetAcceleration( AgentList<Agent> *al, int dim){ UNUSED(dim);

	// RESET ACCELERATION
	for (int i = 0; i < al->size(); i++) {
		al->at(i)->set_acceleration(0.);
		//al->at(i)->pressure() = 0;
	}
}


void addBrownianMotion( AgentList<Agent> *al, int dimensions, float *diffusionCoefficient, float dt){

	float lambda = 1;
	float mass = 1;
	float effectiveDiffusionCoefficient;

	float rnd[dimensions];

	// BROWNIAN MOTION
	//mexPrintf("BROWNIAN MOTION\n");
	for (int i = 0; i < al->size(); i++)
	if (al->at(i)->type() != Agent::substrate){

		normrnd( rnd, dimensions, 0, 1);
		effectiveDiffusionCoefficient = lambda/mass*sqrt(2*diffusionCoefficient[(int)al->at(i)->type()]/dt);
		for( int d=0; d<dimensions; d++)
			al->at(i)->acceleration()[d] += - lambda/mass*al->at(i)->speed()[d] + effectiveDiffusionCoefficient*rnd[d];
			//al->at(i)->acceleration()[d] += lambda/mass *
			//	( - al->at(i)->speed()[d] + RANDNORM( 0, 1) * sqrt(2*diffusionCoefficient[al->at(i)->type()]/dt));
	}

}

void addBrownianMotion( AgentList<Agent> *al, int dimy, float *randVelocity){ UNUSED(dimy);

	// BROWNIAN MOTION
	//mexPrintf("BROWNIAN MOTION\n");
	for (int i = 0; i < al->size(); i++)
	if (al->at(i)->type() != Agent::substrate){
		//mexPrintf("cell %i: set speed\n", i);
		randSphere3( al->at(i)->acceleration(), randVelocity[(int)al->at(i)->type()]);
			//for (int d = 0; d < dimy; d++)
			//	al->at(i)->acceleration()[d] = randVelocity * (RAND01 - 0.5);

	}

}

void addSpringForces( AgentList<Agent> *al, int dimy, Spring<int>** springs, int countSprings){
	// SPRING FORCES

	for (int s = 0; s < countSprings; s++) {
		int i = springs[s]->get(0);
		int j = springs[s]->get(1);

		float dist0 = springs[s]->distance();
		float dist1 = 0;
		for (int d = 0; d < dimy; d++)
			dist1 += pow(
					al->at(j)->position()[d] - al->at(i)->position()[d], 2);
		dist1 = sqrt(dist1);

		// acceleration
		float acceleration[3];
		//float springForce[3];
		//for( int d=0; d<dimy; d++)
		//	 springForce[d] = (al->at(j)->position()[d]-al->at(i)->position()[d]) * (dist1-dist0);

		al->Acceleration(acceleration, i, j, -(dist1 - dist0) * 100000.);

		for (int d = 0; d < dimy; d++) {
			al->at(i)->acceleration()[d] += acceleration[d];
			al->at(j)->acceleration()[d] -= acceleration[d];
		}

		//fprintf( stderr, "(%i, %i) d = %lf -> a (%lf, %lf)\n", i, j, dist1, acceleration[0], acceleration[1]);

		/*if (dist1 > al->at(i)->radius() + dist1 > al->at(j)->radius() && t > 10) {
			countSprings--;
			springs[s] = springs[countSprings];
			s--;
		}*/
	}

}

void addCellCellInteractionForces( AgentList<Agent> *al, int dimy, Boxes<Agent*> *box, float **adhesionEnergies){

	// CELL-CELL-INTERACTION
	//mexPrintf("CELL-CELL-INTERACTION\n");
	BoxIterator<Agent*> boxIt( box, 1);

	for (int i = 0; i < al->size(); i++){

		float tmp_pressure = 0;
		al->at(i)->_substrateContacts = 0;
		if(al->at(i)->type() != Agent::substrate)
		{
			boxIt.init( al->at(i)->position());
			for( ; !boxIt.end(); boxIt++){
				Agent *a = *boxIt;
				int j = a->index();
				if (j != i && al->inContact(i,j)){
					al->ExtendedHertzForce( i, j, adhesionEnergies, al->at(i)->acceleration(), &tmp_pressure);
					if( al->at(j)->type() == Agent::substrate)
						al->at(i)->_substrateContacts++;
				}
			}
		}

		al->at(i)->pressure() = tmp_pressure;
	}
}

void addChemotacticForces( AgentList<Agent> *al){
	for (int i = 0; i < al->size(); i++)
		if(al->at(i)->type() != Agent::substrate)
			al->AccelerationChemotaxis( al->at(i)->acceleration(), i);
}

void updatePressure( AgentList<Agent> *al, Boxes<Agent*> *box, float **adhesionEnergies){ UNUSED(adhesionEnergies);

	BoxIterator<Agent*> boxIt( box, 1);

	for (int i = 0; i < al->size(); i++){

		if(al->at(i)->type() != Agent::substrate){
			float tmp_pressure = 0;
			boxIt.init( al->at(i)->position());
			for( ; !boxIt.end(); boxIt++){
				int j = (*boxIt)->index();
				if (j != i && al->inContact(i,j))
					tmp_pressure += al->Pressure(i, j);
			}
			al->at(i)->pressure() = tmp_pressure;
		}
	}
}

void addCellMediumFriction( AgentList<Agent> *al, int dimy, float dt){
	// CELL-MEDIUM-FRICTION
	//mexPrintf("CELL-MEDIUM-FRICTION\n");

	for (int i = 0; i < al->size(); i++)
		if (al->at(i)->type() != Agent::substrate){
		//mexPrintf("cell %i: set speed\n", i);
		for (int d = 0; d < dimy; d++)
			al->at(i)->speed()[d] *= (1 - dt * 0);
	}
}

void updatePositions(  AgentList<> *al, int dimy, Boxes<Agent*> *box, float dt, int kineticType){
	// UPDATE POSITIONS
	for (int i = 0; i < al->size(); i++) {
		if (al->at(i)->type() != Agent::substrate){

			box->remove(al->at(i)->position(), al->at(i));

			for (int d = 0; d < dimy; d++){

				/*if( fabs(dt*al->at(i)->speed()[d] + 0.5*dt*dt*al->at(i)->acceleration()[d]) > maxDisplacement){
				 mexPrintf("too large displacement =  %e (>%e)\n", fabs(dt*al->at(i)->speed()[d] + 0.5*dt*dt*al->at(i)->acceleration()[d]), maxDisplacement);
				 }*/
				switch( kineticType){
				case OverDamped:
//#ifdef OVERDAMPED

				al->at(i)->position()[d] +=
						+ 0.5 * dt * dt * al->at(i)->acceleration()[d];


				al->at(i)->speed()[d] = 0;
				break;
//#else
				case Exact:
				// update positions
				//box->remove(al->at(i)->position(), i);
				al->at(i)->position()[d] += dt * al->at(i)->speed()[d]
						+ 0.5 * dt * dt * al->at(i)->acceleration()[d];
				//box->add(al->at(i)->position(), i);

				// update speed
				al->at(i)->speed()[d] += dt * al->at(i)->acceleration()[d];
				break;
//#endif
				}
			}



			// Remove Cells outside Domain
			if( !box->inside(al->at(i)->position())){
				fprintf(stderr, "ERASE!\n");
				al->erase(al->at(i));
				i--;
//				box->add(al->at(i)->position(), al->at(i));
			}else
				box->add(al->at(i)->position(), al->at(i));

		}
	}


}


void performCellDivisionAlongSubstrate( AgentList<Agent> *al, int dimy, Boxes<Agent*> *box, float dt, float minR, float maxR, float *cycleTime, bool WNTDependence)
{
	// DIVISION along SUBSTRATE
	//fprintf(stderr, "DIVISION along SUBSTRATE\n");
	for (int i = 0; i < al->size(); i++) {
		float beta = 0.4;

		//printf( "pressure = %e\n", al->at(i)->pressure());
		//if (al->at(i)->type() != Agent::substrate)
		//if (al->at(i)->type() == Agent::StemCell || al->at(i)->type() == Agent::Absorptive || al->at(i)->type() == Agent::Secretory )
		if ( cycleTime[(int)al->at(i)->type()] > 0. )
		{
			// GROW RADIUS
			if (al->at(i)->pressure() < al->pressureThreshold){ //4e2//8e2
				if( WNTDependence)
					al->at(i)->radius() += dt * (maxR - minR) / cycleTime[(int)al->at(i)->type()] * al->at(i)->wnt();
				else
					al->at(i)->radius() += dt * (maxR - minR) / cycleTime[(int)al->at(i)->type()];
			}

			//BoxIterator *boxIt = new BoxIterator( box, al->at(i)->position(), 1);
			//BoxIterator boxIt( box, al->at(i)->position(), 0);
			//fprintf(stderr,"WK(%i) = %e\n", i, al->CompressionEnergy( i, boxIt));
			//delete boxIt;
			//printf("RADIUS=%lf\n", al->at(i)->radius() );
			//mexPrintf("RADIUS=%lf\n", al->at(i)->radius() );
			if (al->at(i)->radius() >= maxR) {

				/// BOXES

				// mean substrate
				double **A = allocMatrix<double>(dimy, 100);
				float substrate[dimy];
				for (int d = 0; d < dimy; d++)
					substrate[d] = 0.;
				int nsubstrate = 0;

				//fprintf(stderr, "count neighbors in contact\n");
				BoxIterator<Agent*> boxIt( box, al->at(i)->position(), 1);
				for( ; !boxIt.end(); boxIt++){
					//count++;
					Agent *a = *boxIt;
					int j = a->index();
					if (j != i && al->at(j)->type() == Agent::substrate && al->inContact( i, j)) {
						/// BOXES
						for (int d = 0; d < dimy; d++) {
							substrate[d] += al->at(j)->position()[d];
							A[d][nsubstrate] = al->at(j)->position()[d];
						}
						nsubstrate += 1;
					}
				}

				for (int d = 0; d < dimy; d++)
					substrate[d] /= nsubstrate;

				//fprintf(stderr, "finished counting neighbors in contact\n");

				/*mexPrintf("cell position = (%e, %e)\n",
						al->at(i)->position()[0], al->at(i)->position()[1]);
				mexPrintf("center of mass = (%e, %e)\n", substrate[0],
						substrate[1]);

				mexPrintf("DIVISION\n");*/

				Agent *newa;

				if (nsubstrate == 0)
				{
					//fprintf(stderr, "Cell Division (Random)\n");

					//mexPrintf("CELL DETACHED!\n");
					//t = tend;
					switch(dimy){
					case 2:
					newa = new Agent(
							al->at(i)->position()[0]+0.1*RAND01,
							al->at(i)->position()[1]+0.1*RAND01);
					break;
					case 3:
					newa = new Agent(
							al->at(i)->position()[0]+0.1*RAND01,
							al->at(i)->position()[1]+0.1*RAND01,
							al->at(i)->position()[2]+0.1*RAND01);
					break;
					}
					/*this->terminate();
					this->quit();
					this->exit();
					return;*/
				}

				if( nsubstrate>1){
					//fprintf(stderr, "Cell Division (Eigen Vectors)\n");

					// EIGENVECTORS

					//printMatrix<double>(A, dimy, nsubstrate);

					// Allocate Memory
					gsl_matrix *covA = gsl_matrix_alloc(dimy, dimy);
					gsl_matrix *evec = gsl_matrix_alloc(dimy, dimy);
					gsl_vector *eval = gsl_vector_alloc(dimy);

					// Calculate Covariance Matrix
					gsl_covariance_matrix(A, dimy, nsubstrate, covA->data);
					//printMatrix<double>(covA->data, dimy, dimy);

					// Calculate Eigen Values & Vectors
					gsl_eigen_symmv_workspace * w = gsl_eigen_symmv_alloc(dimy);
					gsl_eigen_symmv(covA, eval, evec, w);
					gsl_eigen_symmv_free(w);

					// Sort for highest Eigen Values
					gsl_eigen_symmv_sort(eval, evec, GSL_EIGEN_SORT_ABS_DESC);
					//printMatrix<double>(evec->data, dimy, dimy);

					// EIGENVECTORS

					switch(dimy){
					case 3:
					newa = new Agent(
							(al->at(i)->position()[0] + evec->data[0*dimy+0]*minR*0.01)*(1-beta) + beta*substrate[0],
							(al->at(i)->position()[1] + evec->data[1*dimy+0]*minR*0.01)*(1-beta) + beta*substrate[1],
							(al->at(i)->position()[2] + evec->data[2*dimy+0]*minR*0.01)*(1-beta) + beta*substrate[2]);
					break;

					case 2:
					newa = new Agent(
							(al->at(i)->position()[0] + evec->data[0*dimy+0]*minR*0.01)*(1-beta) + beta*substrate[0],
							(al->at(i)->position()[1] + evec->data[1*dimy+0]*minR*0.01)*(1-beta) + beta*substrate[1]);
					break;

					}

					// Free Memory
					gsl_matrix_free(covA);
					gsl_matrix_free(evec);
					gsl_vector_free(eval);


					box->remove( al->at(i)->position(), al->at(i));
					for( int d=0; d<dimy; d++)
						al->at(i)->position()[d] = al->at(i)->position()[d]*(1-beta) + beta*substrate[d];
					box->add( al->at(i)->position(), al->at(i));

					//fprintf(stderr, "finished Cell Division (Eigen Vectors)\n");
				}

				if( nsubstrate==1){
					//fprintf(stderr, "Cell Division (Normal), D=%i\n", dimy);

					// normal
					float normal[dimy];
					float abs_normal = 0.;

					switch(dimy){
					case 2:

						for (int d = 0; d < dimy; d++){
							normal[d] = (d==0?+1:-1) * (substrate[d] - al->at(i)->position()[d]);
							abs_normal += normal[d]*normal[d];
						}
						if(abs_normal < 1e-6)
						{
							fprintf(stderr, "\nNormal Vector %i is singular\n", i);
							exit(0);
							return;
						}
						abs_normal = sqrt(abs_normal);
						for (int d = 0; d < dimy; d++)
							normal[d] /= abs_normal;
						break;

					case 3:
						float temp[] = {RAND01, RAND01, RAND01};
						// cross-product
						normal[0] = substrate[1]*temp[2] - substrate[2]*temp[1];
						normal[1] = substrate[2]*temp[0] - substrate[0]*temp[2];
						normal[2] = substrate[0]*temp[1] - substrate[1]*temp[0];
						for (int d = 0; d < dimy; d++){
							abs_normal += normal[d]*normal[d];
						}
						abs_normal = sqrt(abs_normal);
						for (int d = 0; d < dimy; d++)
							normal[d] /= abs_normal;
						break;
					}

					//mexPrintf("normal vector => (%e, %e) abs = %e\n", normal[0], normal[1], abs_normal);

					newa = new Agent(
							(al->at(i)->position()[0] + normal[1]*minR*0.01)*(1-beta) + beta*substrate[0],
							(al->at(i)->position()[1] + normal[0]*minR*0.01)*(1-beta) + beta*substrate[1],
							al->at(i)->position()[2]);

					box->remove( al->at(i)->position(), al->at(i));
					for( int d=0; d<dimy; d++)
						al->at(i)->position()[d] = al->at(i)->position()[d]*(1-beta) + beta*substrate[d];
					box->add( al->at(i)->position(), al->at(i));
				}



				/*mexPrintf("(%e, %e) -- (%e, %e)\n",
						al->at(i)->position()[0], al->at(i)->position()[1],
						newa->position()[0], newa->position()[1]);*/
				newa->set_type( al->at(i)->type());

				//printf("Adding cell\n");
				if( box->inside( newa->position())){
					box->add(newa->position(), newa);
					al->insert(newa);
					newa->radius() = minR;
				}else
					delete newa;
				al->at(i)->radius() = minR;
				//printf("finished adding cell\n");

				/*if(RAND01 < 0.5)
					newa->differentiate( al, box);
				else
					al->at(i)->differentiate( al, box);*/
				// ATTENTION!!!!!!!!!
				//t=tend;

				freeMatrix(dimy, A);

			}
		}

		//al->at(i)->differentiate( al, box);

	}
}


void performCellDivision( AgentList<Agent> *al, int dimy, Boxes<Agent*> *box, float dt, float minR, float maxR, float *cycleTime, bool WNTDependence)
{
	int delta = 0;
	int notch = 1;

	// DIVISION along SUBSTRATE
	//fprintf(stderr, "DIVISION along SUBSTRATE\n");
	for (int i = 0; i < al->size(); i++) {
		//float beta = 0.4;

		//if (al->at(i)->type() == Agent::StemCell || al->at(i)->type() == Agent::Absorptive || al->at(i)->type() == Agent::Secretory )
		if ( cycleTime[(int)al->at(i)->type()] > 0. )
		{

			// GROW RADIUS
			if (al->at(i)->pressure() < al->pressureThreshold){ //4e2//8e2
				//fprintf(stderr, "GROW RADIUS\n");
				if( WNTDependence)
					al->at(i)->radius() += dt * (maxR - minR) / cycleTime[(int)al->at(i)->type()] * al->at(i)->wnt();
				else
					al->at(i)->radius() += dt * (maxR - minR) / cycleTime[(int)al->at(i)->type()];
			}
			//else
				//fprintf(stderr, "DONT GROW RADIUS: pressure %e > %e\n", al->at(i)->pressure(), al->pressureThreshold);

			// DIVIDE
			if (al->at(i)->radius() >= maxR) {
				Agent *newa;

				// Position of new Agent
				switch(dimy){
				case 2:
				newa = new Agent(
						al->at(i)->position()[0]+0.1*(2*RAND01-1.),
						al->at(i)->position()[1]+0.1*(2*RAND01-1.));
				break;
				case 3:
				newa = new Agent(
						al->at(i)->position()[0]+0.1*(2*RAND01-1.),
						al->at(i)->position()[1]+0.1*(2*RAND01-1.),
						al->at(i)->position()[2]+0.1*(2*RAND01-1.));
				break;
				}

				// Check Boundary Constraint
				if( box->inside(newa->position()) ){
					newa->set_type( al->at(i)->type());
					box->add(newa->position(), newa);
					al->insert(newa);
					al->at(i)->radius() = minR;
					newa->radius() = minR;

					for( int d=0; d<al->nbmolecules; d++)
						al->molecules[al->nbmolecules*newa->index()+d] = al->molecules[al->nbmolecules*i+d];
				}
				else
					delete newa;
			}
		}
	}
}
