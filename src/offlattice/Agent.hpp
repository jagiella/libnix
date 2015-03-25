/*
 * Agent.hpp
 *
 *  Created on: 12.06.2013
 *      Author: jagiella
 */

#ifndef AGENT_HPP_
#define AGENT_HPP_

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "Boxes.hpp"
#include "Macros.hpp"

//#define PI 3.142

template <class AgentType>
class AgentList;

class Agent {
	template <class AgentType> friend class AgentList;
protected:
//	static const float TPwnt;
//	static const float TPnotch;
//	static const  float TDwnt;

	int   _index;
	float _radius;
	char  _type;
	float *_position;
	float *_speed;
	float *_acceleration;
	float _pressure;
public:
	static float TPwnt;
	static float TPnotch;
	static float TDwnt;

	int   _substrateContacts;

	static float wntmin;
	static float wntmax;

	static char NotchKnockout;
	static char WNTKnockout;
	static char EphrinB12Knockout;
	static char EphB2Knockout;
	static char EphB3Knockout;
	enum types{ StemCell,
		        Paneth, Goblet, Enteroendocrine, Enterocyte,
		        Absorptive, Secretory,
		        substrate,
		        TypeSize};

	Agent( float x=0, float y=0, float z=0): _pressure(0) {
		_radius = 3.;
		_position 		= (float*) malloc( sizeof(float)*3);
		_speed 			= (float*) malloc( sizeof(float)*3);
		_acceleration 	= (float*) malloc( sizeof(float)*3);
		set_position( x, y, z);
		set_speed( 0);
		set_acceleration( 0);
	};
	~Agent(){
		free(_position);
		free(_speed);
		free(_acceleration);
	};

	int index() {return _index;}
	void set_index( int index) {_index = index;}
	void set_position    ( float x, float y=0, float z=0) { _position[0]=x;     _position[1]=y;     _position[2]=z; };
	void set_speed       ( float x, float y=0, float z=0) { _speed[0]=x;        _speed[1]=y;        _speed[2]=z; };
	void set_acceleration( float x, float y=0, float z=0) { _acceleration[0]=x; _acceleration[1]=y; _acceleration[2]=z; };
	//void set_radius( float r) { _radius = r;};
	void set_type( char t) {
		_type = t;
		switch(_type){
		case Agent::StemCell:
			break;
		case Agent::substrate:
			_radius = 1;
			break;
		}
	};
	template <class AgentType>
	void differentiate( AgentList<AgentType> *al, Boxes<AgentType*> *boxes)
	{
		//fprintf(stderr, "BEGIN\n");
		switch(_type){
		case Agent::StemCell:
			if( wnt() <= TPwnt){
				// Priming
				if( notch( al, boxes) < TPnotch){
					_type = Agent::Secretory;
				}else{
					_type = Agent::Absorptive;
				}
			}else{
				// Terminal Paneth
				if( notch( al, boxes) < TPnotch)
					_type = Agent::Paneth;
			}
			break;

		case Agent::Secretory:
			// Terminal Goblet?
			if( wnt() <= TDwnt){
				if( RAND01 < 0.5)
					_type = Agent::Goblet;
				else
					_type = Agent::Enteroendocrine;
			}

			// Reverse priming
			else if( wnt() > TPwnt)
				_type = Agent::StemCell;

			// Delta-Notch-Signaling
			else if( TPwnt > wnt() && wnt() >= TDwnt && notch( al, boxes) >= TPnotch)
				_type = Agent::Absorptive;
			break;

		case Agent::Absorptive:
			// Reverse priming
			if( wnt() > TPwnt)
				_type = Agent::StemCell;

			// Terminal Enterocyte?
			else if( wnt() <= TDwnt)
				_type = Agent::Enterocyte;

			// Delta-Notch-Signaling
			else if( notch( al, boxes) < TPnotch)
				_type = Agent::Secretory;
			break;

			default:
				break;
		}
		//fprintf(stderr, "END\n");

	}
	template <class AgentType>
	void differentiate4( AgentList<AgentType> *al, Boxes<AgentType*> *boxes)
	{
		//fprintf(stderr, "BEGIN\n");
		switch(_type){
		case Agent::StemCell:
			if( notch( al, boxes) <= TPnotch)
				// Asymmetric Division -> Secretory
				_type = Agent::Secretory;
			if( notch( al, boxes) >= TPnotch && wnt() <= TPwnt)
				// Asymmetric Division -> Absorptive
				_type = Agent::Absorptive;
			break;

		case Agent::Secretory:
			// Terminal Goblet?
			if( wnt() <= TDwnt){
				if( RAND01 < 0.85){
					_type = Agent::Goblet;
				}
				else
					_type = Agent::Enteroendocrine;
			}else
			// Terminal Paneth?
				if( wnt() >=  TPwnt && radius() == 5){
					_type = Agent::Paneth;
				}else{
					_type = Agent::StemCell;
					differentiate(al, boxes);
				}
			break;

		case Agent::Absorptive:

			//
			if( wnt() <= TDwnt){
				_type = Agent::Enterocyte;
			}else{
				_type = Agent::StemCell;
				differentiate(al, boxes);
			}
			break;

			default:
				break;
		}
		//fprintf(stderr, "END\n");

	}
	template <class AgentType>
	void differentiate2( AgentList<AgentType> *al, Boxes<AgentType*> *boxes)
	{
		//fprintf(stderr, "BEGIN\n");
		switch(_type){
		case Agent::StemCell:{
			float rnd = RAND01;
			if( rnd <= 0.5)
				// Asymmetric Division -> Secretory
				_type = Agent::Secretory;
			if( rnd >= 0.5 && wnt() <= TPwnt)
				// Asymmetric Division -> Absorptive
				_type = Agent::Absorptive;}
			break;

		case Agent::Secretory:
			// Terminal Goblet?
			if( wnt() <= TDwnt){
				if( RAND01 < 0.5)
					_type = Agent::Goblet;
				else
					_type = Agent::Enteroendocrine;
			}else
			// Terminal Paneth?
				if( wnt() >=  TPwnt && radius() == 5){
					_type = Agent::Paneth;
				}else{
					_type = Agent::StemCell;
					differentiate(al, boxes);
				}
			break;

		case Agent::Absorptive:

			//
			if( wnt() <= TDwnt){
				_type = Agent::Enterocyte;
			}else{
				_type = Agent::StemCell;
				differentiate(al, boxes);
			}
			break;

			default:
				break;
		}
		//fprintf(stderr, "END\n");

	}
	template <class AgentType>
	void differentiate3( AgentList<AgentType> *al, Boxes<AgentType*> *boxes)
	{
		//fprintf(stderr, "BEGIN\n");
		switch(_type){
		case Agent::StemCell:{
			float rnd = RAND01;
			if( rnd <= 0.5)
				// Asymmetric Division -> Secretory
				_type = Agent::Secretory;
			if( rnd >= 0.5 && wnt() <= TPwnt)
				// Asymmetric Division -> Absorptive
				_type = Agent::Absorptive;}
			break;

		case Agent::Secretory:
			// Terminal Goblet?
			if( wnt() <= TDwnt){
				if( RAND01 < 0.85)
					_type = Agent::Goblet;
				else
					_type = Agent::Enteroendocrine;
			}else
			// Terminal Paneth?
				if( wnt() >=  TPwnt && radius() == 5){
					_type = Agent::Paneth;
				}else{
					_type = Agent::StemCell;
					differentiate(al, boxes);
				}
			break;

		case Agent::Absorptive:

			//
			if( wnt() <= TDwnt){
				_type = Agent::Enterocyte;
			}else{
				_type = Agent::StemCell;
				differentiate(al, boxes);
			}
			break;

			default:
				break;
		}
		//fprintf(stderr, "END\n");

	}

	template <class AgentType>
	void differentiateProbabilistic( AgentList<AgentType> *al, Boxes<AgentType*> *boxes)
	{
		if( radius() == 5){ // Differentiation after division

			switch(_type){
			case Agent::StemCell:
				if( wnt() <= TPwnt ){ // Priming
					if( RAND01 < 0.5)
						_type = Agent::Absorptive;
					else
						_type = Agent::Secretory;
				}
			break;

			case Agent::Absorptive:
				if( wnt() <= TDwnt){
					_type = Agent::Enterocyte;
				}
				break;

			case Agent::Secretory:
				if( wnt() <= TDwnt){
					_type = Agent::Goblet;
				}else
				// Terminal Paneth?
				if( wnt() >=  TPwnt && radius() == 5){
					_type = Agent::Paneth;
				}else{
					_type = Agent::StemCell;
					differentiate(al, boxes);
				}
				break;
			}

		}
	}

	float &pressure()		{ return _pressure; }
	float &radius()      	{ return _radius; }
	float* &(position())    	{ return _position;};
	float* (speed())       	{ return _speed;};
	float*& (acceleration())	{ return _acceleration;};
	char type()        		{ return _type;};
	float wnt()				{ /*fprintf(stderr, "[WNT]=%lf\n", -_position[1]+7*1000.);*/
		if(WNTKnockout == -1)
			return 0;
		else
		if(WNTKnockout == +1)
			return 1;
		else{
			//fprintf(stderr, "[WNT]=%lf\n", (_position[1] - Agent::wntmin)/(Agent::wntmax-Agent::wntmin));
			return (_position[1] - Agent::wntmin)/(Agent::wntmax-Agent::wntmin);
		}
			//return 1.-(7*1000.-_position[1])/150.;

	}
	static float wnt( float height){
		return (height - Agent::wntmin)/(Agent::wntmax-Agent::wntmin);
	}
	float foxa2(){
		switch(_type ){
		case Agent::Secretory:
		case Agent::Enteroendocrine:
		case Agent::Goblet:
			return 1;

		case Agent::Paneth:
		case Agent::StemCell:
			return 0.5;

		case Agent::Enterocyte:
			return 0;
		}

		return 0.5/(Agent::TPwnt-Agent::TDwnt)*(wnt(_position[1])-Agent::TDwnt);
	}
	//static float wnt( float *&pos)	{ return 1.-(7*1000.-pos[1])/150.; }
	static float gradEphrinB( float *pos, int d)	{
		UNUSED(pos);

		//float dx = 0.1;
		if(d==1 && !EphrinB12Knockout)
			return 1./(Agent::wntmax-Agent::wntmin);
			//return (pos[1]+dx)/150./dx;
		else
			return 0;
	}
	template <class AgentType>
	float notch( AgentList<AgentType> *al, Boxes<Agent*> *boxes)			{

		if( Agent::NotchKnockout == +1)
			return 10;

		if( Agent::NotchKnockout == -1)
			return 0;

		int i=this->_index;
		float notch=0;
		BoxIterator<Agent*> boxIt( boxes, al->at(i)->position(), 1);
		for( ; !boxIt.end(); boxIt++){
			//count++;
			Agent *a = *boxIt;
			int j = a->index();
			if( j!=i && al->inContact( i,j) && (
					al->at(j)->type() == Agent::Secretory ||
					al->at(j)->type() == Agent::Paneth ||
					al->at(j)->type() == Agent::Goblet ||
					al->at(j)->type() == Agent::Enteroendocrine)){
				notch++;
				//fprintf(stderr, "\nPANETH CELL FOUND!!!\n");
			}
		}

		return notch;
	}
};

template <class AgentType = Agent>
class AgentList {
private:
	AgentType **_agents;
	int     	_size;
	int     	_sizeLimit;
	int	   		_sizeIncrement;
	int 		_dimensions;
public:
	//double y[10000];
	double *molecules;
	int   nbmolecules;

	AgentList( int size = 0, int dimensions = 3) : _agents(0), _size(0), _sizeLimit(size), _sizeIncrement(10), _dimensions(dimensions), nbmolecules(0), molecules(0)
	{
		if( _sizeLimit)
			_agents = (AgentType**) malloc( _sizeLimit * sizeof(AgentType*));
	}
	~AgentList()
	{

		free( _agents);
	}

	static  float pressureThreshold;

	inline void setDimensions( int d){
		_dimensions = d;
	}

	// ELASTIC PROPERTIES

	inline float getElasticModulus( int i){
		// Young/Elastic Modulus (Stiffness): 0...inf Pa
		switch(_agents[i]->_type){
		case Agent::StemCell: // CELL (300...1000)
			return 300;
		case Agent::substrate: // BARRIER
			return 480;
			//return 48000;
		default:
			return 300;
		}
		return 0;
	}
	inline float getPoissionRatio( int i){
		// Poisson Ratio (compressibility): 0=compressible ... 0.5=incompressible
		switch(_agents[i]->_type){
		case Agent::StemCell: // CELL
			return 0.33;
		case Agent::substrate: // BARRIER
			return 0.33;
		default:
			return 0.33;
		}
		return 0;
	}
    inline float BulkModulus( int i){
    	float E = getElasticModulus(i);
    	float eta = getPoissionRatio(i);

    	return E / (3*(1-2*eta));
    }



    // AD-/COHESION PROPERTIES

	inline float getAdhesion( int i, int j){
		// Poisson Ratio (compressibility): 0=incompressible ... 0.5=compressible

//		float sc = 3e2;
//		float sc = 1e2; 2D
		float sc = 5e1*0;
		//float sp = sc;
		float cc = 1e1*0;

		/*float A[Agent::TypeSize][Agent::TypeSize];
		for(int i=0;i<Agent::TypeSize;i++) for(int j=0;j<Agent::TypeSize;j++) A[i][j] = 0;

		A[Agent::StemCell][Agent::StemCell] = 1e1;
		A[Agent::StemCell][Agent::Paneth] = A[Agent::Paneth][Agent::StemCell] = 1e2;
		A[Agent::StemCell][Agent::substrate] = A[Agent::substrate][Agent::StemCell] =  sc;
		A[Agent::Tranient][Agent::substrate] = A[Agent::substrate][Agent::Tranient] =  sc;
		A[Agent::Paneth]  [Agent::substrate] = A[Agent::substrate][Agent::Paneth]   =  sc*2;


		return A[_agents[i]->_type][_agents[j]->_type];*/

		/*switch( _agents[i]->_type){

		case Agent::substrate: switch(_agents[j]->_type){
			case Agent::Paneth:	return sp;
			default:			return sc;
			}break;

		case Agent::Paneth: switch(_agents[j]->_type){
			case Agent::substrate:	return sp;
			default:				return cc;
			}break;

		default: switch( _agents[j]->_type){
			case Agent::substrate:	return sc;
			default:				return cc;
			}
		}*/
		if(_agents[i]->_type == Agent::StemCell  && _agents[j]->_type == Agent::StemCell)
			return cc;
		if(_agents[i]->_type == Agent::substrate  || _agents[j]->_type == Agent::substrate)
			return sc;
		return 0;
	}



	bool inContact(  int i, int j)
	{
       	float sqr_distance =  SqrDistance(i,j);
       	float sum_of_radii = _agents[i]->_radius + _agents[j]->_radius;

		return sum_of_radii*sum_of_radii >= sqr_distance;
		//return sum_of_radii >= Distance(i,j);
	}

	int inContact(  int i, char type, Boxes<AgentType*> *boxes)
	{
		int count=0;
		BoxIterator<AgentType*> boxIt( boxes, at(i)->position(), 1);
		for( ; !boxIt.end(); boxIt++){
			Agent *a = *boxIt;
			int j = a->index();
			if( j!=i && inContact( i,j) && this->at(j)->type() == type ){
				count++;
			}
		}

		return count;
	}

    inline void ExtendedHertzForce(
    		int i, int j, float **adhesionEnergies,
    		float *&F, float *pressure = 0){

    	float Fabs;
       	float distance = Distance(i,j);


       	ExtendedHertzForce( i, j, adhesionEnergies, Fabs, pressure);
       	for( int d=0; d<_dimensions; d++)
       		F[d] += - Fabs * (_agents[j]->_position[d]-_agents[i]->_position[d]) / distance;

    }

    inline void ExtendedHertzForce( int i, int j, float **adhesionEnergies, float &F, float *pressure=0){

    	// Derjaguin-Muller-Toporov (DMT) model of elastic contact

    	// Young/Elastic Modulus (Stiffness)
    	// default = 480 Pa
    	float E1 = getElasticModulus(i),
    		  E2 = getElasticModulus(j);

    	// Poisson Ratio (compressibility): 0=incompressible ... 0.5=compressible
    	// default = 0.5
    	float eta1 = getPoissionRatio(i),
    		  eta2 = getPoissionRatio(j);

    	float Estar = 1. / ((1. - eta1*eta1)/E1 + (1. - eta2*eta2)/E2);
    	float Rstar = 1. / ((1./_agents[i]->_radius)   + (1./_agents[j]->_radius)  );
       	float distance = Distance(i ,j);
    	float overlapp = _agents[i]->_radius + _agents[j]->_radius - distance;

    	if(overlapp<=0){
       		F = 0;
       		if( pressure)
       		*pressure += 0;
    	}else{

    		// REPULSION (Hertz)
    		float Frep = 4./3. *Estar * sqrt(Rstar) * pow(overlapp, 3./2.);

    		// ADHESION  (Extension: Derjaguin-Muller-Toporov)
    		float Fadh = -2*M_PI*Rstar*adhesionEnergies[(int)_agents[i]->type()][(int)_agents[j]->type()];

    		if(isnan(Frep)){
    			fprintf(stderr, "Frep=%e\n", Frep);
    			fprintf(stderr, "Estar=%e\n", Estar);
    			fprintf(stderr, "Rstar=%e (R(%i)=%e, R(%i)=%e)\n", Rstar, i, _agents[i]->_radius, j, _agents[j]->_radius);
    			fprintf(stderr, "overlapp=%e\n", overlapp);
    			fprintf(stderr, "getAdhesion(i, j)=%e\n", getAdhesion(i, j));
    		}

    		// TOTAL FORCE
    		F = Frep + Fadh;

    		// CONTACT PRESSURE
    		// http://en.wikipedia.org/wiki/Contact_mechanics#Contact_between_two_spheres
    		if( pressure)
    		*pressure += pow( 6.*Frep, 1./3.) * pow( Estar/Rstar, 2./3.) / M_PI;
    	}

    }

    inline float Pressure( int i, int j)
    {

    	float E1 = getElasticModulus(i),
    		  E2 = getElasticModulus(j);
    	float eta1 = getPoissionRatio(i),
    		  eta2 = getPoissionRatio(j);

    	float Estar = 1. / ((1. - eta1*eta1)/E1 + (1. - eta2*eta2)/E2);
    	float Rstar = 1. / ((1./_agents[i]->_radius)   + (1./_agents[j]->_radius)  );

       	float distance = Distance(i ,j);
    	float overlapp = _agents[i]->_radius + _agents[j]->_radius - distance;

    	if( overlapp > 0 ){
       		float F = 4./3. *Estar * sqrt(Rstar) * pow(overlapp, 3./2.);

       		return 1./M_PI * pow( F*Estar*Estar/Rstar/Rstar, 1./3.);
    	}else
    		return 0;
    }


    inline float Pressure( int i, Boxes<Agent*> *boxes)
    {
    	float pressure = 0;

    	BoxIterator<Agent*> boxIt( boxes, this->at(i)->position(), 1);
    	for( ; !boxIt.end(); boxIt++){
    		int j = (*boxIt)->index();
    		if( j!=i){
    			pressure += Pressure( i, j);
    	    }
    	}

    	return pressure;
    }


    void AccelerationChemotaxis( float *acceleration, int i)
    {
    	// Gradient of WNT ~ Ephrin-B1 ligand
    	float alpha = 1*10000;
    	float beta  = 0.5*10000;

    	switch( at(i)->type()){

    	case Agent::Paneth:
    		if( Agent::EphB3Knockout == 0){
    			for( int d=0; d<_dimensions; d++)
    				acceleration[d] +=  alpha * Agent::gradEphrinB(at(i)->position(), d) ;
    		}
    		if( Agent::EphB3Knockout == +1){
    			for( int d=0; d<_dimensions; d++)
    				acceleration[d] +=  2*alpha * Agent::gradEphrinB(at(i)->position(), d) ;
    		}
			break;

    	//case Agent::StemCell:
    	case Agent::Secretory:
    	case Agent::Absorptive:
    	case Agent::Goblet:
    	case Agent::Enterocyte:
    	case Agent::Enteroendocrine:
    		if( Agent::EphB2Knockout == 0){
    			for( int d=0; d<_dimensions; d++)
    				acceleration[d] -=  beta * Agent::gradEphrinB(at(i)->position(), d) ;
    		}
    		if( Agent::EphB2Knockout == +1){
    			for( int d=0; d<_dimensions; d++)
    				acceleration[d] -=  2*beta * Agent::gradEphrinB(at(i)->position(), d) ;
    		}

			break;
    	}

    	/*fprintf( stderr, "(%e, %e, %e)\n",
    			Agent::gradEphrinB(at(i)->position(), 0),
    			Agent::gradEphrinB(at(i)->position(), 1),
    			Agent::gradEphrinB(at(i)->position(), 2));*/
    	return ;
    }

    void Acceleration( float *acceleration, int i, int j, float **adhesionEnergies)
    {
    	float F;
    	ExtendedHertzForce( i, j, adhesionEnergies, F);
       	float distance = Distance(i,j);
       	float mass = 1;

       	for( int d=0; d<_dimensions; d++)
       		acceleration[d] = - F/mass * (_agents[j]->_position[d]-_agents[i]->_position[d]) / distance;

       	if( isnan(acceleration[0]) || isnan( acceleration[1]) || isinf(acceleration[0]) || isinf( acceleration[1])){
       		//fprintf( stderr, "positions %i,%i = (%e,%e), (%e,%e) \n", i, j, position[i][0],position[i][1],position[j][0],position[j][1]);
       		fprintf( stderr, "distance |%i(%i),%i(%i)| = %e\n", i, _agents[i]->type(), j, _agents[j]->type(), distance);
       		fprintf( stderr, "(Hertz) force = %e\n", F);
       		fprintf( stderr, "acceleration = (%e, %e)\n", acceleration[0], acceleration[1]);
       		fprintf( stderr, "position(%i) = (%e, %e)\n", i, _agents[i]->_position[0], _agents[i]->_position[1]);
       		fprintf( stderr, "position(%i) = (%e, %e)\n", j, _agents[j]->_position[0], _agents[j]->_position[1]);
       	}
    	return ;
    }

    void Acceleration( float *acceleration, int i, int j, float F)
    {
       	float distance = Distance(i,j);
       	float mass = 1;

       	if(distance==0){
       		acceleration[0] = 0;
       		acceleration[1] = 0;
       		acceleration[2] = 0;
       		return;
       	}


       	for( int d=0; d<_dimensions; d++)
       		acceleration[d] = - F/mass * (_agents[j]->_position[d]-_agents[i]->_position[d]) / distance;

       	if( isnan(acceleration[0]) || isnan( acceleration[1]) || isinf(acceleration[0]) || isinf( acceleration[1])){
       		//fprintf( stderr, "positions %i,%i = (%e,%e), (%e,%e) \n", i, j, position[i][0],position[i][1],position[j][0],position[j][1]);
       		fprintf( stderr, "distance |%i,%i| = %e\n", i, j, distance);
       		fprintf( stderr, "force = %e\n", F);
       		fprintf( stderr, "acceleration = (%e, %e)\n", acceleration[0], acceleration[1]);
       	}
    	return ;
    }

    inline float Distance( int i, int j)
    {
    	float temp = 0.;

    	for( int d=0; d<_dimensions; d++)
    		temp += pow(_agents[j]->_position[d] - _agents[i]->_position[d], 2);

       	return sqrt( temp);
    }

    inline float SqrDistance( int i, int j)
    {
    	float temp = 0.;

    	for( int d=0; d<_dimensions; d++)
    		temp += pow(_agents[j]->_position[d] - _agents[i]->_position[d], 2);

       	return temp;
    }

    inline float ContactAreaDistance(  int i, int j){
    	float dij = Distance( i, j);
    	float Ri = _agents[i]->radius();
    	float Rj = _agents[j]->radius();
    	return ( pow( Ri, 2) - pow( Rj, 2) + pow( dij, 2)) / ( 2*dij );
    }
    /*inline float dContactAreaDistance_dDistance(  int i, int j){
    	float dij = Distance( i, j);
    	float Ri = _agents[i]->radius();
    	float Rj = _agents[j]->radius();

    	float a = pow( Ri, 2) - pow( Rj, 2);
    	float x = dij;
    	//return ( a + pow( x, 2)) / ( 2*x );
    	return 0.5 * (1. - a / pow( x, 2));
    }*/

    inline float ContactArea(  int i, int j)
    {
    	float Ri = _agents[i]->radius();
    	float xij = ContactAreaDistance( i, j);
    	return M_PI * ( pow(Ri, 2) - pow( xij, 2));
    }
    /*inline float dContactArea_dDistance(  int i, int j)
    {
    	float Ri = _agents[i]->radius();
    	float xij = dContactAreaDistance_dDistance( i, j);
    	return M_PI * ( pow(Ri, 2) - pow( xij, 2));
    }*/

    inline float AdhesiveInteractionEnergy( int i, int j)
    {
    	float adhesionEnergyPerUnitContactArea = 200; // 10^-6 N / m

    	return adhesionEnergyPerUnitContactArea * ContactArea( i, j);
    }
    inline float dAdhesiveInteractionEnergy_dDistance( int i, int j)
    {
    	float dij = Distance( i, j);
    	float Ri = _agents[i]->radius();
    	float Rj = _agents[j]->radius();

    	float adhesionEnergyPerUnitContactArea = 20; // 10^-6 N / m

    	if( _agents[j]->type() == Agent::substrate)
    		adhesionEnergyPerUnitContactArea = 1000;

    	if( Ri + Rj > dij)
    	return adhesionEnergyPerUnitContactArea * M_PI
    			* (pow( Rj, 2) - pow( dij, 2) + 1.)
    			* (pow( Rj, 2) + pow( dij, 2) + 1.)
    			/ (2. * pow( dij, 3));
    	else return 0;
    }
    inline float dAdhesiveInteractionEnergyBasalMembrane_dDistance( int i, int j)
    {
    	float dij = Distance( i, j);
    	float Ri = _agents[i]->radius();
    	float Rj = _agents[j]->radius();

    	float adhesionEnergyPerUnitContactArea = 5.5 * 1e-12; // Nm

    	if( Ri + Rj > dij)
    	return -adhesionEnergyPerUnitContactArea * (0.95/dij + Ri);
    	else return 0;
    }

    inline float DeformationEnergy( int i, int j)
    {
    	float E1 = getElasticModulus(i),
    		  E2 = getElasticModulus(j);
    	float eta1 = getPoissionRatio(i),
    		  eta2 = getPoissionRatio(j);
    	float Estar = 1. / ((1. - eta1*eta1)/E1 + (1. - eta2*eta2)/E2);

    	float D = 3./2. * 1. / Estar;

    	float dij = Distance( i, j);
    	float Ri = _agents[i]->radius();
    	float Rj = _agents[j]->radius();

    	return 2./5.*pow(Ri + Rj - dij, 5./2.) / D * sqrt( Ri*Rj/(Ri+Rj));
    }
    inline float dDeformationEnergy_dDistance( int i, int j)
    {
    	float E1 = getElasticModulus(i),
    		  E2 = getElasticModulus(j);
    	float eta1 = getPoissionRatio(i),
    		  eta2 = getPoissionRatio(j);
    	float Estar = 1. / ((1. - eta1*eta1)/E1 + (1. - eta2*eta2)/E2);

    	float D = 3./2. * 1. / Estar;

    	float dij = Distance( i, j);
    	float Ri = _agents[i]->radius();
    	float Rj = _agents[j]->radius();
//fprintf(stderr, "Ri(%e) + Rj(%e) > dij(%e)\n", Ri,Rj,dij);
    	//return 2./5.*pow(Ri + Rj - dij, 5./2.) / D * sqrt( Ri*Rj/(Ri+Rj));
    	if( Ri + Rj > dij)
    		return - pow(Ri + Rj - dij, 3./2.) / D * sqrt( Ri*Rj/(Ri+Rj));
    	else
    		return 0;
    }


    inline float CompressionEnergy( int i, BoxIterator<AgentType*> &boxIt)
    {
    	float Ri = _agents[i]->radius();

    	// Volume of Agent if isolated
    	float targetVolume = 4./3. * M_PI * pow(Ri, 3);

    	// Volume Overlap with other objects
    	float overlappedVolume = 0.;
    	for( ; !boxIt.end(); boxIt++){
    		int j = *boxIt;
    		if( j!=i){
    			fprintf(stderr, "neighbor %i\n", j);
    			overlappedVolume += M_PI / 3. * pow(Ri - ContactAreaDistance(i,j), 2) * (2*Ri - ContactAreaDistance(i,j));
    		}
    	}

    	return BulkModulus( i) / (2.*targetVolume) * pow(overlappedVolume, 2);
    }
    inline float dCompressionEnergy_dDistance( int i, int j)
    {
    	float Ri = _agents[i]->radius();

    	// Volume of Agent if isolated
    	float targetVolume = 4./3. * M_PI * pow(Ri, 3);

    	// Volume Overlap with other objects
    	float overlappedVolume = 0.;
    	//for( ; !boxIt.end(); boxIt++){
    	//	int j = *boxIt;
    	//	if( j!=i){
    	    	float dij = Distance( i, j);
    	    	float Rj = _agents[j]->radius();

    			//fprintf(stderr, "neighbor %i\n", j);
    			overlappedVolume += M_PI / 24.
    					* (3.*pow(Rj,2)-3.*dij+10.*Ri*dij+3)
    					* (Rj*Rj + dij*dij + 1)
    					* (Rj*Rj - pow(dij-Ri,2))
    					/ pow(dij,4);
    	//	}
    	//}

    	if( Ri + Rj > dij)
    	return BulkModulus( i) / (2.*targetVolume) * pow(overlappedVolume, 2);
    	else return 0;
    }



	/* CAPACITY */
	bool empty() {return (_size==0);}
	int  size()  {return  _size;}


	/* MODIFIERS */
	void clear()
	{
		for( int a=0; a<_size; a++)
			delete _agents[a];
	}

	AgentType *insert( AgentType *agent){
		if( _size == _sizeLimit){
			_sizeLimit += _sizeIncrement;
			_agents = (AgentType**) realloc( _agents, _sizeLimit * sizeof(AgentType*));
		}
		_agents[_size] = agent;
		_agents[_size]->set_index(_size);
		_size++;

		// molecules

		return agent;
	};
	void erase( int index) {
		delete _agents[index];
		_size--;
		_agents[index] = _agents[_size];
		_agents[index]->set_index(index);
	}
	void erase( AgentType *agent) {
		_size--;
		_agents[_size]->_index  = agent->_index;
		_agents[agent->_index] = _agents[_size];
		for( int d=0; d<nbmolecules; d++)
			molecules[agent->_index*nbmolecules+d]   = molecules[_size*nbmolecules+d];
		//_agents[agent->_index]->_index = agent->_index;
		delete agent;
	}


	/* ELEMENT ACCESS */
	AgentType* at( int index)         {return _agents[index];}
	AgentType* operator[]( int index) {return _agents[index];}
};

template <class AgentType>
//const float AgentList<AgentType>::pressureThreshold = 1e4; // no inhibition
//const float AgentList<AgentType>::pressureThreshold = 6e2; // contact-inhibition
 float AgentList<AgentType>::pressureThreshold = 3e2; // strong contact-inhibition

#endif /* AGENT_HPP_ */
