/*
 * PharmacoKinetics.cpp
 *
 *  Created on: 07.02.2010
 *      Author: jagiella
 */

#include "PharmacoKinetics.h"
//#include "SparseMatrix.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

PharmacoKinetics::PharmacoKinetics() {
	// TODO Auto-generated constructor stub
	fprintf( stderr, "new PharmacoKinetics instance\n");
}

PharmacoKinetics::PharmacoKinetics( VoronoiDiagram *vd) {
	// TODO Auto-generated constructor stub
	fprintf( stderr, "new PharmacoKinetics instance related to Voronoi Diagram (size:%i)\n", vd->countVoronoiCells);

	//this->ACPT11in  = SparseMatrix::newSparseMatrix( vd->countVoronoiCells, vd->countVoronoiCells);
	//this->ACPT11out = SparseMatrix::newSparseMatrix( vd->countVoronoiCells, vd->countVoronoiCells);
	this->A = SparseMatrix::newSparseMatrix( vd->countVoronoiCells, vd->countVoronoiCells);
	this->b  = (float*) malloc( sizeof(float) * vd->countVoronoiCells);
	this->Ain = SparseMatrix::newSparseMatrix( vd->countVoronoiCells, vd->countVoronoiCells);
	this->bin  = (float*) malloc( sizeof(float) * vd->countVoronoiCells);

	this->CPT11in  = (float*) malloc( sizeof(float) * vd->countVoronoiCells);
	this->CPT11out = (float*) malloc( sizeof(float) * vd->countVoronoiCells);

	V_in = 1; // micrometer^3
	V_out = 1./6.; // micrometer^3
	V_blood = 1.; // micrometer^3
	//V_out = 0.2; //312.5; // micrometer^3
	k_uptake_CPT = 5.93;//0.0231 * 312.5; //= 7.21875 // hour^-1
	V_eff_CPT = 1200; //1157; // micromolar * hour^-1
	K_eff_CPT = 70; //71; // micromolar

	diffusionCoefficientCPT11out = 1000000;

	spaceStep = 13.; // micrometer
	timeStep = 0.1; // hour


}

void PharmacoKinetics::initCPT11in( VoronoiDiagram *vd, float value){
	for( int c=0; c<vd->countVoronoiCells; c++){
		vd->voronoiCells[c]->CPT11in = value;
		this->CPT11in[c] = value;
	}
}

void PharmacoKinetics::initCPT11out( VoronoiDiagram *vd, float value){
	for( int c=0; c<vd->countVoronoiCells; c++){
		vd->voronoiCells[c]->CPT11out = value;
		this->CPT11out[c] = value;
	}
}

void PharmacoKinetics::setBorderConditionCPT11in( float value){
	this->CPT11in_Border = value;
}

void PharmacoKinetics::setBorderConditionCPT11out( float value){
	this->CPT11out_Border = value;
}

PharmacoKinetics::~PharmacoKinetics() {
	// TODO Auto-generated destructor stub
}

void PharmacoKinetics::setSystemExplicitCPT11out( VoronoiDiagram *vd, SparseMatrix *_A, float *_b, float *_CPT11in, float *_CPT11out, float spaceStep, float timeStep)
{
	int di   = 1;
	int dii  = vd->xN[0];
	int diii = vd->xN[0]*vd->xN[1];

	double r = diffusionCoefficientCPT11out * timeStep / (spaceStep*spaceStep);

	//fprintf( stderr, "diffusionCoefficientCPT11out is %lf\n", diffusionCoefficientCPT11out);
	//fprintf( stderr, "r is %lf\n", r);

	int m=0;

	for( int iii=0; iii<vd->xN[2]; iii++)
	for( int ii=0; ii<vd->xN[1]; ii++)
	for( int i=0; i<vd->xN[0]; i++)
	{
		_A->resetRow( m);

		double _q;

		// sink/source of [CPT11]out
		if( (i %10==0 && ii %10==0) || (i %10==0 && iii%10==0) || (ii%10==0 && iii%10==0)){
			// vessel compartment
			//double A=27.03;	// microM
			//double B=0.2828;// hour^-1
			double k_exchange = 1.5;
			//_q = - V_blood/V_out * k_exchange*2.5 * _CPT11out[m]
			//     + V_blood/V_out * k_exchange * this->CPT11out_Border /*A*exp(-B*time)*/ ;
			_q = - V_blood/V_out * k_exchange * _CPT11out[m]
			     +                 k_exchange * this->CPT11out_Border /*A*exp(-B*time)*/ ;
		}else{
			// cell compartment
			_q = - V_in/V_out * k_uptake_CPT*_CPT11out[m]
			     + V_in/V_out * V_eff_CPT*_CPT11in[m] / (K_eff_CPT + _CPT11in[m]);
		}


		// vessels
		/*if( (i %10==0 && ii %10==0) ||
			(i %10==0 && iii%10==0) ||
			(ii%10==0 && iii%10==0) )
		{
			_A->setLast(m, m, 1);
			_b[m] = this->CPT11out_Border;
		}

		// element inside domain
		else*/
		if( i>0   && i<vd->xN[0]-1 &&
			ii>0  && ii<vd->xN[1]-1 &&
			iii>0 && iii<vd->xN[2]-1){
			_A->setLast(m, m, 1);
			//double _q = - V_in/V_out * k_uptake_CPT*_CPT11out[m]
			//            + V_in/V_out * V_eff_CPT*_CPT11in[m] / (K_eff_CPT + _CPT11in[m]);

			_b[m] = _CPT11out[m]
			        + timeStep*_q
					- r*6*_CPT11out[m]
					+ r*( _CPT11out[m-diii] + _CPT11out[m-dii] + _CPT11out[m-di] + _CPT11out[m+di] + _CPT11out[m+dii] + _CPT11out[m+diii]);
		}

		// border condition
		else{
			_A->setLast(m, m, 1);

			// DIRICHELET
			//_b[m] = this->CPT11out_Border;

			// NEUMANN
			//double _q = - V_in/V_out * k_uptake_CPT*_CPT11out[m]
			//            + V_in/V_out * V_eff_CPT*_CPT11in[m] / (K_eff_CPT + _CPT11in[m]);
			_b[m] = _CPT11out[m]
			        + timeStep*_q
					- r*6*_CPT11out[m]
					+ r*((i>0?_CPT11out[m-di]:_CPT11out[m]) + (ii>0?_CPT11out[m-dii]:_CPT11out[m]) + (iii>0?_CPT11out[m-diii]:_CPT11out[m]) +
							(i<vd->xN[0]-1?_CPT11out[m+di]:_CPT11out[m]) + (ii<vd->xN[1]-1?_CPT11out[m+dii]:_CPT11out[m]) + (iii<vd->xN[2]-1?_CPT11out[m+diii]:_CPT11out[m]));

		}

		m++;
	}
}

void PharmacoKinetics::setSystemExplicitCPT11in( VoronoiDiagram *vd, SparseMatrix *A, float *b, float *CPT11in, float *CPT11out, float spaceStep, float timeStep)
{
	int m=0;

	for( int iii=0; iii<vd->xN[2]; iii++)
	for( int ii=0; ii<vd->xN[1]; ii++)
	for( int i=0; i<vd->xN[0]; i++)
	{
		A->resetRow( m);

		// element inside domain
		//if( i>0 && i<vd->xN[0]-1 && ii>0 && ii<vd->xN[1]-1 && iii>0 && iii<vd->xN[2]-1){
			A->setLast(m, m, 1);
			double q = k_uptake_CPT * CPT11out[m]
			         - V_eff_CPT * CPT11in[m] / (K_eff_CPT + CPT11in[m]);
			b[m] = CPT11in[m] + timeStep*q;
		//}

		// border condition
		/*else{
			A->setLast(m, m, 1);
			b[m] = this->CPT11in_Border;
		}*/
		m++;
	}
}

float PharmacoKinetics::evolveSystem( VoronoiDiagram *vd, float time, float timeEnd, float timeStep)
{
	//double time;

	fprintf( stderr, "evolveSystem\n");

	float timeInterval = timeEnd - time;

	for( ; time+timeStep <= timeEnd; time += timeStep){
		fprintf(stderr, "\r%lf\%% \b",(timeInterval - timeEnd + time)/timeInterval*100. );
		//fprintf( stderr, "time: %lf\n", time);

		//this->CPT11out_Border =
		//float C0 = 0;
				/*( time < 0.5 ?
						  27.*time/0.5
						: 27.*exp(-time*0.3)); //30.*exp(-time*0.38);
*/
		//float v1 = 1.;		// blood volume
		//float v2 = 1./3.;	// extracellular volume

		//float k = 0.5;		// diffusion coefficient: blood <-> extra-cellular medium
		float A = 27.03;
		float B = 0.2828;
		//float E = 1./(v1/v2 - B/k);

		//this->CPT11out_Border = (C0 - A*E)*exp( -k*v1/v2*time) + E*A*exp(-B*time);

		// Explicite
		this->CPT11out_Border = A*exp(-B*time);

		// Implicite
		this->CPT11out_Border = CPT11out_Border - B * CPT11out_Border * timeStep;

		//this->initCPT11out( vd, this->CPT11out_Border);

		// set up system
		setSystemExplicitCPT11out( vd, this->A, this->b, this->CPT11in, this->CPT11out, spaceStep, timeStep);
		setSystemExplicitCPT11in(  vd, this->Ain, this->bin, this->CPT11in, this->CPT11out, spaceStep, timeStep);

		// update
		SolveExplicit( this->A, this->b, this->CPT11out);
		SolveExplicit( this->Ain, this->bin, this->CPT11in);

		// copy
		/*for( int c=0; c<vd->countVoronoiCells; c++){
			vd->voronoiCells[c]->CPT11out = this->CPT11out[c];
			if( isnan(this->CPT11out[c])){
				fprintf( stderr, "this->CPT11out[%i] is nan\n", c);
				//fprintf( stderr, "m=%i\n", m);
				//fprintf( stderr, "CPT11out[%i]=%lf\n", m, _CPT11out[m]);
				//fprintf( stderr, "CPT11in[%i]=%lf\n", m, _CPT11in[m]);
				exit(0);
			}
			if( isinf(this->CPT11out[c])){
				fprintf( stderr, "this->CPT11out[%i] is inf\n", c);
				exit(0);
			}

			//fprintf(stderr, " %lf\n", this->CPT11out[c]);
		}*/

		if( ceil(time*10.) != ceil((time+timeStep)*10.)){

			// mean [CPT11out]
			float meanCPT11out = 0;
			float meanCPT11in = 0;

//#pragma omp parallel for
			for( int c=0; c<vd->countVoronoiCells; c++){
				meanCPT11in  += this->CPT11in[c];
				meanCPT11out += this->CPT11out[c];
			}
			meanCPT11in  /= vd->countVoronoiCells;
			meanCPT11out /= vd->countVoronoiCells;


			// output
			printf( "%lf %lf %lf %lf\n",
				time,
				meanCPT11out, //this->CPT11out[(int)vd->xN[0]*vd->xN[1]*vd->xN[2]/2 + (int)vd->xN[0]*vd->xN[1]/2 + (int)vd->xN[0]/2],
				meanCPT11in, //this->CPT11in[(int)vd->xN[0]*vd->xN[1]*vd->xN[2]/2 + (int)vd->xN[0]*vd->xN[1]/2 + (int)vd->xN[0]/2],
				this->CPT11out_Border);
		}
	}

#pragma omp parallel for
	for( int c=0; c<vd->countVoronoiCells; c++){
		vd->voronoiCells[c]->CPT11in = this->CPT11in[c];
		vd->voronoiCells[c]->CPT11out = this->CPT11out[c];
	}

	return timeEnd - time;
}
