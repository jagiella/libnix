/*
 * PharmacoKinetics.h
 *
 *  Created on: 07.02.2010
 *      Author: jagiella
 */

#ifndef PHARMACOKINETICS_H_
#define PHARMACOKINETICS_H_

#include "SparseMatrix.h"
#include "VoronoiDiagramExtended.h"

class PharmacoKinetics {
private:
	//SparseMatrix *ACPT11in;
	//SparseMatrix *ACPT11out;
	SparseMatrix *A;
	float *b;
	SparseMatrix *Ain;
	float *bin;

	float *CPT11in;
	float *CPT11out;

	// border conditions
	float CPT11in_Border;
	float CPT11out_Border;

	float timeStep;
	float spaceStep;
	float diffusionCoefficientCPT11in;
	float diffusionCoefficientCPT11out;

	float V_in;
	float V_out;
	float V_blood;
	float k_uptake_CPT;
	float V_eff_CPT;
	float K_eff_CPT;


public:
	PharmacoKinetics();
	PharmacoKinetics(VoronoiDiagram *vd);
	virtual ~PharmacoKinetics();

	void initCPT11in( VoronoiDiagram *vd, float value);
	void initCPT11out( VoronoiDiagram *vd, float value);

	void setBorderConditionCPT11in( float value);
	void setBorderConditionCPT11out( float value);

	void setSystemExplicitCPT11in( VoronoiDiagram *vd, SparseMatrix *A, float *b, float *CPT11in, float *CPT11out, float spaceStep, float timeStep);
	void setSystemExplicitCPT11out( VoronoiDiagram *vd, SparseMatrix *A, float *b, float *CPT11in, float *CPT11out, float spaceStep, float timeStep);

	float evolveSystem( VoronoiDiagram *vd, float time, float timeEnd, float timeStep);
};

#endif /* PHARMACOKINETICS_H_ */
