/*
 * model.hpp
 *
 *  Created on: Aug 13, 2014
 *      Author: jagiella
 */

#ifndef MODEL_HPP_
#define MODEL_HPP_

using namespace std;


typedef struct {
	float k_div;
	float k_nec;
	float k_re;
	int   delta_L;

	bool  USE_ECM;
	float pecm;
	float qecm;
	float ecm_min;

	bool  USE_ATP;
	float atp_min;
} model_input;

double* model( int parc, double *parv);
double compare( int parc, double *parv1, double *parv2, int n1, int n2);

typedef struct{
	bool   H;
	double pValue;
	double KSstatistic;
} KSTestResult;

KSTestResult KolmogorovSmirnoffTest( int n1, double *x1, int n2, double *x2);

#endif /* MODEL_HPP_ */
