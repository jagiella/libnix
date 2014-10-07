/*
 * model.cpp
 *
 *  Created on: Aug 13, 2014
 *      Author: jagiella
 */




/*
 * main.cpp
 *
 *  Created on: 29.04.2014
 *      Author: jagiella
 */
#include <functional>
#include <iostream>

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <vector>

#include "hash.hpp"
#include "statistics.hpp"
#include "model.hpp"

enum CellPhenotypes{ FREE=0, DIVIDING=2, QUESCENT=3, DEAD=1};
enum CellProcessTypes{ DIVISION=0, NECROSIS=1, REENTER=2};
enum StatType{ Proliferation, ECM};

#define myrand( p_rseed) (rand_r(p_rseed)/(RAND_MAX+1.))
#define dmax(a,b) (a>b?a:b)
#define NO_FILE_OUTPUT





class Cell;
class Process {
public:
	int    _type;
	double _rate;
	Cell  *_cell;
	Process( Cell* cell, int type, double rate=0) : _type(type), _rate(rate), _cell(cell) {};
};

class Cell {
public:
	int _index;
	int _type;
	int _pos;
	Process *processes[3];
	Cell( int pos, int type, model_input minp, int idx=0) : _index(idx), _type(type), _pos(pos) {
		processes[0] = new Process( this, DIVISION, minp.k_div);
		processes[1] = new Process( this, NECROSIS, minp.k_nec);
		processes[2] = new Process( this, REENTER,  minp.k_re);
	};
	~Cell(){
		free( processes[0]);
		free( processes[1]);
		free( processes[2]);
	};
	static void add( Cell* cell, Cell* cells[], int n_cells, Cell* lattice[], int n_lattice){

	}
};


#include <stdint.h>




void doThatStat( Cell* lattice[], double molecule[], int size, double contribution, double stat[], int length, double &progress, const char* filename, char stat_type, double *stat_out){
	// find outer border
	int offset=size-1;
	while( lattice[offset]==0) offset--;

	progress += contribution;
	double partial_contribution = ( progress < 1. ? contribution : contribution - (progress-1) );

	for( int i=0; i<length; i++)
		if( i<=offset && lattice[offset-i])
		switch(stat_type){
		case Proliferation:
			stat[i] += partial_contribution*(lattice[offset-i]->_type == DIVIDING); break;
		case ECM:
			stat[i] += partial_contribution*molecule[offset-i]; break;
		}

	// output and reset
	if( progress >= 1.){
#ifndef NO_FILE_OUTPUT
		FILE *fp = fopen( filename, "w+");
#endif

		for( int i=0; i<length; i++){
			// write
#ifndef NO_FILE_OUTPUT
			fprintf( fp, "%e\n", stat[i]);
#endif
			if(stat_out)
				stat_out[i] = stat[i];

			// reset
			if( i<=offset && lattice[offset-i])
			switch(stat_type){
				case Proliferation:
				stat[i] = (contribution-partial_contribution)*(lattice[offset-i]->_type == DIVIDING); break;
				case ECM:
				stat[i] = (contribution-partial_contribution)*molecule[offset-i]; break;
			}else
				stat[i] = 0;
		}
#ifndef NO_FILE_OUTPUT
		fclose( fp);
#endif

		progress -= floor(progress);

	}
}

double oxyCon( double g, double o) {	double qmax=22, qmin=10, ko=0.03, kg=0.1;
	return qmax * o/(ko+o) * (1-(1-qmin/qmax) * g/(kg+g));
}

double glcCon( double g, double o) {	double qmax=54, qmin=15, ko=0.03, kg=0.07;
	return qmax * g/(kg+g) * (1-(1-qmin/qmax) * o/(ko+o));
}

double atpPro( double g, double o) {
	return 2*glcCon(g,o) + 17./3.*oxyCon(g,o);
}

double updateRates(
		Cell* lattice[], int N, std::vector<Cell*> cells, model_input minp,  // input
		double oxy[], double glc[], // multiscale
		Process* processes[], int &M // output
		){
	//printf( "update rates\n");
	M=0;
	double sum=0;
	for( unsigned int i=0; i<cells.size(); i++){
		switch( cells[i]->_type){

		case DIVIDING:
			if( atpPro( glc[cells[i]->_pos], oxy[cells[i]->_pos]) >= minp.atp_min){
				processes[M] = cells[i]->processes[DIVISION];
				sum += processes[M]->_rate;
				M++;
			}else{
				cells[i]->_type = QUESCENT;
				processes[M] = cells[i]->processes[NECROSIS];
				sum += processes[M]->_rate;
				M++;
			}
			break;
		case QUESCENT:
			if( atpPro( glc[cells[i]->_pos], oxy[cells[i]->_pos]) < minp.atp_min){
				//printf( "add NECROSIS\n");
				processes[M] = cells[i]->processes[NECROSIS];
				sum += processes[M]->_rate;
				M++;
			}else{
				processes[M] = cells[i]->processes[REENTER];
				sum += processes[M]->_rate;
				M++;
			}

			break;
		}
	}
	//printf( "...finished\n");

	return sum;
}

void findRates(
		Process* processes[], int &M, double partial_sum, // input
		Process* &process // output
		){

	double sum=0;

	for(int i=0; i<M; i++){
		sum += processes[i]->_rate;
		if( sum >= partial_sum ){
			process = processes[i];
			return;
		}
	}

	fprintf( stderr, "error findRates(): M=%i (sum=%e, parsum=%e)\n", M, sum, partial_sum);
}

void closestFree(
		Cell* lattice[], int N, int idx1,
		int &idx2, int &dist
		){
	for( int i=1; i<N; i++){
		int iminus = fmax( idx1-i, 0);
		int iplus  = fmin( idx1+i, N-1);
		if( lattice[iminus] == 0){
			idx2=iminus; dist=i; return;
		}
		if( lattice[iplus] == 0){
			idx2=iplus; dist=i; return;
		}
		if( iminus==0 && iplus == N-1){
			idx2 = -1; dist=N; return;
		}
	}
}

void shift(
		Cell* lattice[], int N, int origin, int target
		){

    if( origin < target)
        for( int i=target;i>origin; i--){
            lattice[i]=lattice[i-1];
            lattice[i]->_pos = i;
        }
    else
        for( int i=target;i<origin; i++){
            lattice[i]=lattice[i+1];
            lattice[i]->_pos = i;
        }
    lattice[ origin] = 0;
}

void perform( Cell* lattice[], int N, std::vector<Cell*> &cells,
		double glc[], double oxy[], double ecm[],
		Process* process, model_input minp, unsigned int *p_rseed ){

	int origin = process->_cell->_pos;

	switch( process->_type){
	case DIVISION:{
		//fprintf( stderr, "DIVISION: c %i\n", origin);
		int target=0, dist=0;
		closestFree(
				lattice, N, origin,
				target, dist
				);
		//printf( "closestFree(%i)=%i\n", origin, target);

		if( exp(-dist*16.8/minp.delta_L) < myrand( p_rseed) || ecm[origin] < minp.ecm_min ){
			lattice[ origin ]->_type = QUESCENT;
		}
		if( target!=-1){
			// shift all cells away from origin
			shift( lattice, N, origin, target);

			Cell* daughter;
			//if( dist > delta_L || ecm[origin] < ecm_min )
			if( exp(-dist*16.8/minp.delta_L) < myrand( p_rseed) || ecm[origin] < minp.ecm_min )
				daughter = new Cell( origin, QUESCENT, minp, cells.size());
			else
				daughter = new Cell( origin, DIVIDING, minp, cells.size());

			lattice[origin] = daughter;
			cells.push_back( daughter);

		}else
			lattice[ origin ]->_type = QUESCENT;
	}break;

	case NECROSIS:
	{
		//fprintf( stderr, "NECROSIS: c %i\n", origin);
		// detach from lattice
		lattice[ origin ] = 0;

		// erase from cell list
		Cell* remove = process->_cell;
		cells[ remove->_index ]         = cells[ cells.size()-1 ];
		cells[ remove->_index ]->_index = remove->_index;
		cells.pop_back();
		delete( remove );

		//fprintf( stderr, "NECROSIS: l %i -> c %i/%i\n", origin, idx, cells.size());
	}break;

	case REENTER:
	{
		int target = 0, dist = 0;
		closestFree(
				lattice, N, origin,
				target, dist
		);
		if( exp(-dist*16.8/minp.delta_L) < myrand( p_rseed) || ecm[origin] < minp.ecm_min )
			lattice[ origin ]->_type = QUESCENT;
		else
			lattice[ origin ]->_type = DIVIDING;
	}break;

	}

}


int countCells( Cell* lattice[], int N)
{
	int last_cellCount=0;
	for( int i=0; i<N; i++)
		if( lattice[i] && (lattice[i]->_type == DIVIDING || lattice[i]->_type == QUESCENT))
			//last_cellCount++;
			last_cellCount=i;
	return last_cellCount;
}

/*void setModelOutput( model_output &mout){
	mout.growth_curve = mout.all;
	mout.ECM_17   = &mout.all[600];
	mout.ECM_24   = &mout.all[600+300];
	mout.KI67_17  = &mout.all[600+300+300];
	mout.KI67_24  = &mout.all[600+300+300+300];
	mout.TUNEL_17 = &mout.all[600+300+300+300+300];
	mout.TUNEL_24 = &mout.all[600+300+300+300+300+300];

}*/

double* model( int parc, double *parv)
{

	// PARAMETERS
	// initial condition
	int InitialRadius = 50;
	float InitialQuiescent = 0;

	// cell kinetics
	model_input minp;
	minp.delta_L = 100;
	minp.k_div   = 1./24.;

	// ecm
	minp.USE_ECM = false;
	minp.pecm = 0;
	minp.qecm = 0;
	minp.ecm_min = 0;

	// atp
	minp.USE_ATP = false;
	minp.atp_min = 0;

	minp.k_nec   = 0.01;
	minp.k_re    = 0.01;

	unsigned int rseed = 0;
	for( int i=0; i<parc; i++)
	switch( i+1){
	case 9:
		minp.atp_min = (parv[8]); 	if(minp.atp_min) minp.USE_ATP = true; break;
	case 8:
		minp.ecm_min = (parv[7]); 	if(minp.ecm_min) minp.USE_ECM = true; break;
	case 7:
		minp.qecm    = (parv[5]); 	if(minp.qecm) minp.USE_ECM = true; break;
	case 6:
		minp.pecm    = (parv[6]); 	if(minp.pecm) minp.USE_ECM = true; break;
	case 5:
	 	InitialQuiescent =   (parv[4]); break;
	case 4:
	 	InitialRadius = fmax( (int)(parv[3]), 0); break;
	case 3:
		minp.delta_L=        fmax( (int)(parv[2]), 0); break;
	case 2:
		minp.k_div  =        fmax( (parv[1]), 0); break;
	case 1:
		rseed = (unsigned int) parv[0]; break;
	 	//srand( (int)( parv[0]));
	}
	//fprintf(stderr, "seed = %i\n", rseed);
	//srand( rseed);
	//srand( time(NULL));

 	//fprintf(stderr, "%.3e %i %i %.3e %.3e %.3e %.3e %.3e\n", minp.k_div, minp.delta_L, InitialRadius, InitialQuiescent, minp.pecm,minp.qecm,minp.ecm_min,minp.atp_min);

	// declare variables
	int N=1000;
	Cell* lattice[N];
	std::vector<Cell*> cells;

	// initial condition
	for( int i=0; i<InitialRadius && i<N; i++){
		Cell* cell = new Cell(i, DIVIDING, minp, cells.size());
		if( myrand( &rseed) < InitialQuiescent)
			cell->_type = QUESCENT;
		cells.push_back(cell);
		lattice[i]= cell;
	}

	for( int i=InitialRadius; i<N; i++)
		lattice[i] = 0;

	// stat
	int last=0;
	int last_cellCount   = InitialRadius;
	bool STAT_OVER_TIME  = true;
	bool STAT_OVER_CELLS = !STAT_OVER_TIME;

	int    stat_size = 300;
	double stat[stat_size];    for( int i=0; i<stat_size; i++) stat[i] = 0;
	double stat_progress = 0;
	double statECM[stat_size]; for( int i=0; i<stat_size; i++) statECM[i] = 0;
	double statECM_progress = 0;

	double *mout = (double*) malloc( 2400 * sizeof(double));
	double *growth_curve = mout;
	double *KI67_17  = &mout[600];
	double *KI67_24  = &mout[600+300];
	double *ECM_17   = &mout[600+300+300];
	double *ECM_24   = &mout[600+300+300+300];
	//double *TUNEL_17 = &mout[600+300+300+300+300];
	//double *TUNEL_24 = &mout[600+300+300+300+300+300];
	for( int i=0; i<2400; i++)
		mout[i] = 0;

	growth_curve[0] = InitialRadius;

#ifndef	NO_FILE_OUTPUT
	FILE *fp_growthcurve = fopen( "tmp.dat", "w+");
#endif

	// simulate
	float t=0, dt;
	float tend=600;

	// molecules
	double oxy[N], doxy[N], Doxy=6300000/pow(16.8,2), oxyB=0.28; for( int i=0; i<N; i++) oxy[i] = oxyB;
	double glc[N], dglc[N], Dglc=378000 /pow(16.8,2), glcB=25;   for( int i=0; i<N; i++) glc[i] = glcB;
	double ecm[N], decm[N], Decm=100 /pow(16.8,2),    ecmB=0.;  for( int i=0; i<N; i++) ecm[i] = 0.;
	float last_update = 0;
	float update_interval = 0.1;//0.5 * fmin( 1/Doxy, 1/Dglc);
	//fprintf( stderr, "update_interval = %e\n", update_interval);

	// rates
	double k_sum;
	int   processcount=0;
	Process* processes[10000];
	k_sum = updateRates( lattice, N, cells, minp, // input
			oxy, glc,
			processes, processcount);

	while( (STAT_OVER_TIME  && t<tend) ||
		   (STAT_OVER_CELLS && last_cellCount<N))
	{
		//fprintf( stderr, "k_sum = %e, cells=%i\n", k_sum, cells.size());
		//for( int c=0; )

	   // gillespie
		if( processcount){
			//fprintf( stderr, "findRates\n ");
		   double k_partial_sum = myrand( &rseed)*k_sum;
		   Process* process=0;
		   findRates(
				   processes, processcount, k_partial_sum, // input
				   process // output
				);

		   // perform process
		   //fprintf( stderr, "perform process\n ");
		   perform( lattice, N, cells, glc, oxy, ecm, process, minp, &rseed);

		   for( unsigned int c=0; c<cells.size(); c++){
			   if( lattice[ cells[c]->_pos ] == 0){ // all cells non-free?
				   fprintf( stderr, "ERROR!!!");
			   }
		   }
		   unsigned int count_nonfree=0;
		   for( int i=0; i<N; i++){
			   if( lattice[i] != 0)
				   count_nonfree++;
		   }
		   if(count_nonfree!= cells.size())
			   fprintf( stderr, "ERROR2: too many non-free (%i != %i)!!!", count_nonfree, (int)cells.size());

		   // update time
		   dt = fmin( 1./k_sum * log(1./(1.-myrand( &rseed)) ), tend-t);
		   if( k_sum==0) dt = tend-t;
		   t = t + dt;

	   }else{
		   dt = 1;
		   t = t + dt;
	   }

	   // update rates
	   k_sum = updateRates( lattice, N, cells, minp, // input
			   	   oxy, glc,
			   	   processes, processcount);


	   if( 16< t/24 && t/24<17.9 ){
		   doThatStat( lattice, ecm, N, dt/24., statECM, stat_size, statECM_progress, "radial17ECM.dat", ECM, ECM_17);
		   doThatStat( lattice, 0,   N, dt/24., stat,    stat_size, stat_progress,    "radial17.dat", Proliferation, KI67_17);
	   }
	   if( 23< t/24 && t/24<24.9 ){
		   doThatStat( lattice, ecm, N, dt/24., statECM, stat_size, statECM_progress, "radial24ECM.dat", ECM, ECM_24);
	   	   doThatStat( lattice, 0,   N, dt/24., stat,    stat_size, stat_progress,    "radial24.dat", Proliferation, KI67_24);
	   }

	   // update molecule conc.
	   if(t >= last_update + update_interval){

		   if( minp.USE_ECM)
		   while(t >= last_update + update_interval){
			   // diffusion
			   decm[0] = 0;
			   for(int i=0; i<N-1; i++){
				   decm[i] += (ecm[i+1]-ecm[i]) * Decm - minp.qecm*ecm[i]; decm[i+1] = (ecm[i]-ecm[i+1]) * Decm;
			   }
			   decm[N-1] += (ecmB-ecm[N-1]) * Decm - minp.qecm*ecm[N-1];

			   // reaction
			   for( unsigned int c=0; c<cells.size(); c++){ int alive = (cells.at(c)->_type==DIVIDING || cells.at(c)->_type==QUESCENT);
				   int i = cells.at(c)->_pos;
				   decm[i] += minp.pecm*alive;
			   }

			   for(int i=0; i<N; i++)
				   ecm[i] += decm[i] * update_interval;
			   last_update += update_interval;
		   }

		   // STEADY STATE
		   double max_residual = 1e-5;
		   double residual = max_residual;
		   if( minp.USE_ATP)
		   while(residual >= max_residual){
			   residual=0;
			   // diffusion
			   int i=0; int consuming = (lattice[i] && (lattice[i]->_type==DIVIDING || lattice[i]->_type==QUESCENT));
			   doxy[i] = (oxy[i+1])*(consuming*10+1)*Doxy / (consuming*oxyCon( glc[i], oxy[i])/oxy[i] + (consuming*10+1)*Doxy) - oxy[i];
			   dglc[i] = (glc[i+1])*(consuming*10+1)*Dglc / (consuming*glcCon( glc[i], oxy[i])/glc[i] + (consuming*10+1)*Dglc) - glc[i];

			   for( i=1; i<N-1; i++){consuming = (lattice[i] && (lattice[i]->_type==DIVIDING || lattice[i]->_type==QUESCENT));
				   doxy[i] = (oxy[i-1]+oxy[i+1])*(consuming*10+1)*Doxy / (consuming*oxyCon( glc[i], oxy[i])/oxy[i] + 2*(consuming*10+1)*Doxy) - oxy[i];
				   dglc[i] = (glc[i-1]+glc[i+1])*(consuming*10+1)*Dglc / (consuming*glcCon( glc[i], oxy[i])/glc[i] + 2*(consuming*10+1)*Dglc) - glc[i];
			   }

			   i=N-1;consuming = (lattice[i] && (lattice[i]->_type==DIVIDING || lattice[i]->_type==QUESCENT));
			   doxy[i] = (oxyB)*(consuming*10+1)*Doxy / (consuming*oxyCon( glc[i], oxy[i])/oxy[i] + (consuming*10+1)*Doxy) - oxy[i];
			   dglc[i] = (glcB)*(consuming*10+1)*Dglc / (consuming*glcCon( glc[i], oxy[i])/glc[i] + (consuming*10+1)*Dglc) - glc[i];

			   for(i=0; i<N; i++){
				   oxy[i] += doxy[i];
				   glc[i] += dglc[i];

				   residual = fmax( residual, fabs(doxy[i]));
			   }

			   //printf( "residual=%e\n ", residual);
		   }
		   //last_update = t;
	   }

	   // statistics
	   //fprintf( stderr, "statistics\n ");
	   if( STAT_OVER_TIME && ceil(t-dt) != ceil(t)){
		   for( int i=last+1; i<(int)ceil(t); i++ ){
#ifndef	NO_FILE_OUTPUT
			   fprintf( fp_growthcurve, "%i\n ", last_cellCount);
#endif
			   growth_curve[i] = last_cellCount;
		   }
		   last=ceil(t);
		   last_cellCount=countCells( lattice, N);
#ifndef	NO_FILE_OUTPUT
		   fprintf( fp_growthcurve, "%i\n ", last_cellCount);
#endif
		   growth_curve[last] = last_cellCount;
	   }

	   if( STAT_OVER_CELLS && last_cellCount != countCells( lattice, N))
	   {
		   last_cellCount = countCells( lattice, N);
#ifndef	NO_FILE_OUTPUT
		   fprintf( fp_growthcurve, "%e\n ", t);
#endif
	   }

#ifndef	NO_FILE_OUTPUT
	   {
		   FILE *fp = fopen( "arrangement.dat", "w+");
		   for( int i=0; i<N; i++)
			   fprintf( fp, "%i %i %e %e %e %e\n", i, (lattice[i]==0 ? 0 : lattice[i]->_type), glc[i], oxy[i], atpPro(glc[i], oxy[i]), ecm[i]);
		   fclose( fp);
	   }
#endif
	}
#ifndef	NO_FILE_OUTPUT
	fprintf( fp_growthcurve, "\n");
	fclose( fp_growthcurve);
#endif
	//fprintf( stderr, "finished\n ");

	// FREE MEMORY
	while( !cells.empty()){
		delete( cells.at( cells.size()-1));
		cells.pop_back();
	}

	return mout;
}





/*void average( double **x, int n, int m, double *m){
	for( int j=0; j<m; j++){
		m[j] = 0;
		for( int i=0; i<n; i++)
			m[j] += x[i][j];
		m[j] /= n;
	}
}*/






double compare( int parc, double *parv1, double *parv2, int n1, int n2){
	char filename[1024];
	FILE *fp;
	int length = 600+300+300+300+300+300+300;

	double **m_out1 = (double**) malloc( n1 * sizeof(double*));
	double **m_out2 = (double**) malloc( n2 * sizeof(double*));

	double *m_av1, *s2_av1, *m_av2, *s2_av2;


	unsigned int seed = clock();//time(NULL);

	//int test=0;

	// DATA 1
	sprintf( filename, "%u.bin", (unsigned int) hash_array( parv1, parc));
	fp = fopen( filename, "rb");
	if( fp){
		for( int i=0; i<n1; i++){
			m_out1[i] = (double*) malloc( length * sizeof(double));
			fread( m_out1[i], sizeof(double), length, fp);
		}
	}else{
		for( int i=0; i<n1; i++){
			//seed = time(NULL);
			parv1[0] = rand_r(&seed);
			m_out1[i] = model( parc, parv1);
		}
		fp = fopen( filename, "wb");
		for( int i=0; i<n1; i++)
			fwrite( m_out1[i], sizeof(double), length, fp);
	}
	fclose(fp);


//	m_av1 = m_out[0];



	// DATA 2

	for( int i=0; i<n2; i++){
		//seed = time(NULL);
		parv2[0] = rand_r(&seed);
		m_out2[i] = model( parc, parv2);
	}



//	m_av2 = m_out[0];
	// COMPARISON

	//KolmogorovSmirnoffTest( 300, m_av1.growth_curve, 300, m_av2.growth_curve);

	m_av1  = mean( m_out1, n1, 2400, rowwise);
	s2_av1 = std2( m_out1, m_av1, n1, 2400, rowwise);
	m_av2 = mean( m_out2, n2, 2400, rowwise);
	s2_av2 = std2( m_out2, m_av2, n2, 2400, rowwise);
	double negloglik =	- logLikelihood( (double*)m_av1, (double*)m_av2, (double*)s2_av1, (double*)s2_av2, length);


	fp = fopen( "tumordata.dat", "w+");
	for( int i=0; i<length; i++)
		fprintf( fp, "%e %e\n", m_av2[i], s2_av2[i]);
	fclose( fp);

	/*double negloglik = 0;
	double x1[n1], x2[n2];
	for( int i=0; i<2400; i++){
		for( int j=0; j<n1; j++)
			x1[j] = m_out1[j][i];
		for( int j=0; j<n2; j++)
			x2[j] = m_out2[j][i];
		KSTestResult result = KolmogorovSmirnoffTest( n1, x1, n2, x2);
		negloglik += - log( result.pValue + 1e-20 );
	}*/

	for( int i=0; i<n1; i++)
		free(m_out1[i]);
	free( m_out1);
	for( int i=0; i<n2; i++)
		free(m_out2[i]);
	free( m_out2);
	free( m_av1);
	free( s2_av1);
	free( m_av2);
	free( s2_av2);

	//fprintf( stderr, "[ %e ]\n", negloglik);

	return negloglik;
}
