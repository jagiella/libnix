/****************************************************************************
 * Includes                                                                 *
 ****************************************************************************/
#include "Montecarlo.h"

#include <stdio.h>            // Standard-Ein- und Ausgabe
#include <stdlib.h>
#include <stdarg.h>
#include <unistd.h>
#include <string.h>
#include <limits.h>
#include <float.h>
#include <time.h>
#include <math.h>
#include <sys/stat.h>

#include "VoronoiDiagramExtended.h"
//#include "VoronoiDiagram.h"
#include "SquareLattice2D.h"
#include "HexagonalLattice2D.h"
#include "Agent.h"
#include "CellProcesses.h"
#include "Mathematix.h"
//#include "CoarseGraining.h"
#include "Substrate.h"
#include "VesselNetwork.h"
#include "finiteDifferences.h"
#include "SparseMatrix.h"
#include "PharmacoKinetics.h"
#include "Interpolation.h"
#include "EpsIO.h"

#include "DiffusionReactionEquation.hpp"

#include "../matrix.hpp"
#include "../statistics.hpp"
#include "AsciiIO.h"

#define MAX_CELLS 200000
//#define TEST_STATIC_PROFIL

//#define PROB_REENTERING_CELL_CYCLE 1.//0.6//1.//0.7//0.66//1.//0.75
//#define PROB_REENTERING_CELL_CYCLE( DELTA_L) 1.//( exp( -DELTA_L/150.)) // micrometers
//#define PROB_REENTERING_CELL_CYCLE( DELTA_L) ( exp( -DELTA_L/100.)) // micrometers

//#define PROB_REENTERING_CELL_CYCLE( DELTA_L) ( pow(130.,1.)/(pow(130.,1.) + pow(DELTA_L,1.))* 0.8)//0.85 )

//#define PROB_REENTERING_CELL_CYCLE( DELTA_L) ( pow(200.,4.)/(pow(200.,4.) + pow(DELTA_L,4.))*1.0 )
//#define PROB_REENTERING_CELL_CYCLE( DELTA_L) (DELTA_L < 200. ? 0.9 : 0.)//( exp( -DELTA_L/150.)) // micrometers
//#define PROB_REENTERING_CELL_CYCLE( DELTA_L) (DELTA_L < 200. ? 0.7 : 0.)//( exp( -DELTA_L/150.)) // micrometers
//#define PROB_REENTERING_CELL_CYCLE( DELTA_L) ( 130./(130.+DELTA_L) - exp(-DELTA_L/5))

//#define PROB_REENTERING_CELL_CYCLE( DELTA_L) ( 100./(100.+DELTA_L) - 0.7*exp(-DELTA_L/1.0) )
//#define PROB_REENTERING_CELL_CYCLE( DELTA_L) (1-exp(-DELTA_L/1.0)) * exp(-DELTA_L/130.)

//#define PROB_REENTERING_CELL_CYCLE( DELTA_L) exp(-DELTA_L/130.)
//#define PROB_REENTERING_CELL_CYCLE( DELTA_L) 1
#define PROB_REENTERING_CELL_CYCLE( DELTA_L) ( ReentranceProbabilityFunction == HEAVYSIDE ? (DELTA_L <= ReentranceProbabilityReferenceLength ? 1 : 0) : exp(-DELTA_L/ReentranceProbabilityReferenceLength))

double TempTime = 0;

// TEXT STYLE
#define RESET		0
#define BRIGHT 		1
#define DIM			2
#define UNDERLINE 	3
#define BLINK		4
#define REVERSE		7
#define HIDDEN		8

// TEXT COLORS
#define BLACK 		0
#define RED			1
#define GREEN		2
#define YELLOW		3
#define BLUE		4
#define MAGENTA		5
#define CYAN		6
#define	WHITE		7


void fprinttextcolor(int fg) {
	/* Command is the control command to the terminal */
	//sprintf(command, "%c[%d;%d;%dm", 0x1B, attr, fg + 30, bg + 40);
	fprintf(stderr, "%c[;%d;m", 0x1B, fg + 30);
}

void fprinttextattribute(int attr) {
	/* Command is the control command to the terminal */
	//sprintf(command, "%c[%d;%d;%dm", 0x1B, attr, fg + 30, bg + 40);
	fprintf(stderr, "%c[%d;;m", 0x1B, attr);
}

void fprintbackgroundcolor(int bg) {

	/* Command is the control command to the terminal */
	fprintf(stderr, "%c[;;%dm", 0x1B, bg + 40);
}

int Case = 1;

typedef struct _Statisticx {
	int *count_cells;
	int *count_cell_volume;
	int *count_expanded_cells;
	double *count_divided_cells;
	int *count_vessel_cells;
	int *count_necrotic_cells;
	int *count_necrotic_volume;
	int *count_inner_empty_volume;
	int *count_lysed_cells;
	int *count_lysed_volume;
} Statistix;


void updateHistogram( AgentList *agentArray, VoronoiDiagram *vd, int* histogram, int* histogramDividing, int* histogramNecrotic, int* histogramFree, double* histogramECM)
{
	//fprintf(stderr, "[Print Cells to Povray File]\n");
	//Time = EndTime;
//	sprintf(outfilename, "%s/cellsR300.pov", dirname);
//	agentArray->printToPovray(outfilename, voronoiDiagram);
	//exit(0);

	// GET ALL BORDER CELLS
	int countBorderCells = 0;
	int maxBorderCells = 100000;

	VoronoiCell *borderCells[100000];
	//fprintf(stderr, "[Estimate Border]\n");

	for (int a = 0; countBorderCells!=maxBorderCells && a < agentArray->countActiveAgents; a++)
		for (int l = 0; countBorderCells!=maxBorderCells && l < agentArray->agents[a]->countLocations; l++){
			//bool unoccupiedNeighbor = false;
			for( int n=0; countBorderCells!=maxBorderCells && n<agentArray->agents[a]->location[l]->countNeighborCells /*&& !unoccupiedNeighbor*/; n++)
				if( !agentArray->agents[a]->location[l]->neighborCells[n]->agent
						|| agentArray->agents[a]->location[l]->neighborCells[n]->index < 0
						|| agentArray->agents[a]->location[l]->isDomainBorder(vd)){

					//unoccupiedNeighbor=true;
					//if( countBorderCells==100000)
					//	exit(0);
					borderCells[countBorderCells++] = agentArray->agents[a]->location[l]->neighborCells[n];
				}
		}
	//fprintf(stderr, "[Border Estimated]\n");

	for (int a = 0; a < agentArray->countActiveAgents; a++)
	//if( agentArray->agents[a]->state != FREE)
			{
		//fprintf(stderr, "\r[Agent %i]: State                 \b", a );
		// DIVIDING?
		bool dividing = false;
		bool necrotic = false;
		bool free = false;
		for (int l = 0;
				l < agentArray->agents[a]->countLocations;
				l++) {
			if (agentArray->agents[a]->state == FREE)
				free = true;
			if (agentArray->agents[a]->state == NECROTIC)
				necrotic = true;
			if (agentArray->agents[a]->state == ACTIVE
					&& agentArray->agents[a]->divide == 1
					&& !IsQuiescent(agentArray->agents[a]))
				dividing = true;
		}

		// distance to border
		VoronoiCell *border = 0;
		//if( VoronoiCell::SHIFT_TO_UNOCCUPIED)
		/*if( VoronoiCell::USE_SYMBOLIC_EXTENDED_NEIGHBORHOOD)
			border = voronoiDiagram->searchClosestUnoccupiedVoronoiCell(
				agentArray->agents[a]->location[0],
				VoronoiCell::extendedNeighborDepth,
				VoronoiCell::symbolicExtendedNeighborhoodSize,
				VoronoiCell::symbolicExtendedNeighborhood);
		else*/

		//fprintf(stderr, "\r[Agent %i]: Dist to Border                 \b", a );
		if(countBorderCells){
			border = borderCells[0];
			double distClosestBorder = border->getDistanceTo( agentArray->agents[a]->location[0]);
			for( int b=1; b<countBorderCells; b++){
				double dist = borderCells[b]->getDistanceTo( agentArray->agents[a]->location[0]);
				if( distClosestBorder>dist){
					border = borderCells[b];
					distClosestBorder=dist;
				}
			}
		}else{
			//borderCells[]
		}
			//border =
				//// SLOWEST!
				////voronoiDiagram->searchClosestUnoccupiedVoronoiCellNAIVE( agentArray->agents[a]->location[0]);
				// SLOW
				//voronoiDiagram->searchClosestUnoccupiedVoronoiCell(	agentArray->agents[a]->location[0], 100);
		//else
		//	border = voronoiDiagram->searchClosestFreeVoronoiCell( agentArray->agents[a]->location[0], 100);

			//fprintf(stderr, "\r[Agent %i]: Upd. Stats                 \b", a );
		double dist = 1000000;
		if (border) {
			dist = border->getDistanceTo(
					agentArray->agents[a]->location[0])
					* SPATIAL_UNIT;
		} else {
			//fprintf(stderr, "No border found!\n");
			// shortest distance to domain border
			for( int i=0; i<DIMENSIONS; i++){
				dist = fmin( dist, agentArray->agents[a]->location[0]->position[i]);
				dist = fmin( dist, 100 /*pow(agentArray->countAgents, 1./DIMENSIONS)*/ - agentArray->agents[a]->location[0]->position[i]);
			}
			dist *= SPATIAL_UNIT;
		}

		for (int i = (dist - AGENT_DIAMETER > 0. ? (int) (dist - AGENT_DIAMETER) : 0);
				i < (int) (dist); i++) {
			histogram[i]++;
			if (dividing)
				histogramDividing[i]++;
			if (necrotic)
				histogramNecrotic[i]++;
			if (free)
				histogramFree[i]++;
			histogramECM[i] +=
					agentArray->agents[a]->location[0]->ecm;
		}

	}
	//fprintf(stderr, "[finished]\n");

}

double GetDistanceToClosestFreeNeighbor(VoronoiDiagram *voronoiDiagram,
		VoronoiCell *voronoiCell) {
	int dist = 0;
	int min_dist = 0;
	for (int d = 0; d < DIMENSIONS; d++)
		min_dist += voronoiDiagram->xN[d];

	if (VoronoiCell::USE_SYMBOLIC_EXTENDED_NEIGHBORHOOD) {
		VoronoiCell *closestVC;
		if (VoronoiCell::SHIFT_TO_UNOCCUPIED)
			closestVC =
					voronoiDiagram->searchClosestUnoccupiedVoronoiCell(
							voronoiCell, VoronoiCell::extendedNeighborDepth,
							VoronoiCell::symbolicExtendedNeighborhoodSize,
							VoronoiCell::symbolicExtendedNeighborhood);
		else
			closestVC =
					voronoiDiagram->searchClosestFreeVoronoiCell(voronoiCell,
							VoronoiCell::extendedNeighborDepth,
							VoronoiCell::symbolicExtendedNeighborhoodSize,
							VoronoiCell::symbolicExtendedNeighborhood);
		if (closestVC)
			min_dist = closestVC->getDistanceSquareTo(voronoiCell);
	} else if (VoronoiCell::USE_DYNAMIC_EXTENDED_NEIGHBORHOOD) {
		VoronoiCell *closestVC;
		if (VoronoiCell::SHIFT_TO_UNOCCUPIED)
			closestVC =
					voronoiDiagram->searchClosestUnoccupiedVoronoiCell(
							voronoiCell, VoronoiCell::extendedNeighborDepth);
		else
			closestVC =
					voronoiDiagram->searchClosestFreeVoronoiCell(voronoiCell,
							VoronoiCell::extendedNeighborDepth);

		if (closestVC)
			min_dist = closestVC->getDistanceSquareTo(voronoiCell);
	} else {
		if (voronoiCell->countFreeNeighborCells)
			for (int i = 0; i < voronoiCell->countNeighborCells; i++) {
				dist = 0;
				for (int d = 0; d < DIMENSIONS; d++)
					dist +=
							pow(
									voronoiCell->position[d]
											- voronoiCell->neighborCells[i]->position[d],
									2);
				if (min_dist > dist)
					min_dist = dist;
			}

		if (voronoiCell->countFreeExtendedNeighborCells)
			for (int i = 0; i < voronoiCell->countExtendedNeighborCells; i++) {
				dist = 0;
				for (int d = 0; d < DIMENSIONS; d++)
					dist +=
							pow(
									voronoiCell->position[d]
											- voronoiCell->extendedNeighborhood[i]->position[d],
									2);
				if (min_dist > dist)
					min_dist = dist;
			}
	}

	return sqrt(min_dist);
}

void refineSurrounding(VoronoiDiagram *voronoiDiagram, VoronoiCell *cell,
		ActionTree *actionList, AgentList *agentArray, int scale) {
	int index; // = floor(cell->position[0])*voronoiDiagram->xN[1] + floor(cell->position[1]);
	int x = (int) floor(cell->position[0]);
	int y = (int) floor(cell->position[1]);
#if DIMENSIONS == 3
	int z = (int)floor(cell->position[2]);
	//int dz = 1;
	int dy = voronoiDiagram->xN[2];
	int dx = voronoiDiagram->xN[1] * dy;
	int minz = MAX(0,z-1);
	int maxz = MIN(voronoiDiagram->xN[2]-1,z+1);
#elif DIMENSIONS == 2
	int dy = 1;
	int dx = voronoiDiagram->xN[1];
#endif
	int n = 0;
	//for( int ix=MAX(0,x-1); ix<MIN(voronoiDiagram->xN[0],x+2); ix++)
	//for( int iy=MAX(0,y-1); iy<MIN(voronoiDiagram->xN[1],y+2); iy++)
	int minx = MAX(0,x-1);
	int miny = MAX(0,y-1);
	int maxx = MIN(voronoiDiagram->xN[0]-1,x+1);
	int maxy = MIN(voronoiDiagram->xN[1]-1,y+1);

	for (int ix = minx; ix <= maxx; ix++)
		for (int iy = miny; iy <= maxy; iy++)
#if DIMENSIONS == 2
				{
			index = ix * dx + iy;
#elif DIMENSIONS == 3
			for( int iz=minz; iz<=maxz; iz++)
			//for( int iz=MAX(0,z-1); iz<MIN(voronoiDiagram->xN[2],z+2); iz++)
			{
				index = ix*dx +iy*dy +iz;
#endif
			//{
			//index = ix*dx +iy*dy;
			if (voronoiDiagram->voronoiCells[index]->refined == false
					&& voronoiDiagram->voronoiCells[index]->isFree()) {

				//long passedTime = clock();
				VoronoiCell* comp = voronoiDiagram->voronoiCells[index];
				voronoiDiagram->refine(voronoiDiagram->voronoiCells[index],
						scale, actionList);
				//fprintf(stderr, "...finished ( %.3lfsec)\n", (double)(clock() - passedTime)/(double)CLOCKS_PER_SEC);

				/*int n=0;
				 for( ; n<GetVoronoiCell( selected_action->originalCell)->countNeighborCells; )
				 {
				 if( GetVoronoiCell( selected_action->originalCell)->neighborCells[n]->index >= 0 &&
				 GetVoronoiCell( selected_action->originalCell)->neighborCells[n]->refined == false
				 && GetVoronoiCell( selected_action->originalCell)->neighborCells[n]->isFree()){
				 //printTriangulation( voronoiDiagram, "beforeRef.eps", 1);
				 passedTime = clock();
				 VoronoiCell* comp = GetVoronoiCell( selected_action->originalCell)->neighborCells[n];
				 voronoiDiagram->refine( GetVoronoiCell( selected_action->originalCell)->neighborCells[n], CountCellsPerVoronoiCell, actionList);
				 fprintf(stderr, "...finished ( %.3lfsec)\n", (double)(clock() - passedTime)/(double)CLOCKS_PER_SEC);
				 //printTriangulation( voronoiDiagram, "afterRef.eps", 1);
				 //voronoiDiagram->coarsen( GetVoronoiCell( selected_action->originalCell)->neighborCells[n], 10);
				 n=0;*/

				//TEST
				if (comp->agent == NULL)
				//exit(0);
				{
					//fprintf(stderr, "Add Compartment Agent\n");

					Agent *compAgent = agentArray->activateAgent();
					compAgent->attach(comp);
					compAgent->state = COMPARTMENT; //NONACTIVE;
					compAgent->cellCount = 0;
					compAgent->maxCellCount = scale * scale
#if DIMENSIONS == 3
					*scale
#endif
					;
					compAgent->cellCount = compAgent->maxCellCount;
					compAgent->countFree = compAgent->maxCellCount;
					compAgent->countActive = 0;
					compAgent->countNonactive = 0;
					//for( int n=0; n<comp->countNeighborCells; n++)
					//	comp->neighborCells[n]->countFreeNeighborCells--;
				}
				//TEST END
			} else {
				n++;
			}
#if DIMENSIONS >= 2
		}
#endif
}

void refineNeighborhood(VoronoiDiagram *voronoiDiagram, VoronoiCell *cell,
		ActionTree *actionList, AgentList *agentArray, int scale) {
	/*		int index;// = floor(cell->position[0])*voronoiDiagram->xN[1] + floor(cell->position[1]);
	 int x = (int)floor(cell->position[0]);
	 int y = (int)floor(cell->position[1]);
	 #if DIMENSIONS == 3
	 int z = (int)floor(cell->position[2]);
	 //int dz = 1;
	 int dy = voronoiDiagram->xN[2];
	 int dx = voronoiDiagram->xN[1] * dy;
	 int minz = MAX(0,z-1);
	 int maxz = MIN(voronoiDiagram->xN[2]-1,z+1);
	 #elif DIMENSIONS == 2
	 int dy = 1;
	 int dx = voronoiDiagram->xN[1];
	 #endif
	 int n=0;
	 //for( int ix=MAX(0,x-1); ix<MIN(voronoiDiagram->xN[0],x+2); ix++)
	 //for( int iy=MAX(0,y-1); iy<MIN(voronoiDiagram->xN[1],y+2); iy++)
	 int minx = MAX(0,x-1);
	 int miny = MAX(0,y-1);
	 int maxx = MIN(voronoiDiagram->xN[0]-1,x+1);
	 int maxy = MIN(voronoiDiagram->xN[1]-1,y+1);

	 for( int ix=minx; ix<=maxx; ix++)
	 for( int iy=miny; iy<=maxy; iy++)
	 #if DIMENSIONS == 2
	 {
	 index = ix*dx +iy;
	 #elif DIMENSIONS == 3
	 for( int iz=minz; iz<=maxz; iz++)
	 //for( int iz=MAX(0,z-1); iz<MIN(voronoiDiagram->xN[2],z+2); iz++)
	 {
	 index = ix*dx +iy*dy +iz;
	 #endif
	 //{
	 //index = ix*dx +iy*dy;
	 if( voronoiDiagram->voronoiCells[index]->refined == false
	 && voronoiDiagram->voronoiCells[index]->isFree()){

	 long passedTime = clock();
	 VoronoiCell* comp = voronoiDiagram->voronoiCells[index];
	 voronoiDiagram->refine( voronoiDiagram->voronoiCells[index], scale, actionList);
	 fprintf(stderr, "...finished ( %.3lfsec)\n", (double)(clock() - passedTime)/(double)CLOCKS_PER_SEC);
	 */
	//for( int n=0 ; n<cell->countNeighborCells; n++)
	//	fprintf( stderr, "%i, ", cell->neighborCells[n]->index);
	//fprintf( stderr, "\n");
	int n = 0;
	for (; n < cell->countNeighborCells;) {
		if (cell->neighborCells[n]->index >= 0
				&& cell->neighborCells[n]->refined == false
				&& cell->neighborCells[n]->isFree()) {
			//printTriangulation( voronoiDiagram, "beforeRef.eps", 1);
			//long passedTime = clock();
			VoronoiCell* comp = cell->neighborCells[n];
			voronoiDiagram->refine(cell->neighborCells[n],
					CountCellsPerVoronoiCell, actionList);
			//fprintf(stderr, "...finished ( %.3lfsec)\n", (double)(clock() - passedTime)/(double)CLOCKS_PER_SEC);
			//printTriangulation( voronoiDiagram, "afterRef.eps", 1);
			//voronoiDiagram->coarsen( GetVoronoiCell( selected_action->originalCell)->neighborCells[n], 10);
			n = 0;

			//TEST
			if (comp->agent == NULL)
			//exit(0);
			{
				//fprintf(stderr, "Add Compartment Agent\n");

				Agent *compAgent = agentArray->activateAgent();
				compAgent->attach(comp);
				compAgent->state = COMPARTMENT; //NONACTIVE;
				compAgent->cellCount = 0;
				compAgent->maxCellCount = scale * scale
#if DIMENSIONS == 3
				*scale
#endif
				;
				compAgent->cellCount = compAgent->maxCellCount;
				compAgent->countFree = compAgent->maxCellCount;
				compAgent->countActive = 0;
				compAgent->countNonactive = 0;
				//for( int n=0; n<comp->countNeighborCells; n++)
				//	comp->neighborCells[n]->countFreeNeighborCells--;
			}
			//TEST END
		} else {
			n++;
		}
	}
}

void printTriangulation(VoronoiDiagram* voronoiDiagram, const char * filename,
		bool plotIndex) {
	std::fstream fs;
	fs.open(filename, std::fstream::out);
	EpsIO::PSwriteHeader(&fs, voronoiDiagram->xMin[0] * 10,
			voronoiDiagram->xMax[0] * 10, voronoiDiagram->xMin[1] * 10,
			voronoiDiagram->xMax[1] * 10);
	for (int c = 0; c < voronoiDiagram->countVoronoiCells; c++) {
		for (int n = 0; n < voronoiDiagram->voronoiCells[c]->countNeighborCells;
				n++) {
			EpsIO::PSwriteLine(
					&fs,
					voronoiDiagram->voronoiCells[c]->position[0] * 10,
					voronoiDiagram->voronoiCells[c]->position[1] * 10,
					//voronoiDiagram->voronoiCells[c]->neighborCells[n]->position[0]*10,
					//voronoiDiagram->voronoiCells[c]->neighborCells[n]->position[1]*10,
					(voronoiDiagram->voronoiCells[c]->position[0]
							+ voronoiDiagram->voronoiCells[c]->neighborCells[n]->position[0])
							/ 2. * 10,
					(voronoiDiagram->voronoiCells[c]->position[1]
							+ voronoiDiagram->voronoiCells[c]->neighborCells[n]->position[1])
							/ 2. * 10, 0.1, (char*) " 0 0 0 setrgbcolor");
			/*if( voronoiDiagram->voronoiCells[c]==centralCell)
			 EpsIO::PSwriteCircle( &fs,
			 voronoiDiagram->voronoiCells[c]->position[0]*10,
			 voronoiDiagram->voronoiCells[c]->position[1]*10,
			 1,
			 1,
			 " 0 0 1 setrgbcolor");
			 if( c==voronoiDiagram->countVoronoiCells-1)
			 EpsIO::PSwriteCircle( &fs,
			 voronoiDiagram->voronoiCells[c]->position[0]*10,
			 voronoiDiagram->voronoiCells[c]->position[1]*10,
			 1,
			 1,
			 " 1 0 0 setrgbcolor");*/
		}
		// /Times-Roman 270 selectfont
		if (plotIndex && voronoiDiagram->voronoiCells[c]->countNeighborCells > 0) {
			fs << "/Times-Roman 4 selectfont" << endl;
			fs << "0.5 setgray "
					<< voronoiDiagram->voronoiCells[c]->position[0] * 10 << " "
					<< voronoiDiagram->voronoiCells[c]->position[1] * 10
					<< " moveto (" << voronoiDiagram->voronoiCells[c]->index
					<< ") show" << endl;
		}
	}
	fs.close();

}

void printVoronoiDiagram(VoronoiDiagram* voronoiDiagram, const char * filename,
		bool plotIndex) {
	std::fstream fs;
	fs.open(filename, std::fstream::out);
	EpsIO::PSwriteHeader(&fs, voronoiDiagram->xMin[0] * 10,
			voronoiDiagram->xMax[0] * 10, voronoiDiagram->xMin[1] * 10,
			voronoiDiagram->xMax[1] * 10);
	for (int t = 0; t < voronoiDiagram->countTetrahedra; t++) {
		double vp1[2], vp2[2];
		getCircumCenter(voronoiDiagram->tetrahedra[t], vp1);
		for (int n = 0;
				n < voronoiDiagram->tetrahedra[t]->countNeighborTetrahedra;
				n++) {
			getCircumCenter(
					voronoiDiagram->tetrahedra[t]->neighborTetrahedra[n], vp2);
			EpsIO::PSwriteLine(&fs, vp1[0] * 10,
					vp1[1] * 10,
					//voronoiDiagram->voronoiCells[c]->neighborCells[n]->position[0]*10,
					//voronoiDiagram->voronoiCells[c]->neighborCells[n]->position[1]*10,
					(vp1[0] + vp2[0]) / 2. * 10, (vp1[1] + vp2[1]) / 2. * 10,
					0.1, (char*) (" 0 0 0 setrgbcolor"));
			/*if( voronoiDiagram->voronoiCells[c]==centralCell)
			 EpsIO::PSwriteCircle( &fs,
			 voronoiDiagram->voronoiCells[c]->position[0]*10,
			 voronoiDiagram->voronoiCells[c]->position[1]*10,
			 1,
			 1,
			 " 0 0 1 setrgbcolor");
			 if( c==voronoiDiagram->countVoronoiCells-1)
			 EpsIO::PSwriteCircle( &fs,
			 voronoiDiagram->voronoiCells[c]->position[0]*10,
			 voronoiDiagram->voronoiCells[c]->position[1]*10,
			 1,
			 1,
			 " 1 0 0 setrgbcolor");*/
		}
		// /Times-Roman 270 selectfont
		/*if( plotIndex && voronoiDiagram->voronoiCells[c]->countNeighborCells>0){
		 fs << "/Times-Roman 4 selectfont" << endl;
		 fs << "0.5 setgray " << voronoiDiagram->voronoiCells[c]->position[0]*10 << " " << voronoiDiagram->voronoiCells[c]->position[1]*10 << " moveto ("<< voronoiDiagram->voronoiCells[c]->index <<") show" << endl;
		 }*/
	}

	for (int c = 0; c < voronoiDiagram->countVoronoiCells; c++) {
		if (voronoiDiagram->voronoiCells[c]->countNeighborCells > 0) {
			char color[200];
			float radius;

			if (voronoiDiagram->voronoiCells[c]->refined) {
				radius = 0.5;
			} else {
				radius = 2.;
			}

			switch (voronoiDiagram->voronoiCells[c]->getState()) {
			case FREE:
				sprintf(color, " 0 0 0 setrgbcolor");
				break;
			case ACTIVE:
				sprintf(color, " 1 0 0 setrgbcolor");
				break;
			case NONACTIVE:
				sprintf(color, " 0 1 0 setrgbcolor");
				break;
			case COMPARTMENT:
				if (GetAgent(voronoiDiagram->voronoiCells[c])->countFree
						== GetAgent(voronoiDiagram->voronoiCells[c])->maxCellCount)
					sprintf(color, " 0 0 0 setrgbcolor");
				else if (GetAgent(voronoiDiagram->voronoiCells[c])->countActive
						== GetAgent(voronoiDiagram->voronoiCells[c])->maxCellCount)
					sprintf(color, " 1 0 0 setrgbcolor");
				else if (GetAgent(voronoiDiagram->voronoiCells[c])->countNonactive
						== GetAgent(voronoiDiagram->voronoiCells[c])->maxCellCount)
					sprintf(color, " 0 1 0 setrgbcolor");
				else
					sprintf(color, " 1 0 1 setrgbcolor");
				radius = 2.;
				break;
			}

			//if(voronoiDiagram->voronoiCells[c]->refined)
			if (voronoiDiagram->voronoiCells[c]->getState() != FREE)
				EpsIO::PSwriteCircle(&fs,
						voronoiDiagram->voronoiCells[c]->position[0] * 10,
						voronoiDiagram->voronoiCells[c]->position[1] * 10,
						radius, 0., color); //" 1 0 0 setrgbcolor");
			/*else
			 EpsIO::PSwriteCircle( &fs,
			 voronoiDiagram->voronoiCells[c]->position[0]*10,
			 voronoiDiagram->voronoiCells[c]->position[1]*10,
			 2,
			 2,
			 " 0 0 1 setrgbcolor");
			 */
		}
	}

	fs.close();

}

void performMigration(VoronoiDiagram *voronoiDiagram, AgentList *agentArray,
		ActionTree *actionList, Action *selected_action,
		Statistix *myStatistics, double &Time, double &EndTime);

void performChemotacticMigration(VoronoiDiagram *voronoiDiagram,
		AgentList *agentArray, ActionTree *actionList, Action *selected_action,
		Statistix *myStatistics, double &Time, double &EndTime);

VoronoiCell* performGrowth(VoronoiDiagram *voronoiDiagram,
		AgentList *agentArray, ActionTree *actionList, Action *selected_action,
		Statistix *myStatistics, double &Time, double &EndTime);
void performDivision(VoronoiDiagram *voronoiDiagram, AgentList *agentArray,
		ActionTree *actionList, Action *selected_action,
		Statistix *myStatistics, double &Time, double &EndTime);
VoronoiCell* performGrowthAndDivision(VoronoiDiagram *voronoiDiagram,
		AgentList *agentArray, ActionTree *actionList, Action *selected_action,
		Statistix *myStatistics, double &Time, double &EndTime);
void performFreeMigration(VoronoiDiagram *voronoiDiagram, AgentList *agentArray,
		ActionTree *actionList, Action *selected_action,
		Statistix *myStatistics, double &Time, double &EndTime);

double getAvgGlucose(AgentList *agentList) {
	double avg = 0.;
	double N = 0;

	for (int i = 0; i < agentList->countActiveAgents; i++) {
		if (agentList->agents[i]->growingTumorCellCount > 0){
			for (int ii = 0; ii < agentList->agents[i]->countLocations; ii++){
				if (agentList->agents[i]->state == COMPARTMENT) {
					avg +=
							agentList->agents[i]->location[ii]->glucose
									* agentList->agents[i]->growingTumorCellCount;
					N += agentList->agents[i]->growingTumorCellCount;
				} else {
					avg += agentList->agents[i]->location[ii]->glucose;
					N++;
				}
			}
		}
	}

	return (N > 0. ? avg / N : 0.);
}

double getAvgOxygen(AgentList *agentList) {
	double avg = 0.;
	double N = 0;

	for (int i = 0; i < agentList->countActiveAgents; i++) {
		if (agentList->agents[i]->growingTumorCellCount > 0){
			for (int ii = 0; ii < agentList->agents[i]->countLocations; ii++){
				if (agentList->agents[i]->state == COMPARTMENT) {
					if (N == 0
							|| avg
									> agentList->agents[i]->location[ii]->glucose
											* agentList->agents[i]->location[ii]->oxygen) {
						avg =
								agentList->agents[i]->location[ii]->glucose
										* agentList->agents[i]->location[ii]->oxygen;
						//			if( N==0 || avg > agentList->agents[i]->location[ii]->glucose){
						//				avg = agentList->agents[i]->location[ii]->glucose;
						N = 1;
					}
					//			avg += agentList->agents[i]->location[ii]->oxygen * agentList->agents[i]->growingTumorCellCount;
					//			N += agentList->agents[i]->growingTumorCellCount;
				} else {
					if (N == 0
							|| avg
									> agentList->agents[i]->location[ii]->glucose
											* agentList->agents[i]->location[ii]->oxygen) {
						avg =
								agentList->agents[i]->location[ii]->glucose
										* agentList->agents[i]->location[ii]->oxygen;
						//			if( N==0 || avg > agentList->agents[i]->location[ii]->glucose){
						//				avg = agentList->agents[i]->location[ii]->glucose;
						N = 1;
					}
					//			avg += agentList->agents[i]->location[ii]->oxygen;
					//			N ++;
				}
			}
		}
	}

	return (N > 0. ? avg / N : 0.);
}

VoronoiCell **getOuterBorder(AgentList *agentList, int &borderSize) {
	int max_queueSize = agentList->countActiveAgents
			* (MAX_SUBCELLULAR_COMPONENTS) * 2;
	int queueSize = 0;
	queueSize = 0;
	VoronoiCell **queue = (VoronoiCell**) calloc(max_queueSize,
			sizeof(VoronoiCell*));

	//fprintf( stderr, "ENTER!\n");

	//fprintf( stderr, "agentList->countActiveAgents: %i\n", agentList->countActiveAgents);

	// find first border location
	for (int i = 0; !queueSize && i < agentList->countActiveAgents; i++) {
		//	fprintf( stderr, "agent: %i\n", agentList->agents[i]->index);
		if (agentList->agents[i]->state != FREE)
			for (int ii = 0;
					!queueSize && ii < agentList->agents[i]->countLocations;
					ii++) {
				//		fprintf( stderr, "-->location: %i\n", agentList->agents[i]->location[ii]->index);
				for (int iii = 0;
						!queueSize
								&& iii
										< agentList->agents[i]->location[ii]->countNeighborCells;
						iii++) {
					//			fprintf( stderr, "---->neighbor: %i\n", agentList->agents[i]->location[ii]->neighborCells[iii]->index);
					if (GetAgent( agentList->agents[i]->location[ii]->neighborCells[iii])
							== NULL) {
						if (max_queueSize == queueSize) {
							fprintf(
									stderr,
									"WARNING: realloc queue: %i->%i\n",
									max_queueSize,
									max_queueSize + agentList->countActiveAgents);
							max_queueSize += agentList->countActiveAgents;
							queue =
									(VoronoiCell**) realloc(
											queue,
											max_queueSize
													* sizeof(VoronoiCell *));
						}
						queue[queueSize++] =
								agentList->agents[i]->location[ii]->neighborCells[iii];
						//				fprintf( stderr, "    >ADD TO QUEUE!\n");
					}
				}
			}
	}
	//fprintf( stderr, "AHA?!\n");
	//fprintf( stderr, "queueSize: %i, queue[0]: %i\n", queueSize, queue[0]->index);
	if (queueSize == 0)
		return NULL;
	//fprintf( stderr, "AHA!!!!\n");

	// collect all border points
	int max_borderSize = agentList->countActiveAgents
			* (MAX_SUBCELLULAR_COMPONENTS) * 2;
	borderSize = 0;
	VoronoiCell **border = (VoronoiCell**) malloc(
			max_borderSize * sizeof(VoronoiCell*));
	do {
		// take next element from queue
		VoronoiCell *actual = queue[--queueSize];

		// test if containes empty neighbors
		//fprintf( stderr, "-->next: %i\n", actual->index);
		char occupiedNeighborFound = FALSE;
		for (int iii = 0;
				!occupiedNeighborFound && iii < actual->countNeighborCells;
				iii++) {
			//	fprintf( stderr, "---->neighbor: %i\n", actual->neighborCells[iii]->index);
			// test if containes empty neighbors
			if (GetAgent( actual->neighborCells[iii]) != NULL)
				occupiedNeighborFound = TRUE;
		}

		if (occupiedNeighborFound) {
			// add neighbors to queue
			for (int iii = 0; iii < actual->countNeighborCells; iii++) {
				if (GetAgent( actual->neighborCells[iii]) == NULL
						&& !isElementOf(border, borderSize,
								actual->neighborCells[iii])
						&& !isElementOf(queue, queueSize,
								actual->neighborCells[iii])) {

					if (max_queueSize == queueSize) {
						fprintf(stderr, "WARNING: realloc queue: %i->%i\n",
								max_queueSize,
								max_queueSize + agentList->countActiveAgents);
						max_queueSize += agentList->countActiveAgents;
						queue =
								(VoronoiCell**) realloc(queue,
										max_queueSize * sizeof(VoronoiCell *));
					}

					queue[queueSize++] = actual->neighborCells[iii];
					//	fprintf( stderr, "    >ADD TO QUEUE!\n");
				}
			}

			// add to border
			//if( max_borderSize==borderSize)
			//	exit( 0);
			if (max_borderSize == borderSize) {
				fprintf(stderr, "WARNING: realloc border: %i->%i\n",
						max_borderSize,
						max_borderSize + agentList->countActiveAgents);
				max_borderSize += agentList->countActiveAgents;
				border =
						(VoronoiCell**) realloc(border,
								max_borderSize * sizeof(VoronoiCell *));
			}
			border[borderSize++] = actual;
		}
	} while (queueSize > 0);

	free(queue);

	//fprintf( stderr, "border points:");		
	for (int i = 0; i < borderSize; i++) {
		//fprintf( stderr, " %i", border[i]->index);		
	}
	//fprintf( stderr, "\n");
	//fprintf( stderr, "LEAVE!\n");

	return border;

	free(border);
	borderSize = 0;
	return (VoronoiCell**) malloc(1 * sizeof(VoronoiCell*));
}

VoronoiCell **getOuterBorder2(AgentList *agentList, int &borderSize) {
	int max_borderSize = agentList->countActiveAgents
			* (MAX_SUBCELLULAR_COMPONENTS) * 2;
	borderSize = 0;
	VoronoiCell **border = (VoronoiCell**) malloc(
			max_borderSize * sizeof(VoronoiCell*));

	//fprintf( stderr, "ENTER!\n");

	//fprintf( stderr, "agentList->countActiveAgents: %i\n", agentList->countActiveAgents);

	// find first border location
	for (int i = 0; i < agentList->countActiveAgents; i++) {
		//	fprintf( stderr, "agent: %i\n", agentList->agents[i]->index);
		if (agentList->agents[i]->state == ACTIVE)
			for (int ii = 0; ii < agentList->agents[i]->countLocations; ii++) {
				//		fprintf( stderr, "-->location: %i\n", agentList->agents[i]->location[ii]->index);
				for (int iii = 0;
						iii
								< agentList->agents[i]->location[ii]->countNeighborCells;
						iii++) {
					//			fprintf( stderr, "---->neighbor: %i\n", agentList->agents[i]->location[ii]->neighborCells[iii]->index);
					if (GetAgent( agentList->agents[i]->location[ii]->neighborCells[iii])
							== NULL
							&& !isElementOf(
									border,
									borderSize,
									agentList->agents[i]->location[ii]->neighborCells[iii])) {
						if (max_borderSize == borderSize) {
							fprintf(
									stderr,
									"WARNING: realloc queue: %i->%i\n",
									max_borderSize,
									max_borderSize
											+ agentList->countActiveAgents);
							max_borderSize += agentList->countActiveAgents;
							border =
									(VoronoiCell**) realloc(
											border,
											max_borderSize
													* sizeof(VoronoiCell *));
						}
						border[borderSize++] =
								agentList->agents[i]->location[ii]->neighborCells[iii];
						//				fprintf( stderr, "    >ADD TO QUEUE!\n");
					}
				}
			}
	}
	//fprintf( stderr, "AHA?!\n");
	//fprintf( stderr, "queueSize: %i, queue[0]: %i\n", queueSize, queue[0]->index);
	//if( queueSize == 0)
	//	return NULL;
	//fprintf( stderr, "AHA!!!!\n");

	return border;
}

double getGyrationRadius(AgentList *agentList) {
	// center of mass
	//fprintf( stderr, "center of mass:");
	double centerMass[DIMENSIONS];
	int count=agentList->countActiveAgents;
	for (int i = 0; i < DIMENSIONS; i++) {
		centerMass[i] = 0.;
		count=0;
		for (int ii = 0; ii < agentList->countActiveAgents; ii++)
			if( agentList->agents[ii]->state != FREE){
			centerMass[i] += agentList->agents[ii]->location[0]->position[i];
			count++;
		}
		centerMass[i] /= (double) count;
		//	fprintf( stderr, " %lf", centerMass[i]);
	}
	//fprintf( stderr, "\n");

	// radius of gyration
	double gyrRadius = 0.;
	for (int ii = 0; ii < agentList->countActiveAgents; ii++)
		if( agentList->agents[ii]->state != FREE){
		for (int i = 0; i < DIMENSIONS; i++) {
			gyrRadius += pow(agentList->agents[ii]->location[0]->position[i] - centerMass[i], 2.);
		}
	}
	gyrRadius /= (double) count;
	//fprintf( stderr, "radius of gyration: %lf\n", gyrRadius);

	return gyrRadius;
}

double getGyrationRadiusOfBorder(AgentList *agentList, VoronoiDiagram *vd) {

	//fprintf( stderr, "Gyr Rad\n");
	/*int max_queueSize = agentList->countActiveAgents*MAX_SUBCELLULAR_COMPONENTS;
	 int queueSize = 0;
	 VoronoiCell **queue = (VoronoiCell**) calloc( max_queueSize, sizeof(VoronoiCell*));

	 //fprintf( stderr, "agentList->countActiveAgents: %i\n", agentList->countActiveAgents);

	 // find first border location
	 for( int i=0; !queueSize && i<agentList->countActiveAgents; i++){
	 //	fprintf( stderr, "agent: %i\n", agentList->agents[i]->index);
	 for( int ii=0; !queueSize && ii<agentList->agents[i]->countLocations; ii++){
	 //		fprintf( stderr, "-->location: %i\n", agentList->agents[i]->location[ii]->index);
	 for( int iii=0; !queueSize && iii<agentList->agents[i]->location[ii]->countNeighborCells; iii++){
	 //			fprintf( stderr, "---->neighbor: %i\n", agentList->agents[i]->location[ii]->neighborCells[iii]->index);
	 if( GetAgent( agentList->agents[i]->location[ii]->neighborCells[iii]) == NULL){
	 if( max_queueSize==queueSize){
	 queue = (VoronoiCell**) realloc( queue, (max_queueSize + agentList->countActiveAgents) * sizeof(VoronoiCell *));
	 }
	 queue[queueSize++] = agentList->agents[i]->location[ii];
	 //				fprintf( stderr, "    >ADD TO QUEUE!\n");
	 }
	 }
	 }
	 }
	 //fprintf( stderr, "queueSize: %i, queue[0]: %i\n", queueSize, queue[0]->index);
	 if( queueSize == 0)
	 return 0.;


	 // collect all border points
	 int max_borderSize = agentList->countActiveAgents*MAX_SUBCELLULAR_COMPONENTS;
	 int borderSize = 0;
	 VoronoiCell *border[max_borderSize];
	 do{
	 // take next element from queue
	 VoronoiCell *actual = queue[--queueSize];

	 // test if containes empty neighbors
	 //	fprintf( stderr, "-->next: %i\n", actual->index);
	 char emptyNeighborFound = FALSE;
	 for( int iii=0; !emptyNeighborFound && iii<actual->countNeighborCells; iii++){
	 //		fprintf( stderr, "---->neighbor: %i\n", actual->neighborCells[iii]->index);
	 // test if containes empty neighbors
	 if( GetAgent( actual->neighborCells[iii]) == NULL)
	 emptyNeighborFound = TRUE;
	 }

	 if( emptyNeighborFound){
	 // add neighbors to queue
	 for( int iii=0; iii<actual->countNeighborCells; iii++){
	 if( GetAgent( actual->neighborCells[iii]) != NULL &&
	 !isElementOf( border, borderSize, actual->neighborCells[iii]) &&
	 !isElementOf( queue,  queueSize, actual->neighborCells[iii])){
	 queue[queueSize++] = actual->neighborCells[iii];
	 //				fprintf( stderr, "    >ADD TO QUEUE!\n");
	 }
	 }

	 // add to border
	 if( max_borderSize==borderSize)
	 exit( 0);
	 border[borderSize++] = actual;
	 }
	 }while(queueSize>0);

	 free( queue);

	 //fprintf( stderr, "border points:");
	 for( int i=0; i<borderSize; i++){
	 //	fprintf( stderr, " %i", border[i]->index);
	 }
	 //fprintf( stderr, "\n");
	 */

	// collect border points
	int max_borderSize = agentList->countActiveAgents * MAX_SUBCELLULAR_COMPONENTS;
	int borderSize = 0;
	//VoronoiCell *border[max_borderSize];
	VoronoiCell **border = (VoronoiCell**) malloc(sizeof(VoronoiCell*) * max_borderSize);
	//fprintf(stderr, "Call\n");
	for (int i = 0; i < agentList->countActiveAgents; i++) {
		//	fprintf( stderr, "agent: %i\n", agentList->agents[i]->index);
		if (agentList->agents[i]->state != VESSEL)
			for (int ii = 0; ii < agentList->agents[i]->countLocations; ii++) {
				//		fprintf( stderr, "-->location: %i\n", agentList->agents[i]->location[ii]->index);

				char foundEmptyNeighbor = FALSE;
				// check for border
				if( agentList->agents[i]->location[ii]->isDomainBorder(vd)){
					if (borderSize == max_borderSize) {
						max_borderSize += agentList->countActiveAgents;
						border = (VoronoiCell**) realloc( border, max_borderSize * sizeof(VoronoiCell *));
					}
					//fprintf(stderr, "Found\n");
					border[borderSize++] = agentList->agents[i]->location[ii];

					foundEmptyNeighbor = TRUE;
				}

				// check for empty neighbors
				for (int iii = 0; !foundEmptyNeighbor && iii < agentList->agents[i]->location[ii]->countNeighborCells; iii++) {
					//			fprintf( stderr, "---->neighbor: %i\n", agentList->agents[i]->location[ii]->neighborCells[iii]->index);
					if ( GetAgent( agentList->agents[i]->location[ii]->neighborCells[iii]) == NULL) {
						if (borderSize == max_borderSize) {
							max_borderSize += agentList->countActiveAgents;
							border = (VoronoiCell**) realloc( border, max_borderSize * sizeof(VoronoiCell *));
						}

						border[borderSize++] = agentList->agents[i]->location[ii];
						foundEmptyNeighbor = TRUE;
						//				fprintf( stderr, "    >ADD TO QUEUE!\n");
						if (agentList->agents[i]->location[ii] == NULL) {
							fprintf(
									stderr,
									"agentList->agents[%i]->location[%i / %i] == NULL\n",
									i, ii,
									agentList->agents[i]->countLocations);
							exit(0);
						}
						if (agentList->agents[i]->location[ii]->agent == NULL) {
							fprintf(
									stderr,
									"agentList->agents[%i]->location[%i]->agent == NULL\n",
									i, ii);
							exit(0);
						}
						if (agentList->agents[i]->location[ii]->agent
								!= agentList->agents[i]) {
							fprintf(
									stderr,
									"agentList->agents[%i]->location[%i](%i)->agent(%i) != agentList->agents[%i](%i)\n",
									i,
									ii,
									agentList->agents[i]->location[ii]->index,
									GetAgent(agentList->agents[i]->location[ii])->index,
									i, agentList->agents[i]->index);
							exit(0);
						}
						/*if( GetAgent( border[ii]) == NULL){
						 fprintf( stderr, "border element %i == NULL => borderSize=%i\n", ii, borderSize);
						 exit(0);
						 }*/
					}
				}
			}
	}

	if (borderSize > max_borderSize) {
		fprintf(stderr, "border size to big: %i > %i\n", borderSize,
				max_borderSize);
		exit(0);
	}

	// center of mass
	//fprintf( stderr, "center of mass:");		
	double centerMass[DIMENSIONS];
	for (int i = 0; i < DIMENSIONS; i++) {
		centerMass[i] = 0.;
		for (int ii = 0; ii < borderSize; ii++) {
			if (GetAgent( border[ii]) == NULL) {
				fprintf(stderr, "border element %i == NULL => borderSize=%i\n",
						ii, borderSize);
				exit(0);
			}
			//fprintf( stderr, "border %i == (%lf,%lf,%lf)\n", ii, border[ii]->position[0], border[ii]->position[1], border[ii]->position[2]);
			centerMass[i] += border[ii]->position[i];
		}
		centerMass[i] /= (double) borderSize;
		//	fprintf( stderr, " %lf", centerMass[i]);
	}
	//fprintf( stderr, "\n");		

	// radius of gyration
	double gyrRadius = 0.;
	for (int ii = 0; ii < borderSize; ii++) {
		for (int i = 0; i < DIMENSIONS; i++) {
			gyrRadius += pow(border[ii]->position[i] - centerMass[i], 2.);
		}
	}
	gyrRadius /= (double) borderSize;
	//fprintf( stderr, "radius of gyration: %lf\n", gyrRadius);		

	free(border);

	//fprintf( stderr, "Gyr Rad End\n");


	return gyrRadius;
	return 0.;
}

Action *selectActionFixTimeStep(AgentList *agentArray, double &Time) {

	int i;
	// get oldest generation of dividing cells
	Agent *growingAgents[2 * (int) CellDivisionDepth];
	int countGrowingAgents = 0;
	int oldestGeneration = -1;
	int youngestGeneration = -1;
	for (i = 0; i < agentArray->countActiveAgents; i++) {
		if (agentArray->agents[i]->actions[INDEX_GROWTH]->top != NULL) {
			if (youngestGeneration < 0
					|| agentArray->agents[i]->generation > youngestGeneration) {
				youngestGeneration = agentArray->agents[i]->generation;
			}
			if (oldestGeneration < 0
					|| agentArray->agents[i]->generation < oldestGeneration) {
				oldestGeneration = agentArray->agents[i]->generation;
				countGrowingAgents = 0;
				//fprintf( stderr, "*");
			}
			//fprintf( stderr, "%i ", agentArray->agents[i]->generation);
			if (agentArray->agents[i]->generation == oldestGeneration)
				growingAgents[countGrowingAgents++] = agentArray->agents[i];
		}
	}
	//fprintf( stderr, " (%i x %i, %i)\n", countGrowingAgents, oldestGeneration, youngestGeneration);

	// chose random cell to divide
	int random = (int) (myRand() * countGrowingAgents);
	Action *selected_action = growingAgents[random]->actions[INDEX_GROWTH];

	// update time
	if (oldestGeneration == youngestGeneration)
		Time += 1. / selected_action->rate;

	return selected_action;
}

Action *selectActionVariablTimeStep(AgentList *agentArray, double &Time) {

	int i;
	// get oldest generation of dividing cells
	Agent *growingAgents[2 * (int) CellDivisionDepth];
	int countGrowingAgents = 0;
	int oldestGeneration = -1;
	int youngestGeneration = -1;
	for (i = 0; i < agentArray->countActiveAgents; i++) {
		if (agentArray->agents[i]->actions[INDEX_GROWTH]->top != NULL) {
			if (youngestGeneration < 0
					|| agentArray->agents[i]->generation > youngestGeneration) {
				youngestGeneration = agentArray->agents[i]->generation;
			}
			if (oldestGeneration < 0
					|| agentArray->agents[i]->generation < oldestGeneration) {
				oldestGeneration = agentArray->agents[i]->generation;
				//countGrowingAgents = 0;
				//fprintf( stderr, "*");
			}
			//fprintf( stderr, "%i ", agentArray->agents[i]->generation);
			//if( agentArray->agents[i]->generation==oldestGeneration)
			growingAgents[countGrowingAgents++] = agentArray->agents[i];
		}
	}
	//fprintf( stderr, " (%i x %i, %i)\n", countGrowingAgents, oldestGeneration, youngestGeneration);

	// sum rates
	double rateSum = 0.;
	for (i = 0; i < countGrowingAgents; i++)
		rateSum += growingAgents[i]->actions[INDEX_GROWTH]->rate;

	// chose random cell to divide
	double randomRateSum = myRand() * rateSum;
	double tempRateSum = 0.;
	//Action *selected_action = growingAgents[random]->actions[INDEX_GROWTH];
	for (i = 0; i < countGrowingAgents; i++) {
		tempRateSum += growingAgents[i]->actions[INDEX_GROWTH]->rate;
		if (tempRateSum >= randomRateSum) {

			// output
			//int pos;
			//if( growingAgents[i]->actions[INDEX_GROWTH]->originalCell->freeNeighborCells == 0)
			//	pos =
			//printf();

			// update time
			Time += 1. / rateSum;
			return growingAgents[i]->actions[INDEX_GROWTH];
		}
	}

	return NULL;
}

/*double getRadius( ActionList *actionList)
 {
 double center[DIMENSIONS];
 int i, ii, dim;

 // center
 Action *action = actionList->head;
 int countInvolvedVoronoiCells = 0;
 for( i=0; i<actionList->length; i++){
 for( ii=0; ii<action->originalCell->countLocations; ii++){
 for( dim=0; dim<DIMENSIONS; dim++){
 center[dim] += action->originalCell->location[ii]->position[dim];
 }
 countInvolvedVoronoiCells++;
 }
 action = action->next;
 }
 for( dim=0; dim<DIMENSIONS; dim++){
 center[dim] /= countInvolvedVoronoiCells;
 }

 // find first border cell
 int maxStackLength = actionList->length;
 int stackLength = 0;
 VoronoiCell **stack = (VoronoiCell**) calloc( maxStackLength, sizeof(VoronoiCell*));
 for( i=0; i<5; i++){
 action = actionList->head;
 for( ii=0; ii<myRand()*actionList->length; ii++)
 action = action->next;
 for( ii=0; ii<action->originalCell->countLocations; ii++){
 stack[stackLength++] = action->originalCell->location[ii];
 action
 action = actionList->
 }
 }*/

void writeVoronoiCellInformation(VoronoiCell *vc, FILE *fp) {
	double cATP = GiveMeTheATP(vc);
	double alpha =
			MIN( vc->glucose * GiveMeTheGlucoseRate(vc), vc->oxygen * GiveMeTheOxygenRate(vc) / 6.)
					/ vc->glucose / GiveMeTheGlucoseRate(vc);
	if (GiveMeTheGlucoseRate(vc) == 0.)
		alpha = 0.;

	fprintf(
			fp,
			"%lf %lf %lf %lf %lf %lf %lf %lf %lf %i %i %i %i %e %lf %lf %lf %lf %lf %e %e %e\n",

			// 2D-coordinates
			vc->position[0], vc->position[1], vc->position[DIMENSIONS - 1],

			// nutrients & metabolites
			vc->glucose,
			vc->oxygen,
			vc->ecm,

			vc->lactate,
			vc->morphogen,
			// change of nutrients & metabolites
			//0.,//voronoiDiagram->voronoiCells[index]->dglucose,
			//0.,//voronoiDiagram->voronoiCells[index]->doxygen,

			// cell state
			(float) (vc->getState() == COMPARTMENT ? (float) GetAgent(vc)->cellCount
					/ (float) GetAgent(vc)->maxCellCount
					: vc->getState() == ACTIVE),
//					(voronoiDiagram->voronoiCells[index]->getState() == ACTIVE)
			//&& cATP > THRESHOLD_QUIESCENCE_ATP
			/* &&
			 voronoiDiagram->voronoiCells[index]->glucose>THRESHOLD_QUIESCENCE_GLUCOSE &&
			 cATP>THRESHOLD_QUIESCENCE_ATP  ? 1 : 0),*/

			10 * (vc->getState() == NONACTIVE),

			11 * (vc->getState() == NECROTIC), //||

			12 * ((vc->getState() == FREE && GetAgent(vc) != NULL) ? 1 : 0),

			13 * (vc->getState() == VESSEL ? 1 : 0),

			// ATP
			//cATP
			cATP,
			//α = min{qG ,qO /6}/qG ,
			alpha,

			vc->CPT11in, vc->CPT11out,

			vc->CPT11in * 8000., vc->CPT11out * 1000.,

			vc->waste,

			(GetAgent(vc) == NULL ? 0 : GetAgent(vc)->waste),

			( GetAgent(vc) ? GiveMeTheATP( GetAgent(vc)) : 0)

			);

}

void writeSliceInformation(VoronoiDiagram* voronoiDiagram, AgentList *agentList,
		char *dirname) {

	//int dimSize = (int) (pow(voronoiDiagram->countVoronoiCells, 1. / DIMENSIONS) + 0.5);
	int dimSize = voronoiDiagram->xN[0];

	//int dimSize = voronoiDiagram->countVoronoiCells;
	//printf( "dimSize:%i\n", dimSize);
	int i, ii, iii; // = dimSize/2;
	//printf( "indexBase:%i\n", indexBase);

	/*if( agentList->countActiveAgents){

	 // center of mass of agents
	 double centerOfMass[DIMENSIONS];
	 for( k=0; k<DIMENSIONS; k++)
	 centerOfMass[k] = 0.0;
	 for( i=0; i<agentList->countActiveAgents; i++)
	 for( ii=0; ii<agentList->agents[i]->countLocations; ii++)
	 centerOfMass[0] += agentList->agents[i]->location[ii]->position[k];
	 //for( k=0; k<DIMENSIONS; k++)
	 //	centerOfMass[k] += agentList->agents[i]->location[ii]->position[k];

	 iii = (int)(centerOfMass[0]);
	 }else*/

	// open file to write
	FILE *fp;
	char outfilename[512];
	sprintf(outfilename, "%s/slice.dat", dirname);
	if ((fp = fopen(outfilename, "w")) == NULL) {
		fprintf(stderr, "Error opening file %s for writing!\n", outfilename);
	}
	//fprintf( stderr, "DIMENSIONS:%i\n", DIMENSIONS);
	// write slice data
#if DIMENSIONS == 1

	ii = 0;
	iii = 0;

	for (i = 0; i < voronoiDiagram->countVoronoiCells; i++) {
		//if( GetAgent( voronoiDiagram->voronoiCells[i]) != NULL)
		//fprintf( stderr, "%s \n", cellTypeToString( GetAgent( voronoiDiagram->voronoiCells[i])->state));

		writeVoronoiCellInformation( voronoiDiagram->voronoiCells[i], fp);
	}

#elif DIMENSIONS == 2

	/*for (ii = 0; ii < voronoiDiagram->xN[0]; ii++)
	 for (i = 0; i < voronoiDiagram->xN[1]; i++) {

	 int index = ii * dimSize + i;

	 writeVoronoiCellInformation( voronoiDiagram->voronoiCells[index], fp);
	 }*/

	for (i = 0; i < voronoiDiagram->countVoronoiCells; i++) {
		//if( GetAgent( voronoiDiagram->voronoiCells[i]) != NULL)
		//fprintf( stderr, "%s \n", cellTypeToString( GetAgent( voronoiDiagram->voronoiCells[i])->state));
		if (i % voronoiDiagram->xN[0] == 0)
			fprintf(fp, "\n");
		writeVoronoiCellInformation(voronoiDiagram->voronoiCells[i], fp);
	}

#elif DIMENSIONS == 3

	iii = dimSize / 2;
	//int centerOfMassX = 0;
	//int count = 0;
	/*for( i=0; i<agentList->countActiveAgents; i++)
	 for( ii=0; ii<agentList->agents[i]->countLocations; ii++){
	 centerOfMassX += agentList->agents[i]->location[ii]->position[0];
	 count++;
	 }
	 iii = centerOfMassX/count;*/
	for (ii = 0; ii < dimSize; ii++) {
		for (i = 0; i < dimSize; i++) {

			int index = iii * dimSize * dimSize + ii * dimSize + i;

			writeVoronoiCellInformation( voronoiDiagram->voronoiCells[index], fp);
		}
		fprintf(fp, "\n");
	}

#endif

	// close file
	fclose(fp);
}

void statisticsOfSpheroid(VoronoiDiagram* voronoiDiagram,
		double** global_oxygen_concentration,
		double** global_glucose_concentration,
		double** global_active_cell_density,
		double** global_nonactive_cell_density,
		double** global_necrotic_cell_density, double radius,
		int radiusIntervals, int timeIndex) {

	int i, k, kk;
	double radiusStep = radius / (double) radiusIntervals;
	//double timeStep =   time / (double) timeIntervals;

	/*double** global_oxygen_concentration = ( double** ) calloc ( EndTime_int, sizeof( double*));
	 for(k = 0;k < EndTime_int; k++){
	 global_oxygen_concentration[k] = ( double* ) calloc ( EndRadius_int, sizeof( double));
	 for(kk = 0;kk < EndRadius_int; kk++){
	 global_oxygen_concentration[k][kk] = 0.0;
	 }
	 }*/

	// center of mass
	double centerOfMass[DIMENSIONS];
	int countFoundCells = 0;
	for (k = 0; k < DIMENSIONS; k++)
		centerOfMass[k] = 0.0;

	for (i = 0; i < voronoiDiagram->countVoronoiCells; i++) {
		if (voronoiDiagram->voronoiCells[i] != NULL) {
			//if( !voronoiDiagram->voronoiCells[i]->isFree()){
			for (k = 0; k < DIMENSIONS; k++)
				centerOfMass[k] += voronoiDiagram->voronoiCells[i]->position[k];
			countFoundCells++;
		}
	}
	if (countFoundCells)
		for (k = 0; k < DIMENSIONS; k++) {
			centerOfMass[k] /= (double) countFoundCells;
			//fprintf( stderr, "%lf\n", centerOfMass[k]);
		}
	else
		for (k = 0; k < DIMENSIONS; k++) {
			centerOfMass[k] =
					(voronoiDiagram->xMax[k] + voronoiDiagram->xMin[k]) / 2.;
			//fprintf( stderr, "%lf\n", centerOfMass[k]);
		}
	//exit( 0);

	// scan radial oxygen concentration
	//k = timeIndex;
	double distanceToMassCenter;
	int countCellsInInterval[radiusIntervals];
	for (i = 0; i < radiusIntervals; i++)
		countCellsInInterval[i] = 0;
	for (i = 0; i < voronoiDiagram->countVoronoiCells; i++) {
		//if( GetAgent(voronoiDiagram->voronoiCells[i])->state != FREE){
		// determine distance to mass center
		distanceToMassCenter = 0.;
		for (k = 0; k < DIMENSIONS; k++)
			distanceToMassCenter +=
					pow(
							centerOfMass[k]
									- voronoiDiagram->voronoiCells[i]->position[k],
							2);
		distanceToMassCenter = sqrt(distanceToMassCenter);

		//
		kk = (int) (distanceToMassCenter / radiusStep);
		if (kk < radiusIntervals)
			countCellsInInterval[kk]++;
		//}
	}

	for (i = 0; i < voronoiDiagram->countVoronoiCells; i++) {
		//if( GetAgent(voronoiDiagram->voronoiCells[i])->state != FREE){
		// determine distance to mass center
		distanceToMassCenter = 0.;
		for (k = 0; k < DIMENSIONS; k++)
			distanceToMassCenter +=
					pow(
							centerOfMass[k]
									- voronoiDiagram->voronoiCells[i]->position[k],
							2.);
		distanceToMassCenter = sqrt(distanceToMassCenter);

		//
		kk = (int) (distanceToMassCenter / radiusStep);
		if (kk >= radiusIntervals) {
			fprintf(stderr, "Index too big: %i (< %i)!!! <= %lf\n", kk,
					radiusIntervals, distanceToMassCenter);

			exit(0);
		}
		if (kk < radiusIntervals) {
			global_oxygen_concentration[timeIndex][kk] +=
					voronoiDiagram->voronoiCells[i]->oxygen
							/ (double) countCellsInInterval[kk];
			global_glucose_concentration[timeIndex][kk] +=
					voronoiDiagram->voronoiCells[i]->glucose
							/ (double) countCellsInInterval[kk];
		}
		if (voronoiDiagram->voronoiCells[i]->getState() == ACTIVE)
			global_active_cell_density[timeIndex][kk] +=
					1. / (double) countCellsInInterval[kk];
		if (voronoiDiagram->voronoiCells[i]->getState() == NONACTIVE)
			global_nonactive_cell_density[timeIndex][kk] +=
					1. / (double) countCellsInInterval[kk];
		if (voronoiDiagram->voronoiCells[i]->getState() == NECROTIC)
			global_necrotic_cell_density[timeIndex][kk] +=
					1. / (double) countCellsInInterval[kk];

		//}
	}
}

double InitialOxygenConcentration = 1.;
double InitialGlucoseConcentration = 1.;

#ifndef __WithGUI__
double montecarlo(int argc, char **argv)
#else
#include <gui/window.h>
		int montecarlo( int argc, char **argv, Window *gui)
#endif
		{
/*	fprintf(stderr, "DIMENSIONS = %i\n", DIMENSIONS);

	fprinttextcolor((USE_DIVISION ? GREEN : RED));
	fprintf(stderr, "DIVISION [%s]\n", (USE_DIVISION ? "ON" : "OFF"));
	fprinttextcolor((USE_MIGRATION ? GREEN : RED));
	fprintf(stderr, "MIGRATION [%s]\n", (USE_MIGRATION ? "ON" : "OFF"));
	fprinttextcolor((USE_APOPTOSIS ? GREEN : RED));
	fprintf(stderr, "APOPTOSIS [%s]\n", (USE_APOPTOSIS ? "ON" : "OFF"));
	fprinttextcolor((USE_NECROSIS ? GREEN : RED));
	fprintf(stderr, "NECROSIS [%s]\n", (USE_NECROSIS ? "ON" : "OFF"));
	fprinttextcolor((USE_LYSIS ? GREEN : RED));
	fprintf(stderr, "LYSIS [%s]\n", (USE_LYSIS ? "ON" : "OFF"));
	fprinttextcolor((USE_VESSEL_GROWTH ? GREEN : RED));
	fprintf(stderr, "VESSEL_GROWTH [%s]\n", (USE_VESSEL_GROWTH ? "ON" : "OFF"));
	fprinttextcolor((USE_GROWTH ? GREEN : RED));
	fprintf(stderr, "GROWTH [%s]\n", (USE_GROWTH ? "ON" : "OFF"));
	fprinttextattribute(RESET);*/

	int histogram[10000];
	int histogramDividing[10000];
	int histogramNecrotic[10000];
	int histogramFree[10000];
	double histogramECM[10000];

	// AVERAGES & STD DERIVATION
	int histogramSquaresDividing[10000];
	int histogramSquaresNecrotic[10000];
	int histogramSquaresDividingFraction[10000];
	int histogramSquaresNecroticFraction[10000];
	int histogramSquaresFree[10000];
	double histogramSquaresECM[10000];

	int histogramSumDividing[10000];
	int histogramSumNecrotic[10000];
	int histogramSumDividingFraction[10000];
	int histogramSumNecroticFraction[10000];
	int histogramSumFree[10000];
	double histogramSumECM[10000];
	// -

	for (int i = 0; i < 10000; i++) {
		histogram[i] = 0;

		histogramDividing[i] = 0;
		histogramNecrotic[i] = 0;
		histogramFree[i] = 0;
		histogramECM[i] = 0.;

		histogramSquaresDividing[i] = 0;
		histogramSquaresNecrotic[i] = 0;
		histogramSquaresDividingFraction[i] = 0;
		histogramSquaresNecroticFraction[i] = 0;
		histogramSquaresFree[i] = 0;
		histogramSquaresECM[i] = 0.;

		histogramSumDividing[i] = 0;
		histogramSumNecrotic[i] = 0;
		histogramSumDividingFraction[i] = 0;
		histogramSumNecroticFraction[i] = 0;
		histogramSumFree[i] = 0;
		histogramSumECM[i] = 0.;
	}

	int i, k, l;
	//int j;
	unsigned short seed_start = 0;
	int c;
	//Grid *voronoiGrid;    // Grid structure
#if( DETERMINE_COARSE_GRAINED_BEHAVIOR)
	Grid *coarseVoronoiGrid; // Grid structure for the coarse grained model
#endif
	ActionTree *actionList; // probability list structure
	Action *selected_action;

	// File and Directory Variables
	char parameter_file[FILENAMESIZE];
	char dirname[FILENAMESIZE];
	char outfilename[FILENAMESIZE];
	char parameterfilename[FILENAMESIZE];
	//char init_filename[ FILENAMESIZE ]; 

	char filename_fineGrid[FILENAMESIZE] = "";
	char filename_coarseGrid[FILENAMESIZE];

	//char suffix[FILENAMESIZE],seedstring[FILENAMESIZE];//,string[FILENAMESIZE];
	char pnm_suffix[FILENAMESIZE + 10], pnmstring[FILENAMESIZE];
	parameter_file[0] = '\0';

	// Filepointer
	FILE *fp;
	FILE *fp_param;
#if SINGLE_OUTPUT 
	FILE *fp_div_time;
	FILE *fp_cells_vs_time;
	FILE *fp4;
	FILE *fp5;
#endif

	char boundaryCondition[FILENAMESIZE] = "";

	// Statistic Variables
	int count_inoculated_cells = 0;
	double count_divided_cells = 0;
	int count_expanded_cells = 0;
	int count_cells = 0;
	int count_cell_volume = 0;
	int count_vessel_cells = 0;
	int count_necrotic_cells = 0;
	int count_necrotic_volume = 0;
	int count_inner_empty_volume = 0;
	int count_lysed_cells = 0;
	int count_lysed_volume = 0;
	double gyrRadius = 0.;
	double prob_X_equal_k = 0.;
	double prob_X_smaller_k = 0.;
	double prob_X_bigger_k = 0.;

	Statistix myStatistics;

	myStatistics.count_cells = &count_cells;
	myStatistics.count_cell_volume = &count_cell_volume;
	myStatistics.count_expanded_cells = &count_expanded_cells;
	myStatistics.count_divided_cells = &count_divided_cells;
	myStatistics.count_vessel_cells = &count_vessel_cells;
	myStatistics.count_necrotic_cells = &count_necrotic_cells;
	myStatistics.count_necrotic_volume = &count_necrotic_volume;
	myStatistics.count_inner_empty_volume = &count_inner_empty_volume;
	myStatistics.count_lysed_cells = &count_lysed_cells;
	myStatistics.count_lysed_volume = &count_lysed_volume;

	//int global_count_divided_cells = 0;
	//int global_count_necrotic_cells = 0;
	//int global_count_inoculated_cells = 0;
	//int nocim;

	int Last_NumberOfCells = 0;
	double Last_gyrRadius = 0;

	// Time
	double Time, Last_Time, EndEndTime = 120, EndTime = 120;
	//double BeginningTime = 0.001;
	int AdaptEndTime = FALSE;

	// Scaling Variables
	int EndTime_int = 0;
	double OutputRate = 1.;

	// Averaging Arrays
	double *global_cells;
	double *global_volume;
	double *global_inner_empty_volume;
	double *global_necrotic_cells;
	double *global_necrotic_volume;
	double *global_vessel_cells;
	double *global_prob_X_equal_k;
	double *global_prob_X_smaller_k;
	double *global_prob_X_bigger_k;

	// Time & Timer Variables
	time_t timer;
	time_t temp_t;
	char local_time[256];
	struct tm *loctime = (struct tm *) malloc(sizeof(struct tm));
	long passedTime;
	long benchmarkTime = 0;

	int resolution = 10;

	// VESSEL STUFF
	double DistanceToInitialCell = 10.0;
	int CountVessel = 1;
	double BranchingProbability = 0.;
	int BranchingLength = 20;

	int Averages = 0;
	int NumberOfCellsPerSite = 2;
	//int M_div = 0;

	bool NoRadialProfiles = false;
	bool NoSliceOutput = false;
	bool PovrayOutput = false;

	int     RadialProfilesCount = 0;
	double *RadialProfilesTime  = 0;


#if( MONOD_KINETICS)
	//double **OxygenConcentrations;
	//double **GlucoseConcentrations;
	//double **GrowthFactorConcentrations;
	//double InitialOxygenConcentration = 0.;
	//double InitialGlucoseConcentration = 0.;
	double WilliamsCodeParameter = 0.;
	//double InitialGrowthFactorConcentration = 0.;
#endif

	// OUTPUTS
	int OutputAnimation = FALSE;
	char CustomDirectoryName = FALSE;
	int rotations = 0;
	double setoffX = 30, setoffY = 30, setoffZ = 30;
	double startAngle = -0.5 * PI;
	//int OutputPovrayFiles = FALSE;

	//double timeDifference = 0.;
	double customTimestep = 0.1;
	double customDiffusionCoefficientFactor = 1.;

	double custom_Waste_Diffusion = Waste_Diffusion;

	double InitialRadius = 1;
	double InitialQuiescentFraction = 0.;

	// DATA
	comparison_t data_growthcurve = create_comparison();
	comparison_t data_KI67_17 = create_comparison();
	comparison_t data_ECM_17  = create_comparison();
	comparison_t data_TUNEL_17= create_comparison();
	comparison_t data_KI67_24 = create_comparison();
	comparison_t data_ECM_24  = create_comparison();
	comparison_t data_TUNEL_24= create_comparison();
	//comparison_t sim_growthcurve = create_comparison();
	//sim_growthcurve.x = (double*) malloc( sizeof(double) * 1000);
	//sim_growthcurve.m = (double*) malloc( sizeof(double) * 1000);
	double maxRadius  = 0;
	double maxEpsilon = DBL_MAX;
	double cumEpsilon = 0;
	int    idxEpsilon = 0;
	double measurement_error = 0;
	double measurement_error_gc   = 0;
	double measurement_error_ki67 = 0;
	double measurement_error_ecm  = 0;
	double measurement_error_TUNEL= 0;
	{FILE *fp_raw = fopen( "raw_gc.dat", "w"); fprintf(fp_raw, ""); fclose(fp_raw); }
	//{FILE *fp_raw = fopen( "raw_ki67.dat", "w"); fprintf(fp_raw, ""); fclose(fp_raw); }


	// PARSING COMMAND LINE ARGUEMENTS
	while ((c =
			getopt(
					argc,
					argv,
					"1:2:3:4:5:6:7:8:9:R:a:b:c:d:e:f:g:i:jo:l:n:M:m:p:s:v:F:t:k:x:y:r:w:z:Z:D:A:NE:B:I:J:K:C:L:S:M:G:O:V:P:Y:W:h"))
			!= EOF) {
		switch (c) {

		case '1': {
			// Read Growth Curve Data
			   //fprintf( stderr, "Read data from %s\n", optarg) ;
			data_growthcurve.dim = readFileColumn( optarg, data_growthcurve.x, 0);
			                       readFileColumn( optarg, data_growthcurve.m, 1);
			                       readFileColumn( optarg, data_growthcurve.s, 2);
		   for( int j=0; j<data_growthcurve.dim; j++){
			   data_growthcurve.x[j] *= 24.;
			   maxRadius = max<double>( maxRadius, data_growthcurve.m[j] + data_growthcurve.s[j]*1);
		   }
		   //fprintf(stderr, "(( max radius = %e ))\n", maxRadius);

		} break;

		case '2': {
			// Read Proliferation Profile Data
			   //fprintf( stderr, "Read data from %s\n", optarg) ;
			   data_KI67_17.dim = readFileColumn( optarg, data_KI67_17.x, 0);
			   				      readFileColumn( optarg, data_KI67_17.m, 1);
			   				      readFileColumn( optarg, data_KI67_17.s, 2);
		} break;

		case '3': {
			// Read Proliferation Profile Data
			   //fprintf( stderr, "Read data from %s\n", optarg) ;
			   data_ECM_17.dim = readFileColumn( optarg, data_ECM_17.x, 0);
			   				     readFileColumn( optarg, data_ECM_17.m, 1);
			   				     readFileColumn( optarg, data_ECM_17.s, 2);
		} break;

		case '4': {
			// Read Proliferation Profile Data
			   fprintf( stderr, "Read data from %s\n", optarg) ;
			   data_TUNEL_17.dim = readFileColumn( optarg, data_TUNEL_17.x, 0);
			   				       readFileColumn( optarg, data_TUNEL_17.m, 1);
			   				       readFileColumn( optarg, data_TUNEL_17.s, 2);
		} break;

		case '5': {
			// Read Proliferation Profile Data
			   //fprintf( stderr, "Read data from %s\n", optarg) ;
			   data_KI67_24.dim = readFileColumn( optarg, data_KI67_24.x, 0);
			   				      readFileColumn( optarg, data_KI67_24.m, 1);
			   				      readFileColumn( optarg, data_KI67_24.s, 2);
		} break;

		case '6': {
			// Read Proliferation Profile Data
			   //fprintf( stderr, "Read data from %s\n", optarg) ;
			   data_ECM_24.dim = readFileColumn( optarg, data_ECM_24.x, 0);
			   				     readFileColumn( optarg, data_ECM_24.m, 1);
			   				     readFileColumn( optarg, data_ECM_24.s, 2);
		} break;

		case '7': {
			// Read Proliferation Profile Data
			   //fprintf( stderr, "Read data from %s\n", optarg) ;
			   data_TUNEL_24.dim = readFileColumn( optarg, data_TUNEL_24.x, 0);
			   				       readFileColumn( optarg, data_TUNEL_24.m, 1);
			   				       readFileColumn( optarg, data_TUNEL_24.s, 2);
		} break;

		case '8': {
			maxEpsilon = atof( optarg);
			//fprintf( stderr, "Pass epsilon limit %e\n", maxEpsilon) ;
		} break;

		case '9': {
			measurement_error = atof( optarg);
			//fprintf( stderr, "Pass epsilon limit %e\n", maxEpsilon) ;
		} break;

		case 'R': {

			// Parameter Value
			//double value = atof(optarg);
			double value = atof(&optarg[strcspn(optarg, "0123456789.")]);
			//fprintf( stderr, "optarg = %s to %f\n", optarg, value);
			if (strstr(optarg, "alpha") != 0) {
				alpha_ref = value;
				Agent::USE_GRAVITY = true;
			} else if (strstr(optarg, "InitialRadius") != 0) {
				InitialRadius = value;
			} else if (strstr(optarg, "InitialQuiescentFraction") != 0) {
				InitialQuiescentFraction = value;
			} else if (strstr(optarg, "WasteDiffusion") != 0) {
				custom_Waste_Diffusion = value;
			} else if (strstr(optarg, "ExponentialReentranceProbability") != 0) {
				ReentranceProbabilityFunction = EXPONENTIAL;
			} else if (strstr(optarg, "ReentranceProbabilityLength") != 0) {
				ReentranceProbabilityReferenceLength = value;
			} else if (strstr(optarg, "ShiftToUnoccupied") != 0) {
				VoronoiCell::SHIFT_TO_UNOCCUPIED = true;
			} else if (strstr(optarg, "ECMThresholdQuiescence") != 0) {
				//fprintf( stderr, "set ECM_THRESHOLD_QUIESCENCE to %f\n", value);
				VoronoiCell::ECM_THRESHOLD_QUIESCENCE = value;
			} else if (strstr(optarg, "ECMProductionRate") != 0) {
				VoronoiCell::ECM_PRODUCTION_RATE = value;
			} else if (strstr(optarg, "ECMDegradationRate") != 0) {
				VoronoiCell::ECM_DEGRADATION_RATE = value;
			} else if (strstr(optarg, "LactateThresholdQuiescence") != 0) {
				VoronoiCell::LACTATE_THRESHOLD_QUIESCENCE = value;
			} else if (strstr(optarg, "LactateThresholdDeath") != 0) {
				VoronoiCell::LACTATE_THRESHOLD_DEATH = value;
			} else if (strstr(optarg, "WasteThresholdQuiescence") != 0) {
				VoronoiCell::WASTE_THRESHOLD_QUIESCENCE = value;
			} else if (strstr(optarg, "WasteThresholdDeath") != 0) {
				VoronoiCell::WASTE_THRESHOLD_DEATH = value;
			} else if (strstr(optarg, "WasteUptake") != 0) {
				Agent::WASTE_UPTAKE = value;
			} else if (strstr(optarg, "WasteThresholdQuiescenceIntracellular") != 0) {
					Agent::WASTE_THRESHOLD_QUIESCENCE = value;
			} else if (strstr(optarg, "WasteThresholdSlowedGrowth") != 0) {
					Agent::WASTE_THRESHOLD_SLOWED_GROWTH = value;
			} else if (strstr(optarg, "WasteIntoxicatedCellCycles") != 0) {
				Agent::WASTE_INTOXICATED_CELL_CYCLES = (int)value;
			} else if (strstr(optarg, "ATPThresholdQuiescence") != 0) {
				VoronoiCell::ATP_THRESHOLD_QUIESCENCE = value;
			} else if (strstr(optarg, "ATPThresholdDeath") != 0) {
				VoronoiCell::ATP_THRESHOLD_DEATH = value;
			} else if (strstr(optarg, "GlucoseOxygenProductThresholdDeath") != 0) {
				VoronoiCell::GLUCOSE_OXYGEN_PRODUCT_THRESHOLD_DEATH = value;
			} else if (strstr(optarg, "GlucoseOxygenProductThresholdQuiescence")!= 0) {
				VoronoiCell::GLUCOSE_OXYGEN_PRODUCT_THRESHOLD_QUIESCENCE = value;
			} else if (strstr(optarg, "DynamicExtendedNeighborhood") != 0) {
				VoronoiCell::USE_DYNAMIC_EXTENDED_NEIGHBORHOOD = true;
			} else if (strstr(optarg, "SymbolicExtendedNeighborhood") != 0) {
				VoronoiCell::USE_SYMBOLIC_EXTENDED_NEIGHBORHOOD = true;
			} else if (strstr(optarg, "Lactate") != 0) {
				VoronoiCell::USE_LACTATE = true;
			} else if (strstr(optarg, "Waste") != 0) {
				VoronoiCell::USE_WASTE = true;
			} else if (strstr(optarg, "Morphogen") != 0) {
				VoronoiCell::USE_MORPHOGEN = true;
			} else if (strstr(optarg, "MaxGlucoseConsumption") != 0) {
				Agent::NICK_G_MAX = value * 1e-17;
			} else if (strstr(optarg, "PovrayOutput") != 0) {
				PovrayOutput = true;
			} else if (strstr(optarg, "NoSliceOutput") != 0) {
				NoSliceOutput = true;
			} else if (strstr(optarg, "NoRadialProfiles") != 0) {
				NoRadialProfiles = true;
			} else if (strstr(optarg, "RadialProfilesTime") != 0) {

				RadialProfilesCount++;
				RadialProfilesTime = (double*) realloc( RadialProfilesTime, sizeof(double) * RadialProfilesCount);
				RadialProfilesTime[RadialProfilesCount-1] = (int)value;

			} else if (strstr(optarg, "apt") != 0) {
				CellApoptosisRate = value;
			} else if (strstr(optarg, "nec") != 0) {
				CellNecrosisRate = value;
			} else if (strstr(optarg, "lys") != 0) {
				CellLysisRate = value;
			} else if (strstr(optarg, "mi") != 0) {
				CellMigrationRate = value;
			} else if (strstr(optarg, "re") != 0) {
				Agent::ReentranceRate = value;
			} else {
				fprintf(stderr, "\nWARNING: Unknown Parameter: %s\n", optarg);
				fprintf(stderr, "USAGE: -R Option\n");
				fprintf(stderr, "apt [float]\n");
				fprintf(stderr, "nec [float]\n");
				fprintf(stderr, "InitialRadius [float]\n");
				fprintf(stderr, "ExponentialReentranceProbability\n");
				fprintf(stderr, "ShiftToUnoccupied\n");
				fprintf(stderr, "ATPThresholdQiescence [float]\n");
				fprintf(stderr, "ATPThresholdDeath [float]\n");
				fprintf(stderr, "GlucoseOxygenProductThresholdDeath [float]\n");
				fprintf(stderr,
						"GlucoseOxygenProductThresholdQuiescence [float]\n");
				exit(0);
			}

			break;
		}

		case 'j':
			AdaptEndTime = TRUE;
			break;
		case 'Y':
			DistanceToInitialCell = atof(optarg);
			break;
		case 'W':
			CountVessel = atoi(optarg);
			break;
		case 'P':
			BranchingProbability = atof(optarg);
			break;
		case 'b':
			BranchingLength = atoi(optarg);
			break;
		case 'd':
			strncpy(dirname, optarg, FILENAMESIZE);
			CustomDirectoryName = TRUE;
			//fprintf( stderr, "Create dir %s\n", optarg) ;
			break;
		case 'c':
			NumberOfCellsPerSite = atoi(optarg);
			break;
		case 'p':
			strncpy(parameter_file, optarg, FILENAMESIZE);
			break;
		case 'i':
			strncpy(filename_fineGrid, optarg, FILENAMESIZE);
			break;
		case 'g':
			strncpy(filename_coarseGrid, optarg, FILENAMESIZE);
			break;
		case 'o':
			OutputRate = atof(optarg);
			break;
		case 'm':
			M_div = atoi(optarg) - 1;
			break;
		case 'M':
			M_gro = atoi(optarg) - 1;
			break;
		case 'n':
			CountCellsPerVoronoiCell = atoi(optarg);
			break;

		case 's':
			seed_start = atoi(optarg);
			break;
		case 'v':
			MaxCellDivisionRate = atof(optarg);
			break;
		case 'f':
			CellMigrationRate = atof(optarg);
			break;
		case 'z':
			CellApoptosisRate = atof(optarg);
			break;
		case 'V':
			VesselGrowthRate = atof(optarg);
			break;
		case 'k':
			CellDivisionDepth = atof(optarg);
			break;
		case 'l':
			CellLysisRate = atof(optarg);
			break;
		case 't':
			customTimestep = atof(optarg);
			break;
		case 'B':
			strncpy(boundaryCondition, optarg, FILENAMESIZE);
			break;
		case 'x':
			Averages = atoi(optarg);
			//Averages_Meanfield = Averages;
			break;
		case 'y':
			EndEndTime = EndTime = atof(optarg);
			break;
		case 'r':
			resolution = atoi(optarg);
			break;
		//case 'A':
		//	OutputAnimation = atoi(optarg);
		//	break;
		case 'D':
			customDiffusionCoefficientFactor = atof(optarg);
			Oxygen_Diffusion *= customDiffusionCoefficientFactor;
			Glucose_Diffusion *= customDiffusionCoefficientFactor;
			break;
		case 'N':
			//numericalApproache( argc, argv);
			break;
			/*case 'D':
			 Type_of_diffusion=atoi(optarg);
			 break;*/
			/*case 'N':
			 nu=atof(optarg);
			 break;
			 case 'E':
			 E_0=atof(optarg);
			 break;
			 case 'B':
			 E_b=atof(optarg);
			 break;*/
#if MONOD_KINETICS
		case 'O':
			InitialOxygenConcentration = atof(optarg);
			break;
		case 'G':
			InitialGlucoseConcentration = atof(optarg);
			break;
		case 'w':
			WilliamsCodeParameter = atof(optarg);
			break;
#endif
		case 'C':
			Case = atoi(optarg);
			break;
		case 'h':
		default:
			fprintf(stdout, "usage: %s [options]\n", argv[0]);
			fprintf(stdout, "\toptions are:\n");
			fprintf(stdout,
					"\t-d dirname (Creates directory to store files in)\n");
			fprintf(stdout, "\t-i initfilename (init with this file)\n");
			fprintf(stdout, "\t-p parameter-file\n");
			fprintf(stdout,
					"\t-o output rate ( time step between measuring points)\n");
			fprintf(stdout, "\t-m number of interstates M\n");
			fprintf(stdout,
					"\t-s seed for random generator (not implemented yet!)\n");
			fprintf(stdout, "\t-v rate of division\n");
			fprintf(stdout, "\t-f rate of diffusion\n");
			fprintf(stdout, "\t-t rate of attachment\n");
			fprintf(stdout, "\t-e rate of detachment\n");
			fprintf(stdout, "\t-z rate of necrosis (only on grid)\n");
			fprintf(stdout, "\t-Z rate of necrosis in suspension\n");
			fprintf(stdout, "\t-k radius k of possible division and shift\n");
			fprintf(stdout, "\t-a mutation ratio in mitosis as an quotient\n");
			fprintf(stdout, "\t-x number of averages\n");
			fprintf(stdout, "\t-y time for process\n");
			fprintf(
					stdout,
					"\t-D type of diffusion (0, free, 1, limited to border, 2, energyminimizing\n");
			fprintf(stdout,
					"\t-A animation (1 for animation, default FALSE)\n");
			fprintf(stdout, "\t-Y DistanceToInitialCell\n");
			fprintf(stdout, "\t-V VesselGrowthRate\n");
			fprintf(stdout, "\t-W CountVessel\n");
			fprintf(stdout, "\t-P BranchingProbability\n");
			fprintf(stdout, "\t-b BranchingLength\n");
			fprintf(stdout, "\t-h this help\n");
			exit(0);
		}
	}
	srand(0);
	//customDiffusionCoefficientFactor /= pow( (double)CountCellsPerVoronoiCell,2./3.);
	Oxygen_Diffusion /= pow((double) CountCellsPerVoronoiCell, 2. / DIMENSIONS);
	Glucose_Diffusion /= pow((double) CountCellsPerVoronoiCell, 2. / DIMENSIONS);


	// ADDITIVE MEASUREMENT ERROR
	measurement_error_gc = 0;
	for( int j=0; j<data_growthcurve.dim; j++)
		measurement_error_gc = max<double>( measurement_error_gc, data_growthcurve.m[j]);
	measurement_error_gc *= measurement_error;

	measurement_error_ki67 = 0;
	for( int j=0; j<data_KI67_17.dim; j++)
		measurement_error_ki67 = max<double>( measurement_error_ki67, data_KI67_17.m[j]);
	measurement_error_ki67 *= measurement_error;

	measurement_error_ecm = 0;
	for( int j=0; j<data_ECM_17.dim; j++)
		measurement_error_ecm = max<double>( measurement_error_ecm, data_ECM_17.m[j]);
	measurement_error_ecm *= measurement_error;

	measurement_error_TUNEL = 0;
	for( int j=0; j<data_TUNEL_17.dim; j++)
		measurement_error_TUNEL = max<double>( measurement_error_TUNEL, data_TUNEL_17.m[j]);
	measurement_error_TUNEL *= measurement_error;


	// READ DATA





	// FILE AND DIRECTORY NAMES

	if (!CustomDirectoryName) {
#ifdef SQUARE_LATTICE
		sprintf( dirname, "C%iD%iM%im%ik%in%iSL", Case, DIMENSIONS, M_gro+1, M_div+1, (int)CellDivisionDepth, CountCellsPerVoronoiCell );
#else
		sprintf(dirname, "C%iD%iM%im%ik%in%i", Case, DIMENSIONS, M_gro + 1,
				M_div + 1, (int) CellDivisionDepth, CountCellsPerVoronoiCell);
#endif
	}

	mkdir(dirname MODUS);
	//if (mkdir(dirname MODUS))
	//	fprintf(stderr, "WARNING: Error creating directory %s. Exiting\n", dirname);
	//else
	//	fprintf(stderr, "INFO: Created directory %s\n", dirname);

	// INITIALIZE VORONOI GRID

	// read voronoi grid
	passedTime = clock();


#ifdef SQUARE_LATTICE

	fprintf( stderr, "Create Square Lattice... \n");

#if DIMENSIONS == 1
	// TODO:
	//int countPoints[DIMENSIONS] = {1720}; 
	int countPoints[DIMENSIONS] = {200./CountCellsPerVoronoiCell+2};
	int periodic[DIMENSIONS] = {0};

#elif DIMENSIONS == 2

	//int countPoints[DIMENSIONS+1] = {100, 100, 1};
	//int countPoints[DIMENSIONS+1] = {40, 40, 1};
	int countPoints[DIMENSIONS+1] = {2+100/pow(CountCellsPerVoronoiCell,1./DIMENSIONS), 2+100/pow(CountCellsPerVoronoiCell,1./DIMENSIONS), 1};
	//int countPoints[DIMENSIONS+1] = {2, 2, 1};
	//int countPoints[DIMENSIONS] = {40, 40};
	int periodic[DIMENSIONS] = {0, 0};

#elif DIMENSIONS == 3

	//int countPoints[DIMENSIONS] = {3, 3, 3};
	//int countPoints[DIMENSIONS] = {10, 10, 10};
	int countPoints[DIMENSIONS] = {100, 100, 100};
	//int countPoints[DIMENSIONS] = {100, 100, 100};

	//int countPoints[DIMENSIONS] = {2+60/pow(CountCellsPerVoronoiCell,1./3.), 2+60/pow(CountCellsPerVoronoiCell,1./3.), 2+60/pow(CountCellsPerVoronoiCell,1./3.)};
	//int countPoints[DIMENSIONS] = { 80, 80, 80}; 
	//int countPoints[DIMENSIONS] = { 60, 60, 60}; 
	//int countPoints[DIMENSIONS] = { 50, 50, 50}; 
	//int countPoints[DIMENSIONS] = { 30, 30, 30}; 

	//int countPoints[DIMENSIONS] = { 1000, 1, 1}; 
	//int countPoints[DIMENSIONS] = { 20, 20, 20};

	//int countPoints[DIMENSIONS] = { 10, 10, 10};
	//int countPoints[DIMENSIONS] = { 3, 3, 3}; 
	int periodic[DIMENSIONS] = {0, 0, 0};

#endif

	VoronoiDiagram * voronoiDiagram = newSquareLattice( countPoints, periodic);
	//VoronoiDiagram * voronoiDiagram = VoronoiDiagram::newVoronoiDiagram( countPoints[0], countPoints[1], countPoints[2]);

#else

	VoronoiDiagram * voronoiDiagram;
	int countPoints[DIMENSIONS];

	if( strlen(filename_fineGrid)){

		//fprintf(stderr, "Reading Voronoi Lattice from File... \n");
		voronoiDiagram= new VoronoiDiagram(filename_fineGrid);


		for (i = 0; i < DIMENSIONS; i++)
			countPoints[i] = (int) (pow(voronoiDiagram->countVoronoiCells, 1. / DIMENSIONS)	+ 0.5);
	}else{
		for (i = 0; i < DIMENSIONS; i++)
			countPoints[i] = 100;

		char latticeFilename[512], testFilename[512];
		int  latticeSize = 100;
		sprintf(latticeFilename, "%ipow%i", latticeSize, DIMENSIONS);
		sprintf(testFilename,    "%ipow%i.qele", latticeSize, DIMENSIONS);

		if( access( testFilename, F_OK ) != -1 ) {
		    // file exists
			//fprintf(stderr, "Reading Voronoi Lattice from File... \n");
			voronoiDiagram= new VoronoiDiagram(latticeFilename);
		} else {
		    // file doesn't exist
			//fprintf(stderr, "Create new Voronoi Lattice \n");
			voronoiDiagram = VoronoiDiagram::newVoronoiDiagram( latticeSize, latticeSize, (DIMENSIONS == 2 ? 1 : latticeSize));
			voronoiDiagram->writeToFile( latticeFilename);
		}
		//exit( 0);
	}
	voronoiDiagram->boundaryThickness = 1.;


	// STATISTICS
	/*double meanNeighborDistance = 0;
	int meanNeighborDistanceCount = 0;
	for( int i=0; i<voronoiDiagram->countVoronoiCells; i++)
		for( int n=0; n<voronoiDiagram->voronoiCells[i]->countNeighborCells; n++)
		if( voronoiDiagram->voronoiCells[i]->position[0]>2. && voronoiDiagram->voronoiCells[i]->position[0] < 58. &&
			voronoiDiagram->voronoiCells[i]->position[1]>2. && voronoiDiagram->voronoiCells[i]->position[1] < 58. &&
			voronoiDiagram->voronoiCells[i]->position[2]>2. && voronoiDiagram->voronoiCells[i]->position[2] < 58. 	){
		meanNeighborDistance += voronoiDiagram->voronoiCells[i]->getDistanceTo(
				voronoiDiagram->voronoiCells[i]->neighborCells[n]);
		meanNeighborDistanceCount++;
	}
	fprintf(stderr, "meanNeighborDistance=%lf\n", meanNeighborDistance/meanNeighborDistanceCount);*/
	// END STATISTICS

#endif

	voronoiDiagram->setDomain();
	//fprintf(stderr, "...finished ( %lisec)\n",
	//		(clock() - passedTime) / CLOCKS_PER_SEC);

	/*passedTime = clock();
	 fprintf( stderr, "Set Voronoi Grid... \n");
	 voronoiDiagram->setVoronoiGrid();
	 fprintf( stderr, "...finished ( %lisec)\n", (clock() - passedTime)/CLOCKS_PER_SEC);
	 */

	/*for( int t=1; t<voronoiDiagram->countTetrahedra; t++){
	 for( int i=0; i<DIMENSIONS+1; i++){
	 fprintf( stderr, "tet->vertice[%i] (%lf, %lf, %lf)\n", i, voronoiDiagram->tetrahedra[t]->vertices[i]->position[0], voronoiDiagram->tetrahedra[t]->vertices[i]->position[1], voronoiDiagram->tetrahedra[t]->vertices[i]->position[2]);
	 }
	 fprintf( stderr, "center (%lf, %lf, %lf)\n",voronoiDiagram->voronoiGridPoints[t].position[0],voronoiDiagram->voronoiGridPoints[t].position[1],voronoiDiagram->voronoiGridPoints[t].position[2]);
	 }
	 exit( 0);*/

	// alternative methode!
	/*passedTime = clock();
	 fprintf( stderr, "alternative methode!... \n");
	 int ii,j;
	 //voronoiDiagram->countVoronoiCells = 100;
	 char ***matrixA = (char***) calloc( voronoiDiagram->countVoronoiCells, sizeof( char**));
	 for( k=0; k<2; k++){
	 matrixA[k] = (char**) calloc( voronoiDiagram->countVoronoiCells, sizeof( char*));
	 for( i=0; i<voronoiDiagram->countVoronoiCells; i++){
	 //fprintf( stderr, "\r%i \b");
	 matrixA[k][i] = (char*) calloc( voronoiDiagram->countVoronoiCells, sizeof( char));
	 for( ii=0; ii<voronoiDiagram->countVoronoiCells; ii++)
	 matrixA[k][i][ii] = 0;
	 }
	 }
	 // init first adjazenz matrix
	 for( i=0; i<voronoiDiagram->countVoronoiCells; i++){
	 matrixA[0][i][i] = 1;
	 for( ii=i+1; ii<voronoiDiagram->voronoiCells[i]->countNeighborCells; ii++){
	 matrixA[0][i][voronoiDiagram->voronoiCells[i]->neighborCells[ii]->index] = 1;
	 }
	 }

	 for( i=0; i<voronoiDiagram->countVoronoiCells; i++){
	 fprintf( stderr, "\rposition [%i, ?]            \b", i);
	 for( ii=i; ii<voronoiDiagram->countVoronoiCells; ii++){
	 //fprintf( stderr, "\rposition [%i, %i]            \b", i, ii);
	 if( matrixA[0][i][ii] == 1)
	 matrixA[1][i][ii] = 1;
	 else
	 for( j=0; j<voronoiDiagram->countVoronoiCells; j++){
	 //fprintf( stderr, "\rposition [%i, %i] -> %i           \b", i, ii, j);
	 //matrixA[1][i][ii] += matrixA[0][i][j]&&matrixA[0][j][ii];
	 if( (j>=i ? matrixA[0][i][j] : matrixA[0][j][i])&&(ii>=j ? matrixA[0][j][ii] : matrixA[0][ii][j]))
	 matrixA[1][i][ii] = 1;
	 }
	 }
	 }



	 fprintf( stderr, "...finished ( %lisec)\n", (clock() - passedTime)/CLOCKS_PER_SEC);
	 */

	// set agent array
//	passedTime = clock();
//	fprintf(stderr, "Initializing Agents...\n");
	//AgentList* agentArray;
	AgentList* agentArray = AgentList::newAgentList(
			voronoiDiagram->countVoronoiCells);
//	fprintf(stderr, "...finished ( %lisec)\n",
//			(clock() - passedTime) / CLOCKS_PER_SEC);

	// set extended neighborhood
//	passedTime = clock();
//	fprintf(stderr, "Initialize Shifting Neighborhoods of Agents...\n");
	VoronoiCell::extendedNeighborDepth = CellDivisionDepth;
	if (!(VoronoiCell::USE_DYNAMIC_EXTENDED_NEIGHBORHOOD
			|| VoronoiCell::USE_SYMBOLIC_EXTENDED_NEIGHBORHOOD)
			&& CellDivisionDepth > 1) {
//	if ( !(VoronoiCell::extendedNeighborhoodType==VoronoiCell::DYNAMIC||VoronoiCell::extendedNeighborhoodType==VoronoiCell::SYMBOLIC) && CellDivisionDepth > 1) {

#ifdef SQUARE_LATTICE

		// SQUARE LATTICE
		double innerDomainRadius = voronoiDiagram->xMax[0] - voronoiDiagram->xMin[0];
		for( i=1; i<DIMENSIONS; i++) {
			//			if( innerDomainRadius > voronoiDiagram->xMax[i] - voronoiDiagram->xMin[i])
			// largest inner domaine radius
			if( innerDomainRadius < voronoiDiagram->xMax[i] - voronoiDiagram->xMin[i])
			innerDomainRadius = voronoiDiagram->xMax[i] - voronoiDiagram->xMin[i];
		}
		innerDomainRadius /= 2.;
		fprintf( stderr, "innerDomainRadius = %lf, CellDivisionDepth = %lf, CentralCell = %p\n", innerDomainRadius, CellDivisionDepth, getCentralCell( voronoiDiagram));
			//voronoiDiagram->setExtendedNeighborhoodAroundVoronoiCell( (int)CellDivisionDepth, (int)innerDomainRadius, getCentralCell( voronoiDiagram));
			//voronoiDiagram->setExtendedNeighborhood( CellDivisionDepth);
		voronoiDiagram->NEWsetExtendedNeighborhood((int) CellDivisionDepth);

#else

		// VORONOI DIAGRAM
		if (!voronoiDiagram->readExtendedNeighborhoodToFile(filename_fineGrid,
				(int) CellDivisionDepth)) {
			//fprintf( stderr, "Couldn't read Extended Neighborhood from file -> Calculate my self!\n");
			double innerDomainRadius = voronoiDiagram->xMax[0]
					- voronoiDiagram->xMin[0];
			for (i = 1; i < DIMENSIONS; i++) {
				if (innerDomainRadius
						> voronoiDiagram->xMax[i] - voronoiDiagram->xMin[i])
					innerDomainRadius =
							voronoiDiagram->xMax[i] - voronoiDiagram->xMin[i];
			}
			innerDomainRadius /= 2.;
			fprintf(
					stderr,
					"innerDomainRadius = %lf, CellDivisionDepth = %lf, CentralCell = %p\n",
					innerDomainRadius, CellDivisionDepth,
					getCentralCell(voronoiDiagram));
			//voronoiDiagram->setExtendedNeighborhoodAroundVoronoiCell( (int)CellDivisionDepth, (int)innerDomainRadius, getCentralCell( voronoiDiagram));
			voronoiDiagram->NEWsetExtendedNeighborhood((int) CellDivisionDepth);
			//fprintf( stderr, "Write Extended Neighborhood to file\n");
			voronoiDiagram->writeExtendedNeighborhoodToFile(filename_fineGrid,
					(int) CellDivisionDepth);
			//exit( 0);
		}

#endif

	}
//	fprintf(stderr, "...finished ( %lisec)\n",
//			(clock() - passedTime) / CLOCKS_PER_SEC);

	// initialize diffusion/kinetics

#if MONOD_KINETICS
//	passedTime = clock();
//	fprintf(stderr, "Initialize Diffusion...\n");

#ifndef __WithGUI__

	Substrate substrate(countPoints, Case, voronoiDiagram, agentArray,
			InitialOxygenConcentration, InitialGlucoseConcentration,
			WilliamsCodeParameter, dirname, boundaryCondition, customTimestep);
	//Substrate substrate;

	// INIT SOLVER
	/*v0 = (float*) malloc( sizeof(float) * voronoiDiagram->countVoronoiCells);
	 v1 = (float*) malloc( sizeof(float) * voronoiDiagram->countVoronoiCells);
	 v2 = (float*) malloc( sizeof(float) * voronoiDiagram->countVoronoiCells);
	 v3 = (float*) malloc( sizeof(float) * voronoiDiagram->countVoronoiCells);
	 v4 = (float*) malloc( sizeof(float) * voronoiDiagram->countVoronoiCells);
	 v5 = (float*) malloc( sizeof(float) * voronoiDiagram->countVoronoiCells);
	 v6 = (float*) malloc( sizeof(float) * voronoiDiagram->countVoronoiCells);*/

	DiffusionReactionEquation *lactateDynamics = 0;
	float *lactate = 0;
	float *lactateProduction = 0;
	float *lactateDiffusion = 0;

	if (VoronoiCell::USE_LACTATE) {
		fprintf(stderr, "Initialize Lactate Dynamics\n");
		lactateDynamics =
				new DiffusionReactionEquation(voronoiDiagram,
						DiffusionReactionEquation::STEADYSTATE, //time
						DiffusionReactionEquation::EXPLICIT, //reaction
						DiffusionReactionEquation::IMPLICIT //diffusion
						);
		lactateDynamics->setSpaceStep(SPATIAL_UNIT);
		lactate =
				(float*) malloc(
						sizeof(float) * voronoiDiagram->countVoronoiCells);
		lactateProduction =
				(float*) malloc(
						sizeof(float) * voronoiDiagram->countVoronoiCells);
		lactateDiffusion =
				(float*) malloc(
						sizeof(float) * voronoiDiagram->countVoronoiCells);

		lactateDynamics->initSystem(lactate);
		lactateDynamics->initReaction(lactateProduction);
		lactateDynamics->initDiffusion(lactateDiffusion);

		for (int i = 0; i < voronoiDiagram->countVoronoiCells; i++) {
			lactate[i] = 0.;
			lactateDiffusion[i] = Lactate_Diffusion;
		}
	}

	DiffusionReactionEquation *wasteDynamics = 0;
	float *waste = 0;
	float *wasteProduction = 0;
	float *wasteUptake = 0;
	float *wasteDiffusion = 0;

	if (VoronoiCell::USE_WASTE) {
		fprintf(stderr, "Initialize Waste Dynamics\n");
		wasteDynamics =
				new DiffusionReactionEquation(voronoiDiagram,
						DiffusionReactionEquation::IMPLICIT, //time
						DiffusionReactionEquation::IMPLICIT, //reaction
						DiffusionReactionEquation::IMPLICIT //diffusion
						);
		wasteDynamics->setSpaceStep(SPATIAL_UNIT);
		waste =	(float*) malloc(sizeof(float) * voronoiDiagram->countVoronoiCells);
		wasteProduction = (float*) malloc(sizeof(float) * voronoiDiagram->countVoronoiCells);
		wasteUptake = (float*) malloc(sizeof(float) * voronoiDiagram->countVoronoiCells);
		wasteDiffusion = (float*) malloc( sizeof(float) * voronoiDiagram->countVoronoiCells);

		wasteDynamics->initSystem(waste);
		wasteDynamics->initReaction(wasteProduction);
		wasteDynamics->initReactionDerivative(wasteUptake);
		wasteDynamics->initDiffusion(wasteDiffusion);

		for (int i = 0; i < voronoiDiagram->countVoronoiCells; i++) {
			waste[i] = 0.;
			wasteDiffusion[i] = custom_Waste_Diffusion;
		}
	}

	DiffusionReactionEquation *morphogenDynamics = 0;
	float *morphogen = 0;
	float *morphogenProduction = 0;
	float *morphogenDecay = 0;
	float *morphogenDiffusion = 0;

	if (VoronoiCell::USE_MORPHOGEN) {
		fprintf(stderr, "Initialize Morphogen Dynamics\n");
		morphogenDynamics =
				new DiffusionReactionEquation(voronoiDiagram,
						DiffusionReactionEquation::IMPLICIT, //time
						DiffusionReactionEquation::IMPLICIT,//EXPLICIT, //reaction
						DiffusionReactionEquation::IMPLICIT //diffusion
						);
		morphogen =
				(float*) malloc(
						sizeof(float) * voronoiDiagram->countVoronoiCells);
		morphogenProduction =
				(float*) malloc(
						sizeof(float) * voronoiDiagram->countVoronoiCells);
		morphogenDecay =
				(float*) malloc(
						sizeof(float) * voronoiDiagram->countVoronoiCells);
		morphogenDiffusion =
				(float*) malloc(
						sizeof(float) * voronoiDiagram->countVoronoiCells);

		morphogenDynamics->initSystem(morphogen);
		morphogenDynamics->initReaction(morphogenProduction);
		morphogenDynamics->initReactionDerivative(morphogenDecay);
		morphogenDynamics->initDiffusion(morphogenDiffusion);
		morphogenDynamics->setTimeStep(customTimestep);

		for (int i = 0; i < voronoiDiagram->countVoronoiCells; i++) {
			morphogen[i] = 0.;
			morphogenDiffusion[i] = Lactate_Diffusion;
			morphogenProduction[i] = 0.;//-300;
		}
	}

	double *dECM = 0;
	if (VoronoiCell::ECM_THRESHOLD_QUIESCENCE != 0) {

		dECM = (double*)malloc( voronoiDiagram->countVoronoiCells*sizeof(double));
		for( int v=0; v<voronoiDiagram->countVoronoiCells; v++)
			dECM[v]=0;
	}

	// update reaction
	/*for( int i=0; i<voronoiDiagram->countVoronoiCells; i++){
	 lactateProduction[i] = -GetLactateProductionRate( voronoiDiagram->voronoiCells[i], voronoiDiagram->voronoiCells[i]->glucose, voronoiDiagram->voronoiCells[i]->oxygen);
	 }

	 // update lactate
	 lactateDynamics->update();
	 morphogenDynamics->update( 0.0002);

	 // hard copy
	 for( int i=0; i<voronoiDiagram->countVoronoiCells; i++){
	 voronoiDiagram->voronoiCells[i]->lactate   = lactate[i];
	 voronoiDiagram->voronoiCells[i]->morphogen = morphogen[i];
	 / *fprintf( stderr, "g:%e, o:%e -> (dl:%e, dg:%e, do:%e)\n",
	 voronoiDiagram->voronoiCells[i]->glucose,
	 voronoiDiagram->voronoiCells[i]->oxygen,
	 GetLactateProductionRate( voronoiDiagram->voronoiCells[i], voronoiDiagram->voronoiCells[i]->glucose, voronoiDiagram->voronoiCells[i]->oxygen),
	 GetGlucoseConsumptionRate(voronoiDiagram->voronoiCells[i], voronoiDiagram->voronoiCells[i]->glucose, voronoiDiagram->voronoiCells[i]->oxygen),
	 GetOxygenConsumptionRate( voronoiDiagram->voronoiCells[i], voronoiDiagram->voronoiCells[i]->glucose, voronoiDiagram->voronoiCells[i]->oxygen));* /
	 }

	 writeSliceInformation( voronoiDiagram, agentArray, dirname);
	 exit(0);*/

#else
	Substrate substrate( countPoints, Case, voronoiDiagram, agentArray, InitialOxygenConcentration, InitialGlucoseConcentration, WilliamsCodeParameter, dirname, boundaryCondition, customTimestep, gui);
#endif

//	fprintf(stderr, "...finished ( %lisec)\n",
//			(clock() - passedTime) / CLOCKS_PER_SEC);
#endif

	// INITIALIZE STATISTICAL ARRAYS

	// calculate number of discrete measure points 
	EndTime_int = nrOfSteps( BeginningTime, EndTime, OutputRate);

	// allocate memory for arrays
	global_cells = (double*) calloc(EndTime_int + 1, sizeof(double));
	for (k = 0; k < EndTime_int+1; k++)
		global_cells[k] = 0.0;

	global_necrotic_cells = (double*) calloc(EndTime_int + 1, sizeof(double));
	for (k = 0; k < EndTime_int+1; k++)
		global_necrotic_cells[k] = 0.0;

	global_necrotic_volume = (double*) calloc(EndTime_int + 1, sizeof(double));
	for (k = 0; k < EndTime_int+1; k++)
		global_necrotic_volume[k] = 0.0;

	global_vessel_cells = (double*) calloc(EndTime_int + 1, sizeof(double));
	for (k = 0; k < EndTime_int+1; k++)
		global_vessel_cells[k] = 0.0;

	global_volume = (double*) calloc(EndTime_int + 1, sizeof(double));
	for (k = 0; k < EndTime_int+1; k++)
		global_volume[k] = 0.0;

	global_inner_empty_volume =
			(double*) calloc(EndTime_int + 1, sizeof(double));
	for (k = 0; k < EndTime_int+1; k++)
		global_inner_empty_volume[k] = 0.0;

	double *global_gyrRadius = (double*) calloc(EndTime_int + 1, sizeof(double));
	for (k = 0; k < EndTime_int+1; k++)
		global_gyrRadius[k] = 0.0;

	double *global_gyrRadiusSquare = (double*) calloc(EndTime_int + 1,
			sizeof(double));
	for (k = 0; k < EndTime_int+1; k++)
		global_gyrRadiusSquare[k] = 0.0;

	double *global_avgOxygen = (double*) calloc(EndTime_int + 1, sizeof(double));
	for (k = 0; k < EndTime_int+1; k++)
		global_avgOxygen[k] = 0.0;

	double *global_avgGlucose = (double*) calloc(EndTime_int + 1,
			sizeof(double));
	for (k = 0; k < EndTime_int+1; k++)
		global_avgGlucose[k] = 0.0;

#if PROB_OUTPUT
	int MaxCells = 100000;
	int **global_prob = ( int** ) calloc ( EndTime_int, sizeof( int*));
	for(int t = 0;t < EndTime_int; t++) {
		global_prob[t] = ( int* ) calloc ( MaxCells, sizeof( int*));
		for(int c = 0; c < MaxCells; c++) {
			global_prob[t][c] = 0;
		}
	}
#endif

#if RADIAL_OUTPUT
	int kk;
	//double EndRadius = 60.;
	double EndRadius = (voronoiDiagram->xMax[0]-voronoiDiagram->xMin[0]);
	double OutputRateRadius = 1.;
	int EndRadius_int = (int) (EndRadius/OutputRateRadius);

	double** global_oxygen_concentration = ( double** ) calloc ( EndTime_int, sizeof( double*));
	double** global_glucose_concentration = ( double** ) calloc ( EndTime_int, sizeof( double*));
	double** global_active_cell_density = ( double** ) calloc ( EndTime_int, sizeof( double*));
	double** global_nonactive_cell_density = ( double** ) calloc ( EndTime_int, sizeof( double*));
	double** global_necrotic_cell_density = ( double** ) calloc ( EndTime_int, sizeof( double*));
	for(k = 0;k < EndTime_int; k++) {
		global_oxygen_concentration[k] = ( double* ) calloc ( EndRadius_int, sizeof( double));
		global_glucose_concentration[k] = ( double* ) calloc ( EndRadius_int, sizeof( double));
		global_active_cell_density[k] = ( double* ) calloc ( EndRadius_int, sizeof( double));
		global_nonactive_cell_density[k] = ( double* ) calloc ( EndRadius_int, sizeof( double));
		global_necrotic_cell_density[k] = ( double* ) calloc ( EndRadius_int, sizeof( double));
		for(kk = 0;kk < EndRadius_int; kk++) {
			global_oxygen_concentration[k][kk] = 0.0;
			global_glucose_concentration[k][kk] = 0.0;
			global_active_cell_density[k][kk] = 0.0;
			global_nonactive_cell_density[k][kk] = 0.0;
			global_necrotic_cell_density[k][kk] = 0.0;
		}
	}
#endif

	// WRITE ALL PARAMETERS TO FILE
	sprintf(parameterfilename, "%s%cparameter.dat", dirname, SEPERATOR);
	if ((fp_param = fopen(parameterfilename, "w")) == NULL) {
		fprintf(stderr, "Error opening file %s for writing!\n",
				parameterfilename);
		exit(0);
	}

	fprintf(fp_param, "Command Line:\n");
	for (i = 0; i < argc; i++) {
		fprintf(fp_param, "%s ", argv[i]);
	}
	fprintf(fp_param, "\n");

	fprintf(fp_param, "******************************** \n");
	fprintf(fp_param, "\nDefined Parameters: \n");
	fprintf(fp_param, "#define MIN_SPECIFIC_DEATH_RATE %lf\n",
			MIN_SPECIFIC_DEATH_RATE);
	fprintf(fp_param, "#define MAX_SPECIFIC_DEATH_RATE %lf\n",
			MAX_SPECIFIC_DEATH_RATE);
	fprintf(fp_param, "#define MONOD_PARAMETER_DEATH %lf\n",
			MONOD_PARAMETER_DEATH);
	fprintf(fp_param, "#define MAX_SPECIFIC_GROWTH_RATE %lf\n",
			MAX_SPECIFIC_GROWTH_RATE);
	fprintf(fp_param, "#define MONOD_PARAMETER_GROWTH %lf\n",
			MONOD_PARAMETER_GROWTH);
	fprintf(fp_param, "#define CELL_GROWTH_YIELD_OXYGEN %lf\n",
	CELL_GROWTH_YIELD_OXYGEN);
	fprintf(fp_param, "#define CELL_GROWTH_YIELD_GLUCOSE %lf\n",
	CELL_GROWTH_YIELD_GLUCOSE);
	fprintf(fp_param, "#define MAINTENANCE_ENERGY_REQUIREMENT_OXYGEN %lf\n",
	MAINTENANCE_ENERGY_REQUIREMENT_OXYGEN);
	fprintf(fp_param, "#define MAINTENANCE_ENERGY_REQUIREMENT_GLUCOSE %lf\n",
	MAINTENANCE_ENERGY_REQUIREMENT_GLUCOSE);
	fprintf(fp_param, "#define MONOD_KINETICS %i\n", MONOD_KINETICS);
	fprintf(fp_param, "#define USE_DIVISION %i\n", USE_DIVISION);
	fprintf(fp_param, "#define USE_MIGRATION %i\n", USE_MIGRATION);
	fprintf(fp_param, "#define USE_APOPTOSIS %i\n", USE_APOPTOSIS);
	fprintf(fp_param, "#define USE_NECROSIS %i\n", USE_NECROSIS);
	fprintf(fp_param, "#define USE_LYSIS %i\n", USE_LYSIS);
	fprintf(fp_param, "#define MONOD_KINETICS %i\n", MONOD_KINETICS);
	fprintf(fp_param, "#define MULTISCALE %i\n", MULTISCALE);
	fprintf(fp_param, "#define AGENT_DIAMETER %lf\n", AGENT_DIAMETER);
	fprintf(fp_param, "#define SPATIAL_UNIT %lf\n", SPATIAL_UNIT);

	fclose(fp_param);

	// PRINT ALL PARAMETERS TO SCREEN
#if VERBOSE > 1
	fprintf(stderr,"******************************** \n");
	fprintf(stderr,"\nDefined Parameters: \n");
	fprintf(stderr,"#define MIN_SPECIFIC_DEATH_RATE %lf\n", MIN_SPECIFIC_DEATH_RATE);
	fprintf(stderr,"#define MAX_SPECIFIC_DEATH_RATE %lf\n", MAX_SPECIFIC_DEATH_RATE);
	fprintf(stderr,"#define MONOD_PARAMETER_DEATH %lf\n", MONOD_PARAMETER_DEATH);
	fprintf(stderr,"#define MAX_SPECIFIC_GROWTH_RATE %lf\n", MAX_SPECIFIC_GROWTH_RATE);
	fprintf(stderr,"#define MONOD_PARAMETER_GROWTH %lf\n", MONOD_PARAMETER_GROWTH);
	fprintf(stderr,"#define CELL_GROWTH_YIELD_OXYGEN %lf\n", CELL_GROWTH_YIELD_OXYGEN);
	fprintf(stderr,"#define CELL_GROWTH_YIELD_GLUCOSE %lf\n", CELL_GROWTH_YIELD_GLUCOSE);
	fprintf(stderr,"#define MAINTENANCE_ENERGY_REQUIREMENT_OXYGEN %lf\n", MAINTENANCE_ENERGY_REQUIREMENT_OXYGEN);
	fprintf(stderr,"#define MAINTENANCE_ENERGY_REQUIREMENT_GLUCOSE %lf\n", MAINTENANCE_ENERGY_REQUIREMENT_GLUCOSE);
	fprintf(stderr,"#define MONOD_KINETICS %i\n", MONOD_KINETICS);
	fprintf(stderr,"#define USE_DIVISION %i\n", USE_DIVISION);
	fprintf(stderr,"#define USE_MIGRATION %i\n", USE_MIGRATION);
	fprintf(stderr,"#define USE_APOPTOSIS %i\n", USE_APOPTOSIS);
	fprintf(stderr,"#define USE_NECROSIS %i\n", USE_NECROSIS);
	fprintf(stderr,"#define USE_LYSIS %i\n", USE_LYSIS);
	fprintf(stderr,"#define MONOD_KINETICS %i\n", MONOD_KINETICS);
	fprintf(stderr,"#define MULTISCALE %i\n", MULTISCALE);
	fprintf(stderr,"#define AGENT_DIAMETER %lf\n", AGENT_DIAMETER);
	fprintf(stderr,"#define SPATIAL_UNIT %lf\n", SPATIAL_UNIT);
	fprintf( stderr, "Averages = %d \n\n",Averages);
	fprintf( stderr, "****************************************** \n");
#endif   

	/** STOCHASTIC APPROACHE *************************************************/

	double timeActualization = 0, timeSelect = 0, timeExecution = 0;

#define USE_SPARSE_MATRIX

	double timeStep = customTimestep; //0.1;
	double timeDifference = 0.;
	if (timeStep > timeDifference)
		timeDifference = timeStep;

#ifdef USE_SPARSE_MATRIX

	int N = 2 * voronoiDiagram->xN[0]
#if DIMENSIONS >= 2
			* voronoiDiagram->xN[1]
#endif
#if DIMENSIONS >= 3
	* voronoiDiagram->xN[2]
#endif
	;
	//int N = voronoiDiagram->xN[0]*voronoiDiagram->xN[1]*voronoiDiagram->xN[2];
	//fprintf(stderr, "N=%i\n", N);
	if (Case > 1) {
		sA = SparseMatrix::newSparseMatrix(N, N);
		//M = SparseMatrix::newSparseMatrix( N, N);
		//J = SparseMatrix::newSparseMatrix( N, N);
		b = (float *) malloc(N * sizeof(float));
		x = (float *) malloc(N * sizeof(float));
		v0 = (float *) malloc(N * sizeof(float));
		v1 = (float *) malloc(N * sizeof(float));
		v2 = (float *) malloc(N * sizeof(float));
		v3 = (float *) malloc(N * sizeof(float));
		v4 = (float *) malloc(N * sizeof(float));
		v5 = (float *) malloc(N * sizeof(float));
		v6 = (float *) malloc(N * sizeof(float));
		v7 = (float *) malloc(N * sizeof(float));
		v8 = (float *) malloc(N * sizeof(float));
	}
	/*		for( double t=0; t<timeDifference; t+=timeStep){
	 UpdateSystemImplicitSparse( voronoiDiagram, timeStep, timeStep);
	 writeSliceInformation( voronoiDiagram, agentArray, dirname);
	 }
	 */
#else

	initMatrixImplicit( voronoiDiagram);
	for( double t=0; t<timeDifference; t+=timeStep) {
		UpdateSystemImplicit( voronoiDiagram, timeStep, timeStep);
		writeSliceInformation( voronoiDiagram, agentArray, dirname);
	}

#endif	

	// Interpolation
	Interpolation::voronoiDiagram = voronoiDiagram;

	// PHARMACO KINETICS
//	PharmacoKinetics *pharmaco = new PharmacoKinetics( voronoiDiagram);
	//pharmaco->setBorderConditionCPT11in(0.);
//	pharmaco->setBorderConditionCPT11out( 27.03);
	//pharmaco->setBorderConditionCPT11out(0.);
	//pharmaco->initCPT11in( voronoiDiagram, 11.);
	//pharmaco->initCPT11out( voronoiDiagram, 11.);

	//srand( 0);

	for (k = 0; k < Averages; k++) {

//		pharmaco->initCPT11in( voronoiDiagram, 0.);
//		pharmaco->initCPT11out( voronoiDiagram, 0.);

		srand(seed_start + k);

		Last_Time = Time = TempTime = 0;
		gyrRadius = 0;
		count_inoculated_cells = 0;
		count_divided_cells = 0;
		count_expanded_cells = 0;
		count_cells = 0;
		count_cell_volume = 0;
		count_vessel_cells = 0;
		count_necrotic_cells = 0;
		count_necrotic_volume = 0;
		count_inner_empty_volume = 0;
		count_lysed_cells = 0;
		count_lysed_volume = 0;
		prob_X_equal_k = 0;
		prob_X_smaller_k = 0;
		prob_X_bigger_k = 0;
		//#if REFILL_MEDIUM
		double global_oxygen_concentration = InitialOxygenConcentration;
		double global_glucose_concentration = InitialGlucoseConcentration;
		//#endif

		//fprintf( stderr, "INFO: Start cell added to voronoi point %i: cellCount %i, type %s\n", centralCell->nr, centralCell->cellCount, cellTypeToString( centralCell->state));
		timer = time(NULL);

		// INITIALIZE RANDOM GENERATOR
		//srand( seed_start );
		//srand( timer );
		//srand( 0 );

		// INITIALIZE PROBABILITY LIST

		actionList = ActionTree::newActionTree();
		//ActionTree *actionTree = ActionTree::newActionTree();
		//  test_voronoi_grid( voronoiGrid);
		//initProbabilityElements( voronoiGrid);

		// SET INITIAL VESSEL NETWORK

		//fprintf( stderr, "CENTRAL CELL\n");
		VoronoiCell *centralCell = getCentralCell(voronoiDiagram);
		centralCell =
				centralCell->neighborCells[(int) (centralCell->countNeighborCells
						* myRand())];
		//fprintf( stderr, "central cell: ( %lf, %lf)\n", centralCell->position[0], centralCell->position[1]);
		//exit( 0);

		if( VoronoiCell::USE_SYMBOLIC_EXTENDED_NEIGHBORHOOD){
		// Init Symbolic Extended Neighborhood
			// extended neighborship
			VoronoiCell ***extendedNeighborhood;
			int extendedNeighborhoodSize[(int) CellDivisionDepth + 1];
			int extendedNeighborhoodMax[(int) CellDivisionDepth + 1];
			extendedNeighborhood =
					(VoronoiCell ***) malloc(
							(int) (CellDivisionDepth + 1)
									* sizeof(VoronoiCell **));

			// 1st gen
			extendedNeighborhood[0] =
					(VoronoiCell **) malloc(1 * sizeof(VoronoiCell *));
			extendedNeighborhood[0][0] = centralCell;
			extendedNeighborhoodSize[0] = 1;

			// 2nd gen
			extendedNeighborhood[1] =
					(VoronoiCell **) malloc(
							centralCell->countNeighborCells
									* sizeof(VoronoiCell *));
			for (int n = 0; n < centralCell->countNeighborCells; n++)
				extendedNeighborhood[1][n] = centralCell->neighborCells[n];
			extendedNeighborhoodSize[1] = centralCell->countNeighborCells;

			// following gen
			for (int g = 2; g <= (int) CellDivisionDepth; g++) {
				// alloc memory
				extendedNeighborhood[g] =
						(VoronoiCell **) malloc(1000 * sizeof(VoronoiCell *));
				extendedNeighborhoodSize[g] = 0;
				extendedNeighborhoodMax[g] = 1000;

				//fprintf( stderr, "gen:%i\n", g);

				// look for all points in old generation
				for (int i = 0; i < extendedNeighborhoodSize[g - 1]; i++)
					// look at surrounding points
					for (int n = 0;
							n
									< extendedNeighborhood[g - 1][i]->countNeighborCells;
							n++) {
						// if not added already
						//fprintf( stderr, "TEST cell:%i\n", extendedNeighborhood[g-1][i]->neighborCells[n]->index);
						if (!isElementOf(
								extendedNeighborhood[g - 2],
								extendedNeighborhoodSize[g - 2],
								extendedNeighborhood[g - 1][i]->neighborCells[n])
								&& !isElementOf(
										extendedNeighborhood[g - 1],
										extendedNeighborhoodSize[g - 1],
										extendedNeighborhood[g - 1][i]->neighborCells[n])
								&& !isElementOf(
										extendedNeighborhood[g],
										extendedNeighborhoodSize[g],
										extendedNeighborhood[g - 1][i]->neighborCells[n])) {
							// realloc
							if (extendedNeighborhoodSize[g]
									== extendedNeighborhoodMax[g]) {
								extendedNeighborhoodMax[g] += 1000;
								extendedNeighborhood[g] =
										(VoronoiCell **) realloc(
												extendedNeighborhood[g],
												extendedNeighborhoodMax[g]
														* sizeof(VoronoiCell *));
							}
							// add
							extendedNeighborhood[g][extendedNeighborhoodSize[g]] =
									extendedNeighborhood[g - 1][i]->neighborCells[n];
							extendedNeighborhoodSize[g]++;
						}
					}
			}

			// symbolic extended (with direct neighbors)
			VoronoiCell::symbolicExtendedNeighborhood =
					(int**) malloc((int) (CellDivisionDepth) * sizeof(int *));
			VoronoiCell::symbolicExtendedNeighborhoodSize =
					(int*) malloc((int) (CellDivisionDepth) * sizeof(int));
			for (int i = 0; i < (int) CellDivisionDepth; i++) {
				//fprintf(stderr, "generation %i: ", i);
				VoronoiCell::symbolicExtendedNeighborhoodSize[i] =
						extendedNeighborhoodSize[i + 1];
				VoronoiCell::symbolicExtendedNeighborhood[i] =
						(int*) malloc(
								VoronoiCell::symbolicExtendedNeighborhoodSize[i]
										* sizeof(int));
				for (int j = 0;
						j < VoronoiCell::symbolicExtendedNeighborhoodSize[i];
						j++) {
					/*fprintf(
							stderr,
							" %i,",
							extendedNeighborhood[i + 1][j]->index
									- centralCell->index);*/
					VoronoiCell::symbolicExtendedNeighborhood[i][j] =
							extendedNeighborhood[i + 1][j]->index
									- centralCell->index;
				}
				//fprintf(stderr, "\n");
				free(extendedNeighborhood[i + 1]);
			}
			free(extendedNeighborhood[0]);
			//free( extendedNeighborhood[1]);
			free(extendedNeighborhood);
		}

		// TEST

#ifdef TEST_STATIC_PROFIL
		// set central cells
		int size = 0;
		int N = 3;//CountCellsPerVoronoiCell;
		int k = 5;
		int iy=0, iz=0;
		int index = floor(voronoiDiagram->xN[0]/2)
#if DIMENSIONS>=2
		+ voronoiDiagram->xN[0]*floor(voronoiDiagram->xN[1]/2)
		+ voronoiDiagram->xN[0]*voronoiDiagram->xN[1]*floor(voronoiDiagram->xN[2]/2)
#endif
		;
		fprintf( stderr, "domain (%i, %i, %i)\n", voronoiDiagram->xN[0], voronoiDiagram->xN[1], voronoiDiagram->xN[2]);
		int dx = 1;
		int dy = dx*voronoiDiagram->xN[0];
		int dz = dy*voronoiDiagram->xN[1];
		if( CountCellsPerVoronoiCell==1) {
//			for( int ix=-dx*(1+size*N) - ((N+1)%2); ix<=dx*(1+size*N); ix+=dx)
			for( int ix=-dx*(1+size*N) - ((N+1)%2)+1; ix<=dx*(1+size*N) + 1; ix+=dx)
#if DIMENSIONS >= 2
			for( int iy=-dy*(1+size*N); iy<=dy*(1+size*N); iy+=dy)
#if DIMENSIONS == 3
			for( int iz=-dz*(1+size*N); iz<=dz*(1+size*N); iz+=dz)
#endif
#endif
			{
				int i = index + ix
#if DIMENSIONS >= 2
				+ iy
#endif
#if DIMENSIONS == 3
				+ iz
#endif
				;
				fprintf( stderr, "set cell: %i (%i, %i, %i), (%lf, %lf, %lf)\n", index, ix, iy, iz,
						voronoiDiagram->voronoiCells[i]->position[0], voronoiDiagram->voronoiCells[i]->position[1], voronoiDiagram->voronoiCells[i]->position[2]);
				agentArray->activateAgent()->attach( voronoiDiagram->voronoiCells[i]);
				GetAgent( voronoiDiagram->voronoiCells[i])->growingTumorCellCount = 1;
				GetAgent( voronoiDiagram->voronoiCells[i])->dividingTumorCellCount = 0;
				GetAgent( voronoiDiagram->voronoiCells[i])->cellCount = 1;
				GetAgent( voronoiDiagram->voronoiCells[i])->state = ACTIVE;
			}
		} else {
			for( int ix=-dx*size; ix<=dx*size; ix+=dx)
#if DIMENSIONS>=2
			for( int iy=-dy*size; iy<=dy*size; iy+=dy)
#if DIMENSIONS == 3
			for( int iz=-dz*size; iz<=dz*size; iz+=dz)
#endif
#endif
			{
				int i = index + ix
#if DIMENSIONS>=2
				+ iy
#endif
#if DIMENSIONS==3
				+ iz
#endif
				;
				fprintf( stderr, "set cell: %i (%i, %i, %i), (%lf, %lf, %lf)\n", index, ix, iy, iz,
						voronoiDiagram->voronoiCells[i]->position[0], voronoiDiagram->voronoiCells[i]->position[1], voronoiDiagram->voronoiCells[i]->position[2]);
				agentArray->activateAgent()->attach( voronoiDiagram->voronoiCells[i]);
				GetAgent( voronoiDiagram->voronoiCells[i])->growingTumorCellCount = pow(N,DIMENSIONS);
				GetAgent( voronoiDiagram->voronoiCells[i])->dividingTumorCellCount = 0;
				GetAgent( voronoiDiagram->voronoiCells[i])->cellCount = pow(N,DIMENSIONS);
				GetAgent( voronoiDiagram->voronoiCells[i])->maxCellCount = pow(N,DIMENSIONS);
				GetAgent( voronoiDiagram->voronoiCells[i])->state = COMPARTMENT;
			}
		}

		// set initial conc.
		for( i=0; i<voronoiDiagram->countVoronoiCells; i++) {
			voronoiDiagram->voronoiCells[i]->oxygen = InitialOxygenConcentration;
			voronoiDiagram->voronoiCells[i]->glucose = InitialGlucoseConcentration;
		}

		timeDifference = UpdateSystemNonLinearCGSparse(
				voronoiDiagram,
				0., EndTime,
				customTimestep);

		writeSliceInformation( voronoiDiagram, agentArray, dirname);

		for( i=0; i<agentArray->countActiveAgents; i++) {
			initCellActions( agentArray->agents[i]);
			fprintf( stderr, "Growth rate = %lf, conc = %lf -> neighbors = (%lf <> %lf) -> interpol = %lf -> spline = %lf\n",
					GetGrowthRate( agentArray->agents[i]),
					agentArray->agents[i]->location[0]->glucose * agentArray->agents[i]->location[0]->oxygen,
					agentArray->agents[i]->location[0]->neighborCells[0]->glucose * agentArray->agents[i]->location[0]->neighborCells[0]->oxygen,
					agentArray->agents[i]->location[0]->neighborCells[1]->glucose * agentArray->agents[i]->location[0]->neighborCells[1]->oxygen,
					Interpolation::getDiscritizedVolumeAboveThreshold1D( agentArray->agents[i]->location[0]->index, 0, 0.025),
					Interpolation::getDiscritizedVolumeAboveThreshold1DCubicSpline( agentArray->agents[i]->location[0]->index, 0, 0.025));
		}

		exit( 0);
#endif
		// END TEST

		double DistanceToInitialCellArray[DIMENSIONS];
		//DistanceToInitialCellArray[0] = DistanceToInitialCell;
		//for( i=1; i<DIMENSIONS; i++)
		for (i = 0; i < DIMENSIONS; i++)
			DistanceToInitialCellArray[i] = DistanceToInitialCell; // default value

		switch (Case) {
		case 1:
		case 2:
			break;
		case 3:
			break;
		case 4:
		case 5:
			fprintf(stderr, "Initialize Vessel Network...\n");
			//SetInitialVesselNetwork( voronoiDiagram, agentArray, actionList, centralCell, DistanceToInitialCell, CountVessel, BranchingProbability, BranchingLength);
			SetRegularInitialVesselNetwork(voronoiDiagram, agentArray,
					actionList, centralCell, DistanceToInitialCellArray,
					BranchingProbability, BranchingLength);
			//( VoronoiDiagram *voronoiDiagram, AgentList* agentArray, ActionTree *probTree, VoronoiCell *initialCell, int distanceBetweenVessels, double branchingProbability, int branchingLength)
			fprintf(stderr, "...finished\n");
			break;
		}

#ifndef __WithGUI__
		//substrate.SwitchOnVisualization();
#endif
		substrate.Initialization(Case);

		// SET FIRST CELL
		//fprintf( stderr, "INITIALIZE FIRST CELL\n");

		//fprintf( stderr, "INFO: Start cell added to voronoi point %i: cellCount %i, type %s\n", centralCell->index, GetAgent(centralCell)->cellCount, cellTypeToString( GetAgent(centralCell)->state));

		/*int i = (int)(myRand() * (centralCell->countNeighborCells + 1.) - 1.);
		 if(i >= 0)
		 centralCell = centralCell->neighborCells[i];*/

		/*std::fstream fs;
		 fs.open("testBefore.eps", std::fstream::out);
		 EpsIO::PSwriteHeader( &fs,
		 voronoiDiagram->xMin[0]*10,
		 voronoiDiagram->xMax[0]*10,
		 voronoiDiagram->xMin[1]*10,
		 voronoiDiagram->xMax[1]*10);
		 for( int c=0; c<voronoiDiagram->countVoronoiCells; c++){
		 for( int n=0; n<voronoiDiagram->voronoiCells[c]->countNeighborCells; n++){
		 EpsIO::PSwriteLine( &fs,
		 voronoiDiagram->voronoiCells[c]->position[0]*10,
		 voronoiDiagram->voronoiCells[c]->position[1]*10,
		 //voronoiDiagram->voronoiCells[c]->neighborCells[n]->position[0]*10,
		 //voronoiDiagram->voronoiCells[c]->neighborCells[n]->position[1]*10,
		 (voronoiDiagram->voronoiCells[c]->position[0]+voronoiDiagram->voronoiCells[c]->neighborCells[n]->position[0])/2.*10,
		 (voronoiDiagram->voronoiCells[c]->position[1]+voronoiDiagram->voronoiCells[c]->neighborCells[n]->position[1])/2.*10,
		 0.1, " 0 0 0 setrgbcolor");
		 if( voronoiDiagram->voronoiCells[c]==centralCell)
		 EpsIO::PSwriteCircle( &fs,
		 voronoiDiagram->voronoiCells[c]->position[0]*10,
		 voronoiDiagram->voronoiCells[c]->position[1]*10,
		 1,
		 1,
		 " 0 0 1 setrgbcolor");
		 if( c==voronoiDiagram->countVoronoiCells-1)
		 EpsIO::PSwriteCircle( &fs,
		 voronoiDiagram->voronoiCells[c]->position[0]*10,
		 voronoiDiagram->voronoiCells[c]->position[1]*10,
		 1,
		 1,
		 " 1 0 0 setrgbcolor");
		 }
		 }
		 fs.close();*/
		//printTriangulation( voronoiDiagram, "testBefore.eps");
		for (int f = 0; f < voronoiDiagram->countFramePoints; f++) {
			agentArray->activateAgent()->attach(voronoiDiagram->framePoints[f]);
			GetAgent(voronoiDiagram->framePoints[f])->state = FRAME;
		}

#ifdef REFINE
		//fprintf( stderr, "countVoronoiCells = %i\n", voronoiDiagram->countVoronoiCells);
		voronoiDiagram->refine( centralCell, CountCellsPerVoronoiCell, actionList);
		//fprintf( stderr, "countVoronoiCells = %i\n", voronoiDiagram->countVoronoiCells);
		//TEST

		VoronoiCell* comp = centralCell;
		if( comp->agent == NULL)
		//exit(0);
		{
			//fprintf(stderr, "Add Compartment Agent\n");
			Agent *compAgent = agentArray->activateAgent();
			compAgent->attach(comp);
			compAgent->state = COMPARTMENT;//NONACTIVE;
			compAgent->cellCount = 0;
			compAgent->maxCellCount = CountCellsPerVoronoiCell*CountCellsPerVoronoiCell
#if DIMENSIONS == 3
			*CountCellsPerVoronoiCell
#endif
			;
			compAgent->cellCount = compAgent->maxCellCount;
			compAgent->countFree = compAgent->maxCellCount;
			compAgent->countActive = 0;
			compAgent->countNonactive = 0;
			//for( int n=0; n<comp->countNeighborCells; n++)
			//	comp->neighborCells[n]->countFreeNeighborCells--;
		}
			//#ifdef REFINE
		centralCell = voronoiDiagram->voronoiCells[voronoiDiagram->countVoronoiCells-1];
			//#endif
			//TEST END
		refineSurrounding( voronoiDiagram, centralCell, actionList, agentArray, CountCellsPerVoronoiCell);
		refineNeighborhood( voronoiDiagram, centralCell, actionList, agentArray, CountCellsPerVoronoiCell);

#endif

		/*fs.open("testAfter.eps", std::fstream::out);
		 EpsIO::PSwriteHeader( &fs,
		 voronoiDiagram->xMin[0]*10,
		 voronoiDiagram->xMax[0]*10,
		 voronoiDiagram->xMin[1]*10,
		 voronoiDiagram->xMax[1]*10);
		 for( int c=0; c<voronoiDiagram->countVoronoiCells; c++){
		 for( int n=0; n<voronoiDiagram->voronoiCells[c]->countNeighborCells; n++){
		 EpsIO::PSwriteLine( &fs,
		 voronoiDiagram->voronoiCells[c]->position[0]*10,
		 voronoiDiagram->voronoiCells[c]->position[1]*10,
		 voronoiDiagram->voronoiCells[c]->neighborCells[n]->position[0]*10,
		 voronoiDiagram->voronoiCells[c]->neighborCells[n]->position[1]*10,
		 0.1, " 0 0 0 setrgbcolor");
		 if( voronoiDiagram->voronoiCells[c]==centralCell)
		 EpsIO::PSwriteCircle( &fs,
		 voronoiDiagram->voronoiCells[c]->position[0]*10,
		 voronoiDiagram->voronoiCells[c]->position[1]*10,
		 1,
		 1,
		 " 0 0 1 setrgbcolor");
		 if( c==voronoiDiagram->countVoronoiCells-1)
		 EpsIO::PSwriteCircle( &fs,
		 voronoiDiagram->voronoiCells[c]->position[0]*10,
		 voronoiDiagram->voronoiCells[c]->position[1]*10,
		 1,
		 1,
		 " 1 0 0 setrgbcolor");
		 }
		 }
		 fs.close();*/
		//printTriangulation( voronoiDiagram, "testAfter.eps");
		//voronoiDiagram->coarsen( voronoiDiagram->voronoiCells[voronoiDiagram->countVoronoiCells-1], 10);
		/*fs.open("testAfterCoarsen.eps", std::fstream::out);
		 EpsIO::PSwriteHeader( &fs,
		 voronoiDiagram->xMin[0]*10,
		 voronoiDiagram->xMax[0]*10,
		 voronoiDiagram->xMin[1]*10,
		 voronoiDiagram->xMax[1]*10);
		 for( int c=0; c<voronoiDiagram->countVoronoiCells; c++){
		 for( int n=0; n<voronoiDiagram->voronoiCells[c]->countNeighborCells; n++){
		 EpsIO::PSwriteLine( &fs,
		 voronoiDiagram->voronoiCells[c]->position[0]*10,
		 voronoiDiagram->voronoiCells[c]->position[1]*10,
		 voronoiDiagram->voronoiCells[c]->neighborCells[n]->position[0]*10,
		 voronoiDiagram->voronoiCells[c]->neighborCells[n]->position[1]*10,
		 0.1, " 0 0 0 setrgbcolor");
		 if( voronoiDiagram->voronoiCells[c]==centralCell)
		 EpsIO::PSwriteCircle( &fs,
		 voronoiDiagram->voronoiCells[c]->position[0]*10,
		 voronoiDiagram->voronoiCells[c]->position[1]*10,
		 1,
		 1,
		 " 0 0 1 setrgbcolor");
		 if( c==voronoiDiagram->countVoronoiCells-1)
		 EpsIO::PSwriteCircle( &fs,
		 voronoiDiagram->voronoiCells[c]->position[0]*10,
		 voronoiDiagram->voronoiCells[c]->position[1]*10,
		 1,
		 1,
		 " 1 0 0 setrgbcolor");
		 }
		 }
		 fs.close();*/
		//exit(0);
		//printTriangulation( voronoiDiagram, "testAfterCoarsen.eps");

		// START WITH ONE CELL
		/*if (!centralCell->isFree()) {
		 //fprintf( stderr, "!centralCell->isFree()\n");
		 for (i = 0; i < centralCell->countNeighborCells
		 && !centralCell->neighborCells[i]->isFree(); i++)
		 ;
		 if (i < centralCell->countNeighborCells) {
		 centralCell = centralCell->neighborCells[i];
		 } else {
		 fprintf(
		 stderr,
		 "ERROR: All grid sites in central area occupied by vessels! Can't find space to put first tumor cell!\n");
		 exit(0);
		 }
		 }

		 //fprintf( stderr, "INFO: central cell %i is declared as %s!\n", centralCell->index, cellTypeToString(centralCell->getState()));

		 //centralCell = voronoiDiagram->voronoiCells[0];
		 // set agent
		 agentArray->activateAgent()->attach(centralCell);
		 GetAgent( centralCell)->growingTumorCellCount = 1;
		 GetAgent( centralCell)->dividingTumorCellCount = 0;
		 GetAgent( centralCell)->cellCount = 1;
		 //GetAgent( centralCell)->state = ACTIVE;

		 // update the surrounding cells
		 //fprintf( stderr, "cellcount: %i\n", GetAgent( centralCell)->cellCount);
		 if (CountCellsPerVoronoiCell == 1 || centralCell->refined) {
		 fprintf( stderr, "Set Initial Cell: %i\n", centralCell->index);
		 GetAgent( centralCell)->state = ACTIVE;//COMPARTMENT;
		 update_surrounding_added_cell(actionList, centralCell, voronoiDiagram);
		 } else {
		 fprintf( stderr, "Set Initial Compartment: %i\n", centralCell->index);
		 GetAgent( centralCell)->state = COMPARTMENT;
		 GetAgent( centralCell)->actualize(actionList);
		 }
		 GetAgent( centralCell)->actions[INDEX_GROWTH]->internalStateM[0] = 1;
		 count_inoculated_cells++;
		 count_cells++;
		 count_cell_volume++;
		 */

		//fprintf(stderr, "Init cells\n");
		switch( 0){
		case 0:{
			// START WITH BALL OF CELLS
			//double InitialRadius = 10;
			for (int v = 0; v < voronoiDiagram->countVoronoiCells; v++) {
				double dist = centralCell->getDistanceTo(
						voronoiDiagram->voronoiCells[v]);
				if (dist < InitialRadius) {

					// AGENT
					agentArray->activateAgent()->attach(
							voronoiDiagram->voronoiCells[v]);
					GetAgent( voronoiDiagram->voronoiCells[v])->growingTumorCellCount = 1;
					GetAgent( voronoiDiagram->voronoiCells[v])->dividingTumorCellCount = 0;
					GetAgent( voronoiDiagram->voronoiCells[v])->cellCount = 1;
					if( CountCellsPerVoronoiCell == 1){
						GetAgent( voronoiDiagram->voronoiCells[v])->state = ACTIVE; //COMPARTMENT;
						update_surrounding_added_cell(actionList, voronoiDiagram->voronoiCells[v], voronoiDiagram);
						//initCellActions( GetAgent( voronoiDiagram->voronoiCells[v]));
						//GetAgent( voronoiDiagram->voronoiCells[v])->actualize( actionList);
					}else{
						GetAgent( voronoiDiagram->voronoiCells[v])->state = COMPARTMENT;
						GetAgent( voronoiDiagram->voronoiCells[v])->maxCellCount = CountCellsPerVoronoiCell;
						GetAgent( voronoiDiagram->voronoiCells[v])->actualize( actionList);
						//update_surrounding_added_cell(actionList, voronoiDiagram->voronoiCells[v], voronoiDiagram);
					}


					// STATS
					count_inoculated_cells++;
					count_cells++;
					count_cell_volume++;
				}
			}
		}	break;

		case 1:{
			// START with a DENSITY OF CELLS in a SPHERE
			//double InitialRadius = 10;
			double center[3]={50,50,50};
			for (int v = 0; v < voronoiDiagram->countVoronoiCells; v++)
			if(voronoiDiagram->voronoiCells[v]->getDistanceTo(center)<50.){
				if (0.1 > myRand()) {

					// AGENT
					agentArray->activateAgent()->attach(
							voronoiDiagram->voronoiCells[v]);
					GetAgent( voronoiDiagram->voronoiCells[v])->growingTumorCellCount =
							1;
					GetAgent( voronoiDiagram->voronoiCells[v])->dividingTumorCellCount =
							0;
					GetAgent( voronoiDiagram->voronoiCells[v])->cellCount = 1;
					GetAgent( voronoiDiagram->voronoiCells[v])->state = ACTIVE; //COMPARTMENT;
					update_surrounding_added_cell(actionList,
							voronoiDiagram->voronoiCells[v], voronoiDiagram);

					// STATS
					count_inoculated_cells++;
					count_cells++;
					count_cell_volume++;
				}
			}
		}	break;
		}

		// SET 3/4 of the CELLS QUIESCENT
		for (int a = 0; a < agentArray->countActiveAgents; a++) {
			agentArray->agents[a]->divide = //1;
					(myRand() < (1.-InitialQuiescentFraction) * PROB_REENTERING_CELL_CYCLE( SPATIAL_UNIT * GetDistanceToClosestFreeNeighbor( voronoiDiagram, agentArray->agents[a]->location[0])) ? 1: 0);
		}

		//fprintf(stderr, "...finished\n");

		/*for(int a=0; a<agentArray->countActiveAgents; a++){
		 double dist = agentArray->agents[a]->location[0]->getDistanceTo( voronoiDiagram->searchClosestFreeVoronoiCell( agentArray->agents[a]->location[0], VoronoiCell::extendedNeighborDepth, VoronoiCell::symbolicExtendedNeighborhoodSize, VoronoiCell::symbolicExtendedNeighborhood));
		 agentArray->agents[a]->divide = (myRand() < PROB_REENTERING_CELL_CYCLE( sqrt((double)dist) * SPATIAL_UNIT) ? 1 : 0);
		 }*/

		/*fprintf( stderr, "PRINT CELLS\n");
		 sprintf(outfilename, "%s/cellsInit.pov", dirname);
		 agentArray->printToPovray( outfilename, voronoiDiagram);
		 fprintf( stderr, "FINISHED\n");
		 exit(0);*/

		// FILL WHOLE DOMAIN WITH CELLS
		/*fprintf( stderr, "FILL WHOLE DOMAIN WITH CELLS\n");
		 for(int v=0; v<voronoiDiagram->countVoronoiCells; v++){
		 //fprintf( stderr, "Cell %i\n", v);
		 centralCell = voronoiDiagram->voronoiCells[v];
		 // set agent
		 agentArray->activateAgent()->attach(centralCell);
		 GetAgent( centralCell)->growingTumorCellCount = CountCellsPerVoronoiCell;
		 GetAgent( centralCell)->dividingTumorCellCount = 0;
		 GetAgent( centralCell)->cellCount = CountCellsPerVoronoiCell;
		 //GetAgent( centralCell)->state = ACTIVE;

		 // update the surrounding cells
		 //fprintf( stderr, "cellcount: %i\n", GetAgent( centralCell)->cellCount);
		 if (CountCellsPerVoronoiCell == 1) {
		 GetAgent( centralCell)->state = ACTIVE;//COMPARTMENT;
		 update_surrounding_added_cell(actionList, centralCell,
		 voronoiDiagram);
		 } else {
		 GetAgent( centralCell)->state = COMPARTMENT;
		 GetAgent( centralCell)->actualize(actionList);
		 }
		 //GetAgent( centralCell)->actions[INDEX_GROWTH]->internalStateM[0] = 1;
		 count_inoculated_cells+=CountCellsPerVoronoiCell;
		 count_cells+=CountCellsPerVoronoiCell;
		 count_cell_volume+=CountCellsPerVoronoiCell;
		 }
		 fprintf( stderr, "...finished\n");
		 // END: FILL WHOLE DOMAIN WITH CELLS
		 */

		//////////

		// GROWTH PART
		/*for( i=0; i<centralCell->countNeighborCells; i++)
		 if(centralCell->neighborCells[i]->isFree()){

		 VoronoiCell * newLocation = centralCell->neighborCells[i];

		 // create 2nd daughter cell
		 Agent* daughterCell = agentArray->activateAgent();
		 daughterCell->attach( newLocation);
		 daughterCell->state = ACTIVE;
		 initCellActions( daughterCell);
		 update_surrounding_expanded_cell( actionList, newLocation, voronoiDiagram);

		 daughterCell->actualize( actionList);
		 count_inoculated_cells++;
		 count_cells++;
		 count_cell_volume++;
		 //GetAgent(centralCell)->actualize( actionList);
		 }*/

		//////////

		/*for( i=0; i<centralCell->countNeighborCells; i++)
		 if(centralCell->neighborCells[i]->isFree()){

		 VoronoiCell * newLocation = centralCell->neighborCells[i];

		 // create 2nd daughter cell
		 GetAgent( centralCell)->attach( newLocation);
		 update_surrounding_expanded_cell( actionList, newLocation, voronoiDiagram);

		 GetAgent( centralCell)->actualize( actionList);
		 count_cell_volume++;
		 }*/

		/// SPARSE MATRIX TEST ///
		/*{	
		 int N = 10;

		 SparseMatrix *A = SparseMatrix::newSparseMatrix( N, N);
		 float **B = newMatrix( N, N);
		 float b[N];
		 float x[N];
		 fprintf( stderr, "A->set( 0, 1, 0.1);\n");
		 A->set( 0, 1, 0.1);
		 fprintf( stderr, "A->set( 0, 1, 0.2);\n");
		 A->set( 0, 1, 0.2);
		 fprintf( stderr, "A->set( 0, 0, 0.05);\n");
		 A->set( 0, 0, 0.05);
		 for( i=0; i<N; i++){
		 for( int j=0; j<N; j++)
		 if(j>i-2 && j<i+2)
		 {
		 fprintf( stderr, "A->set( %i, %i, %lf);\n",  i, j, (float)(i*j));
		 A->set( i, j, (float)(i*j));
		 B[i][j] =  (float)(i*j);
		 }
		 b[i] = i;
		 }
		 fprintf( stderr, "====================\n");
		 for( i=0; i<N; i++){
		 for( int j=0; j<N; j++){
		 fprintf( stderr, "%.3lf ", A->get(i, j));
		 }
		 fprintf( stderr, "\n");
		 }
		 fprintf( stderr, "====================\n");
		 for( i=0; i<N; i++){
		 for( int j=0; j<A->sizeA[i]; j++){
		 fprintf( stderr, "%.0lf(%i) ", A->A[i][j], A->JA[i][j]);
		 }
		 fprintf( stderr, "\n");
		 }
		 fprintf( stderr, "====================\n");
		 for( i=0; i<N; i++){
		 for( int j=0; j<N; j++)
		 if(j>i-2 && j<i+2)
		 {
		 fprintf( stderr, "A->set( %i, %i, %lf);\n",  i, j, (float)(i*j));
		 A->set( i, j, (float)(i*j));
		 B[i][j] =  (float)(i*j);
		 }
		 b[i] = i;
		 }
		 for( i=0; i<N; i++){
		 for( int j=0; j<N; j++){
		 fprintf( stderr, "%.3lf ", A->get(i, j));
		 }
		 fprintf( stderr, "\n");
		 }
		 fprintf( stderr, "====================\n");
		 for( i=0; i<N; i++){
		 for( int j=0; j<A->sizeA[i]; j++){
		 fprintf( stderr, "%.0lf(%i) ", A->A[i][j], A->JA[i][j]);
		 }
		 fprintf( stderr, "\n");
		 }

		 matrixVectorProduct( B, b, x, N);
		 for( i=0; i<10; i++){
		 fprintf( stderr, "%lf ", x[i]);
		 }
		 fprintf( stderr, "\n");
		 fprintf( stderr, "====================\n");

		 sparseMatrixVectorProduct( A, b, x, N);
		 for( i=0; i<10; i++){
		 fprintf( stderr, "%lf ", x[i]);
		 }
		 fprintf( stderr, "\n");

		 }
		 exit( 0);*/

		// ADI WILLIAM //////////////////////////
		/*double timeStep = customTimestep;//0.1;
		 double timeDifference = 1.; if(timeStep>timeDifference) timeDifference=timeStep;
		 for( double t=0; t<timeDifference; t+=timeStep){
		 substrate.Update_All( Case, timeStep, actionList);
		 writeSliceInformation( voronoiDiagram, agentArray, dirname);
		 }
		 exit( 0);*/
		// END ADI WILLIAM //////////////////////
		/*int M=4;
		 float A[1][M];
		 float B[M][1];
		 float C[M][M];

		 float D[M][M];

		 for( int m=0; m<M; m++){
		 A[0][m]	= M*myRand();
		 B[m][0] = A[0][m];
		 for( int mm=0; mm<M; mm++)
		 C[m][mm] = M*myRand();
		 }
		 fprintf( stderr, "A = B':\n");
		 for( int m=0; m<M; m++)
		 fprintf( stderr, "%lf ", A[0][m]);
		 fprintf( stderr, "\n=================\n");

		 fprintf( stderr, "C:\n");
		 for( int m=0; m<M; m++){
		 for( int mm=0; mm<M; mm++)
		 fprintf( stderr, "%lf ", C[m][mm]);
		 fprintf( stderr, "\n");
		 }
		 fprintf( stderr, "\n");

		 // C*b = x
		 matrixProduct( C, B, D, M, M, 1);

		 // a*C = x
		 matrixProduct( A, C, D, 1, M, M);

		 exit( 0);	*/

		// IMPLICIT //////////////////////////
		/*timeDifference = 0.;
		 char filename[512];

		 sprintf( filename,"%s%cconc.dat",dirname, SEPERATOR);
		 if(( fp = fopen(filename, "w"))==NULL){
		 fprintf( stderr, "Error opening file %s for writing!\n", filename);
		 exit(0);
		 }

		 for( double t=0; t<EndTime; t+=timeStep){
		 UpdateSystemNonLinearCGSparse( voronoiDiagram, t, t+customTimestep, customTimestep);

		 fprintf( fp, "%lf %e %e\n", t+customTimestep, centralCell->glucose, centralCell->oxygen);
		 }
		 writeSliceInformation( voronoiDiagram, agentArray, dirname);

		 fclose( fp);

		 exit( 0);*/

		// END IMPLICIT ////////////////////
		// EXPLICIT & ADI //////////////////
		/*initMatrix( voronoiDiagram);

		 fprintf( stderr, "(%i, %i, %i)\n", voronoiDiagram->xN[0], voronoiDiagram->xN[1], voronoiDiagram->xN[2]);

		 sprintf(outfilename,"%s/diffusion.dat",dirname);
		 if((fp=fopen(outfilename, "w"))==NULL){
		 fprintf(stderr, "Error opening file %s for writing!\n",outfilename);
		 exit( 0);
		 }

		 double timeStep = customTimestep;
		 double output_interval = 0.001;
		 double timeDifference = 1.;
		 for( double t=0; t<timeDifference; t+=timeStep){
		 UpdateSystem( voronoiDiagram, timeStep, timeStep);
		 //printf( "t=%lf: [O2]=%lf, d[O2] = %lf = %lf * (%lf + %lf)\n", t, centralCell->oxygen, centralCell->doxygen * timeStep, timeStep, - centralCell->oxygen * GiveMeTheOxygenRate( centralCell), Oxygen_Diffusion * secondDerivationOxygen( centralCell));
		 //printf( "t=%lf: [O2]=%lf, d[O2] = %lf\n", t, centralCell->oxygen, centralCell->doxygen * timeStep);
		 if( (int)(t/output_interval) != (int)((t+timeStep)/output_interval)){
		 //printf( "write output\n");
		 writeSliceInformation( voronoiDiagram, agentArray, dirname);
		 fprintf( fp,"%lf %lf %lf %e %e \b\n", t, centralCell->oxygen, centralCell->glucose, centralCell->doxygen, centralCell->dglucose);
		 fprintf( stderr,"t:%lf => [O2]=%lf, d[G]=%lf (d[O2]=%e, d[G]=%e)\n", t, centralCell->oxygen, centralCell->glucose, centralCell->doxygen, centralCell->dglucose);
		 }
		 }
		 //UpdateSystem( voronoiDiagram, timeStep, timeDifference);

		 fclose(fp);

		 exit( 0);*/

		// END EXPLICIT & ADI

		// 2nd cell
		/*do{
		 centralCell = centralCell->neighborCells[(int)(centralCell->countNeighborCells*myRand())];
		 }while( centralCell->getState() != FREE);
		 agentArray->activateAgent()->attach( centralCell);
		 GetAgent( centralCell)->state = NONACTIVE;//COMPARTMENT;
		 update_surrounding_added_cell( actionList, centralCell, voronoiDiagram);

		 for( i=0; i<centralCell->countNeighborCells; i++)
		 if(centralCell->neighborCells[i]->isFree()){

		 VoronoiCell * newLocation = centralCell->neighborCells[i];

		 GetAgent( centralCell)->attach( newLocation);
		 update_surrounding_expanded_cell( actionList, newLocation, voronoiDiagram);

		 GetAgent( centralCell)->actualize( actionList);
		 count_cell_volume++;
		 }*/

		//////////
		//for( i=0; i<centralCell->countNeighborCells; i++){
		/*i=0;
		 if( centralCell->neighborCells[i]->isFree()){
		 // set agent
		 agentArray->activateAgent()->attach( centralCell->neighborCells[i]);
		 GetAgent( centralCell->neighborCells[i])->state = ACTIVE;//COMPARTMENT;
		 //GetAgent( centralCell)->state = COMPARTMENT;
		 GetAgent( centralCell->neighborCells[i])->growingTumorCellCount = 1;
		 GetAgent( centralCell->neighborCells[i])->dividingTumorCellCount = 0;
		 GetAgent( centralCell->neighborCells[i])->cellCount = 1;
		 //GetAgent( centralCell)->state = ACTIVE;

		 // update the surrounding cells
		 //fprintf( stderr, "cellcount: %i\n", GetAgent( centralCell)->cellCount);
		 update_surrounding_added_cell( actionList, centralCell->neighborCells[i], voronoiDiagram);
		 count_inoculated_cells++;
		 }
		 //fprintf( stderr, "cellcount: %i\n", count_inoculated_cells);
		 //}*/

		// initialize all possible processes of the new cell
		//initCellActions( GetAgent( centralCell));
		//fprintf( stderr, "maxCellCount: %i\ncellcount: %i\ntumorCellCount: %i\ntumorVolume: %i\nnecroticCellCount: %i\n", GetAgent( centralCell)->maxCellCount, GetAgent( centralCell)->cellCount, GetAgent( centralCell)->tumorCellCount, GetAgent( centralCell)->tumorVolume, GetAgent( centralCell)->necroticCellCount);
		//actionList->getDepth( actionList->root);
		// ACTION TREE
		/*Action* actualAction[10000];
		 for( i=0; i<10000; i++){
		 actualAction[i] = newAction( i, i*0.25);
		 }

		 for( i=0; i<100000; i++){
		 int whichDelete = (int)(10000.*myRand());
		 //fprintf( stderr, "Delete %i\n", whichDelete);
		 actionTree->deleteAction( actualAction[whichDelete]);
		 //actionTree->getDepth( actionTree->root);
		 //actionTree->print();
		 //fprintf( stderr, "\rDELETE (%i): Size = %3i, Depth = %3i (>= %.3lf) \n",
		 //       whichDelete, actionTree->size, actionTree->getDepth( actionTree->root), log( actionTree->size-1.)/log(2) - 1);


		 int whichAdd = (int)(10000.*myRand());
		 //fprintf( stderr, "Add %i\n", whichAdd);
		 actionTree->addAction( actualAction[whichAdd]);
		 //actionTree->getDepth( actionTree->root);
		 //actionTree->print();
		 fprintf( stderr, "\rADD (%6i):    Size = %3i, Depth = %3i (>= %.3lf), Rate = %5.3lf \r",
		 whichAdd, actionTree->size, actionTree->getDepth( actionTree->root), log( actionTree->size-1.)/log(2) - 1, actionTree->rateSum);
		 //

		 actionTree->actualizeRate( actualAction[whichAdd], actualAction[whichAdd]->rate*myRand()*2.);
		 }
		 fprintf( stderr, "Depth = %i\n", actionTree->getDepth( actionTree->root));
		 */
		//actionList->print();
		//global_count_inoculated_cells++;
		//printToPovray( "anim_0.pov", voronoiDiagram, startAngle, setoffX, setoffY, setoffZ);
		//fprintf( stderr, "FIRST CELL INITIALIZED\n");
		//fprintf( stderr, ":::: %lf, %i (%lf, %lf)\n", GetAgent( centralCell)->oxygen, actionList->length, actionList->head->rate, actionList->head->next->rate);

		// SET ALL CELLS
		/*for( i=0; i<voronoiDiagram->countVoronoiCells; i++){
		 if( GetAgent( voronoiDiagram->voronoiCells[i])->state != VESSEL){
		 // update the surrounding cells
		 update_surrounding_added_cell( actionList, GetAgent( voronoiDiagram->voronoiCells[i]), voronoiDiagram);
		 GetAgent( voronoiDiagram->voronoiCells[i])->timeOfCreation = 0.;
		 GetAgent( voronoiDiagram->voronoiCells[i])->tumorCellCount = 1;
		 count_inoculated_cells++;
		 global_count_inoculated_cells++;
		 }
		 }*/

		// SET FIRST CELLS
		/*for( i=0; i<countPointsX; i++){
		 // initialize all possible processes of the new cell
		 if( !GetAgent( voronoiDiagram->voronoiCells[i])->actionsInitialized)
		 initCellActions( GetAgent( voronoiDiagram->voronoiCells[i]));

		 // update the surrounding cells
		 update_surrounding_added_cell( actionList, GetAgent( voronoiDiagram->voronoiCells[i]), voronoiDiagram);
		 GetAgent( voronoiDiagram->voronoiCells[i])->tumorCellCount = (int) sqrt( CountCellsPerVoronoiCell); //1;
		 GetAgent( voronoiDiagram->voronoiCells[i])->timeOfCreation = 0.;
		 count_inoculated_cells += (int) sqrt( CountCellsPerVoronoiCell); //++;
		 global_count_inoculated_cells += (int) sqrt( CountCellsPerVoronoiCell); //++;
		 }*/

		// INITIALIZE STATISTICAL VARIABLES

		timer = time(NULL);
		Last_Time = Time = 0.0;
		//gyrRadius = getGyrationRadius(agentArray);
		Last_gyrRadius = gyrRadius = getGyrationRadiusOfBorder(agentArray, voronoiDiagram);

		l = (int) Time;
		global_cells[l] += count_cells;
		global_volume[l] += count_cell_volume;

		global_necrotic_cells[l] += count_necrotic_cells;
		global_necrotic_volume[l] += count_necrotic_volume;

		global_inner_empty_volume[l] += count_inner_empty_volume;

		global_vessel_cells[l] += count_vessel_cells;

		global_gyrRadius[l] += sqrt(gyrRadius);

		global_gyrRadiusSquare[l] += gyrRadius;

		global_avgGlucose[l] += getAvgGlucose(agentArray);
		global_avgOxygen[l] += getAvgOxygen(agentArray);

		//central_glucose[l] +=
#if RADIAL_OUTPUT
		statisticsOfSpheroid( voronoiDiagram, global_oxygen_concentration, global_glucose_concentration, global_active_cell_density, global_nonactive_cell_density, global_necrotic_cell_density, EndRadius, EndRadius_int, l);
		statisticsOfSpheroid( voronoiDiagram, global_oxygen_concentration, global_glucose_concentration, global_active_cell_density, global_nonactive_cell_density, global_necrotic_cell_density, EndRadius, EndRadius_int, l);
#endif

		// IMMIDIATE OUTPUT

#if ON_THE_FLY_OUTPUT

		sprintf(outfilename, "%s/rel%i.data.dat", dirname, k);
		fp = fopen(outfilename, "w+");
		fclose(fp);

		sprintf(outfilename, "%s/dataOnTheFly.dat", dirname);
		if ((fp = fopen(outfilename, "w")) == NULL) {
			fprintf(stderr, "Error opening file %s for writing!\n",
					outfilename);
			//exit(0);
		}
		fprintf(fp, "# %i realisations\n", k + 1);
		fprintf(fp,	"# rows: time, tumor cells, necrotic cells, vessel segments\n");
		fclose(fp);
#endif
		// END IMMIDIATE OUTPUT

		actionList->addAction( newAction( GAP, 1./OutputRate));
		//actionList->addAction( newAction( GAP, 1./TIME_STEP));

		/** SIMULATION **/

		for (i = 0; i < voronoiDiagram->countVoronoiCells; i++)
			voronoiDiagram->voronoiCells[i]->ecm = 0.;

		Last_Time = -1.;

		int flag_actualizeProcessRates = TRUE;

		do {

			/*VoronoiCell* closestFreeCell = voronoiDiagram->searchClosestFreeVoronoiCell( centralCell);
			 if( closestFreeCell==0){
			 fprintf( stderr, "WARNING: No free Voronoi cell found\n");
			 exit(0);
			 }
			 else{
			 fprintf( stderr, "INFO: Clostest Voronoi cell to cell %i is cell %i\n", centralCell->index, closestFreeCell->index);
			 }*/

			//printf( "count_cells: %i == %i (active agents)\n", count_cells, agentArray->countActiveAgents);
			// INITIALIZE STATISTICAL VARIABLES
			Last_Time = Time;
			Last_NumberOfCells =
					count_inoculated_cells + (int) count_divided_cells
							- count_necrotic_cells - count_lysed_cells;
			//printf("Time = %lf, Last_Time = %lf, actionList->length = %d, count_divided_cells = %d\n", Time, Last_Time, actionList->length, count_divided_cells);

			// =================== TESTS ==================

			// STEADY STATE & TIME STEP ANALYSIS

			/*double allowed_error = 1e-10,
			 error_oxygen = allowed_error,
			 error_glucose = allowed_error;
			 double t_End = 1,
			 t_Now = 0.,
			 t_Step = TIME_STEP;
			 int steps = 100;
			 #ifdef SCHALLER
			 sprintf( outfilename, "schaller_ts%lf_%i.dat", TIME_STEP, Case);
			 #else
			 sprintf( outfilename, "adi_ts%lf_%i.dat", TIME_STEP, Case);
			 #endif
			 if((fp=fopen(outfilename, "w+"))==NULL){
			 fprintf(stderr, "Error opening file %s for writing!\n",outfilename);
			 //exit(0);
			 }*/

			// REACHE STEADY STATE FOR DIFFERENT TIME STEPS
			//substrate.timestep
			/*int sign = 0;
			 for(int j=0; j<10; j++){
			 fprintf( stderr, "time step %lf\n", substrate.timestep);
			 for( i=0; i<1000 ; i++){
			 double old_oxygen  = centralCell->oxygen;
			 double old_glucose = centralCell->glucose;
			 //t_Now += t_End/(double)steps;
			 t_Now += substrate.timestep;
			 substrate.Update_All( Case ,substrate.timestep);
			 fprintf( fp, "%lf %e %e %e %e\n", t_Now, centralCell->glucose, centralCell->oxygen, fabs(old_oxygen - centralCell->oxygen), fabs(old_glucose - centralCell->glucose));
			 }
			 if( sign == 1){
			 substrate.timestep /=2.;
			 //TIME_STEP /= 2.;
			 sign--;
			 }else{
			 substrate.timestep *=2.;
			 //TIME_STEP *= 2.;
			 sign++;
			 }
			 substrate.RebuiltMatrix();
			 }
			 exit( 0);*/

			//fprintf( fp, "%lf %e %e %e %e\n", t_Now, centralCell->glucose, centralCell->oxygen, 0., 0. );
			/*double steadyStateOxygen  = centralCell->oxygen;
			 double steadyStateGlucose = centralCell->glucose;
			 do{
			 steadyStateOxygen  = centralCell->oxygen;
			 steadyStateGlucose = centralCell->glucose;
			 fprintf( stderr, "steady state: %lf, %lf\n", steadyStateOxygen, steadyStateGlucose);

			 //double old_oxygen  = centralCell->oxygen;
			 //double old_glucose = centralCell->glucose;
			 // look for steady state
			 error_oxygen = allowed_error;
			 error_glucose = allowed_error;
			 for( i=0; i<steps && (error_oxygen>=allowed_error || error_glucose>=allowed_error); i++){
			 double old_oxygen  = centralCell->oxygen;
			 double old_glucose = centralCell->glucose;
			 //t_Now += t_End/(double)steps;
			 t_Now += substrate.timestep;
			 substrate.Update_All( Case ,substrate.timestep);
			 //fprintf( fp, "%lf %e %e %e %e\n", t_Now, centralCell->glucose, centralCell->oxygen, fabs(old_oxygen - centralCell->oxygen), fabs(old_glucose - centralCell->glucose));

			 // ozillation check
			 / *if( i && fabs(old_oxygen - centralCell->oxygen) > error_oxygen){
			 fprintf( stderr, "OZILLATION DETECTED: %lf\n", t_Now);
			 if( fabs(old_oxygen - centralCell->oxygen) - error_oxygen > allowed_error){
			 fprintf( stderr, "Ozillation amplitude %e is bigger than allowed error %e\nChose smaller time step: %e\n", fabs(old_oxygen - centralCell->oxygen) - error_oxygen, allowed_error, substrate.timestep / 2.);
			 substrate.timestep /= 2.;
			 substrate.RebuiltMatrix();
			 //fclose( fp);
			 //exit( 0);
			 }
			 }* /
			 error_oxygen = fabs(old_oxygen - centralCell->oxygen);
			 error_glucose = fabs(old_glucose - centralCell->glucose);
			 }

			 fprintf( fp, "%e %lf %lf\n", substrate.timestep, centralCell->glucose, centralCell->oxygen);

			 // decrease time step
			 substrate.timestep /= 2.;
			 substrate.RebuildMatrix();

			 }while( fabs(steadyStateOxygen - centralCell->oxygen)>allowed_error || fabs(steadyStateGlucose - centralCell->glucose)>allowed_error );

			 // increase time step
			 substrate.timestep *= 2.;
			 fprintf( stderr, "Final Time Step: %e\n", substrate.timestep);

			 fclose( fp);
			 exit( 0);
			 */

			/*for( i=0; i<agentArray->countActiveAgents; i++){
			 int countFreeNeighborCells = 0;
			 for( int ii=0; ii<agentArray->agents[i]->countLocations ; ii++)
			 countFreeNeighborCells += agentArray->agents[i]->location[ii]->countFreeNeighborCells + agentArray->agents[i]->location[ii]->countFreeExtendedNeighborCells;
			 if( countFreeNeighborCells==0 && agentArray->agents[i]->state==ACTIVE ){
			 fprintf( stderr, " agent %i shouldn't be active\n", agentArray->agents[i]->index);
			 exit( 0);
			 }
			 if( countFreeNeighborCells>0 && agentArray->agents[i]->state==NONACTIVE ){
			 fprintf( stderr, " agent %i should be active\n", agentArray->agents[i]->index);
			 exit( 0);
			 }
			 }*/

			// AGENTS CONNECTED?
			/*for( i=0; i<agentArray->countActiveAgents; i++){
			 if( !agentArray->agents[i]->isConnected()){
			 fprintf( stderr, "locations of agent %i are not connected!", agentArray->agents[i]->index);
			 for( int ii=0; ii<agentArray->agents[i]->countLocations; ii++)
			 fprintf( stderr, " %i,", agentArray->agents[i]->location[ii]->index);
			 fprintf( stderr, "\b\n");

			 }
			 }*/
			//fprintf( stderr, "\n");
			// END
			// AGENTS WELL ATTACHED
			/*for( i=0; i<agentArray->countActiveAgents; i++){
			 int ii;
			 for( ii=0; ii<agentArray->agents[i]->countLocations; ii++){
			 if( GetAgent( agentArray->agents[i]->location[ii]) != agentArray->agents[i]){
			 fprintf( stderr, "agent %i isn't well attached to Voronoi cell %i!\n", agentArray->agents[i]->index, agentArray->agents[i]->location[ii]->index);
			 }
			 }
			 }*/

			// END
			/*fprintf( stderr, "\n");
			 for( i=0; i<agentArray->countActiveAgents; i++){
			 fprintf( stderr, "Agent %i (%s) ---[attached to]---> location{ ", agentArray->agents[i]->index, cellTypeToString( agentArray->agents[i]->state));
			 int ii;
			 for(ii=0; ii<agentArray->agents[i]->countLocations; ii++){
			 fprintf( stderr, "%i |", agentArray->agents[i]->location[ii]->index);
			 }
			 fprintf( stderr, "\b}\n");
			 }

			 // TEST VORONOI CELLS
			 for( i=0; i<voronoiDiagram->countVoronoiCells; i++)
			 voronoiDiagram->voronoiCells[i]->validate();*/

			// TEST AGENTS
			/*for( i=0; i<voronoiDiagram->countVoronoiCells; i++)
			 if( voronoiDiagram->voronoiCells[i]->countNeighborCells != 0){
			 //fprintf( stderr, "Test Agent %i\n", agentArray->agents[i]->index);
			 //agentArray->agents[i]->validate();
			 //	if( agentArray->agents[i]->countLocations>2 || agentArray->agents[i]->countLocations<0)
			 //		exit( 0);
			 //	if( agentArray->agents[i]->actions[INDEX_DIVISION]->prev!=NULL && agentArray->agents[i]->actions[INDEX_GROWTH]->prev!=NULL)
			 //		exit(  0);
			 / *	if( agentArray->agents[i]->cellCount > agentArray->agents[i]->maxCellCount ||
			 agentArray->agents[i]->tumorCellCount > agentArray->agents[i]->cellCount ||
			 agentArray->agents[i]->tumorVolume > agentArray->agents[i]->cellCount
			 ){
			 fprintf(stderr,"Agent %i:\n", agentArray->agents[i]->index);
			 fprintf(stderr,"-> (cell:%i, maxCell:%i,tumorCell:%i,tumorVolume%i,necroticCell:%i,freeNeighbors:%i/%i)\n",
			 agentArray->agents[i]->cellCount,
			 agentArray->agents[i]->maxCellCount,
			 agentArray->agents[i]->tumorCellCount,
			 agentArray->agents[i]->tumorVolume,
			 agentArray->agents[i]->necroticCellCount,
			 GetVoronoiCell(agentArray->agents[i])->countFreeNeighborCells, GetVoronoiCell(agentArray->agents[i])->countNeighborCells);

			 exit(0);
			 }
			 if( !agentArray->agents[i]->isGrowing() && agentArray->agents[i]->actions[INDEX_GROWTH]->prev!=NULL){
			 fprintf(stderr,"Agent %i should not be growing!\n", agentArray->agents[i]->index);
			 exit(0);
			 }* /
			 int countFreeNeighbors = 0 ,ii;
			 for( ii=0; ii<voronoiDiagram->voronoiCells[i]->countNeighborCells; ii++){
			 //fprintf( stderr, "Neighbor %i/%i\n", ii+1, GetVoronoiCell(agentArray->agents[i])->countNeighborCells);
			 if( voronoiDiagram->voronoiCells[i]->neighborCells[ii]->isFree())
			 countFreeNeighbors++;
			 else{
			 //fprintf( stderr, "Neighbor Agent %i is not free\n", GetNeighborAgent(agentArray->agents[i], ii)->index);
			 }
			 }
			 if( countFreeNeighbors!=voronoiDiagram->voronoiCells[i]->countFreeNeighborCells){
			 fprintf( stderr, "countFreeNeighbors(%i) != voronoiDiagram->voronoiCells[%i]->countFreeNeighborCells(%i/%i)\n", countFreeNeighbors, i, voronoiDiagram->voronoiCells[i]->countFreeNeighborCells, voronoiDiagram->voronoiCells[i]->countNeighborCells);
			 exit( 0);
			 }

			 countFreeNeighbors = 0;
			 for( ii=0; ii<voronoiDiagram->voronoiCells[i]->countExtendedNeighborCells; ii++){
			 //fprintf( stderr, "Neighbor %i/%i\n", ii+1, GetVoronoiCell(agentArray->agents[i])->countNeighborCells);
			 if( voronoiDiagram->voronoiCells[i]->extendedNeighborhood[ii]->isFree())
			 countFreeNeighbors++;
			 else{
			 //fprintf( stderr, "Neighbor Agent %i is not free\n", GetNeighborAgent(agentArray->agents[i], ii)->index);
			 }
			 }
			 if( countFreeNeighbors!=voronoiDiagram->voronoiCells[i]->countFreeExtendedNeighborCells){
			 fprintf( stderr, "countFreeExtendedNeighborCells(%i) != voronoiDiagram->voronoiCells[%i]->countFreeExtendedNeighborCells(%i/%i)\n", countFreeNeighbors, i, voronoiDiagram->voronoiCells[i]->countFreeExtendedNeighborCells, voronoiDiagram->voronoiCells[i]->countExtendedNeighborCells);
			 for( i=0; i<voronoiDiagram->countVoronoiCells; i++){
			 fprintf( stderr, "cell:%i (%s), %i neighbors [", voronoiDiagram->voronoiCells[i]->index, cellTypeToString(voronoiDiagram->voronoiCells[i]->getState()), voronoiDiagram->voronoiCells[i]->countNeighborCells);
			 for( int ii=0; ii<voronoiDiagram->voronoiCells[i]->countNeighborCells; ii++)
			 fprintf( stderr, " %i,", voronoiDiagram->voronoiCells[i]->neighborCells[ii]->index );
			 fprintf( stderr, "\b], %i extended neighbors [", voronoiDiagram->voronoiCells[i]->countExtendedNeighborCells);
			 for( int ii=0; ii<voronoiDiagram->voronoiCells[i]->countExtendedNeighborCells; ii++)
			 fprintf( stderr, " %i,", voronoiDiagram->voronoiCells[i]->extendedNeighborhood[ii]->index );
			 fprintf( stderr, "\b]\n");
			 }

			 exit( 0);
			 }


			 }*/

			//fprintf( stderr, "voronoiDiagram->voronoiCells[%i]->countFreeNeighborCells(%i/%i)\n", 37, voronoiDiagram->voronoiCells[37]->countFreeExtendedNeighborCells, voronoiDiagram->voronoiCells[i]->countExtendedNeighborCells);

			//fprintf( stderr, "VoronoiCell %i, state: %s, countFreeExtendedNeighbors=%i\n", voronoiDiagram->voronoiCells[2169]->index, cellTypeToString( voronoiDiagram->voronoiCells[2169]->getState()), voronoiDiagram->voronoiCells[2169]->countFreeExtendedNeighborCells);
			//fprintf();
			// =================== TESTS ==================

			if (flag_actualizeProcessRates == TRUE) {
				//fprintf( stderr, "actualizeProcessRates\n");
				actionList->actualizeAllRates(voronoiDiagram);
				flag_actualizeProcessRates = FALSE;
				//actionList->getDepth( actionList->root);
			}

			// SELECT ACTION

			//fprintf( stderr, "actionList->size = %i, actionList->rateSum = %lf\n", actionList->size, actionList->rateSum);
			//actionList->print();
			if (actionList->size != 0 && actionList->rateSum != 0.) {

				bool nullEvent = false;

				//actionList->getDepth( actionList->root);
				//fprintf(stderr, "actualizeProcessRates();\n");
				passedTime = clock();
				//if( Time==0 || indexOfTime( Last_Time, BeginningTime, OutputRate) < indexOfTime( Time, BeginningTime, OutputRate))
				//actionList->print();
				//actionList->print();
				timeActualization +=
						(double) (clock() - passedTime)
								/ (double) CLOCKS_PER_SEC;

				//printf("selectAction();\n");
				passedTime = clock();
				selected_action = actionList->selectAction(&Time);
				//selected_action = selectActionFixTimeStep( agentArray, Time);
				//selected_action = selectActionVariablTimeStep( agentArray, Time);
				timeSelect +=
						(double) (clock() - passedTime)
								/ (double) CLOCKS_PER_SEC;
				//actionList->getDepth( actionList->root);
				//print_prob_list( actionList);
				//actionList->print();
				//printf("selected_action->type = %d, actionList->length = %d (rate sum = %lf)\n", selected_action->type, actionList->size, actionList->rateSum);

				//fprintf(stderr,"Perform %s (rate=%lf, cell nr=%i)\n", actionTypeToString( selected_action->type), selected_action->rate, GetIndex( selected_action->originalCell));
				//print_prob_list( actionList);

				/*fprintf(stderr,"Agent: c(%i) = g(%i) + d(%i) + n(%i)\n", selected_action->originalCell->cellCount, selected_action->originalCell->growingTumorCellCount, selected_action->originalCell->dividingTumorCellCount, selected_action->originalCell->necroticCellCount);
				 if( selected_action->originalCell->cellCount != selected_action->originalCell->growingTumorCellCount + selected_action->originalCell->dividingTumorCellCount + selected_action->originalCell->necroticCellCount)exit(0);*/

				// EXECUTE ACTION
				//printf("EXECUTE ACTION\n");
				if (selected_action->rate != 0 && Time <= EndTime) {
					passedTime = clock();

					switch (selected_action->type) {

					// PROCESSES ON THE MICRO CARRIER SURFACE

#if USE_MIGRATION
					case MIGRATION:
						//fprintf(stderr, "MIGRATION\n");
						//performMigration( voronoiDiagram, agentArray, actionList, selected_action, &myStatistics, Time, EndTime);
						if (0.5 < myRand()
								&& GetVoronoiCell(selected_action->originalCell)->isDomainBorder(
										voronoiDiagram)) {
							//selected_action->originalCell->cellCount = 0;
							//selected_action->originalCell->growingTumorCellCount = 0;
							//selected_action->originalCell->dividingTumorCellCount = 0;
							//selected_action->originalCell->necroticCellCount = 0;
							selected_action->originalCell->state = FREE;
							selected_action->originalCell->actualize(
									actionList);
							selected_action->originalCell->detach();
							fprintf(stderr,
									"WARNING: Agent migrates out of domain!\n");
						} else{


							if( Agent::USE_GRAVITY && Time < 0.2/3600 /* 0.1 sec */){
							//performChemotacticMigration
								//performFreeMigration
								performMigration(voronoiDiagram, agentArray,
										actionList, selected_action,
										&myStatistics, Time, EndTime);

								// UPDATE MIGRATION RATE
								double rho_cell = 1035; // kg/m^3
								double rho_liquid = 992.2; // kg/m^3
								double Cd = 0.1;
								double g = 9.81; // m/s^2
								double d = 20 * pow(10, -6); // m
								double k1 = g * (1 - rho_liquid / rho_cell);
								double k2 = Cd * 6 / 8 / d * rho_liquid
										/ rho_cell;
								double velocity = sqrt(k1 / k2)
										* tanh(Time * 60 * 60 * sqrt(k1 * k2)); // m/s
#define pi 3.14159265
								double lattice_constant = d
										* pow(pi / 6., 1. / 3.); // m
								CellMigrationRate =
										velocity / (pi / 4. * lattice_constant)
												* 60 * 60; // h^-1
								//fprintf(stderr, "Velocity: %e [m/s] => Hopping Rate: %e [mu m/h]\n", velocity, CellMigrationRate);
								flag_actualizeProcessRates = TRUE;
							}else
							{
								if( Agent::USE_GRAVITY){
									alpha_ref = 10000;
									CellMigrationRate = 1.17;
									Agent::USE_GRAVITY = false;
									flag_actualizeProcessRates = TRUE;
								}

								//performChemotacticMigration(
								//performFreeMigration
								performMigration(
									voronoiDiagram, agentArray,
									actionList, selected_action,
									&myStatistics, Time, EndTime);

							}
						}
						//diffuse_cell(actionList, selected_action, voronoiDiagram);
						break;
#endif

						/*case GROWTH:
						 //fprintf(stderr,"GROWTH\n");
						 growCell( selected_action->originalCell);
						 break;*/

#if USE_GROWTH
					case GROWTH: {

						// TODO: Correct Refinement Errors!
						//fprintf( stderr, "GROWTH of agent %i on cell %i\n", selected_action->originalCell->index, GetVoronoiCell( selected_action->originalCell)->index);

						// REFINE NEIGHBORS
#ifdef REFINE/*
						int index = floor(GetVoronoiCell( selected_action->originalCell)->position[0])*voronoiDiagram->xN[1] + floor(GetVoronoiCell( selected_action->originalCell)->position[1]);
						int x = (int)floor(GetVoronoiCell( selected_action->originalCell)->position[0]);
						int y = (int)floor(GetVoronoiCell( selected_action->originalCell)->position[1]);
	#if DIMENSIONS == 3
						int z = (int)floor(GetVoronoiCell( selected_action->originalCell)->position[2]);
						//int dz = 1;
						int dy = voronoiDiagram->xN[2];
						int dx = voronoiDiagram->xN[1] * dy;
						int minx = MAX(0,x-1);
						int miny = MAX(0,y-1);
						int minz = MAX(0,z-1);
						int maxx = MIN(voronoiDiagram->xN[0]-1,x+1);
						int maxy = MIN(voronoiDiagram->xN[1]-1,y+1);
						int maxz = MIN(voronoiDiagram->xN[2]-1,z+1);
	#elif DIMENSIONS == 2
						int dy = 1;
						int dx = voronoiDiagram->xN[1];
	#endif
						int n=0;
						//for( int ix=MAX(0,x-1); ix<MIN(voronoiDiagram->xN[0],x+2); ix++)
						//for( int iy=MAX(0,y-1); iy<MIN(voronoiDiagram->xN[1],y+2); iy++)
						for( int ix=minx; ix<=maxx; ix++)
						for( int iy=miny; iy<=maxy; iy++)
	#if DIMENSIONS == 2
						{
							index = ix*dx +iy;
	#elif DIMENSIONS == 3
						for( int iz=minz; iz<=maxz; iz++)
						//for( int iz=MAX(0,z-1); iz<MIN(voronoiDiagram->xN[2],z+2); iz++)
						{
							index = ix*dx +iy*dy +iz;
	#endif
						//{
							//index = ix*dx +iy*dy;
							if( voronoiDiagram->voronoiCells[index]->refined == false
								&& voronoiDiagram->voronoiCells[index]->isFree()){

								passedTime = clock();
								VoronoiCell* comp = voronoiDiagram->voronoiCells[index];
								voronoiDiagram->refine( voronoiDiagram->voronoiCells[index], CountCellsPerVoronoiCell, actionList);
								fprintf(stderr, "...finished ( %.3lfsec)\n", (double)(clock() - passedTime)/(double)CLOCKS_PER_SEC);

						/ *int n=0;
						for( ; n<GetVoronoiCell( selected_action->originalCell)->countNeighborCells; )
						{
							if( GetVoronoiCell( selected_action->originalCell)->neighborCells[n]->index >= 0 &&
									GetVoronoiCell( selected_action->originalCell)->neighborCells[n]->refined == false
									&& GetVoronoiCell( selected_action->originalCell)->neighborCells[n]->isFree()){
								//printTriangulation( voronoiDiagram, "beforeRef.eps", 1);
								passedTime = clock();
								VoronoiCell* comp = GetVoronoiCell( selected_action->originalCell)->neighborCells[n];
								voronoiDiagram->refine( GetVoronoiCell( selected_action->originalCell)->neighborCells[n], CountCellsPerVoronoiCell, actionList);
								fprintf(stderr, "...finished ( %.3lfsec)\n", (double)(clock() - passedTime)/(double)CLOCKS_PER_SEC);
								//printTriangulation( voronoiDiagram, "afterRef.eps", 1);
								//voronoiDiagram->coarsen( GetVoronoiCell( selected_action->originalCell)->neighborCells[n], 10);
								n=0;* /

								//TEST
								if( comp->agent == NULL)
									//exit(0);
								{
									//fprintf(stderr, "Add Compartment Agent\n");

									Agent *compAgent = agentArray->activateAgent();
									compAgent->attach(comp);
									compAgent->state = COMPARTMENT;//NONACTIVE;
									compAgent->cellCount = 0;
									compAgent->maxCellCount = CountCellsPerVoronoiCell*CountCellsPerVoronoiCell;
									compAgent->cellCount = compAgent->maxCellCount;
									compAgent->countFree = compAgent->maxCellCount;
									compAgent->countActive = 0;
									compAgent->countNonactive = 0;
									//for( int n=0; n<comp->countNeighborCells; n++)
									//	comp->neighborCells[n]->countFreeNeighborCells--;
								}
								//TEST END
							}else{
								n++;
							}
						}
*/
						refineNeighborhood( voronoiDiagram, GetVoronoiCell(selected_action->originalCell), actionList, agentArray, CountCellsPerVoronoiCell);

						if( selected_action->originalCell->state != ACTIVE) {
							fprintf(stderr, "WARNING: Cell Division stopped due to Refinement!\n");
						} else
#endif

						/*int index = (int)floor(GetVoronoiCell( selected_action->originalCell)->position[1])
						 + (int)floor(GetVoronoiCell( selected_action->originalCell)->position[0])*voronoiDiagram->xN[1];
						 for( int dx=-1; dx<=1; dx++)
						 for( int dy=-1; dy<=1; dy++){
						 int nindex = index + dy + dx*voronoiDiagram->xN[1];
						 int x = voronoiDiagram->xN[0];
						 int y = voronoiDiagram->xN[1];

						 fprintf( stderr, "nindex: %i (%i, %i) -> refined? = %i\n", nindex, x, y, voronoiDiagram->voronoiCells[nindex]->refined);
						 if( (int)floor(GetVoronoiCell( selected_action->originalCell)->position[0])+dx >= 0 &&
						 (int)ceil( GetVoronoiCell( selected_action->originalCell)->position[0])+dx <= voronoiDiagram->xN[0] &&
						 (int)floor(GetVoronoiCell( selected_action->originalCell)->position[1])+dy >= 0 &&
						 (int)ceil( GetVoronoiCell( selected_action->originalCell)->position[1])+dy <= voronoiDiagram->xN[1] &&
						 voronoiDiagram->voronoiCells[nindex]->refined)
						 voronoiDiagram->coarsen(voronoiDiagram->voronoiCells[nindex], 1);
						 }*/

						/*for( int vc=0; vc<voronoiDiagram->countVoronoiCells; vc++)
						 if( voronoiDiagram->voronoiCells[vc]->countNeighborCells > 0)
						 voronoiDiagram->voronoiCells[vc]->checkFreeNeighbors();
						 */

						//printTriangulation( voronoiDiagram, "testBetween.eps", false);
						//printVoronoiDiagram( voronoiDiagram, "testVDBetween.eps", false);
						/*n=0;
						 for( ; n<GetVoronoiCell( selected_action->originalCell)->countNeighborCells; )
						 {
						 if( GetVoronoiCell( selected_action->originalCell)->neighborCells[n]->refined == true){
						 //voronoiDiagram->refine( GetVoronoiCell( selected_action->originalCell)->neighborCells[n], 10);
						 voronoiDiagram->coarsen( GetVoronoiCell( selected_action->originalCell)->neighborCells[n], 10);
						 n=0;
						 }else{
						 n++;
						 }
						 }*/
						/*int index = (int)floor(GetVoronoiCell( selected_action->originalCell)->position[1])
						 + (int)floor(GetVoronoiCell( selected_action->originalCell)->position[0])*voronoiDiagram->xN[1];
						 for( int dx=-1; dx<=1; dx++)
						 for( int dy=-1; dy<=1; dy++){
						 int nindex = index + dy + dx*voronoiDiagram->xN[1];
						 int x = voronoiDiagram->xN[0];
						 int y = voronoiDiagram->xN[1];

						 fprintf( stderr, "nindex: %i (%i, %i) -> refined? = %i\n", nindex, x, y, voronoiDiagram->voronoiCells[nindex]->refined);
						 if( (int)floor(GetVoronoiCell( selected_action->originalCell)->position[0])+dx >= 0 &&
						 (int)ceil( GetVoronoiCell( selected_action->originalCell)->position[0])+dx <= voronoiDiagram->xN[0] &&
						 (int)floor(GetVoronoiCell( selected_action->originalCell)->position[1])+dy >= 0 &&
						 (int)ceil( GetVoronoiCell( selected_action->originalCell)->position[1])+dy <= voronoiDiagram->xN[1] &&
						 voronoiDiagram->voronoiCells[nindex]->refined)
						 voronoiDiagram->coarsen(voronoiDiagram->voronoiCells[nindex], 1);
						 }*/

						//fprintf( stderr, "countVoronoiCells = %i\n", voronoiDiagram->countVoronoiCells);

						if (selected_action->originalCell->state == COMPARTMENT) {

							//if(selected_action->originalCell->growingTumorCellCount+selected_action->originalCell->dividingTumorCellCount != selected_action->originalCell->cellCount)
							//fprintf(stderr,"%i == g:%i + d:%i\n", selected_action->originalCell->cellCount, selected_action->originalCell->growingTumorCellCount,selected_action->originalCell->dividingTumorCellCount);
							//if(selected_action->originalCell->growingTumorCellCount == 0)
							//fprintf(stderr,"%i == g:%i + d:%i\n", selected_action->originalCell->cellCount, selected_action->originalCell->growingTumorCellCount,selected_action->originalCell->dividingTumorCellCount);
							//fprintf(stderr, "Compartment Growth\n");
							//exit(0);
							int whichCell =
									(int) myRandE(
											(double) selected_action->originalCell->growingTumorCellCount);
							int sumCells = 0;
							int randM = -1;
							for (int m = 0; m <= M_gro && randM < 0; m++) {
								if (selected_action->internalStateM[m] > 0) {
									sumCells +=
											selected_action->internalStateM[m];
									if (sumCells > whichCell) {
										randM = m;
									}
								}
							}

							//if(selected_action->internalState == M_gro){
							if (randM == M_gro) {
								/*fprintf(stderr,"%3.2lf DIVISION: %i =>", Time, selected_action->originalCell->growingTumorCellCount);
								 for(  i=0; i<=M_gro; i++)
								 if( selected_action->internalStateM[i]>0)
								 fprintf(stderr,"[m:%2i, c:%i]", i, selected_action->internalStateM[i]);
								 fprintf(stderr,"\n");*/

								VoronoiCell *destination =
										performGrowthAndDivision(voronoiDiagram,
												agentArray, actionList,
												selected_action, &myStatistics,
												Time, EndTime);
								if (destination != NULL) {
									//if (count_cells == MAX_CELLS)
									//	EndTime = Time;
									if (GetVoronoiCell(selected_action->originalCell)->isDomainBorder(
											voronoiDiagram)) {
										EndTime = Time;
									}

									//selected_action->internalState = 0;
									selected_action->internalStateM[randM]--;
									selected_action->internalStateM[0]++;
									GetAgent(destination)->actions[INDEX_GROWTH]->internalStateM[0]++;

									/*fprintf(stderr,"AFTER DIVISION: %i =>", selected_action->originalCell->growingTumorCellCount);
									 for(  i=0; i<=M_gro; i++)
									 if( selected_action->internalStateM[i]>0)
									 fprintf(stderr,"[m:%2i, c:%i]", i, selected_action->internalStateM[i]);
									 fprintf(stderr,"\n");*/
								}
							} else {
								//selected_action->internalState++;
								selected_action->internalStateM[randM]--;
								selected_action->internalStateM[randM + 1]++;
							}
						}

						else {

							if( selected_action->originalCell->location[0]->waste > Agent::WASTE_THRESHOLD_SLOWED_GROWTH
									|| selected_action->originalCell->getOxygen() < 0.07)
								//selected_action->originalCell->intoxication++;
								selected_action->originalCell->intoxication += 1./selected_action->rate;
								//selected_action->originalCell->intoxication += 1. / (M_div + M_gro + 2.);

							if(selected_action->originalCell->divide){ // <-- DIVISION / GROWTH


							if (selected_action->internalState == M_gro) {
								VoronoiCell *destination = performGrowth(
										voronoiDiagram, agentArray, actionList,
										selected_action, &myStatistics, Time,
										EndTime);
								if (M_div < 0 && destination != NULL/*&& selected_action->originalCell->countLocations >= 2*MIN_SUBCELLULAR_COMPONENTS*/) {
									performDivision(
											voronoiDiagram,
											agentArray,
											actionList,
											GetAgent(destination)->actions[INDEX_DIVISION],
											&myStatistics, Time, EndTime);
									//if (count_cells == MAX_CELLS)
									//	EndTime = Time;
									if (count_cells == voronoiDiagram->countVoronoiCells)
										Time = EndTime;
									if (AdaptEndTime
											&& GetAgent(destination)
													!= selected_action->originalCell)
										if (GetVoronoiCell(selected_action->originalCell)->isDomainBorder(
												voronoiDiagram)) {
											Time = EndTime;
											//EndTime = Time;
										}
								}
								selected_action->internalState = 0;
							}

							else {
								selected_action->internalState++;
							}

							}else{ /// <-- QUIESCENCE

								VoronoiCell *closestVC = voronoiDiagram->searchClosestUnoccupiedVoronoiCell(
										selected_action->originalCell->location[0],
										VoronoiCell::extendedNeighborDepth,
										VoronoiCell::symbolicExtendedNeighborhoodSize,
										VoronoiCell::symbolicExtendedNeighborhood);


								if (closestVC && selected_action->originalCell->countFreeNeighbors()){
									double min_dist =	closestVC->getDistanceSquareTo(
											selected_action->originalCell->location[0]);

									double prob_reentering_cell_cycle =
													PROB_REENTERING_CELL_CYCLE( sqrt((double)min_dist) * SPATIAL_UNIT);

									// INTOXICATION -> QUIESCENCE
									if( selected_action->originalCell->intoxication > Agent::WASTE_INTOXICATED_CELL_CYCLES/MaxCellDivisionRate)
										prob_reentering_cell_cycle=0;

									if (
													( myRand() <= prob_reentering_cell_cycle
													&& VoronoiCell::ECM_THRESHOLD_QUIESCENCE <= selected_action->originalCell->location[0]->ecm))
										selected_action->originalCell->divide = 1.;

								}

								//if( )
									//	selected_action->originalCell->divide = 1.;
							} /// <--
						}

#ifdef COARSEN
						int type;
						int x = (int)floor(GetVoronoiCell( selected_action->originalCell)->position[0]);
						int y = (int)floor(GetVoronoiCell( selected_action->originalCell)->position[1]);
#if DIMENSIONS == 2
						int dy = 1;
						int dx = voronoiDiagram->xN[1];
#elif DIMENSIONS == 3
						int z = (int)floor(GetVoronoiCell( selected_action->originalCell)->position[2]);
						int dz = 1;
						int dy = voronoiDiagram->xN[1];
						int dx = voronoiDiagram->xN[1]*voronoiDiagram->xN[2];
#endif

						for( int ix=MAX(0,x-1); ix<MIN(voronoiDiagram->xN[0],x+2); ix++)
						for( int iy=MAX(0,y-1); iy<MIN(voronoiDiagram->xN[1],y+2); iy++)
#if DIMENSIONS == 2
						if(voronoiDiagram->voronoiCells[ix*dx +iy*dy]->refined)
						{
							c = voronoiDiagram->voronoiCells[ix*dx +iy*dy]->refinement;
#elif DIMENSIONS == 3
							for( int iz=MAX(0,z-1); iz<MIN(voronoiDiagram->xN[2],z+2); iz++)
							if(voronoiDiagram->voronoiCells[ix*dx +iy*dy +iz*dz]->refined) { //TEST!!!!!
								c = voronoiDiagram->voronoiCells[ix*dx +iy*dy +iz*dz]->refinement;
#endif
								//for( int c=0; c<voronoiDiagram->countVoronoiCells; c++){
								//if( voronoiDiagram->voronoiCells[c]->refined)
								if( voronoiDiagram->voronoiCells[c]->countNeighborCells>0)
								if( voronoiDiagram->isHomogen( voronoiDiagram->voronoiCells[c], CountCellsPerVoronoiCell, type))
								if( type == NONACTIVE) {//NONACTIVE){
									//	if( voronoiDiagram->voronoiCells[c]->refined
									//			&& voronoiDiagram->voronoiCells[c]->countNeighborCells>0
									//			&& voronoiDiagram->isHomogen( voronoiDiagram->voronoiCells[c], CountCellsPerVoronoiCell, type)
									//			&& type == NONACTIVE){//NONACTIVE){
									//fprintf( stderr, "Found coarse & homogen cell\n");
									passedTime = clock();
									//printTriangulation( voronoiDiagram, "beforeCoarsenT.eps", true);
									//if( c == 1046000)
									//	checkDelaunayCondition( voronoiDiagram, 0, 0);
									VoronoiCell* comp = voronoiDiagram->coarsen( voronoiDiagram->voronoiCells[c], CountCellsPerVoronoiCell, actionList, agentArray);
									//fprintf(stderr, "...finished ( %.3lfsec)\n", (double)(clock() - passedTime)/(double)CLOCKS_PER_SEC);

									//fprintf(stderr, "still refined? %i\n", voronoiDiagram->voronoiCells[c]->refined);
									//fprintf(stderr, "neighbors? %i\n", voronoiDiagram->voronoiCells[c]->countNeighborCells);

									/*for( int a=0; a<agentArray->countActiveAgents; a++){
									 if( GetVoronoiCell(agentArray->agents[a])->countNeighborCells > 0 &&
									 GetVoronoiCell(agentArray->agents[a])->countFreeNeighborCells == 0 &&
									 agentArray->agents[a]->actionsInitialized &&
									 agentArray->agents[a]->actions[INDEX_GROWTH]->top!=NULL){
									 fprintf( stderr, "Agent %i on Cell %i shouldn't grow\n",
									 agentArray->agents[a]->index,
									 GetVoronoiCell(agentArray->agents[a])->index);
									 printTriangulation( voronoiDiagram, "afterCoarsenT.eps", true);

									 exit(0);
									 }
									 }*/

									//exit( 0);
									if( comp->agent == NULL)
									//exit(0);
									{
										fprintf( stderr, "After Coarsen no Agent attached Compartment!\n");
										exit(0);
										Agent *compAgent = agentArray->activateAgent();
										compAgent->attach(comp);
										compAgent->state = COMPARTMENT;//NONACTIVE;
										compAgent->cellCount = 1;
										compAgent->maxCellCount = 1;
										for( int n=0; n<comp->countNeighborCells; n++)
										comp->neighborCells[n]->countFreeNeighborCells--;
									}
									c=-1;
									//checkDelaunayCondition( voronoiDiagram, 0, 0);
								}
							}
									//printTriangulation( voronoiDiagram, "afterCoarsenT.eps", false);
									//printVoronoiDiagram( voronoiDiagram, "afterCoarsenVD.eps", false);

							/*for( int c=0; c<voronoiDiagram->countVoronoiCells; c++){
							 if( voronoiDiagram->voronoiCells[c]->countNeighborCells>0 &&
							 voronoiDiagram->voronoiCells[c]->countFreeNeighborCells==0 &&
							 voronoiDiagram->voronoiCells[c]->getState()!=FREE &&
							 GetAgent(voronoiDiagram->voronoiCells[c])!=NULL &&
							 GetAgent(voronoiDiagram->voronoiCells[c])->actionsInitialized &&
							 GetAgent(voronoiDiagram->voronoiCells[c])->actions[INDEX_GROWTH]->top!=NULL){
							 fprintf( stderr, "Agent %i on Cell %i shouldn't grow\n",
							 GetAgent(voronoiDiagram->voronoiCells[c])->index,
							 voronoiDiagram->voronoiCells[c]->index);
							 exit(0);
							 }
							 }*/
#endif

					}

						break;
#endif

#if USE_DIVISION
					case DIVISION:
						//fprintf(stderr, "DIVISION?\n");
						if (selected_action->internalState == M_div) {
							selected_action->internalState = 0;
							performDivision(voronoiDiagram, agentArray,
									actionList, selected_action, &myStatistics,
									Time, EndTime);
						} else {
							selected_action->internalState++;
						}

						if (AdaptEndTime)
							if (!(
									GetPosition( selected_action->originalCell)[0]
									        > voronoiDiagram->xMin[0]
											        + voronoiDiagram->boundaryThickness
									&& GetPosition( selected_action->originalCell)[0]
											< voronoiDiagram->xMax[0]
													- voronoiDiagram->boundaryThickness
									&& GetPosition( selected_action->originalCell)[1]
											> voronoiDiagram->xMin[1]
													+ voronoiDiagram->boundaryThickness
									&& GetPosition( selected_action->originalCell)[1]
											< voronoiDiagram->xMax[1]
#if DIMENSIONS == 3
													- voronoiDiagram->boundaryThickness
									&& GetPosition( selected_action->originalCell)[2]
											> voronoiDiagram->xMin[2]
													+ voronoiDiagram->boundaryThickness
									&& GetPosition( selected_action->originalCell)[2]
											< voronoiDiagram->xMax[2]
													- voronoiDiagram->boundaryThickness
#endif
							)
						) {
								Time = EndTime;
								//EndTime = Time;
							}
						break;
#endif //USE_DIVISION
#if USE_NECROSIS

					case NECROSIS:
						//fprintf(stderr,"NECROSIS: %i\n", GetVoronoiCell( selected_action->originalCell)->index);
						//fprintf(stderr,"NECROSIS: %i(%i-1/%i)\n", selected_action->originalCell->nr, selected_action->originalCell->cellCount, selected_action->originalCell->maxCellCount);
						//selected_action->originalCell->state = FREE;
						//selected_action->originalCell->cellCount--;
						//GrosseFactors_Concentration[ selected_action->originalCell->mx][ selected_action->originalCell->my][ selected_action->originalCell->mz] += GROWTHFACTOR_THRESHOLD;
						//GiveMeAll( (Time - Last_Time), GROWTHFACTORS_MIGRATION);

						//if(agentArray->countActiveAgents > 100 /*&& (selected_action->originalCell->divide==0 || selected_action->originalCell->state == NONACTIVE)*/)
						if (selected_action->internalState == M_nec) {
							//for(int i=0; i<selected_action->originalCell->countLocations; i++)
							//	selected_action->originalCell->location[i]->growthfactors += 1000.;//GROWTHFACTOR_THRESHOLD;
							//update_deleted_cell( actionList, selected_action->originalCell, voronoiDiagram);

							/*
							 if( selected_action->originalCell->state != COMPARTMENT){
							 selected_action->originalCell->state = NECROTIC;
							 selected_action->originalCell->actualize( actionList);
							 }

							 if( selected_action->originalCell->growingTumorCellCount > 0)
							 selected_action->originalCell->growingTumorCellCount--;
							 else
							 selected_action->originalCell->dividingTumorCellCount--;
							 selected_action->originalCell->necroticCellCount++;
							 count_necrotic_cells ++;
							 count_necrotic_volume += selected_action->originalCell->countLocations;
							 count_cells --;
							 count_cell_volume -= selected_action->originalCell->countLocations;
							 */

							if (selected_action->originalCell->state
									== COMPARTMENT) {
								//fprintf(stderr, "NECROSIS\n");//exit( 0);

								// CHOSE DYING CELL
								int whichCell =
										(int) myRandE(
												(double) selected_action->originalCell->growingTumorCellCount);
								int sumCells = 0;
								int randM = -1;
								for (int m = 0; m <= M_gro && randM < 0; m++) {
									if (selected_action->internalStateM[m] > 0) {
										sumCells +=
												selected_action->internalStateM[m];
										if (sumCells > whichCell) {
											randM = m;
										}
									}
								}

								// CHANGE STATE
								selected_action->originalCell->growingTumorCellCount--;
#if USE_GROWTH
								selected_action->originalCell->actions[INDEX_GROWTH]->internalStateM[randM]--;
#endif
								selected_action->originalCell->necroticCellCount++;

								// ACTUALIZE ACTIONS OF CELL: NECROSIS & GROWTH
#if USE_GROWTH
								actionList->actualizeRate(
										selected_action->originalCell->actions[INDEX_GROWTH],
										selected_action->originalCell->actions[INDEX_GROWTH]->getActualRate());
#endif
								actionList->actualizeRate(
										selected_action->originalCell->actions[INDEX_NECROSIS],
										selected_action->originalCell->actions[INDEX_NECROSIS]->getActualRate());
								selected_action->originalCell->actualize(
										actionList);

								// STATISTICS
								count_necrotic_cells++;
								count_necrotic_volume +=
										selected_action->originalCell->countLocations;
								count_cells--;
								count_cell_volume -=
										selected_action->originalCell->countLocations;
							}

							else {
								//fprintf(stderr, "NECROSIS\n");

								// CHANGE STATE
								selected_action->originalCell->state = NECROTIC;

								// ACTUALIZE ACTIONS OF CELL
								if (selected_action->originalCell->growingTumorCellCount
										> 0)
									selected_action->originalCell->growingTumorCellCount--;
								else
									selected_action->originalCell->dividingTumorCellCount--;
								selected_action->originalCell->necroticCellCount++;
								selected_action->originalCell->actualize(
										actionList);

								// STATISTICS
								count_necrotic_cells++;
								count_necrotic_volume +=
										selected_action->originalCell->countLocations;
								count_cells--;
								count_cell_volume -=
										selected_action->originalCell->countLocations;
							}

							//fprintf( stderr, "INFO: new cell state: %s, cell count: %i\n", cellTypeToString(selected_action->originalCell->state), selected_action->originalCell->cellCount);
						} else {
							selected_action->internalState++;
						}

						break;
#endif //USE_NECROSIS
#if USE_LYSIS

					case LYSIS:

						//fprintf(stderr, "LYSIS\n");
						//exit(0);
						//fprintf(stderr,"LYSIS: %i %s (total:%i/%i, tumor:%i, necrotic:%i-1)\n", GetIndex(selected_action->originalCell), cellTypeToString(selected_action->originalCell->state), selected_action->originalCell->cellCount, selected_action->originalCell->maxCellCount, selected_action->originalCell->tumorCellCount, selected_action->originalCell->necroticCellCount);
						count_lysed_volume +=
								selected_action->originalCell->countLocations;
						count_necrotic_volume -=
								selected_action->originalCell->countLocations;
						count_inner_empty_volume +=
								selected_action->originalCell->countLocations;
#if MULTISCALE
						if( selected_action->originalCell->cellCount-- == selected_action->originalCell->maxCellCount) {
							update_surrounding_deleted_cell(actionList, selected_action->originalCell, voronoiDiagram);
							if( selected_action->originalCell->state == NONACTIVE)
							addDivisionAction( actionList, selected_action->originalCell, voronoiDiagram);
							selected_action->originalCell->state = FREE;
						}
						if( selected_action->originalCell->cellCount == 0) {
							update_deleted_cell( actionList, selected_action->originalCell, voronoiDiagram);
						}
#else

						//fprintf( stderr, "INFO: do single based lysis\n");
						if (selected_action->originalCell->state != COMPARTMENT) {

							//update_surrounding_deleted_cell(actionList, selected_action->originalCell, voronoiDiagram);
							update_surrounding_removed_cell(actionList,
									selected_action->originalCell,
									voronoiDiagram);
							//deleteAction( actionList, selected_action);
							selected_action->originalCell->state = FREE;
							selected_action->originalCell->actualize(
									actionList);
						}
						//selected_action->originalCell->detach();
						//agentArray->deactivateAgent( selected_action->originalCell);
#endif

						selected_action->originalCell->necroticCellCount--;
						selected_action->originalCell->cellCount--;
						count_necrotic_cells--;
						count_lysed_cells++;

						//flag_actualizeProcessRates = TRUE;

						break;

#endif //USE_LYSIS
					case GAP:
						//fprintf(stderr, "\n GAP \n");
						//Time = Last_Time + 1. / selected_action->rate;
						break;

					default:
						fprintf(
								stderr,
								"\nError: Selected action (Type: %i, %s) is unknown!\n",
								selected_action->type,
								actionTypeToString(selected_action->type));
						exit(0);
					}

					timeExecution +=
							(double) (clock() - passedTime)
									/ (double) CLOCKS_PER_SEC;

				} else {
					if (selected_action->rate == 0.) {
						//	fprintf( stderr, "selected_action->rate = %lf\n", selected_action->rate);
						Time = Last_Time + TIME_STEP * 1.;
					}
					if (Time > EndTime) {
						//	fprintf( stderr, "Time(%lf)>EndTime(%lf)\n", Time, EndTime);
						// correct time
						Time = EndTime;
					}
					//fprintf( stderr, "selected_action->type = %i\n", selected_action->type);
					//fprintf( stderr, "selected_action->rate = %lf\n", selected_action->rate);
				}
			} else {
				//fprintf( stderr, "actionList->size = %i\n", actionList->size);
				// Empty Probability List
				//Last_Time = Time;
				if (actionList->rateSum == 0.) {
					//fprintf( stderr, "actionList->rateSum = %lf\n", actionList->rateSum);
					Time += TIME_STEP * 10.;
				}
				if (actionList->size == 0) {
					//fprintf( stderr, "actionList->size = %i\n", actionList->size);
					Time = EndTime;
				}
			}

			if (Time > EndTime) {
				//	fprintf( stderr, "Time(%lf)>EndTime(%lf)\n", Time, EndTime);
				// correct time
				Time = EndTime;
			}
			//fprintf( stderr, "Finished Selected Action\n");


			//if( Last_Time + 24 < Time){
			//	fprintf( stderr, "Last_Time = %lf, Time = %lf\n", Last_Time, Time);
			//}
			//printf("substrate.Update_All();\n");
			//if( Case > 1) //update in any case for several reason to ask dial: 5413
			/*for( int a=0; a<agentArray->countActiveAgents; a++){
			 if(agentArray->agents[a]->actionsInitialized)
			 if( agentArray->agents[a]->actions[INDEX_MIGRATION]->active() && agentArray->agents[a]->countFreeNeighbors()==0){
			 fprintf( stderr, "MIGRATION of agent %i shouldn't be active!\n", agentArray->agents[a]->index);

			 VoronoiCell* centralCell = agentArray->agents[a]->location[0];
			 fprintf( stderr, "cells: %i/%i -> g:%i, d:%i, n:%i\n", GetAgent(centralCell)->cellCount, GetAgent(centralCell)->maxCellCount, GetAgent(centralCell)->growingTumorCellCount, GetAgent(centralCell)->dividingTumorCellCount, GetAgent(centralCell)->necroticCellCount );
			 for( int v=0; v<centralCell->countNeighborCells; v++)
			 fprintf( stderr, "%i. neighbor(%i): %i/%i -> g:%i, d:%i, n:%i\n", v+1, GetAgent(centralCell->neighborCells[v])->index, GetAgent(centralCell->neighborCells[v])->cellCount, GetAgent(centralCell->neighborCells[v])->maxCellCount, GetAgent(centralCell->neighborCells[v])->growingTumorCellCount, GetAgent(centralCell->neighborCells[v])->dividingTumorCellCount, GetAgent(centralCell->neighborCells[v])->necroticCellCount );

			 exit(0);
			 }
			 }*/

			//timeDifference = pharmaco->evolveSystem( voronoiDiagram, Last_Time - timeDifference, Time, customTimestep);
			if (Time != EndTime && Case != 1 /*&& count_inoculated_cells + (int)count_divided_cells - count_necrotic_cells - count_lysed_cells!=0*/) {

				//printf("STEPS: %i\n", (int)((Time-Last_Time<0.03 ? Time-Last_Time : 0.03)/TIME_STEP));
				//printf("STEPS: %i\n", (int)((Time-Last_Time<0.03 ? Time-Last_Time : 0.03)/TIME_STEP));
				//if( 0==(int)((float)(Time-Last_Time)/TIME_STEP))
				//	fprintf( stderr, "Difference of Time(%lf) and Last_Time(%lf) too small: %e == %lf?\n", Time, Last_Time, Time-Last_Time, TIME_STEP);

#if REFILL_MEDIUM
				double refillIntervall = 10.; //h
				double refilledPortion = 0.90;
				double totalVolume = 200e12;//micrometer^3
				double minSpheroidNumber = 100.;
				double maxSpheroidNumber = 10000.;
				double duration = 24.*24.;// days*hours

				// linear
				//double countSpheroids = maxSpheroidNumber+(minSpheroidNumber-maxSpheroidNumber)/duration*refillIntervall*(double)((int)(Time/refillIntervall));

				// exponential
				double countSpheroids = maxSpheroidNumber*exp(-((int)(Time/refillIntervall)*refillIntervall)*log(maxSpheroidNumber/minSpheroidNumber)/duration);

				//double countSpheroids =   maxSpheroidNumber*pow( 40.5,(-(double)((int)(Time/refillIntervall))*(14./24.)/10.2))+50.;
				//double countSpheroids = 10000.*pow(-((int)(Time/14.))*(14./24.)/3.2, 2.);
				if( countSpheroids<minSpheroidNumber) countSpheroids=minSpheroidNumber;

				//double totalDiameter = (double)( CUBE_LENGTH_TO_DIAMETER( pow(totalVolume / countSpheroids, 1/3.))/SPATIAL_UNIT);
				//double occupiedDiameter = count_cell_volume + count_necrotic_volume + count_inner_empty_volume;
				//fprintf( stderr, "%lf -> %lf -> occupied: %lf, free: %lf\n", countSpheroids, totalDiameter, occupiedDiameter, totalDiameter - occupiedDiameter);
				//printf( "%lf %lf %lf %lf %lf\n", Time, global_glucose_concentration, centralCell->glucose, global_oxygen_concentration, centralCell->oxygen);
				for( int j=0; j<(int)((float)(Time-Last_Time)/TIME_STEP)+1; j++) {
					double myTimeStep = (j!=(int)((float)(Time-Last_Time)/TIME_STEP) ? TIME_STEP : Time-Last_Time-TIME_STEP*j);
				//for( int j=0; j<(int)((Time-Last_Time<0.03 ? Time-Last_Time : 0.03)/TIME_STEP); j++){
				//printf("substrate.Update_All();\n");
				//printf("dURATION:%lf\n",Time-Last_Time);
					/*double timeStep = Time-Last_Time;
					 for( i=0; i<10; i++)
					 if( substrate.Update_All( Case ,timeStep/10.))
					 flag_actualizeProcessRates = TRUE;*/
					/*double global_oxygen = 0.;
					 double global_glucose = 0.;
					 double mol_oxygen = 0.;
					 double mol_glucose = 0.;
					 int countBorderCells = 0;
					 // find number of molecules outside the tumor before
					 for( i=0; i<voronoiDiagram->countVoronoiCells; i++){
					 if( GetAgent(voronoiDiagram->voronoiCells[i])==NULL){
					 global_oxygen  += voronoiDiagram->voronoiCells[i]->oxygen ;
					 global_glucose += voronoiDiagram->voronoiCells[i]->glucose;
					 countBorderCells++;
					 }
					 }
					 mol_oxygen = global_oxygen*AGENT_VOLUME*CUBICMICROMETER_TO_MILLILITER*1000.;
					 mol_glucose = global_glucose*AGENT_VOLUME*CUBICMICROMETER_TO_MILLILITER*1000.;
					 global_oxygen  /= (double)countBorderCells;
					 global_glucose /= (double)countBorderCells;
					 */

					//if( substrate.Update_All( Case ,(Time-Last_Time<0.03 ? Time-Last_Time : 0.03)))
					//if( substrate.Update_All( Case ,myTimeStep)){
					//substrate.Update_All( Case ,myTimeStep);
					double glucoseCon = 0.;
					double glucoseConDiff = 0.;
					double oxygenCon = 0.;
					double oxygenConDiff = 0.;
					int borderSize = 0;
					/*VoronoiCell **border = getOuterBorder( agentArray, borderSize);
					 for( i=0; i<borderSize; i++){
					 //	fprintf(stderr, "%i ", border[i]->index);
					 glucoseCon += border[i]->glucose;
					 }*/
					/*for( i=0; i<voronoiDiagram->countVoronoiCells; i++){
					 //	fprintf(stderr, "%i ", border[i]->index);
					 if( GetAgent( voronoiDiagram->voronoiCells[i])==NULL){
					 voronoiDiagram->voronoiCells[i]->glucose = global_glucose_concentration;
					 voronoiDiagram->voronoiCells[i]->oxygen  = global_oxygen_concentration;
					 glucoseCon += voronoiDiagram->voronoiCells[i]->glucose;
					 oxygenCon  += voronoiDiagram->voronoiCells[i]->oxygen;
					 borderSize++;
					 }
					 }*/
					//fprintf(stderr, "\n");
					//fprintf(stderr, "BEFORE FREE\n");
					//fprintf(stderr, "borderSize = %i\n", borderSize);
					//fprintf(stderr, "[G] = %lf ==> ", glucoseCon/(double)borderSize);
					//glucoseConDiff = glucoseCon/(double)borderSize;
					//oxygenConDiff  = oxygenCon/(double)borderSize;
					//fprintf(stderr, "AFTER FREE\n");
					substrate.Update_All(Case, myTimeStep, actionList, agentArray, voronoiDiagram);
					flag_actualizeProcessRates = TRUE;
					//fprintf(stderr, "AFTER UPDATE\n");
					//printf("flag_actualizeProcessRates: %i \n", flag_actualizeProcessRates);
					//}
					//borderSize = 0;
					//border = getOuterBorder( agentArray, borderSize);
					/*for( i=0; i<borderSize; i++){
					 //	fprintf(stderr, "%i ", border[i]->index);
					 glucoseCon += border[i]->glucose;
					 }*/
					glucoseCon = 0.;
					oxygenCon = 0.;
					for( i=0; i<voronoiDiagram->countVoronoiCells; i++) {
						//	fprintf(stderr, "%i ", border[i]->index);
						if( GetAgent( voronoiDiagram->voronoiCells[i])==NULL) {
							glucoseCon += voronoiDiagram->voronoiCells[i]->glucose;
							oxygenCon += voronoiDiagram->voronoiCells[i]->oxygen;
							borderSize++;
						}
					}
					//fprintf(stderr, "[G] = %lf\n", glucoseCon/(double)borderSize);
					// delta [G] per site in the simulated section
					glucoseConDiff = glucoseCon/(double)borderSize - global_glucose_concentration;
					oxygenConDiff = oxygenCon/(double)borderSize - global_oxygen_concentration;
					//fprintf(stderr, "d[G] = %lf\n", glucoseDiffMol);

					// delta [G] per site in the volume portion of the spheroid
					glucoseConDiff *= (double)borderSize/(double)(totalVolume/AGENT_VOLUME + borderSize - voronoiDiagram->countVoronoiCells);
					oxygenConDiff *= (double)borderSize/(double)(totalVolume/AGENT_VOLUME + borderSize - voronoiDiagram->countVoronoiCells);
					//fprintf(stderr, "d[G] = %e\n", glucoseConDiff);

					global_glucose_concentration += glucoseConDiff;
					global_oxygen_concentration += oxygenConDiff;
					//fprintf(stderr, "d[G] = %lf\n", global_glucose_concentration);

					//free( border);

					// find number of molecules outside the tumor after
					/*double global_oxygen_diff  = 0.;
					 double global_glucose_diff = 0.;
					 double mol_oxygen_diff  = 0.;
					 double mol_glucose_diff = 0.;
					 for( i=0; i<voronoiDiagram->countVoronoiCells; i++){
					 if( GetAgent(voronoiDiagram->voronoiCells[i])==NULL){
					 global_oxygen_diff  += voronoiDiagram->voronoiCells[i]->oxygen;  // - global_oxygen;
					 global_glucose_diff += voronoiDiagram->voronoiCells[i]->glucose; // - global_glucose;
					 }
					 }
					 mol_oxygen_diff = global_oxygen_diff*AGENT_VOLUME*CUBICMICROMETER_TO_MILLILITER*1000.;
					 mol_glucose_diff = global_glucose_diff*AGENT_VOLUME*CUBICMICROMETER_TO_MILLILITER*1000.;
					 global_oxygen_diff  /= (double)countBorderCells;
					 global_glucose_diff /= (double)countBorderCells;

					 //printf( "[G]: %lf, consumed glucose: %e mmol\n", global_glucose, (global_glucose_diff- global_glucose)*(double)countBorderCells*AGENT_VOLUME*CUBICMICROMETER_TO_MILLILITER*1000.);
					 //printf( "[G]: %lf, consumed glucose: %e mmol\n", global_glucose, (global_glucose_diff- global_glucose)*(double)countBorderCells*AGENT_VOLUME*CUBICMICROMETER_TO_MILLILITER*1000.);
					 //printf( "avg. [G]: %e, avg. nG: %e, consumed glucose: %e mmol\n", global_glucose, mol_glucose, (global_glucose_diff- global_glucose)*(double)countBorderCells*AGENT_VOLUME*CUBICMICROMETER_TO_MILLILITER*1000.);

					 // local difference
					 mol_oxygen_diff = (mol_oxygen_diff - mol_oxygen) / 2.;
					 mol_glucose_diff = (mol_glucose_diff - mol_glucose) / 2.;
					 //printf( "local consumed nG: %e mmol\n", mol_glucose_diff);

					 // global difference
					 mol_oxygen_diff *= PI *occupiedDiameter*occupiedDiameter;
					 mol_glucose_diff*= PI *occupiedDiameter*occupiedDiameter;
					 //printf( "global consumed nG: %e mmol\n", mol_glucose_diff);

					 // resulting global concentration
					 //printf( "new global [G]: %e mmol/L\n", global_glucose + mol_glucose_diff/(PI/6.*(pow(countPoints[0],3.)-pow(occupiedDiameter,3.)))/(AGENT_VOLUME*CUBICMICROMETER_TO_MILLILITER*1000.));


					 //double countSpheroids = 10000.*exp(-((int)(Time/14.))*(14./24.)/4.5);
					 //double countSpheroids = 10000.*pow(-((int)(Time/14.))*(14./24.)/3.2, 2.);
					 //if( countSpheroids<50.)countSpheroids=50.;
					 //printf( "%lf -> %lf\n", countSpheroids, (double)(pow(200e12 / countSpheroids, 1/3.)/0.82/SPATIAL_UNIT));
					 if(global_oxygen_diff!=0.)
					 global_oxygen_diff  = (global_oxygen_diff - global_oxygen)*(double)countBorderCells  / (double)(pow(500e12 / countSpheroids, 1/3.)/0.82/SPATIAL_UNIT + (double)countBorderCells - countPoints[0]);
					 if(global_glucose_diff!=0.)
					 global_glucose_diff = (global_glucose_diff- global_glucose)*(double)countBorderCells / (double)(pow(500e12 / countSpheroids, 1/3.)/0.82/SPATIAL_UNIT + (double)countBorderCells - countPoints[0]);	*/

					// stirre medium
					if( (int)((Last_Time+TIME_STEP*j)/refillIntervall) != (int)((Last_Time+TIME_STEP*(j+1))/refillIntervall)) {
						global_glucose_concentration = InitialGlucoseConcentration*refilledPortion + (1.-refilledPortion)*global_glucose_concentration;
						global_oxygen_concentration = InitialOxygenConcentration *refilledPortion + (1.-refilledPortion)*global_oxygen_concentration;
					}//else{
					 //global_glucose_concentration =
					 //}
					for( i=0; i<voronoiDiagram->countVoronoiCells; i++) {
						if( GetAgent(voronoiDiagram->voronoiCells[i])==NULL) {
							voronoiDiagram->voronoiCells[i]->glucose = global_glucose_concentration;
							voronoiDiagram->voronoiCells[i]->oxygen = global_oxygen_concentration;
							/*else{
							 voronoiDiagram->voronoiCells[i]->oxygen  = global_oxygen + mol_oxygen_diff/(PI/6.*(pow(totalDiameter,3.)-pow(occupiedDiameter,3.)))/(AGENT_VOLUME*CUBICMICROMETER_TO_MILLILITER*1000.);
							 voronoiDiagram->voronoiCells[i]->glucose = global_glucose + mol_glucose_diff/(PI/6.*(pow(totalDiameter,3.)-pow(occupiedDiameter,3.)))/(AGENT_VOLUME*CUBICMICROMETER_TO_MILLILITER*1000.);
							 //voronoiDiagram->voronoiCells[i]->oxygen  = global_oxygen  + global_oxygen_diff;
							 //voronoiDiagram->voronoiCells[i]->glucose = global_glucose + global_glucose_diff;
							 }	*/
						}
					}
					//printf("%lf %lf %lf %lf %lf %lf\n", Last_Time+TIME_STEP*j, global_glucose + global_glucose_diff, global_oxygen  + global_oxygen_diff, Last_Time, Time, EndTime);
					//printf("%lf %lf %lf %lf %lf %lf\n", Last_Time+TIME_STEP*j, voronoiDiagram->voronoiCells[0]->oxygen, voronoiDiagram->voronoiCells[0]->glucose, centralCell->oxygen, centralCell->glucose, GiveMeTheATP(centralCell));//, countSpheroids);

					//printf("before: oxy: %e, glu: %e ", global_oxygen, global_glucose);
					//printf(", change: oxy: %e, glu: %e", global_oxygen_diff, global_glucose_diff);
					//printf(", after: oxy: %e, glu: %e \n", global_oxygen/(double)countBorderCells  + global_oxygen_diff, global_glucose/(double)countBorderCells + global_glucose_diff);

				}
#else
				//if( substrate.Update_All( Case ,Time-Last_Time)){
				/*if( substrate.Update_All(Case, Time-Last_Time, actionList)){
				 flag_actualizeProcessRates = TRUE;
				 }*/
				//timeDifference = UpdateSystem( voronoiDiagram, customTimestep, Time-Last_Time + timeDifference);
				//timeDifference += Time-Last_Time;
				timeDifference += Time - Last_Time;

				// UPDATE MORPHOGEN
				if (timeDifference > customTimestep)
				if (VoronoiCell::USE_MORPHOGEN) {

					//morphogenDynamics->setBoundaryCondition( DiffusionReactionEquation::NEUMANN);
					morphogenDynamics->setBoundaryCondition( DiffusionReactionEquation::DIRICHLET);

					for (int i = 0; i < voronoiDiagram->countVoronoiCells;
							i++) {

						// GENERATION
						if (GetAgent(voronoiDiagram->voronoiCells[i])
								&& voronoiDiagram->voronoiCells[i]->getState() != FREE
								//&& voronoiDiagram->voronoiCells[i]->getState() == NECROTIC
										)
							morphogenProduction[i] = -300;

						// DECAY
						//morphogenProduction[i] += 10000*voronoiDiagram->voronoiCells[i]->morphogen;
						morphogenDecay[i] = 100;

						/*if( voronoiDiagram->voronoiCells[i]->getState() == FREE)
						 morphogenDiffusion[i] = Lactate_Diffusion * 10.;
						 else
						 morphogenDiffusion[i] = Lactate_Diffusion;*/
					}

					// update morphogen
					morphogenDynamics->update(timeDifference);

					// hard copy
					for (int i = 0; i < voronoiDiagram->countVoronoiCells; i++)
						voronoiDiagram->voronoiCells[i]->morphogen =
								morphogen[i];
				}

				// UPDATE GLUCOSE & OXYGEN & LACTATE

				if (timeDifference > customTimestep) {

					// UPDATE GLUCOSE & OXYGEN

					// CG for linear system
					//timeDifference = UpdateSystemImplicitSparse( voronoiDiagram, timeDifference, timeDifference);
					//timeDifference = UpdateSystemNewtonSparse( voronoiDiagram, timeDifference, timeDifference);
					//for(int it=0; it<1; it++)
					passedTime = clock();
					if (Case == 5) {
						// Diffusion Reaction of Growth Factors
						//substrate.Produce_Growthfactors();
						UpdateGrowthFactorsNonLinearCGSparse(voronoiDiagram,
								Time - timeDifference, Time, customTimestep);

						// Update Vessel Network
						if (substrate.vasculature.Update_Network(timeDifference,
								actionList, voronoiDiagram)) {
							// Update Pressure in Vessels
							substrate.vasculature.GiveMeTheBlood();
							// Reset Initial Oxygen Concentration in Vessels
							substrate.Refill_Oxygen(Case);
							// Reset Initial Glucose Concentration in Vessels
							substrate.Refill_Glucose(Case);
						}
					}

					// Diffusion Reaction of Glucose & Oxygen
					timeDifference =
							UpdateSystemNonLinearCGSparse(voronoiDiagram,
									Time - timeDifference, Time, customTimestep);

					//pharmaco->evolveSystem( voronoiDiagram, timeDifference, customTimestep);
					//timeDifference = UpdateSystemNewtonCGSparse( voronoiDiagram, Time - timeDifference, Time, customTimestep);

					/*fprintf(stderr,
							"...finished ( %lisec) -> time rest = %lf\n",
							(clock() - passedTime) / CLOCKS_PER_SEC,
							timeDifference);*/
					benchmarkTime += (clock() - passedTime);
					//	exit(0);
					// Newton for non-linear system
					//timeDifference = UpdateSystemNewtonSparse( voronoiDiagram, timeDifference/10, timeDifference);

					//timeDifference = UpdateSystemImplicitSparse( voronoiDiagram, customTimestep, timeDifference);
					flag_actualizeProcessRates = TRUE;
					//timeDifference = 0.;

					// UPDATE LACTATE
					if (VoronoiCell::USE_LACTATE) {
						for (int i = 0; i < voronoiDiagram->countVoronoiCells;
								i++) {
							lactateProduction[i] =
									-GetLactateProductionRate(
											voronoiDiagram->voronoiCells[i],
											voronoiDiagram->voronoiCells[i]->glucose,
											voronoiDiagram->voronoiCells[i]->oxygen);
							if (voronoiDiagram->voronoiCells[i]->isFree())
								lactateDiffusion[i] = 30 * Lactate_Diffusion;
							else
								lactateDiffusion[i] = Lactate_Diffusion;
						}

						// update lactate
						lactateDynamics->update();

						// hard copy
						for (int i = 0; i < voronoiDiagram->countVoronoiCells;
								i++)
							voronoiDiagram->voronoiCells[i]->lactate =
									lactate[i];

					}

					// UPDATE WASTE
					if (VoronoiCell::USE_WASTE) {
						for (int i = 0; i < voronoiDiagram->countVoronoiCells;
								i++) {

							// waste production
							wasteProduction[i] = (voronoiDiagram->voronoiCells[i]->getState() == NECROTIC ? -10 : 0);
							// waste uptake
							if (!voronoiDiagram->voronoiCells[i]->isFree() && voronoiDiagram->voronoiCells[i]->getState() != NECROTIC){
								wasteUptake[i] = Agent::WASTE_UPTAKE;
								GetAgent( voronoiDiagram->voronoiCells[i])->waste += Agent::WASTE_UPTAKE * waste[i] * customTimestep;
							}else{
								wasteUptake[i] = 0;
							}


							if (voronoiDiagram->voronoiCells[i]->isFree())
								wasteDiffusion[i] = 30 * custom_Waste_Diffusion;
							else{
								wasteDiffusion[i] = custom_Waste_Diffusion;
							}
						}

						// update waste
						wasteDynamics->update();

						// hard copy
						for (int i = 0; i < voronoiDiagram->countVoronoiCells;
								i++)
							voronoiDiagram->voronoiCells[i]->waste = waste[i];

					}

				}

#endif
			}
			/*if( (int)(Last_Time/24.) < (int)(Time/24.))
			 fprintf( stderr, "\r%.3f%% (%d/%d) (%.1lfdays, volume: %d, tumor cells: %d, necrotic: %d, central glu: %e oxy: %e)\t\t",
			 (Time * 100)/EndTime,
			 k+1,
			 Averages,
			 Time/24.,
			 count_cell_volume,
			 count_cells,
			 count_necrotic_cells,
			 centralCell->glucose,
			 centralCell->oxygen//,
			 //actionList->size//,
			 //agentArray->countActiveAgents
			 );*/

			//test_prob_list(actionList, voronoiGrid);
			//test_voronoi_grid( voronoiGrid);
			/*printActionList( actionList);
			 int j;
			 for( j=0; j<voronoiDiagram->countVoronoiCells; j++){
			 int count_free_neighbors = 0;
			 for(i = 0; i < voronoiDiagram->voronoiCells[j]->countNeighborCells; i++)
			 if( GetAgent(voronoiDiagram->voronoiCells[j]->neighborCells[i])->state == FREE)
			 count_free_neighbors ++;
			 if( voronoiDiagram->voronoiCells[j]->countFreeNeighborCells != count_free_neighbors){
			 fprintf(stderr, "END OF THE LOOP:  cell %i (%s): free neighbors = %d ? = %d = count_free_neighbors\n",
			 voronoiDiagram->voronoiCells[j]->index,
			 cellTypeToString( voronoiDiagram->voronoiCells[j]->getState()),
			 voronoiDiagram->voronoiCells[j]->countFreeNeighborCells,
			 count_free_neighbors);
			 exit(0);
			 }
			 if( count_free_neighbors==0 && GetAgent(voronoiDiagram->voronoiCells[j])->state == ACTIVE){
			 fprintf(stderr, "END OF THE LOOP:  cell %i (%s) has %i free neighbors!!!\n",
			 voronoiDiagram->voronoiCells[j]->index,
			 cellTypeToString( GetAgent(voronoiDiagram->voronoiCells[j])->state),
			 voronoiDiagram->voronoiCells[j]->countFreeNeighborCells);
			 for( i=0; i<centralCell->countNeighborCells; i++){
			 fprintf( stderr, "%i. neighbor: %i (%s) -> free neighbors %i/%i\n",
			 i+1, voronoiDiagram->voronoiCells[j]->neighborCells[i]->index,
			 cellTypeToString( voronoiDiagram->voronoiCells[j]->neighborCells[i]->getState()),
			 voronoiDiagram->voronoiCells[j]->neighborCells[i]->countFreeNeighborCells,
			 voronoiDiagram->voronoiCells[j]->neighborCells[i]->countNeighborCells);
			 }

			 exit(0);
			 }
			 }*/

			// check if all cells are represented correctly in the list
			/*int valid = 1;
			 for( j=0; j<voronoiDiagram->countVoronoiCells; j++){
			 int countDivisionAction = 0,
			 countNecrosisAction = 0,
			 countLysisAction = 0,
			 countVesselGrowthAction = 0;
			 Action *tempProbListElem = actionList->head;
			 if(actionList->length!=0)
			 do{
			 if( tempProbListElem->originalCell!=NULL && GetVoronoiCell( tempProbListElem->originalCell)==voronoiDiagram->voronoiCells[j])
			 switch( tempProbListElem->type){
			 case DIVISION:
			 countDivisionAction++;
			 break;

			 case LYSIS:
			 countLysisAction++;
			 break;

			 case NECROSIS:
			 countNecrosisAction++;
			 break;

			 case VESSEL_GROWTH:
			 countVesselGrowthAction++;
			 break;
			 }
			 tempProbListElem = tempProbListElem->next;
			 }while( tempProbListElem != actionList->head);
			 switch( GetAgent(voronoiDiagram->voronoiCells[j])->state){
			 case FREE:
			 if( GetAgent(voronoiDiagram->voronoiCells[j])->cellCount==0){
			 // Empty cell
			 valid *= (!countDivisionAction
			 && !countNecrosisAction && !countLysisAction && !countVesselGrowthAction
			 );

			 }else{
			 // Partly occupied cell
			 valid *= ( countDivisionAction
			 && countNecrosisAction && countLysisAction && !countVesselGrowthAction
			 );
			 }
			 valid *= (!countDivisionAction
			 && !countNecrosisAction && !countLysisAction && !countVesselGrowthAction
			 );
			 break;

			 case ACTIVE:
			 valid *= ( countDivisionAction
			 && countNecrosisAction && !countLysisAction && !countVesselGrowthAction
			 );
			 break;

			 case NONACTIVE:
			 valid *= (!countDivisionAction
			 && countNecrosisAction && !countLysisAction && !countVesselGrowthAction
			 );
			 break;

			 case VESSEL:
			 valid *= (!countDivisionAction
			 && !countNecrosisAction && !countLysisAction && countVesselGrowthAction
			 );
			 break;

			 case NECROTIC:
			 valid *= (!countDivisionAction
			 && !countNecrosisAction && countLysisAction && !countVesselGrowthAction
			 );
			 break;
			 }
			 if(!valid){
			 printActionList( actionList);
			 fprintf( stderr, "ERROR: some actions are missing!\n"\
			                 "%s cell %i (cellCount=%i) \ncountDivisionAction: %i\ncountNecrosisAction: %i\ncountLysisAction: %i\ncountVesselGrowthAction: %i\n", 
			 cellTypeToString( GetAgent(voronoiDiagram->voronoiCells[j])->state),
			 voronoiDiagram->voronoiCells[j]->index,
			 GetAgent(voronoiDiagram->voronoiCells[j])->cellCount,
			 countDivisionAction, countNecrosisAction, countLysisAction, countVesselGrowthAction);
			 exit(0);
			 }
			 }*/

			// CONTINUUM MODEL (ECM)
			if (VoronoiCell::ECM_THRESHOLD_QUIESCENCE != 0) {

				//fprintf( stderr, "CONTINUUM MODEL (ECM)\n");


				//double ECM_PRODUCTION_RATE = 0.0005; // 0.0003;
				//double ECM_DEGRADATION = ECM_PRODUCTION_RATE / 0.15;
				double TimeStep = 0.01;
				//double EXP = 10;
				//double diffusionCoefECM = 0.1;

				//flag_actualizeProcessRates = TRUE;

				//double dECM[voronoiDiagram->countVoronoiCells];

				//fprintf( stderr, "update ECM\n");

				for (; TempTime <= Time; TempTime += TimeStep) {

					for (int v = 0; v < voronoiDiagram->countVoronoiCells;
							v++) {
						// DEGRADATION
						dECM[v] =
								-(TimeStep) * VoronoiCell::ECM_DEGRADATION_RATE
										* voronoiDiagram->voronoiCells[v]->ecm;
						//voronoiDiagram->voronoiCells[v]->ecm -=
						//		(TimeStep) * ECM_DEGRADATION * voronoiDiagram->voronoiCells[v]->ecm;

						// DIFFUSION
						for (int n = 0;
								n
										< voronoiDiagram->voronoiCells[v]->countNeighborCells;
								n++) {
							//	dECM[v] += (TimeStep) * diffusionCoefECM * (voronoiDiagram->voronoiCells[v]->neighborCells[n]->ecm - voronoiDiagram->voronoiCells[v]->ecm);
						}
					}

					for (int a = 0; a < agentArray->countActiveAgents; a++) {
						for (int l = 0;
								l < agentArray->agents[a]->countLocations; l++)
							if (agentArray->agents[a]->state == ACTIVE
									|| agentArray->agents[a]->state == NONACTIVE) {

								// PRODUCTION
								dECM[agentArray->agents[a]->location[l]->index] +=
										(TimeStep) * VoronoiCell::ECM_PRODUCTION_RATE / 2.
								//* (pow(0.1,EXP) - pow(agentArray->agents[a]->location[l]->ecm,EXP))/ pow(0.1,EXP)
								;

								for (int n = 0;
										n
												< agentArray->agents[a]->location[l]->countNeighborCells;
										n++) {
									dECM[agentArray->agents[a]->location[l]->neighborCells[n]->index] +=
											(TimeStep) * VoronoiCell::ECM_PRODUCTION_RATE
													/ 2.
													/ (double) agentArray->agents[a]->location[l]->countNeighborCells;
								}

							}
					}

					for (int v = 0; v < voronoiDiagram->countVoronoiCells;
							v++) {
						// DEGRADATION
						voronoiDiagram->voronoiCells[v]->ecm += dECM[v];
					}
				}

				//fprintf( stderr, "CONTINUUM MODEL (ECM) finished\n");
			}
			//flag_actualizeProcessRates = TRUE;

			// OUTPUT


			if(!NoRadialProfiles ||
					( RadialProfilesCount && inbound<double>( RadialProfilesTime, floor(Last_Time/24.), floor(Time/24.), RadialProfilesCount, 1) ) ){

				//fprintf( stderr, "Radial Profile\n");
				// UPDATE
				if( ceil(Last_Time/1) != ceil(Time/1)){
					//***for( int hour=(int)ceil(Last_Time/1.); hour < (int)ceil(Time/1.); hour++ )
					updateHistogram( agentArray, voronoiDiagram,	histogram, histogramDividing, histogramNecrotic, histogramFree, histogramECM);
					//fprintf(stderr, "[Updated History]\n");
				}


				//***if( ceil(Last_Time/24.) != ceil(Time/24.))
				for( int day=(int)ceil(Last_Time/24.); day < (int)ceil(Time/24.); day++ )
				if(	inbound<double>( RadialProfilesTime, day-1, day-1, RadialProfilesCount, 1))
				{
					fprintf(stderr, "[Compare to Data: %i]\n", (int)ceil(Last_Time/24.));

					// data_KI67
					comparison_t *data_KI67;
					switch( day-1){
					case 17: data_KI67 = &data_KI67_17; break;
					case 24: data_KI67 = &data_KI67_24; break;
					}

					if(true && data_KI67->dim){
						//fprintf(stderr, "[Compare to KI67 Data]\n");

						comparison_t sim_KI67 = create_comparison();
						sim_KI67.dim = (int) ceil( max( data_KI67->x, data_KI67->dim) );
						sim_KI67.x = (double*) malloc( sizeof(double) * sim_KI67.dim);
						sim_KI67.m = (double*) malloc( sizeof(double) * sim_KI67.dim);
						FILE *fp_raw = fopen( "raw_ki67.dat", "w+");
						for( int j=0; j<sim_KI67.dim; j++){
							sim_KI67.x[j] = j;
							sim_KI67.m[j] = ( histogram[j] - histogramFree[j]>0 ? (double) histogramDividing[j] / (double) (histogram[j] - histogramFree[j]) : 0);
							sim_KI67.m[j] = max( 0.,  sim_KI67.m[j] + measurement_error_ki67*normrnd());
							fprintf( fp_raw, "%e %e %e %e\n", data_KI67->x[j], sim_KI67.m[j], data_KI67->m[j], data_KI67->s[j]);
						}
						fclose(fp_raw);

						cumEpsilon += compare( *data_KI67, sim_KI67, mean_vs_mean) / data_KI67->dim;
						//fprintf(stderr, "<< %e >>\n", cumEpsilon);
						if( isnan(cumEpsilon)){
							fprintf(stderr, "{data_KI67}!\n"); exit(0);
						}

						if( cumEpsilon > maxEpsilon)
							{fprintf(stderr, "break ki67 %i!\n", day-1);	return cumEpsilon;}
						fprintf(stderr, "[Compared to KI67 data]\n");
					}

					// data_ECM
					comparison_t *data_ECM;
					switch( day-1){
					case 17: data_ECM = &data_ECM_17; break;
					case 24: data_ECM = &data_ECM_24; break;
					}
					if(true && data_ECM->dim){
						comparison_t sim_ECM = create_comparison();
						sim_ECM.dim = (int) ceil( max( data_ECM->x, data_ECM->dim) );
						sim_ECM.x = (double*) malloc( sizeof(double) * data_ECM->dim);
						sim_ECM.m = (double*) malloc( sizeof(double) * data_ECM->dim);
						FILE *fp_raw = fopen( "raw_ecm.dat", "w+");
						for( int j=0; j<sim_ECM.dim; j++){
							sim_ECM.x[j] = j;
							if( histogram[j])
								sim_ECM.m[j] = (double) histogramECM[j]	/ (double) histogram[j];
							else
								sim_ECM.m[j] = 0;
							sim_ECM.m[j] = max( 0.,  sim_ECM.m[j] + measurement_error_ecm*normrnd());
							fprintf( fp_raw, "%e %e %e %e\n", data_ECM->x[j], sim_ECM.m[j], data_ECM->m[j], data_ECM->s[j]);
						}
						fclose(fp_raw);

						cumEpsilon += compare( *data_ECM, sim_ECM, mean_vs_mean) / data_ECM->dim;
						//fprintf(stderr, "<< %e >>\n", cumEpsilon);

						if( isnan(cumEpsilon)){
							fprintf(stderr, "{data_ECM}!\n"); exit(0);
						}

						if( cumEpsilon > maxEpsilon)
							{fprintf(stderr, "break ecm %i!\n", day-1);	return cumEpsilon;}
						fprintf(stderr, "[Compared to ECM data]\n");
					}

					// data_TUNEL
					comparison_t *data_TUNEL;
					switch( day-1){
					case 17: data_TUNEL = &data_TUNEL_17; break;
					case 24: data_TUNEL = &data_TUNEL_24; break;
					}
					if(true && data_TUNEL->dim){
						comparison_t sim_TUNEL = create_comparison();
						fprintf(stderr, "tunel size: %i\n", data_TUNEL->dim);
						sim_TUNEL.dim = (int) ceil( max( data_TUNEL->x, data_TUNEL->dim) );
						sim_TUNEL.x = (double*) malloc( sizeof(double) * data_TUNEL->dim);
						sim_TUNEL.m = (double*) malloc( sizeof(double) * data_TUNEL->dim);
						fprintf(stderr, "1\n");
						FILE *fp_raw_t = fopen( "raw_TUNEL.dat", "w+");
						perror( "bla");
						if( fp_raw_t==NULL)
							fprintf(stderr, "raw_TUNEL.dat couldn't be opened\n");
						for( int j=0; j<sim_TUNEL.dim-1; j++){
							sim_TUNEL.x[j] = j;
							sim_TUNEL.m[j] = ( histogram[j] - histogramFree[j]>0 ? (double) histogramNecrotic[j] / (double) (histogram[j] - histogramFree[j]) : 0);
							sim_TUNEL.m[j] = max( 0.,  sim_TUNEL.m[j] + measurement_error_TUNEL*normrnd());

							fprintf( stderr, "%e %e %e %e\n", data_TUNEL->x[j], sim_TUNEL.m[j], data_TUNEL->m[j], data_TUNEL->s[j]);
							fprintf( fp_raw_t, "%e %e %e %e\n", data_TUNEL->x[j], sim_TUNEL.m[j], data_TUNEL->m[j], data_TUNEL->s[j]);
						}
						fprintf(stderr, "2\n");

						fclose(fp_raw_t);
						fprintf(stderr, "3\n");

						cumEpsilon += compare( *data_TUNEL, sim_TUNEL, mean_vs_mean) / data_TUNEL->dim;
						//fprintf(stderr, "<< %e >>\n", cumEpsilon);
						fprintf(stderr, "4\n");

						if( isnan(cumEpsilon)){
							fprintf(stderr, "{data_TUNEL}!\n"); exit(0);
						}

						if( cumEpsilon > maxEpsilon)
							{fprintf(stderr, "break tunel %i!\n", day-1);	return cumEpsilon;}

						fprintf(stderr, "[Compared to TUNEL data]\n");
					}


					// WRITE TO FILE
					fprintf(stderr, "[Write Radial Profile]\n");
					//fprintf(stderr, "[Print Cells to Povray File]\n");
					//***sprintf(outfilename, "%s/rel%i.radialProfiles_day%i-%i.pov", dirname, k, (int)floor(Last_Time/24.), (int)ceil(Last_Time/24.));
					sprintf(outfilename, "%s/rel%i.radialProfiles_day%i-%i.pov", dirname, k, day-1, day);

					if(PovrayOutput)
						agentArray->printToPovray(outfilename, voronoiDiagram);
					//fprintf(stderr, "[Wrote Povray File]\n");
					char filename2[512];
					sprintf(filename2, "%s.dat", outfilename);
					fp = fopen(filename2, "w+");
					for (int i = 0; i < 400; i++)
						//if (histogram[i])
						{
							fprintf(
									fp,
									"%i %i %i %lf %i %lf %e %lf %lf %i\n",
									i,
									histogram[i],
									histogramDividing[i],
									(double) histogramDividing[i] / (double) (histogram[i] - histogramFree[i]),
									histogramNecrotic[i],
									(double) histogramNecrotic[i] / (double) (histogram[i] - histogramFree[i]),
									histogramECM[i],
									( histogram[i] == 0 ? 0 : (double) histogramECM[i] / (double) histogram[i]),
									PROB_REENTERING_CELL_CYCLE(i),
									histogramFree[i]);
							//fprintf(fs, "%i %i %i\n", i, histogram[i],histogramDividing[i]);
						}
					fclose(fp);
					fprintf(stderr, "[Wrote Radial Profile]\n");

					// RESET
					for( int i=0; i<10000; i++){
						histogram[i]=0;
						histogramDividing[i]=0;
						histogramNecrotic[i]=0;
						histogramFree[i]=0;
						histogramECM[i]=0;
					}
					fprintf(stderr, "[Reset History]\n");
				}
				fprintf( stderr, "Radial Profile finished\n");
			}



			// actual timestep
			if (indexOfTime( Last_Time, BeginningTime, OutputRate)
					< indexOfTime( Time, BeginningTime, OutputRate)) {

				if(data_growthcurve.dim){
					// passed data points
					for( ; data_growthcurve.x[idxEpsilon] < Time && idxEpsilon<data_growthcurve.dim; idxEpsilon++){
						double radius = ( isnan(gyrRadius) ? sqrt(DIMENSIONS*50*50) : sqrt(gyrRadius) ) * AGENT_DIAMETER * 0.8;

						data_growthcurve.m[idxEpsilon] = max( 0.,  data_growthcurve.m[idxEpsilon] + measurement_error_gc*normrnd());

						cumEpsilon += 0.5 * pow( (data_growthcurve.m[idxEpsilon] - radius)/data_growthcurve.s[idxEpsilon], 2)  / data_growthcurve.dim;
						fprintf(stderr, "[[ %e ]]\n", cumEpsilon);

						FILE *fp_raw = fopen( "raw_gc.dat", "a+");
						fprintf( fp_raw, "%e %e %e %e\n", data_growthcurve.x[idxEpsilon], radius, data_growthcurve.m[idxEpsilon], data_growthcurve.s[idxEpsilon]);
						fclose(fp_raw);

						if( isnan(cumEpsilon)){
							fprintf(stderr, "{GC}!\n"); exit(0);
						}

						if( cumEpsilon > maxEpsilon)
							{fprintf(stderr, "break gc!\n");	return cumEpsilon;}
					}
				}

				Last_gyrRadius = gyrRadius;
				gyrRadius = getGyrationRadiusOfBorder(agentArray, voronoiDiagram);

				{
					int l;

					// SINGLE OUTPUT
					sprintf(outfilename, "%s/rel%i.data.dat", dirname, k);
					fp = fopen(outfilename, "a+");
					for (l = indexOfTime( Last_Time, BeginningTime, OutputRate) + 1;
						 l < indexOfTime( Time, BeginningTime, OutputRate);
						 l++) {
						//fprintf(stderr, "< [%i, %f]\n", l, timeOfIndex ( l, BeginningTime, OutputRate));
						fprintf(	fp, "%e %i %e\n",
							timeOfIndex ( l, BeginningTime, OutputRate),
							Last_NumberOfCells,
							( isnan(Last_gyrRadius) ? sqrt(DIMENSIONS*50*50) : sqrt(Last_gyrRadius) ) * AGENT_DIAMETER * 0.8
							);
					}
					l = indexOfTime( Time, BeginningTime, OutputRate);
					//fprintf(stderr, "= [%i, %f], %f\n", l, timeOfIndex ( l, BeginningTime, OutputRate), Time);

					fprintf(	fp, "%e %i %e\n",
							timeOfIndex ( l, BeginningTime, OutputRate),
							count_cells,
							( isnan(gyrRadius) ? sqrt(DIMENSIONS*50*50) : sqrt(gyrRadius) ) * AGENT_DIAMETER * 0.8
						);
					fclose(fp);
				}



				// fill passed time steps
				for (l = indexOfTime( Last_Time, BeginningTime, OutputRate) + 1;
						l < indexOfTime( Time, BeginningTime, OutputRate);
						l++) {
					//fprintf( stderr, "fill\n");
					//global_prob[t][c]
#if PROB_OUTPUT
					if(Last_NumberOfCells<MaxCells)
					global_prob[l][Last_NumberOfCells]++;
#endif
					if (OutputAnimation == TRUE) {
						*pnm_suffix = '\0';
						sprintf(pnmstring, "%05d", l);
						strcat(pnm_suffix, pnmstring);
						sprintf(outfilename, "%s%canim_", dirname, SEPERATOR);
						strcat(outfilename, pnm_suffix);
						strcat(outfilename, ".pov");
						// print_cell_structure_to_pov_file(outfilename,voronoiGrid, NumberOfVoronoiPoints, radiusOfScene, centerOfScene, radiusOfCell);
						printToPovray(
								outfilename,
								voronoiDiagram,
								2 * PI * (double) l / (double) EndTime_int
										* rotations + startAngle, setoffX,
								setoffY, setoffZ);
					}

					// update arrays
					global_cells[l] += Last_NumberOfCells; //count_cells;
					global_volume[l] += count_cell_volume;

					global_necrotic_cells[l] += count_necrotic_cells;
					global_necrotic_volume[l] += count_necrotic_volume;

					global_inner_empty_volume[l] += count_inner_empty_volume;

					global_vessel_cells[l] += count_vessel_cells;

					global_gyrRadius[l] += sqrt(Last_gyrRadius);

					global_gyrRadiusSquare[l] += Last_gyrRadius;

					global_avgGlucose[l] += getAvgGlucose(agentArray);
					global_avgOxygen[l] += getAvgOxygen(agentArray);

#if RADIAL_OUTPUT
					statisticsOfSpheroid( voronoiDiagram, global_oxygen_concentration, global_glucose_concentration, global_active_cell_density, global_nonactive_cell_density, global_necrotic_cell_density, EndRadius, EndRadius_int, l);
#endif
					/*global_prob_X_equal_k[l]   += (Last_NumberOfCells == CellDivisionDepth);
					 global_prob_X_smaller_k[l] += (Last_NumberOfCells <  CellDivisionDepth);
					 global_prob_X_bigger_k[l]  += (Last_NumberOfCells >  CellDivisionDepth);
					 */

#if ON_THE_FLY_OUTPUT


					/*{
						sim_growthcurve.x[l] = timeOfIndex ( l, BeginningTime, OutputRate);
						sim_growthcurve.m[l] = global_gyrRadius[l] / (k + 1.) * AGENT_DIAMETER * 0.8;
						sim_growthcurve.dim = l+1;
					}*/


					sprintf(outfilename, "%s/dataOnTheFly.dat", dirname);
					if ((fp = fopen(outfilename, "a+")) == NULL) {
						fprintf(stderr, "Error opening file %s for writing!\n",
								outfilename);
						//exit(0);
					}
					fprintf(
							fp,
							"%lf %lf %lf %lf %lf %lf %lf %e %e %lf %lf %e %e %e\n",

							// time
							timeOfIndex ( l, BeginningTime, OutputRate),

							// living cells
							global_cells[l] / (k + 1.),

							// cell volume
							global_volume[l] / (k + 1.),

							// necrotic cells
							global_necrotic_cells[l] / (k + 1.),

							// necrotic cells
							global_necrotic_volume[l] / (k + 1.),

							// inner volume
							global_inner_empty_volume[l] / (k + 1.),

							// vessel segments
							global_vessel_cells[l] / (k + 1.),

							global_avgGlucose[l] / (k + 1.), //centralCell->glucose,

							global_avgOxygen[l] / (k + 1.), //centralCell->oxygen,

							global_gyrRadius[l] / (k + 1.),

							global_gyrRadiusSquare[l] / (k + 1.)
									- pow(global_gyrRadius[l] / (k + 1.), 2.),

							global_glucose_concentration,

							global_oxygen_concentration,

							(double) clock()

							);
					fclose(fp);
#endif
				}

				// calculate array position for actual time   
				l = indexOfTime( Time, BeginningTime, OutputRate);

				//sprintf(outfilename, "%s/VoronoiDiagram%i.eps", dirname, l);
				//printVoronoiDiagram( voronoiDiagram, outfilename, false);

				//printTriangulation( voronoiDiagram, "testTriangulation.eps", false);
				//printVoronoiDiagram( voronoiDiagram, "testVoronoiDiagram.eps", false);

#if PROB_OUTPUT
				if(count_cells<MaxCells)
				global_prob[l][count_cells]++;
#endif

				// write povray file
				if (OutputAnimation == TRUE) {
					*pnm_suffix = '\0';
					sprintf(pnmstring, "%05d", l);
					strcat(pnm_suffix, pnmstring);
					sprintf(outfilename, "%s%canim_", dirname, SEPERATOR);
					strcat(outfilename, pnm_suffix);
					strcat(outfilename, ".pov");
					//print_cell_structure_to_pov_file(outfilename,voronoiGrid, NumberOfVoronoiPoints, radiusOfScene, centerOfScene, radiusOfCell);
					printToPovray(
							outfilename,
							voronoiDiagram,
							2 * PI * (double) l / (double) EndTime_int * rotations + startAngle,
							setoffX,
							setoffY,
							setoffZ);
				}

				// update arrays
				global_cells[l] += count_cells;
				global_volume[l] += count_cell_volume;

				global_necrotic_cells[l] += count_necrotic_cells;
				global_necrotic_volume[l] += count_necrotic_volume;

				global_inner_empty_volume[l] += count_inner_empty_volume;

				global_vessel_cells[l] += count_vessel_cells;

				//gyrRadius = getGyrationRadius(agentArray);
				/*if( data_growthcurve.dim && maxRadius < ( isnan(gyrRadius) ? sqrt(DIMENSIONS)*50 : sqrt(gyrRadius) ) * AGENT_DIAMETER * 0.8){
					//fprintf(stderr, "(( max radius exeeded ))\n");
					//return cumEpsilon + maxEpsilon;
					return maxEpsilon;
				}*/

				if( data_growthcurve.dim && maxRadius < sqrt(gyrRadius) * AGENT_DIAMETER * 0.8)
					return maxEpsilon;

				global_gyrRadius[l] += sqrt(gyrRadius);

				global_gyrRadiusSquare[l] += gyrRadius;

				global_avgGlucose[l] += getAvgGlucose(agentArray);
				global_avgOxygen[l] += getAvgOxygen(agentArray);

#if RADIAL_OUTPUT
				statisticsOfSpheroid( voronoiDiagram, global_oxygen_concentration, global_glucose_concentration, global_active_cell_density, global_nonactive_cell_density, global_necrotic_cell_density, EndRadius, EndRadius_int, l);
#endif
#if SLICE_OUTPUT
				if(!NoSliceOutput)
				writeSliceInformation(voronoiDiagram, agentArray, dirname);

#endif
				/*global_prob_X_equal_k[l]   += (count_inoculated_cells + count_divided_cells - count_necrotic_cells - count_lysed_cells == CellDivisionDepth);
				 global_prob_X_smaller_k[l] += (count_inoculated_cells + count_divided_cells - count_necrotic_cells - count_lysed_cells <  CellDivisionDepth);
				 global_prob_X_bigger_k[l]  += (count_inoculated_cells + count_divided_cells - count_necrotic_cells - count_lysed_cells >  CellDivisionDepth);
				 */

#if ON_THE_FLY_OUTPUT
				/*if( global_gyrRadius[l] / (k + 1.) >= 40){
				 fprintf( stderr, "[Print Cells to Povray File]\n");
				 sprintf(outfilename, "%s/cellsR40.pov", dirname);
				 agentArray->printToPovray( outfilename, voronoiDiagram);
				 exit(0);
				 }*/


				if(false)
				if ((sqrt(gyrRadius) * SPATIAL_UNIT >= 300.	&& sqrt(gyrRadius) * SPATIAL_UNIT <= 350)
						|| (Time > 1000 && Time < 2000	&& agentArray->countActiveAgents < 100000)) {
					fprintf(stderr, "[Print Cells to Povray File]\n");
					//Time = EndTime;
					sprintf(outfilename, "%s/cellsR300.pov", dirname);
					agentArray->printToPovray(outfilename, voronoiDiagram);
					//exit(0);

					updateHistogram( agentArray, voronoiDiagram, histogram, histogramDividing, histogramNecrotic, histogramFree, histogramECM);


					// GET ALL BORDER CELLS
/*					VoronoiCell *borderCells[100000];
					int countBorderCells = 0;
					for (int a = 0; countBorderCells!=100000 && a < agentArray->countActiveAgents; a++)
						for (int l = 0; countBorderCells!=100000 && l < agentArray->agents[a]->countLocations; l++){
							bool unoccupiedNeighbor = false;
							for( int n=0; countBorderCells!=100000 && n<agentArray->agents[a]->location[l]->countNeighborCells / *&& !unoccupiedNeighbor* /; n++)
								if( !agentArray->agents[a]->location[l]->neighborCells[n]->agent){
									unoccupiedNeighbor=true;
									//if( countBorderCells==100000)
									//	exit(0);
									borderCells[countBorderCells++] = agentArray->agents[a]->location[l]->neighborCells[n];
								}
						}

					for (int a = 0; a < agentArray->countActiveAgents; a++)
					//if( agentArray->agents[a]->state != FREE)
							{

						// DIVIDING?
						bool dividing = false;
						bool necrotic = false;
						bool free = false;
						for (int l = 0;
								l < agentArray->agents[a]->countLocations;
								l++) {
							if (agentArray->agents[a]->state == FREE)
								free = true;
							if (agentArray->agents[a]->state == NECROTIC)
								necrotic = true;
							if (agentArray->agents[a]->state == ACTIVE
									&& agentArray->agents[a]->divide == 1
									&& !IsQuiescent(agentArray->agents[a]))
								dividing = true;
						}

						// distance to border
						VoronoiCell *border;
						//if( VoronoiCell::SHIFT_TO_UNOCCUPIED)
						/ *if( VoronoiCell::USE_SYMBOLIC_EXTENDED_NEIGHBORHOOD)
							border = voronoiDiagram->searchClosestUnoccupiedVoronoiCell(
								agentArray->agents[a]->location[0],
								VoronoiCell::extendedNeighborDepth,
								VoronoiCell::symbolicExtendedNeighborhoodSize,
								VoronoiCell::symbolicExtendedNeighborhood);
						else* /

							border = borderCells[0];
							double distClosestBorder = border->getDistanceTo( agentArray->agents[a]->location[0]);
							for( int b=1; b<countBorderCells; b++){
								double dist = borderCells[b]->getDistanceTo( agentArray->agents[a]->location[0]);
								if( distClosestBorder>dist){
									border = borderCells[b];
									distClosestBorder=dist;
								}
							}
							//border =
								//// SLOWEST!
								////voronoiDiagram->searchClosestUnoccupiedVoronoiCellNAIVE( agentArray->agents[a]->location[0]);
								// SLOW
								//voronoiDiagram->searchClosestUnoccupiedVoronoiCell(	agentArray->agents[a]->location[0], 100);
						//else
						//	border = voronoiDiagram->searchClosestFreeVoronoiCell( agentArray->agents[a]->location[0], 100);

						if (border) {
							double dist = border->getDistanceTo(
									agentArray->agents[a]->location[0])
									* SPATIAL_UNIT;
							for (int i = (dist - AGENT_DIAMETER > 0. ? (int) (dist - AGENT_DIAMETER) : 0);
									i < (int) (dist); i++) {
								histogram[i]++;
								if (dividing)
									histogramDividing[i]++;
								if (necrotic)
									histogramNecrotic[i]++;
								if (free)
									histogramFree[i]++;
								histogramECM[i] +=
										agentArray->agents[a]->location[0]->ecm;
							}
						} else {
							fprintf(stderr, "No border found!\n");
							exit(0);
						}
					}
*/
					char filename2[512];
					sprintf(filename2, "%s/cellsR300.pov.dat", dirname);
					fp = fopen(filename2, "w+");
					for (int i = 0; i < 10000; i++)
						if (histogram[i]) {
							fprintf(
									fp,
									"%i %i %i %lf %i %lf %lf %lf %i\n",
									i,
									histogram[i],
									histogramDividing[i],
									(double) histogramDividing[i]
											/ (double) (histogram[i] - histogramFree[i]),
									histogramNecrotic[i],
									(double) histogramNecrotic[i]
									        / (double) (histogram[i] - histogramFree[i]),
									(double) histogramECM[i]
											/ (double) histogram[i],
									PROB_REENTERING_CELL_CYCLE(i),
									histogramFree[i]);
							//fprintf(fs, "%i %i %i\n", i, histogram[i],histogramDividing[i]);
						}
					fclose(fp);

					//exit(0);

				}
				/*if(l%10==0){
				 sprintf(outfilename, "%s/cells%i.pov", dirname, l);
				 agentArray->printToPovray( outfilename, voronoiDiagram);
				 }*/
				/*{
					sim_growthcurve.x[l] = timeOfIndex ( l, BeginningTime, OutputRate);
					sim_growthcurve.m[l] = global_gyrRadius[l] / (k + 1.) * AGENT_DIAMETER * 0.8;
					sim_growthcurve.dim = l+1;

					double x_min = sim_growthcurve.x[0];
					double x_max = sim_growthcurve.x[l];

					int new_dim=0;
					for( ; new_dim<data_growthcurve.dim && data_growthcurve.x[new_dim] <= x_max; new_dim++) ;
					int old_dim = data_growthcurve.dim;
					data_growthcurve.dim = new_dim;

					// evaluate Epsilon
					if(l+1 > 5){
						//fprintf( stderr, "[ Do it: %i ]\n", sim_growthcurve.dim);
						//printVector( sim_growthcurve.m, sim_growthcurve.dim, "%.3e ");
						double negLogLik = compare( data_growthcurve, sim_growthcurve, mean_vs_mean);
						fprintf( stderr, "[ %i, %i, Epsilon = %e ]\n", data_growthcurve.dim, sim_growthcurve.dim, negLogLik);
					}
					data_growthcurve.dim = old_dim;

				}*/

				//fprintf( stderr, "outut!\n");
				sprintf(outfilename, "%s/dataOnTheFly.dat", dirname);
				if ((fp = fopen(outfilename, "a+")) == NULL) {
					fprintf(stderr, "Error opening file %s for writing!\n",
							outfilename);
					//exit(0);
				}
				fprintf(
						fp,
						"%lf %lf %lf %lf %lf %lf %lf %e %e %lf %lf %e %e %e\n",

						// time
						timeOfIndex ( l, BeginningTime, OutputRate),

						// living cells
						global_cells[l] / (k + 1.),

						// cell volume
						global_volume[l] / (k + 1.),

						// necrotic cells
						global_necrotic_cells[l] / (k + 1.),

						// necrotic cells
						global_necrotic_volume[l] / (k + 1.),

						// inner volume
						global_inner_empty_volume[l] / (k + 1.),

						// vessel segments
						global_vessel_cells[l] / (k + 1.),

						global_avgGlucose[l] / (k + 1.), //centralCell->glucose,

						global_avgOxygen[l] / (k + 1.), //centralCell->oxygen,

						global_gyrRadius[l] / (k + 1.),

						global_gyrRadiusSquare[l] / (k + 1.)
								- pow(global_gyrRadius[l] / (k + 1.), 2.),

						global_glucose_concentration,

						global_oxygen_concentration,

						(double) clock()

						);
				fclose(fp);
#endif

				// END OF OUTPUT
#if _COMMENTS_ > 0
				/*fprintf( stderr, "\r%.3f%% (%d/%d) (%.1lfdays, volume: %d, tumor cells: %d, necrotic: %d, central glu: %e oxy: %e)\t\t",
				 (Time * 100)/EndTime,
				 k+1,
				 Averages,
				 Time/24.,
				 count_cell_volume,
				 count_cells,
				 count_necrotic_cells,
				 centralCell->glucose,
				 centralCell->oxygen//,
				 //actionList->size//,
				 //agentArray->countActiveAgents
				 );*/
#endif

				/*global_cells[l]  += count_cells;
				 global_volume[l] += count_cell_volume;

				 global_necrotic_cells[l]  += count_necrotic_cells;
				 global_necrotic_volume[l] += count_necrotic_volume;

				 global_inner_empty_volume[l] += count_inner_empty_volume;*/

			}

			fprintf(stderr, "\r > %d/%d: time %.3lf/%.3lf", k + 1, Averages,
					Time, EndTime);
			fprinttextcolor(GREEN);
			fprintf(stderr, "[CELLS:%i] ", count_cells);
			fprinttextcolor(RED);
			fprintf(stderr, "[DEAD:%i]\t\t", count_necrotic_cells);
			//fprinttextcolor(BLACK);
			fprinttextattribute(RESET);

		} while (Time < EndTime && actionList->size != 0);

		fprintf( stderr, "Simulation finished!\n");

		// RADIAL PROFILE AVERAGING
		if(!NoRadialProfiles){
			char filename2[512];
			sprintf(filename2, "%s/cellsR300.average.dat", dirname);
			fp = fopen(filename2, "w+");

			for (int i = 0; i < 10000; i++) {

				// SUM UP VALUES & SQUARES
				histogramSumDividing[i] += histogramDividing[i] / histogram[i];
				histogramSumNecrotic[i] += histogramNecrotic[i] / histogram[i];
				histogramSumDividingFraction[i] += histogramDividing[i] / (histogram[i] - histogramFree[i]);
				histogramSumNecroticFraction[i] += histogramNecrotic[i] / (histogram[i] - histogramFree[i]);
				histogramSumFree[i]     += histogramFree[i]     / histogram[i];
				histogramSumECM[i]      += histogramECM[i]      / histogram[i];

				histogramSquaresDividing[i] += pow( histogramDividing[i] / histogram[i], 2.);
				histogramSquaresNecrotic[i] += pow( histogramNecrotic[i] / histogram[i], 2.);
				histogramSquaresDividingFraction[i] += pow( histogramDividing[i] / (histogram[i] - histogramFree[i]), 2.);
				histogramSquaresNecroticFraction[i] += pow( histogramNecrotic[i] / (histogram[i] - histogramFree[i]), 2.);
				histogramSquaresFree[i]     += pow( histogramFree[i]     / histogram[i], 2.);
				histogramSquaresECM[i]      += pow( histogramECM[i]      / histogram[i], 2.);

				// RESET
				histogram[i] = 0;

				histogramDividing[i] = 0;
				histogramNecrotic[i] = 0;
				histogramFree[i] = 0;
				histogramECM[i] = 0.;

				if (histogram[i]) {
					fprintf(
							fp,
							"%i %lf %lf %lf %lf %lf %lf\n",
							i,
							histogramSumDividing[i] / (double)k,
							histogramSumNecrotic[i] / (double)k,
							histogramSumECM[i]  / (double)k,
							sqrt( histogramSquaresDividingFraction[i] / (double)k - pow(histogramSumDividingFraction[i] / (double)k, 2.)),
							sqrt( histogramSquaresNecroticFraction[i] / (double)k - pow(histogramSumNecroticFraction[i] / (double)k, 2.)),
							sqrt( histogramSquaresECM[i] / (double)k - pow(histogramSumECM[i] / (double)k, 2.)));
				}
			}

			fclose(fp);
		}


//		fprintf(stderr, "\r > %d/%d\t\t", k + 1, Averages);

		//fprintf(stderr, "\nbefore fill following: Time=%lf, EndTime=%lf\n",
		//		Time, EndTime);
		// fill following time step
		if(false)
		for (l = indexOfTime( Time, BeginningTime, OutputRate) + 1;
				l < nrOfSteps( BeginningTime, EndTime, OutputRate); l++) {

			fprintf(stderr, "fill following: \n");

			if (OutputAnimation == TRUE) {
				*pnm_suffix = '\0';
				sprintf(pnmstring, "%05d", l);
				strcat(pnm_suffix, pnmstring);
				sprintf(outfilename, "%s%canim_", dirname, SEPERATOR);
				strcat(outfilename, pnm_suffix);
				strcat(outfilename, ".pov");
				//print_cell_structure_to_pov_file(outfilename,voronoiGrid, NumberOfVoronoiPoints, radiusOfScene, centerOfScene, radiusOfCell);
				printToPovray(
						outfilename,
						voronoiDiagram,
						2 * PI * (double) l / (double) EndTime_int * rotations
								+ startAngle, setoffX, setoffY, setoffZ);
			}

			// update arrays
			global_cells[l] += count_cells;
			global_volume[l] += count_cell_volume;

			global_necrotic_cells[l] += count_necrotic_cells;
			global_necrotic_volume[l] += count_necrotic_volume;

			global_inner_empty_volume[l] += count_inner_empty_volume;

			global_vessel_cells[l] += count_vessel_cells;

			global_gyrRadius[l] += sqrt(gyrRadius);

			global_gyrRadiusSquare[l] += gyrRadius;

			global_avgGlucose[l] += getAvgGlucose(agentArray);
			global_avgOxygen[l] += getAvgOxygen(agentArray);

#if RADIAL_OUTPUT
			statisticsOfSpheroid( voronoiDiagram, global_oxygen_concentration, global_glucose_concentration, global_active_cell_density, global_nonactive_cell_density, global_necrotic_cell_density, EndRadius, EndRadius_int, l);
#endif
			/*global_prob_X_equal_k[l]   += (count_inoculated_cells + count_divided_cells - count_necrotic_cells - count_lysed_cells == CellDivisionDepth);
			 global_prob_X_smaller_k[l] += (count_inoculated_cells + count_divided_cells - count_necrotic_cells - count_lysed_cells <  CellDivisionDepth);
			 global_prob_X_bigger_k[l]  += (count_inoculated_cells + count_divided_cells - count_necrotic_cells - count_lysed_cells >  CellDivisionDepth);
			 */
		}

		/** END OF SIMULATION **/

#if GLOBAL_OUTPUT

		/** GLOBAL OUTPUT **/
		//sprintf(outfilename,"%s/grid_data_cells_vs_time_average.dat",dirname);
		sprintf(outfilename, "%s/data.dat", dirname);

		if ((fp = fopen(outfilename, "w")) == NULL) {
			fprintf(stderr, "Error opening file %s for writing!\n",
					outfilename);
			//exit(0);
		}
		fprintf(fp, "# %i realisations\n", k + 1);
		fprintf(fp,
				"# rows: time, tumor cells, necrotic cells, vessel segments\n");
		EndTime_int = nrOfSteps( BeginningTime, EndTime, OutputRate);
		for (l = 0; l < EndTime_int; l++) {
			fprintf(
					fp,
					"%lf %lf %lf %lf %lf %lf %lf %e %e %lf %lf\n",

					// time
					timeOfIndex ( l, BeginningTime, OutputRate),

					// living cells
					global_cells[l] / (k + 1.),

					// cell volume
					global_volume[l] / (k + 1.),

					// necrotic cells
					global_necrotic_cells[l] / (k + 1.),

					// necrotic cells
					global_necrotic_volume[l] / (k + 1.),

					// inner volume
					global_inner_empty_volume[l] / (k + 1.),

					// vessel segments
					global_vessel_cells[l] / (k + 1.),

					global_avgGlucose[l] / (k + 1.), //centralCell->glucose,

					global_avgOxygen[l] / (k + 1.), //centralCell->oxygen,

					global_gyrRadius[l] / (k + 1.),

					global_gyrRadiusSquare[l] / (k + 1.)
							- pow(global_gyrRadius[l] / (k + 1.), 2.)

					);
		}
		fclose(fp);

#if RADIAL_OUTPUT
		sprintf(outfilename,"%s/radialData.dat",dirname);

		if((fp=fopen(outfilename, "w"))==NULL) {
			fprintf(stderr, "Error opening file %s for writing!\n",outfilename);
			//exit(0);
		}
		//fprintf( fp,"# %i realisations\n", k+1);
		for(l = 0;l < EndTime_int; l++) {
			for(i = 0;i < EndRadius_int; i++) {
				fprintf( fp,"%lf %lf %lf %lf %lf %lf %lf \n",
						// time
						timeOfIndex ( l, BeginningTime, OutputRate),

						// radius
						EndRadius * (double)i/(double)EndRadius_int,

						global_oxygen_concentration[l][i] / (k+1.),
						global_glucose_concentration[l][i] / (k+1.),
						global_active_cell_density[l][i] / (k+1.),
						global_nonactive_cell_density[l][i] / (k+1.),
						global_necrotic_cell_density[l][i] / (k+1.));
			}
			//fprintf( fp,"\n");
		}
		fclose(fp);
#endif

		/** END GLOBAL OUTPUT **/

#endif // GLOBAL_OUTPUT
		// POST ANALYSIS

		sprintf(outfilename, "%s/post.dat", dirname);

		if ((fp = fopen(outfilename, "w")) == NULL) {
			fprintf(stderr, "Error opening file %s for writing!\n",
					outfilename);
			//exit(0);
		}
#if DIMENSIONS == 1	
		int glidingAvLength = (int)(5000. / OutputRate),
#elif DIMENSIONS == 2
		int glidingAvLength = (int) (500. / OutputRate),
#elif DIMENSIONS == 3
				int glidingAvLength = (int) (500. / OutputRate),
#endif
				glidingAvPos = 0;
		double glidingAv[glidingAvLength];
		for (i = 0; i < glidingAvLength; i++)
			glidingAv[i] = 0.;
		double avRate = 0.;

		EndTime_int = nrOfSteps( BeginningTime, EndTime, OutputRate);
		for (l = 0; l < EndTime_int - 1; l++) {
			// rate

#if DIMENSIONS == 1	
			double rate = ( global_cells[l+1] / 2. / (double)(k+1.) - global_cells[l] / 2. / (double)(k+1.))
			/ ( timeOfIndex( l+1, BeginningTime, OutputRate) - timeOfIndex( l, BeginningTime, OutputRate));
#elif DIMENSIONS == 2
			double rate = (sqrt(global_cells[l + 1] / (double) (k + 1.) / PI)
					- sqrt(global_cells[l] / (double) (k + 1.) / PI))
					/ (timeOfIndex( l+1, BeginningTime, OutputRate)
							- timeOfIndex( l, BeginningTime, OutputRate));
#elif DIMENSIONS == 3
			//double rate = ( pow( global_cells[l+1] / (double)(k+1.) * 3./4./PI, 1./3.)  - pow( global_cells[l] / (double)(k+1.) * 3./4./PI, 1./3.))
			//            / ( timeOfIndex( l+1, BeginningTime, OutputRate)    - timeOfIndex( l, BeginningTime, OutputRate));
			double rate = (global_cells[l + 1] / 2. / countPoints[1]
					/ countPoints[2] / (double) (k + 1.) - global_cells[l] / 2.
					/ countPoints[1] / countPoints[2] / (double) (k + 1.))
			/ (timeOfIndex( l+1, BeginningTime, OutputRate)
					- timeOfIndex( l, BeginningTime, OutputRate));
#endif

			glidingAv[glidingAvPos] = rate;
			glidingAvPos = (glidingAvPos + 1) % glidingAvLength;
			avRate = 0.;
			for (i = 0; i < glidingAvLength; i++)
				avRate += glidingAv[i];
			avRate /= (double) glidingAvLength;

			fprintf(
					fp, "%lf %lf %lf\n",

					// time
					(timeOfIndex( l, BeginningTime, OutputRate)
							+ timeOfIndex( l+1, BeginningTime, OutputRate)) / 2.,

					// growth rate cells
					rate,

					// gliding average growth rate
					avRate);
		}
		fclose(fp);

		sprintf(outfilename, "%s/growthrate.dat", dirname);
		if ((fp = fopen(outfilename, "w")) == NULL)
			fprintf(stderr, "Error opening file %s for writing!\n",
					outfilename);
		fprintf(fp, "%i %i %i %lf\n", M_gro + 1, M_div + 1,
				(int) CellDivisionDepth, avRate);
		fclose(fp);

		// END POST ANALYSIS

		// FREE MEMORY
#ifndef __WithoutGUI__
		//substrate.SwitchOffVisualization();
#endif
		//fprintf(stderr,"DESTROY ACTION LIST\n");
		// DESTROY ACTION LIST
		//destroyActionList( actionList);
		actionList->destroyActionTree();

		//fprintf(stderr,"DETACH, DEACTIVATE AND REINIT AGENTS\n");
		// DETACH, DEACTIVATE AND REINIT AGENTS
		while (agentArray->countActiveAgents != 0) {
			agentArray->agents[0]->detach();
			agentArray->deactivateAgent(agentArray->agents[0]);
		}
		for (i = 0; i < agentArray->countAgents; i++)
			if (agentArray->agents[i]->actionsInitialized)
				for (int ii = 0; ii <= INDEX_GROWTH; ii++)
					agentArray->agents[i]->actions[ii]->internalState = 0;

		// REINIT VORONOI CELLS
		for (i = 0; i < voronoiDiagram->countVoronoiCells; i++) {
			voronoiDiagram->voronoiCells[i]->countFreeNeighborCells =
					voronoiDiagram->voronoiCells[i]->countNeighborCells;
			voronoiDiagram->voronoiCells[i]->countFreeExtendedNeighborCells =
					voronoiDiagram->voronoiCells[i]->countExtendedNeighborCells;
		}

	} // END AVERAGES

	//fprintf(stderr,"END AVERAGES\n");
	//// PROBABILITY

#if PROB_OUTPUT
	sprintf(outfilename,"%s/probData.dat",dirname);

	if((fp=fopen(outfilename, "w"))==NULL) {
		fprintf(stderr, "Error opening file %s for writing!\n",outfilename);
		//exit(0);
	}
	//fprintf( fp,"# %i realisations\n", k+1);
	for(l = 0;l < EndTime_int; l++) {
		for(c = 0;c < MaxCells; c++) {
			if(global_prob[l][c]!=0)
			fprintf( fp,"%lf %i %lf\n",
	// time
					timeOfIndex ( l, BeginningTime, OutputRate),

					// radius
					c,

					(double)global_prob[l][c] / (k+1.));
		}
		//fprintf( fp,"\n");
	}
	fclose(fp);
#endif
	//// END PROBABILITY

	/** END OF STOCHASTIC APPROACHE ******************************************/

	/* Get the current time. */
	temp_t = time(NULL);
	loctime = localtime_r(&temp_t, loctime);
	strftime(local_time, 256, "%A, %B %d. %k:%M:%S ", loctime);
	free(loctime);

#if SINGLE_OUTPUT

	if((fp4=fopen(parameterfilename, "a"))==NULL) {
		fprintf(stderr, "Error opening file %s for writing!\n",parameterfilename);
		exit(0);
	}
	fprintf(fp4,"\n*********************************\n");
	fprintf(fp4,"Simulation ended the: %s",local_time);
	fclose(fp4);

#endif

/*	fprintf(stderr, "\n*********************************\n\n");
	fprintf(stderr, "Simulation ended the: %s\n", local_time);
	fprintf(stderr, "%d cells were divided\n", (int) count_divided_cells);
*/
	//free(global_roughness);

	// FREE MEMORY
	//exit(0);
	{fprintf(stderr, "break !\n");	return cumEpsilon;}

#if MONOD_KINETICS
#ifdef __WithoutGUI__
	//The substrate will be destroy at the end of this scope.
#else
	//Here we need to kill the visualization to prevent any errors
	//substrate.SwitchOffVisualization();
#endif
#endif

	// Agents	

	//fprintf( stderr, "INDEX_MIGRATION + USE_MIGRATION = %i\n", INDEX_MIGRATION + USE_MIGRATION);

	for (i = 0; i < agentArray->countAgents; i++) {
		if (agentArray->agents[i]->actionsInitialized) {
			//for( k=0; k<INDEX_GROWTH+1; k++)
			for (k = 0; k < INDEX_MIGRATION + USE_MIGRATION; k++)
				free(agentArray->agents[i]->actions[k]);
#if USE_MIGRATION
			//for( k=0; k < GetVoronoiCell( agentArray->agents[i])->countNeighborCells; k++)
			//	free( agentArray->agents[i]->actions[k+INDEX_MIGRATION]);
#endif
			free(agentArray->agents[i]->actions);
			agentArray->agents[i]->actionsInitialized = FALSE;
		}
		if (agentArray->agents[i]->location != NULL)
			free(agentArray->agents[i]->location);
		free(agentArray->agents[i]);
	}
	free(agentArray->agents);
	free(agentArray);

	// Voronoi Grids

	for (i = 0; i < voronoiDiagram->countVoronoiCells; i++) {
		//free( GetAgent( voronoiDiagram->voronoiCells[i])->neighbors);
		free(voronoiDiagram->voronoiCells[i]->extendedNeighborhood);
	}
	//free( voronoiDiagram);
	//voronoiDiagram->deleteVoronoiDiagram();

	// Statistical Variables
	free(global_cells);
	free(global_necrotic_cells);
	free(global_necrotic_volume);
	free(global_vessel_cells);
	free(global_volume);
	free(global_inner_empty_volume);
	free(global_gyrRadius);
	free(global_gyrRadiusSquare);
	free(global_avgOxygen);
	free(global_avgGlucose);

#if PROB_OUTPUT
	EndTime_int = nrOfSteps( BeginningTime, EndEndTime, OutputRate);
	for(int t = 0;t < EndTime_int; t++)
	free( global_prob[t]);
	free( global_prob);
#endif

#if RADIAL_OUTPUT
	EndTime_int = nrOfSteps( BeginningTime, EndEndTime, OutputRate);
	for(k = 0;k < EndTime_int; k++) {
		free( global_oxygen_concentration[k]);
		free( global_glucose_concentration[k]);
		free( global_active_cell_density[k]);
		free( global_nonactive_cell_density[k]);
		free( global_necrotic_cell_density[k]);
	}
	free( global_oxygen_concentration);
	free( global_glucose_concentration);
	free( global_active_cell_density);
	free( global_nonactive_cell_density);
	free( global_necrotic_cell_density);
#endif

	// Free Continuum Part
	if (Case > 1) {
		sA->deleteSparseMatrix(sA);
		//J->deleteSparseMatrix(J);
		free(b);
		free(x);
		free(v0);
		free(v1);
		free(v2);
		free(v3);
		free(v4);
		free(v5);
		free(v6);
		free(v7);
		free(v8);
	}

	//This has been put into comment by William because of an error of segmentation fault at the end of the execution.
	//free( global_prob_X_equal_k);
	//free( global_prob_X_smaller_k);
	//free( global_prob_X_bigger_k);

	//free( centerOfScene);

	fprintf(
			stderr,
			"Total Computation Time:\nActualization -> %.3lfsec\nSelect -> %.3lfsec\nExecution -> %.3lfsec\n",
			timeActualization, timeSelect, timeExecution);

	fprintf(stderr, "Total CPU Time for Reaction-Diffusion-Solver ( %lisec)\n",
			(benchmarkTime) / CLOCKS_PER_SEC);

	return 1;
}
/****************************************************************************/

double DESubstrate(double substrate, double cells) {
	double cellGrowthRate = MaxCellDivisionRate, Ks = 0.2, // (mmol / L)
			ms = 0.0056, // (mmol / L * h * 10^5 cells)
			Ys = 1.37; // ()

	double apparentGrowthRate = cellGrowthRate * substrate / (Ks + substrate),
			qs = apparentGrowthRate / Ys + ms;

	return -qs * (double) cells; // / 100000;
	//return -( (double)substrate * cellGrowthRate / (Ks + (double)substrate) / Ys + ms ) * (double)cells;
}
/****************************************************************************/

double rungeKutta(double(*equation)(double, double), double xStart, double xEnd,
		double yStart, double z, int steps) {
	double stepSize = (xEnd - xStart) / (double) steps;
	int i;

	for (i = 0; i < steps; i++) {
		yStart = rungeKuttaStep(equation, xStart, xStart + stepSize, yStart, z);
	}

	return yStart;
}
/****************************************************************************/

double rungeKuttaStep(double(*equation)(double, double), double xStart,
		double xEnd, double yStart, double z) {
	double stepSize = xEnd - xStart;
	double yA, yB, yC, dyStart, dyA, dyB, dyC;

	dyStart = equation(yStart, z);

	yA = yStart + stepSize / 2. * dyStart;
	dyA = equation(xStart + stepSize / 2., yA);

	yB = yStart + stepSize / 2. * dyA;
	dyB = equation(xStart + stepSize / 2., yB);

	yC = yStart + stepSize * dyB;
	dyC = equation(xStart + stepSize, yC);

	return yStart + stepSize / 6. * (dyStart + 2. * (dyA + dyB) + dyC);
}
/****************************************************************************/

/*double getEnergyDifferenceNeighbors( VoronoiCell *oldPosition, VoronoiCell *newPosition)
 {
 double energyBilance = 0;

 // check for each neighbor of old position
 for( int i=0; i<oldPosition->countNeighborCells; i++){
 // loss of contact
 if( oldPosition->neighborCells[i] != newPosition)
 energyBilance -= 1;
 }

 // check for each neighbor of old position
 for( int i=0; i<newPosition->countNeighborCells; i++){
 // win contact
 if( newPosition->neighborCells[i] != oldPosition)
 energyBilance += 1;
 }

 return energyBilance;
 }*/

#define INTRA_AGENT_ENERGY 4
#define INTER_AGENT_ENERGY 2 //0.5
////////////////////////////////////////
double getEnergyDifferenceCadherin(VoronoiCell *oldPosition,
		VoronoiCell *newPosition) {
	double energyBilance = 0;

	// loss of contact
	for (int iii = 0; iii < oldPosition->countNeighborCells; iii++) {
		if (oldPosition->neighborCells[iii]->getState() != FREE) {
			if (GetAgent( oldPosition->neighborCells[iii])
					== GetAgent( oldPosition))
				energyBilance -= INTRA_AGENT_ENERGY;
			else
				energyBilance -= INTER_AGENT_ENERGY;
		}
	}

	// gain new contact
	for (int iii = 0; iii < newPosition->countNeighborCells; iii++) {
		if (newPosition->neighborCells[iii]->getState() != FREE
				&& newPosition->neighborCells[iii] != oldPosition) {
			if (GetAgent( newPosition->neighborCells[iii])
					== GetAgent( oldPosition))
				energyBilance += INTRA_AGENT_ENERGY;
			else
				energyBilance += INTER_AGENT_ENERGY;
		}
	}

	return energyBilance;
}
//////////////////////////////////////

double getEnergyRemove(VoronoiCell *oldPosition) {
	double energy = 0;

	for (int iii = 0; iii < oldPosition->countNeighborCells; iii++) {
		if (oldPosition->neighborCells[iii]->getState() != FREE) {
			if (GetAgent( oldPosition->neighborCells[iii])
					== GetAgent( oldPosition)) {
				// cell deformation
				energy -= INTRA_AGENT_ENERGY;
			} else {
				// loss of cell-cell-contact
				if (GetAgent( oldPosition)->getGlucose()
						* GetAgent( oldPosition)->getOxygen()
						> THRESHOLD_NECROSIS_GLUCOSE_OXYGEN)
					energy -= INTER_AGENT_ENERGY;
			}
		}
	}

	return energy;
}
//////////////////////////////////////

double getEnergyPlaced(VoronoiCell *newPosition, VoronoiCell *oldPosition) {
	double energyBilance = 0;

	//fprintf( stderr, "getEnergyPlaced()\n");

	for (int iii = 0; iii < newPosition->countNeighborCells; iii++) {
		if (newPosition->neighborCells[iii]->getState() != FREE
				&& newPosition->neighborCells[iii] != oldPosition) {
			if (GetAgent( newPosition->neighborCells[iii])
					== GetAgent( oldPosition)) {
				// cell relaxation
				energyBilance += INTRA_AGENT_ENERGY;
			} else {
				// gain new contact
				//if( GetAgent( oldPosition)->getGlucose() * GetAgent( oldPosition)->getOxygen() > THRESHOLD_NECROSIS_GLUCOSE_OXYGEN)
				energyBilance += INTER_AGENT_ENERGY;
			}
		}
	}

	return energyBilance;
}
//////////////////////////////////////

double getEnergyDifference(VoronoiCell *oldPosition, VoronoiCell *newPosition) {
	double energyBilance = 0;

	// loss of contact
	for (int iii = 0; iii < oldPosition->countNeighborCells; iii++) {
		if (oldPosition->neighborCells[iii]->getState() != FREE) {
			if (GetAgent( oldPosition->neighborCells[iii])
					== GetAgent( oldPosition))
				energyBilance -= INTRA_AGENT_ENERGY;
			else
				energyBilance -= INTER_AGENT_ENERGY;
		}
	}

	// gain new contact
	for (int iii = 0; iii < newPosition->countNeighborCells; iii++) {
		if (newPosition->neighborCells[iii]->getState() != FREE
				&& newPosition->neighborCells[iii] != oldPosition) {
			if (GetAgent( newPosition->neighborCells[iii])
					== GetAgent( oldPosition))
				energyBilance += INTRA_AGENT_ENERGY;
			else
				energyBilance += INTER_AGENT_ENERGY;
		}
	}

	return energyBilance;
}

double getConnections(VoronoiCell *start, VoronoiCell *ignore, Agent *agent) {
	int ii;
	int stackLength = 1;
	int stackPosition = 0;
	VoronoiCell *actualLocation; // = agentArray->agents[i]->location[0];
	VoronoiCell *locationStack[100];
	locationStack[0] = start;

	//fprintf( stderr, "getConnections()\n");

	do {
		actualLocation = locationStack[stackPosition];
		stackPosition++;
		//int iii;
		//fprintf( stderr, "getConnections(): stack position %i\n", stackPosition);
		for (ii = 0; ii < actualLocation->countNeighborCells; ii++) {
			// neighbor is occupied by same neighbor?
			if (actualLocation->neighborCells[ii]->getState() != FREE
					&& actualLocation->neighborCells[ii] != ignore
					&& agent == GetAgent( actualLocation->neighborCells[ii])) {
				// is already in stack?
				int iv;
				for (iv = 0;
						iv < stackLength
								&& locationStack[iv]
										!= actualLocation->neighborCells[ii];
						iv++)
					;
				if (iv == stackLength) {
					locationStack[stackLength++] =
							actualLocation->neighborCells[ii];
				}
			}
		}
	} while (stackLength != stackPosition);

	//fprintf( stderr, "getConnections() end\n");
	return stackLength;
}

void performFreeMigration(VoronoiDiagram *voronoiDiagram, AgentList *agentArray,
		ActionTree *actionList, Action *selected_action,
		Statistix *myStatistics, double &Time, double &EndTime) {
	VoronoiCell *centralCell = GetVoronoiCell( selected_action->originalCell);

	if (selected_action->originalCell->state == COMPARTMENT) {
		int countFreeNeighbors =
				selected_action->originalCell->countFreeNeighbors();
		/*for( int i=0; i<centralCell->countNeighborCells; i++)
		 if(centralCell->neighborCells[i]->isFree())
		 countFreeNeighbors++;*/

		if (countFreeNeighbors == 0) {
			fprintf(
					stderr,
					"ERROR in performFreeMigration(): countFreeNeighbors = %i (%i)\n",
					countFreeNeighbors,
					selected_action->originalCell->location[0]->countFreeNeighborCells);
			fprintf(stderr, "cells: %i/%i -> g:%i, d:%i, n:%i\n",
					GetAgent(centralCell)->cellCount,
					GetAgent(centralCell)->maxCellCount,
					GetAgent(centralCell)->growingTumorCellCount,
					GetAgent(centralCell)->dividingTumorCellCount,
					GetAgent(centralCell)->necroticCellCount);
			for (int v = 0; v < centralCell->countNeighborCells; v++)
				fprintf(
						stderr,
						"%i. neighbor(%i): %i/%i -> g:%i, d:%i, n:%i\n",
						v + 1,
						GetAgent(centralCell->neighborCells[v])->index,
						GetAgent(centralCell->neighborCells[v])->cellCount,
						GetAgent(centralCell->neighborCells[v])->maxCellCount,
						GetAgent(centralCell->neighborCells[v])->growingTumorCellCount,
						GetAgent(centralCell->neighborCells[v])->dividingTumorCellCount,
						GetAgent(centralCell->neighborCells[v])->necroticCellCount);

			exit(0);
		}

		// randomly chose neighbor
		int i = (int) myRandE((double) countFreeNeighbors);

		for (int ni = 0; ni < centralCell->countNeighborCells; ni++) {
			if (centralCell->neighborCells[ni]->isFree()) {
				if (i == 0) {
					//fprintf( stderr, "FOUND\n");
					VoronoiCell *oldLocation = centralCell, *newLocation =
							centralCell->neighborCells[ni];
					if (newLocation->agent == NULL) {
						agentArray->activateAgent()->attach(newLocation);
						GetAgent( newLocation)->state = COMPARTMENT;
					}
					//fprintf( stderr, "MIGRATION: %i -> %i\n", GetAgent(oldLocation)->index, GetAgent(newLocation)->index);

					int whichCell =
							(int) myRandE(
									(double) selected_action->originalCell->growingTumorCellCount);
					int sumCells = 0;
					int randM = -1;
					for (int m = 0; m <= M_gro && randM < 0; m++) {
						if (selected_action->originalCell->actions[INDEX_GROWTH]->internalStateM[m]
								> 0) {
							sumCells +=
									selected_action->originalCell->actions[INDEX_GROWTH]->internalStateM[m];
							if (sumCells > whichCell) {
								randM = m;
							}
						}
					}

					GetAgent( oldLocation)->cellCount--;
					GetAgent( oldLocation)->growingTumorCellCount--;
					update_surrounding_reduced_compartment(actionList,
							GetAgent( oldLocation), voronoiDiagram);
					GetAgent( newLocation)->cellCount++;
					GetAgent( newLocation)->growingTumorCellCount++;
					update_surrounding_expanded_compartment(actionList,
							newLocation, voronoiDiagram);
					GetAgent( oldLocation)->actualize(actionList);
					GetAgent( newLocation)->actualize(actionList);
					actionList->actualizeRate(
							GetAgent( oldLocation)->actions[INDEX_GROWTH],
							GetAgent( oldLocation)->actions[INDEX_GROWTH]->getActualRate());
					actionList->actualizeRate(
							GetAgent( newLocation)->actions[INDEX_GROWTH],
							GetAgent( newLocation)->actions[INDEX_GROWTH]->getActualRate());
					actionList->actualizeRate(
							GetAgent( oldLocation)->actions[INDEX_MIGRATION],
							GetAgent( oldLocation)->actions[INDEX_MIGRATION]->getActualRate());
					actionList->actualizeRate(
							GetAgent( newLocation)->actions[INDEX_MIGRATION],
							GetAgent( newLocation)->actions[INDEX_MIGRATION]->getActualRate());

					//actionList->actualizeRate( daughterCell->actions[INDEX_MIGRATION], daughterCell->actions[INDEX_MIGRATION]->getActualRate());

					GetAgent( oldLocation)->actions[INDEX_GROWTH]->internalStateM[randM]--;
					GetAgent( newLocation)->actions[INDEX_GROWTH]->internalStateM[randM]++;
				}
				i--;
			}
		}

		//fprintf( stderr, "ERROR in performFreeMigration(): still no implementation for cell type COMPARTMENT!\n");
		//exit( 0);
		return;
	}

	// randomly chose neighbor
	int i = (int) myRandE((double) centralCell->countFreeNeighborCells);

	for (int ni = 0; ni < centralCell->countNeighborCells; ni++) {
		if (centralCell->neighborCells[ni]->isFree()) {
			if (i == 0) {
				// diffuse
				//fprintf( stderr, "Migrate Single Cell\n");
				VoronoiCell *oldLocation = centralCell, *newLocation =
						centralCell->neighborCells[ni];

				selected_action->originalCell->detach(oldLocation);
				update_surrounding_reduced_cell(actionList, oldLocation,
						voronoiDiagram);

				// exchange FREE and migrating cell
				if (newLocation->agent != NULL) {
					Agent *freeAgent = GetAgent( newLocation);
					freeAgent->detach(newLocation);
					freeAgent->attach(oldLocation);
				}
				// place FREE agent and migrating cell
				else {
					Agent* oldAgent = agentArray->activateAgent();
					oldAgent->attach(oldLocation);
					oldAgent->state = FREE;
				}

				selected_action->originalCell->attach(newLocation);
				update_surrounding_expanded_cell(actionList, newLocation,
						voronoiDiagram);

				// exchange concentrations
				double temp_oxy = newLocation->oxygen;
				double temp_glu = newLocation->glucose;
				newLocation->oxygen = oldLocation->oxygen;
				newLocation->glucose = oldLocation->glucose;
				oldLocation->oxygen = temp_oxy;
				oldLocation->glucose = temp_glu;
			}
			i--;
		}
	}

	return;
}

void performMigration(VoronoiDiagram *voronoiDiagram, AgentList *agentArray,
		ActionTree *actionList, Action *selected_action,
		Statistix *myStatistics, double &Time, double &EndTime) {
	//fprintf( stderr, "MIGRATE\n");

	VoronoiCell *oldLocation = NULL, *newLocation = NULL;

	// ENERGY-MINIMIZING MIGRATION

	// count candidates
	int countFreeNeighborSites = 0;
	for (int i = 0; i < selected_action->originalCell->countLocations; i++)
		for (int ii = 0;
				ii
						< selected_action->originalCell->location[i]->countNeighborCells;
				ii++)
			if (selected_action->originalCell->location[i]->neighborCells[ii]->getState()
					== FREE)
				countFreeNeighborSites++;

	// calculate energies
	double energies[countFreeNeighborSites];
	double energySum = 0;
	countFreeNeighborSites = 0;
	double center[2] = {50,50};
	for (int i = 0; i < selected_action->originalCell->countLocations; i++)
		for (int ii = 0;
				ii
						< selected_action->originalCell->location[i]->countNeighborCells;
				ii++)
			if (selected_action->originalCell->location[i]->neighborCells[ii]->getState() == FREE &&
					selected_action->originalCell->location[i]->neighborCells[ii]->getDistanceTo(center)<50.) {
				// each contact to other cell
				double energyDiff =
						getEnergyDifference(
								selected_action->originalCell->location[i],
								selected_action->originalCell->location[i]->neighborCells[ii]);

				if(Agent::USE_GRAVITY){
					double dx = selected_action->originalCell->location[i]->neighborCells[ii]->position[0] - selected_action->originalCell->location[i]->position[0];
					double dy = selected_action->originalCell->location[i]->neighborCells[ii]->position[1] - selected_action->originalCell->location[i]->position[1];
#if DIMENSIONS == 3
					double dz = selected_action->originalCell->location[i]->neighborCells[ii]->position[2] - selected_action->originalCell->location[i]->position[2];
#else
					double dz=0;
#endif
					double vx=0, vy=-1, vz=0;

					double gravitation = exp( - acos(vy*dy / sqrt(vx*vx+vy*vy+vz*vz) / sqrt(dx*dx+dy*dy+dz*dz)) / alpha_ref);

					energies[countFreeNeighborSites] = gravitation;
				}else{
					// NOTOX
					//double gravitation = selected_action->originalCell->location[i]->neighborCells[ii]->position[1] - selected_action->originalCell->location[i]->position[1];
					double gravitation = 0;
					double temperature = 1e-2;
					energies[countFreeNeighborSites] =
							exp( (-gravitation*10 + energyDiff) / temperature);

				}

				if (getConnections(
						selected_action->originalCell->location[i]->neighborCells[ii],
						selected_action->originalCell->location[i],
						selected_action->originalCell)
						< getConnections(
								selected_action->originalCell->location[i],
								selected_action->originalCell->location[i]->neighborCells[ii],
								selected_action->originalCell))
					energies[countFreeNeighborSites] = 0;
				/*if(isinf( energies[countFreeNeighborSites])){
					printf("ERROR: infinit energy! Decreasing the temperature might help!\n");
					energies[countFreeNeighborSites] = (double)FLT_MAX;
					exit(0);
				}
				if(energies[countFreeNeighborSites]==0){
					printf("ERROR: gravitation=%e, energyDiff=%e, temperature=%e\n", gravitation, energyDiff, temperature);
					printf("ERROR: zero energy! Decreasing the temperature might help!\n");
					energies[countFreeNeighborSites] = 1./(double)FLT_MAX;
					exit(0);
				}*/
				energySum += energies[countFreeNeighborSites];
				//fprintf( stderr, "%i.candidate: i=%i, i+1=%i, energy = %lf, prob = %lf\n", countFreeNeighborSites+1, (int)energyBefore, (int)energyAfter, log( energies[countFreeNeighborSites]), energies[countFreeNeighborSites]);
				countFreeNeighborSites++;
			}

	// chose winner
	double randSum = myRand() * energySum;
	energySum = 0;
	char foundWinner = FALSE;
	countFreeNeighborSites = 0;
	for (int i = 0; i < selected_action->originalCell->countLocations; i++)
		for (int ii = 0;
				ii
						< selected_action->originalCell->location[i]->countNeighborCells;
				ii++)
			if (selected_action->originalCell->location[i]->neighborCells[ii]->getState() == FREE  &&
					selected_action->originalCell->location[i]->neighborCells[ii]->getDistanceTo(center)<50.) {
				// each contact to other cell
				energySum += energies[countFreeNeighborSites];
				if (!foundWinner && energySum > randSum) {
					oldLocation = selected_action->originalCell->location[i];
					newLocation =
							selected_action->originalCell->location[i]->neighborCells[ii];
					foundWinner = TRUE;
				}

				countFreeNeighborSites++;
			}

	// ENERGY-MINIMIZING MIGRATION INCLUDING CASE OF NO MIGRATION

	// Chose Moving Cell Cluster
	// Energies => Probabilities
	//fprintf( stderr, "Chose Moving Cell Cluster: Energies => Probabilities\n");
	/*double oldClusterProb[selected_action->originalCell->countLocations];
	 double probSum = 0.;
	 for( int i=0; i<selected_action->originalCell->countLocations; i++){
	 oldClusterProb[i] = exp( getEnergyRemove( selected_action->originalCell->location[i]));
	 probSum += oldClusterProb[i];
	 }

	 // Choice
	 //fprintf( stderr, "Chose Moving Cell Cluster: Choice\n");
	 double randProbSum = myRand() * probSum;
	 probSum = 0.;
	 int i=0;
	 do{
	 probSum += oldClusterProb[i++];
	 }while( probSum <= randProbSum);
	 oldLocation = selected_action->originalCell->location[i-1];
	 //	fprintf( stderr, "Chose Moving Cell Cluster: %i\n", );

	 // Chose Moving Destination
	 // Candidates
	 //fprintf( stderr, "Chose Moving Destination: Candidates\n");
	 int countFreeNeighborSites = 0;
	 for( i=0; i<selected_action->originalCell->countLocations; i++)
	 for( int ii=0; ii<selected_action->originalCell->location[i]->countNeighborCells; ii++)
	 if( selected_action->originalCell->location[i]->neighborCells[ii]->getState() == FREE)
	 countFreeNeighborSites++;

	 // Energies => Probabilities
	 //fprintf( stderr, "Chose Moving Destination: Energies => Probabilities\n");
	 double newClusterProb[countFreeNeighborSites];
	 probSum = 0.;
	 countFreeNeighborSites = 0;
	 for( i=0; i<selected_action->originalCell->countLocations; i++)
	 for( int ii=0; ii<selected_action->originalCell->location[i]->countNeighborCells; ii++)
	 if( selected_action->originalCell->location[i]->neighborCells[ii]->getState() == FREE){
	 //fprintf( stderr, "Chose Moving Destination: found free \n");

	 if( getConnections( selected_action->originalCell->location[i]->neighborCells[ii], oldLocation, selected_action->originalCell)
	 < getConnections( oldLocation, selected_action->originalCell->location[i]->neighborCells[ii], selected_action->originalCell))
	 newClusterProb[countFreeNeighborSites] = 0;
	 else
	 newClusterProb[countFreeNeighborSites] = exp( getEnergyPlaced( selected_action->originalCell->location[i]->neighborCells[ii], oldLocation));

	 //fprintf( stderr, "Chose Moving Destination: add prob \n");
	 probSum += newClusterProb[countFreeNeighborSites];
	 countFreeNeighborSites++;
	 }

	 // Case of no Migration
	 //fprintf( stderr, "Chose Moving Destination: Case of no Migration\n");
	 //double noChangeProb = exp( getEnergyPlaced( oldLocation, oldLocation));
	 //probSum += noChangeProb;
	 //double totalProbSum = probSum;
	 //double choiceProb = 0.;

	 // Choice
	 //fprintf( stderr, "Chose Moving Destination: Choice\n");
	 randProbSum = myRand() * probSum;
	 probSum = 0.;
	 countFreeNeighborSites = 0;
	 for( i=0; i<selected_action->originalCell->countLocations && newLocation==NULL; i++)
	 for( int ii=0; ii<selected_action->originalCell->location[i]->countNeighborCells && newLocation==NULL; ii++)
	 if( selected_action->originalCell->location[i]->neighborCells[ii]->getState() == FREE){
	 probSum += newClusterProb[countFreeNeighborSites];
	 if(probSum > randProbSum){
	 newLocation = selected_action->originalCell->location[i]->neighborCells[ii];
	 //choiceProb = newClusterProb[countFreeNeighborSites];
	 }
	 countFreeNeighborSites++;
	 }

	 if( newLocation==NULL){
	 //fprintf( stderr, "Don't Migrate: %lf\n", noChangeProb/totalProbSum);
	 return;
	 }*/
	//fprintf( stderr, "Migrate: %lf\n", choiceProb/totalProbSum);
	// Perform Migration
	if (newLocation != NULL) {

		selected_action->originalCell->detach(oldLocation);
		update_surrounding_reduced_cell(actionList, oldLocation,
				voronoiDiagram);

		// exchange FREE and migrating cell
		if (newLocation->agent != NULL) {
			Agent *freeAgent = GetAgent( newLocation);
			freeAgent->detach(newLocation);
			freeAgent->attach(oldLocation);
		}
		// place FREE agent and migrating cell
		else {
			Agent* oldAgent = agentArray->activateAgent();
			oldAgent->attach(oldLocation);
			oldAgent->state = FREE;
		}

		selected_action->originalCell->attach(newLocation);
		update_surrounding_expanded_cell(actionList, newLocation,
				voronoiDiagram);

		// exchange concentrations
		double temp_oxy = newLocation->oxygen;
		double temp_glu = newLocation->glucose;
		newLocation->oxygen = oldLocation->oxygen;
		newLocation->glucose = oldLocation->glucose;
		oldLocation->oxygen = temp_oxy;
		oldLocation->glucose = temp_glu;
	}

	//fprintf( stderr, "END MIGRATE\n");
}
/****************************************************************************/

void performChemotacticMigration(VoronoiDiagram *voronoiDiagram,
		AgentList *agentArray, ActionTree *actionList, Action *selected_action,
		Statistix *myStatistics, double &Time, double &EndTime) {
	//fprintf( stderr, "MIGRATE\n");

	VoronoiCell *oldLocation = NULL, *newLocation = NULL;

	// ENERGY-MINIMIZING MIGRATION

	//LUNGSYS
	//if (selected_action->originalCell->location[0]->morphogen < 1e-9)
	//	return;

	// count candidates
	int countFreeNeighborSites = 0;
	for (int i = 0; i < selected_action->originalCell->countLocations; i++)
		for (int ii = 0; ii < selected_action->originalCell->location[i]->countNeighborCells; ii++)
			if (selected_action->originalCell->location[i]->neighborCells[ii]->getState()== FREE)
				countFreeNeighborSites++;

	// calculate energies
	double energies[countFreeNeighborSites];
	double energySum = 0;
	countFreeNeighborSites = 0;
	for (int i = 0; i < selected_action->originalCell->countLocations; i++)
		for (int ii = 0; ii	< selected_action->originalCell->location[i]->countNeighborCells; ii++)
			if (selected_action->originalCell->location[i]->neighborCells[ii]->getState()== FREE) {
				// each contact to other cell
				/*double concDiff = selected_action->originalCell->location[i]->morphogen
				 - selected_action->originalCell->location[i]->neighborCells[ii]->morphogen;

				 energies[countFreeNeighborSites] = exp( -concDiff);
				 */

				//float beta = 1e6; // LUNGSYS: beta = temerature^-1
				double temperature = 1e-5; // NOTOX

				//energies[countFreeNeighborSites] = exp( beta*selected_action->originalCell->location[i]->neighborCells[ii]->morphogen);
				//energies[countFreeNeighborSites] = exp( -selected_action->originalCell->location[i]->neighborCells[ii]->morphogen);
				energies[countFreeNeighborSites] =
						exp( ((double)selected_action->originalCell->location[i]->neighborCells[ii]->morphogen
							- (double)selected_action->originalCell->location[i]->morphogen) / temperature);
				if(isinf( energies[countFreeNeighborSites])){
					printf("ERROR: infinit energy! Decreasing the temperature might help!\n");
					energies[countFreeNeighborSites] = (double)FLT_MAX;
				}
				if(energies[countFreeNeighborSites]==0){
					printf("ERROR: zero energy! Decreasing the temperature might help!\n");
					energies[countFreeNeighborSites] = 1./(double)FLT_MAX;
				}

				/*if(isinf( energies[countFreeNeighborSites])){
					printf("ERROR: %e, %e => %e => %e ==> %e\n",
							selected_action->originalCell->location[i]->neighborCells[ii]->morphogen,
							selected_action->originalCell->location[i]->morphogen,
							selected_action->originalCell->location[i]->neighborCells[ii]->morphogen
							- selected_action->originalCell->location[i]->morphogen,
							(selected_action->originalCell->location[i]->neighborCells[ii]->morphogen
							- selected_action->originalCell->location[i]->morphogen) / temperature,
							exp(  (selected_action->originalCell->location[i]->neighborCells[ii]->morphogen
								- selected_action->originalCell->location[i]->morphogen) / temperature));
					exit(0);
				}*/
				//energies[countFreeNeighborSites] = exp( selected_action->originalCell->location[i]->morphogen
				//		                              - selected_action->originalCell->location[i]->neighborCells[ii]->morphogen);

				energySum += energies[countFreeNeighborSites];
				//fprintf( stderr, "%i.candidate: i=%i, i+1=%i, energy = %lf, prob = %lf\n", countFreeNeighborSites+1, (int)energyBefore, (int)energyAfter, log( energies[countFreeNeighborSites]), energies[countFreeNeighborSites]);
				countFreeNeighborSites++;
			}

	//LUNGSYS
	//if (fabs(energySum) < 1e-9)
	//	return;

	//fprintf( stderr, "energySum=%lf\n", energySum);
	// chose winner
	double randSum = myRand() * energySum;
	energySum = 0;
	bool foundWinner = false;
	countFreeNeighborSites = 0;
	for (int i = 0; i < selected_action->originalCell->countLocations; i++)
		for (int ii = 0; ii < selected_action->originalCell->location[i]->countNeighborCells; ii++)
			if (selected_action->originalCell->location[i]->neighborCells[ii]->getState() == FREE) {
				// each contact to other cell
				energySum += energies[countFreeNeighborSites];
				if (!foundWinner && energySum >= randSum) {
					oldLocation = selected_action->originalCell->location[i];
					newLocation =
							selected_action->originalCell->location[i]->neighborCells[ii];
					foundWinner = true;
					//if( oldLocation->morphogen > newLocation->morphogen )
					/*{
						fprintf(
								stderr,
								"ERROR: migration against the morphogen gradient: %e -> %e =======> PROBABILITY = %e\n",
								oldLocation->morphogen, newLocation->morphogen,
								energies[countFreeNeighborSites] / energySum);
						//	exit(0);
					}*/
					//fprintf( stderr, "%i.candidate: energy = %lf, prob = %lf\n", countFreeNeighborSites+1, log( energies[countFreeNeighborSites]), energies[countFreeNeighborSites]);
					//exit(0);
				}

				countFreeNeighborSites++;
			}

	if( !foundWinner && countFreeNeighborSites){
		printf("ERROR: energySum=%e, randSum=%e, countFreeNeighborSites=%i\n",energySum,randSum, countFreeNeighborSites);
		exit(0);
	}
	// ENERGY-MINIMIZING MIGRATION INCLUDING CASE OF NO MIGRATION

	// Chose Moving Cell Cluster
	// Energies => Probabilities
	//fprintf( stderr, "Chose Moving Cell Cluster: Energies => Probabilities\n");
	/*double oldClusterProb[selected_action->originalCell->countLocations];
	 double probSum = 0.;
	 for( int i=0; i<selected_action->originalCell->countLocations; i++){
	 oldClusterProb[i] = exp( getEnergyRemove( selected_action->originalCell->location[i]));
	 probSum += oldClusterProb[i];
	 }

	 // Choice
	 //fprintf( stderr, "Chose Moving Cell Cluster: Choice\n");
	 double randProbSum = myRand() * probSum;
	 probSum = 0.;
	 int i=0;
	 do{
	 probSum += oldClusterProb[i++];
	 }while( probSum <= randProbSum);
	 oldLocation = selected_action->originalCell->location[i-1];
	 //	fprintf( stderr, "Chose Moving Cell Cluster: %i\n", );

	 // Chose Moving Destination
	 // Candidates
	 //fprintf( stderr, "Chose Moving Destination: Candidates\n");
	 int countFreeNeighborSites = 0;
	 for( i=0; i<selected_action->originalCell->countLocations; i++)
	 for( int ii=0; ii<selected_action->originalCell->location[i]->countNeighborCells; ii++)
	 if( selected_action->originalCell->location[i]->neighborCells[ii]->getState() == FREE)
	 countFreeNeighborSites++;

	 // Energies => Probabilities
	 //fprintf( stderr, "Chose Moving Destination: Energies => Probabilities\n");
	 double newClusterProb[countFreeNeighborSites];
	 probSum = 0.;
	 countFreeNeighborSites = 0;
	 for( i=0; i<selected_action->originalCell->countLocations; i++)
	 for( int ii=0; ii<selected_action->originalCell->location[i]->countNeighborCells; ii++)
	 if( selected_action->originalCell->location[i]->neighborCells[ii]->getState() == FREE){
	 //fprintf( stderr, "Chose Moving Destination: found free \n");

	 if( getConnections( selected_action->originalCell->location[i]->neighborCells[ii], oldLocation, selected_action->originalCell)
	 < getConnections( oldLocation, selected_action->originalCell->location[i]->neighborCells[ii], selected_action->originalCell))
	 newClusterProb[countFreeNeighborSites] = 0;
	 else
	 newClusterProb[countFreeNeighborSites] = exp( getEnergyPlaced( selected_action->originalCell->location[i]->neighborCells[ii], oldLocation));

	 //fprintf( stderr, "Chose Moving Destination: add prob \n");
	 probSum += newClusterProb[countFreeNeighborSites];
	 countFreeNeighborSites++;
	 }

	 // Case of no Migration
	 //fprintf( stderr, "Chose Moving Destination: Case of no Migration\n");
	 //double noChangeProb = exp( getEnergyPlaced( oldLocation, oldLocation));
	 //probSum += noChangeProb;
	 //double totalProbSum = probSum;
	 //double choiceProb = 0.;

	 // Choice
	 //fprintf( stderr, "Chose Moving Destination: Choice\n");
	 randProbSum = myRand() * probSum;
	 probSum = 0.;
	 countFreeNeighborSites = 0;
	 for( i=0; i<selected_action->originalCell->countLocations && newLocation==NULL; i++)
	 for( int ii=0; ii<selected_action->originalCell->location[i]->countNeighborCells && newLocation==NULL; ii++)
	 if( selected_action->originalCell->location[i]->neighborCells[ii]->getState() == FREE){
	 probSum += newClusterProb[countFreeNeighborSites];
	 if(probSum > randProbSum){
	 newLocation = selected_action->originalCell->location[i]->neighborCells[ii];
	 //choiceProb = newClusterProb[countFreeNeighborSites];
	 }
	 countFreeNeighborSites++;
	 }

	 if( newLocation==NULL){
	 //fprintf( stderr, "Don't Migrate: %lf\n", noChangeProb/totalProbSum);
	 return;
	 }*/
	//fprintf( stderr, "Migrate: %lf\n", choiceProb/totalProbSum);
	// Perform Migration
	if (newLocation != NULL) {

		selected_action->originalCell->detach(oldLocation);
		update_surrounding_reduced_cell(actionList, oldLocation,
				voronoiDiagram);

		// exchange FREE and migrating cell
		if (newLocation->agent != NULL) {
			Agent *freeAgent = GetAgent( newLocation);
			freeAgent->detach(newLocation);
			freeAgent->attach(oldLocation);
		}
		// place FREE agent and migrating cell
		else {
			Agent* oldAgent = agentArray->activateAgent();
			oldAgent->attach(oldLocation);
			oldAgent->state = FREE;
		}

		selected_action->originalCell->attach(newLocation);
		update_surrounding_expanded_cell(actionList, newLocation,
				voronoiDiagram);

		// exchange concentrations
		double temp_oxy = newLocation->oxygen;
		double temp_glu = newLocation->glucose;
		newLocation->oxygen = oldLocation->oxygen;
		newLocation->glucose = oldLocation->glucose;
		oldLocation->oxygen = temp_oxy;
		oldLocation->glucose = temp_glu;
	}

	//fprintf( stderr, "END MIGRATE\n");
}
/****************************************************************************/

VoronoiCell* performGrowthAndDivision(VoronoiDiagram *voronoiDiagram,
		AgentList *agentArray, ActionTree *actionList, Action *selected_action,
		Statistix *myStatistics, double &Time, double &EndTime) {
	if (selected_action->originalCell->state == COMPARTMENT) {
		Agent *daughterCell = NULL;
		//printf("%i %lf\n", selected_action->originalCell->growingTumorCellCount, selected_action->rate);

		// VOLUME EXPANSION
		if (selected_action->originalCell->cellCount
				< selected_action->originalCell->maxCellCount) {
			// into same compartment	
			//fprintf( stderr, "GROW cell %i (%i=%ig+%id/%i)\n", selected_action->originalCell->index, selected_action->originalCell->cellCount, selected_action->originalCell->growingTumorCellCount, selected_action->originalCell->dividingTumorCellCount, selected_action->originalCell->maxCellCount);	
			daughterCell = selected_action->originalCell;

			daughterCell->cellCount++;
			daughterCell->growingTumorCellCount++;
			daughterCell->actualize(actionList);
		} else {
			int i;

			// count free neighbors
			int countFreeNeighbors = 0;
			for (i = 0;
					i
							< GetVoronoiCell(selected_action->originalCell)->countNeighborCells;
					i++) {
				if (GetAgent(GetVoronoiCell(selected_action->originalCell)->neighborCells[i])
						== NULL) {
					Agent *newAgent = agentArray->activateAgent();
					newAgent->state = COMPARTMENT;
					newAgent->attach(
							GetVoronoiCell(selected_action->originalCell)->neighborCells[i]);
					//fprintf( stderr, "ACTIVATE cell %i\n", newAgent->index);
				}

				if (GetAgent(GetVoronoiCell(selected_action->originalCell)->neighborCells[i])->cellCount
						< GetAgent(GetVoronoiCell(selected_action->originalCell)->neighborCells[i])->maxCellCount)
					countFreeNeighbors++;
			}
			if (countFreeNeighbors == 0 /*&& countSlightlyFreeNeighbors == 0*/) {
				fprintf(
						stderr,
						"ERROR in performGrowthAndDivision(): countFreeNeighbors = %i\n",
						countFreeNeighbors);
				//exit(0);
				return NULL;
			}

			// chose free neighbor
			//double rand = myRandIE( 0., (double)countFreeNeighbors);
			//int whichFreeNeighbor = (int)(myRand() * (double)countFreeNeighbors);
			int whichFreeNeighbor = (int) (myRandE((double) countFreeNeighbors));
			/*if( whichFreeNeighbor >= countFreeNeighbors){
			 fprintf( stderr, "ERROR in performGrowthAndDivision(): rand(%lf) * countFreeNeighbors(%i) =: whichFreeNeighbor(%i) >= countFreeNeighbors(%i)\n", rand, countFreeNeighbors, whichFreeNeighbor, countFreeNeighbors);
			 exit( 0);
			 }*/
			//fprintf( stderr, "which neighbor %i/%i\n", whichFreeNeighbor, countFreeNeighbors);
			countFreeNeighbors = 0;
			//for( i=0; i<GetVoronoiCell(selected_action->originalCell)->countNeighborCells && countFreeNeighbors!=whichFreeNeighbor; i++)
			for (i = 0;
					i
							< GetVoronoiCell(selected_action->originalCell)->countNeighborCells
							&& daughterCell == NULL; i++) {
				if (GetAgent(GetVoronoiCell(selected_action->originalCell)->neighborCells[i])->cellCount
						< GetAgent(GetVoronoiCell(selected_action->originalCell)->neighborCells[i])->maxCellCount) {
					if (countFreeNeighbors == whichFreeNeighbor)
						daughterCell =
								GetAgent(GetVoronoiCell(selected_action->originalCell)->neighborCells[i]);
					countFreeNeighbors++;
				}
			}
			//selected_action->originalCell->cellCount--;
			//selected_action->originalCell->growingTumorCellCount--;

			//fprintf( stderr, "GROW cell %i (%i/%i) to %i (%i/%i)\n", 
			//	selected_action->originalCell->index, selected_action->originalCell->cellCount, selected_action->originalCell->maxCellCount, 
			//	daughterCell->index, daughterCell->cellCount, daughterCell->maxCellCount);	

			if (daughterCell == NULL) {
				fprintf(
						stderr,
						"ERROR in performGrowthAndDivision(): daughterCell==NULL\n");
				fprintf(stderr,
						"whichFreeNeighbor(%i) >= countFreeNeighbors(%i) == \n",
						whichFreeNeighbor, countFreeNeighbors);
				exit(0);
			}

			daughterCell->cellCount++;
			daughterCell->growingTumorCellCount++;
			daughterCell->actualize(actionList);

			selected_action->originalCell->actualize(actionList);
		}
		//fprintf( stderr, "DIVISION: %i -> %i\n", selected_action->originalCell->index, daughterCell->index);

		myStatistics->count_expanded_cells[0]++;
		myStatistics->count_divided_cells[0]++;
		myStatistics->count_cell_volume[0]++;
		myStatistics->count_cells[0]++;
		//actionList->actualizeAllRates( voronoiDiagram);
		//actionList->actualizeRate( selected_action, selected_action->getActualRate());
		actionList->actualizeRate(
				selected_action->originalCell->actions[INDEX_GROWTH],
				selected_action->originalCell->actions[INDEX_GROWTH]->getActualRate());
		actionList->actualizeRate(daughterCell->actions[INDEX_GROWTH],
				daughterCell->actions[INDEX_GROWTH]->getActualRate());
#if USE_MIGRATION
		actionList->actualizeRate(daughterCell->actions[INDEX_MIGRATION],
				daughterCell->actions[INDEX_MIGRATION]->getActualRate());
#endif
		update_surrounding_expanded_compartment(actionList,
				GetVoronoiCell(daughterCell), voronoiDiagram);

		//actionList->actualizeRate( daughterCell->actions[INDEX_DIVISION], daughterCell->actions[INDEX_DIVISION]->getActualRate());

		/*for(int i=0; i<GetVoronoiCell(selected_action->originalCell)->countNeighborCells; i++)
		 if(GetVoronoiCell(selected_action->originalCell)->neighborCells[i]->agent!=NULL)
		 GetAgent( GetVoronoiCell(selected_action->originalCell)->neighborCells[i])->actualize( actionList);
		 for(int i=0; i<GetVoronoiCell(daughterCell)->countNeighborCells; i++)
		 if(GetVoronoiCell(daughterCell)->neighborCells[i]->agent!=NULL)
		 GetAgent( GetVoronoiCell(daughterCell)->neighborCells[i])->actualize( actionList);
		 */
		return GetVoronoiCell(daughterCell);
	} else {
		fprintf(
				stderr,
				"NOT YET IMPLIMENTED in performGrowthAndDivision(): case selected_action->originalCell->state != COMPARTMENT\n");
		exit(0);
	}
}

VoronoiCell* performGrowth(VoronoiDiagram *voronoiDiagram,
		AgentList *agentArray, ActionTree *actionList, Action *selected_action,
		Statistix *myStatistics, double &Time, double &EndTime) {
	if (selected_action->originalCell->state == COMPARTMENT) {
		Agent *daughterCell = NULL;
		//printf("%i %lf\n", selected_action->originalCell->growingTumorCellCount, selected_action->rate);
		if (selected_action->originalCell->cellCount
				< selected_action->originalCell->maxCellCount) {
			// into same compartment	
			//fprintf( stderr, "GROW cell %i (%i=%ig+%id/%i)\n", selected_action->originalCell->index, selected_action->originalCell->cellCount, selected_action->originalCell->growingTumorCellCount, selected_action->originalCell->dividingTumorCellCount, selected_action->originalCell->maxCellCount);	
			selected_action->originalCell->cellCount++;
			selected_action->originalCell->growingTumorCellCount--;
			selected_action->originalCell->dividingTumorCellCount++;
			selected_action->originalCell->actualize(actionList);
			daughterCell = selected_action->originalCell;
		} else {
			int i;

			// count free neighbors
			int countFreeNeighbors = 0;
			int countSlightlyFreeNeighbors = 0;
			for (i = 0;
					i
							< GetVoronoiCell(selected_action->originalCell)->countNeighborCells;
					i++) {
				if (GetAgent(GetVoronoiCell(selected_action->originalCell)->neighborCells[i])
						== NULL) {
					Agent *newAgent = agentArray->activateAgent();
					newAgent->state = COMPARTMENT;
					newAgent->attach(
							GetVoronoiCell(selected_action->originalCell)->neighborCells[i]);
					//fprintf( stderr, "ACTIVATE cell %i\n", newAgent->index);
				}

				if (GetAgent(GetVoronoiCell(selected_action->originalCell)->neighborCells[i])->cellCount
						+ 1
						< GetAgent(GetVoronoiCell(selected_action->originalCell)->neighborCells[i])->maxCellCount)
					countFreeNeighbors++;
				if (GetAgent(GetVoronoiCell(selected_action->originalCell)->neighborCells[i])->cellCount
						< GetAgent(GetVoronoiCell(selected_action->originalCell)->neighborCells[i])->maxCellCount)
					countSlightlyFreeNeighbors++;
			}
			if (countFreeNeighbors == 0 /*&& countSlightlyFreeNeighbors == 0*/)
				return NULL;
			if (countFreeNeighbors == 0) {
				int whichSlightlyFreeNeighbor = (int) (myRand()
						* (double) countSlightlyFreeNeighbors);
				countSlightlyFreeNeighbors = 0;
				//for( i=0; i<GetVoronoiCell(selected_action->originalCell)->countNeighborCells && countFreeNeighbors!=whichFreeNeighbor; i++)
				for (i = 0;
						i
								< GetVoronoiCell(selected_action->originalCell)->countNeighborCells
								&& daughterCell == NULL; i++) {
					if (GetAgent(GetVoronoiCell(selected_action->originalCell)->neighborCells[i])->cellCount
							< GetAgent(GetVoronoiCell(selected_action->originalCell)->neighborCells[i])->maxCellCount) {
						if (countSlightlyFreeNeighbors
								== whichSlightlyFreeNeighbor) {
							daughterCell =
									GetAgent(GetVoronoiCell(selected_action->originalCell)->neighborCells[i]);
							daughterCell->maxCellCount++;
						}
						countSlightlyFreeNeighbors++;
					}
				}

			} else {
				//if( countFreeNeighbors == 0){
				//	return NULL;
				//}

				// chose free neighbor
				int whichFreeNeighbor = (int) (myRand()
						* (double) countFreeNeighbors);
				//fprintf( stderr, "which neighbor %i/%i\n", whichFreeNeighbor, countFreeNeighbors);
				countFreeNeighbors = 0;
				//for( i=0; i<GetVoronoiCell(selected_action->originalCell)->countNeighborCells && countFreeNeighbors!=whichFreeNeighbor; i++)
				for (i = 0;
						i
								< GetVoronoiCell(selected_action->originalCell)->countNeighborCells
								&& daughterCell == NULL; i++) {
					if (GetAgent(GetVoronoiCell(selected_action->originalCell)->neighborCells[i])->cellCount
							+ 1
							< GetAgent(GetVoronoiCell(selected_action->originalCell)->neighborCells[i])->maxCellCount) {
						if (countFreeNeighbors == whichFreeNeighbor)
							daughterCell =
									GetAgent(GetVoronoiCell(selected_action->originalCell)->neighborCells[i]);
						countFreeNeighbors++;
					}
				}
			}
			selected_action->originalCell->cellCount--;
			selected_action->originalCell->growingTumorCellCount--;

			//fprintf( stderr, "GROW cell %i (%i/%i) to %i (%i/%i)\n", 
			//	selected_action->originalCell->index, selected_action->originalCell->cellCount, selected_action->originalCell->maxCellCount, 
			//	daughterCell->index, daughterCell->cellCount, daughterCell->maxCellCount);	

			daughterCell->cellCount += 2;
			daughterCell->dividingTumorCellCount++;

			selected_action->originalCell->actualize(actionList);
			daughterCell->actualize(actionList);
		}

		myStatistics->count_expanded_cells[0]++;
		myStatistics->count_cell_volume[0]++;
		//actionList->actualizeAllRates( voronoiDiagram);
		//actionList->actualizeRate( selected_action, selected_action->getActualRate());
		actionList->actualizeRate(
				selected_action->originalCell->actions[INDEX_GROWTH],
				selected_action->originalCell->actions[INDEX_GROWTH]->getActualRate());
		actionList->actualizeRate(daughterCell->actions[INDEX_DIVISION],
				daughterCell->actions[INDEX_DIVISION]->getActualRate());

		//update_surrounding_expanded_cell( actionList, GetVoronoiCell(selected_action->originalCell), voronoiDiagram);
		//update_surrounding_expanded_cell( actionList, GetVoronoiCell(daughterCell), voronoiDiagram);
		/*for(int i=0; i<GetVoronoiCell(selected_action->originalCell)->countNeighborCells; i++)
		 if(GetVoronoiCell(selected_action->originalCell)->neighborCells[i]->agent!=NULL)
		 GetAgent( GetVoronoiCell(selected_action->originalCell)->neighborCells[i])->actualize( actionList);
		 for(int i=0; i<GetVoronoiCell(daughterCell)->countNeighborCells; i++)
		 if(GetVoronoiCell(daughterCell)->neighborCells[i]->agent!=NULL)
		 GetAgent( GetVoronoiCell(daughterCell)->neighborCells[i])->actualize( actionList);
		 */
		return GetVoronoiCell(daughterCell);
	}

	int i;
	int shifted = FALSE;

	//fprintf(stderr,"GROWTH\n");
	VoronoiCell *newLocation = NULL;
	VoronoiCell *sourceLocation;
	VoronoiCell *occupiedLocation = growCell(voronoiDiagram,
			selected_action->originalCell, shifted, &sourceLocation);
	//fprintf(stderr,"CELL GROWTH (a:%i,c1:%i -> c2:%i)\n", selected_action->originalCell->index, selected_action->originalCell->location[0]->index, occupiedLocation->index);

	int shiftPathLength;
	VoronoiCell **_shiftPath = NULL;

	// SHIFT
	if (shifted) {
		/*fprintf(stderr,"SHIFT\n");
		 fprintf(stderr," from %i\n", sourceLocation->index);
		 fprintf(stderr," to %i\n", occupiedLocation->index);*/
		// get shift path
		//fprintf(stderr,"Get shift path!\n");
		//int shiftPathLength;
		//VoronoiCell **_shiftPath = getShiftPath( selected_action->originalCell, occupiedLocation, shiftPathLength);
		//fprintf(stderr,"occupiedLocation: %i!\n", occupiedLocation->index);
		//fprintf(stderr,"sourceLocation: %i (%p)!\n",   sourceLocation_->index,   sourceLocation_);
		_shiftPath =
				getShiftPath(selected_action->originalCell, &occupiedLocation,
						sourceLocation, shiftPathLength);

		//fprintf(stderr,"=> %i - > %i!\n", _shiftPath[0]->index, _shiftPath[shiftPathLength-1]->index);
		//fprintf(stderr,"=> %i - > %i!\n", _shiftPath[0]->index, newLocation->index);
		if (_shiftPath[shiftPathLength - 1]->agent != NULL) {
			if (GetAgent( _shiftPath[shiftPathLength-1])->countLocations == 1)
				agentArray->deactivateAgent(
						GetAgent( _shiftPath[shiftPathLength-1]));
			GetAgent( _shiftPath[shiftPathLength-1])->detach(
					_shiftPath[shiftPathLength - 1]);
			myStatistics->count_inner_empty_volume[0]--;
		}

		// shift
		//fprintf(stderr,"Shift along path!\n");
		shiftPath(_shiftPath, shiftPathLength);

		// actualize shifted agents
		//fprintf(stderr,"Actualize shifted agents!\n");
		//for( i=1; i<shiftPathLength; i++)
		//	GetAgent( _shiftPath[i])->actualize( actionList);

		newLocation = _shiftPath[0];
		//fprintf(stderr,"Done!\n");
		//exit( 0);
	} else
		newLocation = occupiedLocation;

#ifdef REFINE
	if( occupiedLocation->coarseParent==NULL) {
		fprintf( stderr, "ERROR in performDivision(): cell: %i has no coarse parent!\n", occupiedLocation->index);
		fprintf( stderr, "growing agent sits on cell %i\n", selected_action->originalCell->location[0]->index);
		printVoronoiDiagram( voronoiDiagram, "errorVD.eps", true);
		printTriangulation( voronoiDiagram, "errorT.eps", true);
	}
	if( GetAgent(occupiedLocation->coarseParent)->countFree == GetAgent(occupiedLocation->coarseParent)->maxCellCount) {
		refineNeighborhood( voronoiDiagram, occupiedLocation, actionList, agentArray, pow( GetAgent(occupiedLocation->coarseParent)->maxCellCount,1./DIMENSIONS ));
		refineSurrounding( voronoiDiagram, occupiedLocation, actionList, agentArray, pow( GetAgent(occupiedLocation->coarseParent)->maxCellCount,1./DIMENSIONS ) + 0.5);
		refineNeighborhood( voronoiDiagram, occupiedLocation, actionList, agentArray, pow( GetAgent(occupiedLocation->coarseParent)->maxCellCount,1./DIMENSIONS ) + 0.5);
	}
	GetAgent(occupiedLocation->coarseParent)->countActive++;
	GetAgent(occupiedLocation->coarseParent)->countFree--;
#endif

	/*if (newLocation != NULL)
	 printf("1:%i\n", newLocation->extendedNeighborCellsInitialized);
	 printf("2:%i\n", selected_action->originalCell->location[0]->extendedNeighborCellsInitialized);
	 */
	// For DYNAMIC_EXTENDED_NEIGHBORHOOD ignore non-growing cells
	if (newLocation == NULL
			&& (VoronoiCell::USE_DYNAMIC_EXTENDED_NEIGHBORHOOD
					|| VoronoiCell::USE_SYMBOLIC_EXTENDED_NEIGHBORHOOD)) {
		selected_action->originalCell->state = NONACTIVE;
		return NULL;
	}

	if (newLocation == NULL) {
		//|| ((newLocation->extendedNeighborCellsInitialized == FALSE
		//		|| selected_action->originalCell->location[0]->extendedNeighborCellsInitialized == FALSE) && CellDivisionDepth > 1.)) {
		fprintf(stderr, "Cell %i isn't able to divide!\n",
				selected_action->originalCell->index);
		//fprintf(stderr,"newLocation->extendedNeighborCellsInitialized == %i\n", newLocation->extendedNeighborCellsInitialized);
		fprintf(
				stderr,
				"Free Neighbors: %i/%i\n",
				selected_action->originalCell->location[0]->countFreeNeighborCells,
				selected_action->originalCell->location[0]->countNeighborCells);
		fprintf(
				stderr,
				"Free Extended Neighbors: %i/%i\n",
				selected_action->originalCell->location[0]->countFreeExtendedNeighborCells,
				selected_action->originalCell->location[0]->countExtendedNeighborCells);
		int countFreeNeighbors = 0;
		for (i = 0;
				i
						< selected_action->originalCell->location[0]->countExtendedNeighborCells;
				i++) {
			if (selected_action->originalCell->location[0]->extendedNeighborhood[i]->isFree())
				countFreeNeighbors++;
		}
		fprintf(stderr, "Found Free Extended Neighbors: %i\n",
				countFreeNeighbors);
		EndTime = (Time < EndTime ? Time : EndTime);
	} else {
#if _COMMENTS_ > 1
		if( newLocation->countExtendedNeighborCells<1)
		fprintf(stderr,"Cell %i is growing into non-initialized area: VoronoiCell %i!\n", selected_action->originalCell->index, newLocation->index);
#endif
		//fprintf(stderr,"CELL GROWTH\n");
		if (newLocation->agent != NULL) {
			if (GetAgent( newLocation)->countLocations == 1)
				agentArray->deactivateAgent(GetAgent( newLocation));
			GetAgent( newLocation)->detach(newLocation);
			myStatistics->count_inner_empty_volume[0]--;
		}
		selected_action->originalCell->attach(newLocation);
		//fprintf(stderr,"update_surrounding_expanded_cell( actionList, newLocation, voronoiDiagram);\n", selected_action->originalCell->index, newLocation->index);
		update_surrounding_expanded_cell(actionList, occupiedLocation,
				voronoiDiagram);
		//fprintf(stderr,"end update_surrounding_expanded_cell( actionList, newLocation, voronoiDiagram);\n", selected_action->originalCell->index, newLocation->index);
		//update_expanded_cell( actionList, newLocation, voronoiDiagram);
		selected_action->originalCell->actualize(actionList);
	}

	if (shifted) {

		// actualize shifted agents
		//fprintf(stderr,"Actualize shifted agents!\n");
		//fprintf( stderr, "[ %i (%i)", sourceLocation->index, GetAgent(sourceLocation)->index);
		for (i = 0; i < shiftPathLength; i++) {
			GetAgent( _shiftPath[i])->actualize(actionList);
			//	fprintf( stderr, ", %i (%i)", _shiftPath[i]->index, GetAgent(_shiftPath[i])->index);
		}
		//fprintf( stderr, "]\n");
		free(_shiftPath);
	}

	selected_action->internalState = 0;
	myStatistics->count_expanded_cells[0]++;
	myStatistics->count_cell_volume[0]++;

	selected_action->originalCell->growingTumorCellCount--;
	selected_action->originalCell->dividingTumorCellCount++;
	//daughterCell->cellCount++;

	return newLocation;

}

/****************************************************************************/

void performDivision(VoronoiDiagram *voronoiDiagram, AgentList *agentArray,
		ActionTree *actionList, Action *selected_action,
		Statistix *myStatistics, double &Time, double &EndTime) {
	if (selected_action->originalCell->state == COMPARTMENT) {

		//fprintf( stderr, "DIVIDE cell %i (%i/%i)\n", selected_action->originalCell->index, selected_action->originalCell->cellCount, selected_action->originalCell->maxCellCount);
		selected_action->originalCell->growingTumorCellCount += 2;
		selected_action->originalCell->dividingTumorCellCount--;
		selected_action->originalCell->actualize(actionList);

		myStatistics->count_divided_cells[0]++;
		myStatistics->count_cells[0]++;
		//actionList->actualizeAllRates( voronoiDiagram);
		//actionList->actualizeRate( selected_action, selected_action->getActualRate());
		actionList->actualizeRate(
				selected_action->originalCell->actions[INDEX_GROWTH],
				selected_action->originalCell->actions[INDEX_GROWTH]->getActualRate());
		actionList->actualizeRate(
				selected_action->originalCell->actions[INDEX_DIVISION],
				selected_action->originalCell->actions[INDEX_DIVISION]->getActualRate());

		return;
	}
	/*fprintf(stderr,"DIVISION (cell %i, state %i: cells %i of %i?)\n",
	 GetVoronoiCell( selected_action->originalCell)->index,
	 selected_action->originalCell->state,
	 selected_action->originalCell->cellCount,
	 selected_action->originalCell->maxCellCount);*/
	/*if( selected_action->originalCell->state == COMPARTMENT){
	 fprintf(stderr,"COMPARTMENT DIVISION\n");
	 VoronoiCell * newLocation = growCell( selected_action->originalCell);
	 Agent* expandedCell;
	 //if( newLocation==NULL)
	 //	expandedCell = selected_action->originalCell;
	 //else
	 expandedCell = GetAgent( newLocation);
	 if( expandedCell==NULL){
	 // init new compartment
	 expandedCell = agentArray->activateAgent();
	 expandedCell->attach( newLocation);
	 expandedCell->state = COMPARTMENT;
	 initCellActions( expandedCell);
	 }
	 selected_action->originalCell->dividingTumorCellCount --;
	 selected_action->originalCell->growingTumorCellCount ++;
	 //selected_action->originalCell->cellCount --;
	 expandedCell->growingTumorCellCount ++;
	 expandedCell->cellCount ++;
	 count_divided_cells++;
	 update_surrounding_expanded_cell( actionList, newLocation, voronoiDiagram);

	 expandedCell->actualize( actionList);
	 selected_action->originalCell->actualize( actionList);
	 if( selected_action->originalCell->isGrowing())
	 actionList->actualizeRate( selected_action->originalCell->actions[INDEX_GROWTH],   selected_action->originalCell->actions[INDEX_GROWTH]->getActualRate());
	 //actionList->getDepth( actionList->root);

	 if( selected_action->originalCell->isDividing())
	 actionList->actualizeRate( selected_action->originalCell->actions[INDEX_DIVISION], selected_action->originalCell->actions[INDEX_DIVISION]->getActualRate());
	 //actionList->getDepth( actionList->root);
	 //selected_action->originalCell->actions[INDEX_GROWTH]->actualizeRate();
	 //selected_action->originalCell->actions[INDEX_DIVISION]->actualizeRate();
	 }else
	 */

	/*#if SUBCELLULAR_COMPONENTS > 0
	 if( selected_action->originalCell->countLocations == 2*SUBCELLULAR_COMPONENTS){


	 #endif*/
	//sumDivTimes += Time - selected_action->originalCell->timeOfCreation;
	myStatistics->count_divided_cells[0]++;
	myStatistics->count_cells[0]++;
	//global_count_divided_cells++;
	//fprintf( stderr, "DIVISION\n");
#if SUBCELLULAR_COMPONENTS > 0
	{
		// CELL DIVISION

		//fprintf(stderr,"CELL DIVISION\n");
		// STILL WORK IN PROGRESS
		// chose sub-cellular components forming 2nd daughter cell
		int i = selected_action->originalCell->countLocations - 1; //(int)(selected_action->originalCell->countLocations*myRand());
		VoronoiCell * daughterCellLocation =
				selected_action->originalCell->location[i];

		// adapt 1st daughter cell
		selected_action->originalCell->detach(daughterCellLocation);
		//selected_action->originalCell->location[i] = selected_action->originalCell->location[--(selected_action->originalCell->countLocations)];

		// create 2nd daughter cell
		Agent* daughterCell = agentArray->activateAgent();
		daughterCell->attach(daughterCellLocation);
		daughterCell->state = ACTIVE;
		initCellActions(daughterCell);

//		GetAgent(daughterCellLocation->coarseParent)->countActive++;
//		GetAgent(daughterCellLocation->coarseParent)->countFree--;
		/*fprintf( stderr, "daughterCell: active:%i, nonact:%i, free:%i, max:%i\n",
		 GetAgent(daughterCellLocation->coarseParent)->countActive,
		 GetAgent(daughterCellLocation->coarseParent)->countNonactive,
		 GetAgent(daughterCellLocation->coarseParent)->countFree,
		 GetAgent(daughterCellLocation->coarseParent)->maxCellCount);*/

		// Get Distance to closest free Voronoi Cell
		int dist = 0;
		int min_dist = 0;
		for (int d = 0; d < DIMENSIONS; d++)
			min_dist += voronoiDiagram->xN[d];

		if (VoronoiCell::USE_SYMBOLIC_EXTENDED_NEIGHBORHOOD) {
			VoronoiCell *closestVC;
			if (VoronoiCell::SHIFT_TO_UNOCCUPIED)
				closestVC =
						voronoiDiagram->searchClosestUnoccupiedVoronoiCell(
								daughterCell->location[0],
								VoronoiCell::extendedNeighborDepth,
								VoronoiCell::symbolicExtendedNeighborhoodSize,
								VoronoiCell::symbolicExtendedNeighborhood);
			else
				closestVC =
						voronoiDiagram->searchClosestFreeVoronoiCell(
								daughterCell->location[0],
								VoronoiCell::extendedNeighborDepth,
								VoronoiCell::symbolicExtendedNeighborhoodSize,
								VoronoiCell::symbolicExtendedNeighborhood);
			if (closestVC)
				min_dist =
						closestVC->getDistanceSquareTo(
								daughterCell->location[0]);
		} else if (VoronoiCell::USE_DYNAMIC_EXTENDED_NEIGHBORHOOD) {
			VoronoiCell *closestVC;
			if (VoronoiCell::SHIFT_TO_UNOCCUPIED)
				closestVC =
						voronoiDiagram->searchClosestUnoccupiedVoronoiCell(
								daughterCell->location[0],
								VoronoiCell::extendedNeighborDepth);
			else
				closestVC =
						voronoiDiagram->searchClosestFreeVoronoiCell(
								daughterCell->location[0],
								VoronoiCell::extendedNeighborDepth);

			if (closestVC)
				min_dist =
						closestVC->getDistanceSquareTo(
								daughterCell->location[0]);
		} else {
			if (daughterCell->location[0]->countFreeNeighborCells)
				for (int i = 0;
						i < daughterCell->location[0]->countNeighborCells;
						i++) {
					dist = 0;
					for (int d = 0; d < DIMENSIONS; d++)
						dist +=
								pow(
										daughterCell->location[0]->position[d]
												- daughterCell->location[0]->neighborCells[i]->position[d],
										2);
					if (min_dist > dist)
						min_dist = dist;
				}

			if (daughterCell->location[0]->countFreeExtendedNeighborCells)
				for (int i = 0;
						i
								< daughterCell->location[0]->countExtendedNeighborCells;
						i++) {
					dist = 0;
					for (int d = 0; d < DIMENSIONS; d++)
						dist +=
								pow(
										daughterCell->location[0]->position[d]
												- daughterCell->location[0]->extendedNeighborhood[i]->position[d],
										2);
					if (min_dist > dist)
						min_dist = dist;
				}
		}

		double prob_reentering_cell_cycle =
				PROB_REENTERING_CELL_CYCLE( sqrt((double)min_dist) * SPATIAL_UNIT);

		// INTOXICATION -> QUIESCENCE
		if( selected_action->originalCell->intoxication > Agent::WASTE_INTOXICATED_CELL_CYCLES/MaxCellDivisionRate)
		//if( selected_action->originalCell->intoxication > Agent::WASTE_INTOXICATED_CELL_CYCLES)
			prob_reentering_cell_cycle=0;

		//if(agentArray->countActiveAgents >100 && selected_action->originalCell->location[0]->ecm <  )
		//	prob_reentering_cell_cycle = 0.;

		//double threshold = 0.003;

		if ( /*agentArray->countActiveAgents <100 ||*/
				( myRand() <= prob_reentering_cell_cycle
				&& VoronoiCell::ECM_THRESHOLD_QUIESCENCE <= daughterCell->location[0]->ecm))
			daughterCell->divide = 1.;
		else
			daughterCell->divide = 0.;

		if ( /*agentArray->countActiveAgents <100 ||*/
				( myRand() <= prob_reentering_cell_cycle
				&& VoronoiCell::ECM_THRESHOLD_QUIESCENCE <= selected_action->originalCell->location[0]->ecm))
			selected_action->originalCell->divide = 1.;
		else
			selected_action->originalCell->divide = 0.;

		//actionList->print();
		//fprintf( stderr, "daughterCell->actualize( actionList);\n ");
		//daughterCell->tumorCellCount++;
		daughterCell->actualize(actionList);
		//actionList->getDepth( actionList->root);
		//actionList->print();
		//fprintf( stderr, "selected_action->originalCell->actualize( actionList);\n ");
		selected_action->originalCell->actualize(actionList);

		selected_action->originalCell->generation++;
		daughterCell->generation = selected_action->originalCell->generation;
		daughterCell->waste = selected_action->originalCell->waste;
		daughterCell->intoxication = selected_action->originalCell->intoxication;

		selected_action->originalCell->growingTumorCellCount++;
		selected_action->originalCell->dividingTumorCellCount--;
		daughterCell->growingTumorCellCount++;
		daughterCell->cellCount++;

		//actionList->print();
		//actionList->getDepth( actionList->root);
		//fprintf( stderr, "after selected_action->originalCell->actualize( actionList);\n ");
		//fprintf(stderr,"CELL DIVISION (a:%i -> a:%i)\n", selected_action->originalCell->index, daughterCell->index);

		//update_surrounding_added_cell( actionList, daughterCell, voronoiDiagram);

		/*fprintf( stderr, "CELL DIVISION: ");
		 fprintf( stderr, "Agent %i ---[division to]---> agent %i { ",
		 selected_action->originalCell->index, selected_action->originalCell->index);
		 int ii;
		 for(ii=0; ii<selected_action->originalCell->countLocations; ii++){
		 fprintf( stderr, "%i |", selected_action->originalCell->location[ii]->index);
		 }
		 fprintf( stderr, "\b} && agent %i { ",
		 daughterCell->index);

		 for(ii=0; ii<daughterCell->countLocations; ii++){
		 fprintf( stderr, "%i |", daughterCell->location[ii]->index);
		 }
		 fprintf( stderr, "\b}\n");*/

		//flag_actualizeProcessRates = TRUE;
	}

#else // SUBCELLULAR_COMPONENTS == 0
	/*VoronoiCell * daughterCellLocation = growCell( selected_action->originalCell);

	 // create 2nd daughter cell
	 Agent* daughterCell = agentArray->activateAgent();
	 daughterCell->attach( daughterCellLocation);
	 update_surrounding_added_cell( actionList, daughterCell, voronoiDiagram);
	 daughterCell->tumorCellCount++;*/

	// I. CELL EXPANSION
	/*VoronoiCell * daughterCellLocation = growCell( selected_action->originalCell);
	 //selected_action->originalCell->attach( daughterCellLocation);
	 //update_surrounding_expanded_cell( actionList, daughterCellLocation, voronoiDiagram);
	 //selected_action->originalCell->actualize( actionList);

	 // II. CELL DIVISION
	 //int i = (int)(SUBCELLULAR_COMPONENTS*myRand());
	 //VoronoiCell * daughterCellLocation = selected_action->originalCell->location[i];

	 // adapt 1st daughter cell
	 //selected_action->originalCell->detach( daughterCellLocation);


	 // create 2nd daughter cell
	 Agent* daughterCell = agentArray->activateAgent();
	 daughterCell->attach( daughterCellLocation);
	 daughterCell->state = ACTIVE;
	 initCellActions( daughterCell);
	 update_surrounding_expanded_cell( actionList, daughterCellLocation, voronoiDiagram);
	 daughterCell->actualize( actionList);
	 selected_action->originalCell->actualize( actionList);

	 daughterCell->tumorCellCount++;
	 selected_action->internalState=0;
	 */

#endif // SUBCELLULAR_COMPONENTS > 0
	/*#if SUBCELLULAR_COMPONENTS > 0
	 }else{
	 // CELL EXPANSION
	 fprintf(stderr,"EXPANSION\n");
	 exit( 0);
	 //VoronoiCell * newLocation = growCell( selected_action->originalCell);
	 selected_action->originalCell->attach( newLocation);
	 update_surrounding_expanded_cell( actionList, newLocation, voronoiDiagram);
	 //update_expanded_cell( actionList, newLocation, voronoiDiagram);
	 selected_action->originalCell->actualize( actionList);
	 / *fprintf( stderr,"CELL EXPANSION: agent %i ---[expandes to]---> location %i ==> locations {",
	 selected_action->originalCell->index, newLocation->index);
	 for(i=0; i<selected_action->originalCell->countLocations; i++){
	 fprintf( stderr, "%i |", selected_action->originalCell->location[i]->index);
	 }
	 fprintf( stderr, "\b}\n");* /

	 //agentArray->printActiveAgents();
	 }
	 #endif*/

}


