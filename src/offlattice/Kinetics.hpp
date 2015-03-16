/*
 * Kinetics.hpp
 *
 *  Created on: 07.11.2013
 *      Author: jagiella
 */

#ifndef KINETICS_HPP_
#define KINETICS_HPP_

enum KineticTypes {Exact, OverDamped};
enum Kinetics {Reset = 1, BrownianMotion = 2, ExtendedHertz = 4, Chemotaxis = 8, Springs = 16};

float getTimeStep( AgentList<Agent> *al, int dimy, int kineticType);
void resetAcceleration( AgentList<Agent> *al, int dimy);

void addBrownianMotion( AgentList<Agent> *al, int dimy, float *randomVelocityCellTypes);
void addBrownianMotion( AgentList<Agent> *al, int dimy, float *randomVelocityCellTypes, float dt);
void addCellCellInteractionForces( AgentList<Agent> *al, int dimy, Boxes<Agent*> *box, float **adhesionEnergies);
void addCellMediumFriction( AgentList<Agent> *al, int dimy, float dt);
void addSpringForces( AgentList<Agent> *al, int dimy, Spring<int>** springs, int countSprings);
void addChemotacticForces( AgentList<Agent> *al);

void updatePressure( AgentList<Agent> *al, int dimy, Boxes<Agent*> *box, float **adhesionEnergies);
void updatePositions(  AgentList<> *al, int dimy, Boxes<Agent*> *box, float dt, int kineticType);

void performCellDivision( AgentList<Agent> *al, int dimy, Boxes<Agent*> *box, float dt, float minR, float maxR, float *cycleTime, bool WNTDependence = true);
void performCellDivisionAlongSubstrate( AgentList<Agent> *al, int dimy, Boxes<Agent*> *box, float dt, float minR, float maxR, float *cycleTime, bool WNTDependence = true);

void normrnd( float  *vector, int dimensions, float  mean, float  variance);
void normrnd( double *vector, int dimensions, double mean, double variance);

#endif /* KINETICS_HPP_ */
