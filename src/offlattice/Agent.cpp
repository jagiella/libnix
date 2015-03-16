/*
 * Agent.cpp
 *
 *  Created on: 12.06.2013
 *      Author: jagiella
 */

#include "Agent.hpp"

//static const float AgentList::pressureThreshold = 6e2;

// Notch Threshold (Priming Differentiation)
const float Agent::TPnotch = 1;
// WNT Threshold   (Priming Differentiation)
const float Agent::TPwnt = 0.8;
// WNT Threshold   (Terminal Differentiation)
const float Agent::TDwnt = 0.2;

float Agent::wntmin = 0;
float Agent::wntmax = 0;

char Agent::NotchKnockout = 0;
char Agent::WNTKnockout = 0;
char Agent::EphrinB12Knockout = 0;
char Agent::EphB2Knockout = 0;
char Agent::EphB3Knockout = 0;

