/*
 * main.cpp
 *
 *  Created on: Nov 26, 2012
 *      Author: jagiella
 */

#include "SimulationThread.hpp"
#include "Cells.hpp"
#include "Agent.hpp"

#ifndef NOGUI
#include "Window.hpp"
#include "WindowGL.hpp"

#include <QtWidgets/QApplication>
#include <QtWidgets/QWidget>
#include <QtGui/QPainter>
#endif

#include <stdio.h>
#include <locale.h>

#define EndOfLine 10

void readIntParameter( const char *text, int &parameter){
    size_t nbytes = 512;
    char *string = (char *) malloc (nbytes + 1);
    fprintf(stderr, "%s? [%i] : ", text, parameter);
    getline( &string, &nbytes, stdin);
    if( string[0]!=EndOfLine)
    	parameter = atoi(string);
    free(string);
}

void readFloatParameter( const char *text, float &parameter){
    size_t nbytes = 512;
    char *string = (char *) malloc (nbytes + 1);
    fprintf(stderr, "%s? [%f] : ", text, parameter);
    getline( &string, &nbytes, stdin);
    if( string[0]!=EndOfLine)
    	parameter = atof(string);
    free(string);
}


int main(int argc, char **argv)
 {
#ifndef NOGUI
	 setlocale(LC_ALL,"C");

     QApplication app(argc, argv);
     qsrand(0);
#endif

     // default parameters
     float minRadius = 5, maxRadius = 7; UNUSED(minRadius);
     float latticeConstant = maxRadius*2;

     // cell array
     //Cells *cells = 0;//new Cells(Nmax, minRadius, maxRadius);
     AgentList<>*agents = new AgentList<>();


     // boxes
     //Boxes<Agent*> *box = new Boxes<Agent*>( 400.0f, 1000.0f, latticeConstant);
     //Boxes<Agent*> *box = new Boxes<Agent*>( 2000.0f, 2000.0f, 20.0f, latticeConstant);
     Boxes<Agent*> *box = new Boxes<Agent*>( 400.0f, 400.0f, 400.0f, latticeConstant);

     float time;

     // visualisation window
#ifndef NOGUI
     WindowGL *wgl = new WindowGL(agents, box, &time);
     wgl->show();

     // simulation thread
     MyThread *t = new MyThread( agents, box, &time);
     t->start();

     return app.exec();
#else
     simCrypt( agents, box, &time);

#endif
 }
