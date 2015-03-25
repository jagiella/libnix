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
#include <unistd.h>

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
	// PREDEFINED PARAMETERS
	bool verbosity = false;
	bool snapshots = false;
	int num_snapshots = 1;
	float epsilon;
	char  directory[1024];
	double  parameters[100];
	int num_parameters=0;
	char data_type = 0; // 0=cell types, 1=foxa2

	// READ COMMANDLINE OPTIONS
	char c;
	while ((c = getopt (argc, argv, "vsS:e:d:F")) != -1)
		switch (c){
		case 'v':
			verbosity = true;
			break;
		case 's':
			snapshots = true;
			break;
		case 'S':
			num_snapshots = atoi(optarg);
			break;
		case 'e':
			epsilon = atof(optarg);
			break;
		case 'd':
			sprintf( directory, optarg);
			break;
		case 'F':
			data_type = 1; // use foxa2 data
			break;
		}
	for (int index = optind; index < argc; index++){
		parameters[num_parameters++] = atof( argv[index]);
	    //printf ("Non-option argument %s\n", argv[index]);
	}

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
     simCrypt( agents, box, &time, verbosity, snapshots, num_snapshots, parameters, num_parameters, data_type);

#endif
 }
