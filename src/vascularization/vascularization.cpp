/*
 * vascularization.cpp
 *
 *  Created on: May 27, 2015
 *      Author: jagiella
 */

#include <stdio.h>
#include "../tinyxml.h"

//using namespace tinyxml;

int main( int argc, char **argv)
{
	// parse xml
	fprintf(stderr, "Read File %s\n", argv[1]);
	//XMLDocument doc;
	//doc.Parse( argv[1] );
	TiXmlDocument doc( argv[1]);
	doc.LoadFile();

	TiXmlHandle docHandle( &doc );

	//doc.Parse(argv[1]);

	fprintf(stderr, "Get Element\n");
	//XMLElement* vesselGraphElement = doc.FirstChildElement( "Model" )->FirstChildElement( "VesselGraph" );
	TiXmlElement *vesselGraphElement = doc.FirstChildElement( "Model" )->FirstChildElement( "VesselGraph");
	fprintf(stderr, "finished\n");

	float lc = atof(vesselGraphElement->Attribute("latticeConstant"));
	int x = atof(vesselGraphElement->Attribute("x")) / lc;
	int y = atof(vesselGraphElement->Attribute("y")) / lc;
	int z = atof(vesselGraphElement->Attribute("z")) / lc;
	fprintf(stderr, "%i x %i x %i elements\n", x,y,z);

	int N = x*y*z;
	float **nodes  = (float**) malloc(N*sizeof(float*));
	float  *inFlow = (float *) malloc(N*sizeof(float));
	for( int i=0; i<N; i++){
		nodes[i]  = (float*)malloc(3*sizeof(float));
		inFlow[i] = 0;
	}

	fprintf( stderr,  "Read Nodes\n");
	TiXmlElement *vesselNodeElement = vesselGraphElement->FirstChildElement( "Node");
	while( vesselNodeElement){
		int idx = atoi(vesselNodeElement->Attribute( "id"));
		nodes[idx][0] = atof(vesselNodeElement->Attribute( "x"));
		nodes[idx][1] = atof(vesselNodeElement->Attribute( "y"));
		nodes[idx][2] = atof(vesselNodeElement->Attribute( "z"));
		vesselNodeElement = vesselNodeElement->NextSiblingElement( "Node");
	}

	fprintf( stderr,  "Read Segments\n");
	TiXmlElement *vesselSegmentElement = vesselGraphElement->FirstChildElement( "Segment");
	//fprintf( stderr,  "Read Attribute\n");
	while( vesselSegmentElement){
		//fprintf( stderr,  "Read Attribute\n");
		float flow = atof(vesselSegmentElement->Attribute( "flow"));

		int targetNodeIdx;
		if( flow > 0)
			targetNodeIdx = atoi(vesselSegmentElement->Attribute( "node2"));
		else
			targetNodeIdx = atoi(vesselSegmentElement->Attribute( "node1"));
		//fprintf( stderr,  "%i %i\n", targetNodeIdx, N);
		inFlow[targetNodeIdx] += abs(flow);


		vesselSegmentElement = vesselSegmentElement->NextSiblingElement( "Segment");
	}

	float flowCoarsened[x/5][y/5];
	for( int i=0; i<x/5; i++)
		for( int j=0; j<y/5; j++)
			flowCoarsened[i][j] = 0;

	for( int i=0; i<N; i++)
		//if(nodes[i][2] == lc * z/2)
		{
			//fprintf( stdout,  "%f %f %f\n", nodes[i][0],nodes[i][1],inFlow[i]);
			flowCoarsened[(int)(nodes[i][0]/lc/5)][(int)(nodes[i][1]/lc/5)] += inFlow[i];
		}

	for( int i=0; i<x/5; i++)
			for( int j=0; j<y/5; j++)
				fprintf( stdout,  "%i %i %f\n", i,j,flowCoarsened[i][j] / 250.);
	/*for( int i=0; i<x; i++)
		for( int j=0; j<y; j++){
			int idx = i + j*x + z/2*x*y;
			if( inFlow[idx]!=0)
			fprintf( stdout,  "%i %i %f\n", i,j,inFlow[idx]);
		}*/

	return 0;
}
