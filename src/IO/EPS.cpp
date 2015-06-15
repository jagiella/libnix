/*
 * EpsIO.cpp
 *
 *  Created on: Dec 4, 2009
 *      Author: jagiella
 */

#include "EPS.hpp"
#include <fstream>


void EPS::PSwriteHeader(std::fstream *fs, double minx, double maxx, double miny, double maxy) {
	//printf( "\%%!PS-Adobe-3.0 EPSF-3.0\n");
	(*fs) << "%!PS-Adobe-3.0 EPSF-3.0\n";
	//(*fs) << "%%BoundingBox: "<<minx<<" "<<miny<<" "<<maxx<<" "<<maxy<<" "<<"\n";
	char temp[512];
	sprintf( temp, "%%%%BoundingBox: %i %i %i %i\n", (int)(minx+0.5), (int)(miny+0.5), (int)(maxx+0.5), (int)(maxy+0.5));
	(*fs) << temp;
	sprintf( temp, "%%%%HiResBoundingBox: %lf %lf %lf %lf\n", minx, miny, maxx, maxy);
	(*fs) << temp;
}

void EPS::PSwriteSphere(std::fstream *fs, double x, double y, double radius,
		char *color) {
	(*fs) << "0.0 setlinewidth\n ";
	(*fs) << color << "\n ";
	(*fs) << x << " "<< y << " " << radius << " 0 360 arc closepath\n";
	(*fs) << color << " fill\n";
}

void EPS::PSwriteLine(std::fstream *fs, double x1, double y1, double x2, double y2,
		double radius, char *color) {
	//printf( "0.01 setlinewidth\n %s setrgbcolor\n newpath %lf %lf moveto %lf %lf lineto stroke\n", color, x1, y1, x2, y2);
	(*fs) << "newpath\n";
	(*fs) << radius << " setlinewidth\n ";
	(*fs) << x1 << " " << y1 << " moveto\n ";
	(*fs) << x2 << " " << y2 << " lineto\n ";
	(*fs) << color << "\n";
	(*fs) << "stroke\n";
}


void EPS::PSfillPolygon( std::fstream *fs,
		double *x, double *y, int length,
		double radius, char *color) {

	(*fs) << "newpath\n";
	(*fs) << " "<<x[0]<<" "<<y[0]<<" moveto\n";

	for(int i=1; i<length; i++)
	(*fs) << " "<<x[i]<<" "<<y[i]<<" lineto\n";

	(*fs) << " closepath\n " << color << "\n fill\n";
}

void EPS::PSwriteCircle(std::fstream *fs, double x, double y, double radius,
		double lineWidth, char *color) {
	//printf( "0.01 setlinewidth\n %s setrgbcolor %lf %lf %lf 0 360 arc stroke\n", color, x, y, majorRadius);
	(*fs) << lineWidth << " setlinewidth\n ";
	(*fs) << color << "\n ";
	(*fs) << x << " "<< y << " " << radius << " 0 360 arc stroke\n";
	//torus {1.0,0.25 scale <1,1,1> rotate<0,0,0> translate<0,0,0>
}

void PSfillCircle(std::fstream *fs, double x, double y, double majorRadius,
		double minorRadius, char *color) {
	//printf( "0.01 setlinewidth\n %s setrgbcolor %lf %lf %lf 0 360 arc stroke\n", color, x, y, majorRadius);
	(*fs) << x << " "<< y << " " << majorRadius << " 0 360 arc closepath\n";
	(*fs) << color << " fill\n";
}

void EPS::PScolor( char color[], double R, double G, double B) {
	sprintf( color, "%.10lf %.10lf %.10lf setrgbcolor\n", R,G,B);
}

