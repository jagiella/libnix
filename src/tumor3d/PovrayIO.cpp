/*
 * PovrayIO.cpp
 *
 *  Created on: 05.12.2009
 *      Author: jagiella
 */

#include "PovrayIO.h"

void PovrayIO::writePovrayHeader( std::fstream *fs, double minx, double miny, double maxx, double maxy)
{
	(*fs) << "#include \"colors.inc\"\n"
	      << "background { color White }\n"
	      << "camera {";
	(*fs) << "location <"<<(maxx-minx)/2.<<","<<(maxy-miny)/2.<<","<<(maxx-minx)/2.+(maxy-miny)/2.<<">\n";
	(*fs) << "look_at  <"<<(maxx-minx)/2.<<","<<(maxy-miny)/2.<<","<<0.<<">";
	(*fs) << "}\n";
	(*fs) << "light_source { <"<<(maxx-minx)/2.<<","<<(maxy-miny)/2.<<","<<(maxx-minx)/2.+(maxy-miny)/2.<<"> color White}\n";
	(*fs) << "global_settings { ambient_light rgb<0.5, 0.5, 0.5> }\n";
}


void PovrayIO::writeSphere( std::fstream *fs, double x, double y, double radius, const char *color)
{
	(*fs) << "sphere { <"<<x<<","<<y<<">,"<<radius<<" texture {pigment { color "<<color<<"}} finish{ambient 1}}\n";
}

void PovrayIO::writeSphere( std::fstream *fs, double x, double y, double z, double radius, const char *color)
{
	(*fs) << "sphere { <"<<x<<","<<y<<","<<z<<">,"<<radius<<" texture {pigment { color "<<color<<"}} finish{ambient 1}}\n";
}

void PovrayIO::writeCylinder( std::fstream *fs, double x1, double y1, double x2, double y2, double radius, const char *color)
{
	(*fs) << "cylinder{ <"<<x1<<","<<y1<<">,<"<<x2<<","<<y2<<">,"<<radius<<" pigment{color "<<color<<"}}\n";
}

void PovrayIO::writeCircle( std::fstream *fs, double x, double y, double majorRadius, double minorRadius, const char *color)
{
	(*fs) << "torus { "<<majorRadius<<","<<minorRadius<<" rotate<90,0,0> translate<"<<x<<","<<y<<"> pigment{"<<color<<"}}\n";
	//torus {1.0,0.25 scale <1,1,1> rotate<0,0,0> translate<0,0,0>
}
