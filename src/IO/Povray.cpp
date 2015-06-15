/*
 * PovrayIO.cpp
 *
 *  Created on: 05.12.2009
 *      Author: jagiella
 */

#include "Povray.hpp"

void PovrayIO::writePovrayHeader( std::fstream *fs, double minx, double miny, double maxx, double maxy)
{
	(*fs) << "#include \"colors.inc\"\n"
	      << "background { color White }\n"
	      << "camera {";
	(*fs) << "location <"<<(maxx-minx)/2.<<","<<(maxy-miny)/2.<<","<<(maxx-minx)/2.+(maxy-miny)/2.<<">\n";
	(*fs) << "look_at  <"<<(maxx-minx)/2.<<","<<(maxy-miny)/2.<<","<<0.<<">";
	(*fs) << "}\n";
	(*fs) << "light_source { <"<<(maxx-minx)/2.<<","<<(maxy-miny)/2.<<","<<(maxx-minx)/2.+(maxy-miny)/2.<<"> color White}\n";
	(*fs) << "global_settings { ambient_light rgb<0.4, 0.4, 0.4> }\n";
}


void PovrayIO::writePovrayHeader( std::fstream *fs, double minx, double miny, double maxx, double maxy, double z)
{
	(*fs) << "#include \"colors.inc\"\n"
	      << "background { color White }\n"
	      << "camera {";
	(*fs) << "location <"<<(maxx-minx)/2.<<","<<(maxy-miny)/2.<<","<<z<<">\n";
	(*fs) << "look_at  <"<<(maxx-minx)/2.<<","<<(maxy-miny)/2.<<","<<0.<<">";
	(*fs) << "}\n";
	(*fs) << "light_source { <"<<(maxx-minx)/2.<<","<<(maxy-miny)/2.<<","<<(maxx-minx)/2.+(maxy-miny)/2.<<"> color White}\n";
	(*fs) << "global_settings { ambient_light rgb<0.4, 0.4, 0.4> }\n";
}


void PovrayIO::writeSphere( std::fstream *fs, double x, double y, double z, double radius, const char *color)
{
	(*fs) << "sphere { <"<<x<<","<<y<<","<<z<<">,"<<radius<<" texture {pigment { color "<<color<<"}} finish{ambient 1}}\n";
}

void PovrayIO::writeCylinder( std::fstream *fs, double x1, double y1, double x2, double y2, double radius, const char *color)
{
	(*fs) << "cylinder{ <"<<x1<<","<<y1<<">,<"<<x2<<","<<y2<<">,"<<radius<<" pigment{color "<<color<<"} finish{ambient 1}}\n";
}

void PovrayIO::writeCylinder( std::fstream *fs, double x1, double y1, double z1, double x2, double y2, double z2, double radius, const char *color)
{
	(*fs) << "cylinder{ <"<<x1<<","<<y1<<","<<z1<<">,<"<<x2<<","<<y2<<","<<z2<<">,"<<radius<<" pigment{color "<<color<<"} finish{ambient 1}}\n";
}

void PovrayIO::writeCube( std::fstream *fs, double x1, double y1, double z1, double x2, double y2, double z2, const char *color)
{
	(*fs) << "box{ <"<<x1<<","<<y1<<","<<z1<<">,<"<<x2<<","<<y2<<","<<z2<<">, pigment{color "<<color<<"} finish{ambient 1}}\n";
}

void PovrayIO::writeCircle( std::fstream *fs, double x, double y, double majorRadius, double minorRadius, const char *color)
{
	(*fs) << "torus { "<<majorRadius<<","<<minorRadius<<" rotate<90,0,0> translate<"<<x<<","<<y<<"> pigment{"<<color<<"} finish{ambient 1}}\n";
	//torus {1.0,0.25 scale <1,1,1> rotate<0,0,0> translate<0,0,0>
}


void PovrayIO::writePrism( std::fstream *fs,
		double height1, double height2,
		int numberVertex, double *x, double *y, double *z,
		const char *color)
{
	int i=0;
	(*fs) << "prism {"<<height1<<","<<height2<<","<< numberVertex+1;
	while(i<numberVertex)
	{
		(*fs) << ",<" << x[i] << "," << y[i] << "," << z[i] << ">";
		++i;
	}
	(*fs) << ",<" << x[0] << "," << y[0] << "," << z[0] << ">";
	(*fs) << "pigment{color "<<color<<"} finish{ambient 1}}\n";
}

void PovrayIO::beginIntersectionWithBox( std::fstream *fs, double x1, double y1, double z1, double x2, double y2, double z2)
{
	(*fs) << "intersection{\n";
	(*fs) << "box{<"<<x1<<","<<y1<<","<<z1<<">,<"<<x2<<","<<y2<<","<<z2<<">}\n";
	(*fs) << "union{\n";
}

void PovrayIO::endIntersectionWithBox( std::fstream *fs)
{
	(*fs) << "}}\n";
}

