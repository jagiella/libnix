/*
 * PovrayIO.h
 *
 *  Created on: 05.12.2009
 *      Author: jagiella
 */

#ifndef POVRAYIO_H_
#define POVRAYIO_H_

#include <fstream>

class PovrayIO {
public:
	static void writePovrayHeader( std::fstream *fs, double minx, double miny, double maxx, double maxy);
	static void writeSphere( std::fstream *fs, double x, double y, double radius, const char *color);
	static void writeSphere( std::fstream *fs, double x, double y, double z, double radius, const char *color);
	static void writeCylinder( std::fstream *fs, double x1, double y1, double x2, double y2, double radius, const char *color);
	static void writeCircle( std::fstream *fs, double x, double y, double majorRadius, double minorRadius, const char *color);
};

#endif /* POVRAYIO_H_ */
