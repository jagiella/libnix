/*
 * PovrayIO.h
 *
 *  Created on: 05.12.2009
 *      Author: jagiella
 */

#ifndef POVRAYIO_H_
#define POVRAYIO_H_

#include <fstream>

/** \addtogroup IO
 *  @{
 */

/** is a static class providing a number of functions to write Povray files */
class PovrayIO {
private:
	PovrayIO();

public:
	static void writePovrayHeader( std::fstream *fs, double minx, double miny, double maxx, double maxy);
	static void writePovrayHeader( std::fstream *fs, double minx, double miny, double maxx, double maxy, double z);
	static void writeSphere( std::fstream *fs, double x, double y, double radius, const char *color);
	static void writeSphere( std::fstream *fs, double x, double y, double z, double radius, const char *color);
	static void writeCylinder( std::fstream *fs, double x1, double y1, double x2, double y2, double radius, const char *color);
	static void writeCylinder( std::fstream *fs, double x1, double y1, double z1, double x2, double y2, double z2, double radius, const char *color);
	static void writeCube( std::fstream *fs, double x1, double y1, double z1, double x2, double y2, double z2, const char *color);
	static void writeCircle( std::fstream *fs, double x, double y, double majorRadius, double minorRadius, const char *color);
	template <class T>
	static void writePolygon( std::fstream *fs, int numberVertex, T *x, T *y, T *z, const char *color);
	static void writePrism( std::fstream *fs,
			double height1, double height2,
			int numberVertex, double *x, double *y, double *z,
			const char *color);

	static void beginIntersectionWithBox( std::fstream *fs, double x1, double y1, double z1, double x2, double y2, double z2);
	static void endIntersectionWithBox( std::fstream *fs);
};

/** @}*/

#include "Povray.ipp"

#endif /* POVRAYIO_H_ */
