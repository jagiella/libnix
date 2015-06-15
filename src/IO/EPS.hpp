/*
 * EpsIO.h
 *
 *  Created on: Dec 4, 2009
 *      Author: jagiella
 */

#ifndef EPSIO_H_
#define EPSIO_H_

#include <fstream>

/** \addtogroup IO
 *  @{
 */

class EPS {
private:
	EPS();

public:
	static void PSwriteHeader( std::fstream *fs, double minx, double miny, double maxx, double maxy);
	static void PSwriteSphere( std::fstream *fs, double x, double y, double radius, char *color);
	static void PSwriteLine( std::fstream *fs, double x1, double y1, double x2, double y2, double radius, char *color);
	static void PSwriteCircle( std::fstream *fs, double x, double y, double radius, double lineWidth, char *color);
	static void PSfillCircle( std::fstream *fs, double x, double y, double majorRadius, double minorRadius, char *color);
	template <typename T>
	static void PSwritePolygon( std::fstream *fs,
			T *x, T *y, int length,
			double radius, char *color);

	static void PSfillPolygon( std::fstream *fs,
			double *x, double *y, int length,
			double radius, char *color);

	static void PScolor( char color[], double R, double G, double B);
};

/** @}*/

#include "EPS.ipp"

#endif /* EPSIO_H_ */
