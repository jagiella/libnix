/*
 * Povray.ipp
 *
 *  Created on: 30.04.2012
 *      Author: jagiella
 */

#ifndef POVRAY_IPP_
#define POVRAY_IPP_

template <class T>
void PovrayIO::writePolygon( std::fstream *fs, int numberVertex, T *x, T *y, T *z, const char *color)
{
	int i=0;
	(*fs) << "polygon {"<< numberVertex+1;
	while(i<numberVertex)
	{
		(*fs) << ",<" << x[i] << "," << y[i] << "," << z[i] << ">";
		++i;
	}
	(*fs) << ",<" << x[0] << "," << y[0] << "," << z[0] << ">";
	(*fs) << "pigment{color "<<color<<"} finish{ambient 1}}\n";
}

#endif /* POVRAY_IPP_ */
