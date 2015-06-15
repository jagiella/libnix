/*
 * EPS.ipp
 *
 *  Created on: 30.04.2012
 *      Author: jagiella
 */

#ifndef EPS_IPP_
#define EPS_IPP_

template <typename T>
void EPS::PSwritePolygon( std::fstream *fs,
		T *x, T *y, int length,
		double radius, char *color) {

	(*fs) << "newpath\n";
	(*fs) << "0.01 setlinewidth\n ";
	(*fs) << " "<<x[0]<<" "<<y[0]<<" moveto\n";

	for(int i=1; i<length; i++)
	(*fs) << " "<<x[i]<<" "<<y[i]<<" lineto\n";

	(*fs) << " closepath\n " << color << "\n";
}


#endif /* EPS_IPP_ */
