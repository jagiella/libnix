/*
 * Color.cpp
 *
 *  Created on: 17.10.2013
 *      Author: jagiella
 */

#include <math.h>

#define MINMAX(  a,  b,  c) \
	( b<a?a: (b>c?c : b))

double rgbformulae(char formulae, double x) {
	double ret = 0;
	double pi = 3.14159265;
	switch ( (int)fabs(formulae)) {
	case 0:
		ret = 0.;
		break;
	case 1:
		ret = 0.5;
		break;
	case 2:
		ret = 1.;
		break;
	case 3:
		ret = x;
		break;
	case 4:
		ret = pow(x, 2.);
		break;
	case 5:
		ret = pow(x, 3.);
		break;
	case 6:
		ret = pow(x, 4.);
		break;
	case 7:
		ret = sqrt(x);
		break;
	case 8:
		ret = sqrt(sqrt(x));
		break;
	case 9:
		ret = sin(M_PI / 2. * x);
		break;
	case 10:
		ret = cos(M_PI / 2. * x);
		break;
	case 11:
		ret = fabs(x - 0.5);
		break;
	case 12:
		ret = pow(2. * x - 1., 2.);
		break;
	case 13:
		ret = sin(pi * x);
		break;
	case 14:
		ret = fabs(cos(pi * x));
		break;
	case 15:
		ret = sin(2. * pi * x);
		break;
	case 16:
		ret = cos(2. * pi * x);
		break;
	case 17:
		ret = fabs(sin(2. * pi * x));
		break;
	case 18:
		ret = fabs(cos(2. * pi * x));
		break;
	case 19:
		ret = fabs(sin(4. * pi * x));
		break;
	case 20:
		ret = fabs(cos(4. * pi * x));
		break;
	case 21:
		ret = 3. * x;
		break;
	case 22:
		ret = 3. * x - 1.;
		break;
	case 23:
		ret = 3. * x - 2.;
		break;
	case 24:
		ret = fabs(3. * x - 1.);
		break;
	case 25:
		ret = fabs(3. * x - 2.);
		break;
	case 26:
		ret = (3. * x - 1.) / 2.;
		break;
	case 27:
		ret = (3. * x - 2.) / 2.;
		break;
	case 28:
		ret = fabs((3. * x - 1.) / 2.);
		break;
	case 29:
		ret = fabs((3. * x - 2.) / 2.);
		break;
	case 30:
		ret = x / 0.32 - 0.78125;
		break;
	case 31:
		ret = 2. * x - 0.84;
		break;
		//case 32: ret = 4x;1;-2x+1.84;x/0.08-11.5; break;
	case 33:
		ret = fabs(2. * x - 0.5);
		break;
	case 34:
		ret = 2. * x;
		break;
	case 35:
		ret = 2. * x - 0.5;
		break;
	case 36:
		ret = 2. * x - 1.;
		break;
	}

	if (formulae < 0)
		return MINMAX(0., -ret, 1.);
	else
		return MINMAX(0., ret, 1.);
}

/*
color_t rgbformulaeMapping(char rformulae, char gformulae, char bformulae,
		double x, double min, double max) {
	color_t color = 0;
	char *formulae[3];
	formulae[0] = &rformulae;
	formulae[1] = &gformulae;
	formulae[2] = &bformulae;

	for (int c = 0; c < 3; c++) {
		//fprintf(stderr, "%lf\n", rgbformulae( *formulae[c], (x-min) / (max-min)));
		((channel_t*) (&color))[c] = (channel_t) 255. * rgbformulae(
				*formulae[c], (x - min) / (max - min));
	}

	return color;
}*/
