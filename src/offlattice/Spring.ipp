/*
 * Spring.cpp
 *
 *  Created on: 10.06.2013
 *      Author: jagiella
 */

#include "Spring.hpp"

/*template <class elemType>
Spring<elemType>::Spring( elemType a = 0, elemType b = 0, float distance = -1) :
	a(a), b(b), _distance(distance)
{
}*/

template <class elemType>
Spring<elemType>::Spring( elemType a, elemType b, float distance, float stiffness) :
	a(a), b(b), _distance(distance), _stiffness(stiffness)
{
}


template <class elemType>
Spring<elemType>::~Spring() {
	// TODO Auto-generated destructor stub
}

template <class elemType>
elemType Spring<elemType>::get( int i)
{
	if( i==0)
		return a;
	else
		return b;
}

template <class elemType>
float Spring<elemType>::distance()
{
	return _distance;
}
