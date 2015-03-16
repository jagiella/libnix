/*
 * Spring.h
 *
 *  Created on: 10.06.2013
 *      Author: jagiella
 */

#ifndef SPRING_H_
#define SPRING_H_

template <class elemType>
class Spring {
private:
	elemType a;
	elemType b;
	float _distance;
	float _stiffness;
public:
	Spring( elemType a = 0, elemType b = 0, float distance = 0, float stiffness = 0);
	elemType get( int i);
	float distance();
	float &stiffness() { return _stiffness;};
	virtual ~Spring();
};


#include "Spring.ipp"

#endif /* SPRING_H_ */
