/*
 * DataVertex.hpp
 *
 *  Created on: 09.05.2012
 *      Author: jagiella
 */

#ifndef DATAVERTEX_HPP_
#define DATAVERTEX_HPP_

#include "Triangulation.hpp"

template <class DataType>
class DataContainer {
private:
	DataType *data;
public:
	DataContainer() : data(0) {};
	DataContainer(DataType *data) : data(data) {};
	void set( DataType *data);
	DataType* get();
};

template <class DataContainerType>
class Data {
private:
	DataContainerType **associatedVertices;
	int countAssociatedVertices;
public:
	Data() : associatedVertices(0),countAssociatedVertices(0) {};
	void attach( DataContainerType *vertex);
	void detach( DataContainerType *vertex);
	void detach();
	DataContainerType* get( int i);
	int count();
};

#include "DataVertex.ipp"

#endif /* DATAVERTEX_HPP_ */
