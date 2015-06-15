/*
 * DataVertex.ipp
 *
 *  Created on: 09.05.2012
 *      Author: jagiella
 */

#include <typeinfo>

template <class DataType>
void DataContainer<DataType>::set( DataType *data)
{
	this->data = data;
}

template <class DataType>
DataType* DataContainer<DataType>::get()
{
	return data;
}

template <class DataContainerType>
void Data<DataContainerType>::attach( DataContainerType *vertex)
{
	// link vertex -> data
	vertex->set( this);

	// link data -> vertex
	associatedVertices = (DataContainerType**) realloc(
			associatedVertices,
			++(countAssociatedVertices) * sizeof(DataContainerType*));
	this->associatedVertices[countAssociatedVertices-1] = vertex;
}

template <class DataContainerType>
void Data<DataContainerType>::detach( DataContainerType *vertex)
{
	// link vertex -> data
	vertex->set(0);

	// link data -> vertex
	for( int i=0; i<countAssociatedVertices; i++)
		if( this->associatedVertices[i] == vertex){
			countAssociatedVertices--;
			this->associatedVertices[i] = this->associatedVertices[countAssociatedVertices];
			return;
		}
	return;
}

template <class DataContainerType>
void Data<DataContainerType>::detach()
{
	for( int i=0; i<countAssociatedVertices; i++){
		// link vertex -> data
		associatedVertices[i]->set(0);
		// link data -> vertex
	}
	countAssociatedVertices=0;
}


template <class DataContainerType>
DataContainerType* Data<DataContainerType>::get( int i)
{
	return associatedVertices[i];
}

template <class DataContainerType>
int Data<DataContainerType>::count()
{
	return countAssociatedVertices;
}
