/*
 * Quad.cpp
 *
 *  Created on: 22/03/2014
 *      Author: juan
 */

#include "Quad.h"
#include <iostream>

Quad::Quad(const Point::Vector& quadPoints):id(0),type(0){
	points[0] = quadPoints[0];
	points[1] = quadPoints[1];
	points[2] = quadPoints[2];
	points[3] = quadPoints[3];
}

Quad::Quad():id(0),type(0){

}

Quad::Quad(const Point& point1,const Point& point2,const Point& point3,const Point& point4,ids_t id):type(0){
	points[0] = point1;
	points[1] = point2;
	points[2] = point3;
	points[3] = point4;
	this->id = id;
}

Quad::~Quad() {
}

coordinates_t Quad::getX(unsigned pos)const{
	if(validPos(pos)){
		return points[pos].getZ();
	}else{
		return 0;
	}

}
coordinates_t Quad::getY(unsigned pos)const{
	if(validPos(pos)){
		return points[pos].getR();
	}else{
		return 0;
	}
}

Point Quad::getPoint(unsigned pos)const{
	if(validPos(pos)){
		return points[pos];
	}else{
		return Point();
	}
}

bool Quad::validPos(unsigned pos)const{
	return (pos <= 3);
}

void Quad::setType(ids_t value){
	type = value;
}
ids_t Quad::getType()const{
	return type;
}
void Quad::setID(ids_t value){
	id = value;
}
ids_t Quad::getID()const{
	return id;
}

coordinates_t Quad::getSurface(){
	coordinates_t surface = 0;

	for(index_t i = 0; i < 4 ; ++i)
	{
		surface += (this->points[(i+1)%4].getZ() - this->points[i].getZ())*
				   (this->points[(i+1)%4].getR() + this->points[i].getR());
	}

	return surface;
}
