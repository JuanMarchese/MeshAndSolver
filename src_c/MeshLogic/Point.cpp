/*
 * Point.cpp
 *
 *  Created on: 22/03/2014
 *      Author: juan
 */

#include "Point.h"

Point::Point():z(0),r(0),id(-1){
}
Point::Point(coordinates_t z, coordinates_t r,bool special):z(z),r(r),id(-1){
}
Point::Point(const Point& other){
	z = other.z;
	r = other.r;
	id = other.id;
}
Point& Point::operator=( const Point& other ){
	z = other.z;
	r = other.r;
	id = other.id;
	return *this;
}
Point Point::operator+( const Point& other )const{
	return Point(z+other.z,r+other.r);
}
Point Point::operator-( const Point& other )const{
	return Point(z-other.z,r-other.r);
}

bool Point::operator<(const Point& other) const{
	if (z == other.z){
		return r < other.r;
	}else{
		return z < other.z;
	}
}
bool Point::operator==(const Point& other) const{
	return (z == other.z) && (r == other.r);
}

Point::~Point() {
}

coordinates_t Point::getZ()const{
	return z;
}
coordinates_t Point::getR()const{
	return r;
}

void Point::setZ(coordinates_t z){
	this->z = z;
}
void Point::setR(coordinates_t r){
	this->r = r;
}
Point& Point::setID(ids_t value){
	id = value;
	return *this;
}
ids_t Point::getID()const{
	return id;
}

