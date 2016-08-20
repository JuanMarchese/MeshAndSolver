/*
 * Point.cpp
 *
 *  Created on: 22/03/2014
 *      Author: juan
 */

#include "Point.h"

Point::Point():x(0),y(0),id(-1){
}
Point::Point(coordinates_t x, coordinates_t y,bool special):x(x),y(y),id(-1){
}
Point::Point(const Point& other){
	x = other.x;
	y = other.y;
	id = other.id;
}
Point& Point::operator=( const Point& other ){
	x = other.x;
	y = other.y;
	id = other.id;
	return *this;
}
Point Point::operator+( const Point& other )const{
	return Point(x+other.x,y+other.y);
}
Point Point::operator-( const Point& other )const{
	return Point(x-other.x,y-other.y);
}

bool Point::operator<(const Point& other) const{
	if (x == other.x){
		return y < other.y;
	}else{
		return x < other.x;
	}
}
bool Point::operator==(const Point& other) const{
	return (x == other.x) && (y == other.y);
}

Point::~Point() {
}

coordinates_t Point::getX()const{
	return x;
}
coordinates_t Point::getY()const{
	return y;
}

void Point::setX(coordinates_t x){
	this->x = x;
}
void Point::setY(coordinates_t y){
	this->y = y;
}
Point& Point::setID(ids_t value){
	id = value;
	return *this;
}
ids_t Point::getID()const{
	return id;
}

