/*
 * Point.h
 *
 *  Created on: 22/03/2014
 *      Author: juan
 */

#ifndef POINT_H_
#define POINT_H_

#include <set>
#include <vector>
#include "../Problem/Types.h"

class Point {
public:

	typedef std::vector<Point> Vector;

	Point();
	Point(const Point& other);
	Point& operator=(const Point& other);
	Point(coordinates_t x, coordinates_t y,bool special = false);

	Point operator+( const Point& other )const;
	Point operator-( const Point& other )const;
	bool operator<(const Point& v) const;
	bool operator==(const Point& v) const;

	~Point();

	coordinates_t getX()const;
	coordinates_t getY()const;

	void setX(coordinates_t x);
	void setY(coordinates_t y);

	Point& setID(ids_t value);
	ids_t getID()const;


private:
	coordinates_t x;
	coordinates_t y;
	ids_t id;


};

#endif /* POINT_H_ */
