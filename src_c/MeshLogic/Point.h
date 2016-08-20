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
	Point(coordinates_t z, coordinates_t r,bool special = false);

	Point operator+( const Point& other )const;
	Point operator-( const Point& other )const;
	bool operator<(const Point& v) const;
	bool operator==(const Point& v) const;

	~Point();

	coordinates_t getZ()const;
	coordinates_t getR()const;

	void setZ(coordinates_t z);
	void setR(coordinates_t r);

	Point& setID(ids_t value);
	ids_t getID()const;


private:
	coordinates_t z;
	coordinates_t r;
	ids_t id;


};

#endif /* POINT_H_ */
