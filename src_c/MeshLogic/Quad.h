/*
 * Quad.h
 *
 *  Created on: 22/03/2014
 *      Author: juan
 */

#ifndef QUAD_H_
#define QUAD_H_

#include <vector>

#include "Point.h"
#include "../Problem/Types.h"

class Quad {
public:

	typedef std::vector<Quad> Vector;

	Quad(const Point::Vector& quadPoints);
	Quad();
	Quad(const Point& point1,const Point& point2,const Point& point3,const Point& point4,ids_t id=0);
	virtual ~Quad();

	coordinates_t getX(unsigned pos)const;
	coordinates_t getY(unsigned pos)const;
	Point getPoint(unsigned pos)const;

	void setType(ids_t value);
	ids_t getType()const;

	void setID(ids_t value);
	ids_t getID()const;

	coordinates_t getSurface();

private:
	bool validPos(unsigned pos)const;
	Point points[4];
	ids_t id;
	ids_t type;

};

#endif /* QUAD_H_ */
