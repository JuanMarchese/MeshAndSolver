/*
 * Context.cpp
 *
 *  Created on: 26/04/2014
 *      Author: juan
 */

#include "Context.h"
#include <stdlib.h>

Context::Context(unsigned points_count) {

	points_size = points_count*2;//two coordinates, x and y

	// T = 2N - n - 2
	// T -> triangles
	// N -> points
	// n -> outer points (min is 3)

	triangles_size = (points_count*2 - 3 - 2) * 3;//Each triangle has

	points = (coordinates_t*) malloc(sizeof(coordinates_t) * points_size);
	triangles = (index_t*) malloc(sizeof(index_t) * points_size);

	points_instanced = 0;
	triangles_instanced = 0;
}

Context::~Context() {
	free(points);
	free(triangles);
}

coordinates_t* Context::getPointsVector(){
	return points;
}
index_t* Context::getTrianglesVector(){
	return triangles;
}

index_t Context::getNextPointID(){
	return points_instanced++;
}
index_t Context::getNextTriangleID(){
	return triangles_instanced++;
}

index_t Context::getPointsSize(){
	return points_size;
}
index_t Context::getTrianglesSize(){
	return triangles_size;
}
