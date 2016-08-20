/*
 * Context.h
 */

#ifndef CONTEXT_H_
#define CONTEXT_H_

#include "Types.h"

class Context {
public:
	Context(unsigned points_count);
	~Context();

	coordinates_t* getPointsVector();
	index_t* getTrianglesVector();

	index_t getPointsSize();
	index_t getTrianglesSize();

	index_t getNextPointID();
	index_t getNextTriangleID();

private:
	index_t points_size;
	index_t triangles_size;

	coordinates_t* points;
	index_t* triangles;

	index_t points_instanced;
	index_t triangles_instanced;


};

#endif /* CONTEXT_H_ */
