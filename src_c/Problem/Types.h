/*
 * Types.h
 *
 */

#ifndef TYPES_H_
#define TYPES_H_

typedef double coordinates_t;
typedef unsigned index_t;
typedef long long int ids_t;

#define COORDINATE_PRECISION	0.0000001
#define SQR_TWO_OVER_TWO		0.707106781
#define COST_COEFICIENT			3.414213562
#define PI						3.14159265358979323846

#define INNER_CELL_MATERIAL			1
#define IN_MEMBRANE_MATERIAL		2
#define OUTER_CELL_MATERIAL			3
#define BAD_MATERIAL				-1

#endif /* TYPES_H_ */
