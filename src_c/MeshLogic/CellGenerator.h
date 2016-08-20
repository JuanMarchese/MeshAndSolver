/*
 * CellGenerator.h
 *
 *  Created on: 13/7/2014
 *      Author: juan
 */

#ifndef CELLGENERATOR_H_
#define CELLGENERATOR_H_

#include <vector>
#include <map>
#include "Data.h"

class CellGenerator {

private:
	std::vector<Point> firstInnerPoints;
	std::vector<Point> realOuterPoints;
	std::vector<Point> realInnerPoints;
	std::vector<Point> innerPoints;
	std::vector<Point> outerPoints;

	coordinates_t innerRadius;
	coordinates_t outerRadius;

	coordinates_t cellZCenter;
	coordinates_t cellRCenter;

	index_t cellLayers;
	index_t cellIterations;

	index_t boundaryLayers;
	coordinates_t boundaryLinesSeparation;

	index_t innerLayers;
	coordinates_t innerLinesSeparation;

	Data* data;

	std::vector<coordinates_t> angles;
	std::vector<coordinates_t> radius;
	std::vector<coordinates_t> outerRadiusVector;
	std::vector<coordinates_t> innerRadiusVector;

	ids_t inCellType;
	ids_t outerCellType;
	ids_t innerCellType;

public:
	CellGenerator();
	virtual ~CellGenerator();

	void initData(Data* data);
	void initCellGeometrics(coordinates_t innerRadius,
						coordinates_t outerRadius);
	void initCellParameters(index_t cellLayers,
							index_t cellIterations);
	void initBoundaryParameters(index_t boundaryLayers,
								coordinates_t boundaryLinesSeparation);

	void initInnerParameters(index_t boundaryLayers,
							 coordinates_t boundaryLinesSeparation);

	void initTypes(ids_t inCell,ids_t outerCell,ids_t innerCell);

	void process();
	std::vector<Point>& getOuterPoints();
	std::vector<Point>& getInnerPoints();

	bool circlesOverlap();

	coordinates_t getInnerRadius();
	coordinates_t getOuterRadius();
	index_t getCellIterations();

private:
	void generateRadiusVector();
	void generateOuterRadiusVector();
	void generateCellCirclePoints(std::vector<Point>& vector,index_t iteration);
	void generateBoundaryCirclePoints(std::vector<Point>& vector,index_t iteration);
	void generateQuads(ids_t type = 0,bool saveIDs = false,ids_t layer = 0);
	void generateAngles();
	void generateFinalVector(index_t start,
							 index_t increment,
							 index_t iterations,
							 std::vector<Point>& vector,
							 std::vector<coordinates_t> vectorZ,
							 std::vector<coordinates_t> vectorR,
							 coordinates_t zSign,
							 coordinates_t rSign,
							 bool addToData = true);

	void generateInnerRadiusVector();
	void generateInnerCirclePoints(std::vector<Point>& vector,index_t iteration,bool addPointsToData = true);

	void addPointsToData(std::vector<Point>& newPoints);
};

#endif /* CELLGENERATOR_H_ */

