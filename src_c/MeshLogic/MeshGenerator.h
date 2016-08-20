/*
 * MeshGenerator.h
 *
 *  Created on: 9/7/2014
 *      Author: juan
 */

#ifndef MESHGENERATOR_H_
#define MESHGENERATOR_H_

#include <vector>
#include <map>
#include "Data.h"

class MeshGenerator {
private:
	std::vector<Point> originalInnerPoints;
	std::vector<Point> innerPoints;
	std::vector<Point> outerPoints;

	std::vector<Point> inCellPoints;

	index_t outerRowsCount;
	index_t outerColumsCount;
	index_t innerRowsCount;
	index_t innerColumsCount;

	coordinates_t zCenter;
	coordinates_t rCenter;

	coordinates_t blockSize;
	index_t blockZIterations;
	index_t blockRIterations;
	index_t numberOfLayers;

	index_t numberOfFixingLayers;
	std::vector<coordinates_t> fixingIncremenst;
	std::vector<coordinates_t> sizes;

	Data* data;

	ids_t outerCellType;
	ids_t innerCellType;

public:
	MeshGenerator();
	virtual ~MeshGenerator();

	void initData(Data* data);
	void initBlocksSizes(coordinates_t blockSize,index_t blockXIterations,index_t blockYIterations,index_t numberOfLayers);
	void initCenterPosition(coordinates_t zCenter,coordinates_t rCenter);
	void initTransitionData(index_t numberOfFixingLayers);
	void initTypes(ids_t outerCell,ids_t innerCell);
	void process();
	void setInnerVector(std::vector<Point>& innerPoints);
	void setInCellVector(std::vector<Point>& inCellPoints);

private:
	coordinates_t calculateAverageSize(std::vector<Point>& vector);
	coordinates_t getMinSize(std::vector<Point>& vector);
	void processBoundary(coordinates_t blockSize,index_t iterations);
	void processInnerCell(index_t iterations);
	bool addBlocks(coordinates_t blockSize,ids_t type = 0);
	void generateQuads(ids_t type = 0);
	void calculateSizes();
	void calculateFixIncrements();
	void calculateBlockSize();
	void calculateFixingLayers();
	void generateFixPoints();
	void generatePoints(std::vector<Point>& vector,
					    Point start,
					    coordinates_t incrementZ,
					    coordinates_t incrementR,
					    index_t iterations,
					    coordinates_t fixFirstIncrementZ = 0.0,
					    coordinates_t fixFirstIncrementR = 0.0);
	void generateFixedPoints(std::vector<Point>& vector,
							 Point start,
							 coordinates_t increment,
							 coordinates_t mainIncrementMultiplier,
							 bool fixX,
							 index_t iterations,
							 index_t globalIteration,
							 index_t& fixedPoints);
	void generateFixedPointsInverseOrder(std::vector<Point>& vector,
							 	 	 	 Point start,
							 	 	 	 coordinates_t increment,
							 	 	 	 coordinates_t mainIncrementMultiplier,
							 	 	 	 bool fixX,
							 	 	 	 index_t iterations,
							 	 	 	 index_t globalIteration,
							 	 	 	 index_t& fixedPoints);

	void addDoubleSizeExternalTransition();
	void addTripleSizeExternalTransition();
	void addDoubleSizeInternalTransition();
	void addTripleSizeInternalTransition();

	void verticalEndQuads();
	void horizontalEndQuads();

	void addRow(std::vector<Point>& points, coordinates_t increment,ids_t type = 0);
	void addColumn(std::vector<Point>& points, coordinates_t increment,bool left,ids_t type = 0);
	void addInnerColumn(std::vector<Point>& points, coordinates_t increment,bool left,ids_t type = 0);
	void generateTransitionTwoBlocs(std::vector<Point>& points,
								    coordinates_t sizeZ,
								    coordinates_t sizeR,
								    std::vector<Point>& outPointsVector,
								    index_t start,
								    index_t end,
								    Point firstPoint,
								    ids_t type = 0,
								    bool inverseOrder = false);

	void addTowBlocksTransition(std::vector<Point>& points,coordinates_t size);
	index_t findTransition(std::vector<Point>& points,index_t start);

	void addInnerLayer(std::vector<Point>& points,coordinates_t size,bool hasPadding);
	void addOuterLayer(std::vector<Point>& points,coordinates_t size);
	void addInnerTransition(std::vector<Point>& points,coordinates_t size,std::vector<Point>& outPoints);

	void processOuterFixingLayer(coordinates_t zStart,coordinates_t rStart);
	bool thereIsPlaceToContinue(coordinates_t size);
	index_t possibleIterations(coordinates_t size);
	void finishBoundary(coordinates_t size);
	void setUpBoundaries();
};

#endif /* MESHGENERATOR_H_ */
