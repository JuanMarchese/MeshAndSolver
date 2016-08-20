/*
 * MeshGenerator.cpp
 *
 *  Created on: 9/7/2014
 *      Author: juan
 */

#include "MeshGenerator.h"
#include <iostream>
#include <cmath>

MeshGenerator::MeshGenerator():	outerRowsCount(0),
								outerColumsCount(0),
								innerRowsCount(0),
								innerColumsCount(0),
								xCenter(0),
								yCenter(0),
							    blockSize(0),
							    blockXIterations(0),
							    blockYIterations(0),
							    numberOfLayers(0),
							    numberOfFixingLayers(0),
							    data(0),
							    outerCellType(0),
							    innerCellType(0) {
}

MeshGenerator::~MeshGenerator() {
}

void MeshGenerator::initBlocksSizes(coordinates_t blockSize,
									index_t blockXIterations,
									index_t blockYIterations,
									index_t numberOfLayers){
	this->blockSize = blockSize;
	this->blockXIterations = blockXIterations;
	this->blockYIterations = blockYIterations;
	this->numberOfLayers = numberOfLayers;
}
void MeshGenerator::initCenterPosition(coordinates_t xCenter,coordinates_t yCenter){
	this->xCenter = xCenter;
	this->yCenter = yCenter;
}
void MeshGenerator::setInnerVector(std::vector<Point>& innerPoints){


	this->innerPoints = innerPoints;
	this->originalInnerPoints = innerPoints;
}

void MeshGenerator::setInCellVector(std::vector<Point>& inCellPoints){
	this->inCellPoints = inCellPoints;
}

void MeshGenerator::process(){

	std::vector<Point> innerPointsTemp;
	coordinates_t x;
	coordinates_t y;

	/*Init the paremeters*/

	calculateBlockSize();
	calculateSizes();
	calculateFixingLayers();
	calculateFixIncrements();

	/*Hack to inverse the order in the points*/
	for(index_t i = innerPoints.size(); i > 0;--i){
		innerPointsTemp.push_back(innerPoints[i-1]);
	}

	innerPoints = innerPointsTemp;
	x = xCenter - blockSize*(blockXIterations/2+1);
	y = yCenter;

	/*Generate the outer fixing layers points*/
	processOuterFixingLayer(x,y);

	/*Generate the outer layers until the end*/

	coordinates_t size = blockSize;
	while(thereIsPlaceToContinue(size)){
		processBoundary(size,4);
		size = size * 2;
	}
	/*Is not size/2 because processBoundary ends whit the double of the size*/
	finishBoundary(size);

	setUpBoundaries();

	innerPointsTemp.clear();
	for(index_t i = inCellPoints.size(); i > 0;--i){
		innerPointsTemp.push_back(inCellPoints[i-1]);
	}
	inCellPoints = innerPointsTemp;
	this->innerPoints = inCellPoints;

	coordinates_t sizeTemp = calculateAverageSize(this->innerPoints);
	addInnerLayer(this->innerPoints,sizeTemp,true);



	/*MIN for transition: 8 up; 4 sides*/
	/*In this way we eliminate the rest needed for the padding*/
	index_t side = findTransition(this->innerPoints,0);
	index_t up = findTransition(this->innerPoints,side+1) - side;

	index_t sideEfective =  (side / 2)*2;
	index_t upEfective = (up / 4)*4;

	/*UP + 4  because there is a chance that we will need padding*/
	/*SIDE + 2  because there is a chance that we will need padding*/

	/*+2 for the iterations pases to the function*/
	while((sideEfective > (4+2)) && (upEfective > (8+2))){
		processInnerCell(2);
		side = findTransition(this->innerPoints,0);
		up = findTransition(this->innerPoints,side+1) - side;
		sideEfective =  (side / 2)*2;
		upEfective = (up / 4)*4;
	}

	index_t endIterations;

	if(side > up)
		endIterations = up - 1;
	else
		endIterations = side - 1;



	for(index_t i = 0; i < endIterations; ++i)
	{
		addInnerLayer(this->innerPoints,this->innerPoints[1].getY()- this->innerPoints[0].getY(),false);
	}

	side = findTransition(this->innerPoints,0);
	up = findTransition(this->innerPoints,side+1) - side;

	if(side > up){
		horizontalEndQuads();
	}else{
		verticalEndQuads();
	}

}

void MeshGenerator::processOuterFixingLayer(coordinates_t xStart,coordinates_t yStart){

	index_t fixedPoints;

	for(index_t i = 0; i < numberOfFixingLayers ; ++i){
		this->outerPoints.clear();
		fixedPoints = 0;

		generateFixedPoints(this->outerPoints,Point(xStart,yStart),
					   	    0.0,1.0,false,blockYIterations,i,fixedPoints);
		generatePoints(this->outerPoints,
					   this->outerPoints[0] + Point(0.0,blockSize*blockYIterations),
					   0.0,blockSize,i + 1);
		generatePoints(this->outerPoints,
					   this->outerPoints[this->outerPoints.size()-1] + Point(blockSize,0.0),
					   blockSize,0.0,i);
		generateFixedPointsInverseOrder(this->outerPoints,Point(xCenter,(blockYIterations+i)*blockSize),
										0.0,1.0,true,blockXIterations/2+1,i,fixedPoints);

		Point point(xCenter,(blockYIterations+i)*blockSize);
		point.setID(data->addPoint(point));
		this->outerPoints.push_back(point);
		--fixedPoints;

		generateFixedPoints(this->outerPoints,Point(xCenter,(blockYIterations+i)*blockSize),
					   	    0.0,1.0,true,blockXIterations/2,i,fixedPoints);
		generatePoints(this->outerPoints,
					   Point(xCenter+blockSize*(blockXIterations/2+1),this->outerPoints[this->outerPoints.size()-1].getY()),
					   blockSize,0.0,i+1);
		generatePoints(this->outerPoints,
					   this->outerPoints[this->outerPoints.size()-1] + Point(0.0,-1.0*blockSize),
					   0.0,-1.0*blockSize,i);
		generateFixedPointsInverseOrder(this->outerPoints,
										Point(this->outerPoints[this->outerPoints.size()-1].getX(),yCenter),
										0.0,-1.0,false,blockYIterations+1,i,fixedPoints);
		xStart -= blockSize;

		generateQuads(outerCellType);
		this->innerPoints = this->outerPoints;
	}
}


void MeshGenerator::processInnerCell(index_t iterations){

	std::vector<Point> innerPointsTemp;
	index_t firstPointChange;
	index_t middlePointChange;


	firstPointChange = findTransition(this->innerPoints,0);
	middlePointChange = findTransition(this->innerPoints,firstPointChange+1);

	index_t rest = (middlePointChange - firstPointChange) % 4;

	coordinates_t size = calculateAverageSize(this->innerPoints);


	/*Test if the top is 4 multiple*/
	if(rest != 0)
	{
		if(rest > 1){
			/*The inner layer extracts two blocks to the top*/
			addInnerLayer(this->innerPoints,size,true);
		}

		if(rest % 2 != 0)
			addInnerColumn(this->innerPoints, this->innerPoints[1].getY()- this->innerPoints[0].getY(),true,innerCellType);

		firstPointChange = findTransition(this->innerPoints,0);
		middlePointChange = findTransition(this->innerPoints,firstPointChange+1);
	}


	/*Test if the sides are 2 multiple*/
	if(firstPointChange % 2 != 0)
	{
		addRow(this->innerPoints, this->innerPoints[0].getY()- this->innerPoints[1].getY(),innerCellType);
		firstPointChange = findTransition(this->innerPoints,0);
		middlePointChange = findTransition(this->innerPoints,firstPointChange+1);
	}

	innerPointsTemp.clear();



	/*Add the transition*/
	addInnerTransition(this->innerPoints,innerPoints[1].getY()- innerPoints[0].getY(),innerPointsTemp);

	this->innerPoints = innerPointsTemp;

	/*Add an inner layer*/
	for(index_t i = 0; i < iterations;++i)
	{
		addInnerLayer(this->innerPoints,this->innerPoints[1].getY()- this->innerPoints[0].getY(),false);
	}

}

bool MeshGenerator::addBlocks(coordinates_t blockSize,ids_t type){
	index_t firstPointChange;
	index_t middlePointChange;
	coordinates_t increment;
	coordinates_t oldIncrement;
	bool retVal = true;

	this->outerPoints.clear();

	firstPointChange = findTransition(this->innerPoints,0);
	middlePointChange = findTransition(this->innerPoints,firstPointChange+1);


	if((this->innerPoints[0].getX() + blockSize)<(this->data->getMinX())){
		retVal = false;
		increment = this->data->getMinX() - this->innerPoints[0].getX();
	}else{
		increment = -1.0 * blockSize;
	}

	generatePoints(this->outerPoints,
				   this->innerPoints[0] + Point(increment,0.0),
				   0.0,blockSize,firstPointChange+1);


	oldIncrement = increment;
	if((this->innerPoints[firstPointChange].getY() + blockSize)>(this->data->getMaxY())){
		retVal = false;
		increment = this->data->getMaxY() - this->innerPoints[firstPointChange].getY();
	}else{
		increment = blockSize;
	}

	generatePoints(this->outerPoints,
				   this->outerPoints[firstPointChange] + Point(0.0,increment),
				   blockSize,0.0,middlePointChange-firstPointChange+2,
				   oldIncrement + blockSize,0.0);


	oldIncrement = increment;
	if((this->innerPoints[middlePointChange+2].getX() + blockSize)>(this->data->getMaxX())){
			retVal = false;
			increment = this->data->getMaxX() - this->innerPoints[middlePointChange+2].getX();
		}else{
			increment = blockSize;
		}

	generatePoints(this->outerPoints,
				   this->outerPoints[middlePointChange+2] + Point(increment,0.0),
				   0.0,-1.0*blockSize,
				   firstPointChange+2,
				   0.0,blockSize - oldIncrement);

	generateQuads(type);

	return retVal;
}

index_t MeshGenerator::possibleIterations(coordinates_t size){
	index_t iterations = 0;

	coordinates_t startX;
	coordinates_t endX;
	coordinates_t endY;

	index_t iterationsStartX;
	index_t iterationsEndX;
	index_t iterationsEndY;

	index_t firstPointChange;
	index_t middlePointChange;

	firstPointChange = findTransition(this->outerPoints,0);
	middlePointChange = findTransition(this->outerPoints,firstPointChange+1);

	startX = this->outerPoints[0].getX() - this->data->getMinX();
	endX = this->data->getMaxX() - this->outerPoints[middlePointChange].getX();
	endY = this->data->getMaxY() - this->outerPoints[firstPointChange].getY();

	iterationsStartX = (index_t)(startX / size);
	iterationsEndX = (index_t)(endX / size);
	iterationsEndY = (index_t)(endY / size);

	if(iterationsStartX < iterationsEndX){
		iterations = iterationsStartX;
	}else{
		iterations = iterationsEndX;
	}

	if(iterationsEndY < iterations){
		iterations = iterationsEndY;
	}

	return iterations;
}

void MeshGenerator::finishBoundary(coordinates_t size){
	coordinates_t startX;
	coordinates_t endX;
	coordinates_t endY;
	coordinates_t endIncrementXleft;
	coordinates_t endIncrementXright;
	coordinates_t endIncrementY;

	int iterationsStartX;
	int iterationsEndX;
	int iterationsEndY;

	int firstPointChange;
	int middlePointChange;

	firstPointChange = findTransition(this->innerPoints,0);
	middlePointChange = findTransition(this->innerPoints,firstPointChange+1);

	startX = this->innerPoints[0].getX() - this->data->getMinX();
	endX = this->data->getMaxX() - this->innerPoints[middlePointChange].getX();
	endY = this->data->getMaxY() - this->innerPoints[firstPointChange].getY();

	iterationsStartX = (int)(startX / size);
	iterationsEndX = (int)(endX / size);
	iterationsEndY = (int)(endY / size);

	endIncrementXleft = (this->innerPoints[0].getX() - (size*iterationsStartX)) - this->data->getMinX();
	endIncrementXright = this->data->getMaxX() - (this->innerPoints[middlePointChange].getX() + (iterationsEndX*size));
	endIncrementY = this->data->getMaxY() - (this->innerPoints[firstPointChange].getY() + (iterationsEndY*size));

	for(int i = 0 ; i < iterationsStartX; i++){
		addColumn(this->innerPoints,size*-1.0,true,outerCellType);
	}

	if(endIncrementXleft > 0.0)
		addColumn(this->innerPoints,endIncrementXleft*-1.0,true,outerCellType);


	for(int i = 0 ; i < iterationsEndX; i++){
		addColumn(this->innerPoints,size,false,outerCellType);
	}
	if(endIncrementXright > 0.0)
		addColumn(this->innerPoints,endIncrementXright,false,outerCellType);


	for(int i = 0 ; i < iterationsEndY; i++){
		addRow(this->innerPoints,size,outerCellType);
	}
	if(endIncrementY > 0.0)
		addRow(this->innerPoints,endIncrementY,outerCellType);
}

/*We need a minimum of 10 blockSize to work properly*/
bool MeshGenerator::thereIsPlaceToContinue(coordinates_t size){
	bool retVal = false;
	coordinates_t startX;
	coordinates_t endX;
	coordinates_t endY;

	index_t firstPointChange;
	index_t middlePointChange;

	firstPointChange = findTransition(this->outerPoints,0);
	middlePointChange = findTransition(this->outerPoints,firstPointChange+1);

	startX = this->outerPoints[0].getX() - (size * 6.0);
	endX = this->outerPoints[middlePointChange].getX() + (size * 6.0);
	endY = this->outerPoints[firstPointChange].getY() + (size * 6.0);

	retVal = (startX > this->data->getMinX()) &&
			 (endX < this->data->getMaxX()) &&
			 (endY < this->data->getMaxY());

	return retVal;
}

void MeshGenerator::processBoundary(coordinates_t blockSize,index_t iterations){


	index_t firstPointChange;
	index_t middlePointChange;
	index_t rest;

	this->innerPoints = this->outerPoints;
	firstPointChange = findTransition(this->innerPoints,0);
	middlePointChange = findTransition(this->innerPoints,firstPointChange+1);
	rest = firstPointChange % 4;

	/* Before add the transition we must make the amount of quads a multiple
	 * of the transition*/

	/*First horizontally*/
	if(rest != 0)
	{
		for(index_t i = 0; i < (4-rest);++i)
		{
			addBlocks(blockSize,outerCellType);
			this->innerPoints = this->outerPoints;
		}

		firstPointChange = findTransition(this->innerPoints,0);
		middlePointChange = findTransition(this->innerPoints,firstPointChange+1);
	}

	rest = (middlePointChange-firstPointChange) % 4;

	/*Then vertically*/
	if(rest != 0)
	{
		bool left = true;
		for(index_t i = 0; i < (4-rest);++i)
		{
			addColumn(this->innerPoints,-1.0*blockSize,left,outerCellType);
			left = !left;
		}
	}


	addTowBlocksTransition(this->outerPoints,blockSize);
	this->innerPoints = this->outerPoints;


	index_t maxIterations = possibleIterations(blockSize*2);

	if(maxIterations < iterations){
		iterations = maxIterations;
	}

	for(index_t i = 0; i < iterations; ++i)
	{
		addBlocks(blockSize*2,outerCellType);
		this->innerPoints = this->outerPoints;
	}

}

void MeshGenerator::initData(Data* data){
	this->data = data;
}

void MeshGenerator::generateQuads(ids_t type){

	index_t innerIterator;
	index_t outerIterator;

	bool verticalIncrement = std::abs(innerPoints[1].getY() - innerPoints[0].getY())
							 >std::abs(innerPoints[1].getX() - innerPoints[0].getX());

	bool newVerticalIncrement = verticalIncrement;

	outerIterator = 0;
	for(innerIterator = 0; innerIterator < innerPoints.size() - 1;++innerIterator){

		newVerticalIncrement =  std::abs(innerPoints[innerIterator+1].getY() - innerPoints[innerIterator].getY()) >
	 	 	 	 	 	 	    std::abs(innerPoints[innerIterator+1].getX() - innerPoints[innerIterator].getX());

		if(verticalIncrement == newVerticalIncrement){
			data->addQuad(innerPoints[innerIterator],
						  innerPoints[innerIterator+1],
						  outerPoints[outerIterator+1],
						  outerPoints[outerIterator],
						  type);

			++outerIterator;

		}else{
			/*has change*/
			verticalIncrement = newVerticalIncrement;

			data->addQuad(innerPoints[innerIterator],
						  outerPoints[outerIterator+2],
					 	  outerPoints[outerIterator+1],
					 	  outerPoints[outerIterator],
					 	  type);

			outerIterator += 2;

			data->addQuad(innerPoints[innerIterator],
					 	  innerPoints[innerIterator+1],
					 	  outerPoints[outerIterator+1],
					 	  outerPoints[outerIterator],
					 	  type);

			++outerIterator;
		}

	}

}
void MeshGenerator::generatePoints(std::vector<Point>& vector,
								   Point start,
								   coordinates_t incrementX,
								   coordinates_t incrementY,
								   index_t iterations,
								   coordinates_t fixFirstIncrementX,
								   coordinates_t fixFirstIncrementY){

	coordinates_t X = start.getX();
	coordinates_t Y = start.getY();

	for(index_t i = 0; i < iterations; ++i){
		Point point(X,Y);
		point.setID(data->addPoint(point));
		vector.push_back(point);

		if(i == 0){
			X += incrementX + fixFirstIncrementX;
			Y += incrementY + fixFirstIncrementY;
		}else{
			X += incrementX;
			Y += incrementY;
		}
	}
}


void MeshGenerator::generateFixedPoints(std::vector<Point>& vector,
										Point start,
			   	   	    	 	 	 	coordinates_t increment,
			   	   	    	 	 	 	coordinates_t mainIncrementMultiplier,
			   	   	    	 	 	 	bool fixX,
			   	   	    	 	 	 	index_t iterations,
			   	   	    	 	 	 	index_t globalIteration,
			   	   	    	 	 	 	index_t& fixedPoints){

	coordinates_t X = start.getX();
	coordinates_t Y = start.getY();


	for(index_t i = 0; i < iterations; ++i){

		if(fixX){
			X += mainIncrementMultiplier*(sizes[fixedPoints] + fixingIncremenst[fixedPoints] * ((coordinates_t)globalIteration+1));
			Y += increment;
		}else{
			X += increment;
			Y += mainIncrementMultiplier*(sizes[fixedPoints] + fixingIncremenst[fixedPoints] * ((coordinates_t)globalIteration+1));
		}

		Point point(X,Y);
		point.setID(data->addPoint(point));

		vector.push_back(point);

		++fixedPoints;
	}
}
void MeshGenerator::generateFixedPointsInverseOrder(std::vector<Point>& vector,
						 	 	 	 	 	 	   	Point start,
						 	 	 	 	 	 	   	coordinates_t increment,
						 	 	 	 	 	 	   	coordinates_t mainIncrementMultiplier,
						 	 	 	 	 	 	   	bool fixX,
						 	 	 	 	 	 	   	index_t iterations,
						 	 	 	 	 	 	   	index_t globalIteration,
						 	 	 	 	 	 	   	index_t& fixedPoints){
	coordinates_t X = start.getX();
	coordinates_t Y = start.getY();
	std::vector<Point> temporalStore;

	for(index_t i = iterations-1; i > 0; --i){

		if(fixX){
			X -= mainIncrementMultiplier*(sizes[fixedPoints + i - 1] + fixingIncremenst[fixedPoints + i - 1] * ((coordinates_t)globalIteration+1));
			Y += increment;
		}else{
			Y -= mainIncrementMultiplier*(sizes[fixedPoints + i - 1] + fixingIncremenst[fixedPoints + i - 1] * ((coordinates_t)globalIteration+1));
			X += increment;
		}

		Point point(X,Y);
		point.setID(data->addPoint(point));

		temporalStore.push_back(point);

	}

	for(index_t j = temporalStore.size(); j > 0;--j){
		vector.push_back(temporalStore[j-1]);
	}

	fixedPoints += iterations;
}


void MeshGenerator::initTransitionData(index_t numberOfFixingLayers){
	this->numberOfFixingLayers = numberOfFixingLayers;
}

void MeshGenerator::calculateSizes(){
	coordinates_t xDifference;
	coordinates_t yDifference;

	sizes.clear();
	sizes.push_back(0.0);

	for(index_t  i = 1; i < originalInnerPoints.size(); ++i){
		xDifference = std::abs(originalInnerPoints[i].getX() - originalInnerPoints[i-1].getX());
		yDifference = std::abs(originalInnerPoints[i].getY() - originalInnerPoints[i-1].getY());

		if(xDifference > yDifference){
			sizes.push_back(xDifference);
		}else{
			sizes.push_back(yDifference);
		}
	}
}

void MeshGenerator::calculateBlockSize(){

	coordinates_t xDifference;
	coordinates_t yDifference;
	coordinates_t maxDifference;

	blockSize = 0.0;

	for(index_t  i = 1; i < originalInnerPoints.size(); ++i){
		xDifference = std::abs(originalInnerPoints[i].getX() - originalInnerPoints[i-1].getX());
		yDifference = std::abs(originalInnerPoints[i].getY() - originalInnerPoints[i-1].getY());

		if(xDifference > yDifference){
			maxDifference = xDifference;
		}else{
			maxDifference = yDifference;
		}
		if(blockSize < maxDifference)
			blockSize = maxDifference;
	}

}
void MeshGenerator::calculateFixingLayers(){
	index_t baseFixingLayers = numberOfFixingLayers;

	coordinates_t max_x = 0;
	coordinates_t max_y = 0;

	index_t max_fixing_x = 0;
	index_t max_fixing_y = 0;

	for(index_t  i = 1; i < this->innerPoints.size(); ++i){
		if(max_x < this->innerPoints[i].getX()){
			max_x = this->innerPoints[i].getX();
		}

		if(max_y < this->innerPoints[i].getY()){
			max_y = this->innerPoints[i].getY();
		}
	}

	max_fixing_x = (index_t) ((this->data->getMaxX() - max_x) / blockSize);
	max_fixing_y = (index_t) ((this->data->getMaxY() - max_y) / blockSize);

	if(max_fixing_x < numberOfFixingLayers)
		numberOfFixingLayers = max_fixing_x;

	if(max_fixing_y < numberOfFixingLayers)
		numberOfFixingLayers = max_fixing_y;

}

void MeshGenerator::calculateFixIncrements(){

	coordinates_t fix;
	fixingIncremenst.clear();

	fixingIncremenst.push_back(0.0);

	for(index_t  i = 1; i < sizes.size(); ++i){

		fix = (blockSize - sizes[i])/ ((coordinates_t)numberOfFixingLayers);
		fixingIncremenst.push_back(fix);

	}

	fixingIncremenst.push_back(0.0);
}

index_t MeshGenerator::findTransition(std::vector<Point>& points,index_t start){
	index_t i = start;
	bool changeFound = false;
	bool XgraterThanY = false;

	coordinates_t deltaX;
	coordinates_t deltaY;

	deltaX = points[i+1].getX() - points[i].getX();
	deltaY = points[i+1].getY() - points[i].getY();

	XgraterThanY = std::abs(deltaY) < std::abs(deltaX);

	while(!changeFound && ((i+1) < points.size() )){

		deltaX = points[i+1].getX() - points[i].getX();
		deltaY = points[i+1].getY() - points[i].getY();

		++i;
		changeFound = XgraterThanY ^ (std::abs(deltaY) < std::abs(deltaX));
	}

	return i - 1;
}


void MeshGenerator::addRow(std::vector<Point>& points, coordinates_t increment,ids_t type){

	std::vector<Point> newPoints;
	index_t firstPoint;
	index_t lastPoint;
	index_t iterations;
	index_t start;

	firstPoint = findTransition(points,0);
	lastPoint = findTransition(points,firstPoint+1);

	iterations = firstPoint;
	if(increment < 0.0)
		iterations-=2;


	for(index_t j = 0; j <= iterations;++j){
		newPoints.push_back(points[j]);
	}

	for(index_t j = firstPoint; j <= lastPoint;++j){
		Point newPoint(points[j].getX(),points[j].getY()+increment);
		newPoint.setID(data->addPoint(newPoint));
		newPoints.push_back(newPoint);
	}

	start = lastPoint;
	if(increment < 0.0)
		start+=2;

	for(index_t j = start; j < points.size();++j){
		newPoints.push_back(points[j]);
	}

	/*Build quads*/
	for(index_t j = firstPoint; j < lastPoint;++j){

		if(increment > 0.0){
			data->addQuad(points[j],
						  points[j+1],
						  newPoints[j+2],
						  newPoints[j+1],
						  type);
		}else{
			data->addQuad(newPoints[j+1-2],
					  	  newPoints[j+2-2],
					  	  points[j+1],
						  points[j],
						  type);
		}
	}

	points = newPoints;
}

void MeshGenerator::addInnerColumn(std::vector<Point>& points, coordinates_t increment,bool left,ids_t type){
	std::vector<Point> newPoints;
	index_t firstPoint;
	index_t lastPoint;
	index_t newPointsCounter;
	index_t newPointsOffset;

	firstPoint = findTransition(points,0);
	newPointsCounter = firstPoint;

	if(left){
		lastPoint = firstPoint+1;
		firstPoint = 0;
		newPointsOffset = 0;
	}else{
		firstPoint = findTransition(points,firstPoint+1)+1;
		newPointsOffset = firstPoint-1;
		lastPoint = points.size();
	}

	for(index_t j = 0; j < firstPoint;++j){
		newPoints.push_back(points[j]);
	}

	for(index_t j = 0; j < newPointsCounter;++j){
		Point newPoint(points[j+newPointsOffset].getX()+increment,points[j+newPointsOffset].getY());
		newPoint.setID(data->addPoint(newPoint));
		newPoints.push_back(newPoint);
	}

	for(index_t j = lastPoint; j < points.size();++j){
		newPoints.push_back(points[j]);
	}


	/*Build quads*/
	for(index_t j = 0; j < newPointsCounter;++j){

		if(increment > 0.0){
			data->addQuad(newPoints[j+newPointsOffset],
						  newPoints[j+1+newPointsOffset],
						  points[j+1+newPointsOffset],
						  points[j+newPointsOffset],
						  type);
		}else{
			data->addQuad(newPoints[j+1+newPointsOffset],
						  newPoints[j+newPointsOffset],
					  	  points[j+newPointsOffset],
					  	  points[j+1+newPointsOffset],
						  type);
		}
	}

	points = newPoints;
}

void MeshGenerator::addColumn(std::vector<Point>& points, coordinates_t increment,bool left,ids_t type){
	std::vector<Point> newPoints;
	index_t firstPoint;
	index_t lastPoint;
	index_t newPointsCounter;
	index_t newPointsOffset;

	firstPoint = findTransition(points,0);
	newPointsCounter = firstPoint;

	if(left){
		lastPoint = firstPoint;
		firstPoint = 0;
		newPointsOffset = 0;
	}else{
		firstPoint = findTransition(points,firstPoint+1)+1;
		newPointsOffset = firstPoint-1;
		lastPoint = points.size();
	}

	for(index_t j = 0; j < firstPoint;++j){
		newPoints.push_back(points[j]);
	}

	for(index_t j = 0; j <= newPointsCounter;++j){
		Point newPoint(points[j+newPointsOffset].getX()+increment,points[j+newPointsOffset].getY());
		newPoint.setID(data->addPoint(newPoint));
		newPoints.push_back(newPoint);
	}

	for(index_t j = lastPoint; j < points.size();++j){
		newPoints.push_back(points[j]);
	}


	/*Build quads*/
	for(index_t j = 0; j < newPointsCounter;++j){

		if(increment > 0.0){
			data->addQuad(points[j+newPointsOffset],
						  points[j+1+newPointsOffset],
						  newPoints[j+1+newPointsOffset+1],
						  newPoints[j+newPointsOffset+1],
						  type);
		}else{
			data->addQuad(points[j+newPointsOffset],
						  points[j+1+newPointsOffset],
						  newPoints[j+1+newPointsOffset],
					  	  newPoints[j+newPointsOffset],
						  type);
		}
	}

	points = newPoints;
}

void MeshGenerator::generateTransitionTwoBlocs(std::vector<Point>& points,
											   coordinates_t sizeX,
											   coordinates_t sizeY,
											   std::vector<Point>& outPointsVector,
											   index_t start,
											   index_t end,
											   Point firstPoint,
											   ids_t type,
											   bool inverseOrder){

	Point outerDownPoint;
	Point outerMiddlePoint;
	Point outerUpPoint;

	Point innerDownPoint;
	Point innerMiddlePoint;
	Point innerUpPoint;

	index_t i = start;
	outerDownPoint = firstPoint;

	while(i < end){
		outerMiddlePoint = Point(points[i+2].getX() + sizeX*2.0,points[i+2].getY() + sizeY*2.0);
		outerUpPoint = Point(points[i+4].getX() + sizeX*2.0,points[i+4].getY() + sizeY*2.0);

		innerDownPoint = Point(points[i+1].getX() + sizeX,points[i+1].getY() + sizeY);
		innerMiddlePoint = Point(points[i+2].getX() + sizeX,points[i+2].getY() + sizeY);
		innerUpPoint = Point(points[i+3].getX() + sizeX,points[i+3].getY() + sizeY);


		outerMiddlePoint.setID(this->data->addPoint(outerMiddlePoint));
		outerUpPoint.setID(this->data->addPoint(outerUpPoint));
		innerDownPoint.setID(this->data->addPoint(innerDownPoint));
		innerMiddlePoint.setID(this->data->addPoint(innerMiddlePoint));
		innerUpPoint.setID(this->data->addPoint(innerUpPoint));

		outPointsVector.push_back(outerMiddlePoint);
		outPointsVector.push_back(outerUpPoint);

		if(!inverseOrder){
			this->data->addQuad(outerDownPoint,points[i],points[i+1],innerDownPoint,type);
			this->data->addQuad(points[i+1],points[i+2],innerMiddlePoint,innerDownPoint,type);
			this->data->addQuad(outerDownPoint,innerDownPoint,innerMiddlePoint,outerMiddlePoint,type);
			this->data->addQuad(points[i+2],points[i+3],innerUpPoint,innerMiddlePoint,type);
			this->data->addQuad(outerMiddlePoint,innerMiddlePoint,innerUpPoint,outerUpPoint,type);
			this->data->addQuad(outerUpPoint,innerUpPoint,points[i+3],points[i+4],type);
		}else{
			this->data->addQuad(innerDownPoint,points[i+1],points[i],outerDownPoint,type);
			this->data->addQuad(innerDownPoint,innerMiddlePoint,points[i+2],points[i+1],type);
			this->data->addQuad(outerMiddlePoint,innerMiddlePoint,innerDownPoint,outerDownPoint,type);
			this->data->addQuad(innerMiddlePoint,innerUpPoint,points[i+3],points[i+2],type);
			this->data->addQuad(outerUpPoint,innerUpPoint,innerMiddlePoint,outerMiddlePoint,type);
			this->data->addQuad(points[i+4],points[i+3],innerUpPoint,outerUpPoint,type);
		}
		outerDownPoint = outerUpPoint;
		i += 4;
	}

}


void MeshGenerator::addTowBlocksTransition(std::vector<Point>& points,coordinates_t size){

	index_t firstPointChange;
	index_t middlePointChange;
	index_t lastPointChange;
	std::vector<Point> newPoints;

	Point firstPoint;
	Point lastPoint;

	firstPointChange = findTransition(points,0);
	middlePointChange = findTransition(points,firstPointChange+1);
	lastPointChange = findTransition(points,middlePointChange+1);

	/*Add checks??*/

	/*Left blocks*/
	firstPoint = Point(points[0].getX() - size*2,points[0].getY());
	firstPoint.setID(this->data->addPoint(firstPoint));

	newPoints.push_back(firstPoint);

	generateTransitionTwoBlocs(points,
							   -1.0*size,
							   0.0,
							   newPoints,
							   0,
							   firstPointChange,
							   firstPoint,
							   outerCellType);

	/*Corner block*/

	lastPoint = Point(points[firstPointChange].getX() - size*2,points[firstPointChange].getY()+size*2);
	lastPoint.setID(this->data->addPoint(lastPoint));
	firstPoint = Point(points[firstPointChange].getX(),points[firstPointChange].getY()+size*2);
	firstPoint.setID(this->data->addPoint(firstPoint));

	this->data->addQuad(lastPoint,newPoints[newPoints.size()-1],points[firstPointChange],firstPoint,outerCellType);
	newPoints.push_back(lastPoint);
	newPoints.push_back(firstPoint);

	/*Upper blocks*/

	generateTransitionTwoBlocs(points,
							   0.0,
							   1.0*size,
							   newPoints,
							   firstPointChange,
							   middlePointChange,
							   firstPoint,
							   outerCellType);

	/*Corner block*/

	lastPoint = Point(points[middlePointChange].getX() + size*2,points[middlePointChange].getY()+size*2);
	lastPoint.setID(this->data->addPoint(lastPoint));
	firstPoint = Point(points[middlePointChange].getX() + size*2,points[middlePointChange].getY());
	firstPoint.setID(this->data->addPoint(firstPoint));

	this->data->addQuad(lastPoint,newPoints[newPoints.size()-1],points[middlePointChange],firstPoint,outerCellType);
	newPoints.push_back(lastPoint);
	newPoints.push_back(firstPoint);

	/*Right blocks*/

	generateTransitionTwoBlocs(points,
							   1.0*size,
							   0.0,
							   newPoints,
							   middlePointChange,
							   lastPointChange,
							   firstPoint,
							   outerCellType);

	points = newPoints;
}


void MeshGenerator::addInnerLayer(std::vector<Point>& points,coordinates_t size,bool hasPadding){

	index_t firstPointChange;
	index_t middlePointChange;
	index_t pointsSize;
	coordinates_t padding;

	firstPointChange = findTransition(points,0);
	middlePointChange = findTransition(points,firstPointChange+1);
	pointsSize = points.size();

	this->outerPoints.clear();
	this->outerPoints = points;
	this->innerPoints.clear();

	if(hasPadding){
		padding = this->outerPoints[this->outerPoints.size()-1].getX() - this->outerPoints[0].getX();
		padding = (padding - (((middlePointChange-firstPointChange)-2)*size))/2;
	}else{
		padding = size;
	}


	generatePoints(this->innerPoints,
				   this->outerPoints[0] + Point(padding,0.0),
				   0.0,
				   size,
				   firstPointChange);

	generatePoints(this->innerPoints,
				   this->innerPoints[this->innerPoints.size()-1] + Point(1.0*size,0.0),
				   size,
				   0.0,
				   (middlePointChange-firstPointChange)-2);

	generatePoints(this->innerPoints,
				   this->innerPoints[this->innerPoints.size()-1] + Point(0.0,-1.0*size),
				   0.0,
				   -1.0*size,
				   (pointsSize - middlePointChange)-2);

	generateQuads(innerCellType);


}
void MeshGenerator::addOuterLayer(std::vector<Point>& points,coordinates_t size){

}


void MeshGenerator::addInnerTransition(std::vector<Point>& points,coordinates_t size,std::vector<Point>& outPoints){

	index_t firstPointChange;
	index_t middlePointChange;

	index_t offset = 0;

	Point cornerPointOne;
	Point cornerPointTwo;
	Point downPointOne;
	Point downPointTwo;
	Point middlePoint;
	Point upperPoint;

	firstPointChange = findTransition(points,0);
	middlePointChange = findTransition(points,firstPointChange+1);
	outPoints.clear();

	if(((firstPointChange-2)%4)>0){
		offset = 2;
		downPointOne = points[0] + Point(size,0.0);
		downPointTwo = points[0] + Point(2.0*size,0.0);
		middlePoint = points[1] + Point(size,0.0);
		upperPoint = points[2] + Point(2.0*size,0.0);

		downPointOne.setID(this->data->addPoint(downPointOne));
		downPointTwo.setID(this->data->addPoint(downPointTwo));
		middlePoint.setID(this->data->addPoint(middlePoint));
		upperPoint.setID(this->data->addPoint(upperPoint));

		outPoints.push_back(downPointTwo);
		outPoints.push_back(upperPoint);

		this->data->addQuad(points[0],downPointOne,middlePoint,points[1],innerCellType);
		this->data->addQuad(downPointOne,downPointTwo,upperPoint,middlePoint,innerCellType);
		this->data->addQuad(points[1],middlePoint,upperPoint,points[2],innerCellType);
	}else{
		upperPoint = points[0] + Point(2.0*size,0.0);
		upperPoint.setID(this->data->addPoint(upperPoint));
		outPoints.push_back(upperPoint);
	}

	/*left colum*/
	generateTransitionTwoBlocs(points,
							   1.0*size,
							   0.0,
							   outPoints,
							   offset,
							   firstPointChange-2,
							   upperPoint,
							   innerCellType,true);


	/*left corner*/
	cornerPointOne = points[firstPointChange] + Point(size,-1.0*size);
	cornerPointOne.setID(this->data->addPoint(cornerPointOne));
	cornerPointTwo = outPoints[outPoints.size()-1];



	this->data->addQuad(points[firstPointChange],points[firstPointChange-1],cornerPointOne,points[firstPointChange+1],innerCellType);
	this->data->addQuad(points[firstPointChange-1],points[firstPointChange-2],cornerPointTwo,cornerPointOne,innerCellType);
	this->data->addQuad(points[firstPointChange+2],points[firstPointChange+1],cornerPointOne,cornerPointTwo,innerCellType);



	/*upper row*/
	generateTransitionTwoBlocs(points,
							   0.0,
							   -1.0*size,
							   outPoints,
							   firstPointChange+2,
							   middlePointChange-2,
							   cornerPointTwo,
							   innerCellType,true);


	/*right corner*/
	cornerPointOne = points[middlePointChange] + Point(-1.0*size,-1.0*size);
	cornerPointTwo = outPoints[outPoints.size()-1];
	cornerPointOne.setID(this->data->addPoint(cornerPointOne));
	//The data already has this point
	//cornerPointTwo.setID(this->data->addPoint(cornerPointTwo));



	this->data->addQuad(points[middlePointChange],points[middlePointChange-1],cornerPointOne,points[middlePointChange+1],innerCellType);
	this->data->addQuad(points[middlePointChange-1],points[middlePointChange-2],cornerPointTwo,cornerPointOne,innerCellType);
	this->data->addQuad(points[middlePointChange+2],points[middlePointChange+1],cornerPointOne,cornerPointTwo,innerCellType);

	/*right  column*/
	generateTransitionTwoBlocs(points,
							   -1.0*size,
							   0.0,
							   outPoints,
							   middlePointChange+2,
							   points.size()-1-offset,
							   cornerPointTwo,
							   innerCellType,true);

	if(((middlePointChange-2)%4)>0){
		downPointOne = points[points.size()-1] + Point(-1.0*size,0.0);
		downPointTwo = points[points.size()-1] + Point(-2.0*size,0.0);
		middlePoint = points[points.size()-2] + Point(-1.0*size,0.0);
		upperPoint = outPoints[outPoints.size()-1];

		downPointOne.setID(this->data->addPoint(downPointOne));
		downPointTwo.setID(this->data->addPoint(downPointTwo));
		middlePoint.setID(this->data->addPoint(middlePoint));

		outPoints.push_back(downPointTwo);

		this->data->addQuad(points[points.size()-1],points[points.size()-2],middlePoint,downPointOne,innerCellType);
		this->data->addQuad(downPointTwo,downPointOne,middlePoint,upperPoint,innerCellType);
		this->data->addQuad(points[points.size()-3],upperPoint,middlePoint,points[points.size()-2],innerCellType);
	}
}

coordinates_t MeshGenerator::calculateAverageSize(std::vector<Point>& vector){

	coordinates_t sum = 0.0;
	coordinates_t count = (coordinates_t)(vector.size()-1);
	coordinates_t result = 0.0;
	coordinates_t diferenceX;
	coordinates_t diferenceY;

	if(count > 0){

		for(index_t i = 1; i < vector.size();++i)
		{
			diferenceX = std::abs(vector[i].getX()-vector[i-1].getX());
			diferenceY = std::abs(vector[i].getY()-vector[i-1].getY());

			if(diferenceX > diferenceY){
				sum += diferenceX;
			}else{
				sum += diferenceY;
			}
		}

		result = sum / count;
	}

	return result;
}

coordinates_t MeshGenerator::getMinSize(std::vector<Point>& vector){
	coordinates_t result = 0.0;
	bool firstIteration = true;
	coordinates_t diferenceX;
	coordinates_t diferenceY;
	coordinates_t diference;

	for(index_t i = 1; i < vector.size();++i)
	{
		diferenceX = std::abs(vector[i].getX()-vector[i-1].getX());
		diferenceY = std::abs(vector[i].getY()-vector[i-1].getY());

		if(diferenceX > diferenceY){
			diference = diferenceX;
		}else{
			diference = diferenceY;
		}

		if((diference < result)||firstIteration){
			firstIteration = false;
			result = diference;
		}
	}

	return result;
}

void MeshGenerator::verticalEndQuads(){
	index_t verticalPoints;
	index_t middleQuadsStart;
	coordinates_t size;
	verticalPoints = (this->innerPoints.size()-1)/2;

	this->outerPoints = this->innerPoints;
	this->innerPoints.clear();

	size = this->outerPoints[1].getY() - this->outerPoints[0].getY();

	/*Generate the points in the middle*/
	for(index_t i = 0; i < verticalPoints-1;++i)
	{
		Point newPoint;
		newPoint = this->outerPoints[i] + Point(size,0.0);
		newPoint.setID(this->data->addPoint(newPoint));

		this->innerPoints.push_back(newPoint);
	}

	middleQuadsStart = 0;
	/*Generate the left quads */
	for(index_t i = 0; i < verticalPoints-2;++i)
	{
		this->data->addQuad(this->outerPoints[i+1],this->outerPoints[i],this->innerPoints[i],this->innerPoints[i+1],outerCellType);
		middleQuadsStart = i + 1;
	}
	/*Generate the middle quads*/
	this->data->addQuad(this->outerPoints[middleQuadsStart+2],
						this->outerPoints[middleQuadsStart+1],
						this->outerPoints[middleQuadsStart],
						this->innerPoints[middleQuadsStart],
						innerCellType);

	this->data->addQuad(this->innerPoints[middleQuadsStart],
						this->outerPoints[middleQuadsStart+4],
						this->outerPoints[middleQuadsStart+3],
						this->outerPoints[middleQuadsStart+2],
						innerCellType);

	for(index_t i = 0; i < verticalPoints-2;++i)
	{
		this->data->addQuad(this->innerPoints[middleQuadsStart+4+i],
							this->outerPoints[middleQuadsStart-i],
							this->outerPoints[middleQuadsStart-i-1],
							this->innerPoints[middleQuadsStart+4+i+1],
							innerCellType);
	}
}

void MeshGenerator::horizontalEndQuads(){
	index_t horizontalPoints;
	coordinates_t size;

	/*-4 because the corners, only count the middle points*/
	horizontalPoints = (this->innerPoints.size()-4);

	this->outerPoints = this->innerPoints;
	this->innerPoints.clear();

	size = this->outerPoints[1].getY() - this->outerPoints[0].getY();

	/*Generate the points in the middle*/
	for(index_t i = 0; i < horizontalPoints;++i)
	{
		Point newPoint;
		newPoint = this->outerPoints[i] - Point(0.0,size);
		newPoint.setID(this->data->addPoint(newPoint));

		this->innerPoints.push_back(newPoint);
	}

	/*Generate the left quad*/
	this->data->addQuad(this->innerPoints[0],
						this->outerPoints[2],
						this->outerPoints[1],
						this->outerPoints[0],
						innerCellType);


	/*Generate the middle quads */
	for(index_t i = 0; i < horizontalPoints-1;++i)
	{
		this->data->addQuad(this->innerPoints[i],
							this->innerPoints[i+1],
							this->outerPoints[i+1],
							this->outerPoints[i],
							innerCellType);
	}

	/*Generate the right quad*/
	this->data->addQuad(this->innerPoints[horizontalPoints-1],
						this->outerPoints[this->outerPoints.size()-1],
						this->outerPoints[this->outerPoints.size()-2],
						this->outerPoints[this->outerPoints.size()-3],
						innerCellType);

}

void MeshGenerator::initTypes(ids_t outerCell,ids_t innerCell){
	this->outerCellType = outerCell;
	this->innerCellType = innerCell;
}

void MeshGenerator::setUpBoundaries(){
	index_t firstPointChange;
	index_t middlePointChange;

	firstPointChange = findTransition(this->innerPoints,0);
	middlePointChange = findTransition(this->innerPoints,firstPointChange+1);

	for(index_t i = 0; i <= firstPointChange; ++i){
		this->data->addLeftBoundaryPoint(this->innerPoints[i].getID());
	}

	for(index_t i = middlePointChange; i < this->innerPoints.size(); ++i){
		this->data->addRigthBoundaryPoints(this->innerPoints[i].getID());
	}

}
