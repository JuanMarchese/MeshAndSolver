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
								zCenter(0),
								rCenter(0),
							    blockSize(0),
							    blockZIterations(0),
							    blockRIterations(0),
							    numberOfLayers(0),
							    numberOfFixingLayers(0),
							    data(0),
							    outerCellType(0),
							    innerCellType(0) {
}

MeshGenerator::~MeshGenerator() {
}

void MeshGenerator::initBlocksSizes(coordinates_t blockSize,
									index_t blockZIterations,
									index_t blockRIterations,
									index_t numberOfLayers){
	this->blockSize = blockSize;
	this->blockZIterations = blockZIterations;
	this->blockRIterations = blockRIterations;
	this->numberOfLayers = numberOfLayers;
}
void MeshGenerator::initCenterPosition(coordinates_t zCenter,coordinates_t rCenter){
	this->zCenter = zCenter;
	this->rCenter = rCenter;
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
	coordinates_t z;
	coordinates_t r;

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
	z = zCenter - blockSize*(blockZIterations/2+1);
	r = rCenter;

	/*Generate the outer fixing layers points*/
	processOuterFixingLayer(z,r);

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
		addInnerLayer(this->innerPoints,this->innerPoints[1].getR()- this->innerPoints[0].getR(),false);
	}

	side = findTransition(this->innerPoints,0);
	up = findTransition(this->innerPoints,side+1) - side;

	if(side > up){
		horizontalEndQuads();
	}else{
		verticalEndQuads();
	}

}

void MeshGenerator::processOuterFixingLayer(coordinates_t zStart,coordinates_t rStart){

	index_t fixedPoints;

	for(index_t i = 0; i < numberOfFixingLayers ; ++i){
		this->outerPoints.clear();
		fixedPoints = 0;

		generateFixedPoints(this->outerPoints,Point(zStart,rStart),
					   	    0.0,1.0,false,blockRIterations,i,fixedPoints);
		generatePoints(this->outerPoints,
					   this->outerPoints[0] + Point(0.0,blockSize*blockRIterations),
					   0.0,blockSize,i + 1);
		generatePoints(this->outerPoints,
					   this->outerPoints[this->outerPoints.size()-1] + Point(blockSize,0.0),
					   blockSize,0.0,i);
		generateFixedPointsInverseOrder(this->outerPoints,Point(zCenter,(blockRIterations+i)*blockSize),
										0.0,1.0,true,blockZIterations/2+1,i,fixedPoints);

		Point point(zCenter,(blockRIterations+i)*blockSize);
		point.setID(data->addPoint(point));
		this->outerPoints.push_back(point);
		--fixedPoints;

		generateFixedPoints(this->outerPoints,Point(zCenter,(blockRIterations+i)*blockSize),
					   	    0.0,1.0,true,blockZIterations/2,i,fixedPoints);
		generatePoints(this->outerPoints,
					   Point(zCenter+blockSize*(blockZIterations/2+1),this->outerPoints[this->outerPoints.size()-1].getR()),
					   blockSize,0.0,i+1);
		generatePoints(this->outerPoints,
					   this->outerPoints[this->outerPoints.size()-1] + Point(0.0,-1.0*blockSize),
					   0.0,-1.0*blockSize,i);
		generateFixedPointsInverseOrder(this->outerPoints,
										Point(this->outerPoints[this->outerPoints.size()-1].getZ(),rCenter),
										0.0,-1.0,false,blockRIterations+1,i,fixedPoints);
		zStart -= blockSize;

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
			addInnerColumn(this->innerPoints, this->innerPoints[1].getR()- this->innerPoints[0].getR(),true,innerCellType);

		firstPointChange = findTransition(this->innerPoints,0);
		middlePointChange = findTransition(this->innerPoints,firstPointChange+1);
	}


	/*Test if the sides are 2 multiple*/
	if(firstPointChange % 2 != 0)
	{
		addRow(this->innerPoints, this->innerPoints[0].getR()- this->innerPoints[1].getR(),innerCellType);
		firstPointChange = findTransition(this->innerPoints,0);
		middlePointChange = findTransition(this->innerPoints,firstPointChange+1);
	}

	innerPointsTemp.clear();



	/*Add the transition*/
	addInnerTransition(this->innerPoints,innerPoints[1].getR()- innerPoints[0].getR(),innerPointsTemp);

	this->innerPoints = innerPointsTemp;

	/*Add an inner layer*/
	for(index_t i = 0; i < iterations;++i)
	{
		addInnerLayer(this->innerPoints,this->innerPoints[1].getR()- this->innerPoints[0].getR(),false);
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


	if((this->innerPoints[0].getZ() + blockSize)<(this->data->getMinZ())){
		retVal = false;
		increment = this->data->getMinZ() - this->innerPoints[0].getZ();
	}else{
		increment = -1.0 * blockSize;
	}

	generatePoints(this->outerPoints,
				   this->innerPoints[0] + Point(increment,0.0),
				   0.0,blockSize,firstPointChange+1);


	oldIncrement = increment;
	if((this->innerPoints[firstPointChange].getR() + blockSize)>(this->data->getMaxR())){
		retVal = false;
		increment = this->data->getMaxR() - this->innerPoints[firstPointChange].getR();
	}else{
		increment = blockSize;
	}

	generatePoints(this->outerPoints,
				   this->outerPoints[firstPointChange] + Point(0.0,increment),
				   blockSize,0.0,middlePointChange-firstPointChange+2,
				   oldIncrement + blockSize,0.0);


	oldIncrement = increment;
	if((this->innerPoints[middlePointChange+2].getZ() + blockSize)>(this->data->getMaxZ())){
			retVal = false;
			increment = this->data->getMaxZ() - this->innerPoints[middlePointChange+2].getZ();
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

	coordinates_t startZ;
	coordinates_t endZ;
	coordinates_t endR;

	index_t iterationsStartZ;
	index_t iterationsEndZ;
	index_t iterationsEndR;

	index_t firstPointChange;
	index_t middlePointChange;

	firstPointChange = findTransition(this->outerPoints,0);
	middlePointChange = findTransition(this->outerPoints,firstPointChange+1);

	startZ = this->outerPoints[0].getZ() - this->data->getMinZ();
	endZ = this->data->getMaxZ() - this->outerPoints[middlePointChange].getZ();
	endR = this->data->getMaxR() - this->outerPoints[firstPointChange].getR();

	iterationsStartZ = (index_t)(startZ / size);
	iterationsEndZ = (index_t)(endZ / size);
	iterationsEndR = (index_t)(endR / size);

	if(iterationsStartZ < iterationsEndZ){
		iterations = iterationsStartZ;
	}else{
		iterations = iterationsEndZ;
	}

	if(iterationsEndR < iterations){
		iterations = iterationsEndR;
	}

	return iterations;
}

void MeshGenerator::finishBoundary(coordinates_t size){
	coordinates_t startZ;
	coordinates_t endZ;
	coordinates_t endR;
	coordinates_t endIncrementZleft;
	coordinates_t endIncrementZright;
	coordinates_t endIncrementR;

	int iterationsStartZ;
	int iterationsEndZ;
	int iterationsEndR;

	int firstPointChange;
	int middlePointChange;

	firstPointChange = findTransition(this->innerPoints,0);
	middlePointChange = findTransition(this->innerPoints,firstPointChange+1);

	startZ = this->innerPoints[0].getZ() - this->data->getMinZ();
	endZ = this->data->getMaxZ() - this->innerPoints[middlePointChange].getZ();
	endR = this->data->getMaxR() - this->innerPoints[firstPointChange].getR();

	iterationsStartZ = (int)(startZ / size);
	iterationsEndZ = (int)(endZ / size);
	iterationsEndR = (int)(endR / size);

	endIncrementZleft = (this->innerPoints[0].getZ() - (size*iterationsStartZ)) - this->data->getMinZ();
	endIncrementZright = this->data->getMaxZ() - (this->innerPoints[middlePointChange].getZ() + (iterationsEndZ*size));
	endIncrementR = this->data->getMaxR() - (this->innerPoints[firstPointChange].getR() + (iterationsEndR*size));

	for(int i = 0 ; i < iterationsStartZ; i++){
		addColumn(this->innerPoints,size*-1.0,true,outerCellType);
	}

	if(endIncrementZleft > 0.0)
		addColumn(this->innerPoints,endIncrementZleft*-1.0,true,outerCellType);


	for(int i = 0 ; i < iterationsEndZ; i++){
		addColumn(this->innerPoints,size,false,outerCellType);
	}
	if(endIncrementZright > 0.0)
		addColumn(this->innerPoints,endIncrementZright,false,outerCellType);


	for(int i = 0 ; i < iterationsEndR; i++){
		addRow(this->innerPoints,size,outerCellType);
	}
	if(endIncrementR > 0.0)
		addRow(this->innerPoints,endIncrementR,outerCellType);
}

/*We need a minimum of 10 blockSize to work properly*/
bool MeshGenerator::thereIsPlaceToContinue(coordinates_t size){
	bool retVal = false;
	coordinates_t startZ;
	coordinates_t endZ;
	coordinates_t endR;

	index_t firstPointChange;
	index_t middlePointChange;

	firstPointChange = findTransition(this->outerPoints,0);
	middlePointChange = findTransition(this->outerPoints,firstPointChange+1);

	startZ = this->outerPoints[0].getZ() - (size * 6.0);
	endZ = this->outerPoints[middlePointChange].getZ() + (size * 6.0);
	endR = this->outerPoints[firstPointChange].getR() + (size * 6.0);

	retVal = (startZ > this->data->getMinZ()) &&
			 (endZ < this->data->getMaxZ()) &&
			 (endR < this->data->getMaxR());

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

	bool verticalIncrement = std::abs(innerPoints[1].getR() - innerPoints[0].getR())
							 >std::abs(innerPoints[1].getZ() - innerPoints[0].getZ());

	bool newVerticalIncrement = verticalIncrement;

	outerIterator = 0;
	for(innerIterator = 0; innerIterator < innerPoints.size() - 1;++innerIterator){

		newVerticalIncrement =  std::abs(innerPoints[innerIterator+1].getR() - innerPoints[innerIterator].getR()) >
	 	 	 	 	 	 	    std::abs(innerPoints[innerIterator+1].getZ() - innerPoints[innerIterator].getZ());

		if(verticalIncrement == newVerticalIncrement){
			data->addQuad(innerPoints[innerIterator+1],
						  innerPoints[innerIterator],
						  outerPoints[outerIterator],
						  outerPoints[outerIterator+1],
						  type);

			++outerIterator;

		}else{
			/*has change*/
			verticalIncrement = newVerticalIncrement;

			data->addQuad(innerPoints[innerIterator],
						  outerPoints[outerIterator],
					 	  outerPoints[outerIterator+1],
					 	  outerPoints[outerIterator+2],
					 	  type);

			outerIterator += 2;

			data->addQuad(innerPoints[innerIterator+1],
					 	  innerPoints[innerIterator],
					 	  outerPoints[outerIterator],
					 	  outerPoints[outerIterator+1],
					 	  type);

			++outerIterator;
		}

	}

}
void MeshGenerator::generatePoints(std::vector<Point>& vector,
								   Point start,
								   coordinates_t incrementZ,
								   coordinates_t incrementR,
								   index_t iterations,
								   coordinates_t fixFirstIncrementZ,
								   coordinates_t fixFirstIncrementR){

	coordinates_t Z = start.getZ();
	coordinates_t R = start.getR();

	for(index_t i = 0; i < iterations; ++i){
		Point point(Z,R);
		point.setID(data->addPoint(point));
		vector.push_back(point);

		if(i == 0){
			Z += incrementZ + fixFirstIncrementZ;
			R += incrementR + fixFirstIncrementR;
		}else{
			Z += incrementZ;
			R += incrementR;
		}
	}
}


void MeshGenerator::generateFixedPoints(std::vector<Point>& vector,
										Point start,
			   	   	    	 	 	 	coordinates_t increment,
			   	   	    	 	 	 	coordinates_t mainIncrementMultiplier,
			   	   	    	 	 	 	bool fixZ,
			   	   	    	 	 	 	index_t iterations,
			   	   	    	 	 	 	index_t globalIteration,
			   	   	    	 	 	 	index_t& fixedPoints){

	coordinates_t Z = start.getZ();
	coordinates_t R = start.getR();


	for(index_t i = 0; i < iterations; ++i){

		if(fixZ){
			Z += mainIncrementMultiplier*(sizes[fixedPoints] + fixingIncremenst[fixedPoints] * ((coordinates_t)globalIteration+1));
			R += increment;
		}else{
			Z += increment;
			R += mainIncrementMultiplier*(sizes[fixedPoints] + fixingIncremenst[fixedPoints] * ((coordinates_t)globalIteration+1));
		}

		Point point(Z,R);
		point.setID(data->addPoint(point));

		vector.push_back(point);

		++fixedPoints;
	}
}
void MeshGenerator::generateFixedPointsInverseOrder(std::vector<Point>& vector,
						 	 	 	 	 	 	   	Point start,
						 	 	 	 	 	 	   	coordinates_t increment,
						 	 	 	 	 	 	   	coordinates_t mainIncrementMultiplier,
						 	 	 	 	 	 	   	bool fixZ,
						 	 	 	 	 	 	   	index_t iterations,
						 	 	 	 	 	 	   	index_t globalIteration,
						 	 	 	 	 	 	   	index_t& fixedPoints){
	coordinates_t Z = start.getZ();
	coordinates_t R = start.getR();
	std::vector<Point> temporalStore;

	for(index_t i = iterations-1; i > 0; --i){

		if(fixZ){
			Z -= mainIncrementMultiplier*(sizes[fixedPoints + i - 1] + fixingIncremenst[fixedPoints + i - 1] * ((coordinates_t)globalIteration+1));
			R += increment;
		}else{
			R -= mainIncrementMultiplier*(sizes[fixedPoints + i - 1] + fixingIncremenst[fixedPoints + i - 1] * ((coordinates_t)globalIteration+1));
			Z += increment;
		}

		Point point(Z,R);
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
	coordinates_t zDifference;
	coordinates_t rDifference;

	sizes.clear();
	sizes.push_back(0.0);

	for(index_t  i = 1; i < originalInnerPoints.size(); ++i){
		zDifference = std::abs(originalInnerPoints[i].getZ() - originalInnerPoints[i-1].getZ());
		rDifference = std::abs(originalInnerPoints[i].getR() - originalInnerPoints[i-1].getR());

		if(zDifference > rDifference){
			sizes.push_back(zDifference);
		}else{
			sizes.push_back(rDifference);
		}
	}
}

void MeshGenerator::calculateBlockSize(){

	coordinates_t zDifference;
	coordinates_t rDifference;
	coordinates_t maxDifference;

	blockSize = 0.0;

	for(index_t  i = 1; i < originalInnerPoints.size(); ++i){
		zDifference = std::abs(originalInnerPoints[i].getZ() - originalInnerPoints[i-1].getZ());
		rDifference = std::abs(originalInnerPoints[i].getR() - originalInnerPoints[i-1].getR());

		if(zDifference > rDifference){
			maxDifference = zDifference;
		}else{
			maxDifference = rDifference;
		}
		if(blockSize < maxDifference)
			blockSize = maxDifference;
	}

}
void MeshGenerator::calculateFixingLayers(){
	index_t baseFixingLayers = numberOfFixingLayers;

	coordinates_t max_z = 0;
	coordinates_t max_r = 0;

	index_t max_fixing_z = 0;
	index_t max_fixing_r = 0;

	for(index_t  i = 1; i < this->innerPoints.size(); ++i){
		if(max_z < this->innerPoints[i].getZ()){
			max_z = this->innerPoints[i].getZ();
		}

		if(max_r < this->innerPoints[i].getR()){
			max_r = this->innerPoints[i].getR();
		}
	}

	max_fixing_z = (index_t) ((this->data->getMaxZ() - max_z) / blockSize);
	max_fixing_r = (index_t) ((this->data->getMaxR() - max_r) / blockSize);

	if(max_fixing_z < numberOfFixingLayers)
		numberOfFixingLayers = max_fixing_z;

	if(max_fixing_r < numberOfFixingLayers)
		numberOfFixingLayers = max_fixing_r;

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
	bool ZgraterThanR = false;

	coordinates_t deltaZ;
	coordinates_t deltaR;

	deltaZ = points[i+1].getZ() - points[i].getZ();
	deltaR = points[i+1].getR() - points[i].getR();

	ZgraterThanR = std::abs(deltaR) < std::abs(deltaZ);

	while(!changeFound && ((i+1) < points.size() )){

		deltaZ = points[i+1].getZ() - points[i].getZ();
		deltaR = points[i+1].getR() - points[i].getR();

		++i;
		changeFound = ZgraterThanR ^ (std::abs(deltaR) < std::abs(deltaZ));
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
		Point newPoint(points[j].getZ(),points[j].getR()+increment);
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
			data->addQuad(points[j+1],
						  points[j],
						  newPoints[j+1],
						  newPoints[j+2],
						  type);
		}else{
			data->addQuad(newPoints[j+2-2],
					  	  newPoints[j+1-2],
					  	  points[j],
						  points[j+1],
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
		Point newPoint(points[j+newPointsOffset].getZ()+increment,points[j+newPointsOffset].getR());
		newPoint.setID(data->addPoint(newPoint));
		newPoints.push_back(newPoint);
	}

	for(index_t j = lastPoint; j < points.size();++j){
		newPoints.push_back(points[j]);
	}


	/*Build quads*/
	for(index_t j = 0; j < newPointsCounter;++j){

		if(increment > 0.0){
			data->addQuad(newPoints[j+1+newPointsOffset],
						  newPoints[j+newPointsOffset],
						  points[j+newPointsOffset],
						  points[j+1+newPointsOffset],
						  type);
		}else{
			data->addQuad(newPoints[j+newPointsOffset],
						  newPoints[j+1+newPointsOffset],
					  	  points[j+1+newPointsOffset],
					  	  points[j+newPointsOffset],
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
		Point newPoint(points[j+newPointsOffset].getZ()+increment,points[j+newPointsOffset].getR());
		newPoint.setID(data->addPoint(newPoint));
		newPoints.push_back(newPoint);
	}

	for(index_t j = lastPoint; j < points.size();++j){
		newPoints.push_back(points[j]);
	}


	/*Build quads*/
	for(index_t j = 0; j < newPointsCounter;++j){

		if(increment > 0.0){
			data->addQuad(points[j+1+newPointsOffset],
						  points[j+newPointsOffset],
						  newPoints[j+newPointsOffset+1],
						  newPoints[j+1+newPointsOffset+1],
						  type);
		}else{
			data->addQuad(points[j+1+newPointsOffset],
						  points[j+newPointsOffset],
						  newPoints[j+newPointsOffset],
					  	  newPoints[j+1+newPointsOffset],
						  type);
		}
	}

	points = newPoints;
}

void MeshGenerator::generateTransitionTwoBlocs(std::vector<Point>& points,
											   coordinates_t sizeZ,
											   coordinates_t sizeR,
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
		outerMiddlePoint = Point(points[i+2].getZ() + sizeZ*2.0,points[i+2].getR() + sizeR*2.0);
		outerUpPoint = Point(points[i+4].getZ() + sizeZ*2.0,points[i+4].getR() + sizeR*2.0);

		innerDownPoint = Point(points[i+1].getZ() + sizeZ,points[i+1].getR() + sizeR);
		innerMiddlePoint = Point(points[i+2].getZ() + sizeZ,points[i+2].getR() + sizeR);
		innerUpPoint = Point(points[i+3].getZ() + sizeZ,points[i+3].getR() + sizeR);


		outerMiddlePoint.setID(this->data->addPoint(outerMiddlePoint));
		outerUpPoint.setID(this->data->addPoint(outerUpPoint));
		innerDownPoint.setID(this->data->addPoint(innerDownPoint));
		innerMiddlePoint.setID(this->data->addPoint(innerMiddlePoint));
		innerUpPoint.setID(this->data->addPoint(innerUpPoint));

		outPointsVector.push_back(outerMiddlePoint);
		outPointsVector.push_back(outerUpPoint);

		if(inverseOrder){
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
	firstPoint = Point(points[0].getZ() - size*2,points[0].getR());
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

	lastPoint = Point(points[firstPointChange].getZ() - size*2,points[firstPointChange].getR()+size*2);
	lastPoint.setID(this->data->addPoint(lastPoint));
	firstPoint = Point(points[firstPointChange].getZ(),points[firstPointChange].getR()+size*2);
	firstPoint.setID(this->data->addPoint(firstPoint));

	this->data->addQuad(firstPoint,points[firstPointChange],newPoints[newPoints.size()-1],lastPoint,outerCellType);
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

	lastPoint = Point(points[middlePointChange].getZ() + size*2,points[middlePointChange].getR()+size*2);
	lastPoint.setID(this->data->addPoint(lastPoint));
	firstPoint = Point(points[middlePointChange].getZ() + size*2,points[middlePointChange].getR());
	firstPoint.setID(this->data->addPoint(firstPoint));

	this->data->addQuad(firstPoint,points[middlePointChange],newPoints[newPoints.size()-1],lastPoint,outerCellType);
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
		padding = this->outerPoints[this->outerPoints.size()-1].getZ() - this->outerPoints[0].getZ();
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

		this->data->addQuad(points[1],middlePoint,downPointOne,points[0],innerCellType);
		this->data->addQuad(middlePoint,upperPoint,downPointTwo,downPointOne,innerCellType);
		this->data->addQuad(points[2],upperPoint,middlePoint,points[1],innerCellType);

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

	this->data->addQuad(points[firstPointChange+1],cornerPointOne,points[firstPointChange-1],points[firstPointChange],innerCellType);
	this->data->addQuad(cornerPointOne,cornerPointTwo,points[firstPointChange-2],points[firstPointChange-1],innerCellType);
	this->data->addQuad(cornerPointTwo,cornerPointOne,points[firstPointChange+1],points[firstPointChange+2],innerCellType);



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


	this->data->addQuad(points[middlePointChange+1],cornerPointOne,points[middlePointChange-1],points[middlePointChange],innerCellType);
	this->data->addQuad(cornerPointOne,cornerPointTwo,points[middlePointChange-2],points[middlePointChange-1],innerCellType);
	this->data->addQuad(cornerPointTwo,cornerPointOne,points[middlePointChange+1],points[middlePointChange+2],innerCellType);

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

		this->data->addQuad(downPointOne,middlePoint,points[points.size()-2],points[points.size()-1],innerCellType);
		this->data->addQuad(upperPoint,middlePoint,downPointOne,downPointTwo,innerCellType);
		this->data->addQuad(points[points.size()-2],middlePoint,upperPoint,points[points.size()-3],innerCellType);
	}
}

coordinates_t MeshGenerator::calculateAverageSize(std::vector<Point>& vector){

	coordinates_t sum = 0.0;
	coordinates_t count = (coordinates_t)(vector.size()-1);
	coordinates_t result = 0.0;
	coordinates_t diferenceZ;
	coordinates_t diferenceR;

	if(count > 0){

		for(index_t i = 1; i < vector.size();++i)
		{
			diferenceZ = std::abs(vector[i].getZ()-vector[i-1].getZ());
			diferenceR = std::abs(vector[i].getR()-vector[i-1].getR());

			if(diferenceZ > diferenceR){
				sum += diferenceZ;
			}else{
				sum += diferenceR;
			}
		}

		result = sum / count;
	}

	return result;
}

coordinates_t MeshGenerator::getMinSize(std::vector<Point>& vector){
	coordinates_t result = 0.0;
	bool firstIteration = true;
	coordinates_t diferenceZ;
	coordinates_t diferenceR;
	coordinates_t diference;

	for(index_t i = 1; i < vector.size();++i)
	{
		diferenceZ = std::abs(vector[i].getZ()-vector[i-1].getZ());
		diferenceR = std::abs(vector[i].getR()-vector[i-1].getR());

		if(diferenceZ > diferenceR){
			diference = diferenceZ;
		}else{
			diference = diferenceR;
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

	size = this->outerPoints[1].getR() - this->outerPoints[0].getR();

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
		this->data->addQuad(this->outerPoints[i],
							this->outerPoints[i+1],
							this->innerPoints[i+1],
							this->innerPoints[i],
							outerCellType);
		middleQuadsStart = i + 1;
	}
	/*Generate the middle quads*/
	this->data->addQuad(this->outerPoints[middleQuadsStart],
						this->outerPoints[middleQuadsStart+1],
						this->outerPoints[middleQuadsStart+2],
						this->innerPoints[middleQuadsStart],
						innerCellType);

	this->data->addQuad(this->innerPoints[middleQuadsStart],
						this->outerPoints[middleQuadsStart+2],
						this->outerPoints[middleQuadsStart+3],
						this->outerPoints[middleQuadsStart+4],
						innerCellType);

	for(index_t i = 0; i < verticalPoints-2;++i)
	{
		this->data->addQuad(this->innerPoints[middleQuadsStart+4+i+1],
							this->outerPoints[middleQuadsStart-i-1],
							this->outerPoints[middleQuadsStart-i],
							this->innerPoints[middleQuadsStart+4+i],
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

	size = this->outerPoints[1].getR() - this->outerPoints[0].getR();

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
						this->outerPoints[0],
						this->outerPoints[1],
						this->outerPoints[2],
						innerCellType);


	/*Generate the middle quads */
	for(index_t i = 0; i < horizontalPoints-1;++i)
	{
		this->data->addQuad(this->innerPoints[i+1],
							this->innerPoints[i],
							this->outerPoints[i],
							this->outerPoints[i+1],
							innerCellType);
	}

	/*Generate the right quad*/
	this->data->addQuad(this->innerPoints[horizontalPoints-1],
						this->outerPoints[this->outerPoints.size()-3],
						this->outerPoints[this->outerPoints.size()-2],
						this->outerPoints[this->outerPoints.size()-1],
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
