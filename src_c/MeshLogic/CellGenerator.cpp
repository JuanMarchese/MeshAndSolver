/*
 * CellGenerator.cpp
 *
 *  Created on: 13/7/2014
 *      Author: juan
 */

#include "CellGenerator.h"
#include <cmath>

CellGenerator::CellGenerator():innerRadius(0),
							   outerRadius(0),
							   cellZCenter(0),
							   cellRCenter(0),
							   cellLayers(0),
							   cellIterations(0),
							   boundaryLayers(0),
							   boundaryLinesSeparation(0),
							   innerLayers(0),
							   innerLinesSeparation(0),
							   data(0),
							   inCellType(0),
							   outerCellType(0),
							   innerCellType(0) {
}

CellGenerator::~CellGenerator() {
}

void CellGenerator::initData(Data* data){
	this->data=data;
}
void CellGenerator::initCellGeometrics(coordinates_t innerRadius,
									   coordinates_t outerRadius){
	this->innerRadius=innerRadius;
	this->outerRadius=outerRadius;
}
void CellGenerator::initCellParameters(index_t cellLayers,
									   index_t cellIterations){
	this->cellLayers=cellLayers;
	this->cellIterations=cellIterations;
}
void CellGenerator::initBoundaryParameters(index_t boundaryLayers,
										   coordinates_t boundaryLinesSeparation){
	this->boundaryLayers=boundaryLayers;
	this->boundaryLinesSeparation=boundaryLinesSeparation;
}
void CellGenerator::initInnerParameters(index_t innerLayers,
										coordinates_t innerLinesSeparation){
	this->innerLayers=innerLayers;
	this->innerLinesSeparation=innerLinesSeparation;
}
std::vector<Point>& CellGenerator::getOuterPoints(){
	return realOuterPoints;
}
std::vector<Point>& CellGenerator::getInnerPoints(){
	return realInnerPoints;
}
void CellGenerator::generateRadiusVector(){

	coordinates_t actualRadius = innerRadius;
	coordinates_t radiusIncrement = (outerRadius - innerRadius)/(cellLayers);

	radius.push_back(actualRadius);

	for(index_t i = 0;i <=cellLayers;++i){

		actualRadius += radiusIncrement;
		radius.push_back(actualRadius);
	}

}

void CellGenerator::generateOuterRadiusVector(){

	coordinates_t coordinateRadius = outerRadius + boundaryLinesSeparation;
	coordinates_t actualRadius;
	coordinates_t expConstant;

	coordinates_t endRealRadiusMult = 20.0;

	expConstant = log(endRealRadiusMult-1) / (((coordinates_t)boundaryLayers));


	for(int i = -1.0*boundaryLayers;i <=(int)boundaryLayers;i+=2){
		actualRadius = coordinateRadius * (1+exp(((coordinates_t)(i))*expConstant));

		outerRadiusVector.push_back(actualRadius);

		coordinateRadius += boundaryLinesSeparation;
	}

}


void CellGenerator::process(){

	cellZCenter = data->getCellZCenter();
	cellRCenter = data->getCellRCenter();

	generateAngles();
	generateRadiusVector();
	generateOuterRadiusVector();
	generateInnerRadiusVector();

	generateCellCirclePoints(innerPoints,0);
	firstInnerPoints = innerPoints;


	data->addMembraneLayer(firstInnerPoints);

	for(index_t i = 1;i <=cellLayers;++i){
		generateCellCirclePoints(outerPoints,i);
		generateQuads(inCellType,true,i-1);
		innerPoints = outerPoints;
		data->addMembraneLayer(outerPoints);
	}

	for(index_t i = 1;i <=boundaryLayers;++i){
		generateBoundaryCirclePoints(outerPoints,i);
		generateQuads(outerCellType);
		innerPoints = outerPoints;
	}

	realOuterPoints = outerPoints;
	outerPoints = firstInnerPoints;

	for(index_t i = 1;i <=innerLayers;++i){
		generateInnerCirclePoints(innerPoints,i,false);

		if(!circlesOverlap()){
			addPointsToData(innerPoints);
			generateQuads(innerCellType);
			outerPoints = innerPoints;
		}

	}
	realInnerPoints = innerPoints;

	for(index_t i = 0; i < realOuterPoints.size() ; ++i)
	{
		data->AddMembraneOuterPoints(realOuterPoints[i].getID());
	}

	for(index_t i = 0; i < realInnerPoints.size() ; ++i)
	{
		data->AddMembraneInnerPoints(realInnerPoints[i].getID());
	}

}
void CellGenerator::generateCellCirclePoints(std::vector<Point>& vector,index_t iteration){

	coordinates_t angle = 0.0;

	coordinates_t Z;
	coordinates_t R;

	coordinates_t radius = this->radius[iteration];
	vector.clear();

	for(index_t i = 0; i <= cellIterations; ++i){

		angle = angles[i];
		Z = cellZCenter + cos(angle)*radius;
		R = cellRCenter + sin(angle)*radius;

		Point point(Z,R);

		point.setID(data->addPoint(point));
		vector.push_back(point);
	}

}

void CellGenerator::generateFinalVector(index_t start,
							 	 	 	index_t increment,
							 	 	 	index_t iterations,
										std::vector<Point>& vector,
										std::vector<coordinates_t> vectorZ,
										std::vector<coordinates_t> vectorR,
										coordinates_t zSign,
										coordinates_t rSign,
										bool addToData){

	coordinates_t Z;
	coordinates_t R;
	index_t i = start;

	for(index_t j = 0;j < iterations; ++j){
		Z = cellZCenter + zSign * vectorZ[i];
		R = cellRCenter + rSign * vectorR[i];
		Point point(Z,R);

		if(addToData)
			point.setID(data->addPoint(point));

		vector.push_back(point);
		i += increment;
	}

}

void CellGenerator::generateBoundaryCirclePoints(std::vector<Point>& vector,index_t iteration){

	coordinates_t distanceFromCenter;
	coordinates_t currentRealRadius;

	coordinates_t coordinateRadius;

	coordinates_t angle = 0.0;

	coordinates_t betaAngle;

	coordinates_t Z;
	coordinates_t R;

	std::vector<coordinates_t> Zoriginal;
	std::vector<coordinates_t> Roriginal;


	coordinateRadius = outerRadius + boundaryLinesSeparation * (iteration);
	currentRealRadius = outerRadiusVector[iteration];
	distanceFromCenter = coordinateRadius - currentRealRadius;


	vector.clear();

	for(index_t i=0.0; i <= (cellIterations/4);++i){
		angle = angles[i];
		betaAngle = 3.14159265 - angle
				  - asin( sin ( (3.14159265 - angle))
					* distanceFromCenter / currentRealRadius);

		Z = currentRealRadius * cos(betaAngle) - distanceFromCenter;
		R = currentRealRadius * sin(betaAngle);

		Zoriginal.push_back(Z);
		Roriginal.push_back(R);
	}

	index_t size = Zoriginal.size();

	generateFinalVector(0,1,size,vector,Zoriginal,Roriginal,-1.0,1.0);
	generateFinalVector(size-2,-1,size-1,vector,Roriginal,Zoriginal,1.0,-1.0);
	generateFinalVector(1,1,size-1,vector,Roriginal,Zoriginal,-1.0,-1.0);
	generateFinalVector(size-2,-1,size-1,vector,Zoriginal,Roriginal,1.0,1.0);


}

void CellGenerator::generateAngles(){

	coordinates_t distanceFromCenter;
	coordinates_t perimeter;
	coordinates_t blockSize;
	coordinates_t currentAnge = 0.0;
	coordinates_t edge;
	index_t size;

	angles.clear();
	distanceFromCenter = outerRadius + (boundaryLinesSeparation * boundaryLayers);

	/*left + 2 times up + right*/
	perimeter = distanceFromCenter * 4;

	blockSize = perimeter / (coordinates_t)cellIterations;/*cellIterations is multiple of 4 */

	edge = 0.0;
	for(index_t i = 0; i <= (cellIterations/4);++i)
	{
		angles.push_back(currentAnge);
		edge += blockSize;

		currentAnge = asin(edge / sqrt(edge*edge + distanceFromCenter*distanceFromCenter));
	}
	currentAnge = 3.14159265/2.0;
	size = angles.size();
	for(index_t i = size-1; i > 0;--i)
	{
		angles.push_back(currentAnge-angles[i-1]);
	}

	size = angles.size();
	currentAnge = 3.14159265;
	for(index_t i = size-1; i > 0;--i)
	{
		angles.push_back(currentAnge-angles[i-1]);
	}

}

void CellGenerator::generateQuads(ids_t type,bool saveIDs,ids_t layer){

	for(unsigned i = 0; i < innerPoints.size()-1;++i){
		Quad newQuad(outerPoints[i+1],
					 outerPoints[i],
					 innerPoints[i],
					 innerPoints[i+1]);

		ids_t id = data->addQuad(newQuad,type);

		if(saveIDs){
			data->addMembraneQuadsIDs(id,layer);
		}

	}
}


//TODO cambiar esto para que termine mas derecho
void CellGenerator::generateInnerRadiusVector(){
	coordinates_t coordinateRadius = innerRadius - innerLinesSeparation;
	coordinates_t actualRadius;

	coordinates_t startRealRadiusMult = 1.0;
	coordinates_t endRealRadiusMult = 2.0;

	coordinates_t realRadiusMultiplier = startRealRadiusMult;
	coordinates_t realRadiusIncrement = (endRealRadiusMult - startRealRadiusMult) / innerLayers;


	for(index_t i = 0;i <=innerLayers;++i){
		actualRadius = innerRadius * realRadiusMultiplier;

		innerRadiusVector.push_back(actualRadius);

		realRadiusMultiplier += realRadiusIncrement;
		coordinateRadius -= innerLinesSeparation;
	}
}

void CellGenerator::generateInnerCirclePoints(std::vector<Point>& vector,index_t iteration,bool addPointsToData){
	coordinates_t distanceFromCenter;
	coordinates_t currentRealRadius;

	coordinates_t coordinateRadius;

	coordinates_t angle = 0.0;
	//coordinates_t angleIncrement = 3.14159265/cellIterations;

	coordinates_t betaAngle;

	coordinates_t Z;
	coordinates_t R;

	std::vector<coordinates_t> Zoriginal;
	std::vector<coordinates_t> Roriginal;


	coordinateRadius = innerRadius -(innerLinesSeparation * (iteration));
	currentRealRadius = innerRadiusVector[iteration];
	distanceFromCenter = coordinateRadius - currentRealRadius;


	vector.clear();

	//angle = 0.0;
	for(index_t i = 0; i <= (cellIterations/4.0);++i){
		angle = angles[i];
		betaAngle = 3.14159265 - angle
				  - asin( sin ( (3.14159265 - angle))
					* distanceFromCenter / currentRealRadius);

		Z = currentRealRadius * cos(betaAngle) - distanceFromCenter;
		R = currentRealRadius * sin(betaAngle);

		Zoriginal.push_back(Z);
		Roriginal.push_back(R);
		//angle += angleIncrement;
	}

	index_t size = Zoriginal.size();

	generateFinalVector(0,1,size,vector,Zoriginal,Roriginal,-1.0,1.0,addPointsToData);
	generateFinalVector(size-2,-1,size-1,vector,Roriginal,Zoriginal,1.0,-1.0,addPointsToData);
	generateFinalVector(1,1,size-1,vector,Roriginal,Zoriginal,-1.0,-1.0,addPointsToData);
	generateFinalVector(size-2,-1,size-1,vector,Zoriginal,Roriginal,1.0,1.0,addPointsToData);
}

void CellGenerator::initTypes(ids_t inCell,ids_t outerCell,ids_t innerCell){
	this->inCellType = inCell;
	this->outerCellType = outerCell;
	this->innerCellType = innerCell;
}

bool CellGenerator::circlesOverlap(){
	bool overlap = false;
	index_t size = innerPoints.size();
	index_t i = 0;


	while(!overlap && (i < size)){
		coordinates_t innerRadius = (innerPoints[i].getZ() - cellZCenter)*(innerPoints[i].getZ() - cellZCenter);
		innerRadius += (innerPoints[i].getR() - cellRCenter)*(innerPoints[i].getR() - cellRCenter);

		coordinates_t outerRadius = (outerPoints[i].getZ() - cellZCenter)*(outerPoints[i].getZ() - cellZCenter);
		outerRadius += (outerPoints[i].getR() - cellRCenter)*(outerPoints[i].getR() - cellRCenter);

		overlap = innerRadius >= outerRadius;
		++i;
	}

	return overlap;
}

void CellGenerator::addPointsToData(std::vector<Point>& newPoints){
	for(std::vector<Point>::iterator it = newPoints.begin();it != newPoints.end();++it){
		it->setID(data->addPoint(*it));
	}
}
coordinates_t CellGenerator::getInnerRadius(){
	return innerRadius;
}
coordinates_t CellGenerator::getOuterRadius(){
	return outerRadius;
}
index_t CellGenerator::getCellIterations(){
	return cellIterations;
}
