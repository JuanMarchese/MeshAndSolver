/*
 * Data.cpp
 *
 *  Created on: 9/7/2014
 *      Author: juan
 */

#include "Data.h"
#include "../Problem/Types.h"
#include <cmath>
#include <fstream>



Data::Data():pointID(0),
			 quadID(0),
			 sigma_inner_cell(0),
			 sigma_outer_cell(0),
			 sigma_membrane(0),
			 minimunDifference(100),
			 type_inner_cell(0),
			 type_outer_cell(0),
			 type_membrane(0),
			 potential_density(0),
			 maxZ(0),
			 maxR(0),
			 minZ(0),
			 minR(0),
			 centerZ(0),
			 centerR(0),
			 funyoung_in(0),
			 funyoung_membrane(0),
			 funyoung_out(0),
			 funposs_in(0),
			 funposs_membrane(0),
			 funposs_out(0),
		     deformation_rep(0),
			 inner_radius(0),
			 outer_radius(0){
}

Data::~Data() {
}

std::map<ids_t,Point>& Data::getPoints(){
	return points;
}
std::map<ids_t,Quad>& Data::getQuads(){
	return quads;
}

std::map<ids_t,std::vector<ids_t> >& Data::getPointToQuads(){
	return pointToQuads;
}

void Data::setPath(const std::string& value){
	path = value;
}
std::string Data::getPath(){
	return this->path;
}

coordinates_t Data::getMaxZ(){
	return maxZ;
}
coordinates_t Data::getMaxR(){
	return maxR;
}
coordinates_t Data::getMinZ(){
	return minZ;
}
coordinates_t Data::getMinR(){
	return minR;
}

void Data::setFileSufix(const std::string& sufix){
	this->fileSufix = sufix;
}

std::string Data::getFileSufix(){
	return this->fileSufix;
}

ids_t Data::addPoint(Point point){
	point.setID(pointID++);
	points[point.getID()] = point;


	coordinates_t z = point.getZ();
	coordinates_t r = point.getR();

	if(z < minZ)
		minZ = z;

	if(z > maxZ)
		maxZ = z;

	if(r < minR)
		minR = r;

	if(r > maxR)
		maxR = r;

	return point.getID();
}
ids_t Data::addQuad(Quad quad,ids_t type){
	quad.setID(quadID++);
	quad.setType(type);
	quads[quad.getID()] = quad;

	updateMinimunDifference(quad);

	return quad.getID();
}

ids_t Data::addQuad(const Point& point1,const Point& point2,const Point& point3,const Point& point4,ids_t type){
	return addQuad(Quad(point1,point2,point3,point4),type);
}

void Data::setBoundaries(coordinates_t minZ,coordinates_t maxZ,coordinates_t minR,coordinates_t maxR){
	this->minZ = minZ;
	this->minR = minR;
	this->maxZ = maxZ;
	this->maxR = maxR;

	this->centerZ = maxZ / 2.0;
	this->centerR = 0.0;
}

void Data::saveInFile(const std::string& Zfile,
					  const std::string& Rfile,
					  const std::string& QUADfile,
					  const std::string& MaterialFile,
					  const std::string& BoundaryFile,
					  const std::string& SolutionFile,
					  const std::string& MembranePotential){

	std::ofstream z;
	std::ofstream r;
	std::ofstream quad;
	std::ofstream material;
	std::ofstream boundary;
	std::ofstream solution;
	std::ofstream membrane;
	ids_t maxID = 0;
	std::map<ids_t,coordinates_t> zMap;
	std::map<ids_t,coordinates_t> rMap;

	z.open(addPath(Zfile).c_str(), std::ios::out | std::ios::trunc);
	r.open(addPath(Rfile).c_str(), std::ios::out | std::ios::trunc);
	quad.open(addPath(QUADfile).c_str(), std::ios::out | std::ios::trunc);
	material.open(addPath(MaterialFile).c_str(), std::ios::out | std::ios::trunc);
	boundary.open(addPath(BoundaryFile).c_str(), std::ios::out | std::ios::trunc);
	solution.open(addPath(SolutionFile).c_str(), std::ios::out | std::ios::trunc);
	membrane.open(addPath(MembranePotential).c_str(), std::ios::out | std::ios::trunc);

	z.precision(15);
	r.precision(15);

	for(std::map<ids_t,Point>::iterator it = this->points.begin();it != this->points.end();++it){
		zMap[it->first] = it->second.getZ();
		rMap[it->first] = it->second.getR();

		if(it->second.getID() > maxID){
			maxID = it->second.getID();
		}
	}

	for(ids_t i = 0; i <= maxID;++i){
		z << zMap[i] << std::endl;
		r << rMap[i] << std::endl;
	}
	maxID = points.size();

	for(ids_t i = 0; i < maxID;++i){
		solution << this->solution[i] << std::endl;
	}

	for(ids_t i = 0;i < (ids_t)this->quads.size();++i){
		quad << this->quads[i].getPoint(0).getID() + 1 << " "
			 << this->quads[i].getPoint(1).getID() + 1 << " "
			 << this->quads[i].getPoint(2).getID() + 1 << " "
			 << this->quads[i].getPoint(3).getID() + 1 << " "
			 << std::endl;

		material << i + 1 << " "
				 << this->quads[i].getType()
				 << std::endl;

	}

	for(std::vector<ids_t>::iterator i = upperBoundaryPoints.begin(); i != upperBoundaryPoints.end(); ++i){
		boundary << *i + 1 << " " << 1 << std::endl;
	}

	for(std::vector<ids_t>::iterator i = downBoundaryPoints.begin(); i != downBoundaryPoints.end(); ++i){
		boundary << *i + 1 << " " << 2 << std::endl;
	}

	for(ids_t i = 0; i < membraneOuterPoints.size();++i){
		membrane << this->solution[i] << std::endl;
	}

}

void Data::saveDeformation(const std::string& file,int iteration){

	if(iteration < 1)
	{
		for(int j = 0; j < membraneLayersPoints[0].size(); ++j)
		{
			outerMembranePoints[membraneLayersPoints[0][j]] = points[membraneLayersPoints[0][j]];
		}
	}
	else
	{
		std::ofstream deformationFile;
		deformationFile.open(addPath(file).c_str(), std::ios::out | std::ios::trunc);
		deformationFile.precision(15);

		for(int j = 0; j < membraneLayersPoints[0].size(); ++j)
		{
			coordinates_t tita;
			coordinates_t z;
			coordinates_t r;

			Point pointTemp;
			pointTemp = outerMembranePoints[membraneLayersPoints[0][j]];
			z = pointTemp.getZ();
			r = pointTemp.getR();

			if(z-centerZ > 0.0)
				tita= atan(r/(z-centerZ));
			else
			{
				if((z - centerZ) == 0)
					tita = PI * 0.5;
				else
					tita = atan(r/(z - centerZ))+PI;

			}

			deformationFile << PI - tita;
			deformationFile << "\t";
			deformationFile << z;
			deformationFile << "\t";
			deformationFile << r;
			deformationFile << "\t";
			deformationFile << points[membraneLayersPoints[0][j]].getZ();
			deformationFile << "\t";
			deformationFile << points[membraneLayersPoints[0][j]].getR();
			deformationFile << std::endl;
		}

		deformationFile.close();
	}

}

void Data::saveTension(const std::string& file,int iteration){
	if(iteration >= 1)
	{
		std::ofstream tensionFile;
		tensionFile.open(addPath(file).c_str(), std::ios::out | std::ios::trunc);
		tensionFile.precision(15);

		std::vector<std::vector<coordinates_t> > outer = tension[0];

		for(int j = 0; j < outer.size(); ++j)
		{
			coordinates_t tita;
			coordinates_t z;
			coordinates_t r;

			r = outer[j][0];
			z = outer[j][1];

			if(z-centerZ > 0.0)
				tita= atan(r/(z-centerZ));
			else
			{
				if((z - centerZ) == 0)
					tita = PI * 0.5;
				else
					tita = atan(r/(z - centerZ))+PI;

			}

			tensionFile << PI - tita;

			tensionFile << "\t";
			tensionFile << outer[j][0];//Z

			tensionFile << "\t";
			tensionFile << outer[j][1];//R

			tensionFile << "\t";
			tensionFile << outer[j][2];

			tensionFile << "\t";
			tensionFile << outer[j][3];

			tensionFile << "\t";
			tensionFile << outer[j][4];

			tensionFile << "\t";
			tensionFile << outer[j][5];

			tensionFile << "\t";
			tensionFile << outer[j][6];

			tensionFile << "\t";
			tensionFile << outer[j][7];
			
			tensionFile << std::endl;
		}

		tensionFile.close();
	}
}

void Data::saveMembraneFile(const std::string& membrane){

	std::ofstream membraneFile;

	membraneFile.open(addPath(membrane).c_str(), std::ios::out | std::ios::trunc);
	membraneFile.precision(15);

	coordinates_t tita;
	coordinates_t z;
	coordinates_t r;
	coordinates_t potential = getPotential();
	Point pointTemp;


	for(int j = 0; j < membraneLayersPoints[0].size(); ++j)
	{
		pointTemp = points[membraneLayersPoints[0][j]];
		z = pointTemp.getZ();
		r = pointTemp.getR();

		if(z-centerZ > 0.0)
			tita= atan(r/(z-centerZ));
		else
		{
			if((z - centerZ) == 0)
				tita = PI * 0.5;
			else
				tita = atan(r/(z - centerZ))+PI;

		}

		membraneFile << PI - tita;
		for(int i = membraneLayersPoints.size() - 1; i >= 0 ; --i)
		{

			membraneFile << "\t";
			membraneFile << this->solution[membraneLayersPoints[i][j]] - (potential/2);
		}
		membraneFile << std::endl;
	}

}



void Data::saveVTK(const std::string& file){

	std::ofstream vtkFile;


	ids_t maxID = 0;
	std::map<ids_t,coordinates_t> zMap;
	std::map<ids_t,coordinates_t> rMap;

	vtkFile.open(addPath(file).c_str(), std::ios::out | std::ios::trunc);
	vtkFile.precision(15);

	vtkFile << "# vtk DataFile Version 3.0" << std::endl;
	vtkFile << "vtk output" << std::endl;
	vtkFile << "ASCII" << std::endl;
	vtkFile << "DATASET POLYDATA" << std::endl;
	vtkFile << "POINTS " << this->points.size() << " float" << std::endl;



	for(std::map<ids_t,Point>::iterator it = this->points.begin();it != this->points.end();++it){
		zMap[it->first] = it->second.getZ();
		rMap[it->first] = it->second.getR();

		if(it->second.getID() > maxID){
			maxID = it->second.getID();
		}
	}

	for(ids_t i = 0; i <= maxID;++i){
		vtkFile << rMap[i];
		vtkFile << " ";
		vtkFile << zMap[i];
		vtkFile << " ";
		vtkFile << 0.0;
		vtkFile << std::endl;
	}
	vtkFile << "POLYGONS " << this->quads.size() << " " << this->quads.size() * 5 << std::endl;


	for(ids_t i = 0;i < (ids_t)this->quads.size();++i){
		vtkFile << 4 << " "
				<< this->quads[i].getPoint(0).getID() << " "
				<< this->quads[i].getPoint(1).getID() << " "
				<< this->quads[i].getPoint(2).getID() << " "
				<< this->quads[i].getPoint(3).getID() << " "
				<< std::endl;

	}

	vtkFile << "POINT_DATA " << this->points.size() << std::endl;
	vtkFile << "SCALARS nodal float" << std::endl;
	vtkFile << "LOOKUP_TABLE default" << std::endl;

	for(ids_t i = 0; i <= maxID;++i){
		vtkFile << this->solution[i]<< std::endl;
	}

}

void Data::saveFEM(const std::string& femFile,const std::string& inputFile){

	std::ofstream fem;
	std::ofstream input;
	ids_t id = 0;
	ids_t maxID = 0;
	ids_t type1 = 0;
	ids_t type2 = 0;
	ids_t type3 = 0;


	fem.open(addPath(femFile).c_str(), std::ios::out | std::ios::trunc);
	input.open(addPath(inputFile).c_str(), std::ios::out | std::ios::trunc);

	std::map<ids_t,coordinates_t> zMap;
	std::map<ids_t,coordinates_t> rMap;


	fem << "*COORDINATES" << std::endl;


	for(std::map<ids_t,Point>::iterator it = this->points.begin();it != this->points.end();++it){
		zMap[it->first] = it->second.getZ();
		rMap[it->first] = it->second.getR();

		if(it->second.getID() > maxID){
			maxID = it->second.getID();
		}
	}

	fem << (maxID + 1) << std::endl;

	for(ids_t i = 0; i <= maxID;++i){
		ids_t j = i+1;
		fem << j << "\t";
		fem << zMap[i] << "\t";
		fem << rMap[i] << std::endl;
	}

	fem << "*ELEMENT_GROUPS" << std::endl;
	fem << " 3" << std::endl;

	for(ids_t i = 0;i < (ids_t)this->quads.size();++i){
		if(this->quads[i].getType() == 1){
			type1++;
		}
		if(this->quads[i].getType() == 2){
			type2++;
		}
		if(this->quads[i].getType() == 3){
			type3++;
		}
	}

	fem << "1 "<< type1 <<" Tri3" << std::endl;
	fem << "2 "<< type2 <<" Tri3" << std::endl;
	fem << "3 "<< type3 <<" Tri3" << std::endl;
	fem << "*INCIDENCES" << std::endl;


	for(ids_t i = 0;i < (ids_t)this->quads.size();++i){
		if(this->quads[i].getType() == 1){
			id++;
			fem << id << "\t"
				<< this->quads[i].getPoint(0).getID() + 1 << "\t"
				<< this->quads[i].getPoint(1).getID() + 1 << "\t"
				<< this->quads[i].getPoint(2).getID() + 1 << "\t"
				<< this->quads[i].getPoint(3).getID() + 1 << std::endl;
		}
	}

	for(ids_t i = 0;i < (ids_t)this->quads.size();++i){
		if(this->quads[i].getType() == 2){
			id++;
			fem << id << "\t"
				<< this->quads[i].getPoint(3).getID() + 1 << "\t"
				<< this->quads[i].getPoint(2).getID() + 1 << "\t"
				<< this->quads[i].getPoint(1).getID() + 1 << "\t"
				<< this->quads[i].getPoint(0).getID() + 1 << std::endl;
		}
	}

	for(ids_t i = 0;i < (ids_t)this->quads.size();++i){
		if(this->quads[i].getType() == 3){
			id++;
			fem << id << "\t"
				<< this->quads[i].getPoint(3).getID() + 1 << "\t"
				<< this->quads[i].getPoint(2).getID() + 1 << "\t"
				<< this->quads[i].getPoint(1).getID() + 1 << "\t"
				<< this->quads[i].getPoint(0).getID() + 1 << std::endl;
		}
	}

	input << "Problema: celula" << std::endl;
	input << "data_entrada" << std::endl;
	input << "dimen:  2" << std::endl;
	input << "modo:  4" << std::endl;

	input << "opcion: 1"<< std::endl;
	input << "archivo3: sistema.dat" << std::endl;
	input << "archivo1: " << femFile << std::endl;
	input << "nodpel: 4" << std::endl;

	input << "sigint: " << sigma_inner_cell << std::endl;
	input << "sigext: " << sigma_outer_cell << std::endl;
	input << "sigmem: " << sigma_membrane << std::endl;

	input << "permit: 1.0"<< std::endl;

	input << "Densidad de potencial: " << potential_density << std::endl;
	input << "Frecuencia: 0.0" << std::endl;
	input << "radioExt: " << outer_radius << std::endl;
	input << "radioInt: " << inner_radius << std::endl;
	input << "end_data" << std::endl;



	input << "dirichV: " << upperBoundaryPoints.size() << std::endl;

	for(unsigned i = 0; i < upperBoundaryPoints.size(); ++i)
	{
		input << upperBoundaryPoints[i] + 1 << std::endl;
	}

	input << "dirichT: " << downBoundaryPoints.size() << std::endl;

	for(unsigned i = 0; i < downBoundaryPoints.size(); ++i)
	{
		input << downBoundaryPoints[i] + 1 << std::endl;
	}

	input << "end_dataCC" << std::endl;




}


std::vector<ids_t> Data::getLeftBoundaryPoints(){
	return this->upperBoundaryPoints;
}
std::vector<ids_t> Data::getRigthBoundaryPoints(){
	return this->downBoundaryPoints;
}
void Data::addLeftBoundaryPoint(ids_t id){
	this->upperBoundaryPoints.push_back(id);
}
void Data::addRigthBoundaryPoints(ids_t id){
	this->downBoundaryPoints.push_back(id);
}

ids_t Data::getNodesPerElement(){
	return 4; /*Their are quads!!!*/
}
coordinates_t Data::getSigma_inner_cell(){
	return sigma_inner_cell;
}
coordinates_t Data::getSigma_outer_cell(){
	return sigma_outer_cell;
}
coordinates_t Data::getSigma_membrane(){
	return sigma_membrane;
}
coordinates_t Data::getPotentialDensity(){
	return potential_density;
}
coordinates_t Data::getPotential(){
	//Density is in kV/m,  maxX is in um
	return this->potential_density * this->maxZ * (1000.0 / 1000000.0);
}

void Data::setSigmas(coordinates_t inner_cell,coordinates_t outer_cell,coordinates_t membrane){
	sigma_inner_cell = inner_cell;
	sigma_outer_cell = outer_cell;
	sigma_membrane = membrane;
}
void Data::setPotentialDensity(coordinates_t potential_density){
	this->potential_density = potential_density;
}
void Data::initTypes(ids_t  innerCell,ids_t outerCell,ids_t membrane){
	type_inner_cell = innerCell;
	type_outer_cell = outerCell;
	type_membrane = membrane;
}
ids_t Data::getTypeInnerCell(){
	return type_inner_cell;
}
ids_t Data::getTypeOuterCell(){
	return type_outer_cell;
}
ids_t Data::getTypeMembrane(){
	return type_membrane;
}
void Data::addSolutionPoint(ids_t pointID,coordinates_t value){
	this->solution[pointID] = value;
}
std::map<ids_t,coordinates_t>& Data::getSolution(){
	return solution;
}
std::map<ids_t,std::vector<std::vector<coordinates_t> > >& Data::getTension(){
	return tension;
}
void Data::AddMembraneOuterPoints(ids_t point){
	membraneOuterPoints.push_back(point);
}
void Data::AddMembraneInnerPoints(ids_t point){
	membraneInnerPoints.push_back(point);
}
std::vector<ids_t> Data::getMembraneOuterPoints(){
	return membraneOuterPoints;
}
std::vector<ids_t> Data::getMembraneInnerPoints(){
	return membraneInnerPoints;
}
std::string Data::addPath(const std::string& file){
	std::string newPath = this->path;
	std::string sufixPath = "";
	size_t pointPos;

	bool fileHasDot = file.find('.') != file.npos;

	newPath.append(file);
	sufixPath.append(fileSufix);

	pointPos = newPath.rfind('.');

	if(fileHasDot)
		newPath.insert(pointPos,sufixPath);
	else
		newPath.append(sufixPath);

	return newPath;
}
coordinates_t Data::getMinimunDifference(){
	return minimunDifference;
}
void Data::updateMinimunDifference(Quad quad){
	coordinates_t minZ,maxZ;
	coordinates_t minR,maxR;
	coordinates_t newDifferenceZ;
	coordinates_t newDifferenceR;

	minZ = quad.getPoint(0).getZ();
	maxZ = quad.getPoint(0).getZ();
	minR = quad.getPoint(0).getR();
	maxR = quad.getPoint(0).getR();

	coordinates_t z;
	coordinates_t r;
	for(unsigned i = 1; i < 4 ; ++i)
	{
		z = quad.getPoint(i).getZ();
		r = quad.getPoint(i).getR();

		if(z < minZ)
			minZ = z;

		if(z > maxZ)
			maxZ = z;

		if(r < minR)
			minR = r;

		if(r > maxR)
			maxR = r;
	}

	newDifferenceZ =  fabs(maxZ - minZ);
	newDifferenceR =  fabs(maxR - minR);

	if(newDifferenceZ < minimunDifference)
		minimunDifference = newDifferenceZ;

	if(newDifferenceR < minimunDifference)
		minimunDifference = newDifferenceR;

}

void Data::initDeformationYoung(coordinates_t in,coordinates_t mem,coordinates_t out){
	this->funyoung_in = in;
	this->funyoung_membrane = mem;
	this->funyoung_out = out;
}
void Data::initDeformationPoss(coordinates_t in,coordinates_t mem,coordinates_t out){
	this->funposs_in = in;
	this->funposs_membrane = mem;
	this->funposs_out = out;
}
void Data::initDeformationIterations(ids_t value){
	this->deformation_rep = value;
}

coordinates_t Data::getFunyoung_in(){
	return funyoung_in;
}
coordinates_t Data::getFunyoung_membrane(){
	return funyoung_membrane;
}
coordinates_t Data::getFunyoung_out(){
	return funyoung_out;
}
coordinates_t Data::getFunposs_in(){
	return funposs_in;
}
coordinates_t Data::getFunposs_membrane(){
	return funposs_membrane;
}
coordinates_t Data::getFunposs_out(){
	return funposs_out;
}

ids_t Data::getDeformation_rep(){
	return deformation_rep;
}

void Data::setRadius(coordinates_t inner,coordinates_t outer){
	inner_radius = inner;
	outer_radius = outer;
}

void Data::addMembraneLayer(const std::vector<Point>& points){

	std::vector<ids_t> newVector;

	for(int i = 0; i < points.size();++i){
		newVector.push_back(points[i].getID());
	}

	membraneLayersPoints.push_back(newVector);
}
void Data::addMembraneQuadsIDs(ids_t quad ,ids_t layer){
	membraneLayersQuads[quad] = layer;
}
std::map<ids_t,ids_t>& Data::getMembraneLayersQuads(){
	return membraneLayersQuads;
}

void Data::addMembraneLayerFront(const std::vector<Point>& points){
	std::vector<ids_t> newVector;
	std::vector<std::vector<ids_t> > newVectorContainer;

	for(int i = 0; i < points.size();++i){
		newVector.push_back(points[i].getID());
	}

	newVectorContainer.push_back(newVector);

	for(int i = 0; i < membraneLayersPoints.size();++i){
		newVectorContainer.push_back(membraneLayersPoints[i]);
	}

	membraneLayersPoints = newVectorContainer;

}


coordinates_t Data::getCellZCenter(){
	return centerZ;
}
coordinates_t Data::getCellRCenter(){
	return centerR;
}


void Data::scaleMesh(coordinates_t z,coordinates_t r){

	std::map<ids_t,Point>& points = this->getPoints();

	for(ids_t i = 0;i < points.size();++i )
	{
		points[i].setZ(scale(points[i].getZ()-this->minZ,this->maxZ-this->minZ,z));
		points[i].setR(scale(points[i].getR()-this->minR,this->maxR-this->minR,r));
	}

}

coordinates_t Data::scale(coordinates_t current_value,coordinates_t max_value,coordinates_t scaling_factor){

	coordinates_t scaling = 1 + (pow((1 - (current_value / max_value)),5) * scaling_factor);

    return current_value / scaling;

}
