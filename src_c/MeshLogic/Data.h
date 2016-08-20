/*
 * Data.h
 *
 *  Created on: 9/7/2014
 *      Author: juan
 */

#ifndef DATA_H_
#define DATA_H_

#include <vector>
#include <map>
#include "Point.h"
#include "Quad.h"

class Data {
private:

	ids_t pointID;
	ids_t quadID;

	coordinates_t sigma_inner_cell;
	coordinates_t sigma_outer_cell;
	coordinates_t sigma_membrane;

	coordinates_t minimunDifference;

	ids_t type_inner_cell;
	ids_t type_outer_cell;
	ids_t type_membrane;

	coordinates_t potential_density;

	std::map<ids_t,Point> points;
	std::map<ids_t,Point> outerMembranePoints;

	std::map<ids_t,Quad> quads;

	std::map<ids_t,std::vector<ids_t> > pointToQuads;

	std::map<ids_t,coordinates_t> solution;
	std::map<ids_t,std::vector<std::vector<coordinates_t> > > tension;

	std::vector<ids_t> leftBoundaryPoints;
	std::vector<ids_t> rigthBoundaryPoints;

	std::vector<ids_t> membraneOuterPoints;
	std::vector<ids_t> membraneInnerPoints;

	std::vector<std::vector<ids_t> > membraneLayersPoints;
	std::map<ids_t,ids_t> membraneLayersQuads;

	coordinates_t maxX;
	coordinates_t maxY;
	coordinates_t minX;
	coordinates_t minY;

	std::string path;
	std::string fileSufix;

	coordinates_t centerX;
	coordinates_t centerY;

	coordinates_t funyoung_in;
	coordinates_t funyoung_membrane;
	coordinates_t funyoung_out;
	coordinates_t funposs_in;
	coordinates_t funposs_membrane;
	coordinates_t funposs_out;

	ids_t deformation_rep;

	coordinates_t inner_radius;
	coordinates_t outer_radius;

public:
	Data();
	virtual ~Data();

	void fakeInit();

	void setPath(const std::string& value);
	std::string getPath();

	std::map<ids_t,Point>& getPoints();
	std::map<ids_t,Quad>& getQuads();
	std::map<ids_t,std::vector<ids_t> >& getPointToQuads();

	std::map<ids_t,ids_t>& getMembraneLayersQuads();

	std::vector<ids_t> getLeftBoundaryPoints();
	std::vector<ids_t> getRigthBoundaryPoints();

	void AddMembraneOuterPoints(ids_t point);
	void AddMembraneInnerPoints(ids_t point);

	std::vector<ids_t> getMembraneOuterPoints();
	std::vector<ids_t> getMembraneInnerPoints();

	void addLeftBoundaryPoint(ids_t id);
	void addRigthBoundaryPoints(ids_t id);

	void addSolutionPoint(ids_t pointID,coordinates_t value);
	std::map<ids_t,coordinates_t>& getSolution();
	std::map<ids_t,std::vector<std::vector<coordinates_t> > >& getTension();

	ids_t getNextPointID();
	ids_t getNextQuadID();

	coordinates_t getMaxX();
	coordinates_t getMaxY();
	coordinates_t getMinX();
	coordinates_t getMinY();

	void initTypes(ids_t  innerCell,ids_t outerCell,ids_t membrane);
	ids_t getTypeInnerCell();
	ids_t getTypeOuterCell();
	ids_t getTypeMembrane();


	ids_t getNodesPerElement();
	coordinates_t getSigma_inner_cell();
	coordinates_t getSigma_outer_cell();
	coordinates_t getSigma_membrane();

	coordinates_t getPotentialDensity();
	coordinates_t getPotential();

	void setSigmas(coordinates_t inner_cell,coordinates_t outer_cell,coordinates_t membrane);
	void setPotentialDensity(coordinates_t potential_density);

	void setBoundaries(coordinates_t minX,coordinates_t maxX,coordinates_t minY,coordinates_t maxY);

	ids_t addPoint(Point point);
	ids_t addQuad(Quad quad,ids_t type = 0);
	ids_t addQuad(const Point& point1,const Point& point2,const Point& point3,const Point& point4,ids_t type = 0);

	void saveInFile(const std::string& Xfile,
					const std::string& Yfile,
					const std::string& QUADfile,
					const std::string& MaterialFile,
					const std::string& BoundaryFile,
					const std::string& SolutionFile,
					const std::string& MembranePotential);

	void saveInParavieFile(const std::string& mesh,
						   const std::string& membrane,
						   const std::string& solution);

	void saveVTK(const std::string& file);
	void saveDeformation(const std::string& file,int iteration);
	void saveTension(const std::string& file,int iteration);
	void saveFEM(const std::string& femFile,const std::string& inputFile);
	void saveMembraneFile(const std::string& membrane);


	void setFileSufix(const std::string& sufix);
	std::string getFileSufix();

	coordinates_t getMinimunDifference();

	void initDeformationYoung(coordinates_t in,coordinates_t mem,coordinates_t out);
	void initDeformationPoss(coordinates_t in,coordinates_t mem,coordinates_t out);
	void initDeformationIterations(ids_t value);

	coordinates_t getFunyoung_in();
	coordinates_t getFunyoung_membrane();
	coordinates_t getFunyoung_out();
	coordinates_t getFunposs_in();
	coordinates_t getFunposs_membrane();
	coordinates_t getFunposs_out();

	ids_t getDeformation_rep();
	void setRadius(coordinates_t inner,coordinates_t outer);

	void addMembraneLayer(const std::vector<Point>& points);
	void addMembraneQuadsIDs(ids_t quad,ids_t layer);
	void addMembraneLayerFront(const std::vector<Point>& points);

	coordinates_t getCellXCenter();
	coordinates_t getCellYCenter();

private:
	std::string addPath(const std::string& file);
	void updateMinimunDifference(Quad quad);
};

#endif /* DATA_H_ */
