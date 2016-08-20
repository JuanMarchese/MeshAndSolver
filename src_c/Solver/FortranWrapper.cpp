/*
 * FortranWrapper.cpp
 *
 *  Created on: 16/11/2014
 *      Author: juan
 */

#include "FortranWrapper.h"
#include <stdio.h>
#include <iostream>
#include <cmath>


FortranWrapper::FortranWrapper():data(0),
								 maxIterations(0),
								 tolerance(0),
								 cellLayers(0){
}

FortranWrapper::~FortranWrapper() {
}

void FortranWrapper::setData(Data* data){
	this->data = data;
}

void FortranWrapper::setParameters(ids_t maxIterations,coordinates_t tolerance,ids_t cellLayers){
	this->maxIterations = maxIterations;
	this->tolerance = tolerance;
	this->cellLayers = cellLayers;
}

void FortranWrapper::initSolver(){
	int pointsSize = data->getPoints().size();
	int elemPerNode = 4;
	int quadsSize = data->getQuads().size();

	allocate_memory_(&pointsSize,&elemPerNode,&quadsSize);

	for(ids_t i = 0; i < this->data->getPoints().size()*2; ++i)
	{
		__def_variables_MOD_vec_quieto[i] = -1;
	}

	copyCoordinatesData();
	copyQuadsData();
	copyBoundaryData();
	copyMaterialData();

	__def_variables_MOD_potencial = data->getPotential();
	__def_variables_MOD_tierra = 0.0;

	__def_variables_MOD_sigmaext = data->getSigma_outer_cell();
	__def_variables_MOD_sigmamem = data->getSigma_membrane();
	__def_variables_MOD_sigmaint = data->getSigma_inner_cell();

	__def_variables_MOD_funyoung_in = data->getFunyoung_in();
	__def_variables_MOD_funyoung_membrane = data->getFunyoung_membrane();
	__def_variables_MOD_funyoung_out = data->getFunyoung_out();

	__def_variables_MOD_funposs_in = data->getFunposs_in();
	__def_variables_MOD_funposs_membrane = data->getFunposs_membrane();
	__def_variables_MOD_funposs_out = data->getFunposs_out();


	__def_constantes_MOD_toler = tolerance;
	__def_constantes_MOD_itermax = maxIterations;

}

void FortranWrapper::copyCoordinatesData(){

	std::map<ids_t,Point>& points = this->data->getPoints();
	coordinates_t difference = this->data->getMinimunDifference() / 2.0;

	for(ids_t i = 0;i < points.size();++i )
	{
		__def_variables_MOD_coor_y[i] = points[i].getZ();
		__def_variables_MOD_coor_x[i] = points[i].getR();

	}


}
void FortranWrapper::copyQuadsData(){

	std::map<ids_t,Quad>& quads = this->data->getQuads();
	ids_t size = quads.size();
	for(ids_t i = 0;i < size;++i )
	{
		__def_variables_MOD_conect[0*size + i] = quads[i].getPoint(0).getID() + 1;
		__def_variables_MOD_conect[1*size + i] = quads[i].getPoint(1).getID() + 1;
		__def_variables_MOD_conect[2*size + i] = quads[i].getPoint(2).getID() + 1;
		__def_variables_MOD_conect[3*size + i] = quads[i].getPoint(3).getID() + 1;
	}
}
void FortranWrapper::copyBoundaryData(){
	ids_t size = this->data->getPoints().size();

	for(ids_t i = 0; i < size; ++i)
	{
		__def_variables_MOD_vec_poten[i] = -1;
		__def_variables_MOD_vec_tierra[i] = -1;
	}

	std::vector<ids_t> potential = data->getLeftBoundaryPoints();
	std::vector<ids_t> ground = data->getRigthBoundaryPoints();

	for(ids_t i = 0; i < potential.size();++i){
		__def_variables_MOD_vec_poten[potential[i]] = 2;/*any number different from -1*/
		__def_variables_MOD_vec_quieto[potential[i]*2] = 0;
		__def_variables_MOD_vec_quieto[potential[i]*2+1] = 0;

	}
	for(ids_t i = 0; i < ground.size();++i){
		__def_variables_MOD_vec_tierra[ground[i]] = 2;/*any number different from -1*/
		__def_variables_MOD_vec_quieto[ground[i]*2] = 0;
		__def_variables_MOD_vec_quieto[ground[i]*2+1] = 0;
	}


}
void FortranWrapper::copyMaterialData(){

	for(std::map<ids_t,Quad>::iterator it = this->data->getQuads().begin();
		it != this->data->getQuads().end();
		++it){

			if(it->second.getType() == data->getTypeMembrane()){
				__def_variables_MOD_material[it->first] = 2;
			}else if(it->second.getType() == data->getTypeInnerCell()){
				__def_variables_MOD_material[it->first] = 1;
			}else{/*data->getTypeOuterCell() will be the default*/
				__def_variables_MOD_material[it->first] = 3;
			}
	}
}

void FortranWrapper::solveElectromagnecticField(){

	poisson_();
	//print_memory_();
	for(ids_t i = 0; i < data->getPoints().size();++i){
		data->addSolutionPoint(i,__def_solver_MOD_solucion[i]);
	}

}

void FortranWrapper::solveDeformation(){

	deforma_();
	//print_memory_();
	std::map<ids_t,Point>& points = this->data->getPoints();



	for(ids_t i = 0;i < points.size();++i )
	{
		/*The mesh gener is transposed*/
		points[i].setR((coordinates_t)__def_variables_MOD_coor_x[i]);
		points[i].setZ((coordinates_t)__def_variables_MOD_coor_y[i]);
	}


	std::map<ids_t,Quad>& quads = this->data->getQuads();
	ids_t size = quads.size();
	std::map<ids_t,std::vector<std::vector<coordinates_t> > >& tesnion = this->data->getTension();
	std::map<ids_t,ids_t> membraneLayer = this->data->getMembraneLayersQuads();

	for(int i = 0; i < cellLayers-1;i++)
	{
		std::vector<std::vector<coordinates_t> > empty;
		tesnion[i] = empty;
	}

	for(ids_t i = 0;i < size;++i )
	{
		if(membraneLayer.find(i) != membraneLayer.end())
		{
			std::vector<coordinates_t> temp;

			coordinates_t x = 0.0;
			coordinates_t y = 0.0;
			coordinates_t r = 0.0;
			coordinates_t z = 0.0;
			coordinates_t tita = 0.0;
			coordinates_t nose = 0.0;
			coordinates_t lado = 0.0;

			coordinates_t x_value = quads[i].getPoint(0).getR() - quads[i].getPoint(1).getR();
			coordinates_t y_value = quads[i].getPoint(0).getZ() - quads[i].getPoint(1).getZ();


			lado = sqrt((x_value*x_value)+(y_value*y_value));

			for(ids_t j = 0;j < 4;++j )
			{
				x += quads[i].getPoint(j).getR() / 4.0;
				y += quads[i].getPoint(j).getZ() / 4.0;
				r += __def_variables_MOD_tension[i + j*size + 0*size*4] / 4.0;
				z += __def_variables_MOD_tension[i + j*size + 1*size*4] / 4.0;
				tita += __def_variables_MOD_tension[i + j*size + 2*size*4] / 4.0;
				nose += __def_variables_MOD_tension[i + j*size + 3*size*4] / 4.0;
			}

			temp.push_back(x);
			temp.push_back(y);
			temp.push_back(r);
			temp.push_back(z);
			temp.push_back(tita);
			temp.push_back(nose);
			temp.push_back(lado);

			tesnion[membraneLayer[i]].push_back(temp);

		}


	}

}

void FortranWrapper::destroySolver(){
	de_allocate_memory_();
	de_allocate_external_memory_();
}




void FortranWrapper::solveProblemsIteratively(){

	std::cout << "Iniciando iteraciones, cantidad: " << data->getDeformation_rep() << std::endl;

	char buffer [10];
	std::string originalSufix = data->getFileSufix();

	if(data->getDeformation_rep() < 1){
		data->saveVTK(".vtk");
	}else{
		for(ids_t i = 0; i < data->getDeformation_rep();++i){

				std::string newSufix = originalSufix;

				if(data->getDeformation_rep() > 1){
					snprintf( buffer, 10, "%03d",(int)i);

					if(originalSufix.length() > 0)
						newSufix.append("-");

					newSufix.append(buffer);
				}

				data->setFileSufix(newSufix);

			    this->initSolver();
			    std::cout << "Realizando los calculos";

			    if(data->getDeformation_rep() > 1){
			    	std::cout << " iteracion " << buffer;
			    }

			    std::cout << std::endl;


			    this->solveElectromagnecticField();

			    /*We save the data of the field before the new deformation*/
			    std::cout << "Guardando los resultados";
			    if(data->getDeformation_rep() > 1){
			    	std::cout << " iteracion " << buffer;
			    }
			    std::cout << std::endl;

			    data->saveVTK(".vtk");
			    data->saveDeformation(".dcsv",i);
			    data->saveTension(".tcsv",i);
			    data->saveMembraneFile(".csv");


			    if(i < (data->getDeformation_rep()-1)){

			    	std::cout << "Calculo la deformacion";
					if(data->getDeformation_rep() > 1){
						std::cout << " iteracion " << buffer;
					}
					std::cout << std::endl;

			    	this->solveDeformation();
			    }

			    this->destroySolver();

			}
	}



}


