/*
 * ParameterReader.cpp
 *
 *  Created on: 2/1/2015
 *      Author: juan
 */

#include "ParameterReader.h"
#include <iostream>

int ParameterReader::error = 0;
int ParameterReader::parameterMode = 1;
int ParameterReader::fileMode = 2;
int ParameterReader::helpMode = 3;

std::string ParameterReader::inner_radius_str = "Inner radius";
std::string ParameterReader::outer_radius_str = "Outer radius";
std::string ParameterReader::cell_iterations_str = "Cell iterations";
std::string ParameterReader::sigma_inner_cell_str = "Sigma inner cell";
std::string ParameterReader::sigma_outer_cell_str = "Sigma outer cell";
std::string ParameterReader::sigma_membrane_str = "Sigma membrane";
std::string ParameterReader::potential_density_str = "Potential Density";
std::string ParameterReader::max_iterations_str = "Max iterations";
std::string ParameterReader::tolerance_str = "Tolerance";
std::string ParameterReader::output_path_str = "Output path";
std::string ParameterReader::funyoung_in_str = "Young in";
std::string ParameterReader::funyoung_membrane_str = "Young membrane";
std::string ParameterReader::funyoung_out_str = "Young out";
std::string ParameterReader::funposs_in_str = "Poisson in";
std::string ParameterReader::funposs_membrane_str = "Poisson membrane";
std::string ParameterReader::funposs_out_str = "Poisson out";
std::string ParameterReader::deformation_rep_str = "Deformation repetitions";

std::string ParameterReader::cellX_center_str = "Cell center X";
std::string ParameterReader::cellY_center_str = "Cell center Y";
std::string ParameterReader::membrane_layers_str = "Membrane layers";

std::string ParameterReader::file_sufix_str = "File sufix";

ParameterReader::ParameterReader():argc(0),argv(0) {

	addParameter("-ir",ParameterReader::inner_radius_str,"9.995");
	addParameter("-or",ParameterReader::outer_radius_str,"10.0");
	addParameter("-ci",ParameterReader::cell_iterations_str,"40");
	addParameter("-is",ParameterReader::sigma_inner_cell_str,"0.00002");
	addParameter("-os",ParameterReader::sigma_outer_cell_str,"0.000015");
	addParameter("-ms",ParameterReader::sigma_membrane_str,"0.000000000005");
	addParameter("-v",ParameterReader::potential_density_str,"10.0");
	addParameter("-i",ParameterReader::max_iterations_str,"100000");
	addParameter("-t",ParameterReader::tolerance_str,"0.000001");
	addParameter("-path",ParameterReader::output_path_str,"./output/");
	addParameter("-iy",ParameterReader::funyoung_in_str,"0.0001");
	addParameter("-my",ParameterReader::funyoung_membrane_str,"0.1");
	addParameter("-oy",ParameterReader::funyoung_out_str,"0.0001");
	addParameter("-ip",ParameterReader::funposs_in_str,"0.49");
	addParameter("-mp",ParameterReader::funposs_membrane_str,"0.49");
	addParameter("-op",ParameterReader::funposs_out_str,"0.49");
	addParameter("-dr",ParameterReader::deformation_rep_str,"2");

	addParameter("-ml",ParameterReader::membrane_layers_str,"2");

	addParameter("-sufix",ParameterReader::file_sufix_str,"");

}

ParameterReader::~ParameterReader() {
}

int ParameterReader::process(int argc, char **argv){
	int mode;
	int iteration;
	this->argc = argc;
	this->argv = argv;

	mode = readMode();

	if((mode == ParameterReader::parameterMode) && (argc > 1)){
		iteration = 0;

		while(processParameter(iteration)){
			++iteration;
		}
	}else{
		if((mode == ParameterReader::fileMode) && (argc < 3)){
			mode = ParameterReader::error;
		}
	}

	return mode;
}

int ParameterReader::readMode(){

	int mode = ParameterReader::error;

	if(this->argc > 1){
		std::string param(this->argv[1]);
		if(param.compare("-p") == 0)
			mode = ParameterReader::parameterMode;

		if(param.compare("-f") == 0)
			mode = ParameterReader::fileMode;

		if(param.compare("-h") == 0)
			mode = ParameterReader::helpMode;
	}else{
		mode = ParameterReader::parameterMode;
	}
	return mode;
}

bool ParameterReader::processParameter(int pos){

	bool success = true;
	bool found = false;
	int iteration = 0;

	/* First parameter is the binary name
	 * second parameter is the mode*/

	/*test if the parameter has a value*/
	if(((pos + 1)+2) > argc){
		success = false;
	}else{
		while(!found && (iteration < prefixes.size())){
			found = readParameter(iteration,pos);
			++iteration;
		}
	}

	return success;
}

bool ParameterReader::readParameter(int parameter,int pos){
	bool success = false;

	if(prefixes[parameter].compare(argv[pos + 2]) == 0){
		success = true;
		parameters[savingNames[parameter]] = std::string(argv[pos + 1 + 2]);
	}

	return success;
}

void ParameterReader::addParameter(const std::string& prefix,const std::string& savingName,const std::string& defaultValue){
	prefixes.push_back(prefix);
	savingNames.push_back(savingName);
	parameters[savingName] = defaultValue;
}

void ParameterReader::initData(Data& data,CellGenerator& generator,MeshGenerator& generatorSquare,FortranWrapper& solver){

    data.setPotentialDensity(std::stod(parameters[ParameterReader::potential_density_str]));

    std::cout << parameters[ParameterReader::sigma_inner_cell_str] << " - " << parameters[ParameterReader::sigma_outer_cell_str] << " - " << parameters[ParameterReader::sigma_membrane_str] << std::endl;
    data.setSigmas(std::stod(parameters[ParameterReader::sigma_inner_cell_str]),
    			   std::stod(parameters[ParameterReader::sigma_outer_cell_str]),
    			   std::stod(parameters[ParameterReader::sigma_membrane_str]));

    data.setPath(parameters[ParameterReader::output_path_str]);

    data.setRadius(std::stod(parameters[ParameterReader::inner_radius_str]),
			 	   std::stod(parameters[ParameterReader::outer_radius_str]));

    generator.initCellGeometrics(std::stod(parameters[ParameterReader::inner_radius_str]),
    							 std::stod(parameters[ParameterReader::outer_radius_str]));

    generator.initCellParameters(std::stoi(parameters[ParameterReader::membrane_layers_str]),
    							 std::stoi(parameters[ParameterReader::cell_iterations_str]));

    solver.setParameters(std::stoi(parameters[ParameterReader::max_iterations_str]),
    					 std::stod(parameters[ParameterReader::tolerance_str]),
    					 std::stoi(parameters[ParameterReader::membrane_layers_str]));


    data.initDeformationIterations(std::stoi(parameters[ParameterReader::deformation_rep_str]));

    data.initDeformationPoss(std::stod(parameters[ParameterReader::funposs_in_str]),
    						 std::stod(parameters[ParameterReader::funposs_membrane_str]),
    						 std::stod(parameters[ParameterReader::funposs_out_str]));

    data.initDeformationYoung(std::stod(parameters[ParameterReader::funyoung_in_str]),
    						  std::stod(parameters[ParameterReader::funyoung_membrane_str]),
    						  std::stod(parameters[ParameterReader::funyoung_out_str]));


    data.setFileSufix(parameters[ParameterReader::file_sufix_str]);

}

void ParameterReader::printUsage(){

	for(int i = 0; i < prefixes.size();++i){
		std::cout<< prefixes[i] << " : " << savingNames[i] << std::endl;
	}

}
void ParameterReader::printParametersInUse(){

	for(int i = 0; i < prefixes.size();++i){
		std::cout<< savingNames[i] << " : " << parameters[savingNames[i]] << std::endl;
	}

}
