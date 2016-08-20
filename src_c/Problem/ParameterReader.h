/*
 * ParameterReader.h
 *
 *  Created on: 2/1/2015
 *      Author: juan
 */

#ifndef PARAMETERREADER_H_
#define PARAMETERREADER_H_

#include <string>
#include <vector>
#include <map>

#include "../MeshLogic/Data.h"
#include "../MeshLogic/CellGenerator.h"
#include "../MeshLogic/MeshGenerator.h"
#include "../Solver/FortranWrapper.h"

class ParameterReader {

private:
	std::map<std::string,std::string> parameters;
	std::vector<std::string> prefixes;
	std::vector<std::string> savingNames;

	static std::string inner_radius_str;
	static std::string outer_radius_str;
	static std::string cell_iterations_str;
	static std::string sigma_inner_cell_str;
	static std::string sigma_outer_cell_str;
	static std::string sigma_membrane_str;
	static std::string potential_density_str;
	static std::string max_iterations_str;
	static std::string tolerance_str;
	static std::string output_path_str;
	static std::string funyoung_in_str;
	static std::string funyoung_membrane_str;
	static std::string funyoung_out_str;
	static std::string funposs_in_str;
	static std::string funposs_membrane_str;
	static std::string funposs_out_str;
	static std::string deformation_rep_str;

	static std::string cellX_center_str;
	static std::string cellY_center_str;
	static std::string membrane_layers_str;

	static std::string file_sufix_str;


	int argc;
	char **argv;

public:
	ParameterReader();
	virtual ~ParameterReader();

	int process(int argc, char **argv);
	void initData(Data& data,CellGenerator& generator,MeshGenerator& generatorSquare,FortranWrapper& solver);
	void printUsage();
	void printParametersInUse();

	static int parameterMode;
	static int fileMode;
	static int error;
	static int helpMode;

private:
	int readMode();
	bool processParameter(int pos);
	bool readParameter(int parameter,int pos);
	void addParameter(const std::string& prefix,const std::string& savingName,const std::string& defaultValue);

};

#endif /* PARAMETERREADER_H_ */
